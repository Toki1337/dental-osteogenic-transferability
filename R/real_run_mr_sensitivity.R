# R/real_run_mr_sensitivity.R
# Targeted multi-SNP sensitivity for the handful of genes actually discussed,
# using CORRECT LD-clumping against the FULL 1000G EUR panel (8.55M variants),
# with proper handling: keep only clump-index SNPs; on clump failure fall back to
# the single lead cis-SNP (never pass all correlated SNPs through). Reports IVW,
# MR-Egger intercept (pleiotropy) where >=3 instruments, and per-gene nsnp so the
# instrument counts can be verified as biologically plausible.
.libPaths(c("F:/Rlib", .libPaths()))
suppressMessages({ library(data.table); library(ieugwasr); library(MendelianRandomization) })
setDTthreads(0)
G <- "00_data/genetic"; PLINK <- file.path(G,"plink/plink.exe"); BFILE <- file.path(G,"ld/EUR")
fl <- function(...) { cat(sprintf("[%s] %s\n", format(Sys.time(),"%H:%M:%S"), paste0(...))); flush.console() }

TARGETS <- c("GREM2","SRGN","WNT16","KREMEN1","MGP","SMAD3","HSPG2",
             "TLR4","LOXL2","EDNRB","PDGFD","IRS1")

sig <- fread("03_rra_meta/regeneration_competent_signature.tsv")
fl("reading eQTLGen ..."); eq <- fread(file.path(G,"eqtlgen_cis_sig.txt.gz"))
eq <- eq[GeneSymbol %in% TARGETS & Pvalue < 5e-8]
fl("reading eBMD ..."); eb <- fread(file.path(G,"ebmd_GCST006979.h.tsv.gz"),
  select=c("hm_rsid","hm_effect_allele","hm_other_allele","hm_beta","standard_error","hm_effect_allele_frequency","p_value","n"))
setnames(eb, c("rsid","ea_o","oa_o","beta_o","se_o","eaf_o","p_o","n_o"))
eb <- eb[!is.na(rsid)&!is.na(beta_o)&!is.na(eaf_o)]; setkey(eb, rsid)

zhu  <- function(z,p,N){ d<-2*p*(1-p)*(N+z^2); ifelse(d>0, z/sqrt(d), NA) }
zhse <- function(p,N,z){ d<-2*p*(1-p)*(N+z^2); ifelse(d>0, 1/sqrt(d), NA) }
is_palindromic <- function(a1,a2){ paste0(toupper(a1),toupper(a2)) %in% c("AT","TA","CG","GC") }

run_gene <- function(g){
  e <- eq[GeneSymbol==g]; m <- merge(e, eb, by.x="SNP", by.y="rsid")
  if(!nrow(m)) return(NULL)
  ok <- (toupper(m$AssessedAllele)==toupper(m$ea_o)&toupper(m$OtherAllele)==toupper(m$oa_o)) |
        (toupper(m$AssessedAllele)==toupper(m$oa_o)&toupper(m$OtherAllele)==toupper(m$ea_o))
  palin <- is_palindromic(m$AssessedAllele, m$OtherAllele)
  m <- m[ok & !palin]; if(!nrow(m)) return(NULL)
  flip <- toupper(m$AssessedAllele)!=toupper(m$ea_o); eaf <- ifelse(flip,1-m$eaf_o,m$eaf_o)
  d <- data.frame(rsid=m$SNP, pval=m$Pvalue, bx=zhu(m$Zscore,eaf,m$NrSamples),
                  bxse=zhse(eaf,m$NrSamples,m$Zscore), by=ifelse(flip,-m$beta_o,m$beta_o), byse=m$se_o)
  d <- d[is.finite(d$bx)&is.finite(d$bxse)&is.finite(d$by)&is.finite(d$byse)&d$bxse>0&d$byse>0,]
  if(!nrow(d)) return(NULL)
  n_pre <- nrow(d)
  if(nrow(d)>1){
    cl <- tryCatch(ld_clump(dplyr::tibble(rsid=d$rsid,pval=d$pval,id=g),
            clump_r2=0.001, clump_kb=10000, plink_bin=PLINK, bfile=BFILE), error=function(e) NULL)
    if(!is.null(cl) && nrow(cl)>0){ d <- d[d$rsid %in% cl$rsid,,drop=FALSE]
    } else { d <- d[order(d$pval),][1,,drop=FALSE] }  # CORRECT fallback: lead SNP only
  }
  ns <- nrow(d)
  out <- data.frame(gene=g, nsnp_preclump=n_pre, nsnp=ns, method=NA, b=NA, se=NA, pval=NA, egger_int_p=NA)
  if(ns==1){ out$method<-"Wald"; out$b<-d$by/d$bx; out$se<-abs(d$byse/d$bx); out$pval<-2*pnorm(-abs(out$b/out$se))
  } else { mi<-mr_input(bx=d$bx,bxse=d$bxse,by=d$by,byse=d$byse); iv<-tryCatch(mr_ivw(mi),error=function(e)NULL)
    if(!is.null(iv)){ out$method<-"IVW"; out$b<-iv$Estimate; out$se<-iv$StdError; out$pval<-iv$Pvalue }
    if(ns>=3){ eg<-tryCatch(mr_egger(mi),error=function(e)NULL); if(!is.null(eg)) out$egger_int_p<-eg$Pvalue.Int } }
  out
}
res <- do.call(rbind, lapply(TARGETS, function(g) tryCatch(run_gene(g), error=function(e){fl("ERR ",g,": ",conditionMessage(e)); NULL})))
res$direction_signature <- sig$direction[match(res$gene, sig$gene)]
fwrite(res, "07_mr_coloc/mr_sensitivity_fullpanel_IVW.tsv", sep="\t")
fl("=== targeted full-panel-clumped IVW/Egger sensitivity ===")
print(res)
fl("DONE sensitivity.")
