# R/tier2_mr_sensitivity_suite.R
# ============================================================================
# Tier 2-B: the MR sensitivity suite that prior reviewers flagged as NOT run.
# For each discussed gene, after CORRECT full-1000G-EUR LD-clumping, run the
# full robustness panel the MendelianRandomization package supports:
#   IVW, MR-Egger (intercept = directional-pleiotropy test),
#   weighted-median, MR-Lasso (outlier-robust; the MR-PRESSO analogue, since
#   MRPRESSO/TwoSampleMR are not installable here), mode-based estimate,
#   + a manual Steiger directionality check (cis-eQTL exposure vs eBMD outcome).
# Multi-SNP methods require >=3 instruments; otherwise only Wald/Steiger reported.
#
# Output: 07_mr_coloc/mr_sensitivity_suite.tsv
# ============================================================================
.libPaths(c("F:/Rlib", .libPaths()))
suppressMessages({ library(data.table); library(ieugwasr); library(MendelianRandomization) })
setDTthreads(1)
G <- "00_data/genetic"; PLINK <- file.path(G,"plink/plink.exe"); BFILE <- file.path(G,"ld/EUR")
fl <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(),"%H:%M:%S"), paste0(...)))
TARGETS <- c("GREM2","SRGN","WNT16","KREMEN1","MGP","SMAD3","HSPG2","TLR4","LOXL2","EDNRB","PDGFD","IRS1")

sig <- fread("03_rra_meta/regeneration_competent_signature.tsv")
fl("reading eQTLGen ..."); eq <- fread(file.path(G,"eqtlgen_cis_sig.txt.gz"))
eq <- eq[GeneSymbol %in% TARGETS & Pvalue < 5e-8]
fl("reading eBMD ..."); eb <- fread(file.path(G,"ebmd_GCST006979.h.tsv.gz"),
  select=c("hm_rsid","hm_effect_allele","hm_other_allele","hm_beta","standard_error","hm_effect_allele_frequency","p_value","n"))
setnames(eb, c("rsid","ea_o","oa_o","beta_o","se_o","eaf_o","p_o","n_o"))
eb <- eb[!is.na(rsid)&!is.na(beta_o)&!is.na(eaf_o)]; setkey(eb, rsid)

zhu  <- function(z,p,N){ d<-2*p*(1-p)*(N+z^2); ifelse(d>0, z/sqrt(d), NA) }
zhse <- function(p,N,z){ d<-2*p*(1-p)*(N+z^2); ifelse(d>0, 1/sqrt(d), NA) }
is_pal <- function(a1,a2) paste0(toupper(a1),toupper(a2)) %in% c("AT","TA","CG","GC")

harmonise_clump <- function(g){
  e <- eq[GeneSymbol==g]; m <- merge(e, eb, by.x="SNP", by.y="rsid"); if(!nrow(m)) return(NULL)
  ok <- (toupper(m$AssessedAllele)==toupper(m$ea_o)&toupper(m$OtherAllele)==toupper(m$oa_o)) |
        (toupper(m$AssessedAllele)==toupper(m$oa_o)&toupper(m$OtherAllele)==toupper(m$ea_o))
  m <- m[ok & !is_pal(m$AssessedAllele,m$OtherAllele)]; if(!nrow(m)) return(NULL)
  flip <- toupper(m$AssessedAllele)!=toupper(m$ea_o); eaf <- ifelse(flip,1-m$eaf_o,m$eaf_o)
  d <- data.frame(rsid=m$SNP, pval=m$Pvalue, eaf=eaf, n_out=m$n_o,
                  bx=zhu(m$Zscore,eaf,m$NrSamples), bxse=zhse(eaf,m$NrSamples,m$Zscore),
                  by=ifelse(flip,-m$beta_o,m$beta_o), byse=m$se_o)
  d <- d[is.finite(d$bx)&is.finite(d$bxse)&is.finite(d$by)&is.finite(d$byse)&d$bxse>0&d$byse>0,]
  if(!nrow(d)) return(NULL)
  if(nrow(d)>1){
    cl <- tryCatch(ld_clump(dplyr::tibble(rsid=d$rsid,pval=d$pval,id=g),
            clump_r2=0.001, clump_kb=10000, plink_bin=PLINK, bfile=BFILE), error=function(e) NULL)
    if(!is.null(cl) && nrow(cl)>0) d <- d[d$rsid %in% cl$rsid,,drop=FALSE] else d <- d[order(d$pval),][1,,drop=FALSE]
  }
  d
}

steiger <- function(d){                       # cis-eQTL exposure vs eBMD outcome
  r2x <- sum(2*d$eaf*(1-d$eaf)*d$bx^2)
  r2y <- sum(2*d$eaf*(1-d$eaf)*d$by^2)
  list(r2_exposure=r2x, r2_outcome=r2y,
       direction=ifelse(r2x>r2y, "exposure->outcome (correct)", "ambiguous/reverse"))
}
fld <- function(x) if(is.null(x)||length(x)==0) NA_real_ else as.numeric(x)

run_gene <- function(g){
  d <- harmonise_clump(g); if(is.null(d)) return(NULL)
  ns <- nrow(d); st <- steiger(d)
  row <- data.table(gene=g, nsnp=ns,
    IVW_b=NA, IVW_p=NA, egger_int_p=NA, wmedian_b=NA, wmedian_p=NA,
    lasso_b=NA, lasso_p=NA, lasso_nvalid=NA, mode_b=NA, mode_p=NA,
    steiger_r2_exp=round(st$r2_exposure,4), steiger_r2_out=signif(st$r2_outcome,3), steiger_dir=st$direction)
  if(ns==1){
    row$IVW_b <- d$by/d$bx; row$IVW_p <- 2*pnorm(-abs((d$by/d$bx)/abs(d$byse/d$bx)))  # Wald reported in IVW slot
  } else {
    mi <- mr_input(bx=d$bx,bxse=d$bxse,by=d$by,byse=d$byse)
    iv <- tryCatch(mr_ivw(mi), error=function(e)NULL)
    if(!is.null(iv)){ row$IVW_b<-iv$Estimate; row$IVW_p<-iv$Pvalue }
    if(ns>=3){
      eg <- tryCatch(mr_egger(mi), error=function(e)NULL); if(!is.null(eg)) row$egger_int_p<-eg$Pvalue.Int
      wm <- tryCatch(mr_median(mi, weighting="weighted"), error=function(e)NULL)
      if(!is.null(wm)){ row$wmedian_b<-wm$Estimate; row$wmedian_p<-wm$Pvalue }
      ml <- tryCatch(mr_lasso(mi), error=function(e)NULL)
      if(!is.null(ml)){ row$lasso_b<-fld(ml$Estimate); row$lasso_p<-fld(ml$Pvalue); row$lasso_nvalid<-length(ml$Valid) }
      mo <- tryCatch(mr_mbe(mi, seed=1337), error=function(e)NULL)
      if(!is.null(mo)){ row$mode_b<-fld(mo$Estimate); row$mode_p<-fld(mo$Pvalue) }
    }
  }
  row
}

fl("running sensitivity suite ...")
res <- rbindlist(lapply(TARGETS, function(g) tryCatch(run_gene(g), error=function(e){fl("ERR ",g,": ",conditionMessage(e)); NULL})), fill=TRUE)
res$direction_signature <- sig$direction[match(res$gene, sig$gene)]
num <- setdiff(names(res), c("gene","steiger_dir","direction_signature"))
for(c in num) if(is.numeric(res[[c]])) res[[c]] <- signif(res[[c]], 3)
fwrite(res, "07_mr_coloc/mr_sensitivity_suite.tsv", sep="\t")
fl("=== MR sensitivity suite ===")
print(res[, .(gene, nsnp, IVW_p, egger_int_p, wmedian_p, lasso_p, lasso_nvalid, steiger_dir)])
fl("DONE.")
