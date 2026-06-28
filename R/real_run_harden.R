# R/real_run_harden.R — address reviewer demands: (1) binomial test on the jaw-position
# signal vs chance; (2) proper plink LD-clumping for MR (replacing distance pruning),
# via a pre-extracted small reference panel for speed; (3) FinnGen replication attempt.
.libPaths(c("F:/Rlib", .libPaths()))
options(repos=c(CRAN="https://cloud.r-project.org"))
if(!requireNamespace("R.utils",quietly=TRUE)) install.packages("R.utils")
suppressMessages({ library(data.table); library(ieugwasr); library(MendelianRandomization) })
setDTthreads(0)
G <- "00_data/genetic"; PLINK <- file.path(G,"plink/plink.exe"); BFILE <- file.path(G,"ld/EUR")
fl <- function(...) { cat(sprintf("[%s] %s\n", format(Sys.time(),"%H:%M:%S"), paste0(...))); flush.console() }

# ---------- (1) position signal: binomial test vs 0.5 ----------
jc <- fread("05_projection/jaw_position_context.tsv")
k <- sum(jc$jaw_weak_direction_consistent %in% c(TRUE,"TRUE")); n <- nrow(jc)
bt <- binom.test(k, n, 0.5)
core <- fread("04_modules_wgcna/core_program.tsv")
jc_core <- jc[gene %in% core$gene]
kc <- sum(jc_core$jaw_weak_direction_consistent %in% c(TRUE,"TRUE")); nc <- nrow(jc_core)
btc <- binom.test(kc, nc, 0.5)
pos <- data.frame(set=c("signature","core"), k=c(k,kc), n=c(n,nc),
                  frac=round(c(k/n,kc/nc),3), binom_p=signif(c(bt$p.value,btc$p.value),3))
fwrite(pos, "05_projection/position_binomial_test.tsv", sep="\t")
fl("POSITION binomial: signature ", k,"/",n," (",round(k/n,3),") p=",signif(bt$p.value,3),
   " | core ",kc,"/",nc," p=",signif(btc$p.value,3), "  => NOT significant if p>0.05")

# ---------- (2) MR with proper LD clumping ----------
sig <- fread("03_rra_meta/regeneration_competent_signature.tsv"); cand <- unique(sig$gene)
fl("reading eQTLGen ..."); eq <- fread(file.path(G,"eqtlgen_cis_sig.txt.gz"))
eq <- eq[GeneSymbol %in% cand & Pvalue < 5e-8]
fl("instruments: ", nrow(eq), " across ", length(unique(eq$GeneSymbol)), " genes")
fl("reading eBMD ..."); eb <- fread(file.path(G,"ebmd_GCST006979.h.tsv.gz"),
  select=c("hm_rsid","hm_effect_allele","hm_other_allele","hm_beta","standard_error","hm_effect_allele_frequency","p_value","n"))
setnames(eb, c("rsid","ea_o","oa_o","beta_o","se_o","eaf_o","p_o","n_o")); eb <- eb[!is.na(rsid)&!is.na(beta_o)]; setkey(eb, rsid)

# pre-extract a SMALL reference panel containing only candidate instrument SNPs (fast clumping)
writeLines(unique(eq$SNP), file.path(G,"cand_snps.txt"))
system2(PLINK, c("--bfile", shQuote(BFILE), "--extract", shQuote(file.path(G,"cand_snps.txt")),
                 "--make-bed", "--out", shQuote(file.path(G,"cand_panel"))), stdout=NULL, stderr=NULL)
SMALL <- file.path(G,"cand_panel")
fl("pre-extracted small panel: ", file.exists(paste0(SMALL,".bed")))

zhu  <- function(z,p,N){ d<-2*p*(1-p)*(N+z^2); ifelse(d>0, z/sqrt(d), NA) }
zhse <- function(p,N,z){ d<-2*p*(1-p)*(N+z^2); ifelse(d>0, 1/sqrt(d), NA) }

run_gene <- function(g, outcome=eb, label="eBMD"){
  e <- eq[GeneSymbol==g]; m <- merge(e, outcome, by.x="SNP", by.y="rsid")
  if(!nrow(m)) return(NULL)
  ok <- (toupper(m$AssessedAllele)==toupper(m$ea_o)&toupper(m$OtherAllele)==toupper(m$oa_o)) |
        (toupper(m$AssessedAllele)==toupper(m$oa_o)&toupper(m$OtherAllele)==toupper(m$ea_o))
  m <- m[ok]; if(!nrow(m)) return(NULL)
  flip <- toupper(m$AssessedAllele)!=toupper(m$ea_o); eaf <- ifelse(flip,1-m$eaf_o,m$eaf_o)
  d <- data.frame(rsid=m$SNP, bx=zhu(m$Zscore,eaf,m$NrSamples), bxse=zhse(eaf,m$NrSamples,m$Zscore),
                  by=ifelse(flip,-m$beta_o,m$beta_o), byse=m$se_o, pval=m$Pvalue)
  d <- d[is.finite(d$bx)&is.finite(d$bxse)&is.finite(d$by)&is.finite(d$byse)&d$bxse>0&d$byse>0,]
  if(!nrow(d)) return(NULL)
  if(nrow(d)>1){ cl <- tryCatch(ld_clump(dplyr::tibble(rsid=d$rsid,pval=d$pval,id=g),
                    clump_r2=0.001, clump_kb=10000, plink_bin=PLINK, bfile=SMALL), error=function(e) NULL)
    if(!is.null(cl)) d <- d[d$rsid %in% cl$rsid,,drop=FALSE] }
  ns <- nrow(d); out <- data.frame(gene=g, nsnp=ns, method=NA, b=NA, se=NA, pval=NA, egger_int_p=NA, label=label)
  if(ns==1){ out$method<-"Wald"; out$b<-d$by/d$bx; out$se<-abs(d$byse/d$bx); out$pval<-2*pnorm(-abs(out$b/out$se))
  } else { mi<-mr_input(bx=d$bx,bxse=d$bxse,by=d$by,byse=d$byse); iv<-tryCatch(mr_ivw(mi),error=function(e)NULL)
    if(!is.null(iv)){ out$method<-"IVW"; out$b<-iv$Estimate; out$se<-iv$StdError; out$pval<-iv$Pvalue }
    if(ns>=3){ eg<-tryCatch(mr_egger(mi),error=function(e)NULL); if(!is.null(eg)) out$egger_int_p<-eg$Pvalue.Int } }
  out
}
genes <- unique(eq$GeneSymbol); fl("LD-clumped MR over ", length(genes), " genes ...")
res <- do.call(rbind, lapply(genes, function(g) tryCatch(run_gene(g), error=function(e) NULL)))
res <- res[!is.na(res$b),]; res$padj <- p.adjust(res$pval,"BH")
res$direction_signature <- sig$direction[match(res$gene, sig$gene)]
res <- res[order(res$pval),]
fwrite(res, "07_mr_coloc/mr_results_ldclumped.tsv", sep="\t")
fl("LD-clumped MR: ", nrow(res), " genes; FDR<0.05: ", sum(res$padj<0.05,na.rm=TRUE),
   "; multi-SNP(IVW): ", sum(res$method=="IVW"), "; single(Wald): ", sum(res$method=="Wald"))
fl("top: ", paste(head(res$gene,8), collapse=", "))

# ---------- (3) FinnGen replication attempt (osteoporosis/fracture endpoints) ----------
fg_urls <- c("https://storage.googleapis.com/finngen-public-data-r12/summary_stats/finngen_R12_M13_OSTEOPOROSIS.gz",
             "https://storage.googleapis.com/finngen-public-data-r11/summary_stats/finngen_R11_M13_OSTEOPOROSIS.gz",
             "https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_M13_OSTEOPOROSIS.gz")
fgf <- file.path(G,"finngen_osteoporosis.gz"); got <- FALSE
for(u in fg_urls){ ok <- tryCatch({ download.file(u, fgf, quiet=TRUE, mode="wb"); file.info(fgf)$size > 1e6 }, error=function(e) FALSE)
  if(isTRUE(ok)){ fl("FinnGen got: ", u); got<-TRUE; break } }
if(got){
  fg <- fread(fgf); cn <- colnames(fg)
  rs <- intersect(c("rsids","rsid","SNP"),cn)[1]; bc<-intersect(c("beta","BETA"),cn)[1]; sc<-intersect(c("sebeta","se","SE"),cn)[1]
  ec <- intersect(c("alt","effect_allele","ALT"),cn)[1]; oc<-intersect(c("ref","other_allele","REF"),cn)[1]
  fgo <- data.table(rsid=fg[[rs]], ea_o=toupper(fg[[ec]]), oa_o=toupper(fg[[oc]]), beta_o=fg[[bc]], se_o=fg[[sc]], eaf_o=NA, p_o=NA, n_o=NA)
  setkey(fgo, rsid)
  top <- head(res$gene[res$padj<0.05], 40)
  repl <- do.call(rbind, lapply(top, function(g) tryCatch(run_gene(g, outcome=fgo, label="FinnGen_osteoporosis"), error=function(e) NULL)))
  if(!is.null(repl)){ repl <- repl[!is.na(repl$b),]; fwrite(repl, "07_mr_coloc/mr_finngen_replication.tsv", sep="\t")
    fl("FinnGen replication: ", nrow(repl), " of ", length(top), " top genes tested") }
} else fl("FinnGen download failed for all endpoints; replication not run (URL/endpoint issue)")
fl("DONE harden.")
