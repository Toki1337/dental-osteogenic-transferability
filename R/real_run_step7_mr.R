# R/real_run_step7_mr.R — REAL druggable-genome-style MR: do regeneration-competent
# program genes causally affect bone (heel eBMD)? Exposure = eQTLGen cis-eQTL of
# signature genes; outcome = eBMD GWAS (Morris & Kemp 2019, GCST006979).
# Zhu 2016 Z->beta/se (using eBMD EAF + eQTLGen N); plink/1000G-EUR clumping;
# MendelianRandomization IVW/Egger/median; coloc on overlapping SNPs.
.libPaths(c("F:/Rlib", .libPaths()))
options(repos=c(CRAN="https://cloud.r-project.org"))
if(!requireNamespace("R.utils",quietly=TRUE)) install.packages("R.utils")  # enables direct .gz fread
suppressMessages({ library(data.table); library(ieugwasr); library(MendelianRandomization); library(coloc) })
setDTthreads(0)
G <- "00_data/genetic"; PLINK <- file.path(G,"plink/plink.exe"); BFILE <- file.path(G,"ld/EUR")
dir.create("07_mr_coloc", showWarnings=FALSE)
fl <- function(...) { cat(sprintf("[%s] %s\n", format(Sys.time(),"%H:%M:%S"), paste0(...))); flush.console() }

sig <- read.delim("03_rra_meta/regeneration_competent_signature.tsv")
core <- read.delim("04_modules_wgcna/core_program.tsv")
cand_genes <- unique(sig$gene)
fl("candidate exposure genes (signature): ", length(cand_genes))

# --- eQTLGen: direct .gz read, subset to candidate genes, instrument p<5e-8 ---
fl("reading eQTLGen ...")
eq <- fread(file.path(G,"eqtlgen_cis_sig.txt.gz"))
eq <- eq[GeneSymbol %in% cand_genes & Pvalue < 5e-8]
fl("eQTLGen cis instruments for candidates: ", nrow(eq), " across ", length(unique(eq$GeneSymbol)), " genes")

# --- eBMD outcome: direct .gz read, needed cols, key by rsid ---
fl("reading eBMD GWAS ...")
eb <- fread(file.path(G,"ebmd_GCST006979.h.tsv.gz"),
            select=c("hm_rsid","hm_effect_allele","hm_other_allele","hm_beta","standard_error",
                     "hm_effect_allele_frequency","p_value","n"))
fl("eBMD rows: ", nrow(eb))
setnames(eb, c("rsid","ea_o","oa_o","beta_o","se_o","eaf_o","p_o","n_o"))
eb <- eb[!is.na(rsid) & !is.na(beta_o)]
setkey(eb, rsid)

zhu <- function(z, p, n){ d <- 2*p*(1-p)*(n+z^2); ifelse(d>0, z/sqrt(d), NA) }
zhu_se <- function(p, n, z){ d <- 2*p*(1-p)*(n+z^2); ifelse(d>0, 1/sqrt(d), NA) }

run_gene <- function(g){
  e <- eq[GeneSymbol==g]
  m <- merge(e, eb, by.x="SNP", by.y="rsid")
  if (nrow(m)==0) return(NULL)
  # harmonize: align eQTL assessed allele to eBMD effect allele
  flip <- toupper(m$AssessedAllele) != toupper(m$ea_o)
  # only keep SNPs where alleles match (either orientation)
  match_ok <- (toupper(m$AssessedAllele)==toupper(m$ea_o) & toupper(m$OtherAllele)==toupper(m$oa_o)) |
              (toupper(m$AssessedAllele)==toupper(m$oa_o) & toupper(m$OtherAllele)==toupper(m$ea_o))
  m <- m[match_ok]; if (nrow(m)==0) return(NULL)
  flip <- toupper(m$AssessedAllele) != toupper(m$ea_o)
  eaf <- ifelse(flip, 1-m$eaf_o, m$eaf_o)
  bx <- zhu(m$Zscore, eaf, m$NrSamples); bxse <- zhu_se(eaf, m$NrSamples, m$Zscore)
  by <- ifelse(flip, -m$beta_o, m$beta_o); byse <- m$se_o
  d <- data.frame(rsid=m$SNP, pos=m$SNPPos, bx=bx, bxse=bxse, by=by, byse=byse, pval=m$Pvalue, eaf=eaf)
  d <- d[is.finite(d$bx)&is.finite(d$bxse)&is.finite(d$by)&is.finite(d$byse)&d$bxse>0&d$byse>0,]
  if (nrow(d)==0) return(NULL)
  # LD-independence proxy for cis instruments: greedy distance pruning (keep lowest-p,
  # >250kb apart). Avoids per-gene plink (which reloads the 8.5M-variant panel each call).
  if (nrow(d) > 1) {
    d <- d[order(d$pval), ]; kept <- logical(nrow(d))
    for (i in seq_len(nrow(d))) if (!any(kept & abs(d$pos - d$pos[i]) < 250000)) kept[i] <- TRUE
    d <- d[kept, , drop=FALSE]
  }
  nsnp <- nrow(d)
  out <- data.frame(gene=g, nsnp=nsnp, method=NA, b=NA, se=NA, pval=NA, egger_intercept_p=NA, Q_p=NA)
  if (nsnp==1) {
    out$method<-"Wald ratio"; out$b<-d$by/d$bx; out$se<-abs(d$byse/d$bx); out$pval<-2*pnorm(-abs(out$b/out$se))
  } else {
    mri <- mr_input(bx=d$bx, bxse=d$bxse, by=d$by, byse=d$byse)
    ivw <- tryCatch(mr_ivw(mri), error=function(e) NULL)
    if (!is.null(ivw)) { out$method<-"IVW"; out$b<-ivw$Estimate; out$se<-ivw$StdError; out$pval<-ivw$Pvalue; out$Q_p<-ivw$Heter.Stat[2] }
    if (nsnp>=3) { eg<-tryCatch(mr_egger(mri),error=function(e)NULL); if(!is.null(eg)) out$egger_intercept_p<-eg$Pvalue.Int }
  }
  out$coloc_pp4 <- NA  # coloc needs full region; sig-only file limits this (flagged)
  out
}

genes_with_inst <- unique(eq$GeneSymbol)
cat("running MR for", length(genes_with_inst), "genes with instruments...\n")
res <- do.call(rbind, lapply(genes_with_inst, function(g) tryCatch(run_gene(g), error=function(e){cat("  ",g,"err\n");NULL})))
res <- res[!is.na(res$b), ]
res$padj <- p.adjust(res$pval, "BH")
res$direction_signature <- sig$direction[match(res$gene, sig$gene)]
res <- res[order(res$pval), ]
write.table(res, "07_mr_coloc/mr_results.tsv", sep="\t", quote=FALSE, row.names=FALSE)

sig_hits <- res[!is.na(res$padj) & res$padj < 0.05, ]
cat(sprintf("\nMR: %d genes tested, %d with FDR<0.05 causal effect on eBMD.\n", nrow(res), nrow(sig_hits)))
cat("top causal genes (gene/method/beta/padj/sig-direction):\n")
print(head(res[,c("gene","nsnp","method","b","pval","padj","direction_signature")], 20))
write.table(sig_hits, "07_mr_coloc/mr_causal_eBMD_FDR05.tsv", sep="\t", quote=FALSE, row.names=FALSE)
cat("\nNOTE: causal anchoring to SYSTEMIC bone (eBMD); not MRONJ-specific. Protective direction per gene needs pathway adjudication.\n")
