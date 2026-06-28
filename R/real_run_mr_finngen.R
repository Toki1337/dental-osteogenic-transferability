# R/real_run_mr_finngen.R
# Context-matched / independent-cohort MR replication (referee R4 action B):
# re-run the corrected lead-cis-SNP Wald MR with eQTLGen exposures against an
# INDEPENDENT outcome cohort and a CLINICAL endpoint: FinnGen osteoporosis
# (M13_OSTEOPOROSIS), instead of UK-Biobank heel-eBMD. Tests whether the bone
# genes recovered on eBMD replicate on a clinical osteoporosis diagnosis in a
# separate population. Same lead-SNP design (no LD reference, palindromes dropped).
.libPaths(c("F:/Rlib", .libPaths()))
suppressMessages({ library(data.table) })
setDTthreads(0)
G <- "00_data/genetic"
fl <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(),"%H:%M:%S"), paste0(...)))

sig <- fread("03_rra_meta/regeneration_competent_signature.tsv"); cand <- unique(sig$gene)
fl("reading eQTLGen ..."); eq <- fread(file.path(G,"eqtlgen_cis_sig.txt.gz"))
eq <- eq[GeneSymbol %in% cand & Pvalue < 5e-8]

fl("reading FinnGen osteoporosis ...")
fg <- fread(file.path(G,"finngen_osteoporosis.gz"),
            select=c("rsids","ref","alt","beta","sebeta","af_alt","pval"))
# FinnGen: beta is for ALT allele; effect=alt, other=ref, eaf=af_alt
setnames(fg, c("rsid","oa_o","ea_o","beta_o","se_o","eaf_o","p_o"))
fg <- fg[!is.na(rsid) & rsid!="" & !is.na(beta_o) & !is.na(eaf_o)]
# rsids field may contain multiple comma-separated ids; keep first
fg[, rsid := tstrsplit(rsid, ",", keep=1L)]
setkey(fg, rsid)
fl("FinnGen variants: ", nrow(fg))

zhu  <- function(z,p,N){ d<-2*p*(1-p)*(N+z^2); ifelse(d>0, z/sqrt(d), NA) }
zhse <- function(p,N,z){ d<-2*p*(1-p)*(N+z^2); ifelse(d>0, 1/sqrt(d), NA) }
is_palin <- function(a1,a2) paste0(toupper(a1),toupper(a2)) %in% c("AT","TA","CG","GC")

run_gene <- function(g){
  e <- eq[GeneSymbol==g]; m <- merge(e, fg, by.x="SNP", by.y="rsid")
  if(!nrow(m)) return(NULL)
  ok <- (toupper(m$AssessedAllele)==toupper(m$ea_o)&toupper(m$OtherAllele)==toupper(m$oa_o)) |
        (toupper(m$AssessedAllele)==toupper(m$oa_o)&toupper(m$OtherAllele)==toupper(m$ea_o))
  palin <- is_palin(m$AssessedAllele, m$OtherAllele)
  m <- m[ok & !palin]; if(!nrow(m)) return(NULL)
  flip <- toupper(m$AssessedAllele)!=toupper(m$ea_o); eaf <- ifelse(flip,1-m$eaf_o,m$eaf_o)
  d <- data.frame(rsid=m$SNP, pval=m$Pvalue, bx=zhu(m$Zscore,eaf,m$NrSamples),
                  bxse=zhse(eaf,m$NrSamples,m$Zscore), by=ifelse(flip,-m$beta_o,m$beta_o), byse=m$se_o)
  d <- d[is.finite(d$bx)&is.finite(d$bxse)&is.finite(d$by)&is.finite(d$byse)&d$bxse>0&d$byse>0,]
  if(!nrow(d)) return(NULL)
  d <- d[order(d$pval),][1,,drop=FALSE]            # lead cis-SNP
  b <- d$by/d$bx; se <- abs(d$byse/d$bx)
  data.frame(gene=g, lead_snp=d$rsid, b_osteoporosis=b, se=se, pval=2*pnorm(-abs(b/se)))
}
genes <- unique(eq$GeneSymbol)
res <- do.call(rbind, lapply(genes, function(g) tryCatch(run_gene(g), error=function(e) NULL)))
res <- res[!is.na(res$b_osteoporosis),]
res$padj <- p.adjust(res$pval, "BH")
res$direction_signature <- sig$direction[match(res$gene, sig$gene)]
res <- res[order(res$pval),]
fwrite(res, "07_mr_coloc/mr_finngen_osteoporosis.tsv", sep="\t")
fl("FinnGen-osteoporosis lead-SNP MR: ", nrow(res), " genes; FDR<0.05: ", sum(res$padj<0.05,na.rm=TRUE))

# ---- concordance with eBMD for positive controls + druggable core ----
eb <- fread("07_mr_coloc/mr_results_leadSNP.tsv")  # eBMD lead-SNP
mg <- merge(res[,c("gene","b_osteoporosis","pval","padj")],
            eb[,c("gene","b","padj")], by="gene", suffixes=c("_op",""))
setnames(mg, c("b","padj"), c("b_eBMD","padj_eBMD"))
# osteoporosis is a DISEASE (higher beta = higher risk); eBMD higher = more bone.
# Concordant biology = OPPOSITE signs (more bone-expression that RAISES eBMD should LOWER osteoporosis risk).
mg$concordant_opposite_sign <- sign(mg$b_eBMD) != sign(mg$b_osteoporosis)
known <- c("WNT16","KREMEN1","MGP","SMAD3","HSPG2")
fl("=== positive controls: eBMD vs FinnGen osteoporosis ===")
print(mg[mg$gene %in% known, ])
fl("=== druggable double-evidence candidates ===")
print(mg[mg$gene %in% c("GREM2","SRGN"), ])
# overall concordance among eBMD-FDR-significant genes
sigset <- mg[!is.na(mg$padj_eBMD) & mg$padj_eBMD<0.05, ]
fl("among eBMD-FDR-sig genes also tested in FinnGen (n=", nrow(sigset),
   "): opposite-sign concordance = ", round(100*mean(sigset$concordant_opposite_sign),1), "%; ",
   "FinnGen-FDR-sig too: ", sum(sigset$padj<0.05, na.rm=TRUE))
fwrite(mg, "07_mr_coloc/mr_eBMD_vs_finngen_concordance.tsv", sep="\t")
fl("DONE FinnGen replication.")
