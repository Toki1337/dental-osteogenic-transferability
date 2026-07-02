# R/real_run_step4_success.R — REAL Step4 success-side validation:
# project the regeneration-competent program onto GSE104473 RNA-seq (sorted
# skeletal stem/progenitor cells from regenerating MANDIBLE; SSC->BCSP->OP
# osteogenic lineage; Distraction vs Fracture; POD5/10/15). Tests whether the
# dental-derived program marks in-vivo jaw osteogenic commitment.
.libPaths(c("F:/Rlib", .libPaths()))
suppressMessages({ library(GEOquery); library(data.table); library(babelgene); library(ggplot2) })
options(timeout = 2400)
acc <- "GSE104473"; d <- file.path("00_data/success_jaw", acc)
g <- getGEO(acc, destdir="00_data/meta", GSEMatrix=TRUE, getGPL=FALSE); es <- if(is.list(g)) g[[1]] else g
pd <- Biobase::pData(es)
# per-sample featureCounts files: GeneID(mouse symbol) Length counts rpk fpkm tpm perc
ff <- list.files(d, pattern="GSM30845.*counts.*\\.txt\\.gz$", full.names=TRUE)
cat("RNA-seq sample files:", length(ff), "\n")
read1 <- function(f){ dt <- fread(f); setNames(dt$tpm, dt$GeneID) }
mats <- lapply(ff, read1)
genes <- Reduce(intersect, lapply(mats, names))
expr <- sapply(mats, function(v) v[genes]); rownames(expr) <- genes
# column names from filenames
colnames(expr) <- sub("^GSM\\d+_", "", sub("\\.counts.*$", "", basename(ff)))
expr <- log2(expr + 1)
ms <- rownames(expr)
# mouse symbol -> human ortholog. Prefer babelgene 1:1; fall back to uppercase
# (valid for the ~85% of orthologs that share the symbol) when babelgene maps too few.
hs <- tryCatch({
  ort <- babelgene::orthologs(genes = ms, species = "mouse", human = TRUE)
  m <- setNames(ort$human_symbol, ort$symbol)[ms]; m
}, error = function(e) rep(NA_character_, length(ms)))
if (mean(!is.na(hs)) < 0.4) { cat("babelgene mapped", sum(!is.na(hs)), "-> using uppercase ortholog fallback\n"); hs <- toupper(ms) }
keep <- !is.na(hs) & hs != ""
expr <- expr[keep,,drop=FALSE]; hs <- hs[keep]
o <- order(rowMeans(expr), decreasing=TRUE); expr <- expr[o,]; hs <- hs[o]
expr <- expr[!duplicated(hs),,drop=FALSE]; rownames(expr) <- hs[!duplicated(hs)]
cat("human-mapped expr:", nrow(expr), "x", ncol(expr), "\n")

# signature score
sig <- read.delim("03_rra_meta/regeneration_competent_signature.tsv")
up <- intersect(sig$gene[sig$direction=="up"], rownames(expr)); dn <- intersect(sig$gene[sig$direction=="down"], rownames(expr))
cat("signature mapped:", length(up), "up,", length(dn), "down\n")
z <- t(scale(t(expr))); z[is.na(z)] <- 0
score <- colMeans(z[up,,drop=FALSE]) - colMeans(z[dn,,drop=FALSE])

# annotate samples from column names (RNAseq_POD10_DO_BCSP1 etc.)
sn <- colnames(expr)
celltype <- ifelse(grepl("OP", sn), "OP", ifelse(grepl("BCSP", sn), "BCSP", ifelse(grepl("SSC", sn), "SSC", NA)))
cond <- ifelse(grepl("DO", sn), "Distraction", ifelse(grepl("Fx|Fracture", sn), "Fracture", NA))
pod <- as.integer(sub(".*POD([0-9]+).*", "\\1", sn))
res <- data.frame(sample=sn, celltype=factor(celltype, levels=c("SSC","BCSP","OP")), condition=cond, POD=pod, score=round(score,3))
res <- res[!is.na(res$celltype), ]
dir.create("05_projection", showWarnings=FALSE)
write.table(res, "05_projection/success_GSE104473_scores.tsv", sep="\t", quote=FALSE, row.names=FALSE)
print(res[order(res$celltype), ])

# test: score increases along osteogenic lineage SSC<BCSP<OP (Jonckheere-style: spearman on ordinal)
res$lin <- as.integer(res$celltype)
rho <- cor(res$lin, res$score, method="spearman"); p <- cor.test(res$lin, res$score, method="spearman")$p.value
cat(sprintf("\nLineage trend (SSC<BCSP<OP) vs reg-competent score: spearman rho=%.2f, p=%.2e\n", rho, p))
by_ct <- aggregate(score~celltype, res, mean); print(by_ct)

p1 <- ggplot(res, aes(celltype, score, fill=celltype)) + geom_boxplot(alpha=.6, outlier.shape=NA) +
  geom_jitter(aes(shape=condition), width=.15, size=2) +
  labs(title="Regeneration-competent program rises along the in-vivo jaw osteogenic lineage",
       subtitle=sprintf("GSE104473 mandibular regeneration (sorted SSC->BCSP->OP); spearman rho=%.2f, p=%.1e", rho, p),
       y="reg-competent score", x="sorted skeletal cell type") + theme_bw(base_size=11)
ggsave("figures/Fig3c_success_GSE104473_lineage.pdf", p1, width=6, height=4)
ggsave("figures/Fig3c_success_GSE104473_lineage.png", p1, width=6, height=4, dpi=130)
cat("DONE.\n")
