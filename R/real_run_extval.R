# R/real_run_extval.R — TRUE external validation of the regeneration-competent
# signature on GSE316449 (PDLF, HOMER RNA-seq; NOT used in RRA training).
# Contrast: day0 (undifferentiated, RH1-2) vs day6 osteogenic (RH3-8).
.libPaths(c("F:/Rlib", .libPaths()))
suppressMessages({ library(data.table); library(pROC); library(ggplot2) })
f <- "00_data/success_dental/GSE316449/GSE316449_RH_diffOutput_condenseGenes.txt.gz"
d <- fread(f)
expr <- as.matrix(d[, 9:16]); colnames(expr) <- paste0("RH", 1:8)
sym <- sub("\\|.*$", "", d[[8]]); sym <- sub(" .*$", "", sym)   # HOMER Annotation/Divergence -> symbol
keep <- !is.na(sym) & sym != "" & sym != "NA"
expr <- expr[keep, ]; sym <- sym[keep]
expr <- log2(expr + 1)
# collapse to symbol by max mean
o <- order(rowMeans(expr), decreasing=TRUE); expr <- expr[o,]; sym <- sym[o]
expr <- expr[!duplicated(sym), ]; rownames(expr) <- sym[!duplicated(sym)]
cat("GSE316449 expr:", nrow(expr), "genes x", ncol(expr), "samples\n")
cat("symbol sanity (should be gene names):", paste(head(rownames(expr),8), collapse=", "), "\n")

group <- c("control","control","osteo","osteo","osteo","osteo","osteo","osteo")  # day0 vs day6
sig <- read.delim("03_rra_meta/regeneration_competent_signature.tsv")
up <- intersect(sig$gene[sig$direction=="up"], rownames(expr))
dn <- intersect(sig$gene[sig$direction=="down"], rownames(expr))
cat(sprintf("signature genes mapped: %d/%d up, %d/%d down\n",
            length(up), sum(sig$direction=="up"), length(dn), sum(sig$direction=="down")))
z <- t(scale(t(expr))); z[is.na(z)] <- 0
score <- colMeans(z[up,,drop=FALSE]) - colMeans(z[dn,,drop=FALSE])
res <- data.frame(sample=colnames(expr), group=group, osteo_score=round(as.numeric(score),3))
print(res)
auc <- as.numeric(pROC::auc(response=group, predictor=score, levels=c("control","osteo"), direction="<"))
p_w <- wilcox.test(score~group)$p.value
cat(sprintf("\nEXTERNAL validation GSE316449: AUC=%.3f, wilcoxon p=%.3f (day0 vs day6)\n", auc, p_w))
dir.create("05_projection", showWarnings=FALSE)
write.table(res, "05_projection/external_validation_GSE316449.tsv", sep="\t", quote=FALSE, row.names=FALSE)

p <- ggplot(res, aes(group, osteo_score, fill=group)) +
  geom_boxplot(alpha=.6, outlier.shape=NA) + geom_jitter(width=.12, size=2.4) +
  scale_fill_manual(values=c(control="#4575b4", osteo="#d73027")) +
  labs(title="External validation: regeneration-competent signature on GSE316449 (PDLF, not in training)",
       subtitle=sprintf("day0 (undiff) vs day6 osteogenic; AUC=%.2f", auc),
       y="osteogenic score", x=NULL) + theme_bw(base_size=11) + theme(legend.position="none")
ggsave("figures/Fig3b_external_validation_GSE316449.pdf", p, width=4.2, height=3.8)
ggsave("figures/Fig3b_external_validation_GSE316449.png", p, width=4.2, height=3.8, dpi=130)
cat("DONE.\n")
