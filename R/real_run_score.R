# R/real_run_score.R — REAL osteogenic-score (regeneration-competent signature) +
# within-dataset separability (osteo vs control) and AUC. Score = mean z(up) - mean z(down)
# per sample (a standard signature-score; GSVA ssGSEA is the pipeline default when its
# full dependency chain is available — KEGGREST/GSVA binary was unavailable here).
.libPaths(c("F:/Rlib", .libPaths()))
suppressMessages({ library(pROC); library(ggplot2) })
ex <- readRDS("03_rra_meta/exprs_cache.rds")
ph <- readRDS("03_rra_meta/pheno_cache.rds")
sig <- read.delim("03_rra_meta/regeneration_competent_signature.tsv")
up <- sig$gene[sig$direction=="up"]; dn <- sig$gene[sig$direction=="down"]

zscore <- function(m) { m <- t(scale(t(m))); m[is.na(m)] <- 0; m }
rows <- list()
for (id in names(ex)) {
  m <- ex[[id]]; z <- zscore(m)
  su <- intersect(up, rownames(z)); sd <- intersect(dn, rownames(z))
  score <- colMeans(z[su,,drop=FALSE]) - colMeans(z[sd,,drop=FALSE])
  g <- ph[[id]]$group[match(colnames(m), ph[[id]]$sample)]
  rows[[id]] <- data.frame(dataset=id, sample=colnames(m), group=g, osteo_score=as.numeric(score))
}
sc <- do.call(rbind, rows)
write.table(sc, "05_projection/osteogenic_score_insample.tsv", sep="\t", quote=FALSE, row.names=FALSE)

# per-dataset AUC + wilcoxon
auc_tab <- do.call(rbind, lapply(split(sc, sc$dataset), function(d){
  if(length(unique(d$group))<2) return(NULL)
  a <- tryCatch(as.numeric(pROC::auc(response=d$group, predictor=d$osteo_score, levels=c("control","osteo"), direction="<")), error=function(e) NA)
  p <- tryCatch(wilcox.test(osteo_score~group, d)$p.value, error=function(e) NA)
  data.frame(dataset=d$dataset[1], n=nrow(d), AUC=round(a,3), wilcox_p=signif(p,3))
}))
print(auc_tab)
write.table(auc_tab, "05_projection/osteogenic_score_AUC.tsv", sep="\t", quote=FALSE, row.names=FALSE)
overall_auc <- as.numeric(pROC::auc(response=sc$group, predictor=sc$osteo_score, levels=c("control","osteo"), direction="<"))
cat(sprintf("\nPooled AUC (osteo vs control, all %d samples): %.3f\n", nrow(sc), overall_auc))

# Fig3: boxplot of osteogenic score by group, faceted by dataset
dir.create("figures", showWarnings=FALSE)
p <- ggplot(sc, aes(group, osteo_score, fill=group)) +
  geom_boxplot(outlier.shape=NA, alpha=.6) + geom_jitter(width=.15, size=1.6) +
  facet_wrap(~dataset, scales="free_y", nrow=1) +
  scale_fill_manual(values=c(control="#4575b4", osteo="#d73027")) +
  labs(title="In-vitro osteogenic signature score separates osteo-induced vs control",
       subtitle=sprintf("z-score signature (350 up / 394 down LOO-stable genes); pooled AUC=%.3f", overall_auc),
       y="osteogenic score (mean z up - mean z down)", x=NULL) +
  theme_bw(base_size=11) + theme(legend.position="none")
ggsave("figures/Fig3_osteogenic_score_separation.pdf", p, width=11, height=3.6)
ggsave("figures/Fig3_osteogenic_score_separation.png", p, width=11, height=3.6, dpi=130)
cat("DONE. Fig3 + AUC table written.\n")
