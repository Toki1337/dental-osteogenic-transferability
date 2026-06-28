# R/real_run_dedup_signature.R — referee R5 M3 fix.
# The direction-stratified RRA let 24 genes enter BOTH the up and down signatures,
# so the 768-row signature contains only 744 UNIQUE genes (double counting the
# headline "370 up / 398 down"). Resolve each gene to ONE direction by its stronger
# (smaller) RRA score = net RRA direction; report the true unique counts.
.libPaths(c("F:/Rlib", .libPaths()))
suppressMessages(library(data.table))
sig <- fread("03_rra_meta/regeneration_competent_signature.tsv")
cat("rows:", nrow(sig), "  unique genes:", uniqueN(sig$gene), "  bidirectional:", sum(table(sig$gene)>1), "\n")

# keep, per gene, the direction with the smaller rra_score (stronger evidence)
setorder(sig, gene, rra_score)
dd <- sig[, .SD[1], by=gene]                 # first row per gene = strongest direction
setorder(dd, direction, rra_score)
fwrite(dd, "03_rra_meta/regeneration_competent_signature_dedup.tsv", sep="\t")
cat(sprintf("deduplicated: %d unique genes (%d up, %d down)\n",
            nrow(dd), sum(dd$direction=="up"), sum(dd$direction=="down")))

# how were the 24 conflicted genes resolved?
conf <- sig[gene %in% sig[, .N, by=gene][N>1, gene]]
res  <- dd[gene %in% conf$gene, .(gene, kept_direction=direction, rra_score)]
fwrite(res, "03_rra_meta/dedup_conflicted_genes.tsv", sep="\t")
cat("conflicted genes resolved to:", paste(table(res$kept_direction), names(table(res$kept_direction)), collapse="  "), "\n")
cat("DONE dedup.\n")
