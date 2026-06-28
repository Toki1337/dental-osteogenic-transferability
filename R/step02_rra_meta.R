# R/step02_rra_meta.R
# Step 2 — Pan-dental-source RRA meta-signature (L1 "regeneration-competent program").
# Robust Rank Aggregation across per-dataset ranked lists -> direction-consistent
# osteogenic signature, with leave-one-dataset-out (LOO) stability.
# Output:
#   03_rra_meta/rra_up.tsv, rra_down.tsv         (meta-DEG with RRA Score = corrected p)
#   03_rra_meta/loo_stability.tsv                (per-gene fraction of LOO runs retained)
#   03_rra_meta/regeneration_competent_signature.tsv (final stable, direction-consistent)
# Innovation guard: signature must be ROBUST to dropping any single dataset (LOO),
# directly answering the "cross-platform drift" reviewer concern.

source("R/utils.R"); init_pipeline()
suppressMessages(library(RobustRankAggreg))

ranked_files <- list.files(PARAMS$dir_de, pattern = "_ranked\\.tsv$", full.names = TRUE)
stopifnot(length(ranked_files) >= PARAMS$rra_min_datasets)
lists <- lapply(ranked_files, function(f) {
  d <- utils::read.delim(f, stringsAsFactors = FALSE)
  list(acc = sub("_ranked\\.tsv$", "", basename(f)), d = d)
})
names(lists) <- vapply(lists, function(x) x$acc, character(1))
.log("RRA over ", length(lists), " datasets: ", paste(names(lists), collapse = ", "))

# Build top-fraction up / down gene rankings per dataset (by stat).
make_glists <- function(lst, direction = c("up", "down")) {
  direction <- match.arg(direction)
  lapply(lst, function(x) {
    d <- x$d[!is.na(x$d$stat), ]
    d <- d[order(if (direction == "up") -d$stat else d$stat), ]
    n <- max(50, ceiling(PARAMS$rra_topN * nrow(d)))
    head(d$gene, n)
  })
}

# Background gene universe = union of ALL measured genes across the per-dataset
# ranked.tsv (the untruncated tested universe), and per-gene coverage = number
# of datasets that actually measured the gene. Both are computed from `lst` so
# they shrink correctly inside leave-one-out runs.
measured_universe <- function(lst) unique(unlist(lapply(lst, function(x) x$d$gene)))
gene_coverage <- function(lst) table(unlist(lapply(lst, function(x) unique(x$d$gene))))

run_rra <- function(glists, lst) {
  # N must be the real measured background (untruncated universe across datasets),
  # NOT the size of the truncated top-fraction union.
  N <- length(measured_universe(lst))
  agg <- aggregateRanks(glist = glists, N = N, method = PARAMS$rra_method)
  # The RRA Score returned by aggregateRanks(method="RRA") is already a
  # multiple-testing-corrected significance value (Bonferroni-on-beta-order
  # statistic; Kõlde et al. 2012, Bioinformatics). Do NOT apply a second BH on
  # top of it -- that would be double correction. Use Score directly as padj.
  agg$padj <- agg$Score
  # coverage: number of datasets in which the gene was actually MEASURED
  # (tested universe), used for eligibility per Methods.
  cov <- gene_coverage(lst)
  agg$coverage <- as.integer(cov[agg$Name])
  # top-fraction co-occurrence count (how many top lists the gene appears in),
  # kept for transparency but NOT used for eligibility.
  freq <- table(unlist(glists))
  agg$n_topfrac <- as.integer(freq[agg$Name])
  agg[order(agg$Score), ]
}

up   <- run_rra(make_glists(lists, "up"),   lists)
down <- run_rra(make_glists(lists, "down"), lists)
sig_up   <- up[up$padj < PARAMS$rra_fdr & up$coverage >= PARAMS$rra_min_datasets, ]
sig_down <- down[down$padj < PARAMS$rra_fdr & down$coverage >= PARAMS$rra_min_datasets, ]
save_tsv(up,   file.path(PARAMS$dir_rra, "rra_up.tsv"))
save_tsv(down, file.path(PARAMS$dir_rra, "rra_down.tsv"))

# ---- leave-one-dataset-out stability ----
loo_hits <- list()
for (drop in names(lists)) {
  sub <- lists[setdiff(names(lists), drop)]
  u <- run_rra(make_glists(sub, "up"),   sub); u <- u$Name[u$padj < PARAMS$rra_fdr & u$coverage >= PARAMS$rra_min_datasets]
  d <- run_rra(make_glists(sub, "down"), sub); d <- d$Name[d$padj < PARAMS$rra_fdr & d$coverage >= PARAMS$rra_min_datasets]
  loo_hits[[drop]] <- list(up = u, down = d)
}
stability <- function(genes, dir) {
  frac <- vapply(genes, function(g) mean(vapply(loo_hits, function(h) g %in% h[[dir]], logical(1))), numeric(1))
  data.frame(gene = genes, direction = dir, loo_fraction = frac, stringsAsFactors = FALSE)
}
loo_tab <- rbind(stability(sig_up$Name, "up"), stability(sig_down$Name, "down"))
save_tsv(loo_tab, file.path(PARAMS$dir_rra, "loo_stability.tsv"))

final <- loo_tab[loo_tab$loo_fraction >= PARAMS$loo_stability_min, ]
# Merge by (gene, direction): a gene could appear in both up and down across
# directions, so joining on "gene" alone would Cartesian-explode the rows.
final <- merge(final, rbind(
  data.frame(gene = sig_up$Name,   direction = "up",   rra_score = sig_up$Score,   padj = sig_up$padj,   coverage = sig_up$coverage,   n_topfrac = sig_up$n_topfrac,   stringsAsFactors = FALSE),
  data.frame(gene = sig_down$Name, direction = "down", rra_score = sig_down$Score, padj = sig_down$padj, coverage = sig_down$coverage, n_topfrac = sig_down$n_topfrac, stringsAsFactors = FALSE)
), by = c("gene", "direction"))
final <- final[order(final$direction, final$rra_score), ]
save_tsv(final, file.path(PARAMS$dir_rra, "regeneration_competent_signature.tsv"))

cat(sprintf("\nRRA meta-signature: %d up, %d down (RRA Score<%.2f, measured in >=%d datasets); LOO-stable (>=%.0f%%): %d\n",
            nrow(sig_up), nrow(sig_down), PARAMS$rra_fdr, PARAMS$rra_min_datasets,
            100 * PARAMS$loo_stability_min, nrow(final)))
dump_session()
