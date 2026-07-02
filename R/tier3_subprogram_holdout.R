# R/tier3_subprogram_holdout.R
# ============================================================================
# Non-circular check for the 27-gene "in-vivo-correspondent subprogram" (defined in
# R/tier3_dfc_deepdive.R on the three DISEASE atlases): does it also localise to the
# osteoprogenitor compartment in the HELD-OUT homeostatic atlas GSE316924, which was
# NOT used to define it? Answers the circularity objection.
#   Rscript R/tier3_subprogram_holdout.R
# Requires: 06_scrna/GSE316924_seurat.rds (from R/tier3_add_homeostatic_GSE316924.R)
#           10_transferability/rra744_correspondent_subprogram_genes.txt
# ============================================================================
.libPaths(c(Sys.getenv("R_TRANSFERAUDIT_LIB", unset = .libPaths()[1]), .libPaths()))
suppressWarnings(suppressMessages({ library(SeuratObject); library(Matrix) }))
source("R/utils.R"); source("config/params.R"); set.seed(PARAMS$seed); B <- 2000L
sub <- readLines("10_transferability/rra744_correspondent_subprogram_genes.txt")
obj <- readRDS("06_scrna/GSE316924_seurat.rds")
dat <- GetAssayData(obj, assay = "RNA", layer = "data")
hs <- to_human_symbols(rownames(dat), from = "mouse", strict = TRUE)
keep <- !is.na(hs) & !duplicated(hs); dat <- dat[keep, , drop = FALSE]; rownames(dat) <- hs[keep]
comp <- as.character(obj$compartment); cl <- sort(unique(comp))
mu <- Matrix::rowMeans(dat); sdg <- sqrt(pmax(Matrix::rowMeans(dat * dat) - mu^2, 0))
M <- (sapply(cl, function(c) Matrix::rowMeans(dat[, comp == c, drop = FALSE])) - mu) / sdg
M <- M[is.finite(rowSums(M)) & sdg > 0, , drop = FALSE]; univ <- rownames(M); osteo <- "Osteo_MSC"
up <- intersect(sub, univ); sc <- colMeans(M[up, , drop = FALSE])
L <- (sc[osteo] - mean(sc)) / sd(sc); rk <- rank(-sc)[osteo]
Ln <- numeric(B); for (b in seq_len(B)) { x <- colMeans(M[sample(univ, length(up)), , drop = FALSE]); Ln[b] <- (x[osteo] - mean(x)) / sd(x) }
cat(sprintf("27-gene correspondent subprogram HELD-OUT on GSE316924 (healthy, not used to define it): n_used=%d osteo_rank=%d z=%.2f p=%.4f localises=%s\n",
  length(up), rk, (L - mean(Ln)) / sd(Ln), (1 + sum(Ln >= L)) / (B + 1), (rk == 1 & (1 + sum(Ln >= L)) / (B + 1) < 0.05)))
