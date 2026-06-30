# R/tier1_project_invivo.R
# ============================================================================
# Tier 1-B: project EVERY signature in the panel onto the in-vivo single-cell
# atlases and ask, with a calibrated permutation null, whether each signature
# localises to the in-vivo osteoprogenitor/mesenchymal compartment MORE than a
# size-matched random gene set would.
#
# Efficiency / exactness trick
# ----------------------------
# The per-compartment MEAN of a per-cell z-scored signature score equals the
# signature-averaged per-compartment mean of each gene's z-scored expression.
# So we precompute a small genes x compartments matrix
#     M[g,c] = (mean_expr(g in compartment c) - mean_expr(g)) / sd_expr(g)
# from the sparse matrix (no densification). Then for any gene set,
#     score_c = mean_{up} M[,c]  -  mean_{down} M[,c]
# and a 2000-fold size-matched permutation null is just random row subsets.
#
# Localisation statistic L = z of the osteoprogenitor compartment's score among
# all compartments. Transferability z = (L_obs - mean(L_null)) / sd(L_null);
# one-sided empirical p tests "localises to osteoprogenitor MORE than random".
# A genuine in-vivo osteoprogenitor program (positive control RUNX2/SP7/COL1A1)
# should sit far in the right tail; a non-transferring in-vitro program should not.
#
# Outputs (10_transferability/):
#   invivo_localization.tsv        - one row per (atlas, signature): L, z, p, CI, rank
#   invivo_compartment_scores.tsv  - per-compartment standardized score per signature
#   invivo_method_check.tsv        - M-based vs stored AddModuleScore concordance
# ============================================================================
.libPaths(c("F:/Rlib", .libPaths()))
suppressWarnings(suppressMessages({
  library(Seurat); library(Matrix); library(data.table)
}))
setDTthreads(1)
source("R/utils.R")
source("config/params.R")
OUT <- "10_transferability"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(PARAMS$seed)
B <- 2000L  # permutations

sigs <- readRDS(file.path(OUT, "signatures.rds"))
# positive control: canonical osteoblast / osteoprogenitor identity markers
POS <- c("RUNX2","SP7","COL1A1","COL1A2","BGLAP","ALPL","IBSP","SPP1","SPARC","DLX5")
sigs[["POSctrl_osteo_markers"]] <- list(name="POSctrl_osteo_markers",
  source_type="positive_control", organism="human", up=POS, down=character(0),
  citation="canonical osteoblast markers", notes="should localise to osteoprogenitor")

## ---- helpers ---------------------------------------------------------------
# genes x compartments standardized-mean matrix from a sparse data matrix
build_M <- function(dat, comp) {
  comp <- as.character(comp)
  cl <- sort(unique(comp))
  ncell <- ncol(dat)
  mu  <- Matrix::rowMeans(dat)
  m2  <- Matrix::rowMeans(dat * dat)           # E[x^2]
  sdg <- sqrt(pmax(m2 - mu^2, 0))
  cmean <- sapply(cl, function(c) Matrix::rowMeans(dat[, comp == c, drop = FALSE]))
  M <- (cmean - mu) / sdg
  M <- M[is.finite(rowSums(M)) & sdg > 0, , drop = FALSE]
  list(M = M, compartments = cl)
}
# per-compartment raw signature score, then standardized across compartments
score_comps <- function(M, up, dn) {
  up <- intersect(up, rownames(M)); dn <- intersect(dn, rownames(M))
  if (length(up) < 3) return(NULL)
  s <- colMeans(M[up, , drop = FALSE])
  if (length(dn) >= 3) s <- s - colMeans(M[dn, , drop = FALSE])
  s
}
loc_stat <- function(s, osteo) {                # z of osteo compartment among compartments
  z <- (s - mean(s)) / sd(s)
  list(L = unname(z[osteo]), rank = unname(rank(-s)[osteo]), z = z, raw = s)
}

analyze_atlas <- function(tag, rds, osteo_comp, comp_col, organism) {
  .log("=== atlas ", tag, " (", organism, ") ===")
  obj <- readRDS(rds)
  dat <- GetAssayData(obj, assay = "RNA", layer = "data")        # log-normalized, sparse
  comp <- obj@meta.data[[comp_col]]
  # mouse -> human symbol on the matrix rows (strict 1:1; quantitative-grade)
  if (organism == "mouse") {
    hs <- to_human_symbols(rownames(dat), from = "mouse", strict = TRUE)
    keep <- !is.na(hs) & !duplicated(hs)
    dat <- dat[keep, , drop = FALSE]; rownames(dat) <- hs[keep]
  }
  stopifnot(osteo_comp %in% comp)
  MM <- build_M(dat, comp); M <- MM$M
  universe <- rownames(M)
  .log(tag, ": ", ncol(dat), " cells, ", length(universe), " usable genes, compartments: ",
       paste(MM$compartments, collapse = "/"))

  res <- list(); comp_tab <- list()
  for (nm in names(sigs)) {
    s <- sigs[[nm]]
    sc <- score_comps(M, s$up, s$down)
    if (is.null(sc)) next
    ls <- loc_stat(sc, osteo_comp)
    nu <- length(intersect(s$up, universe)); nd <- length(intersect(s$down, universe))
    # size-matched permutation null on L
    Lnull <- numeric(B)
    for (b in seq_len(B)) {
      up_b <- sample(universe, nu)
      sc_b <- colMeans(M[up_b, , drop = FALSE])
      if (nd >= 3) sc_b <- sc_b - colMeans(M[sample(universe, nd), , drop = FALSE])
      Lnull[b] <- (sc_b[osteo_comp] - mean(sc_b)) / sd(sc_b)
    }
    zt <- (ls$L - mean(Lnull)) / sd(Lnull)
    p_hi <- (1 + sum(Lnull >= ls$L)) / (B + 1)            # localises MORE than random
    p_lo <- (1 + sum(Lnull <= ls$L)) / (B + 1)            # localises LESS than random
    res[[nm]] <- data.table(atlas = tag, organism = organism, signature = nm,
      source_type = s$source_type, n_up_used = nu, n_down_used = nd,
      osteo_compartment = osteo_comp, osteo_score_z = round(ls$L, 3),
      osteo_rank = ls$rank, n_compartments = length(sc),
      top_compartment = MM$compartments[which.max(sc)],
      transferability_z = round(zt, 3),
      p_localizes_more = signif(p_hi, 3), p_localizes_less = signif(p_lo, 3),
      null_mean = round(mean(Lnull), 3), null_sd = round(sd(Lnull), 3))
    comp_tab[[nm]] <- data.table(atlas = tag, signature = nm,
      compartment = names(sc), score_z = round(as.numeric(scale(sc)), 3))
  }
  # method check: M-based RRA744 per-compartment vs stored AddModuleScore (prog_net), if present
  mc <- NULL
  if ("prog_net" %in% colnames(obj@meta.data) && "RRA744_pan_dental" %in% names(res)) {
    ams <- tapply(obj@meta.data$prog_net, comp, mean)
    mb  <- score_comps(M, sigs$RRA744_pan_dental$up, sigs$RRA744_pan_dental$down)
    common <- intersect(names(ams), names(mb))
    mc <- data.table(atlas = tag, metric = "RRA744 per-compartment: M-based vs AddModuleScore",
      spearman_rho = round(cor(ams[common], mb[common], method = "spearman"), 3))
  }
  rm(obj, dat, M); gc()
  list(res = rbindlist(res), comp = rbindlist(comp_tab), mc = mc)
}

atlases <- list(
  list("GSE303003_human_MRONJ", "06_scrna/GSE303003_seurat.rds", "Mesenchymal_osteo", "compartment", "human"),
  list("GSE295106_mouse_BRONJ", "06_scrna/GSE295106_seurat.rds", "Osteo_MSC",          "cell_type",    "mouse"),
  list("GSE269255_mouse_ORNJ",  "06_scrna/GSE269255_seurat.rds", "Mesenchymal_osteo", "compartment",  "mouse"))

all_res <- list(); all_comp <- list(); all_mc <- list()
for (a in atlases) {
  out <- analyze_atlas(a[[1]], a[[2]], a[[3]], a[[4]], a[[5]])
  all_res[[a[[1]]]] <- out$res; all_comp[[a[[1]]]] <- out$comp; all_mc[[a[[1]]]] <- out$mc
}
res <- rbindlist(all_res)
save_tsv(as.data.frame(res), file.path(OUT, "invivo_localization.tsv"))
save_tsv(as.data.frame(rbindlist(all_comp)), file.path(OUT, "invivo_compartment_scores.tsv"))
save_tsv(as.data.frame(rbindlist(all_mc[!sapply(all_mc, is.null)])), file.path(OUT, "invivo_method_check.tsv"))

cat("\n================ in-vivo osteoprogenitor localisation ================\n")
print(res[, .(atlas, signature, source_type, osteo_score_z, osteo_rank, top_compartment,
              transferability_z, p_localizes_more)])
