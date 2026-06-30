# R/tier1_project_lineage.R
# ============================================================================
# Tier 1-B (bulk arm): project every panel signature onto the in-vivo sorted
# skeletal-progenitor lineage GSE104473 (mandibular distraction/fracture;
# SSC -> BCSP). Readout = Spearman rho of per-sample signature score along the
# osteogenic lineage rank, with a size-matched permutation null (transferability
# z + empirical p that the signature RISES along the lineage more than random).
#
# Also re-derives the RRA-744 lineage rho under CORRECT babelgene mouse->human
# mapping (utils::to_human_symbols, now human=FALSE) to confirm the published
# negative result (rho = -0.14) is robust to the ortholog-mapping fix.
#
# Output: 10_transferability/invivo_lineage_GSE104473.tsv
# ============================================================================
.libPaths(c("F:/Rlib", .libPaths()))
suppressWarnings(suppressMessages({ library(data.table) }))
setDTthreads(1)
source("R/utils.R")
source("config/params.R")
OUT <- "10_transferability"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(PARAMS$seed); B <- 2000L

d <- "00_data/success_jaw/GSE104473"
ff <- list.files(d, pattern = "GSM30845.*counts.*\\.txt\\.gz$", full.names = TRUE)
read1 <- function(f){ dt <- fread(f); setNames(dt$tpm, dt$GeneID) }
mats <- lapply(ff, read1)
genes <- Reduce(intersect, lapply(mats, names))
expr <- sapply(mats, function(v) v[genes]); rownames(expr) <- genes
colnames(expr) <- sub("^GSM\\d+_", "", sub("\\.counts.*$", "", basename(ff)))
expr <- log2(expr + 1)

# mouse -> human (CORRECT 1:1 babelgene mapping via fixed utils::to_human_symbols)
hs <- to_human_symbols(rownames(expr), from = "mouse", strict = TRUE)
keep <- !is.na(hs); expr <- expr[keep, , drop = FALSE]; hs <- hs[keep]
o <- order(rowMeans(expr), decreasing = TRUE); expr <- expr[o, ]; hs <- hs[o]
expr <- expr[!duplicated(hs), , drop = FALSE]; rownames(expr) <- hs[!duplicated(hs)]
.log("GSE104473 human-mapped expr: ", nrow(expr), " x ", ncol(expr), " (correct babelgene 1:1)")

sn <- colnames(expr)
celltype <- ifelse(grepl("BCSP", sn), "BCSP", ifelse(grepl("SSC", sn), "SSC", NA))
lin <- ifelse(celltype == "SSC", 1L, ifelse(celltype == "BCSP", 2L, NA))
keepS <- !is.na(lin); exprS <- expr[, keepS, drop = FALSE]; lin <- lin[keepS]
.log("usable lineage samples: ", length(lin), " (SSC=", sum(lin==1), ", BCSP=", sum(lin==2), ")")

z <- t(scale(t(exprS))); z[is.na(z)] <- 0
universe <- rownames(z)
sig_score <- function(up, dn){
  up <- intersect(up, universe); dn <- intersect(dn, universe)
  if (length(up) < 3) return(NULL)
  s <- colMeans(z[up, , drop = FALSE]); if (length(dn) >= 3) s <- s - colMeans(z[dn, , drop = FALSE]); s
}

sigs <- readRDS(file.path(OUT, "signatures.rds"))
sigs[["POSctrl_osteo_markers"]] <- list(name="POSctrl_osteo_markers", source_type="positive_control",
  organism="human", up=c("RUNX2","SP7","COL1A1","COL1A2","BGLAP","ALPL","IBSP","SPP1","SPARC","DLX5"),
  down=character(0), citation="canonical osteoblast markers", notes="")

res <- list()
for (nm in names(sigs)) {
  s <- sigs[[nm]]; sc <- sig_score(s$up, s$down); if (is.null(sc)) next
  rho <- cor(lin, sc, method = "spearman")
  nu <- length(intersect(s$up, universe)); nd <- length(intersect(s$down, universe))
  rn <- numeric(B)
  for (b in seq_len(B)) {
    sb <- colMeans(z[sample(universe, nu), , drop = FALSE])
    if (nd >= 3) sb <- sb - colMeans(z[sample(universe, nd), , drop = FALSE])
    rn[b] <- cor(lin, sb, method = "spearman")
  }
  res[[nm]] <- data.table(atlas = "GSE104473_mouse_lineage", signature = nm, source_type = s$source_type,
    n_up_used = nu, n_down_used = nd, lineage_rho = round(rho, 3),
    transferability_z = round((rho - mean(rn)) / sd(rn), 3),
    p_rises_more = signif((1 + sum(rn >= rho)) / (B + 1), 3),
    null_mean = round(mean(rn), 3), null_sd = round(sd(rn), 3))
}
res <- rbindlist(res)
save_tsv(as.data.frame(res), file.path(OUT, "invivo_lineage_GSE104473.tsv"))
cat("\n========== GSE104473 lineage transfer (SSC->BCSP) ==========\n")
print(res[, .(signature, source_type, lineage_rho, transferability_z, p_rises_more)])
cat(sprintf("\nRRA-744 lineage rho under correct mapping = %.3f (published toupper-based value: -0.14)\n",
            res[signature=="RRA744_pan_dental", lineage_rho]))
