# env/install_packages.R — one-shot dependency installer for the pipeline.
# Run once:  Rscript env/install_packages.R
# Tested target: R >= 4.3, Bioconductor >= 3.18.

options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

cran_pkgs <- c(
  "data.table", "ggplot2", "ggrepel", "pheatmap", "RColorBrewer", "patchwork",
  "RobustRankAggreg", "WGCNA", "igraph", "ggraph",
  "babelgene", "matrixStats", "Matrix",
  "Seurat", "harmony",
  "remotes", "yaml", "jsonlite", "ggpubr"
)

bioc_pkgs <- c(
  "GEOquery", "Biobase", "limma", "edgeR", "DESeq2",
  "GSVA", "GSEABase", "clusterProfiler", "org.Hs.eg.db", "org.Mm.eg.db",
  "STRINGdb", "ComplexHeatmap", "fgsea", "AnnotationDbi",
  "TwoSampleMR", "MendelianRandomization", "coloc"
)

github_pkgs <- c(
  "jinworks/CellChat",          # cell-cell communication
  "MRCIEU/TwoSampleMR",         # in case CRAN/Bioc copy lags
  "saezlab/decoupleR"           # optional regulon/pathway activity
)

inst <- function(p, fn) {
  miss <- p[!vapply(sub("^.*/", "", p), function(x) requireNamespace(x, quietly = TRUE), logical(1))]
  if (length(miss)) { message("Installing: ", paste(miss, collapse = ", ")); fn(miss) }
}

inst(cran_pkgs, function(m) install.packages(m))
inst(bioc_pkgs, function(m) BiocManager::install(m, update = FALSE, ask = FALSE))
for (g in github_pkgs) {
  nm <- sub("^.*/", "", g)
  if (!requireNamespace(nm, quietly = TRUE)) try(remotes::install_github(g, upgrade = "never"))
}

message("\nDone. Verify with: Rscript -e 'source(\"R/utils.R\"); init_pipeline()'")
message("NOTE: clue.io (CMap/L1000) needs a free account; L1000CDS2/L1000FWD have open APIs (see R/step09_cmap_repurposing.R).")
