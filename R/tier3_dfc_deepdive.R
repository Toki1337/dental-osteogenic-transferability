# R/tier3_dfc_deepdive.R
# ============================================================================
# Item A (reviewer request): why does the dental-follicle (DFC) DE signature
# transfer to the in-vivo osteoprogenitor compartment in all three scRNA atlases
# when the RRA-744 meta-signature does not?
#
# A1 (already computed WITHOUT R, from in-hand files, see
#     10_transferability/dfc_vs_rra_contrast.txt):
#   - DFC (394 genes) shares only 193/394 (49%) with RRA-744; 201 DFC-unique.
#   - DFC osteo-identity content (5.3%) ~ RRA-744 (5.8%)  -> NOT the discriminator.
#   - DFC does NOT overlap the curated GO osteo sets more than RRA (Jaccard 0.015 vs 0.028).
#   - DFC carries LESS culture contamination (interferon 2.8% vs 5.9%; cell-cycle
#     0.25% vs 1.6%), but culture load does not predict TS across the panel.
#   => the answer must be WHICH genes DFC contains and whether they are expressed
#      in the in-vivo osteoprogenitor compartment. That needs the atlas (A2/A3 below).
#
# A2  leading-edge: per gene, its standardized osteoprogenitor-compartment expression
#     M[g, osteo] in each atlas; rank DFC vs RRA-744 up-genes -> which genes drive
#     (or fail to drive) osteoprogenitor localisation.
# A3  define an "in-vivo-correspondent subprogram" = signature genes whose mean
#     osteoprogenitor z across the 3 atlases exceeds the positive-control genes'
#     median, and re-project it with the calibrated Transferability Score to test
#     whether such a subprogram, extracted from RRA-744, localises (i.e. whether a
#     transferable core is diluted inside the full meta-signature).
#
# Reuses the exact machinery of R/tier1_project_invivo.R (build_M, permutation null).
# Run AFTER tier1 (needs 06_scrna/*_seurat.rds and 10_transferability/signatures.rds).
#   Rscript R/tier3_dfc_deepdive.R
# ============================================================================
.libPaths(c(Sys.getenv("R_TRANSFERAUDIT_LIB", unset = .libPaths()[1]), .libPaths()))
# reading the atlas matrices only needs SeuratObject (GetAssayData), not full Seurat
suppressWarnings(suppressMessages({ library(SeuratObject); library(Matrix); library(data.table) }))
setDTthreads(1)
source("R/utils.R"); source("config/params.R")
OUT <- "10_transferability"; set.seed(PARAMS$seed); B <- 2000L

sigs <- readRDS(file.path(OUT, "signatures.rds"))
POS  <- c("RUNX2","SP7","COL1A1","COL1A2","BGLAP","ALPL","IBSP","SPP1","SPARC","DLX5")

# --- same standardized genes x compartments matrix as tier1 ------------------
build_M <- function(dat, comp) {
  comp <- as.character(comp); cl <- sort(unique(comp))
  mu <- Matrix::rowMeans(dat); m2 <- Matrix::rowMeans(dat * dat)
  sdg <- sqrt(pmax(m2 - mu^2, 0))
  cmean <- sapply(cl, function(c) Matrix::rowMeans(dat[, comp == c, drop = FALSE]))
  M <- (cmean - mu) / sdg
  M[is.finite(rowSums(M)) & sdg > 0, , drop = FALSE]
}
score_comps <- function(M, up, dn) {
  up <- intersect(up, rownames(M)); dn <- intersect(dn, rownames(M))
  if (length(up) < 3) return(NULL)
  s <- colMeans(M[up, , drop = FALSE]); if (length(dn) >= 3) s <- s - colMeans(M[dn, , drop = FALSE]); s
}
transfer_z <- function(M, up, dn, osteo, univ) {              # calibrated Transferability z + p
  sc <- score_comps(M, up, dn); if (is.null(sc)) return(NULL)
  L <- (sc[osteo] - mean(sc)) / sd(sc); rk <- rank(-sc)[osteo]
  nu <- length(intersect(up, univ)); nd <- length(intersect(dn, univ))
  Ln <- numeric(B)
  for (b in seq_len(B)) { s <- colMeans(M[sample(univ, nu), , drop = FALSE])
    if (nd >= 3) s <- s - colMeans(M[sample(univ, nd), , drop = FALSE]); Ln[b] <- (s[osteo]-mean(s))/sd(s) }
  list(L = unname(L), rank = unname(rk), z = (L-mean(Ln))/sd(Ln), p = (1+sum(Ln>=L))/(B+1))
}

atlases <- list(
  c("GSE303003_human_MRONJ","06_scrna/GSE303003_seurat.rds","Mesenchymal_osteo","compartment","human"),
  c("GSE295106_mouse_BRONJ","06_scrna/GSE295106_seurat.rds","Osteo_MSC","cell_type","mouse"),
  c("GSE269255_mouse_ORNJ","06_scrna/GSE269255_seurat.rds","Mesenchymal_osteo","compartment","mouse"))

per_gene <- list()   # A2: per-gene osteoprogenitor z per atlas
Ms <- list()         # keep M per atlas for A3 re-projection
for (a in atlases) {
  tag <- a[1]; obj <- readRDS(a[2]); dat <- GetAssayData(obj, assay="RNA", layer="data")
  comp <- obj@meta.data[[a[4]]]
  if (a[5] == "mouse") { hs <- to_human_symbols(rownames(dat), from="mouse", strict=TRUE)
    keep <- !is.na(hs) & !duplicated(hs); dat <- dat[keep,,drop=FALSE]; rownames(dat) <- hs[keep] }
  M <- build_M(dat, comp); Ms[[tag]] <- list(M = M, osteo = a[3])
  per_gene[[tag]] <- data.table(atlas = tag, gene = rownames(M), osteo_z = round(M[, a[3]], 3))
  rm(obj, dat); gc()
}
pg <- rbindlist(per_gene)
# mean osteoprogenitor z across the 3 atlases, per gene
meanz <- pg[, .(mean_osteo_z = round(mean(osteo_z), 3), n_atlas = .N), by = gene]

# ---- A2: leading-edge of DFC vs RRA-744 (up genes) --------------------------
le <- function(nm) meanz[gene %in% sigs[[nm]]$up][order(-mean_osteo_z)][, signature := nm]
le_tab <- rbindlist(list(le("DE_DFC_GSE49007"), le("RRA744_pan_dental"), le("POSctrl_osteo_markers")), fill = TRUE)
save_tsv(as.data.frame(le_tab), file.path(OUT, "dfc_rra_leading_edge.tsv"))

# positive-control threshold: median osteo-z of canonical markers (present in all atlases)
pos_thr <- median(meanz[gene %in% POS & n_atlas == length(atlases)]$mean_osteo_z)
cat(sprintf("positive-control median osteo-z (correspondence threshold) = %.3f\n", pos_thr))

# ---- A3: in-vivo-correspondent subprogram from RRA-744 up-genes --------------
# genes in RRA-744 that ARE expressed in the in-vivo osteoprogenitor compartment
sub_up <- meanz[gene %in% sigs$RRA744_pan_dental$up & mean_osteo_z >= pos_thr]$gene
cat(sprintf("RRA-744 up-genes above correspondence threshold: %d / %d\n",
            length(sub_up), length(sigs$RRA744_pan_dental$up)))
# re-project the subprogram (up-only) and, for reference, the full RRA-744, DFC, positive control
proj <- list()
for (tag in names(Ms)) { M <- Ms[[tag]]$M; os <- Ms[[tag]]$osteo; univ <- rownames(M)
  for (set in list(list("RRA744_correspondent_subprogram", sub_up, character(0)),
                   list("RRA744_full", sigs$RRA744_pan_dental$up, sigs$RRA744_pan_dental$down),
                   list("DFC_full", sigs$DE_DFC_GSE49007$up, sigs$DE_DFC_GSE49007$down),
                   list("POSctrl", POS, character(0)))) {
    r <- transfer_z(M, set[[2]], set[[3]], os, univ); if (is.null(r)) next
    proj[[paste(tag,set[[1]])]] <- data.table(atlas=tag, set=set[[1]], n_up=length(intersect(set[[2]],univ)),
      osteo_rank=r$rank, transferability_z=round(r$z,3), p=signif(r$p,3), localises = r$rank==1 & r$p<0.05) }
}
proj <- rbindlist(proj)
save_tsv(as.data.frame(proj), file.path(OUT, "rra_correspondent_subprogram_projection.tsv"))
writeLines(sort(sub_up), file.path(OUT, "rra744_correspondent_subprogram_genes.txt"))
cat("\n== subprogram re-projection (does an osteoprogenitor-expressed core inside RRA-744 transfer?) ==\n")
print(proj)
cat("\n[wrote dfc_rra_leading_edge.tsv, rra_correspondent_subprogram_projection.tsv, rra744_correspondent_subprogram_genes.txt]\n")
