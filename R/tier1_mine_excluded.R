# R/tier1_mine_excluded.R
# ============================================================================
# Tier 1-A (extension): re-mine independent osteogenic-CONTEXT signatures from
# series that comparison-arm governance EXCLUDED from the RRA training pool
# (perturbation x osteogenesis designs). These are the hardest test of
# systematicity: signatures distilled from entirely separate labs / inducers.
# Appends to 10_transferability/signatures.rds and rewrites the registry.
#
# Datasets mined (osteogenic-context contrast, treated vs control):
#   GSE286540  SCAP   simvastatin-induced osteogenesis (ready DE table)
#   GSE271641  DPSC   cannabidiol odonto/osteogenic    (FPKM, T vs N)
#   GSE236009  PDLSC  CALD1-knockdown osteogenesis      (FPKM, shCALD1 vs shNC)
# ============================================================================
.libPaths(c("F:/Rlib", .libPaths()))
suppressWarnings(suppressMessages({ library(data.table); library(limma) }))
setDTthreads(1)
source("R/utils.R")
OUT <- "10_transferability"
CAP <- 250L; MINDIR <- 20L
sigs <- readRDS(file.path(OUT, "signatures.rds"))

# limma signature from a log2 expression matrix (rows=symbol) + 2 groups
mk_limma <- function(mat, treated, control, name, src, note) {
  mat <- mat[rowSums(is.finite(mat)) == ncol(mat), , drop = FALSE]
  grp <- factor(c(rep("trt", length(treated)), rep("ctl", length(control))), levels = c("ctl","trt"))
  X <- mat[, c(treated, control), drop = FALSE]
  des <- model.matrix(~ grp)
  fit <- eBayes(lmFit(X, des), trend = TRUE, robust = TRUE)
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "t")
  tt$gene <- rownames(tt)
  up <- tt$gene[tt$adj.P.Val < 0.05 & tt$logFC >  1]
  dn <- tt$gene[tt$adj.P.Val < 0.05 & tt$logFC < -1]
  if (length(up) < MINDIR) up <- tt$gene[order(-tt$t)][seq_len(min(MINDIR, nrow(tt)))]
  if (length(dn) < MINDIR) dn <- tt$gene[order(tt$t)][seq_len(min(MINDIR, nrow(tt)))]
  list(name = name, source_type = "excluded_perturbation_DE", organism = "human",
       up = head(unique(up), CAP), down = head(unique(dn), CAP), citation = src, notes = note)
}
collapse_sym <- function(dt, sym_col, val_cols) {
  m <- as.matrix(dt[, ..val_cols]); suppressWarnings(storage.mode(m) <- "double")
  sym <- as.character(dt[[sym_col]]); ok <- !is.na(sym) & nzchar(sym) & rowSums(is.finite(m)) == ncol(m)
  m <- log2(m[ok, , drop = FALSE] + 1); sym <- sym[ok]
  o <- order(rowMeans(m), decreasing = TRUE); m <- m[o, ]; sym <- sym[o]
  m <- m[!duplicated(sym), , drop = FALSE]; rownames(m) <- sym[!duplicated(sym)]; m
}

## --- GSE286540 SCAP simvastatin: ready DE table (GeneSymbol/log2FC/padj) -----
try({
  de <- fread("00_data/success_dental/GSE286540/GSE286540_Group_SK-VS-SV_DE_significant_anno.txt.gz")
  de <- de[!is.na(GeneSymbol) & nzchar(GeneSymbol)]
  up <- de[padj < 0.05 & log2FoldChange >  1][order(-log2FoldChange), unique(GeneSymbol)]
  dn <- de[padj < 0.05 & log2FoldChange < -1][order( log2FoldChange), unique(GeneSymbol)]
  if (length(up) < MINDIR) up <- de[order(-log2FoldChange), unique(GeneSymbol)][seq_len(MINDIR)]
  if (length(dn) < MINDIR) dn <- de[order( log2FoldChange), unique(GeneSymbol)][seq_len(MINDIR)]
  sigs[["EX_SCAP_simvastatin_GSE286540"]] <- list(name="EX_SCAP_simvastatin_GSE286540",
    source_type="excluded_perturbation_DE", organism="human",
    up=head(up,CAP), down=head(dn,CAP), citation="GSE286540 SCAP simvastatin-osteo (excluded)",
    notes="re-mined from author DE table; treated(sim) vs control")
  .log("mined GSE286540: ", length(up), " up / ", length(dn), " down")
})

## --- GSE271641 DPSC cannabidiol: FPKM T vs N --------------------------------
try({
  dt <- fread("00_data/success_dental/GSE271641/GSE271641_all.fpkm_anno.txt.gz")
  fk <- grep("_FPKM$", names(dt), value = TRUE)
  trt <- grep("^T[0-9].*_FPKM$", fk, value = TRUE); ctl <- grep("Control.*_FPKM$", fk, value = TRUE)
  m <- collapse_sym(dt, "GeneSymbol", c(trt, ctl))
  sigs[["EX_DPSC_cannabidiol_GSE271641"]] <- mk_limma(m, trt, ctl,
    "EX_DPSC_cannabidiol_GSE271641", "GSE271641 DPSC CBD odonto/osteo (excluded)", "FPKM limma, CBD vs control")
  .log("mined GSE271641: ", length(sigs[["EX_DPSC_cannabidiol_GSE271641"]]$up), " up / ",
       length(sigs[["EX_DPSC_cannabidiol_GSE271641"]]$down), " down")
})

## --- GSE236009 PDLSC CALD1 knockdown: FPKM shCALD1 vs shNC ------------------
try({
  dt <- fread("00_data/success_dental/GSE236009/GSE236009_ExpressionProfile.txt.gz")
  trt <- c("shCALD1_1","shCALD1_2","shCALD1_3"); ctl <- c("shNC_1","shNC_2","shNC_3")
  m <- collapse_sym(dt, "gene_name", c(trt, ctl))
  sigs[["EX_PDLSC_CALD1kd_GSE236009"]] <- mk_limma(m, trt, ctl,
    "EX_PDLSC_CALD1kd_GSE236009", "GSE236009 PDLSC CALD1-KD osteo (excluded)", "FPKM limma, shCALD1 vs shNC")
  .log("mined GSE236009: ", length(sigs[["EX_PDLSC_CALD1kd_GSE236009"]]$up), " up / ",
       length(sigs[["EX_PDLSC_CALD1kd_GSE236009"]]$down), " down")
})

## --- rewrite outputs --------------------------------------------------------
save_rds(sigs, file.path(OUT, "signatures.rds"))
reg <- rbindlist(lapply(sigs, function(s) data.table(signature=s$name, source_type=s$source_type,
  organism=s$organism, n_up=length(s$up), n_down=length(s$down), n_total=length(s$up)+length(s$down),
  citation=s$citation, notes=s$notes)))
save_tsv(as.data.frame(reg), file.path(OUT, "signature_registry.tsv"))
long <- rbindlist(lapply(sigs, function(s) rbindlist(list(
  if (length(s$up))   data.table(signature=s$name, direction="up",   gene=s$up)   else NULL,
  if (length(s$down)) data.table(signature=s$name, direction="down", gene=s$down) else NULL))))
save_tsv(as.data.frame(long), file.path(OUT, "signatures_long.tsv"))
.log("panel now has ", length(sigs), " signatures")
print(reg[, .(signature, source_type, n_up, n_down)])
