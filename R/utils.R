# R/utils.R — shared helpers for the maxillary-regeneration dry-lab pipeline.
# Source this at the top of every Step script: source("R/utils.R")

suppressWarnings(suppressMessages({
  needed <- c("data.table")
  for (p in needed) if (!requireNamespace(p, quietly = TRUE))
    message("[utils] NOTE: package '", p, "' not installed; see env/install_packages.R")
}))

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a

## ---- logging ----
.log <- function(..., level = "INFO") {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), " ", level, "] ", paste0(..., collapse = ""))
  message(msg)
}
log_step <- function(step, msg) .log(sprintf("Step%s | %s", step, msg))

## ---- registry ----
load_registry <- function(path = "config/datasets.tsv") {
  stopifnot(file.exists(path))
  reg <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  reg$include_in_rra <- tolower(reg$include_in_rra) %in% c("yes", "true", "1")
  reg$is_external_validation <- tolower(reg$is_external_validation) %in% c("yes", "true", "1")
  reg
}

rra_training_sets <- function(reg) reg[reg$include_in_rra & reg$modality == "mrna", , drop = FALSE]

## ---- GEO download (cached) ----
# Returns list(eset=ExpressionSet or matrix, meta=pdata). Caches raw download under 00_data/<role>/.
fetch_geo <- function(accession, role = "misc", destdir = "00_data", force = FALSE) {
  if (!requireNamespace("GEOquery", quietly = TRUE))
    stop("GEOquery required. Run env/install_packages.R")
  dir <- file.path(destdir, role)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  cache <- file.path(dir, paste0(accession, "_eset.rds"))
  if (file.exists(cache) && !force) {
    .log("cache hit: ", accession)
    return(readRDS(cache))
  }
  .log("downloading ", accession, " -> ", dir)
  g <- GEOquery::getGEO(accession, destdir = dir, GSEMatrix = TRUE, AnnotGPL = TRUE)
  eset <- if (is.list(g)) g[[1]] else g
  saveRDS(eset, cache)
  eset
}

## ---- expression matrix + symbol collapse ----
# Detects whether values are in log space; collapses probes->symbol by max mean.
get_expr_symbol <- function(eset) {
  if (!requireNamespace("Biobase", quietly = TRUE)) stop("Biobase required.")
  m <- Biobase::exprs(eset)
  fd <- Biobase::fData(eset)
  # find a gene-symbol column
  sym_col <- grep("symbol|gene.symbol|gene_assignment|GENE_SYMBOL", colnames(fd), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(sym_col)) {
    .log("no symbol column in fData; returning probe-level matrix", level = "WARN")
    return(m)
  }
  sym <- as.character(fd[[sym_col]])
  sym <- sub(" ?//.*$", "", sym)         # array_comprehensive gene_assignment cleanup
  keep <- !is.na(sym) & sym != "" & sym != "---"
  m <- m[keep, , drop = FALSE]; sym <- sym[keep]
  # log-transform if needed
  qx <- as.numeric(stats::quantile(m, c(0, .25, .5, .75, .99, 1), na.rm = TRUE))
  logged <- (qx[5] - qx[1] < 50) || (qx[6] - qx[2] < 100 && qx[2] > 0)
  if (!logged && max(m, na.rm = TRUE) > 100) { m[m < 0] <- 0; m <- log2(m + 1); .log("applied log2 transform") }
  # collapse by max-mean probe per symbol
  o <- order(rowMeans(m, na.rm = TRUE), decreasing = TRUE)
  m <- m[o, , drop = FALSE]; sym <- sym[o]
  m <- m[!duplicated(sym), , drop = FALSE]
  rownames(m) <- sym[!duplicated(sym)]
  m
}

## ---- ortholog mapping (mouse <-> human) ----
# strict: quantitative steps (Step4-7) SHOULD call with strict=TRUE so a missing
#   babelgene install hard-fails instead of silently degrading to a toupper()
#   pseudo-mapping. Default strict=FALSE preserves the original signature/behaviour;
#   the toupper fallback is acceptable ONLY for illustrative scRNA visualisation.
# Mouse->human orthology is not 1:1; babelgene can return multiple human symbols
#   per mouse gene (and vice versa). We enforce a deterministic 1:1 map: drop
#   non-1:1 mouse genes entirely (returned as NA) rather than silently keeping the
#   first hit, so a downstream quantitative join never depends on row order.
to_human_symbols <- function(genes, from = c("mouse", "human"), strict = FALSE) {
  from <- match.arg(from)
  if (from == "human") return(genes)
  if (!requireNamespace("babelgene", quietly = TRUE)) {
    if (strict)
      stop("babelgene required for ortholog mapping (strict=TRUE). Run env/install_packages.R")
    .log("babelgene not installed; returning toupper() fallback (NOT publication-grade; ",
         "illustrative only — quantitative steps must call strict=TRUE)", level = "WARN")
    return(toupper(genes))
  }
  # NOTE: babelgene's `human` flag declares whether the INPUT genes are human.
  # The input here is MOUSE, so human = FALSE (human = TRUE silently maps only the
  # handful of mouse symbols that coincide with a human symbol, ~15 of ~23k).
  ort <- babelgene::orthologs(genes = genes, species = "mouse", human = FALSE)
  ort <- ort[!is.na(ort$symbol) & !is.na(ort$human_symbol), , drop = FALSE]
  # Drop ambiguous mappings: keep only mouse symbols with exactly one distinct
  # human symbol (and human symbols claimed by exactly one mouse symbol). Among
  # exact duplicate rows order alphabetically first for determinism.
  ort <- ort[order(ort$symbol, ort$human_symbol), , drop = FALSE]
  ort <- ort[!duplicated(ort[, c("symbol", "human_symbol")]), , drop = FALSE]
  n_human_per_mouse <- tapply(ort$human_symbol, ort$symbol, function(x) length(unique(x)))
  n_mouse_per_human <- tapply(ort$symbol, ort$human_symbol, function(x) length(unique(x)))
  ort_11 <- ort[n_human_per_mouse[ort$symbol] == 1L & n_mouse_per_human[ort$human_symbol] == 1L, , drop = FALSE]
  n_dropped <- length(unique(genes)) - length(unique(ort_11$symbol))
  if (n_dropped > 0)
    .log("to_human_symbols: ", n_dropped, " mouse gene(s) had non-1:1 or no human ortholog; returned NA",
         level = "WARN")
  unname(stats::setNames(ort_11$human_symbol, ort_11$symbol)[genes])
}

## ---- io ----
save_rds <- function(obj, path) { dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE); saveRDS(obj, path); .log("saved ", path) }
save_tsv <- function(df, path) { dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE); .log("saved ", path) }

dump_session <- function(path = "supp/sessionInfo.txt") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(capture.output(utils::sessionInfo()), path)
}

# Convenience: load params + registry once.
init_pipeline <- function() {
  source("config/params.R", local = FALSE)
  assign("REG", load_registry(), envir = .GlobalEnv)
  .log("pipeline initialized; ", nrow(REG), " datasets in registry; ",
       nrow(rra_training_sets(REG)), " in RRA training pool")
  invisible(TRUE)
}
