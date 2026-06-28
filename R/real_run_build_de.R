# R/real_run_build_de.R  — REAL DATA ingestion + per-dataset DE for the curated
# clean osteogenic-induction pool (Step0+Step1 executed on actual GEO data).
# Produces 02_per_dataset_DE/<ID>_ranked.tsv and caches log-expr matrices for WGCNA.
.libPaths(c("F:/Rlib", .libPaths()))
suppressMessages({ library(data.table); library(limma); library(DESeq2); library(GEOquery) })
options(timeout = 1800)
SD <- "00_data/success_dental"; MD <- "00_data/meta"
dir.create("02_per_dataset_DE", showWarnings = FALSE); dir.create("03_rra_meta", showWarnings = FALSE)
log <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n")

# ---- ENSG -> symbol map, built from GSE226347 annotated files (have both) ----
ens2sym <- local({
  f <- list.files(file.path(SD, "GSE226347"), pattern = "N_DPSCs1.*\\.txt\\.gz$", full.names = TRUE)[1]
  d <- fread(f); setNames(d$GeneSymbol, sub("\\..*$", "", d$gene_id))
})
map_sym <- function(ensg) unname(ens2sym[sub("\\..*$", "", ensg)])

collapse_symbol <- function(mat, sym) {
  keep <- !is.na(sym) & sym != "" & sym != "---"
  mat <- mat[keep, , drop = FALSE]; sym <- sym[keep]
  o <- order(rowMeans(mat, na.rm = TRUE), decreasing = TRUE)
  mat <- mat[o, , drop = FALSE]; sym <- sym[o]
  mat <- mat[!duplicated(sym), , drop = FALSE]; rownames(mat) <- sym[!duplicated(sym)]
  mat
}

# ---------- per-dataset loaders: return list(counts|expr, group, kind) ----------
loaders <- list()

loaders$GSE99958 <- function() {
  d <- fread(file.path(SD, "GSE99958", "GSE99958_all.counts.exp.txt.gz"))
  sm <- d$GeneID; cnt <- as.matrix(d[, c("PDLC", "PDLC-D3", "PDLC-D7", "PDLC-D14")])
  rownames(cnt) <- sm
  cnt <- collapse_symbol(cnt, sm)
  list(counts = round(cnt), group = factor(c("control","osteo","osteo","osteo"), c("control","osteo")), kind = "count")
}

load_226347 <- function(which) {  # which = "DPSCs" or "GSCs"
  fs <- list.files(file.path(SD, "GSE226347"), pattern = paste0("_(N|OM)_", which, "\\d.*\\.txt\\.gz$"), full.names = TRUE)
  fs <- sort(fs)
  mats <- lapply(fs, function(f) { d <- fread(f); setNames(d[[3]], sub("\\..*$","",d$gene_id)) })  # col3 = raw count
  genes <- Reduce(intersect, lapply(mats, names))
  cnt <- sapply(mats, function(v) v[genes]); rownames(cnt) <- genes
  grp <- ifelse(grepl("_OM_", fs), "osteo", "control")
  sym <- map_sym(genes)
  cnt <- collapse_symbol(cnt, sym)
  list(counts = round(cnt), group = factor(grp, c("control","osteo")), kind = "count")
}
loaders$GSE226347_DPSC <- function() load_226347("DPSCs")
loaders$GSE226347_GSC  <- function() load_226347("GSCs")

loaders$GSE266257 <- function() {
  fs <- sort(list.files(file.path(SD, "GSE266257"), pattern = "counts\\.txt\\.gz$", full.names = TRUE))
  mats <- lapply(fs, function(f) { d <- fread(f, skip = 1); setNames(d[[ncol(d)]], sub("\\..*$","",d$Geneid)) })
  genes <- Reduce(intersect, lapply(mats, names))
  cnt <- sapply(mats, function(v) v[genes]); rownames(cnt) <- genes
  grp <- ifelse(grepl("_S4_", fs), "osteo", "control")   # S1=control, S4=Pi-treated(osteo)
  cnt <- collapse_symbol(cnt, map_sym(genes))
  list(counts = round(cnt), group = factor(grp, c("control","osteo")), kind = "count")
}

load_array <- function(acc, ctrl_kw, osteo_kw) {
  g <- getGEO(acc, destdir = MD, GSEMatrix = TRUE, AnnotGPL = TRUE, getGPL = TRUE)
  es <- if (is.list(g)) g[[1]] else g
  ex <- Biobase::exprs(es); fd <- Biobase::fData(es); tt <- tolower(Biobase::pData(es)$title)
  sc <- grep("^gene.?symbol$", colnames(fd), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(sc)) sc <- grep("symbol", colnames(fd), ignore.case = TRUE, value = TRUE)[1]
  sym <- as.character(fd[[sc]]); sym <- sub(" ?//.*$", "", sym)
  if (max(ex, na.rm = TRUE) > 100) ex <- log2(ex + 1)
  ex <- collapse_symbol(ex, sym)
  grp <- rep(NA_character_, length(tt))
  grp[grepl(osteo_kw, tt)] <- "osteo"; grp[grepl(ctrl_kw, tt)] <- "control"
  keep <- !is.na(grp); ex <- ex[, keep, drop = FALSE]; grp <- grp[keep]
  list(expr = ex, group = factor(grp, c("control","osteo")), kind = "array")
}
loaders$GSE49007  <- function() load_array("GSE49007",  "control",            "differentiat")
loaders$GSE159507 <- function() load_array("GSE159507", "cultured",           "osteogenic|differentiat")

# ---------------- DE per dataset ----------------
de_count <- function(cnt, grp) {
  cnt <- cnt[rowSums(cnt) >= 10, ]
  dds <- DESeqDataSetFromMatrix(cnt, data.frame(group = grp), ~group)
  dds <- DESeq(dds, quiet = TRUE)
  r <- as.data.frame(results(dds, contrast = c("group","osteo","control")))
  list(de = data.frame(gene = rownames(r), log2FC = r$log2FoldChange, stat = r$stat,
                       pval = r$pvalue, padj = r$padj), vst = assay(vst(dds, blind = TRUE)))
}
de_array <- function(ex, grp) {
  design <- model.matrix(~grp)
  fit <- eBayes(lmFit(ex, design), trend = TRUE, robust = TRUE)
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  list(de = data.frame(gene = rownames(tt), log2FC = tt$logFC, stat = tt$t,
                       pval = tt$P.Value, padj = tt$adj.P.Val), vst = ex)
}

ids <- c("GSE99958","GSE226347_DPSC","GSE226347_GSC","GSE266257","GSE49007","GSE159507")
exprs_cache <- list(); pheno_cache <- list(); summ <- list()
for (id in ids) {
  log("==== ", id, " ====")
  dat <- tryCatch(loaders[[id]](), error = function(e) { log("LOAD FAIL: ", conditionMessage(e)); NULL })
  if (is.null(dat)) next
  res <- tryCatch({
    if (dat$kind == "count") de_count(dat$counts, dat$group) else de_array(dat$expr, dat$group)
  }, error = function(e) { log("DE FAIL: ", conditionMessage(e)); NULL })
  if (is.null(res)) next
  de <- res$de[!is.na(res$de$stat) & res$de$gene != "", ]
  de <- de[order(-de$stat), ]; de$rank <- seq_len(nrow(de))
  write.table(de, file.path("02_per_dataset_DE", paste0(id, "_ranked.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  exprs_cache[[id]] <- res$vst
  pheno_cache[[id]] <- data.frame(sample = colnames(res$vst), dataset = id,
                                  group = as.character(dat$group)[seq_len(ncol(res$vst))], stringsAsFactors = FALSE)
  ng <- sum(de$padj < 0.05 & abs(de$log2FC) > 1, na.rm = TRUE)
  summ[[id]] <- data.frame(id = id, n_ctrl = sum(dat$group=="control"), n_osteo = sum(dat$group=="osteo"),
                           n_genes = nrow(de), n_DEG = ng, kind = dat$kind)
  log(id, ": ", nrow(de), " genes, ", ng, " DEG (FDR<.05,|lfc|>1)")
}
saveRDS(exprs_cache, "03_rra_meta/exprs_cache.rds")
saveRDS(pheno_cache, "03_rra_meta/pheno_cache.rds")
st <- do.call(rbind, summ); print(st)
write.table(st, "02_per_dataset_DE/_dataset_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
log("DONE. ranked files + exprs_cache.rds written.")
