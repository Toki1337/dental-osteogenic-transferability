# R/step01_per_dataset_DE.R
# Step 1 — Per-dataset differential expression -> ranked gene lists for RRA.
# RNA-seq with raw counts  -> DESeq2 ; otherwise (arrays / normalized matrices) -> limma.
# Output per usable dataset: 02_per_dataset_DE/<ACC>_ranked.tsv  (gene, log2FC, stat, pval, padj, rank)
# Ranking direction: osteo vs control (positive log2FC = up in osteogenic induction).

source("R/utils.R"); init_pipeline()

arms_path <- file.path(PARAMS$dir_curation, "sample_arms_final.tsv")
if (!file.exists(arms_path)) stop("Run Step0 and produce sample_arms_final.tsv first.")
arms <- utils::read.delim(arms_path, stringsAsFactors = FALSE)
# Consume the human purity gate: keep only samples whose arm_clean is yes/true.
# Fail-safe: blank/NA arm_clean is treated as NOT clean (excluded).
arms <- arms[arms$final_group %in% c("osteo", "control") &
               (tolower(arms$arm_clean) %in% c("yes", "y", "true", "1")), ]
usable <- read.delim(file.path(PARAMS$dir_curation, "rra_usable_datasets.tsv"))$accession

# Try to locate raw-count supplementary for a GSE (returns matrix or NULL).
load_counts <- function(acc, role) {
  files <- list.files(file.path(PARAMS$dir_data, role), pattern = paste0(acc, ".*(count|raw|featurecounts).*\\.(txt|tsv|csv|gz)$"),
                      full.names = TRUE, ignore.case = TRUE)
  if (!length(files)) return(NULL)
  m <- tryCatch(as.matrix(data.table::fread(files[1]), rownames = 1), error = function(e) NULL)
  if (is.null(m) || !is.numeric(m[1, 1])) return(NULL)
  m
}

de_deseq2 <- function(counts, grp) {
  suppressMessages(library(DESeq2))
  counts <- round(counts[rowSums(counts) >= PARAMS$min_count, ])
  cd <- data.frame(group = factor(grp, levels = c("control", "osteo")))
  dds <- DESeqDataSetFromMatrix(counts, cd, ~group)
  dds <- DESeq(dds, quiet = TRUE)
  res <- as.data.frame(results(dds, contrast = c("group", "osteo", "control")))
  data.frame(gene = rownames(res), log2FC = res$log2FoldChange, stat = res$stat,
             pval = res$pvalue, padj = res$padj, stringsAsFactors = FALSE)
}

de_limma <- function(expr, grp) {
  suppressMessages(library(limma))
  grp <- factor(grp, levels = c("control", "osteo"))
  design <- model.matrix(~grp)
  fit <- eBayes(lmFit(expr, design), trend = TRUE, robust = TRUE)
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  data.frame(gene = rownames(tt), log2FC = tt$logFC, stat = tt$t,
             pval = tt$P.Value, padj = tt$adj.P.Val, stringsAsFactors = FALSE)
}

for (acc in usable) {
  ri <- which(REG$accession == acc); role <- REG$role[ri]; plat <- REG$platform[ri]
  a <- arms[arms$accession == acc, ]
  eset <- fetch_geo(acc, role)
  pd_ids <- rownames(Biobase::pData(eset))
  grp <- a$final_group[match(pd_ids, a$gsm)]
  keep <- !is.na(grp); grp <- grp[keep]
  if (length(unique(grp)) < 2) { .log("skip ", acc, " (need both arms)", level = "WARN"); next }

  counts <- if (grepl("rnaseq", plat)) load_counts(acc, role) else NULL
  if (!is.null(counts)) {
    counts <- counts[, pd_ids[keep], drop = FALSE]
    res <- tryCatch(de_deseq2(counts, grp), error = function(e) { .log(acc, " DESeq2 failed -> limma: ", conditionMessage(e), level = "WARN"); NULL })
  } else res <- NULL
  if (is.null(res)) {
    expr <- get_expr_symbol(eset)[, pd_ids[keep], drop = FALSE]
    res <- de_limma(expr, grp)
  }
  res <- res[!is.na(res$stat) & !is.na(res$gene) & res$gene != "", ]
  res <- res[order(-res$stat), ]
  res$rank <- seq_len(nrow(res))
  save_tsv(res, file.path(PARAMS$dir_de, paste0(acc, "_ranked.tsv")))
  .log(acc, ": ", nrow(res), " genes ranked (",
       sum(res$padj < PARAMS$de_fdr & abs(res$log2FC) > PARAMS$de_lfc, na.rm = TRUE), " DEG)")
}
dump_session()
.log("Step1 done. Ranked lists in ", PARAMS$dir_de)
