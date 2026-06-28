# R/step08_position_context.R
# Step 8 — Jaw position-specific context: does the regeneration-competent program
# sit in a baseline position that explains the jaw's osteogenic disadvantage?
#
# Idea: jaw (neural-crest) bone is documented to be osteogenically weaker than
# appendicular bone. For each core-program gene we ask whether its baseline jaw-vs-
# other-bone direction is CONSISTENT with that weakness:
#   core osteo-UP gene that is LOWER in jaw   -> consistent (weakness explained)
#   core osteo-DOWN gene that is HIGHER in jaw -> consistent
# Magnitude = |log2FC(jaw vs other bone)|.
#
# Data: GSE58474 (human mandible vs iliac osteoblasts), GSE30167 (mouse jaw/alveolar
# vs long bone; ortholog-mapped). These are baseline comparisons, NOT osteo-induction.
#
# Output (consumed by Step10 jaw_context dimension):
#   05_projection/jaw_position_context.tsv
#     columns: gene, jaw_weak_direction_consistent(logical), abs_jaw_effect,
#              jaw_log2FC_human, jaw_log2FC_mouse, n_sources
#
# Honesty: baseline position is a SUSCEPTIBILITY context, not a regeneration assay;
# treated as supporting evidence, not a standalone claim.

source("R/utils.R"); init_pipeline()
suppressMessages(library(limma))

core <- utils::read.delim(file.path(PARAMS$dir_wgcna, "core_program.tsv"), stringsAsFactors = FALSE)
stopifnot(all(c("gene", "direction") %in% colnames(core)))
succ_dir <- stats::setNames(core$direction, core$gene)

# datasets: jaw label vs other-bone label
POS <- data.frame(
  accession = c("GSE58474",        "GSE30167"),
  organism  = c("human",           "mouse"),
  jaw_kw    = c("mandib|jaw|alveol","jaw|mandib|alveol"),
  other_kw  = c("iliac|ilium|hip",  "long|femur|tibia|limb|appendic"),
  stringsAsFactors = FALSE
)

sample_text <- function(eset) {
  pd <- Biobase::pData(eset)
  cols <- grep("characteristics|title|source_name|description|:ch1|tissue|site",
               colnames(pd), ignore.case = TRUE, value = TRUE)
  if (!length(cols)) cols <- colnames(pd)
  desc <- apply(pd[, cols, drop = FALSE], 1, function(r) paste(na.omit(as.character(r)), collapse = " | "))
  stats::setNames(tolower(desc), rownames(pd))
}

# returns named vector of log2FC(jaw - other) per gene (human symbols), or NULL
jaw_logfc <- function(acc, organism, jaw_kw, other_kw) {
  ri <- which(REG$accession == acc); role <- if (length(ri)) REG$role[ri] else "position"
  eset <- tryCatch(fetch_geo(acc, role), error = function(e) { .log("fetch ", acc, " failed: ", conditionMessage(e), level = "WARN"); NULL })
  if (is.null(eset)) return(NULL)
  expr <- tryCatch(get_expr_symbol(eset), error = function(e) NULL)
  if (is.null(expr) || ncol(expr) < 2) { .log(acc, ": no usable matrix", level = "WARN"); return(NULL) }
  txt <- sample_text(eset)[colnames(expr)]
  grp <- rep(NA_character_, length(txt))
  grp[grepl(jaw_kw, txt)]   <- "jaw"
  grp[grepl(other_kw, txt)] <- "other"
  if (sum(grp == "jaw", na.rm = TRUE) < 1 || sum(grp == "other", na.rm = TRUE) < 1) {
    .log(acc, ": cannot auto-assign jaw vs other (manual confirm needed) -> SKIP", level = "WARN")
    return(NULL)
  }
  keep <- !is.na(grp); expr <- expr[, keep, drop = FALSE]; grp <- grp[keep]
  design <- model.matrix(~factor(grp, levels = c("other", "jaw")))
  fit <- eBayes(lmFit(expr, design), trend = TRUE, robust = TRUE)
  lfc <- topTable(fit, coef = 2, number = Inf, sort.by = "none")$logFC
  names(lfc) <- rownames(expr)
  if (organism == "mouse") {
    hs <- to_human_symbols(names(lfc), from = "mouse")
    keep <- !is.na(hs) & hs != ""
    lfc <- lfc[keep]; names(lfc) <- hs[keep]
    lfc <- lfc[!duplicated(names(lfc))]
  }
  .log(acc, ": jaw-vs-other log2FC for ", length(lfc), " genes")
  lfc
}

h <- jaw_logfc("GSE58474", "human", POS$jaw_kw[1], POS$other_kw[1])
m <- jaw_logfc("GSE30167", "mouse", POS$jaw_kw[2], POS$other_kw[2])

genes <- unique(core$gene)
hf <- if (!is.null(h)) h[genes] else stats::setNames(rep(NA_real_, length(genes)), genes)
mf <- if (!is.null(m)) m[genes] else stats::setNames(rep(NA_real_, length(genes)), genes)

# combined jaw effect: prefer human, fall back to mouse; consistency uses available sign
combined <- ifelse(!is.na(hf), hf, mf)
sd <- succ_dir[genes]
# consistent with jaw osteogenic weakness:
#   up-gene lower in jaw (combined<0)  OR  down-gene higher in jaw (combined>0)
consistent <- (sd == "up" & combined < 0) | (sd == "down" & combined > 0)
consistent[is.na(consistent)] <- FALSE

out <- data.frame(
  gene = genes,
  jaw_weak_direction_consistent = consistent,
  abs_jaw_effect = ifelse(is.na(combined), 0, abs(combined)),
  jaw_log2FC_human = unname(hf),
  jaw_log2FC_mouse = unname(mf),
  n_sources = (!is.na(hf)) + (!is.na(mf)),
  stringsAsFactors = FALSE
)
out <- out[order(-out$jaw_weak_direction_consistent, -out$abs_jaw_effect), ]
save_tsv(out, file.path(PARAMS$dir_proj, "jaw_position_context.tsv"))

cat(sprintf("\nStep8 position context: %d core genes; %d direction-consistent with jaw osteogenic weakness (sources: human=%s, mouse=%s).\n",
            nrow(out), sum(out$jaw_weak_direction_consistent),
            !is.null(h), !is.null(m)))
cat("NOTE: baseline jaw-vs-other-bone position = susceptibility context (supporting evidence), not a regeneration assay.\n")
dump_session()
