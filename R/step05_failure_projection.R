# R/step05_failure_projection.R
# Step 5 — Failure-side projection + the "failure-specific" dysregulated program.
#
# IMPORTANT positioning (read before interpreting any number):
#   The failure side is DATA-THIN, especially in humans. The mouse mandibular axis
#   (GSE295106 control-vs-BRONJ scRNA; GSE269255 irradiated ORNJ stroma) is the main
#   evidence; human GSE7116 (MM-ONJ blood, NO healthy ctrl, MM/chemo confound) and
#   GSE303003 (n=2 MRONJ scRNA, scoop risk) are SUPPORTIVE / illustrative only and are
#   NOT used for quantitative cohort inference. ALL Step5 conclusions are
#   HYPOTHESIS-GENERATING, not confirmatory.
#
# What we do:
#   1) Build pseudobulk (or use bulk) per failure dataset, homogenize to human symbols.
#   2) ssGSEA-score four functional modules (osteoblast diff / osteoclast diff /
#      inflammatory response / angiogenesis from MSigDB GOBP via msigdbr) AND the
#      core up/down sets, to contrast success vs failure module activity.
#   3) Identify "failure-specific" core genes: DOWN in failure but UP in success
#      (direction reversal); mouse results ortholog-mapped to human and intersected
#      with any human-side support for a conservative call.
#
# Output:
#   05_projection/failure_module_scores.tsv     (per-sample/pseudobulk module + core scores)
#   05_projection/failure_specific_program.tsv  (gene, dir_success, dir_failure, reversed, conserved_human_mouse)

source("R/utils.R"); init_pipeline()
suppressMessages({ library(GSVA) })

# ---- core program (success-side directions) ---------------------------------
core <- utils::read.delim(file.path(PARAMS$dir_wgcna, "core_program.tsv"), stringsAsFactors = FALSE)
stopifnot(all(c("gene", "direction") %in% colnames(core)))
core$dir_success <- core$direction
core_sets_human <- list(
  core_up   = unique(core$gene[core$direction == "up"]),
  core_down = unique(core$gene[core$direction == "down"])
)

# ---- functional module gene sets (MSigDB GOBP via msigdbr) ------------------
# We pull human GOBP gene sets; failure expression is homogenized to human first.
module_terms <- c(
  osteoblast_diff = "GOBP_OSTEOBLAST_DIFFERENTIATION",
  osteoclast_diff = "GOBP_OSTEOCLAST_DIFFERENTIATION",
  inflammatory    = "GOBP_INFLAMMATORY_RESPONSE",
  angiogenesis    = "GOBP_ANGIOGENESIS"
)
module_sets <- tryCatch({
  if (!requireNamespace("msigdbr", quietly = TRUE)) stop("msigdbr not installed")
  m <- msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
  out <- lapply(module_terms, function(g) unique(m$gene_symbol[m$gs_name == g]))
  names(out) <- names(module_terms)
  out <- out[vapply(out, length, integer(1)) > 0]
  out
}, error = function(e) {
  # TODO(manual): if msigdbr unavailable, supply these GOBP sets manually (gene symbols):
  #   GOBP_OSTEOBLAST_DIFFERENTIATION, GOBP_OSTEOCLAST_DIFFERENTIATION,
  #   GOBP_INFLAMMATORY_RESPONSE, GOBP_ANGIOGENESIS  (download from MSigDB, human symbols)
  .log("msigdbr unavailable (", conditionMessage(e), "); module sets EMPTY -> see TODO. ",
       "Core up/down still scored.", level = "WARN")
  list()
})
gene_sets <- c(core_sets_human, module_sets)
.log("scoring ", length(gene_sets), " gene sets (", paste(names(gene_sets), collapse = ", "), ")")

# ---- failure datasets -------------------------------------------------------
# tier: "primary" = mouse mandibular axis (quantitative-ish, still hypothesis-generating);
#       "support" = human thin/confounded (illustrative only, flagged, NOT in reversal call by default).
failure_targets <- data.frame(
  accession = c("GSE295106",          "GSE269255",        "GSE7116",          "GSE303003"),
  organism  = c("mouse",              "mouse",            "human",            "human"),
  modality  = c("scrnaseq",           "scrnaseq",         "array",            "scrnaseq"),
  context   = c("BRONJ_mandible_sc",  "ORNJ_mandible_sc", "MM_ONJ_blood",     "MRONJ_granulation_sc"),
  tier      = c("primary",            "primary",          "support",          "support"),
  stringsAsFactors = FALSE
)

# disease-vs-control inference from pData text
disease_kw <- "bronj|mronj|ornj|onj|bisphosphon|zoledron|alendron|irradiat|radiation|necros|fail|disease|granulation|patient"
control_kw <- "control|ctrl|sham|healthy|normal|vehicle|untreated|wild[- ]?type|\\bWT\\b|non[- ]?irradiat|baseline|cyst"
guess_disease_group <- function(text) {
  t <- tolower(text)
  has_d <- grepl(disease_kw, t, perl = TRUE)
  has_c <- grepl(control_kw, t, perl = TRUE)
  if (has_d && !has_c) return("failure")
  if (has_c && !has_d) return("control")
  NA_character_
}
sample_desc <- function(eset) {
  pd <- Biobase::pData(eset)
  cols <- grep("characteristics|title|source_name|description|:ch1|treatment|group|genotype|agent|disease",
               colnames(pd), ignore.case = TRUE, value = TRUE)
  if (!length(cols)) cols <- colnames(pd)
  desc <- apply(pd[, cols, drop = FALSE], 1, function(r) paste(na.omit(as.character(r)), collapse = " | "))
  stats::setNames(desc, rownames(pd))
}

# Homogenize an expression matrix (symbol x sample) to human symbols.
to_human_matrix <- function(expr, organism) {
  if (organism == "mouse") {
    hs <- to_human_symbols(rownames(expr), from = "mouse")
    keep <- !is.na(hs) & hs != ""
    expr <- expr[keep, , drop = FALSE]; rownames(expr) <- hs[keep]
    expr <- expr[!duplicated(rownames(expr)), , drop = FALSE]
  }
  expr
}

# Build a per-group pseudobulk / bulk matrix for one failure dataset.
# For scRNA we cannot assume the Step6 Seurat object exists here, so we degrade
# gracefully: if get_expr_symbol yields a usable bulk-like matrix (e.g. authors
# deposited a processed matrix), we use it; otherwise we emit a TODO and skip.
load_failure_expr <- function(acc, role, organism) {
  eset <- tryCatch(fetch_geo(acc, role), error = function(e) { .log("fetch fail ", acc, ": ", conditionMessage(e), level = "WARN"); NULL })
  if (is.null(eset)) return(NULL)
  expr <- tryCatch(get_expr_symbol(eset), error = function(e) NULL)
  desc <- if (!is.null(eset)) sample_desc(eset) else NULL
  if (is.null(expr) || !nrow(expr) || ncol(expr) < 2) {
    .log(acc, ": no usable bulk/processed matrix from GEO series (likely 10x raw only). ",
         "TODO: pseudobulk from Step6 Seurat object (06_scrna/", acc, "_seurat.rds) by group, then re-run.",
         level = "WARN")
    return(NULL)
  }
  expr <- to_human_matrix(expr, organism)
  list(expr = expr, desc = desc[colnames(expr)])
}

score_failure <- function(i) {
  acc <- failure_targets$accession[i]; organism <- failure_targets$organism[i]
  context <- failure_targets$context[i]; tier <- failure_targets$tier[i]
  ri <- which(REG$accession == acc); role <- if (length(ri)) REG$role[ri] else "failure"
  fx <- load_failure_expr(acc, role, organism)
  if (is.null(fx)) return(NULL)

  par <- ssgseaParam(as.matrix(fx$expr), gene_sets, normalize = TRUE)
  es  <- tryCatch(as.matrix(gsva(par)),
                  error = function(e) as.matrix(gsva(as.matrix(fx$expr), gene_sets,
                            method = PARAMS$ssgsea_method, ssgsea.norm = TRUE, verbose = FALSE)))
  grp <- vapply(fx$desc %||% rep(NA_character_, ncol(fx$expr)), guess_disease_group, character(1))

  long <- do.call(rbind, lapply(rownames(es), function(set) data.frame(
    accession = acc, context = context, tier = tier,
    sample = colnames(es), gene_set = set, score = es[set, ],
    inferred_group = ifelse(is.na(grp), "unassigned", grp),
    group_confirmed = FALSE,
    stringsAsFactors = FALSE)))
  long
}

mod_scores <- do.call(rbind, lapply(seq_len(nrow(failure_targets)), score_failure))
if (is.null(mod_scores)) mod_scores <- data.frame()
save_tsv(mod_scores, file.path(PARAMS$dir_proj, "failure_module_scores.tsv"))

# ---- failure-specific program: reversal (success UP, failure DOWN) ----------
# Per primary (mouse) failure dataset, compute the per-gene failure direction of the
# core genes as sign(mean(failure) - mean(control)) on human-homogenized expression.
# A core gene is "failure-specific" if its success direction is UP but failure direction
# is DOWN (reversal). "conserved_human_mouse" = the reversal also has any human support.
failure_dir_for <- function(i) {
  acc <- failure_targets$accession[i]; organism <- failure_targets$organism[i]
  ri <- which(REG$accession == acc); role <- if (length(ri)) REG$role[ri] else "failure"
  fx <- load_failure_expr(acc, role, organism)
  if (is.null(fx)) return(NULL)
  grp <- vapply(fx$desc %||% rep(NA_character_, ncol(fx$expr)), guess_disease_group, character(1))
  if (sum(grp == "failure", na.rm = TRUE) < 1 || sum(grp == "control", na.rm = TRUE) < 1) {
    .log(acc, ": cannot infer both failure & control arms -> no direction (manual confirm needed)", level = "WARN")
    return(NULL)
  }
  genes <- intersect(rownames(fx$expr), core$gene)
  if (!length(genes)) return(NULL)
  fmean <- rowMeans(fx$expr[genes, grp == "failure", drop = FALSE], na.rm = TRUE)
  cmean <- rowMeans(fx$expr[genes, grp == "control", drop = FALSE], na.rm = TRUE)
  d <- fmean - cmean
  data.frame(accession = acc, organism = organism, gene = genes,
             dir_failure = ifelse(d > 0, "up", "down"), delta = d,
             stringsAsFactors = FALSE)
}

dir_tabs <- do.call(rbind, lapply(seq_len(nrow(failure_targets)), failure_dir_for))

# success direction per gene (from core program)
succ_dir <- stats::setNames(core$dir_success, core$gene)

# schema for the failure-specific program (Step7/Step10 consume gene, failure_direction, abs_effect, conserved_human_mouse)
build_empty_fsp <- function() data.frame(
  gene = character(0), dir_success = character(0), dir_failure = character(0),
  failure_direction = character(0), abs_effect = numeric(0),
  reversed = logical(0), conserved_human_mouse = logical(0), stringsAsFactors = FALSE)

# NOTE: no quit() here — run_all.R sources steps in one session; an empty result is a
# valid "data not yet downloaded" state, not an error. We end the script naturally so
# the driver proceeds (downstream steps also degrade gracefully on empty inputs).
if (is.null(dir_tabs) || !nrow(dir_tabs)) {
  .log("No failure dataset yielded a usable per-gene direction (scRNA pseudobulk pending Step6). ",
       "Emitting empty failure_specific_program.tsv with schema; re-run after Step6 pseudobulk.", level = "WARN")
  save_tsv(build_empty_fsp(), file.path(PARAMS$dir_proj, "failure_specific_program.tsv"))
} else {
  consensus <- function(v) if (all(v == "down")) "down" else if (all(v == "up")) "up" else "mixed"
  mouse_dirs <- dir_tabs[dir_tabs$organism == "mouse", ]
  human_dirs <- dir_tabs[dir_tabs$organism == "human", ]

  mk_call <- function(d, col) {
    if (!nrow(d)) return(stats::setNames(data.frame(gene = character(0), x = character(0)), c("gene", col)))
    agg <- tapply(d$dir_failure, d$gene, consensus)
    stats::setNames(data.frame(gene = names(agg), x = as.character(agg), stringsAsFactors = FALSE), c("gene", col))
  }
  mouse_call <- mk_call(mouse_dirs, "dir_failure_mouse")
  human_call <- mk_call(human_dirs, "dir_failure_human")

  # magnitude: prefer mouse-primary mean |delta|, fall back to human mean |delta|
  abs_mouse <- if (nrow(mouse_dirs)) tapply(abs(mouse_dirs$delta), mouse_dirs$gene, mean, na.rm = TRUE) else numeric(0)
  abs_human <- if (nrow(human_dirs)) tapply(abs(human_dirs$delta), human_dirs$gene, mean, na.rm = TRUE) else numeric(0)

  all_genes <- unique(c(mouse_call$gene, human_call$gene))
  fsp <- data.frame(
    gene        = all_genes,
    dir_success = unname(succ_dir[all_genes]),
    dir_failure = mouse_call$dir_failure_mouse[match(all_genes, mouse_call$gene)],
    stringsAsFactors = FALSE
  )
  fsp$dir_failure[is.na(fsp$dir_failure)] <-
    human_call$dir_failure_human[match(fsp$gene[is.na(fsp$dir_failure)], human_call$gene)]
  fsp$failure_direction <- fsp$dir_failure                       # alias consumed by Step10
  am <- abs_mouse[fsp$gene]; ah <- abs_human[fsp$gene]
  fsp$abs_effect <- ifelse(!is.na(am), am, ah); fsp$abs_effect[is.na(fsp$abs_effect)] <- 0

  fsp$reversed <- !is.na(fsp$dir_success) & !is.na(fsp$dir_failure) &
                  fsp$dir_success == "up" & fsp$dir_failure == "down"
  hcall <- human_call$dir_failure_human[match(fsp$gene, human_call$gene)]
  mcall <- mouse_call$dir_failure_mouse[match(fsp$gene, mouse_call$gene)]
  fsp$conserved_human_mouse <- fsp$dir_success == "up" &
    !is.na(mcall) & mcall == "down" & !is.na(hcall) & hcall == "down"
  fsp$conserved_human_mouse[is.na(fsp$conserved_human_mouse)] <- FALSE

  fsp <- fsp[order(-fsp$reversed, -fsp$conserved_human_mouse, -fsp$abs_effect, fsp$gene), ]
  save_tsv(fsp, file.path(PARAMS$dir_proj, "failure_specific_program.tsv"))

  cat("\n=============== Step5 failure-side projection ===============\n")
  cat(sprintf("core genes with a failure direction: %d\n", nrow(fsp)))
  cat(sprintf("  reversed (success-up / failure-down): %d\n", sum(fsp$reversed)))
  cat(sprintf("  conserved human+mouse reversal:       %d\n", sum(fsp$conserved_human_mouse)))
  cat("\nHYPOTHESIS-GENERATING: human failure data are thin/confounded (GSE7116 no healthy ctrl;\n",
      "GSE303003 n=2). The mouse mandibular axis is the primary signal; treat the reversed/\n",
      "conserved lists as candidate failure-specific nodes, not validated effects.\n", sep = "")
}

dump_session()
.log("Step5 done. Failure-specific program written to ", PARAMS$dir_proj)
