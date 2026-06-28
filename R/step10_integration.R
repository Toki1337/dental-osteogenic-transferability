# R/step10_integration.R
# Step 10 — Multi-evidence weighted target prioritization (transparent, SHAP-style).
# ----------------------------------------------------------------------------
# For each candidate target gene we assemble seven evidence dimensions, min-max
# normalize each to [0,1], then take a weighted sum using PARAMS$weights. Because
# the score is a linear combination of normalized features, the per-dimension
# WEIGHTED VALUE is exactly that dimension's additive contribution to the total
# (an exact additive decomposition; no surrogate model needed — we report it as a
# SHAP-style contribution table for transparency, not a fitted SHAP run).
#
# Evidence dimensions (-> PARAMS$weights keys):
#   rra_strength      core_program RRA score, normalized (robustness in pan-dental program)
#   failure_specific  failure-specific dysregulation magnitude (Step5 projection)
#   scrna_support     single-cell corroboration of dysregulation (Step6)
#   mr_causal         MR + coloc + HEIDI + FinnGen passing (Step7) — highest weight
#   druggability      Finan druggable-genome tier (1>2>3)
#   jaw_context       jaw position-specific direction-consistent osteo weakness (GSE58474/GSE30167/DRA006607)
#   cmap_connectivity double-evidence (MR INTERSECT CMap) flag (Step9)
#
# Outputs:
#   09_integration/priority_targets.tsv         (full: per-dim score, weighted contribution, total,
#                                                protective direction, matched double-evidence drug,
#                                                blank annotation columns to fill: local delivery /
#                                                safety / clinical pipeline)
#   09_integration/candidate_table_for_paper.tsv (trimmed "pro-jaw-healing target / old-drug" table)
#
# Honesty guard: failure-side human data are thin (GSE7116 MM/N-BP confound; small
# scRNA). Targets resting on thin evidence are flagged hypothesis-generating in the
# evidence_tier column; we do not over-claim MRONJ causality.

source("R/utils.R"); init_pipeline()
set.seed(PARAMS$seed)

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
read_if <- function(path, required = FALSE) {
  if (file.exists(path)) return(utils::read.delim(path, stringsAsFactors = FALSE))
  if (required) stop("Missing required input: ", path)
  .log("optional input absent: ", path, level = "WARN"); NULL
}
# min-max to [0,1]; NA-safe; constant/empty vector -> all 0 (no discriminative signal)
mm01 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) return(rep(0, length(x)))
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(rep(0, length(x)))
  out <- (x - rng[1]) / diff(rng)
  out[is.na(out)] <- 0
  out
}

# ---------------------------------------------------------------------------
# 1. candidate universe = core program genes (anchor on the competent program)
# ---------------------------------------------------------------------------
core <- read_if(file.path(PARAMS$dir_wgcna, "core_program.tsv"), required = TRUE)
stopifnot(all(c("gene", "direction") %in% colnames(core)))
genes <- unique(core$gene)
tgt <- data.frame(gene = genes, stringsAsFactors = FALSE)
tgt$direction <- core$direction[match(tgt$gene, core$gene)]

# ---------------------------------------------------------------------------
# 2. dimension: rra_strength (from core_program rra_score; smaller RRA = stronger)
# ---------------------------------------------------------------------------
if ("rra_score" %in% colnames(core)) {
  rs <- core$rra_score[match(tgt$gene, core$gene)]
  tgt$rra_strength <- mm01(-rs)            # invert: smaller score -> larger strength
} else tgt$rra_strength <- 0

# ---------------------------------------------------------------------------
# 3. dimension: failure_specific (Step5 projection of failure-specific program)
# ---------------------------------------------------------------------------
# Expected schema 05_projection/failure_specific_program.tsv:
#   gene, failure_direction, abs_effect (|ssGSEA/logFC magnitude|), conserved_human_mouse(logical)
fs <- read_if(file.path(PARAMS$dir_proj, "failure_specific_program.tsv"))
tgt$failure_abs_effect    <- if (!is.null(fs)) (fs$abs_effect %||% NA)[match(tgt$gene, fs$gene)] else NA
tgt$failure_conserved     <- if (!is.null(fs)) (fs$conserved_human_mouse %||% NA)[match(tgt$gene, fs$gene)] else NA
tgt$failure_specific      <- mm01(tgt$failure_abs_effect)

# ---------------------------------------------------------------------------
# 4. dimension: scrna_support (Step6 single-cell dysregulation corroboration)
# ---------------------------------------------------------------------------
# Expected schema 06_scrna/scrna_dysregulated_genes.tsv:
#   gene, cell_type, is_dysregulated(logical), abs_avg_log2FC, dataset
sc <- read_if(file.path(PARAMS$dir_scrna, "scrna_dysregulated_genes.tsv"))
if (!is.null(sc)) {
  scl <- sc$abs_avg_log2FC %||% rep(1, nrow(sc))
  dys <- if ("is_dysregulated" %in% colnames(sc)) tolower(sc$is_dysregulated) %in% c("true","yes","1") else TRUE
  agg <- tapply(scl[dys], sc$gene[dys], function(v) max(v, na.rm = TRUE))
  tgt$scrna_abs_log2FC <- as.numeric(agg[tgt$gene])
} else tgt$scrna_abs_log2FC <- NA
tgt$scrna_support <- mm01(tgt$scrna_abs_log2FC)

# ---------------------------------------------------------------------------
# 5. dimension: mr_causal (Step7 MR + coloc + HEIDI + FinnGen)
# ---------------------------------------------------------------------------
# Expected schema 07_mr_coloc/mr_protective_druggable_targets.tsv:
#   gene, protective_direction, mr_method, b, se, pval, padj, coloc_pp4, heidi_p,
#   finngen_replicated(logical), finan_tier, passes(logical)
mr <- read_if(file.path(PARAMS$dir_mr, "mr_protective_druggable_targets.tsv"))
tgt$protective_direction <- NA_character_
tgt$mr_padj <- NA_real_; tgt$coloc_pp4 <- NA_real_; tgt$heidi_p <- NA_real_
tgt$finngen_replicated <- NA; tgt$finan_tier <- NA
if (!is.null(mr)) {
  mi <- match(tgt$gene, mr$gene)
  tgt$protective_direction <- (mr$protective_direction %||% NA)[mi]
  tgt$mr_padj    <- (mr$padj %||% mr$pval %||% NA)[mi]
  tgt$coloc_pp4  <- (mr$coloc_pp4 %||% NA)[mi]
  tgt$heidi_p    <- (mr$heidi_p %||% NA)[mi]
  tgt$finngen_replicated <- (mr$finngen_replicated %||% NA)[mi]
  tgt$finan_tier <- (mr$finan_tier %||% NA)[mi]
  passes <- if ("passes" %in% colnames(mr)) (tolower(mr$passes) %in% c("true","yes","1"))[mi] else !is.na(mi)
  # graded causal score: pass gate, then strengthen by coloc PP4 and -log10(padj)
  passes[is.na(passes)] <- FALSE
  pp4    <- tgt$coloc_pp4; pp4[is.na(pp4)] <- 0
  causal <- -log10(pmax(tgt$mr_padj, .Machine$double.xmin)); causal[is.na(causal)] <- 0
  raw_mr <- ifelse(passes, 0.5 + 0.5 * mm01(causal + pp4), 0)
  tgt$mr_causal <- raw_mr
} else tgt$mr_causal <- 0

# ---------------------------------------------------------------------------
# 6. dimension: druggability (Finan tier: 1 > 2 > 3 > none)
# ---------------------------------------------------------------------------
tier_score <- function(t) {
  ti <- suppressWarnings(as.integer(gsub("[^0-9]", "", as.character(t))))
  s <- c(`1` = 1, `2` = 0.66, `3` = 0.33)[as.character(ti)]
  s <- unname(s)
  if (length(s) != 1 || is.na(s)) 0 else s
}
tgt$druggability <- vapply(tgt$finan_tier, tier_score, numeric(1))

# ---------------------------------------------------------------------------
# 7. dimension: jaw_context (Step8 / position-specific osteo weakness alignment)
# ---------------------------------------------------------------------------
# Expected schema 05_projection/jaw_position_context.tsv (from GSE58474/GSE30167/DRA006607):
#   gene, jaw_weak_direction_consistent(logical), abs_jaw_effect
jc <- read_if(file.path(PARAMS$dir_proj, "jaw_position_context.tsv"))
if (!is.null(jc)) {
  consistent <- if ("jaw_weak_direction_consistent" %in% colnames(jc))
    tolower(jc$jaw_weak_direction_consistent) %in% c("true","yes","1") else TRUE
  eff <- jc$abs_jaw_effect %||% rep(1, nrow(jc))
  eff[!consistent] <- 0
  tgt$jaw_abs_effect <- as.numeric(tapply(eff, jc$gene, max)[tgt$gene])
  tgt$jaw_direction_consistent <- as.logical(tapply(consistent, jc$gene, any)[tgt$gene])
} else { tgt$jaw_abs_effect <- NA; tgt$jaw_direction_consistent <- NA }
tgt$jaw_context <- mm01(tgt$jaw_abs_effect)

# ---------------------------------------------------------------------------
# 8. dimension: cmap_connectivity (Step9 double-evidence flag + matched drug)
# ---------------------------------------------------------------------------
de <- read_if(file.path(PARAMS$dir_cmap, "double_evidence_candidates.tsv"))
tgt$double_evidence_drug <- NA_character_; tgt$cmap_connectivity <- 0
if (!is.null(de) && nrow(de) && "target_gene" %in% colnames(de)) {
  di <- match(tgt$gene, de$target_gene)
  tgt$double_evidence_drug <- (de$drug %||% NA)[di]
  # binary connectivity signal: presence in the double-evidence set => 1, else 0
  tgt$cmap_connectivity <- as.numeric(!is.na(di))
}

# ---------------------------------------------------------------------------
# 9. weighted aggregation + SHAP-style additive decomposition
# ---------------------------------------------------------------------------
w <- PARAMS$weights
dims <- names(w)
stopifnot(all(dims %in% colnames(tgt)))

# per-dimension WEIGHTED contribution = weight * normalized_score (exact additive)
for (d in dims) tgt[[paste0("contrib_", d)]] <- w[[d]] * tgt[[d]]
tgt$priority_score <- rowSums(tgt[, paste0("contrib_", dims), drop = FALSE])

# evidence tier (honesty): downgrade targets leaning on thin failure-side data
tgt$evidence_tier <- ifelse(tgt$mr_causal > 0 & tgt$cmap_connectivity > 0, "double_evidence",
                   ifelse(tgt$mr_causal > 0, "genetic_causal",
                   ifelse(tgt$failure_specific > 0 | tgt$scrna_support > 0,
                          "transcriptomic_hypothesis_generating", "program_member")))

tgt <- tgt[order(-tgt$priority_score), , drop = FALSE]

# blank annotation columns to fill manually (kept explicit & empty)
tgt$local_delivery_feasibility <- ""
tgt$safety_note               <- ""
tgt$clinical_pipeline         <- ""

priority_cols <- c(
  "gene", "direction", "protective_direction", "priority_score", "evidence_tier",
  dims,
  paste0("contrib_", dims),
  "mr_padj", "coloc_pp4", "heidi_p", "finngen_replicated", "finan_tier",
  "failure_abs_effect", "failure_conserved", "scrna_abs_log2FC",
  "jaw_abs_effect", "jaw_direction_consistent",
  "double_evidence_drug",
  "local_delivery_feasibility", "safety_note", "clinical_pipeline")
priority_cols <- intersect(priority_cols, colnames(tgt))
save_tsv(tgt[, priority_cols, drop = FALSE], file.path(PARAMS$dir_integr, "priority_targets.tsv"))

# ---------------------------------------------------------------------------
# 10. trimmed paper-facing candidate table
# ---------------------------------------------------------------------------
paper <- tgt[, c("gene", "direction", "protective_direction", "priority_score",
               "evidence_tier", "mr_causal", "druggability", "cmap_connectivity",
               "double_evidence_drug", "finan_tier",
               "local_delivery_feasibility", "safety_note", "clinical_pipeline")]
colnames(paper) <- c("target_gene", "program_direction", "protective_direction",
                     "priority_score", "evidence_tier", "mr_causal_score",
                     "druggability_score", "cmap_double_evidence", "candidate_old_drug",
                     "finan_tier", "local_delivery_feasibility", "safety_note", "clinical_pipeline")
# paper table: lead with the strongest, most-defensible candidates
paper <- utils::head(paper[order(-paper$priority_score), , drop = FALSE], 50)
save_tsv(paper, file.path(PARAMS$dir_integr, "candidate_table_for_paper.tsv"))

cat(sprintf("\nStep10 integration done: %d candidate targets scored.\n", nrow(tgt)))
cat(sprintf("  double_evidence: %d | genetic_causal: %d | transcriptomic(hyp-gen): %d\n",
            sum(tgt$evidence_tier == "double_evidence"),
            sum(tgt$evidence_tier == "genetic_causal"),
            sum(tgt$evidence_tier == "transcriptomic_hypothesis_generating")))
cat("  weights: ", paste(sprintf("%s=%.2f", dims, unlist(w[dims])), collapse = ", "), "\n", sep = "")
dump_session()
