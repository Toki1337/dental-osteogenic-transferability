# run_all.R — top-level driver for the maxillary-regeneration dry-lab pipeline.
# Sequentially sources Step0 -> Step11. Step0 produces an auto-proposal that a
# HUMAN must curate (sample_arms_final.tsv + rra_usable_datasets.tsv) before the
# quantitative stages can run; an explicit guard enforces that hand-off.
#
# Usage:
#   Rscript run_all.R                 # run the whole pipeline (stops after Step0 if curation missing)
#   Rscript run_all.R --from=2        # resume from Step2 (e.g. after human curation)
#   Rscript run_all.R --from=07       # resume from Step7 (MR + coloc)
#
# Design note: this is a thin orchestrator. ALL thresholds live in config/params.R,
# the dataset registry in config/datasets.tsv. Nothing here computes results;
# every number is produced by the Step scripts on the user's downloaded data.

source("R/utils.R")

## ---- step table (script-number ordering = execution order) -----------------
STEPS <- data.frame(
  num    = sprintf("%02d", 0:11),
  script = c(
    "R/step00_curation.R",            # 0  口径治理 (human-in-the-loop)
    "R/step01_per_dataset_DE.R",      # 1  per-dataset DE -> ranked lists
    "R/step02_rra_meta.R",            # 2  pan-dental RRA meta-signature (+LOO)
    "R/step03_wgcna_modules.R",       # 3  WGCNA module + PPI hub -> core program
    "R/step04_success_projection.R",  # 4  success-side cross-context projection
    "R/step05_failure_projection.R",  # 5  failure-specific program (mouse->human)
    "R/step06_scrna_cellchat.R",      # 6  scRNA success-vs-failure (hypothesis-generating)
    "R/step07_mr_coloc.R",            # 7  druggable-genome MR + coloc + SMR-HEIDI
    "R/step08_position_context.R",    # 8  jaw position-specific context
    "R/step09_cmap_repurposing.R",    # 9  CMap coda (MR-intersected double-evidence)
    "R/step10_integration.R",         # 10 weighted multi-evidence priority score
    "R/step11_limitations.R"          # 11 limitations + assembled report tables
  ),
  label = c(
    "comparison-arm governance", "per-dataset DE", "RRA meta-signature",
    "WGCNA core program", "success-side projection", "failure-specific program",
    "scRNA success-vs-failure", "druggable-genome MR + coloc", "jaw position context",
    "CMap repurposing coda", "integration priority score", "report + limitations"
  ),
  stringsAsFactors = FALSE
)

## ---- parse --from=NN --------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
from_arg <- sub("^--from=", "", grep("^--from=", args, value = TRUE))
from_idx <- if (length(from_arg)) {
  i <- match(sprintf("%02d", as.integer(from_arg)), STEPS$num)
  if (is.na(i)) stop("--from must be an integer 0..11; got: ", from_arg)
  i
} else 1L

## ---- human-curation guard (THE make-or-break hand-off) ----------------------
# Step0 only auto-proposes group labels; a human must confirm them. Steps 1+ are
# blocked until BOTH curated files exist. We re-source params here to read paths.
source("config/params.R", local = TRUE)
curation_ready <- function() {
  all(file.exists(
    file.path(PARAMS$dir_curation, "sample_arms_final.tsv"),
    file.path(PARAMS$dir_curation, "rra_usable_datasets.tsv")
  ))
}

## ---- run -------------------------------------------------------------------
.log("================ run_all: maxillary regeneration pipeline ================")
.log("executing steps ", STEPS$num[from_idx], " .. 11",
     if (length(from_arg)) sprintf("  (--from=%s)", from_arg) else "")

t_all <- Sys.time()
for (i in seq(from_idx, nrow(STEPS))) {
  s <- STEPS[i, ]

  # Hard stop between Step0 and the quantitative stages until human curation lands.
  if (as.integer(s$num) >= 1L && !curation_ready()) {
    stop(paste0(
      "\n=================== ACTION REQUIRED — pipeline paused ===================\n",
      "Step0 produced an AUTO-PROPOSAL only. Before continuing you MUST hand-curate:\n",
      "  1) ", file.path(PARAMS$dir_curation, "sample_arms_final.tsv"), "\n",
      "       (fill final_group = osteo|control|exclude and arm_clean = yes/no per sample;\n",
      "        resolve every 'review' row; drop药物×成骨交互的非纯对照臂)\n",
      "  2) ", file.path(PARAMS$dir_curation, "rra_usable_datasets.tsv"), "\n",
      "       (produced by validate_curation() in Step0 after you edit (1))\n",
      "Then resume with:  Rscript run_all.R --from=", s$num, "\n",
      "Rationale: comparison-arm governance is THE make-or-break step; it cannot be\n",
      "automated without inflating the RRA training pool with impure contrasts.\n",
      "========================================================================\n"
    ), call. = FALSE)
  }

  if (!file.exists(s$script)) stop("missing step script: ", s$script, call. = FALSE)

  .log(sprintf("---- Step%s START : %s  (%s) ----", s$num, s$label, s$script))
  t0 <- Sys.time()
  ok <- tryCatch({ source(s$script, local = new.env()); TRUE },
                 error = function(e) { .log("Step", s$num, " FAILED: ", conditionMessage(e), level = "ERROR"); FALSE })
  if (!ok) stop("aborting at Step", s$num, " — fix the error above, then resume with --from=", s$num, call. = FALSE)
  dt <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  .log(sprintf("---- Step%s DONE  : %s  (%ss) ----", s$num, s$label, dt))
}
dt_all <- round(as.numeric(difftime(Sys.time(), t_all, units = "mins")), 2)
.log("================ pipeline complete in ", dt_all, " min ================")
