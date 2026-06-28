# R/step00_curation.R
# Step 0 — Comparison-arm governance (THE make-or-break step).
# Goal: from each candidate dataset, retain ONLY clean "osteogenic-induction vs
# non-induced/control" arms. Auto-propose a group label per sample from GEO
# metadata, then REQUIRE a human to confirm/override in final_group.
#
# Output:
#   01_curation/sample_arms_template.tsv   (auto-proposal; human fills final_group)
#   01_curation/sample_arms_final.tsv      (human-edited; consumed by Step1)
# A dataset enters the RRA training pool only if it has >=1 'osteo' and >=1 'control'
# clean arm in sample_arms_final.tsv.

source("R/utils.R"); init_pipeline()

osteo_kw   <- "osteogen|mineraliz|alizarin|OM |OIM|osteo[- ]?induc|differentiat.*osteo"
control_kw <- "control|ctrl|undifferentiat|growth medium|GM\\b|basal|non[- ]?induc|day ?0|uninduc|vehicle|WT|wild[- ]?type"
exclude_kw <- "odontogen(?!.*osteo)|adipogen|chondrogen|ETV2|reprogram|knockout|overexpress|siRNA|sh[A-Z]"

guess_group <- function(text) {
  t <- tolower(text)
  if (grepl(exclude_kw, t, perl = TRUE)) return("exclude")
  has_o <- grepl(osteo_kw, t, perl = TRUE)
  has_c <- grepl(control_kw, t, perl = TRUE)
  if (has_o && !has_c) return("osteo")
  if (has_c && !has_o) return("control")
  if (has_o && has_c)  return("review")   # both keywords -> human decides
  "review"
}

reg <- REG[REG$include_in_rra | REG$is_external_validation, , drop = FALSE]
rows <- list()

for (i in seq_len(nrow(reg))) {
  acc <- reg$accession[i]; role <- reg$role[i]
  ok <- tryCatch({ eset <- fetch_geo(acc, role); TRUE }, error = function(e) { .log("skip ", acc, ": ", conditionMessage(e), level = "WARN"); FALSE })
  if (!ok) next
  pd <- Biobase::pData(eset)
  # collapse all characteristics columns into one descriptive string per sample
  char_cols <- grep("characteristics|title|source_name|description|:ch1", colnames(pd), ignore.case = TRUE, value = TRUE)
  desc <- apply(pd[, char_cols, drop = FALSE], 1, function(r) paste(na.omit(as.character(r)), collapse = " | "))
  rows[[acc]] <- data.frame(
    accession   = acc,
    role        = role,
    cell_type   = reg$cell_type[i],
    organism    = reg$organism[i],
    gsm         = rownames(pd),
    title       = pd$title %||% NA,
    description = desc,
    proposed_group = vapply(desc, guess_group, character(1)),
    final_group = "",         # <-- HUMAN FILLS: osteo | control | exclude
    arm_clean   = "",         # <-- HUMAN: yes/no (is this a pure osteo-vs-ctrl contrast?)
    reviewer_note = "",
    stringsAsFactors = FALSE
  )
}

tmpl <- do.call(rbind, rows)
save_tsv(tmpl, file.path(PARAMS$dir_curation, "sample_arms_template.tsv"))

cat("\n================ Step0 summary ================\n")
print(table(tmpl$accession, tmpl$proposed_group))
cat("\nACTION REQUIRED:\n",
    "  1) Open 01_curation/sample_arms_template.tsv\n",
    "  2) Fill final_group (osteo/control/exclude) and arm_clean (yes/no) per sample\n",
    "  3) Resolve every 'review' row; drop药物×成骨交互臂的非纯对照 (e.g. GSE226347 confirm recipient cell)\n",
    "  4) Save as 01_curation/sample_arms_final.tsv\n",
    "  -> The actual number of clean RRA training arms (audit predicts ~11-12) is REPORTED from this file, not assumed.\n", sep = "")

# ---- validator (run after human edit) ----
validate_curation <- function(path = file.path(PARAMS$dir_curation, "sample_arms_final.tsv")) {
  stopifnot(file.exists(path))
  d <- utils::read.delim(path, stringsAsFactors = FALSE)
  # Purity gate: a sample only counts if a human marked arm_clean as yes/true.
  # Fail-safe: blank/NA/anything-else arm_clean is treated as NOT clean (excluded).
  d <- d[d$final_group %in% c("osteo", "control") &
           (tolower(d$arm_clean) %in% c("yes", "y", "true", "1")), , drop = FALSE]
  tab <- table(d$accession, d$final_group)
  # Guard against a missing osteo/control column (table() drops absent levels).
  oc <- if ("osteo" %in% colnames(tab)) tab[, "osteo"] else 0
  cc <- if ("control" %in% colnames(tab)) tab[, "control"] else 0
  usable <- rownames(tab)[oc >= 1 & cc >= 1]
  cat("Clean osteo-vs-control datasets usable for RRA:", length(usable), "\n")
  print(tab[usable, , drop = FALSE])
  save_tsv(data.frame(accession = usable), file.path(PARAMS$dir_curation, "rra_usable_datasets.tsv"))
  invisible(usable)
}
dump_session()
