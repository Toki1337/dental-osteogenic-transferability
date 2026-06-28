# R/step09_cmap_repurposing.R
# Step 9 — CMap / L1000 connectivity repurposing (CODA, NOT the main line).
# ----------------------------------------------------------------------------
# Rationale & framing (write this verbatim into Methods/Discussion):
#   The CMap->osteogenesis->old-drug repurposing chain is an ESTABLISHED paradigm
#   (Brey DM et al., PNAS 2015, parbendazole; PMID 26420877 / DOI 10.1073/pnas.1501597112).
#   To avoid an "incremental" critique we DELIBERATELY DOWNGRADE CMap to a
#   hypothesis-generating coda and keep ONLY candidates that intersect the
#   druggable-genome MR targets (Step7) = "double evidence" (genetic causality
#   AND pharmacological connectivity). No efficacy is claimed; all outputs are
#   computational hypotheses for a dry-lab, no-wet-validation study.
#
# Query construction:
#   The "regeneration-competent program" up/down genes (04_modules_wgcna/core_program.tsv)
#   are used as the query signature. We probe the free, open maayanlab APIs:
#     - L1000CDS2  (POST JSON, mimic AND reverse modes)
#     - L1000FWD   (signature search; reverse/forward endpoints)
#   We seek NEGATIVE connectivity (perturbations that REVERSE the failure state /
#   MIMIC the competent osteogenic program), consistent across cell lines
#   (>= PARAMS$cmap_min_celllines).
#
#   clue.io (CMap/LINCS Touchstone, tau-based) requires a registered account &
#   API key; it is NOT called here. An offline interface stub is provided below
#   (query_clue_io) documenting the exact payload so a user with an API key can
#   drop it in. Our negative-tau threshold lives in PARAMS$cmap_neg_score.
#
# Output:
#   08_cmap/l1000_raw_hits.tsv          (all returned perturbation hits, annotated)
#   08_cmap/cmap_candidates.tsv         (negative-connectivity drugs, cross-cell-line consistent)
#   08_cmap/double_evidence_candidates.tsv  (cmap_candidates  INTERSECT  MR protective druggable targets)

source("R/utils.R"); init_pipeline()
suppressMessages({ library(httr); library(jsonlite) })
set.seed(PARAMS$seed)

# ---------------------------------------------------------------------------
# 0. inputs (defensive: this is a coda; degrade gracefully if upstream missing)
# ---------------------------------------------------------------------------
core_path <- file.path(PARAMS$dir_wgcna, "core_program.tsv")
mr_path   <- file.path(PARAMS$dir_mr, "mr_protective_druggable_targets.tsv")

if (!file.exists(core_path)) stop("Missing ", core_path, " — run Step3 first.")
core <- utils::read.delim(core_path, stringsAsFactors = FALSE)
stopifnot(all(c("gene", "direction") %in% colnames(core)))

# L1000 platform measures ~978 landmark genes; oversized lists are silently
# truncated server-side. Cap query to the strongest nodes (by RRA score) if available.
cap <- 150L
trim <- function(dir) {
  g <- core[core$direction == dir, , drop = FALSE]
  if ("rra_score" %in% colnames(g)) g <- g[order(g$rra_score), ]  # RRA score: smaller = stronger
  utils::head(unique(g$gene), cap)
}
up_genes   <- trim("up")
down_genes <- trim("down")
.log("CMap query signature: ", length(up_genes), " up / ", length(down_genes), " down (core program)")

# ---------------------------------------------------------------------------
# 1. L1000CDS2 — open API (mimic + reverse). https://maayanlab.cloud/L1000CDS2/query
# ---------------------------------------------------------------------------
# config = list(aggravate = FALSE) => REVERSE mode (perturbations that reverse the
#   input signature); aggravate = TRUE => MIMIC mode. We run both and label.
query_l1000cds2 <- function(up, down, aggravate, max_hits = 50L) {
  url <- "https://maayanlab.cloud/L1000CDS2/query"
  payload <- list(
    data = list(upGenes = as.list(up), dnGenes = as.list(down)),
    config = list(aggravate = aggravate, searchMethod = "geneSet",
                  share = FALSE, combination = FALSE, "db-version" = "latest"),
    metadata = list()
  )
  resp <- tryCatch(
    POST(url, body = toJSON(payload, auto_unbox = TRUE), encode = "raw",
         content_type("application/json"), timeout(120)),
    error = function(e) { .log("L1000CDS2 POST failed: ", conditionMessage(e), level = "WARN"); NULL }
  )
  if (is.null(resp) || http_error(resp)) {
    .log("L1000CDS2 returned no usable response (aggravate=", aggravate, ")", level = "WARN")
    return(NULL)
  }
  js <- tryCatch(fromJSON(content(resp, "text", encoding = "UTF-8"), simplifyVector = FALSE),
                 error = function(e) NULL)
  topMeta <- js$topMeta
  if (is.null(topMeta) || !length(topMeta)) return(NULL)
  rows <- lapply(topMeta, function(h) data.frame(
    pert_desc   = (h$pert_desc %||% h$pert_id) %||% NA_character_,
    pert_id     = h$pert_id %||% NA_character_,
    cell_line   = h$cell_id %||% NA_character_,
    dose        = h$pert_dose %||% NA_character_,
    time        = h$pert_time %||% NA_character_,
    score       = as.numeric(h$score %||% NA),     # CDS2 cosine-based; lower=closer reverse
    mode        = if (aggravate) "mimic" else "reverse",
    source_api  = "L1000CDS2",
    stringsAsFactors = FALSE))
  out <- do.call(rbind, rows)
  utils::head(out[order(out$score), , drop = FALSE], max_hits)
}

# ---------------------------------------------------------------------------
# 2. L1000FWD — open API. https://maayanlab.cloud/L1000FWD/
# ---------------------------------------------------------------------------
# Returns 'similar' (mimic) and 'opposite' (reverse) signature ids; we keep opposite.
query_l1000fwd <- function(up, down, max_hits = 50L) {
  base <- "https://maayanlab.cloud/L1000FWD"
  payload <- list(up_genes = as.list(up), down_genes = as.list(down))
  resp <- tryCatch(
    POST(file.path(base, "sig_search"),
         body = toJSON(payload, auto_unbox = TRUE), encode = "raw",
         content_type("application/json"), timeout(120)),
    error = function(e) { .log("L1000FWD sig_search failed: ", conditionMessage(e), level = "WARN"); NULL }
  )
  if (is.null(resp) || http_error(resp)) return(NULL)
  rid <- tryCatch(fromJSON(content(resp, "text", encoding = "UTF-8"))$result_id,
                  error = function(e) NULL)
  if (is.null(rid)) return(NULL)
  res <- tryCatch(
    GET(file.path(base, "result/topn", rid), timeout(120)),
    error = function(e) NULL)
  if (is.null(res) || http_error(res)) return(NULL)
  js <- fromJSON(content(res, "text", encoding = "UTF-8"), simplifyDataFrame = FALSE)
  opp <- js$opposite  # signatures whose effect opposes input (= reverse failure / restore program)
  if (is.null(opp) || !length(opp)) return(NULL)
  meta <- function(s, mode) data.frame(
    pert_desc  = s$pert_desc %||% s$pert_id %||% NA_character_,
    pert_id    = s$pert_id %||% NA_character_,
    cell_line  = s$cell_id %||% NA_character_,
    dose       = NA_character_, time = NA_character_,
    score      = as.numeric(s$scores %||% s$combined_score %||% NA),
    mode       = mode, source_api = "L1000FWD", stringsAsFactors = FALSE)
  out <- do.call(rbind, lapply(opp, meta, mode = "reverse"))
  utils::head(out[order(out$score), , drop = FALSE], max_hits)
}

# ---------------------------------------------------------------------------
# 3. clue.io — OFFLINE STUB (account/API key required; NOT called here).
# ---------------------------------------------------------------------------
# To enable: register at https://clue.io, obtain a user_key, then POST a gene-set
# query to the Query API. Negative connectivity = normalized tau <= PARAMS$cmap_neg_score.
# The returned per-cell-line tau lets you apply the same cross-cell-line consistency
# filter as below. We expose the contract so a key-holder can wire it in unchanged.
query_clue_io <- function(up, down, user_key = Sys.getenv("CLUE_API_KEY")) {
  if (!nzchar(user_key)) {
    .log("clue.io skipped: no CLUE_API_KEY (account-gated; offline by design)", level = "WARN")
    return(NULL)
  }
  # Reference payload (Touchstone gene-set query); enable only with a valid key.
  url <- "https://api.clue.io/api/jobs"
  body <- list(
    name = "maxillary_regen_competent_query",
    uptag  = as.list(up),
    dntag  = as.list(down),
    data_type = "L1000", dataset = "Touchstone"
  )
  resp <- tryCatch(POST(url, add_headers(user_key = user_key),
                        body = toJSON(body, auto_unbox = TRUE), encode = "raw",
                        content_type("application/json"), timeout(180)),
                   error = function(e) { .log("clue.io POST failed: ", conditionMessage(e), level = "WARN"); NULL })
  if (is.null(resp) || http_error(resp)) return(NULL)
  # NOTE: clue jobs are asynchronous; poll job status then fetch ncs/tau matrix.
  # Left as a documented stub; cross-cell-line tau <= PARAMS$cmap_neg_score => keep.
  .log("clue.io job submitted; async result retrieval not automated in this coda", level = "WARN")
  NULL
}

# ---------------------------------------------------------------------------
# 4. run queries, pool raw hits
# ---------------------------------------------------------------------------
raw_list <- list(
  query_l1000cds2(up_genes, down_genes, aggravate = FALSE),  # reverse
  query_l1000cds2(up_genes, down_genes, aggravate = TRUE),   # mimic
  query_l1000fwd(up_genes, down_genes),
  query_clue_io(up_genes, down_genes)
)
raw_list <- Filter(Negate(is.null), raw_list)

# ---------------------------------------------------------------------------
# Drug -> target gene map (for the small-molecule double-evidence intersection).
# Loaded once here so both the empty branch and main body share the same schema.
# ---------------------------------------------------------------------------
DRUG_TARGET_MAP_TSV <- "config/drug_target_map.tsv"
DRUG_TARGET_FALLBACK <- data.frame(
  drug = c("SIMVASTATIN", "ATORVASTATIN", "ALENDRONATE", "ZOLEDRONATE",
           "RISEDRONATE", "DENOSUMAB", "ROMOSOZUMAB", "TERIPARATIDE",
           "ODANACATIB", "PARBENDAZOLE"),
  target_gene = c("HMGCR", "HMGCR", "FDPS", "FDPS", "FDPS", "TNFSF11",
                  "SOST", "PTH1R", "CTSK", "TUBB"),
  stringsAsFactors = FALSE)
load_drug_target_map <- function(path = DRUG_TARGET_MAP_TSV) {
  if (!file.exists(path)) {
    .log("drug-target map not found at ", path, level = "WARN")
    .log("TODO(user): supply a curated drug->target table (DGIdb / ChEMBL mechanism / ",
         "DrugBank) at ", path, " with columns drug<TAB>target_gene. ",
         "Using built-in ", nrow(DRUG_TARGET_FALLBACK), "-pair EXAMPLE map (NOT complete).",
         level = "WARN")
    return(DRUG_TARGET_FALLBACK)
  }
  d <- utils::read.delim(path, stringsAsFactors = FALSE, comment.char = "#", check.names = FALSE)
  if (!all(c("drug", "target_gene") %in% colnames(d))) {
    .log("drug-target map ", path, " missing required columns drug/target_gene; ",
         "falling back to built-in EXAMPLE map.", level = "WARN")
    return(DRUG_TARGET_FALLBACK)
  }
  d <- d[!is.na(d$drug) & d$drug != "" & !is.na(d$target_gene) & d$target_gene != "",
         c("drug", "target_gene"), drop = FALSE]
  d$drug        <- toupper(trimws(d$drug))         # match against uppercased perturbagen names
  d$target_gene <- toupper(trimws(d$target_gene))
  if (nrow(d) < 50)
    .log("drug-target map has only ", nrow(d), " pairs — looks like the placeholder. ",
         "TODO(user): replace with a full DGIdb/ChEMBL/DrugBank export before final results.",
         level = "WARN")
  unique(d)
}

if (length(raw_list)) {

raw <- do.call(rbind, raw_list)
raw$drug <- toupper(trimws(as.character(raw$pert_desc)))
raw <- raw[!is.na(raw$drug) & raw$drug != "" & raw$drug != "NA", , drop = FALSE]
save_tsv(raw, file.path(PARAMS$dir_cmap, "l1000_raw_hits.tsv"))

# ---------------------------------------------------------------------------
# 5. negative-connectivity candidates, cross-cell-line consistent
# ---------------------------------------------------------------------------
# Keep 'reverse'-mode hits (reverse failure / restore competent program). Require
# the SAME drug to appear in >= PARAMS$cmap_min_celllines distinct cell lines.
neg <- raw[raw$mode == "reverse", , drop = FALSE]
agg <- by(neg, neg$drug, function(d) data.frame(
  drug         = d$drug[1],
  n_celllines  = length(unique(d$cell_line[!is.na(d$cell_line)])),
  n_sources    = length(unique(d$source_api)),
  mean_score   = mean(d$score, na.rm = TRUE),
  best_mode    = "reverse",
  cell_lines   = paste(sort(unique(stats::na.omit(d$cell_line))), collapse = ";"),
  apis         = paste(sort(unique(d$source_api)), collapse = ";"),
  stringsAsFactors = FALSE))
cand <- do.call(rbind, agg)
cand <- cand[cand$n_celllines >= PARAMS$cmap_min_celllines, , drop = FALSE]
cand$connectivity <- "negative_reverse"     # reverses failure state / mimics competent program
cand <- cand[order(cand$mean_score), , drop = FALSE]
save_tsv(cand, file.path(PARAMS$dir_cmap, "cmap_candidates.tsv"))
.log("CMap candidates (>= ", PARAMS$cmap_min_celllines, " cell lines, negative connectivity): ", nrow(cand))

# ---------------------------------------------------------------------------
# 6. DOUBLE EVIDENCE = CMap candidates  INTERSECT  MR protective druggable targets
# ---------------------------------------------------------------------------
# This is the only output we let drive the narrative: genetic causality (MR/coloc)
# AND pharmacological connectivity. Two complementary match paths (match_type):
#   drug_target          — the CMap perturbagen is a SMALL-MOLECULE drug whose
#                          annotated target gene (via config/drug_target_map.tsv)
#                          IS an MR protective druggable target. This is the path
#                          that lets compounds (e.g. simvastatin->HMGCR) intersect;
#                          the old "drug name == gene symbol" test could never match
#                          a small molecule and was structurally empty.
#   genetic_perturbation — the perturbagen name IS a gene symbol (CMap genetic
#                          OE/KD/cDNA perturbation) that equals an MR target. Kept
#                          as a complementary supplement to the drug_target path.
# Expected schema of mr_protective_druggable_targets.tsv (from Step7):
#   gene, protective_direction, mr_method, b, se, pval, padj, coloc_pp4, heidi_p,
#   finngen_replicated(logical), finan_tier, passes(logical)
double <- data.frame(drug = character(0), target_gene = character(0),
                     match_type = character(0),
                     n_celllines = integer(0), mean_score = numeric(0),
                     protective_direction = character(0), mr_pval = numeric(0),
                     coloc_pp4 = numeric(0), finan_tier = character(0),
                     stringsAsFactors = FALSE)

if (file.exists(mr_path)) {
  mr <- utils::read.delim(mr_path, stringsAsFactors = FALSE)
  mr_genes <- toupper(mr$gene)
  dt_map   <- load_drug_target_map()   # drug (UPPER) -> target_gene (UPPER)

  # helper: build a double-evidence frame given matched (drug, target_gene) pairs
  # carried alongside the originating CMap candidate row + the MR row index.
  assemble <- function(drug, target_gene, match_type, cand_idx, mr_idx) {
    if (!length(drug)) return(double[0, , drop = FALSE])
    data.frame(
      drug                 = drug,
      target_gene          = target_gene,
      match_type           = match_type,
      n_celllines          = cand$n_celllines[cand_idx],
      mean_score           = cand$mean_score[cand_idx],
      protective_direction = (mr$protective_direction %||% rep(NA, nrow(mr)))[mr_idx],
      mr_pval              = (mr$pval %||% rep(NA, nrow(mr)))[mr_idx],
      coloc_pp4            = (mr$coloc_pp4 %||% rep(NA, nrow(mr)))[mr_idx],
      finan_tier           = (mr$finan_tier %||% rep(NA, nrow(mr)))[mr_idx],
      stringsAsFactors = FALSE)
  }

  # ---- path A: small-molecule drug -> mapped target_gene -> MR protective target
  # join CMap candidates to the drug-target map, then keep pairs whose target is MR.
  jt <- merge(
    data.frame(drug = cand$drug, cand_idx = seq_len(nrow(cand)), stringsAsFactors = FALSE),
    dt_map, by = "drug", all = FALSE)
  jt <- jt[jt$target_gene %in% mr_genes, , drop = FALSE]
  dt_double <- if (nrow(jt))
    assemble(jt$drug, jt$target_gene, "drug_target",
             jt$cand_idx, match(jt$target_gene, mr_genes)) else double[0, , drop = FALSE]

  # ---- path B: perturbagen name IS a gene symbol that equals an MR target
  gp_idx <- which(cand$drug %in% mr_genes)
  gp_double <- if (length(gp_idx))
    assemble(cand$drug[gp_idx], mr$gene[match(cand$drug[gp_idx], mr_genes)],
             "genetic_perturbation", gp_idx, match(cand$drug[gp_idx], mr_genes)) else double[0, , drop = FALSE]

  double <- rbind(dt_double, gp_double)
  if (nrow(double)) {
    double <- double[order(double$mr_pval, double$mean_score), , drop = FALSE]
    rownames(double) <- NULL
  }
  .log("Double-evidence candidates: ", nrow(double),
       " (drug_target=", nrow(dt_double), ", genetic_perturbation=", nrow(gp_double), ")")
} else {
  .log("MR targets file absent (", mr_path, "); double-evidence left empty pending Step7.", level = "WARN")
}
save_tsv(double, file.path(PARAMS$dir_cmap, "double_evidence_candidates.tsv"))

cat(sprintf("\nStep9 CMap coda done: %d raw hits -> %d cross-cell-line candidates -> %d double-evidence.\n",
            nrow(raw), nrow(cand), nrow(double)))
cat("NOTE: hypothesis-generating coda only; no efficacy claimed (cf. Brey 2015 PNAS, PMID 26420877).\n")

} else {
  # -------------------------------------------------------------------------
  # No CMap API returned results (offline/timeout). Write empty coda outputs
  # WITH schema so downstream steps parse cleanly; end naturally (no quit()).
  # -------------------------------------------------------------------------
  .log("No CMap API returned results (offline/timeout). Writing empty coda outputs with schema; rerun when online.", level = "WARN")
  empty_raw <- data.frame(pert_desc = character(0), pert_id = character(0), cell_line = character(0),
                          dose = character(0), time = character(0), score = numeric(0),
                          mode = character(0), source_api = character(0),
                          drug = character(0), stringsAsFactors = FALSE)
  save_tsv(empty_raw, file.path(PARAMS$dir_cmap, "l1000_raw_hits.tsv"))
  empty_cand <- data.frame(drug = character(0), n_celllines = integer(0), n_sources = integer(0),
                           mean_score = numeric(0), best_mode = character(0),
                           cell_lines = character(0), apis = character(0),
                           connectivity = character(0), stringsAsFactors = FALSE)
  save_tsv(empty_cand, file.path(PARAMS$dir_cmap, "cmap_candidates.tsv"))
  empty_double <- data.frame(drug = character(0), target_gene = character(0),
                             match_type = character(0),
                             n_celllines = integer(0), mean_score = numeric(0),
                             protective_direction = character(0), mr_pval = numeric(0),
                             coloc_pp4 = numeric(0), finan_tier = character(0),
                             stringsAsFactors = FALSE)
  save_tsv(empty_double, file.path(PARAMS$dir_cmap, "double_evidence_candidates.tsv"))
}

dump_session()
