# R/step04_success_projection.R
# Step 4 — Success-side cross-context validation of the regeneration-competent core program.
# Question: does the WGCNA/RRA "core program" (04_modules_wgcna/core_program.tsv), split into
# up/down gene sets, ALSO go up when jaw/dental regeneration is driven by an INDEPENDENT
# success mechanism (mechanical strain distraction osteogenesis, socket healing, alpha-KG)?
# We score the program with ssGSEA (GSVA) on three success-side datasets and test the
# osteogenic-score separability "regeneration/osteo vs control" (Wilcoxon + AUC/ROC).
#
# Datasets (success side):
#   GSE104473  mouse mandible distraction osteogenesis (DO)  -> ortholog-map to human first
#   GSE223778  extraction socket healing (human+mouse)        -> mouse part ortholog-mapped
#   GSE316449  PDLF alpha-KG osteogenic, human (EXTERNAL independent validation set)
#
# Output:
#   05_projection/success_scores.tsv   (per-sample up/down/osteo ssGSEA scores + inferred group)
#   05_projection/success_auc.tsv      (per-dataset Wilcoxon p + AUC for osteo-score separability)
#
# Note: where group metadata is not explicit, the group is AUTO-GUESSED from pData and the
# row is flagged group_confirmed=FALSE -> requires manual confirmation before claims are made.

source("R/utils.R"); init_pipeline()
suppressMessages({ library(GSVA); library(pROC) })

# ---- core program gene sets (up / down) -------------------------------------
core <- utils::read.delim(file.path(PARAMS$dir_wgcna, "core_program.tsv"), stringsAsFactors = FALSE)
stopifnot(all(c("gene", "direction") %in% colnames(core)))
core_sets_human <- list(
  core_up   = unique(core$gene[core$direction == "up"]),
  core_down = unique(core$gene[core$direction == "down"])
)
.log("core program: ", length(core_sets_human$core_up), " up / ",
     length(core_sets_human$core_down), " down genes (human symbols)")

# Success-side datasets to project onto, with the organism we must homogenize to human.
success_targets <- data.frame(
  accession = c("GSE104473", "GSE223778", "GSE316449"),
  organism  = c("mouse",     "human_mouse", "human"),
  context   = c("mandible_distraction_osteogenesis", "extraction_socket_healing", "PDLF_alphaKG"),
  stringsAsFactors = FALSE
)

# ---- auto group inference from pData (osteo/regeneration vs control) ---------
# Returns a factor with values in {regen, control, NA}; NA where text is ambiguous.
# "regen" = osteogenic induction / distracted / healing / later timepoint vs day0.
regen_kw   <- "osteogen|mineraliz|distract|regenerat|heal|union|consolidat|alpha[- ]?kg|akg|\\bOM\\b|osteo[- ]?induc|treated|day ?(3|5|7|10|14|21|28)\\b"
control_kw <- "control|ctrl|sham|undifferentiat|growth medium|\\bGM\\b|basal|non[- ]?induc|uninduc|vehicle|untreated|intact|day ?0\\b|\\bD0\\b|wild[- ]?type|\\bWT\\b"

guess_regen_group <- function(text) {
  t <- tolower(text)
  has_r <- grepl(regen_kw, t, perl = TRUE)
  has_c <- grepl(control_kw, t, perl = TRUE)
  if (has_r && !has_c) return("regen")
  if (has_c && !has_r) return("control")
  if (has_r && has_c)  return(NA_character_)   # both keywords -> needs manual confirmation
  NA_character_
}

# Pull a descriptive string per sample from pData characteristics/title/source.
sample_desc <- function(eset) {
  pd <- Biobase::pData(eset)
  cols <- grep("characteristics|title|source_name|description|:ch1|treatment|group|genotype|time|agent",
               colnames(pd), ignore.case = TRUE, value = TRUE)
  if (!length(cols)) cols <- colnames(pd)
  desc <- apply(pd[, cols, drop = FALSE], 1, function(r) paste(na.omit(as.character(r)), collapse = " | "))
  stats::setNames(desc, rownames(pd))
}

# ---- score one dataset with ssGSEA on human-homogenized expression ----------
score_dataset <- function(acc, organism, context) {
  ri <- which(REG$accession == acc)
  if (!length(ri)) { .log("registry miss for ", acc, level = "WARN"); return(NULL) }
  role <- REG$role[ri]
  eset <- tryCatch(fetch_geo(acc, role), error = function(e) { .log("fetch fail ", acc, ": ", conditionMessage(e), level = "WARN"); NULL })
  if (is.null(eset)) return(NULL)

  expr <- tryCatch(get_expr_symbol(eset), error = function(e) { .log("expr fail ", acc, ": ", conditionMessage(e), level = "WARN"); NULL })
  if (is.null(expr) || !nrow(expr)) {
    .log(acc, ": no symbol-level matrix (likely raw rnaseq supp not co-loaded) -> TODO manual quant", level = "WARN")
    return(NULL)
  }

  # Homogenize gene space to human. For mouse / mixed mouse rows, map mouse symbols -> human.
  if (organism == "mouse") {
    hs <- to_human_symbols(rownames(expr), from = "mouse")
    keep <- !is.na(hs) & hs != ""
    expr <- expr[keep, , drop = FALSE]; rownames(expr) <- hs[keep]
    expr <- expr[!duplicated(rownames(expr)), , drop = FALSE]
  } else if (organism == "human_mouse") {
    # mixed series: rows already symbol; ortholog-map any rows that fail to match human core,
    # but since core is human, also try mouse->human and union (keep human if collision).
    hs <- to_human_symbols(rownames(expr), from = "mouse")
    mapped <- ifelse(!is.na(hs) & hs != "", hs, rownames(expr))   # fall back to original (already-human) symbol
    rownames(expr) <- mapped
    expr <- expr[!duplicated(rownames(expr)), , drop = FALSE]
  }
  # human stays as-is

  # ssGSEA scoring of up/down core sets
  par <- ssgseaParam(as.matrix(expr), core_sets_human, normalize = TRUE)
  es  <- tryCatch(gsva(par),
                  error = function(e) {
                    # backward-compat for older GSVA signature
                    gsva(as.matrix(expr), core_sets_human, method = PARAMS$ssgsea_method, ssgsea.norm = TRUE, verbose = FALSE)
                  })
  es <- as.matrix(es)

  # group inference
  desc <- sample_desc(eset)[colnames(expr)]
  grp  <- vapply(desc, guess_regen_group, character(1))
  confirmed <- !is.na(grp)              # auto-guess -> not yet human-confirmed
  osteo_score <- es["core_up", ] - es["core_down", ]   # signed osteogenic activation

  data.frame(
    accession      = acc,
    context        = context,
    sample         = colnames(expr),
    description    = unname(desc),
    inferred_group = ifelse(is.na(grp), "unassigned", grp),
    group_confirmed = FALSE,             # ALWAYS FALSE here: auto-guessed, needs manual confirmation
    score_up       = es["core_up", ],
    score_down     = es["core_down", ],
    osteo_score    = osteo_score,
    stringsAsFactors = FALSE
  )
}

scores <- do.call(rbind, lapply(seq_len(nrow(success_targets)), function(i)
  score_dataset(success_targets$accession[i], success_targets$organism[i], success_targets$context[i])))

# Empty-schema definitions so downstream readers see stable column names either way.
empty_scores_schema <- data.frame(
  accession = character(0), context = character(0), sample = character(0),
  description = character(0), inferred_group = character(0), group_confirmed = logical(0),
  score_up = numeric(0), score_down = numeric(0), osteo_score = numeric(0),
  stringsAsFactors = FALSE)
empty_auc_schema <- data.frame(
  accession = character(0), context = character(0), n_regen = integer(0),
  n_control = integer(0), wilcox_p = numeric(0), auc = numeric(0),
  provisional = logical(0), note = character(0), stringsAsFactors = FALSE)

if (is.null(scores) || !nrow(scores)) {
  .log("No success-side dataset could be scored (expression matrices not available locally). ",
       "Download/quantify the rnaseq supplementaries, then re-run.", level = "WARN")
  save_tsv(empty_scores_schema, file.path(PARAMS$dir_proj, "success_scores.tsv"))
  save_tsv(empty_auc_schema, file.path(PARAMS$dir_proj, "success_auc.tsv"))
} else {
  save_tsv(scores, file.path(PARAMS$dir_proj, "success_scores.tsv"))

  # ---- separability: Wilcoxon + AUC of osteo_score (regen vs control) --------
  # Gate quantitative testing on manually-confirmed arms. When all groups are
  # auto-inferred (group_confirmed all FALSE), still compute separability but
  # flag provisional=TRUE and route to success_auc_PROVISIONAL.tsv so the numbers
  # are NOT mistaken for confirmed results.
  any_confirmed <- any(scores$group_confirmed)
  provisional <- !any_confirmed

  auc_rows <- list()
  for (acc in unique(scores$accession)) {
    s <- scores[scores$accession == acc & scores$inferred_group %in% c("regen", "control"), ]
    ng <- table(s$inferred_group)
    if (length(ng) < 2 || any(ng < 2)) {
      auc_rows[[acc]] <- data.frame(accession = acc,
        context = scores$context[match(acc, scores$accession)],
        n_regen = unname(ng["regen"]) %||% 0L, n_control = unname(ng["control"]) %||% 0L,
        wilcox_p = NA_real_, auc = NA_real_, provisional = provisional,
        note = "insufficient inferred groups; CONFIRM arms manually before testing",
        stringsAsFactors = FALSE)
      next
    }
    wt <- stats::wilcox.test(osteo_score ~ inferred_group, data = s)
    roc_obj <- pROC::roc(response = s$inferred_group, predictor = s$osteo_score,
                         levels = c("control", "regen"), direction = "<", quiet = TRUE)
    auc_rows[[acc]] <- data.frame(accession = acc,
      context  = s$context[1],
      n_regen  = unname(ng["regen"]), n_control = unname(ng["control"]),
      wilcox_p = wt$p.value, auc = as.numeric(pROC::auc(roc_obj)),
      provisional = provisional,
      note = if (provisional)
        "groups AUTO-INFERRED from pData (provisional); hypothesis-generating until arms manually confirmed"
      else
        "groups manually confirmed",
      stringsAsFactors = FALSE)
  }
  auc_tab <- do.call(rbind, auc_rows)

  # Route output by confirmation status: provisional results never land in the
  # canonical success_auc.tsv (write an empty-schema placeholder there instead).
  if (provisional) {
    save_tsv(auc_tab, file.path(PARAMS$dir_proj, "success_auc_PROVISIONAL.tsv"))
    save_tsv(empty_auc_schema, file.path(PARAMS$dir_proj, "success_auc.tsv"))
    .log("All success-side groups auto-inferred (group_confirmed=FALSE): AUC written to ",
         "success_auc_PROVISIONAL.tsv; canonical success_auc.tsv left as empty schema.", level = "WARN")
  } else {
    save_tsv(auc_tab, file.path(PARAMS$dir_proj, "success_auc.tsv"))
  }

  cat("\n=============== Step4 success-side projection ===============\n")
  print(auc_tab[, c("accession", "context", "n_regen", "n_control", "wilcox_p", "auc", "provisional")])
  if (provisional) {
    cat("\nNOTE: groups were auto-inferred from GEO metadata (group_confirmed=FALSE).\n",
        "AUC table is PROVISIONAL (success_auc_PROVISIONAL.tsv). Manually confirm regen/control\n",
        "arms in success_scores.tsv before any quantitative claim.\n", sep = "")
  }
}

dump_session()
.log("Step4 done. Success-side core-program activation written to ", PARAMS$dir_proj)
