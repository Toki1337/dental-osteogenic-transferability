# R/step07_mr_coloc.R
# Step 7 — Druggable-genome Mendelian randomization (MR) + colocalization.
#
# WHAT THIS ANCHORS — read before interpreting any output:
#   This step does NOT and CANNOT demonstrate a causal effect on MRONJ / jaw
#   regeneration failure. There is no MRONJ GWAS of usable size. What MR here
#   establishes is a narrower, defensible claim: the subset of failure-specific
#   candidate genes that are (i) druggable (Finan 2017) AND (ii) have human
#   genetic causal support for SYSTEMIC BONE (heel eBMD, Morris & Kemp 2019;
#   replicated in FinnGen fracture/osteoporosis). Eligibility for therapeutic
#   prioritisation is "druggable + genetically causal for bone", extrapolated to
#   the jaw context as hypothesis-generating, not proof.
#
#   The "protective direction" call (raises bone / opposes failure) is NOT taken
#   mechanically from the MR sign. cis-eQTL MR sign reports the effect of *higher
#   expression* on eBMD; whether higher expression of a given gene is biologically
#   protective in the jaw must be argued per gene from pathway direction in
#   Steps 2-6 (RRA direction, module membership, scRNA cell-type context). The
#   MR sign is recorded; the protective label is flagged for manual pathway
#   adjudication (column needs_pathway_adjudication) — never auto-asserted.
#
# Pipeline:
#   candidates  = 05_projection/failure_specific_program.tsv  ∩  druggable genome
#   instruments = eQTLGen cis-eQTL (p<mr_instrument_p, LD-clumped)
#   outcome     = eBMD GWAS (GWAS Catalog GCST006979)
#   MR          = TwoSampleMR IVW + Egger + weighted median + sensitivity
#   coloc       = coloc::coloc.abf per locus (PP.H4 > coloc_pp4_min); SMR-HEIDI hook
#   replication = FinnGen fracture/osteoporosis IVW
#
# Output:
#   07_mr_coloc/mr_results.tsv                       (all genes, all methods + sensitivity)
#   07_mr_coloc/coloc_results.tsv                    (per-gene coloc PP.H0..H4)
#   07_mr_coloc/mr_protective_druggable_targets.tsv  (genetic causal ∩ coloc ∩ druggable,
#                                                     protective direction FLAGGED for review)
#
# Data the user must place under 00_data/genetic/ (see READ-ME notes in code below):
#   eqtlgen_cis.tsv     eQTLGen cis-eQTL summary stats (full or per-gene)
#   ebmd_gwas.tsv       eBMD GWAS summary stats (GCST006979)
#   finngen_bone.tsv    FinnGen fracture/osteoporosis endpoint summary stats (replication)

source("R/utils.R"); init_pipeline()
suppressMessages({
  library(TwoSampleMR)
  library(coloc)
})

# ----------------------------------------------------------------------------
# 0. Config constants (file path conventions; kept here, not hard-coded in body)
# ----------------------------------------------------------------------------
DRUGGABLE_GENOME_TSV <- "config/druggable_genome.tsv"     # Finan 2017 list (placeholder shipped)
GENETIC_DIR          <- file.path(PARAMS$dir_data, "genetic")
EQTLGEN_FILE         <- file.path(GENETIC_DIR, "eqtlgen_cis.tsv")   # cis-eQTL instruments
EBMD_FILE            <- file.path(GENETIC_DIR, "ebmd_gwas.tsv")     # outcome (GCST006979)
FINNGEN_FILE         <- file.path(GENETIC_DIR, "finngen_bone.tsv")  # replication outcome
CIS_WINDOW_KB        <- 1000L   # +/- window around gene TSS defining the cis region

# Fallback sample sizes used ONLY when the summary stats carry no per-SNP N column.
# Preferred behaviour is to read N from the data; these are last-resort constants
# (logged with a WARN when used) and the N actually used is recorded in coloc output.
EQTLGEN_N_ASSUMED    <- 31684L    # eQTLGen 2019 cis-eQTL discovery N (Vosa et al.)
EBMD_N_ASSUMED       <- 426824L   # Morris & Kemp 2019 heel eBMD discovery N

# Built-in fallback druggable list (illustrative bone/druggable targets only).
# Used ONLY if config/druggable_genome.tsv is missing, so the script still runs.
DRUGGABLE_FALLBACK <- c(
  "TNFRSF11B", "TNFSF11", "SOST", "LRP5", "LRP6", "WNT16", "DKK1",
  "CTSK", "PTH1R", "ALPL", "ACP5", "BMP2", "BMPR1A", "ENPP1",
  "GGPS1", "FDPS"
)

# ----------------------------------------------------------------------------
# 1. Candidate genes = failure-specific program ∩ druggable genome
# ----------------------------------------------------------------------------
load_druggable <- function(path = DRUGGABLE_GENOME_TSV) {
  if (!file.exists(path)) {
    .log("druggable genome file not found at ", path, level = "WARN")
    .log("TODO(user): download Finan 2017 (Sci Transl Med 9:eaag1166) Table S1 and ",
         "write it to ", path, " with columns gene<TAB>tier. ",
         "Using built-in ", length(DRUGGABLE_FALLBACK), "-gene EXAMPLE list (NOT complete).",
         level = "WARN")
    return(data.frame(gene = DRUGGABLE_FALLBACK, tier = NA_integer_, stringsAsFactors = FALSE))
  }
  d <- utils::read.delim(path, stringsAsFactors = FALSE, comment.char = "#", check.names = FALSE)
  stopifnot("gene" %in% colnames(d))
  d$tier <- if ("tier" %in% colnames(d)) suppressWarnings(as.integer(d$tier)) else NA_integer_
  d <- d[!is.na(d$gene) & d$gene != "", c("gene", "tier")]
  if (nrow(d) < 200)
    .log("druggable list has only ", nrow(d), " genes — looks like the placeholder. ",
         "TODO(user): replace with the full Finan 2017 list before final results.",
         level = "WARN")
  unique(d)
}

failure_file <- file.path(PARAMS$dir_proj, "failure_specific_program.tsv")
if (!file.exists(failure_file))
  stop("missing upstream input: ", failure_file,
       " (Step5 failure-specific program). Run Step5 first.")
failure <- utils::read.delim(failure_file, stringsAsFactors = FALSE)
stopifnot("gene" %in% colnames(failure))

druggable <- load_druggable()
candidates <- intersect(unique(failure$gene), druggable$gene)
.log("MR candidate set: ", length(candidates), " genes = failure-specific ∩ druggable ",
     "(", nrow(druggable), " druggable; ", length(unique(failure$gene)), " failure-specific)")
if (length(candidates) == 0)
  stop("no overlap between failure-specific program and druggable genome; nothing to test.")

# carry the upstream failure-side direction (for downstream pathway adjudication)
fail_dir <- if ("direction" %in% colnames(failure))
  stats::setNames(failure$direction, failure$gene)[candidates] else stats::setNames(rep(NA, length(candidates)), candidates)

# ----------------------------------------------------------------------------
# 2. Summary-stat readers / formatters (TwoSampleMR::format_data skeletons)
#    These intentionally fail loud with an actionable message if the user has
#    not yet placed the harmonised summary-stat files. Column maps below assume
#    the public release layouts; adjust col names to match the user's download.
# ----------------------------------------------------------------------------

# ---- eQTLGen cis-eQTL -> TwoSampleMR "exposure" frame, restricted to one gene ----
# Expected eQTLGen full-cis columns (2019 release):
#   SNP, Pvalue, Zscore, AssessedAllele, OtherAllele, GeneSymbol, SNPChr, SNPPos,
#   NrSamples, AlleleB_all (allele frequency of the assessed allele).
#
# eQTLGen reports a Z-score, not a per-allele beta/se. We convert Z -> beta/se on
# the standardised-expression scale using the standard Zhu et al. 2016 (Nat Genet,
# SMR) approximation, which requires allele frequency p and sample size N:
#     beta = Z / sqrt(2*p*(1-p)*(N + Z^2))
#     se   = 1 / sqrt(2*p*(1-p)*(N + Z^2))
# This gives beta/se on a proper per-allele (standardised) scale so that IVW over
# MULTIPLE SNPs is dimensionally consistent and inference (se, CI) is valid.
#
# Only if allele frequency is genuinely missing do we fall back to a SINGLE best
# cis-SNP Wald-ratio mode: keep the top SNP, pass beta=Z, se=1. In that mode the
# point estimate is interpretable as a ratio but its MAGNITUDE is on a
# standardised, NOT per-allele, scale — flagged loudly in the log.
read_eqtlgen_exposure <- function(gene, path = EQTLGEN_FILE) {
  if (!file.exists(path))
    stop("eQTLGen file not found: ", path,
         "\n  TODO(user): download eQTLGen cis-eQTL summary stats from ",
         "https://www.eqtlgen.org/ and save (or per-gene subset) to ", path)
  e <- data.table::fread(path)
  gcol <- intersect(c("GeneSymbol", "Gene", "gene"), colnames(e))[1]
  e <- e[e[[gcol]] == gene, ]
  if (nrow(e) == 0) return(NULL)
  # genome-wide-for-cis instrument threshold
  pcol <- intersect(c("Pvalue", "P", "pval"), colnames(e))[1]
  e <- e[e[[pcol]] < PARAMS$mr_instrument_p, ]
  if (nrow(e) == 0) return(NULL)

  scol  <- intersect(c("SNP", "rsid"), colnames(e))[1]
  zcol  <- intersect(c("Zscore", "Z", "zscore"), colnames(e))[1]
  eacol <- intersect(c("AssessedAllele", "effect_allele", "A1"), colnames(e))[1]
  oacol <- intersect(c("OtherAllele", "other_allele", "A2"), colnames(e))[1]
  ncol  <- intersect(c("NrSamples", "N", "n", "NrSamplesUsed"), colnames(e))[1]
  # frequency of the assessed (effect) allele; p used below is treated as MAF-style
  # in the variance term 2*p*(1-p), which is symmetric in p vs 1-p.
  fcol  <- intersect(c("AlleleB_all", "freq", "MAF", "EAF", "eaf", "AF"), colnames(e))[1]

  has_native_betase <- all(c("beta", "se") %in% colnames(e))
  N_e <- if (!is.na(ncol)) as.numeric(e[[ncol]]) else rep(EQTLGEN_N_ASSUMED, nrow(e))
  if (is.na(ncol))
    .log("eQTLGen has no N column for ", gene, "; using assumed N=", EQTLGEN_N_ASSUMED,
         ", verify", level = "WARN")
  freq <- if (!is.na(fcol)) suppressWarnings(as.numeric(e[[fcol]])) else rep(NA_real_, nrow(e))

  if (has_native_betase) {
    df <- data.frame(
      SNP          = e[[scol]],
      beta         = as.numeric(e$beta),
      se           = as.numeric(e$se),
      effect_allele= e[[eacol]],
      other_allele = e[[oacol]],
      eaf          = freq,
      pval         = e[[pcol]],
      Phenotype    = paste0("eQTL_", gene),
      gene         = gene,
      stringsAsFactors = FALSE
    )
  } else if (!is.na(fcol) && any(is.finite(freq) & freq > 0 & freq < 1)) {
    # Zhu 2016 Z -> beta/se conversion (per-allele standardised scale)
    z   <- as.numeric(e[[zcol]])
    p   <- freq
    den <- 2 * p * (1 - p) * (N_e + z^2)
    den[!is.finite(den) | den <= 0] <- NA_real_
    beta <- z / sqrt(den)
    se   <- 1 / sqrt(den)
    df <- data.frame(
      SNP          = e[[scol]],
      beta         = beta,
      se           = se,
      effect_allele= e[[eacol]],
      other_allele = e[[oacol]],
      eaf          = p,
      pval         = e[[pcol]],
      Phenotype    = paste0("eQTL_", gene),
      gene         = gene,
      stringsAsFactors = FALSE
    )
    df <- df[is.finite(df$beta) & is.finite(df$se) & df$se > 0, , drop = FALSE]
    if (nrow(df) == 0) return(NULL)
  } else {
    # No usable allele frequency: fall back to SINGLE best cis-SNP Wald ratio.
    .log("eQTLGen has no usable allele frequency for ", gene,
         "; falling back to single best cis-SNP Wald-ratio mode (beta=Z, se=1). ",
         "Effect sizes are on a standardized scale, magnitudes are NOT per-allele.",
         level = "WARN")
    e1 <- e[which.min(e[[pcol]]), , drop = FALSE]
    df <- data.frame(
      SNP          = e1[[scol]],
      beta         = as.numeric(e1[[zcol]]),
      se           = 1,
      effect_allele= e1[[eacol]],
      other_allele = e1[[oacol]],
      eaf          = NA_real_,
      pval         = e1[[pcol]],
      Phenotype    = paste0("eQTL_", gene),
      gene         = gene,
      stringsAsFactors = FALSE
    )
  }

  exp <- TwoSampleMR::format_data(df, type = "exposure")
  # LD clump cis instruments (independent signals)
  exp <- tryCatch(
    TwoSampleMR::clump_data(exp, clump_r2 = PARAMS$mr_clump_r2, clump_kb = PARAMS$mr_clump_kb),
    error = function(err) { .log("clump failed for ", gene, " (", conditionMessage(err),
                                 "); using unclumped instruments", level = "WARN"); exp })
  exp
}

# ---- eBMD GWAS (GCST006979) -> TwoSampleMR "outcome" frame for a SNP set ----
# Expected harmonised columns (Morris & Kemp 2019 release / GWAS Catalog):
#   SNP/rsid, BETA, SE, P, EA (effect allele), NEA (non-effect allele), EAF
read_outcome_gwas <- function(snps, path, label) {
  if (!file.exists(path))
    stop(label, " GWAS file not found: ", path,
         "\n  TODO(user): place harmonised summary stats at ", path)
  g <- data.table::fread(path)
  scol <- intersect(c("SNP", "rsid", "variant_id"), colnames(g))[1]
  g <- g[g[[scol]] %in% snps, ]
  if (nrow(g) == 0) return(NULL)
  eaf_col <- intersect(c("EAF", "eaf", "A1FREQ"), colnames(g))[1]   # may be NA if absent
  df <- data.frame(
    SNP          = g[[scol]],
    beta         = g[[intersect(c("BETA", "beta", "Effect"), colnames(g))[1]]],
    se           = g[[intersect(c("SE", "se", "StdErr"), colnames(g))[1]]],
    effect_allele= g[[intersect(c("EA", "effect_allele", "A1", "ALLELE1"), colnames(g))[1]]],
    other_allele = g[[intersect(c("NEA", "other_allele", "A2", "ALLELE0"), colnames(g))[1]]],
    eaf          = if (!is.na(eaf_col)) g[[eaf_col]] else NA_real_,
    pval         = g[[intersect(c("P", "pval", "P_BOLT_LMM"), colnames(g))[1]]],
    Phenotype    = label,
    stringsAsFactors = FALSE
  )
  TwoSampleMR::format_data(df, type = "outcome")
}

# ----------------------------------------------------------------------------
# 3. Per-gene MR (IVW + Egger + weighted median) with sensitivity analyses
# ----------------------------------------------------------------------------
run_mr_gene <- function(gene, outcome_path, outcome_label) {
  exp <- tryCatch(read_eqtlgen_exposure(gene), error = function(e) {
    .log("exposure read failed for ", gene, ": ", conditionMessage(e), level = "WARN"); NULL })
  if (is.null(exp) || nrow(exp) == 0) return(NULL)
  out <- tryCatch(read_outcome_gwas(exp$SNP, outcome_path, outcome_label),
                  error = function(e) { .log(outcome_label, " read failed for ", gene, ": ",
                                             conditionMessage(e), level = "WARN"); NULL })
  if (is.null(out) || nrow(out) == 0) return(NULL)
  # action=2 aligns ambiguous (palindromic) SNPs using EAF; this needs an exposure
  # EAF. If the exposure carries no EAF (single-SNP Wald-ratio fallback, or a feed
  # without frequency), use action=3 which simply DROPS ambiguous palindromes.
  has_exp_eaf <- "eaf.exposure" %in% colnames(exp) && any(!is.na(exp$eaf.exposure))
  harm_action <- if (has_exp_eaf) 2L else 3L
  dat <- TwoSampleMR::harmonise_data(exp, out, action = harm_action)
  n_dropped <- sum(!dat$mr_keep, na.rm = TRUE)
  if (harm_action == 3L)
    .log("harmonise action=3 (no exposure EAF) for ", gene, ": dropped ", n_dropped,
         " ambiguous/unalignable SNP(s) of ", nrow(dat), level = "WARN")
  dat <- dat[dat$mr_keep, , drop = FALSE]
  if (nrow(dat) == 0) return(NULL)

  res <- TwoSampleMR::mr(dat, method_list = PARAMS$mr_methods)
  # With a single instrument TwoSampleMR returns "Wald ratio" instead of IVW.
  # Downstream the primary causal set is filtered on the IVW label, which would
  # silently drop single-cis-SNP genes. Relabel the Wald-ratio row into the IVW
  # slot (it is the single-SNP IVW estimate) and flag it so the instrument-count
  # distribution stays auditable.
  res$mr_method_native <- res$method
  res$is_wald_ratio <- res$method == "Wald ratio"
  res$method[res$is_wald_ratio] <- "Inverse variance weighted"
  res$gene <- gene; res$n_snp_total <- nrow(dat); res$harmonise_action <- harm_action
  res$n_snp_dropped <- n_dropped

  # ---- sensitivity ----
  egger_int <- tryCatch(TwoSampleMR::mr_pleiotropy_test(dat), error = function(e) NULL)
  het       <- tryCatch(TwoSampleMR::mr_heterogeneity(dat), error = function(e) NULL)
  res$egger_intercept   <- if (!is.null(egger_int)) egger_int$egger_intercept else NA_real_
  res$egger_intercept_p <- if (!is.null(egger_int)) egger_int$pval else NA_real_
  res$Q_pval <- if (!is.null(het)) het$Q_pval[match("Inverse variance weighted", het$method)] else NA_real_

  # leave-one-out + single-SNP (kept as side artefacts for the figure step; flag only)
  loo <- tryCatch(TwoSampleMR::mr_leaveoneout(dat), error = function(e) NULL)
  res$loo_robust <- if (!is.null(loo)) {
    ivw_loo <- loo[loo$SNP != "All", ]
    all(sign(ivw_loo$b) == sign(ivw_loo$b[1])) && all(ivw_loo$p < 0.05)
  } else NA
  single <- tryCatch(TwoSampleMR::mr_singlesnp(dat), error = function(e) NULL)
  res$n_single_sig <- if (!is.null(single)) sum(single$p[single$SNP != "All"] < 0.05, na.rm = TRUE) else NA_integer_

  res
}

mr_all <- list()
for (g in candidates) {
  r <- run_mr_gene(g, EBMD_FILE, "eBMD")
  if (!is.null(r)) mr_all[[g]] <- r
}
if (length(mr_all) > 0) {
  nsnp_by_gene <- vapply(mr_all, function(r) r$n_snp_total[1], numeric(1))
  n_single <- sum(nsnp_by_gene == 1)
  .log("MR instrument-count distribution across ", length(mr_all), " genes: ",
       "single-SNP (Wald ratio) = ", n_single, "; multi-SNP (IVW) = ",
       length(mr_all) - n_single, "; median nSNP = ",
       stats::median(nsnp_by_gene), ", range ", min(nsnp_by_gene), "-", max(nsnp_by_gene))
}
if (length(mr_all) == 0) {
  .log("no gene yielded a complete MR (instruments + outcome). ",
       "Check that 00_data/genetic/ files are present and column-mapped. ",
       "Writing empty result stubs.", level = "WARN")
  mr_df <- data.frame()
} else {
  mr_df <- do.call(rbind, mr_all)
  # FDR across IVW estimates only (primary causal test)
  ivw <- mr_df$method == "Inverse variance weighted"
  mr_df$ivw_fdr <- NA_real_
  mr_df$ivw_fdr[ivw] <- stats::p.adjust(mr_df$pval[ivw], method = "BH")
}
save_tsv(mr_df, file.path(PARAMS$dir_mr, "mr_results.tsv"))

# ----------------------------------------------------------------------------
# 4. Colocalization — coloc::coloc.abf per candidate locus (eQTL vs eBMD)
#    SMR-HEIDI hook left as an interface stub (smr_heidi_p) for users with the
#    SMR binary; coloc is the primary single-causal-variant test here.
# ----------------------------------------------------------------------------
# Reads the FULL cis region (not just instruments) for one gene from both
# eQTLGen and eBMD, aligned on SNP, and runs coloc.abf. Returns PP.H0..H4.
run_coloc_gene <- function(gene) {
  if (!file.exists(EQTLGEN_FILE) || !file.exists(EBMD_FILE)) return(NULL)
  e <- data.table::fread(EQTLGEN_FILE)
  gcol <- intersect(c("GeneSymbol", "Gene", "gene"), colnames(e))[1]
  e <- e[e[[gcol]] == gene, ]
  if (nrow(e) < 5) return(NULL)
  scol_e <- intersect(c("SNP", "rsid"), colnames(e))[1]
  pcol_e <- intersect(c("Pvalue", "P", "pval"), colnames(e))[1]
  ncol_e <- intersect(c("NrSamples", "N", "n", "NrSamplesUsed"), colnames(e))[1]
  bcol_e <- intersect(c("beta", "Beta"), colnames(e))[1]
  secol_e<- intersect(c("se", "SE", "StdErr"), colnames(e))[1]
  zcol_e <- intersect(c("Zscore", "Z", "zscore"), colnames(e))[1]
  freq_e <- intersect(c("AlleleB_all", "freq", "MAF", "EAF", "eaf", "AF"), colnames(e))[1]

  o <- data.table::fread(EBMD_FILE)
  scol_o <- intersect(c("SNP", "rsid", "variant_id"), colnames(o))[1]
  o <- o[o[[scol_o]] %in% e[[scol_e]], ]
  if (nrow(o) < 5) return(NULL)
  pcol_o <- intersect(c("P", "pval", "P_BOLT_LMM"), colnames(o))[1]
  bcol_o <- intersect(c("BETA", "beta", "Effect"), colnames(o))[1]
  secol_o<- intersect(c("SE", "se", "StdErr"), colnames(o))[1]
  eaf_o  <- intersect(c("EAF", "eaf", "A1FREQ"), colnames(o))[1]
  ncol_o <- intersect(c("N", "n", "N_total", "Neff"), colnames(o))[1]

  common <- intersect(e[[scol_e]], o[[scol_o]])
  if (length(common) < 5) return(NULL)
  e <- e[match(common, e[[scol_e]]), ]; o <- o[match(common, o[[scol_o]]), ]

  # ---- sample sizes: prefer per-SNP N from data; fall back to assumed constant ----
  if (!is.na(ncol_e)) {
    N_e <- stats::median(as.numeric(e[[ncol_e]]), na.rm = TRUE)
  } else {
    N_e <- EQTLGEN_N_ASSUMED
    .log("coloc ", gene, ": eQTLGen N column absent; using assumed N=", EQTLGEN_N_ASSUMED,
         ", verify", level = "WARN")
  }
  if (!is.na(ncol_o)) {
    N_o <- stats::median(as.numeric(o[[ncol_o]]), na.rm = TRUE)
  } else {
    N_o <- EBMD_N_ASSUMED
    .log("coloc ", gene, ": eBMD N column absent; using assumed N=", EBMD_N_ASSUMED,
         ", verify", level = "WARN")
  }

  # ---- eQTL dataset: prefer beta/varbeta; otherwise pvalues + MAF (never MAF=NULL
  #      with type=quant + pvalues, which coloc cannot handle) ----
  maf_e <- if (!is.na(freq_e)) {
    f <- suppressWarnings(as.numeric(e[[freq_e]])); pmin(f, 1 - f)
  } else NULL
  if (!is.na(bcol_e) && !is.na(secol_e)) {
    D_eqtl <- list(snp = common, type = "quant",
                   beta = as.numeric(e[[bcol_e]]), varbeta = as.numeric(e[[secol_e]])^2,
                   N = N_e, MAF = maf_e)
  } else if (!is.na(zcol_e) && !is.null(maf_e) && any(is.finite(maf_e) & maf_e > 0)) {
    # reconstruct beta/varbeta from Z + MAF + N (Zhu 2016), consistent with the
    # exposure reader, so coloc never runs with MAF=NULL.
    z   <- as.numeric(e[[zcol_e]])
    den <- 2 * maf_e * (1 - maf_e) * (N_e + z^2)
    den[!is.finite(den) | den <= 0] <- NA_real_
    D_eqtl <- list(snp = common, type = "quant",
                   beta = z / sqrt(den), varbeta = (1 / sqrt(den))^2,
                   N = N_e, MAF = maf_e)
  } else if (!is.null(maf_e)) {
    D_eqtl <- list(snp = common, type = "quant",
                   pvalues = as.numeric(e[[pcol_e]]), N = N_e, MAF = maf_e)
  } else {
    .log("coloc ", gene, ": no eQTL MAF and no beta/se -> cannot build a valid quant ",
         "dataset (coloc needs MAF with pvalues); skipping coloc for this gene.",
         level = "WARN")
    return(NULL)
  }

  D_bmd  <- list(snp = common, type = "quant",
                 beta = as.numeric(o[[bcol_o]]), varbeta = as.numeric(o[[secol_o]])^2,
                 MAF = if (!is.na(eaf_o)) { f <- as.numeric(o[[eaf_o]]); pmin(f, 1 - f) } else NULL,
                 N = N_o)

  # validate both datasets before running coloc.abf (catches MAF/varbeta issues)
  ok <- tryCatch({
    coloc::check_dataset(D_eqtl, suffix = "(eQTL)")
    coloc::check_dataset(D_bmd, suffix = "(eBMD)")
    TRUE
  }, error = function(err) {
    .log("coloc ", gene, ": check_dataset failed (", conditionMessage(err),
         "); skipping coloc for this gene", level = "WARN"); FALSE })
  if (!ok) return(NULL)

  cc <- tryCatch(coloc::coloc.abf(D_eqtl, D_bmd), error = function(e) NULL)
  if (is.null(cc)) return(NULL)
  s <- cc$summary
  data.frame(gene = gene, nsnps = unname(s["nsnps"]),
             N_eqtl = N_e, N_bmd = N_o,
             PP_H0 = unname(s["PP.H0.abf"]), PP_H1 = unname(s["PP.H1.abf"]),
             PP_H2 = unname(s["PP.H2.abf"]), PP_H3 = unname(s["PP.H3.abf"]),
             PP_H4 = unname(s["PP.H4.abf"]), stringsAsFactors = FALSE)
}

# SMR-HEIDI interface stub: if the user has the SMR binary + BESD eQTLGen, they
# can run `smr --beqtl-summary eqtlgen --gwas-summary ebmd ...` and filter on
# p_HEIDI > PARAMS$smr_heidi_p (HEIDI not rejected => single shared variant).
# We expose the threshold here and leave the call to the user's environment.
.log("SMR-HEIDI filter threshold available as PARAMS$smr_heidi_p = ", PARAMS$smr_heidi_p,
     " (interface stub; run SMR externally and merge p_HEIDI if available)")

coloc_all <- list()
for (g in unique(c(candidates, names(mr_all)))) {
  cr <- tryCatch(run_coloc_gene(g), error = function(e) NULL)
  if (!is.null(cr)) coloc_all[[g]] <- cr
}
coloc_df <- if (length(coloc_all)) do.call(rbind, coloc_all) else data.frame()
if (nrow(coloc_df)) coloc_df$coloc_pass <- coloc_df$PP_H4 > PARAMS$coloc_pp4_min
save_tsv(coloc_df, file.path(PARAMS$dir_mr, "coloc_results.tsv"))

# ----------------------------------------------------------------------------
# 5. FinnGen replication — re-run IVW on fracture/osteoporosis endpoint
# ----------------------------------------------------------------------------
repl_all <- list()
if (file.exists(FINNGEN_FILE)) {
  for (g in names(mr_all)) {
    r <- tryCatch(run_mr_gene(g, FINNGEN_FILE, "FinnGen_bone"), error = function(e) NULL)
    if (!is.null(r)) repl_all[[g]] <- r[r$method == "Inverse variance weighted", ]
  }
} else {
  .log("FinnGen replication file absent (", FINNGEN_FILE, "); skipping replication. ",
       "TODO(user): add FinnGen fracture/osteoporosis endpoint stats for replication.",
       level = "WARN")
}
repl_df <- if (length(repl_all)) do.call(rbind, repl_all) else data.frame()

# ----------------------------------------------------------------------------
# 6. Final keep-set: genetic causal (IVW FDR) ∩ coloc ∩ druggable.
#    Protective direction is FLAGGED, never auto-asserted (see header).
# ----------------------------------------------------------------------------
build_targets <- function() {
  if (!nrow(mr_df)) return(data.frame())
  ivw <- mr_df[mr_df$method == "Inverse variance weighted", ]
  ivw <- ivw[!is.na(ivw$ivw_fdr) & ivw$ivw_fdr < PARAMS$mr_fdr, , drop = FALSE]
  if (!nrow(ivw)) { .log("no gene passes IVW FDR; keep-set empty", level = "WARN"); return(data.frame()) }

  coloc_pass_genes <- if (nrow(coloc_df)) coloc_df$gene[coloc_df$coloc_pass] else character(0)
  keep <- ivw[ivw$gene %in% coloc_pass_genes, , drop = FALSE]
  if (!nrow(keep)) { .log("genes pass MR but none pass coloc PP.H4>", PARAMS$coloc_pp4_min,
                          "; keep-set empty (report as hypothesis-generating)", level = "WARN")
    return(data.frame()) }

  out <- data.frame(
    gene            = keep$gene,
    mr_beta_eBMD    = keep$b,
    mr_se           = keep$se,
    mr_pval         = keep$pval,
    mr_ivw_fdr      = keep$ivw_fdr,
    n_snp           = keep$nsnp,
    # mr_sign: effect of higher cis-expression on eBMD (+ = raises bone density).
    # This is a CANDIDATE protective signal ONLY — must be reconciled with the
    # failure-side direction and pathway biology before claiming protective.
    mr_sign_raises_bmd = ifelse(keep$b > 0, "yes", "no"),
    failure_side_direction = unname(fail_dir[keep$gene]),
    coloc_PP_H4     = coloc_df$PP_H4[match(keep$gene, coloc_df$gene)],
    druggable_tier  = druggable$tier[match(keep$gene, druggable$gene)],
    stringsAsFactors = FALSE
  )
  # replication annotation
  if (nrow(repl_df)) {
    out$finngen_ivw_p   <- repl_df$pval[match(out$gene, repl_df$gene)]
    out$finngen_concord <- ifelse(!is.na(out$finngen_ivw_p) & out$finngen_ivw_p < 0.05 &
                                    sign(repl_df$b[match(out$gene, repl_df$gene)]) == sign(out$mr_beta_eBMD),
                                  "yes", "no")
  } else { out$finngen_ivw_p <- NA_real_; out$finngen_concord <- NA }

  # Protective label is NOT auto-set: flag for manual pathway adjudication.
  out$needs_pathway_adjudication <- TRUE
  out$protective_call <- "PENDING_MANUAL_PATHWAY_REVIEW"
  out$interpretation_note <-
    "Genetic causal support for SYSTEMIC BONE (eBMD) + coloc + druggable. NOT MRONJ causal proof. Protective direction requires per-gene pathway argument (Steps 2-6)."

  # ---- consumer-facing column ALIASES (Step9/Step10 read these names) ----
  # We keep the honest native columns above AND expose standardized aliases so the
  # downstream integration does not silently read NA. Semantics are unchanged.
  out$protective_direction <- out$protective_call          # PENDING until manual pathway review (honest)
  out$pval        <- out$mr_pval
  out$padj        <- out$mr_ivw_fdr
  out$coloc_pp4   <- out$coloc_PP_H4
  out$heidi_p     <- NA_real_                              # filled if user runs external SMR-HEIDI
  out$finan_tier  <- out$druggable_tier
  out$finngen_replicated <- !is.na(out$finngen_concord) & out$finngen_concord == "yes"
  out$mr_method   <- "Inverse variance weighted"
  out$passes      <- TRUE                                  # all rows here passed IVW-FDR AND coloc gates
  out[order(out$mr_ivw_fdr), ]
}

targets <- build_targets()
save_tsv(targets, file.path(PARAMS$dir_mr, "mr_protective_druggable_targets.tsv"))

cat(sprintf(
  "\nStep7 MR+coloc: %d candidates tested | %d with complete MR | %d coloc-pass | %d final druggable causal targets.\n",
  length(candidates), length(mr_all),
  if (nrow(coloc_df)) sum(coloc_df$coloc_pass, na.rm = TRUE) else 0L,
  nrow(targets)))
cat("CAVEAT: targets are 'druggable + genetically causal for systemic bone (eBMD/FinnGen)',\n",
    "        extrapolated to jaw failure as HYPOTHESIS-GENERATING. Not MRONJ causal proof.\n",
    "        Protective direction flagged PENDING_MANUAL_PATHWAY_REVIEW — adjudicate per gene.\n", sep = "")

dump_session()
.log("Step7 done. See 07_mr_coloc/mr_protective_druggable_targets.tsv (review protective_call manually).")
