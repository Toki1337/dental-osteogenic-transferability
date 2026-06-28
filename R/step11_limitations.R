# R/step11_limitations.R
# Step 11 — Robustness & limitations self-audit (reproducible audit tables).
# ----------------------------------------------------------------------------
# This step does NOT generate new biology; it consolidates the project's own
# weaknesses and mitigations into two auditable, reviewer-facing tables, computed
# from the actual upstream artifacts and the dataset registry (never asserted).
#
# Outputs:
#   supp/limitations_audit.tsv    (per row: risk, evidence, mitigation_taken, residual_limitation)
#   supp/dataset_provenance.tsv   (per dataset: platform, n, role, organism, confound flags, usage)
#
# Honesty mandate: thin failure-side human data (GSE7116 = MM/N-BP confound, no
# healthy control; GSE303003 = n=2 scRNA + scoop risk) are flagged explicitly and
# every failure-leaning claim is labeled hypothesis-generating.

source("R/utils.R"); init_pipeline()

read_if <- function(path) if (file.exists(path)) utils::read.delim(path, stringsAsFactors = FALSE) else NULL
n_or_na <- function(x) if (is.null(x)) NA_integer_ else nrow(x)

# ===========================================================================
# A. dataset provenance table (from REG; objective, no claims)
# ===========================================================================
prov <- data.frame(
  accession   = REG$accession,
  role        = REG$role,
  cell_type   = REG$cell_type,
  organism    = REG$organism,
  platform    = REG$platform,
  modality    = REG$modality,
  n_samples   = REG$n_samples,
  year        = REG$year,
  in_rra_pool = REG$include_in_rra & REG$modality == "mrna",
  is_external_validation = REG$is_external_validation,
  status      = REG$status,
  stringsAsFactors = FALSE
)
# explicit confound / caveat flags (claimed up-front per audit instructions)
prov$confound_flag <- ""
prov$confound_flag[prov$accession == "GSE7116"]   <- "MM patients + N-BP/chemotherapy; NO healthy control; supportive-only, not a quantitative cohort"
prov$confound_flag[prov$accession == "GSE303003"] <- "scRNA n=2; pre-print scoop risk; mechanistic illustration only"
prov$confound_flag[prov$accession == "GSE269255"] <- "published mechanism (taurine axis); corroborative not first-discovery; mandible sourcing must be confirmed"
prov$confound_flag[prov$accession == "GSE295106"] <- "primary mouse failure axis (has control); cross-species -> human via ortholog mapping required"
prov$confound_flag[prov$accession == "GSE104473"] <- "reviewer-reported, not independently re-verified by this cohort; rare large n on success side"
prov$confound_flag[prov$accession == "DRA006607"] <- "DDBJ raw fastq (not GEO matrix); needs alignment/quant before use"
prov$confound_flag[prov$role == "neg_control"]    <- "negative control (lineage reprogramming); excluded from RRA by design"
prov$confound_flag[prov$role == "position"]       <- "position/context only; explains jaw osteogenic weakness, not in RRA training"

# how each dataset is actually used in the pipeline
prov$pipeline_use <- with(prov, ifelse(in_rra_pool, "RRA training (success/competent program)",
                               ifelse(is_external_validation, "external validation (osteogenic score)",
                               ifelse(role == "failure", "failure-side projection / scRNA (hypothesis-generating)",
                               ifelse(role %in% c("success_jaw"), "cross-context success-side reference",
                               ifelse(role == "anchor", "cross-source baseline anchor",
                               ifelse(role == "position", "jaw position-specific context",
                               ifelse(role == "neg_control", "negative control", "support"))))))))
save_tsv(prov, file.path(PARAMS$dir_supp, "dataset_provenance.tsv"))
.log("dataset_provenance.tsv: ", nrow(prov), " datasets catalogued")

# ===========================================================================
# B. quantitative robustness summaries (computed from real upstream artifacts)
# ===========================================================================
# B1. RRA leave-one-dataset-out stability distribution
loo <- read_if(file.path(PARAMS$dir_rra, "loo_stability.tsv"))
if (!is.null(loo) && "loo_fraction" %in% colnames(loo)) {
  loo_q <- stats::quantile(loo$loo_fraction, c(0, .25, .5, .75, 1), na.rm = TRUE)
  loo_stable <- mean(loo$loo_fraction >= PARAMS$loo_stability_min, na.rm = TRUE)
  loo_evidence <- sprintf("LOO fraction median=%.2f (IQR %.2f-%.2f); %.0f%% of signature genes stable at >=%.2f over %d genes",
                          loo_q[3], loo_q[2], loo_q[4], 100 * loo_stable, PARAMS$loo_stability_min, nrow(loo))
} else loo_evidence <- "loo_stability.tsv not yet produced (run Step2)"

# B2. platform / sample-size mix of the RRA training pool
rra_pool <- prov[prov$in_rra_pool, , drop = FALSE]
plat_mix <- paste(sprintf("%s=%d", names(table(rra_pool$platform)), as.integer(table(rra_pool$platform))), collapse = ", ")
n_rra <- nrow(rra_pool)
nmin  <- suppressWarnings(min(as.numeric(rra_pool$n_samples), na.rm = TRUE))
cross_platform_evidence <- sprintf("RRA pool = %d mrna datasets; platforms: %s; per-dataset n as low as %s; RRA chosen for platform-drift robustness + LOO gate",
                                   n_rra, plat_mix, ifelse(is.finite(nmin), nmin, "NA"))

# B3. failure-side sample size + confound load
fail <- prov[prov$role == "failure", , drop = FALSE]
failure_evidence <- sprintf("failure-side datasets = %d (%s); human bulk = GSE7116 (MM/N-BP confound, no healthy ctrl); human scRNA = GSE303003 (n=2, scoop risk); primary failure axis is mouse mandible scRNA (GSE295106 ctrl-vs-BRONJ, GSE269255 ORNJ)",
                            nrow(fail), paste(fail$accession, collapse = ", "))

# B4. MR instrument / sensitivity status (from Step7 artifacts if present)
mr_full <- read_if(file.path(PARAMS$dir_mr, "mr_results.tsv"))            # Step7 writes "mr_results.tsv" (not _full)
mr_tgt  <- read_if(file.path(PARAMS$dir_mr, "mr_protective_druggable_targets.tsv"))
if (!is.null(mr_tgt)) {
  n_pass <- if ("passes" %in% colnames(mr_tgt)) sum(tolower(mr_tgt$passes) %in% c("true","yes","1")) else nrow(mr_tgt)
  # coloc PP4 lives under the alias column "coloc_pp4" (Step7 also keeps native "coloc_PP_H4"); read either.
  coloc_col <- intersect(c("coloc_pp4", "coloc_PP_H4"), colnames(mr_tgt))[1]
  n_coloc <- if (!is.na(coloc_col)) sum(mr_tgt[[coloc_col]] >= PARAMS$coloc_pp4_min, na.rm = TRUE) else NA
  mr_evidence <- sprintf("druggable-genome MR: %d protective druggable targets passing (coloc PP4>=%.2f in %s, HEIDI p>%.2f); methods=%s; eBMD exposure + eQTLGen IVs + FinnGen replication",
                         n_pass, PARAMS$coloc_pp4_min, ifelse(is.na(n_coloc), "NA", n_coloc),
                         PARAMS$smr_heidi_p, paste(PARAMS$mr_methods, collapse = "/"))
} else {
  mr_evidence <- sprintf("MR not yet run (Step7); planned: eBMD (GCST006979) + eQTLGen cis-eQTL IVs, methods=%s, coloc PP4>=%.2f, HEIDI p>%.2f, FinnGen replication",
                         paste(PARAMS$mr_methods, collapse = "/"), PARAMS$coloc_pp4_min, PARAMS$smr_heidi_p)
}

# B5. CMap reproducibility note (coda)
cmap_cand <- read_if(file.path(PARAMS$dir_cmap, "cmap_candidates.tsv"))
double    <- read_if(file.path(PARAMS$dir_cmap, "double_evidence_candidates.tsv"))
cmap_evidence <- sprintf("CMap coda: %s cross-cell-line candidates (>= %d cell lines), %s double-evidence (MR INTERSECT CMap); known L1000 reproducibility limits acknowledged (cf. Brey 2015 paradigm)",
                         ifelse(is.null(cmap_cand), "0", n_or_na(cmap_cand)), PARAMS$cmap_min_celllines,
                         ifelse(is.null(double), "0", n_or_na(double)))

# ===========================================================================
# C. limitations audit table (risk / evidence / mitigation / residual)
# ===========================================================================
audit <- rbind(
  data.frame(risk = "Cross-platform hub drift in meta-analysis",
             evidence = cross_platform_evidence,
             mitigation_taken = "RRA (rank-based, platform-robust) + leave-one-dataset-out gate (>= loo_stability_min) + cell-type-stratified sensitivity; WGCNA only on batch-corrected homogeneous subset",
             residual_limitation = "Residual platform/cell-type heterogeneity; some rare sources (DFC/SHED, 2013 array) carry low cross-platform weight"),
  data.frame(risk = "RRA signature instability to dataset choice",
             evidence = loo_evidence,
             mitigation_taken = "leave-one-dataset-out stability fraction reported per gene; only LOO-stable genes enter the core program",
             residual_limitation = "Small RRA pool means dropping a single influential dataset can still shift borderline genes"),
  data.frame(risk = "Comparison-arm heterogeneity (odontogenic / recipient-cell / drug x osteo interaction arms)",
             evidence = "Step0 manual per-arm curation; sample_arms_final.tsv + rra_usable_datasets.tsv define the actual clean osteo-vs-control arms (reported, not assumed; may be < 14)",
             mitigation_taken = "Per-arm human governance; odontogenic-only and TF-reprogramming (GSE160451) excluded; negative control retained as such",
             residual_limitation = "Inducer chemistry differs across datasets (statin, flavonoid, shear, CBD, phosphate); the program is 'osteogenic-direction-shared', not inducer-identical"),
  data.frame(risk = "Failure-side human data scarce and confounded",
             evidence = failure_evidence,
             mitigation_taken = "GSE7116 MM/N-BP confound claimed explicitly and used supportive-only; primary failure axis = mouse mandible scRNA with controls; human-mouse conserved intersection via to_human_symbols(babelgene)",
             residual_limitation = "No clean human healthy-vs-MRONJ bulk cohort; failure-specific calls are hypothesis-generating"),
  data.frame(risk = "Single-cell small-sample / scoop risk",
             evidence = "GSE303003 scRNA n=2 (pre-print scoop risk); GSE295106/GSE269255 mouse mandible scRNA few replicates",
             mitigation_taken = "scRNA used for mechanistic hypothesis + cell-state corroboration only, explicitly labeled; increment framed as success-vs-failure framework, not re-description of any single scRNA paper",
             residual_limitation = "Proportion/communication shifts are illustrative; not powered for quantitative differential abundance"),
  data.frame(risk = "Cross-species (mouse failure axis) comparability",
             evidence = "Failure axis leans on mouse mandible (GSE295106/GSE269255); success-side mandible DO is mouse (GSE104473)",
             mitigation_taken = "Human/mouse never merged; all cross-species comparisons via ortholog mapping (to_human_symbols) then intersection; conservative conserved set reported",
             residual_limitation = "Ortholog mapping loses non-1:1 genes; mouse-only or human-only nodes are not cross-validated"),
  data.frame(risk = "No jaw-specific GWAS; 'raise bone = anti-MRONJ' direction not automatic",
             evidence = mr_evidence,
             mitigation_taken = "MR used only as a 'druggable + systemic-bone genetic causality' FILTER on eBMD/fracture; protective direction argued per-pathway, not by MR sign alone; jaw specificity supplied by transcriptome + scRNA + position datasets",
             residual_limitation = "Genetic causality is for systemic bone density, not MRONJ per se; MRONJ-specific causality remains unproven"),
  data.frame(risk = "MR pleiotropy / LD-driven false colocalization",
             evidence = sprintf("Sensitivity battery: %s + coloc PP4>=%.2f + SMR-HEIDI p>%.2f; FinnGen replication", paste(PARAMS$mr_methods, collapse = "/"), PARAMS$coloc_pp4_min, PARAMS$smr_heidi_p),
             mitigation_taken = "IVW + MR-Egger intercept + weighted median agreement; coloc + HEIDI to filter LD artifacts; clumping r2<",
             residual_limitation = "cis-eQTL instruments can still tag pleiotropic effects; bidirectional/weak-instrument bias not fully excludable"),
  data.frame(risk = "CMap repurposing reproducibility / incrementalism",
             evidence = cmap_evidence,
             mitigation_taken = "CMap demoted to coda; only MR-intersecting double-evidence retained; cross-cell-line consistency required; framed against Brey 2015 paradigm (PMID 26420877); no efficacy claimed",
             residual_limitation = "L1000 connectivity is noisy and cell-line dependent; candidate drugs are computational hypotheses only"),
  data.frame(risk = "2025-26 new deposits may have incomplete metadata",
             evidence = paste("datasets flagged in registry status:", paste(unique(REG$status[grepl("verify|small|confound|sra", REG$status)]), collapse = ", ")),
             mitigation_taken = "first-week per-accession GEO verification; methodology not dependent on any single dataset; unusable sets demoted to citation or replaced",
             residual_limitation = "Some accessions await full metadata release; provenance must be re-checked at submission"),
  data.frame(risk = "No wet-lab validation",
             evidence = "Dry-lab, resource/hypothesis-generating study; deliverables = reproducible signature + prioritized target table + candidate old-drug list",
             mitigation_taken = "Genre declared up-front; internal cross-validation (LOO + cross-context success-side + external osteogenic set + MR INTERSECT CMap double evidence) substitutes for external wet validation; wet steps listed as future work",
             residual_limitation = "Causal/therapeutic efficacy in jaw regeneration is not demonstrated; conclusions are nominations, not proof"),
  stringsAsFactors = FALSE
)
audit$audit_id <- sprintf("L%02d", seq_len(nrow(audit)))
audit <- audit[, c("audit_id", "risk", "evidence", "mitigation_taken", "residual_limitation")]
save_tsv(audit, file.path(PARAMS$dir_supp, "limitations_audit.tsv"))

cat(sprintf("\nStep11 audit done: %d limitation rows, %d datasets in provenance table.\n",
            nrow(audit), nrow(prov)))
cat("  -> supp/limitations_audit.tsv, supp/dataset_provenance.tsv\n")
dump_session()
