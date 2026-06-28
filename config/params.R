# config/params.R — central tunable parameters for the whole pipeline
# Sourced by every Step script. Keep ALL thresholds here for reproducibility.

PARAMS <- list(

  ## ---- paths ----
  root        = normalizePath(".", mustWork = FALSE),
  dir_data    = "00_data",
  dir_curation= "01_curation",
  dir_de      = "02_per_dataset_DE",
  dir_rra     = "03_rra_meta",
  dir_wgcna   = "04_modules_wgcna",
  dir_proj    = "05_projection",
  dir_scrna   = "06_scrna",
  dir_mr      = "07_mr_coloc",
  dir_cmap    = "08_cmap",
  dir_integr  = "09_integration",
  dir_fig     = "figures",
  dir_supp    = "supp",

  ## ---- differential expression (Step1) ----
  de_fdr        = 0.05,     # adjusted p cutoff
  de_lfc        = 1.0,      # |log2FC| cutoff for "DEG" calls (RRA uses full ranking, not cutoff)
  min_count     = 10,       # DESeq2 prefilter: rowSums >= min_count
  array_logspace= TRUE,     # assume arrays already log2; auto-checked in utils

  ## ---- RRA meta (Step2) ----
  rra_method    = "RRA",    # RobustRankAggreg method
  rra_topN      = 0.05,     # fraction of top genes contributed per dataset to RRA
  rra_fdr       = 0.05,     # RRA score BH cutoff for the meta-signature
  rra_min_datasets = 3,     # gene must rank in >= this many datasets to be eligible
  loo_stability_min = 0.7,  # leave-one-dataset-out: keep hub if stable in >=70% of LOO runs

  ## ---- WGCNA (Step3) ----
  wgcna_power_candidates = 1:20,
  wgcna_rsq_cut  = 0.85,    # scale-free topology R^2 target
  wgcna_minModuleSize = 30,
  wgcna_mergeCutHeight = 0.25,
  wgcna_deepSplit = 2,
  ppi_string_score = 700,   # STRING high-confidence (0-1000)
  hub_top_k      = 30,      # cytoHubba top-k hubs

  ## ---- ssGSEA / projection (Step4-5) ----
  ssgsea_method  = "ssgsea",# GSVA::gsva
  proj_perm      = 1000,    # permutations for module score significance
  homolog_db     = "babelgene", # mouse<->human ortholog mapping

  ## ---- scRNA (Step6) ----
  sc_min_genes   = 200,
  sc_max_mt_pct  = 20,
  sc_n_pcs       = 30,
  sc_cluster_res = 0.5,
  cellchat_db    = "human", # switch per-object

  ## ---- MR + coloc (Step7) ----
  mr_instrument_p = 5e-8,   # genome-wide for primary; relax to 5e-6 sensitivity
  mr_instrument_p_relax = 5e-6,
  mr_clump_r2    = 0.001,
  mr_clump_kb    = 10000,
  mr_methods     = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"),
  mr_fdr         = 0.05,
  coloc_pp4_min  = 0.75,    # PP.H4 posterior for colocalization
  smr_heidi_p    = 0.01,    # HEIDI p > 0.01 => not rejected (no pleiotropy)

  ## ---- CMap / repurposing (Step9) ----
  cmap_neg_score = -90,     # connectivity score threshold (clue.io tau)
  cmap_min_celllines = 3,   # require consistency across >= cell lines

  ## ---- integration scoring (Step10) ----
  # transparent weights for the multi-evidence priority score (sum to 1)
  weights = c(
    rra_strength      = 0.20,  # strength/robustness in pan-dental RRA program
    failure_specific  = 0.20,  # dysregulation magnitude in failure-specific program
    scrna_support     = 0.10,  # single-cell corroboration
    mr_causal         = 0.25,  # genetic causal support (highest weight)
    druggability      = 0.10,  # Finan tier
    jaw_context       = 0.05,  # jaw position-specific relevance
    cmap_connectivity = 0.10   # pharmacological reverse-connectivity (coda)
  ),

  seed = 1337
)

stopifnot(abs(sum(PARAMS$weights) - 1) < 1e-9)
set.seed(PARAMS$seed)
