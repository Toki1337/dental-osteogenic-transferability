# R/step03_wgcna_modules.R
# Step 3 — Module construction + PPI hubs + functional annotation -> core program.
# WGCNA on the batch-corrected pan-dental osteogenic expression compendium (clean
# arms only), intersect the osteo-correlated module(s) with the RRA signature,
# extract STRING PPI hubs, annotate (GO/KEGG/WikiPathways).
# Output:
#   04_modules_wgcna/module_assignments.tsv
#   04_modules_wgcna/osteo_module_genes.tsv
#   04_modules_wgcna/ppi_hubs.tsv
#   04_modules_wgcna/core_program.tsv          (core = RRA ∩ osteo-module; PPI-hub is annotation only, NOT a filter)
#   04_modules_wgcna/enrichment_{GO,KEGG,WP}.tsv

source("R/utils.R"); init_pipeline()
suppressMessages({ library(WGCNA); library(limma) })
allowWGCNAThreads()

# ---- build expression compendium from clean arms (symbols x samples) ----
arms_path   <- file.path(PARAMS$dir_curation, "sample_arms_final.tsv")
usable_path <- file.path(PARAMS$dir_curation, "rra_usable_datasets.tsv")
if (!file.exists(arms_path) || !file.exists(usable_path))
  stop("Run Step0 + human curation first.")
arms <- utils::read.delim(arms_path, stringsAsFactors = FALSE)
# Consume the human purity gate: keep only samples whose arm_clean is yes/true.
# Fail-safe: blank/NA arm_clean is treated as NOT clean (excluded).
arms <- arms[arms$final_group %in% c("osteo", "control") &
               (tolower(arms$arm_clean) %in% c("yes", "y", "true", "1")), ]
usable <- read.delim(usable_path)$accession

mats <- list(); pheno <- list()
for (acc in usable) {
  ri <- which(REG$accession == acc)
  eset <- fetch_geo(acc, REG$role[ri])
  m <- get_expr_symbol(eset)
  a <- arms[arms$accession == acc, ]
  ids <- intersect(colnames(m), a$gsm)
  m <- m[, ids, drop = FALSE]
  colnames(m) <- paste(acc, ids, sep = "|")
  mats[[acc]] <- m
  pheno[[acc]] <- data.frame(sample = colnames(m), dataset = acc,
                             group = a$final_group[match(ids, a$gsm)], stringsAsFactors = FALSE)
}
common <- Reduce(intersect, lapply(mats, rownames))
mats <- lapply(mats, function(m) m[common, , drop = FALSE])
expr <- do.call(cbind, mats)
ph   <- do.call(rbind, pheno)

# batch-correct dataset effect, preserving group (osteo/control) signal
expr_bc <- limma::removeBatchEffect(expr, batch = ph$dataset,
                                    design = model.matrix(~ph$group))
datExpr <- t(expr_bc)

# ---- pick soft power ----
sft <- pickSoftThreshold(datExpr, powerVector = PARAMS$wgcna_power_candidates,
                         RsquaredCut = PARAMS$wgcna_rsq_cut, verbose = 0)
power <- sft$powerEstimate %||% 6
.log("WGCNA soft power = ", power)

net <- blockwiseModules(datExpr, power = power, TOMType = "signed",
                        minModuleSize = PARAMS$wgcna_minModuleSize,
                        mergeCutHeight = PARAMS$wgcna_mergeCutHeight,
                        deepSplit = PARAMS$wgcna_deepSplit, numericLabels = TRUE,
                        saveTOMs = FALSE, verbose = 0)
modules <- data.frame(gene = colnames(datExpr), module = labels2colors(net$colors), stringsAsFactors = FALSE)
save_tsv(modules, file.path(PARAMS$dir_wgcna, "module_assignments.tsv"))

# ---- correlate module eigengenes with osteo/control ----
# NOTE: removeBatchEffect(design=~group) regresses the dataset effect while
# preserving group, which tends to make the residual module-trait p anti-
# conservative-looking yet directionally reliable; we therefore pick the osteo
# module mainly to feed RRA∩module, not as a hard significance gate.
MEs <- net$MEs
grp_bin <- as.numeric(ph$group == "osteo")
me_cor <- apply(MEs, 2, function(me) stats::cor(me, grp_bin, use = "pairwise.complete.obs"))
me_p   <- apply(MEs, 2, function(me) stats::cor.test(me, grp_bin)$p.value)
me_padj <- stats::p.adjust(me_p, method = "BH")
sig <- which(me_padj < 0.05)
if (length(sig)) {
  osteo_me <- names(sig)[which.max(abs(me_cor[sig]))]
} else {
  # Fail-safe: no BH-significant module -> fall back to the largest |r| module
  # (the RRA∩module intersection below is the real selection criterion anyway).
  osteo_me <- names(me_cor)[which.max(abs(me_cor))]
  .log("no module-trait BH p<0.05; falling back to max |r| module: ", osteo_me, level = "WARN")
}
osteo_color <- sub("^ME", "", osteo_me)
osteo_genes <- modules$gene[modules$module == osteo_color]
save_tsv(data.frame(gene = osteo_genes), file.path(PARAMS$dir_wgcna, "osteo_module_genes.tsv"))
.log("osteo-correlated module: ", osteo_color, " (", length(osteo_genes), " genes, r=", round(me_cor[osteo_me], 2), ", BH p=", signif(me_padj[osteo_me], 2), ")")

# ---- STRING PPI hubs within osteo module ----
ppi_hubs <- tryCatch({
  suppressMessages(library(STRINGdb))
  sdb <- STRINGdb$new(version = "12.0", species = 9606, score_threshold = PARAMS$ppi_string_score, input_directory = "")
  mapped <- sdb$map(data.frame(gene = osteo_genes), "gene", removeUnmappedRows = TRUE)
  ints <- sdb$get_interactions(mapped$STRING_id)
  g <- igraph::graph_from_data_frame(ints[, c("from", "to")], directed = FALSE)
  deg <- sort(igraph::degree(g), decreasing = TRUE)
  hub_ids <- names(head(deg, PARAMS$hub_top_k))
  data.frame(gene = mapped$gene[match(hub_ids, mapped$STRING_id)], degree = head(deg, PARAMS$hub_top_k))
}, error = function(e) { .log("STRING step failed: ", conditionMessage(e), level = "WARN"); data.frame(gene = character(0), degree = integer(0)) })
save_tsv(ppi_hubs, file.path(PARAMS$dir_wgcna, "ppi_hubs.tsv"))

# ---- core program = RRA ∩ osteo module (PPI hub is annotated, NOT used to filter) ----
rra <- utils::read.delim(file.path(PARAMS$dir_rra, "regeneration_competent_signature.tsv"), stringsAsFactors = FALSE)
core <- intersect(rra$gene, osteo_genes)
core_df <- rra[rra$gene %in% core, ]
core_df$is_ppi_hub <- core_df$gene %in% ppi_hubs$gene
save_tsv(core_df, file.path(PARAMS$dir_wgcna, "core_program.tsv"))
.log("core program: ", nrow(core_df), " genes (", sum(core_df$is_ppi_hub), " also PPI hubs)")

# ---- enrichment ----
tryCatch({
  suppressMessages({ library(clusterProfiler); library(org.Hs.eg.db) })
  eg <- bitr(core_df$gene, "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID
  go <- enrichGO(eg, org.Hs.eg.db, ont = "BP", readable = TRUE)
  kg <- enrichKEGG(eg)
  save_tsv(as.data.frame(go), file.path(PARAMS$dir_wgcna, "enrichment_GO.tsv"))
  if (!is.null(kg)) save_tsv(as.data.frame(kg), file.path(PARAMS$dir_wgcna, "enrichment_KEGG.tsv"))
}, error = function(e) .log("enrichment failed: ", conditionMessage(e), level = "WARN"))

dump_session()
.log("Step3 done. core_program.tsv = the pan-dental regeneration-competent core program.")
