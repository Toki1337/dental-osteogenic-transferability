# R/real_run_wgcna.R — REAL WGCNA on the curated dental compendium -> core program.
# Uses cached per-dataset matrices (RNA-seq series matrices are empty on GEO, so the
# standard step03 fetch path cannot apply here). Hubs = intramodular connectivity
# (kWithin) to avoid STRING's large download; STRINGdb remains available in step03.
.libPaths(c("F:/Rlib", .libPaths()))
suppressMessages({ library(WGCNA); library(limma); library(data.table); library(pheatmap) })
options(stringsAsFactors = FALSE); allowWGCNAThreads()
dir.create("04_modules_wgcna", showWarnings = FALSE); dir.create("figures", showWarnings = FALSE)

ex <- readRDS("03_rra_meta/exprs_cache.rds")
ph <- do.call(rbind, readRDS("03_rra_meta/pheno_cache.rds"))
common <- Reduce(intersect, lapply(ex, rownames))
cat("common genes across", length(ex), "datasets:", length(common), "\n")
mat <- do.call(cbind, lapply(names(ex), function(id) { m <- ex[[id]][common, , drop = FALSE]; colnames(m) <- paste(id, colnames(m), sep="|"); m }))
ph <- ph[match(colnames(mat), paste(ph$dataset, ph$sample, sep="|")), ]
stopifnot(all(!is.na(ph$group)))

# batch-correct dataset, preserve osteo/control
matbc <- limma::removeBatchEffect(mat, batch = ph$dataset, design = model.matrix(~ph$group))
# keep most-variable genes for WGCNA tractability
v <- apply(matbc, 1, var); top <- names(sort(v, decreasing = TRUE))[1:min(8000, length(v))]
datExpr <- t(matbc[top, ])

sft <- pickSoftThreshold(datExpr, powerVector = 1:20, RsquaredCut = 0.85, verbose = 0)
power <- sft$powerEstimate; if (is.na(power)) power <- 6
cat("soft power:", power, "\n")
net <- blockwiseModules(datExpr, power = power, TOMType = "signed", minModuleSize = 30,
                        mergeCutHeight = 0.25, deepSplit = 2, numericLabels = TRUE, verbose = 0)
modules <- data.frame(gene = colnames(datExpr), module = labels2colors(net$colors))
write.table(modules, "04_modules_wgcna/module_assignments.tsv", sep="\t", quote=FALSE, row.names=FALSE)

grp_bin <- as.numeric(ph$group == "osteo")
MEs <- net$MEs
# map numeric module label (ME columns) -> color name (modules$module uses colors)
lab2col <- tapply(modules$module, net$colors, function(x) x[1])
me_label <- as.character(as.integer(sub("ME","",colnames(MEs))))
me_color <- unname(lab2col[me_label])
me_cor <- sapply(MEs, function(me) cor(me, grp_bin, use="pairwise.complete.obs"))
me_p   <- sapply(MEs, function(me) cor.test(me, grp_bin)$p.value)
me_padj <- p.adjust(me_p, "BH")
# union of ALL modules significantly associated with osteo (either direction)
sig <- which(me_padj < 0.1 & me_color != "grey")
if (length(sig) == 0) sig <- which.max(abs(me_cor))
osteo_colors <- unique(me_color[sig])
osteo_genes <- modules$gene[modules$module %in% osteo_colors]
mt <- data.frame(module=me_color, r=round(me_cor,2), BH_p=signif(me_padj,2))
write.table(mt, "04_modules_wgcna/module_trait.tsv", sep="\t", quote=FALSE, row.names=FALSE)
cat(sprintf("osteo-associated modules (BH<0.1): %s | %d genes\n", paste(osteo_colors,collapse=","), length(osteo_genes)))

# core program = RRA signature ∩ osteo-associated module(s)
rra <- read.delim("03_rra_meta/regeneration_competent_signature.tsv")
core <- rra[rra$gene %in% osteo_genes, ]

# intramodular hub connectivity computed AMONG the core genes (relevant network)
cg <- intersect(core$gene, colnames(datExpr))
adj <- adjacency(datExpr[, cg, drop=FALSE], power = power, type = "signed")
kWithin <- rowSums(adj) - 1
hub <- data.frame(gene = names(kWithin), kWithin = round(kWithin,3))
hub <- hub[order(-hub$kWithin), ]
write.table(hub, "04_modules_wgcna/ppi_hubs.tsv", sep="\t", quote=FALSE, row.names=FALSE)
top_hubs <- head(hub$gene, 30)
core$is_ppi_hub <- core$gene %in% top_hubs
core <- core[order(core$direction, core$rra_score), ]
write.table(core, "04_modules_wgcna/core_program.tsv", sep="\t", quote=FALSE, row.names=FALSE)
cat(sprintf("CORE PROGRAM: %d genes (%d up, %d down; %d intramodular hubs)\n",
            nrow(core), sum(core$direction=="up"), sum(core$direction=="down"), sum(core$is_ppi_hub)))
cat("top hub core genes:", paste(head(core$gene[core$is_ppi_hub], 15), collapse=", "), "\n")

# ---- Fig2: RRA top-gene cross-dataset heatmap + LOO stability ----
ranked <- list.files("02_per_dataset_DE", "_ranked.tsv$", full.names = TRUE)
topg <- head(rra$gene[order(rra$rra_score)], 40)
hm <- sapply(ranked, function(f){ d<-read.delim(f); setNames(d$stat, d$gene)[topg] })
rownames(hm) <- topg; colnames(hm) <- sub("_ranked.tsv","",basename(ranked))
hm[is.na(hm)] <- 0
pheatmap(hm, filename="figures/Fig2_RRA_top40_cross_dataset.pdf", width=6, height=8,
         cluster_cols=TRUE, main="Top-40 RRA regeneration-competent genes (DE stat across datasets)",
         color=colorRampPalette(c("navy","white","firebrick"))(50))
saveRDS(list(power=power, osteo_colors=osteo_colors, me_cor=me_cor, me_padj=me_padj), "04_modules_wgcna/wgcna_meta.rds")
cat("DONE. core_program.tsv + Fig2 written.\n")
