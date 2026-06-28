# R/real_run_mronj_human_projection.R
# HUMAN MRONJ failure-side analysis (referee R4 main upgrade): project the
# pre-validated pan-dental-source L1 osteogenic "regeneration-competent" program
# onto GSE303003 (human: 4 MRONJ vs 4 odontogenic-cyst controls, 10x), and test
# whether the program is SUPPRESSED in the MRONJ MESENCHYMAL/osteoprogenitor
# compartment. Replaces the n=1v1 mouse cartoon with a replicated human contrast.
#
# Design guards (per review):
#  - control = cyst granuloma (inflammatory tissue, NOT healthy bone) -> project
#    ONLY within the mesenchymal/osteoprogenitor compartment to avoid immune dilution;
#  - PRE-REGISTERED whole-program sign test (768 + 52), NOT a curated IL6/IGFBP5 subset;
#  - we do NOT redo the dataset's own atlas / COL27A1 / osteoclast-fibroblast story.
.libPaths(c("F:/Rlib", .libPaths()))
suppressMessages({ library(Seurat); library(Matrix); library(data.table); library(harmony) })
set.seed(1)
D <- "00_data/failure/GSE303003/extracted"
OUT <- "06_scrna"; dir.create(OUT, showWarnings=FALSE)
fl <- function(...) { cat(sprintf("[%s] %s\n", format(Sys.time(),"%H:%M:%S"), paste0(...))); flush.console() }

# ---- samples ----
samples <- list(
  Cyst_1=c("GSM9115504_Cyst_1"), Cyst_2=c("GSM9115504_Cyst_2"),
  Cyst_3=c("GSM9115505_Cyst_3"), Cyst_4=c("GSM9115505_Cyst_4"),
  MRONJ_1=c("GSM9115504_MRONJ_1"), MRONJ_3=c("GSM9115504_MRONJ_3"),
  MRONJ_4=c("GSM9115505_MRONJ_4"), MRONJ_5=c("GSM9115505_MRONJ_5"))

read_one <- function(pfx){
  m <- readMM(file.path(D, paste0(pfx, ".matrix.mtx.gz")))
  ft <- fread(file.path(D, paste0(pfx, ".features.tsv.gz")), header=FALSE)
  bc <- fread(file.path(D, paste0(pfx, ".barcodes.tsv.gz")), header=FALSE)
  rownames(m) <- make.unique(ft$V2); colnames(m) <- bc$V1   # gene SYMBOL x cell
  m
}
objs <- list()
for(s in names(samples)){
  m <- read_one(samples[[s]][1])
  o <- CreateSeuratObject(m, project=s, min.cells=3, min.features=200)
  o$sample <- s; o$condition <- ifelse(grepl("MRONJ", s), "MRONJ", "Cyst")
  o[["percent.mt"]] <- PercentageFeatureSet(o, pattern="^MT-")
  o <- subset(o, subset = nFeature_RNA>200 & nFeature_RNA<7000 & percent.mt<20)
  objs[[s]] <- o
  fl(s, ": ", ncol(o), " cells after QC")
}
sc <- merge(objs[[1]], objs[-1], add.cell.ids=names(objs))
sc <- JoinLayers(sc)
fl("merged: ", ncol(sc), " cells x ", nrow(sc), " genes")

# ---- standard pipeline + Harmony integration by sample (checkpointed) ----
ckpt <- file.path(OUT,"GSE303003_clustered.rds")
if(file.exists(ckpt)){
  sc <- readRDS(ckpt); fl("loaded clustered checkpoint: ", ncol(sc), " cells, ",
                          length(unique(sc$seurat_clusters)), " clusters")
} else {
  sc <- NormalizeData(sc, verbose=FALSE)
  sc <- FindVariableFeatures(sc, nfeatures=2000, verbose=FALSE)
  sc <- ScaleData(sc, verbose=FALSE)
  sc <- RunPCA(sc, npcs=30, verbose=FALSE)
  sc <- RunHarmony(sc, group.by.vars="sample", verbose=FALSE)
  sc <- FindNeighbors(sc, reduction="harmony", dims=1:30, verbose=FALSE)
  sc <- FindClusters(sc, resolution=0.5, verbose=FALSE)
  sc <- RunUMAP(sc, reduction="harmony", dims=1:30, verbose=FALSE)
  saveRDS(sc, ckpt)
  fl("clusters: ", length(unique(sc$seurat_clusters)))
}

# ---- annotate compartments by canonical markers (mean scaled expr per cluster) ----
mk <- list(
  Mesenchymal_osteo = c("COL1A1","COL1A2","LUM","DCN","PDGFRA","PDGFRB","RUNX2","SP7","ALPL","IBSP","BGLAP"),
  Osteoclast        = c("CTSK","ACP5","MMP9","OCSTAMP","DCSTAMP"),
  Myeloid           = c("CD68","LYZ","CD14","FCGR3A","ITGAM"),
  Tcell             = c("CD3D","CD3E","CD2","IL7R"),
  Bcell_plasma      = c("MS4A1","CD79A","MZB1","IGHG1"),
  Endothelial       = c("PECAM1","VWF","CLDN5"),
  Epithelial        = c("KRT14","KRT5","EPCAM"))
present <- function(g) g[g %in% rownames(sc)]
avg <- AverageExpression(sc, group.by="seurat_clusters", assays="RNA", layer="data")$RNA
# per-cluster compartment score = mean (across markers) of each marker's z across clusters
scoreC <- sapply(mk, function(gs){ gs<-present(gs); if(!length(gs)) return(rep(0,ncol(avg)))
  rowMeans(scale(t(avg[gs,,drop=FALSE])), na.rm=TRUE) })   # rows = clusters, cols = compartments
rownames(scoreC) <- sub("^g","", colnames(avg))
cl_label <- colnames(scoreC)[apply(scoreC, 1, which.max)]
names(cl_label) <- rownames(scoreC)
sc$compartment <- unname(cl_label[as.character(sc$seurat_clusters)])
fl("compartment counts:"); print(table(sc$compartment))
ct <- as.data.frame.matrix(table(sc$sample, sc$compartment))
fwrite(data.frame(sample=rownames(ct), ct), file.path(OUT,"GSE303003_compartment_counts.tsv"), sep="\t")

# ---- L1 program scoring (AddModuleScore: up and down separately) ----
sig <- fread("03_rra_meta/regeneration_competent_signature.tsv")
core <- fread("04_modules_wgcna/core_program.tsv")
up768 <- present(sig$gene[sig$direction=="up"]); dn768 <- present(sig$gene[sig$direction=="down"])
upC   <- present(core$gene[core$direction=="up"]); dnC <- present(core$gene[core$direction=="down"])
sc <- AddModuleScore(sc, features=list(up768), name="prog_up_", seed=1)
sc <- AddModuleScore(sc, features=list(dn768), name="prog_dn_", seed=1)
sc <- AddModuleScore(sc, features=list(upC),   name="core_up_", seed=1)
sc <- AddModuleScore(sc, features=list(dnC),   name="core_dn_", seed=1)
sc$prog_net <- sc$prog_up_1 - sc$prog_dn_1     # regeneration-competent program score
sc$core_net <- sc$core_up_1 - sc$core_dn_1
# optional UCell (rank-based) if installed
if(requireNamespace("UCell", quietly=TRUE)){
  suppressMessages(library(UCell))
  sc <- AddModuleScore_UCell(sc, features=list(prog_up=up768, prog_dn=dn768))
  sc$prog_net_ucell <- sc$prog_up_UCell - sc$prog_dn_UCell
  fl("UCell scoring done")
}

# ---- KEY TEST: program suppressed in MRONJ within the MESENCHYMAL compartment? ----
mes <- sc$compartment=="Mesenchymal_osteo"
fl("mesenchymal/osteo cells: ", sum(mes), " (", round(100*mean(mes),1), "%)")
md <- data.frame(sample=sc$sample[mes], condition=sc$condition[mes],
                 prog_net=sc$prog_net[mes], core_net=sc$core_net[mes])
persamp <- aggregate(cbind(prog_net, core_net) ~ sample + condition, md, mean)
fwrite(persamp, file.path(OUT,"GSE303003_mes_program_persample.tsv"), sep="\t")
fl("=== per-sample mean program score in MESENCHYMAL compartment ===")
print(persamp)
wt_prog <- wilcox.test(prog_net ~ condition, persamp)
tt_prog <- t.test(prog_net ~ condition, persamp)
wt_core <- wilcox.test(core_net ~ condition, persamp)
mronj_lo <- mean(persamp$prog_net[persamp$condition=="MRONJ"]) < mean(persamp$prog_net[persamp$condition=="Cyst"])
fl("PROGRAM (768) MRONJ vs Cyst (mesenchymal): Wilcoxon p=", signif(wt_prog$p.value,3),
   "  t p=", signif(tt_prog$p.value,3), "  MRONJ_lower=", mronj_lo)
fl("CORE (52) MRONJ vs Cyst (mesenchymal): Wilcoxon p=", signif(wt_core$p.value,3))

# ---- pseudobulk DE within mesenchymal compartment (4v4, DESeq2) ----
suppressMessages(library(DESeq2))
mcells <- colnames(sc)[mes]
cnt <- GetAssayData(sc, layer="counts")[, mcells, drop=FALSE]
grp <- factor(sc$sample[mes])
pb <- t(rowsum(t(as.matrix(cnt)), grp))                 # gene x sample pseudobulk
cond <- ifelse(grepl("MRONJ", colnames(pb)), "MRONJ", "Cyst")
coldata <- data.frame(row.names=colnames(pb), condition=factor(cond, levels=c("Cyst","MRONJ")))
pb <- pb[rowSums(pb)>=10, ]
dds <- DESeqDataSetFromMatrix(round(pb), coldata, ~condition)
dds <- DESeq(dds, quiet=TRUE)
res <- as.data.frame(results(dds, contrast=c("condition","MRONJ","Cyst")))
res$gene <- rownames(res)
fwrite(res, file.path(OUT,"GSE303003_mes_pseudobulk_DE.tsv"), sep="\t")
fl("pseudobulk DE genes tested: ", nrow(res), "; FDR<0.05: ", sum(res$padj<0.05, na.rm=TRUE))

# ---- PRE-REGISTERED whole-program sign test (no curated subset) ----
# Hypothesis: in MRONJ mesenchymal cells the program is SUPPRESSED, i.e. signature
# UP-genes go DOWN (log2FC<0) and DOWN-genes go UP (log2FC>0) vs cyst control.
sign_test <- function(genes_up, genes_dn, label){
  ru <- res[res$gene %in% genes_up & !is.na(res$log2FoldChange),]
  rd <- res[res$gene %in% genes_dn & !is.na(res$log2FoldChange),]
  suppressed <- c(ru$log2FoldChange < 0, rd$log2FoldChange > 0)   # TRUE = consistent with suppression
  k <- sum(suppressed); n <- length(suppressed)
  bt <- binom.test(k, n, 0.5, alternative="greater")
  fl(label, ": ", k, "/", n, " (", round(100*k/n,1), "%) genes move in the SUPPRESSED direction; binom p(>0.5)=", signif(bt$p.value,3))
  data.frame(set=label, k=k, n=n, frac=round(k/n,3), binom_p_greater=signif(bt$p.value,3))
}
st <- rbind(
  sign_test(up768, dn768, "signature_768"),
  sign_test(upC,   dnC,   "core_52"))
fwrite(st, file.path(OUT,"GSE303003_program_signtest.tsv"), sep="\t")

# ---- compartment proportion test (4v4, with variance) ----
prop <- prop.table(table(sc$sample, sc$compartment), 1)
pdf_ <- as.data.frame.matrix(prop); pdf_$sample <- rownames(pdf_)
pdf_$condition <- ifelse(grepl("MRONJ", pdf_$sample), "MRONJ", "Cyst")
proptests <- lapply(setdiff(colnames(as.data.frame.matrix(prop)),""), function(cp){
  v <- pdf_[[cp]]; if(is.null(v)) return(NULL)
  p <- tryCatch(wilcox.test(v ~ pdf_$condition)$p.value, error=function(e) NA)
  data.frame(compartment=cp, Cyst_mean=mean(v[pdf_$condition=="Cyst"]),
             MRONJ_mean=mean(v[pdf_$condition=="MRONJ"]), wilcox_p=signif(p,3))
})
proptab <- do.call(rbind, proptests)
fwrite(proptab, file.path(OUT,"GSE303003_compartment_proptest.tsv"), sep="\t")
fl("=== compartment proportions MRONJ vs Cyst ==="); print(proptab)

# ---- per-compartment program score (localization: where does the program live?) ----
loc <- aggregate(prog_net ~ compartment, data.frame(prog_net=sc$prog_net, compartment=sc$compartment), mean)
loc <- loc[order(-loc$prog_net),]
fwrite(loc, file.path(OUT,"GSE303003_program_by_compartment.tsv"), sep="\t")
fl("=== program score by compartment (localization) ==="); print(loc)

saveRDS(sc, file.path(OUT,"GSE303003_seurat.rds"))
fl("DONE human MRONJ projection.")
