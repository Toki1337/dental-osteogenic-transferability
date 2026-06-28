# R/real_run_step56_scrna.R — REAL Step5/6 failure-side single cell.
# GSE295106 mandibular marrow Control vs BRONJ (10x). Cluster+annotate, compartment
# shifts, pseudobulk projection of the regeneration-competent program (failure
# suppresses it?), per-celltype DE (scrna_dysregulated_genes.tsv for Step10), and
# CellChat Control-vs-BRONJ communication rewiring (RANKL/TNF/SPP1).
.libPaths(c("F:/Rlib", .libPaths()))
suppressMessages({ library(Seurat); library(Matrix); library(ggplot2); library(data.table) })
set.seed(1337)
D <- "00_data/failure/GSE295106"; dir.create("06_scrna", showWarnings=FALSE); dir.create("05_projection", showWarnings=FALSE)
lg <- function(...) cat(sprintf("[%s] ", format(Sys.time(),"%H:%M:%S")), ..., "\n")

read_sample <- function(prefix, cond){
  m <- ReadMtx(mtx=file.path(D,paste0(prefix,"_matrix.mtx.gz")),
               cells=file.path(D,paste0(prefix,"_barcodes.tsv.gz")),
               features=file.path(D,paste0(prefix,"_features.tsv.gz")))
  o <- CreateSeuratObject(m, project=cond, min.cells=3, min.features=200)
  o$condition <- cond; o$sample <- prefix; o
}
lg("reading 10x (total Control + total BRONJ)...")
ctrl <- read_sample("GSM8942338_Control__total","Control")
bronj<- read_sample("GSM8942339_BRONJ__total","BRONJ")
obj <- merge(ctrl, bronj); obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern="^mt-")
obj <- subset(obj, nFeature_RNA>=200 & nFeature_RNA<7000 & percent.mt<20)
lg("cells after QC:", ncol(obj), " (", paste(names(table(obj$condition)),table(obj$condition),collapse=" "),")")

obj <- NormalizeData(obj, verbose=FALSE) |> FindVariableFeatures(nfeatures=2000, verbose=FALSE) |>
       ScaleData(verbose=FALSE) |> RunPCA(npcs=30, verbose=FALSE)
red <- "pca"
if (requireNamespace("harmony", quietly=TRUE)) {
  obj <- harmony::RunHarmony(obj, group.by.vars="sample", verbose=FALSE); red <- "harmony" }
obj <- FindNeighbors(obj, reduction=red, dims=1:30, verbose=FALSE) |>
       FindClusters(resolution=0.5, verbose=FALSE) |> RunUMAP(reduction=red, dims=1:30, verbose=FALSE)

# annotate (mouse markers)
MK <- list(Osteo_MSC=c("Runx2","Sp7","Col1a1","Lepr","Bglap","Pdgfra"),
           Osteoclast=c("Ctsk","Acp5","Mmp9","Oscar"), Macrophage=c("Cd68","Lyz2","Adgre1","Itgam"),
           Neutrophil=c("S100a8","S100a9","Ly6g"), Endothelial=c("Pecam1","Cdh5","Emcn"),
           Tcell=c("Cd3e","Cd8a","Cd4"), Bcell=c("Cd79a","Ms4a1"), Erythroid=c("Hba-a1","Hbb-bs"))
ad <- AggregateExpression(obj, features=unique(unlist(MK)), group.by="seurat_clusters", assays="RNA")$RNA
adz <- t(scale(t(as.matrix(ad))))
scoreCT <- sapply(names(MK), function(ct){ g<-intersect(MK[[ct]],rownames(adz)); if(!length(g)) return(rep(0,ncol(adz))); colMeans(adz[g,,drop=FALSE]) })
ctlab <- apply(scoreCT,1,function(r) names(MK)[which.max(r)])
names(ctlab) <- sub("^g","",rownames(scoreCT))   # AggregateExpression prefixes numeric clusters with 'g'
obj$cell_type <- unname(ctlab[as.character(obj$seurat_clusters)])
obj$cell_type[is.na(obj$cell_type)] <- "Unassigned"
lg("cell types:", paste(names(table(obj$cell_type)), table(obj$cell_type), collapse=" "))

ggsave("figures/Fig5a_scRNA_umap.png",
  (DimPlot(obj,group.by="cell_type",label=TRUE,repel=TRUE)|DimPlot(obj,group.by="condition")), width=12,height=5, dpi=130)

# compartment proportions Control vs BRONJ
pt <- as.data.frame.matrix(table(obj$cell_type, obj$condition))
pt <- sweep(pt,2,colSums(pt),"/"); pt$cell_type<-rownames(pt)
write.table(pt, "06_scrna/sc_celltype_props.tsv", sep="\t", quote=FALSE, row.names=FALSE)
lg("compartment proportions written")

# pseudobulk per condition -> project regeneration-competent signature
pb <- AggregateExpression(obj, group.by="condition", assays="RNA", slot="counts")$RNA
pb <- as.matrix(pb); pb <- log2(t(t(pb)/colSums(pb))*1e6 + 1)            # logCPM
rownames(pb) <- toupper(rownames(pb))                                    # mouse->human ortholog (uppercase)
sig <- read.delim("03_rra_meta/regeneration_competent_signature.tsv")
up <- intersect(sig$gene[sig$direction=="up"], rownames(pb)); dn <- intersect(sig$gene[sig$direction=="down"], rownames(pb))
z <- t(scale(t(pb))); z[is.na(z)]<-0
relctl <- mean(z[up,"Control"])-mean(z[dn,"Control"]); relbr <- mean(z[up,"BRONJ"])-mean(z[dn,"BRONJ"])
# delta per signature gene (BRONJ - Control); failure-specific = signature-up genes DOWN in BRONJ
delta <- pb[,"BRONJ"]-pb[,"Control"]
fs <- data.frame(gene=names(delta), delta_BRONJ_vs_Control=round(delta,3))
fs$signature_dir <- sig$direction[match(fs$gene, sig$gene)]
fs <- fs[!is.na(fs$signature_dir),]
fs$failure_specific <- (fs$signature_dir=="up" & fs$delta_BRONJ_vs_Control<0) | (fs$signature_dir=="down" & fs$delta_BRONJ_vs_Control>0)
fs$abs_effect <- abs(fs$delta_BRONJ_vs_Control); fs$failure_direction <- ifelse(fs$delta_BRONJ_vs_Control<0,"down","up")
fs <- fs[order(-fs$failure_specific,-fs$abs_effect),]
write.table(fs, "05_projection/failure_specific_program.tsv", sep="\t", quote=FALSE, row.names=FALSE)
lg(sprintf("reg-competent score Control=%.3f BRONJ=%.3f ; failure-specific genes=%d/%d", relctl, relbr, sum(fs$failure_specific), nrow(fs)))

# per-celltype DE (BRONJ vs Control) -> scrna_dysregulated_genes.tsv (human symbols)
Idents(obj) <- "cell_type"; dys <- list()
for (ct in unique(obj$cell_type)) {
  sub <- subset(obj, cell_type==ct)
  if (sum(sub$condition=="BRONJ")<20 || sum(sub$condition=="Control")<20) next
  mk <- tryCatch(FindMarkers(sub, ident.1="BRONJ", ident.2="Control", group.by="condition", logfc.threshold=0.25, min.pct=0.1), error=function(e) NULL)
  if (is.null(mk)||!nrow(mk)) next
  mk$gene <- toupper(rownames(mk)); mk$cell_type<-ct
  dys[[ct]] <- data.frame(gene=mk$gene, cell_type=ct, is_dysregulated=mk$p_val_adj<0.05,
                          abs_avg_log2FC=abs(mk$avg_log2FC), dataset="GSE295106")
}
dysdf <- if(length(dys)) do.call(rbind, dys) else data.frame(gene=character(0),cell_type=character(0),is_dysregulated=logical(0),abs_avg_log2FC=numeric(0),dataset=character(0))
write.table(dysdf, "06_scrna/scrna_dysregulated_genes.tsv", sep="\t", quote=FALSE, row.names=FALSE)
saveRDS(obj, "06_scrna/GSE295106_seurat.rds")
lg("per-celltype DE written:", nrow(dysdf), "rows")

# CellChat Control vs BRONJ
ok <- tryCatch({
  suppressMessages(library(CellChat))
  runcc <- function(cells){
    sub <- subset(obj, cells=cells); din <- GetAssayData(sub, layer="data")
    rownames(din) <- toupper(rownames(din))                  # mouse->human for CellChatDB.human
    din <- din[!duplicated(rownames(din)),]
    cc <- createCellChat(din, meta=data.frame(group=sub$cell_type, row.names=colnames(sub)), group.by="group")
    cc@DB <- CellChatDB.human; cc <- subsetData(cc); cc <- identifyOverExpressedGenes(cc)
    cc <- identifyOverExpressedInteractions(cc); cc <- computeCommunProb(cc, type="triMean")
    cc <- filterCommunication(cc, min.cells=10); cc <- computeCommunProbPathway(cc); aggregateNet(cc)
  }
  cells_c <- colnames(obj)[obj$condition=="Control"]; cells_b <- colnames(obj)[obj$condition=="BRONJ"]
  ccc <- runcc(cells_c); ccb <- runcc(cells_b)
  saveRDS(list(Control=ccc, BRONJ=ccb), "06_scrna/cellchat_objs.rds")
  pw <- function(cc){ if(is.null(cc@netP$pathways)||!length(cc@netP$pathways)) return(setNames(numeric(0),character(0)))
    setNames(sapply(seq_along(cc@netP$pathways), function(k) sum(cc@netP$prob[,,k],na.rm=TRUE)), cc@netP$pathways) }
  sc<-pw(ccc); sb<-pw(ccb); allp<-union(names(sc),names(sb))
  pwd <- data.frame(pathway=allp, strength_Control=as.numeric(sc[allp]), strength_BRONJ=as.numeric(sb[allp]))
  pwd[is.na(pwd)]<-0; pwd$delta_BRONJ_minus_Control<-pwd$strength_BRONJ-pwd$strength_Control
  pwd$is_focus <- pwd$pathway %in% c("RANKL","TNF","SPP1","TGFb","BMP","WNT","PTN","PDGF")
  pwd <- pwd[order(-abs(pwd$delta_BRONJ_minus_Control)),]
  write.table(pwd, "06_scrna/cellchat_pathway_diff.tsv", sep="\t", quote=FALSE, row.names=FALSE)
  lg("CellChat pathway diff written:", nrow(pwd), "pathways"); TRUE
}, error=function(e){ lg("CellChat failed:", conditionMessage(e)); FALSE })
lg("Step5/6 DONE. CellChat ok=", ok)
