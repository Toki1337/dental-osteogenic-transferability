# R/tier1_atlas_GSE269255.R
# ============================================================================
# Process GSE269255 (mouse mandibular irradiation / ORNJ scRNA; Control/Day1/Day7)
# into a compartment-annotated Seurat object, so the transferability benchmark
# gains a THIRD in-vivo jaw-disease context (ORNJ) alongside human MRONJ
# (GSE303003) and mouse BRONJ (GSE295106). Same QC/normalization/clustering as the
# other atlases; compartments assigned by canonical mouse markers.
# Output: 06_scrna/GSE269255_seurat.rds
# ============================================================================
.libPaths(c("F:/Rlib", .libPaths()))
suppressWarnings(suppressMessages({ library(Seurat); library(Matrix); library(data.table) }))
setDTthreads(1); set.seed(1337); options(future.globals.maxSize = 6e9)
DIR <- "00_data/failure/GSE269255/extracted"
samps <- c("GSM8311031_Jaw-Control","GSM8311032_Jaw-Day1","GSM8311033_Jaw-Day7")

objs <- lapply(samps, function(s) {
  m <- ReadMtx(mtx = file.path(DIR, paste0(s, ".matrix.mtx.gz")),
               cells = file.path(DIR, paste0(s, ".barcodes.tsv.gz")),
               features = file.path(DIR, paste0(s, ".features.tsv.gz")), feature.column = 2)
  o <- CreateSeuratObject(m, project = s, min.cells = 3, min.features = 200)
  o$sample <- s; o$condition <- sub(".*Jaw-", "", s); o
})
obj <- merge(objs[[1]], objs[-1], add.cell.ids = sub("_.*","",samps))
obj <- JoinLayers(obj)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 20)
cat("cells after QC:", ncol(obj), "\n")
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)

# canonical mouse compartment markers
mk <- list(
  Mesenchymal_osteo = c("Col1a1","Col1a2","Pdgfra","Lum","Dcn","Runx2","Sp7","Bglap","Alpl","Sparc"),
  Osteoclast        = c("Ctsk","Acp5","Mmp9","Oscar","Dcstamp"),
  Myeloid           = c("Lyz2","Cd68","Itgam","Csf1r","C1qa","C1qb"),
  Tcell             = c("Cd3e","Cd3d","Cd3g","Cd8a","Cd4"),
  Bcell             = c("Cd79a","Cd79b","Ms4a1","Igkc"),
  Endothelial       = c("Pecam1","Cdh5","Flt1","Kdr"),
  Epithelial        = c("Epcam","Krt14","Krt5","Krt8"),
  Erythroid         = c("Hba-a1","Hbb-bs","Alas2"))
avg <- AverageExpression(obj, features = unique(unlist(mk)), group.by = "seurat_clusters",
                         assay = "RNA", layer = "data")$RNA
zz <- t(scale(t(log1p(as.matrix(avg)))))
clcomp <- sapply(colnames(zz), function(cl) {
  sc <- sapply(mk, function(g) { gg <- intersect(g, rownames(zz)); if (!length(gg)) return(-Inf); mean(zz[gg, cl]) })
  names(which.max(sc)) })
names(clcomp) <- sub("^g", "", names(clcomp))   # AverageExpression prefixes numeric cluster names with 'g'
obj$compartment <- unname(clcomp[as.character(obj$seurat_clusters)])
cat("=== cluster -> compartment ===\n"); print(clcomp)
cat("=== compartment cell counts ===\n"); print(table(obj$compartment))
saveRDS(obj, "06_scrna/GSE269255_seurat.rds")
cat("wrote 06_scrna/GSE269255_seurat.rds\n")
