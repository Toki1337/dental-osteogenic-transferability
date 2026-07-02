# R/tier3_add_homeostatic_GSE316924.R
# ============================================================================
# Item B (reviewer request): test whether the in-vitro -> in-vivo non-correspondence
# holds in a HOMEOSTATIC (non-diseased) orofacial-bone reference, not only in the
# three DISEASE atlases (MRONJ/BRONJ/ORNJ).
#
# Dataset (verified to exist, 2026):
#   GSE316924 - "8-week-old mouse mandible scRNA-seq reveals bone marrow
#   microenvironment" (Mus musculus, scRNA-seq). Single sample / ~3568 cells:
#   use for REFERENCE-LOCALISATION ONLY (positive-control-validated Transferability
#   Score), NOT for differential statistics.
#
# This is an ADDITIVE robustness panel: it is scored the same way as the three
# disease atlases but is reported separately and is NOT merged into the 4/4-vs-1/9
# provenance Fisher test.
#
# PREREQUISITES / GATE (checklist B0-B1):
#   1. Download GSE316924 processed matrix + cell annotations from GEO.
#   2. The dataset MUST contain an annotatable osteolineage compartment
#      (EMP/LMP/pre-osteoblast/osteoblast/MSC_OLC). If canonical markers
#      (RUNX2/SP7/COL1A1/BGLAP/ALPL) do NOT localise there (positive control
#      fails), DROP this atlas - it is not method-validated.
#   Edit RAW_PATH and the annotation block below to match the deposited files.
#   Rscript R/tier3_add_homeostatic_GSE316924.R
# ============================================================================
.libPaths(c(Sys.getenv("R_TRANSFERAUDIT_LIB", unset = .libPaths()[1]), .libPaths()))
suppressWarnings(suppressMessages({ library(Seurat); library(Matrix); library(data.table) }))
setDTthreads(1); source("R/utils.R"); source("config/params.R")
OUT <- "10_transferability"; set.seed(PARAMS$seed); B <- 2000L
RAW_PATH <- "00_data/homeostatic/GSE316924"      # <-- put the downloaded 10x/matrix here

## ---- 1. load + minimal QC + cluster (edit to the deposited format) ----------
# If GEO provides a processed .rds/.h5ad with cell-type labels, load it directly
# and skip clustering; otherwise the block below builds one from a 10x matrix.
obj <- Read10X(RAW_PATH) |> CreateSeuratObject(min.cells = 3, min.features = 200)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 20)
obj <- NormalizeData(obj) |> FindVariableFeatures(nfeatures = 2000) |> ScaleData() |>
       RunPCA(npcs = 30) |> FindNeighbors(dims = 1:30) |> FindClusters(resolution = 0.5) |>
       RunUMAP(dims = 1:30)

## ---- 2. annotate an "Osteo_MSC" osteolineage compartment --------------------
# Prefer the authors' deposited labels if present in meta.data; else annotate by
# canonical markers. Mouse symbols are Title-case (Runx2, Sp7, Col1a1, Bglap, Alpl).
mk <- list(
  Osteo_MSC   = c("Runx2","Sp7","Col1a1","Col1a2","Bglap","Alpl","Ibsp","Spp1","Sparc","Lepr","Cxcl12"),
  Osteoclast  = c("Ctsk","Acp5","Mmp9"),
  Myeloid     = c("Lyz2","Cd68","Csf1r"),
  Tcell       = c("Cd3e","Cd3d"), Bcell = c("Cd79a","Ms4a1"),
  Endothelial = c("Pecam1","Cdh5"), Erythroid = c("Hba-a1","Hbb-bs"))
obj <- AddModuleScore(obj, features = mk, name = "cs")
cs <- obj@meta.data[, grep("^cs[0-9]+$", colnames(obj@meta.data))]; colnames(cs) <- names(mk)
cl_lab <- tapply(seq_len(ncol(obj)), obj$seurat_clusters, function(ix)
  names(mk)[which.max(colMeans(cs[ix, , drop = FALSE]))])
obj$compartment <- cl_lab[as.character(obj$seurat_clusters)]
stopifnot("Osteo_MSC" %in% obj$compartment)

## ---- 3. score the panel with the calibrated Transferability Score -----------
build_M <- function(dat, comp) { comp <- as.character(comp); cl <- sort(unique(comp))
  mu <- Matrix::rowMeans(dat); sdg <- sqrt(pmax(Matrix::rowMeans(dat*dat) - mu^2, 0))
  M <- (sapply(cl, function(c) Matrix::rowMeans(dat[, comp==c, drop=FALSE])) - mu)/sdg
  M[is.finite(rowSums(M)) & sdg > 0, , drop = FALSE] }
score_comps <- function(M, up, dn){ up<-intersect(up,rownames(M)); dn<-intersect(dn,rownames(M))
  if(length(up)<3) return(NULL); s<-colMeans(M[up,,drop=FALSE]); if(length(dn)>=3) s<-s-colMeans(M[dn,,drop=FALSE]); s }

sigs <- readRDS(file.path(OUT, "signatures.rds"))
sigs[["POSctrl_osteo_markers"]] <- list(source_type="positive_control",
  up=c("RUNX2","SP7","COL1A1","COL1A2","BGLAP","ALPL","IBSP","SPP1","SPARC","DLX5"), down=character(0))

dat <- GetAssayData(obj, assay="RNA", layer="data")
hs <- to_human_symbols(rownames(dat), from="mouse", strict=TRUE)      # mouse -> human, strict 1:1
keep <- !is.na(hs) & !duplicated(hs); dat <- dat[keep,,drop=FALSE]; rownames(dat) <- hs[keep]
M <- build_M(dat, obj$compartment); univ <- rownames(M); osteo <- "Osteo_MSC"

res <- list()
for (nm in names(sigs)) { s <- sigs[[nm]]; sc <- score_comps(M, s$up, s$down); if (is.null(sc)) next
  L <- (sc[osteo]-mean(sc))/sd(sc); rk <- rank(-sc)[osteo]
  nu <- length(intersect(s$up,univ)); nd <- length(intersect(s$down,univ)); Ln <- numeric(B)
  for (b in seq_len(B)){ x<-colMeans(M[sample(univ,nu),,drop=FALSE]); if(nd>=3) x<-x-colMeans(M[sample(univ,nd),,drop=FALSE]); Ln[b]<-(x[osteo]-mean(x))/sd(x) }
  res[[nm]] <- data.table(atlas="GSE316924_mouse_mandible_homeostatic", signature=nm,
    source_type=s$source_type, osteo_rank=unname(rk), transferability_z=round((L-mean(Ln))/sd(Ln),3),
    p_localizes_more=signif((1+sum(Ln>=L))/(B+1),3), localises = unname(rk)==1 & (1+sum(Ln>=L))/(B+1)<0.05) }
res <- rbindlist(res)
save_tsv(as.data.frame(res), file.path(OUT, "invivo_localization_GSE316924_homeostatic.tsv"))
saveRDS(obj, "06_scrna/GSE316924_seurat.rds")
cat("\n== homeostatic orofacial-bone reference (GSE316924) — POSITIVE CONTROL MUST localise ==\n")
print(res[order(-transferability_z)])
cat("\nGATE: if POSctrl_osteo_markers does not have osteo_rank==1 & p<0.05, this atlas is NOT",
    "method-validated and should be dropped (do not report signature results from it).\n")
