# R/tier3_aucell_check.R
# ============================================================================
# Reviewer (label-transfer/AUCell item): show the Transferability Score's per-compartment
# scoring is not a method artifact by cross-checking it against an INDEPENDENT, established
# signature scorer — AUCell (Aibar et al. 2017). For each atlas and each signature (up-genes),
# we compute AUCell per-cell AUC, average it per compartment, and correlate that per-compartment
# AUC profile with the transferaudit M-based per-compartment score (up-genes). High Spearman ρ
# means the two tools agree on where each signature localises. This complements the in-paper
# AddModuleScore-vs-UCell agreement (ρ = 1.0) with a third, rank-AUC-based method.
#   Rscript R/tier3_aucell_check.R  ->  10_transferability/aucell_consistency.tsv
# ============================================================================
.libPaths(c(Sys.getenv("R_TRANSFERAUDIT_LIB", unset = .libPaths()[1]), .libPaths()))
suppressWarnings(suppressMessages({ library(SeuratObject); library(Matrix); library(data.table); library(AUCell) }))
setDTthreads(1); source("R/utils.R"); source("config/params.R")
OUT <- "10_transferability"; set.seed(PARAMS$seed)

sigs <- readRDS(file.path(OUT, "signatures.rds"))
sigs[["POSctrl_osteo_markers"]] <- list(source_type="positive_control",
  up=c("RUNX2","SP7","COL1A1","COL1A2","BGLAP","ALPL","IBSP","SPP1","SPARC","DLX5"), down=character(0))

atlases <- list(
  c("GSE303003_human_MRONJ","06_scrna/GSE303003_seurat.rds","Mesenchymal_osteo","compartment","human"),
  c("GSE295106_mouse_BRONJ","06_scrna/GSE295106_seurat.rds","Osteo_MSC","cell_type","mouse"),
  c("GSE269255_mouse_ORNJ","06_scrna/GSE269255_seurat.rds","Mesenchymal_osteo","compartment","mouse"))

res <- list()
for(a in atlases){ tag<-a[1]; obj<-readRDS(a[2]); dat<-GetAssayData(obj,assay="RNA",layer="data")
  comp<-as.character(obj@meta.data[[a[4]]])
  if(a[5]=="mouse"){ hs<-to_human_symbols(rownames(dat),from="mouse",strict=TRUE)
    keep<-!is.na(hs)&!duplicated(hs); dat<-dat[keep,,drop=FALSE]; rownames(dat)<-hs[keep] }
  # transferaudit M-based per-compartment score (up-genes only, for apples-to-apples vs AUCell)
  cl<-sort(unique(comp)); mu<-Matrix::rowMeans(dat); sdg<-sqrt(pmax(Matrix::rowMeans(dat*dat)-mu^2,0))
  M<-(sapply(cl,function(c) Matrix::rowMeans(dat[,comp==c,drop=FALSE]))-mu)/sdg
  M<-M[is.finite(rowSums(M))&sdg>0,,drop=FALSE]; univ<-rownames(M)
  # AUCell rankings (once) + AUC per signature up-set
  rk<-AUCell_buildRankings(dat, plotStats=FALSE, verbose=FALSE)   # sparse dgCMatrix input (no densify)
  gsets<-lapply(sigs, function(s) intersect(s$up, univ)); gsets<-gsets[sapply(gsets,length)>=3]
  auc<-AUCell_calcAUC(gsets, rk, verbose=FALSE, aucMaxRank=ceiling(0.05*nrow(rk)))
  aucm<-getAUC(auc)                                        # signatures x cells
  for(nm in names(gsets)){ up<-gsets[[nm]]
    m_score <- colMeans(M[up,,drop=FALSE])                 # transferaudit per-compartment (up-only)
    au_cell <- aucm[nm,]; au_comp <- tapply(au_cell, comp, mean)[cl]   # AUCell per-compartment mean AUC
    rho <- suppressWarnings(cor(m_score, au_comp, method="spearman"))
    res[[paste(tag,nm)]] <- data.table(atlas=tag, signature=nm, source_type=sigs[[nm]]$source_type,
      spearman_rho_M_vs_AUCell=round(rho,3),
      osteo_rank_transferaudit=unname(rank(-m_score)[a[3]]),
      osteo_rank_AUCell=unname(rank(-au_comp)[a[3]])) }
  rm(obj,dat,M,rk,aucm); gc() }
res<-rbindlist(res)
save_tsv(as.data.frame(res), file.path(OUT,"aucell_consistency.tsv"))
cat(sprintf("\nAUCell vs transferaudit per-compartment Spearman rho: median %.2f (IQR %.2f-%.2f), min %.2f\n",
  median(res$spearman_rho_M_vs_AUCell,na.rm=TRUE),
  quantile(res$spearman_rho_M_vs_AUCell,.25,na.rm=TRUE), quantile(res$spearman_rho_M_vs_AUCell,.75,na.rm=TRUE),
  min(res$spearman_rho_M_vs_AUCell,na.rm=TRUE)))
cat(sprintf("osteo-compartment rank agreement (transferaudit vs AUCell identical): %d/%d rows\n",
  sum(res$osteo_rank_transferaudit==res$osteo_rank_AUCell), nrow(res)))
print(res)
cat("[wrote aucell_consistency.tsv]\n")
