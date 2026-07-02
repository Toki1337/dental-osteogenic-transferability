# R/tier3_coexpr_null.R
# ============================================================================
# Reviewer (Major 3): the Transferability Score's size-matched RANDOM-gene-set null
# ignores within-signature co-expression and can under-estimate the null variance,
# inflating |z| / deflating p. Here we recompute the Transferability Score against an
# EXPRESSION-MATCHED (co-expression-preserving) null: null gene sets are drawn to match
# the signature's per-bin mean-expression profile (20 quantile bins) in each atlas, so
# the null sets have similar expression level — and hence similar co-expression / score
# variance — as the real signature. If the prior-knowledge-vs-data-driven split (GO/
# markers localise, data-driven do not) survives this stricter null, the qualitative
# conclusion does not depend on the random-null calibration.
#
# Reuses the exact M-matrix machinery of R/tier1_project_invivo.R. The observed osteo_rank
# is null-independent (unchanged); only z and p change.
#   Rscript R/tier3_coexpr_null.R   ->  10_transferability/invivo_localization_coexprnull.tsv
# ============================================================================
.libPaths(c(Sys.getenv("R_TRANSFERAUDIT_LIB", unset = .libPaths()[1]), .libPaths()))
suppressWarnings(suppressMessages({ library(SeuratObject); library(Matrix); library(data.table) }))
setDTthreads(1); source("R/utils.R"); source("config/params.R")
OUT <- "10_transferability"; set.seed(PARAMS$seed); B <- 2000L; NBIN <- 20L

sigs <- readRDS(file.path(OUT, "signatures.rds"))
sigs[["POSctrl_osteo_markers"]] <- list(source_type="positive_control",
  up=c("RUNX2","SP7","COL1A1","COL1A2","BGLAP","ALPL","IBSP","SPP1","SPARC","DLX5"), down=character(0))

build_M_mu <- function(dat, comp){ comp<-as.character(comp); cl<-sort(unique(comp))
  mu<-Matrix::rowMeans(dat); sdg<-sqrt(pmax(Matrix::rowMeans(dat*dat)-mu^2,0))
  M<-(sapply(cl,function(c) Matrix::rowMeans(dat[,comp==c,drop=FALSE]))-mu)/sdg
  keep<-is.finite(rowSums(M)) & sdg>0; list(M=M[keep,,drop=FALSE], mu=mu[keep]) }
score_comps <- function(M,up,dn){ up<-intersect(up,rownames(M)); dn<-intersect(dn,rownames(M))
  if(length(up)<3) return(NULL); s<-colMeans(M[up,,drop=FALSE]); if(length(dn)>=3) s<-s-colMeans(M[dn,,drop=FALSE]); s }

atlases <- list(
  c("GSE303003_human_MRONJ","06_scrna/GSE303003_seurat.rds","Mesenchymal_osteo","compartment","human"),
  c("GSE295106_mouse_BRONJ","06_scrna/GSE295106_seurat.rds","Osteo_MSC","cell_type","mouse"),
  c("GSE269255_mouse_ORNJ","06_scrna/GSE269255_seurat.rds","Mesenchymal_osteo","compartment","mouse"))

res <- list()
for(a in atlases){ tag<-a[1]; obj<-readRDS(a[2]); dat<-GetAssayData(obj,assay="RNA",layer="data")
  comp<-obj@meta.data[[a[4]]]
  if(a[5]=="mouse"){ hs<-to_human_symbols(rownames(dat),from="mouse",strict=TRUE)
    keep<-!is.na(hs)&!duplicated(hs); dat<-dat[keep,,drop=FALSE]; rownames(dat)<-hs[keep] }
  MM<-build_M_mu(dat,comp); M<-MM$M; osteo<-a[3]; univ<-rownames(M)
  bins<-cut(MM$mu, breaks=quantile(MM$mu,probs=seq(0,1,length.out=NBIN+1)), include.lowest=TRUE, labels=FALSE)
  by_bin<-split(univ, bins)
  # expression-matched draw: replicate the query's per-bin composition
  draw_matched <- function(genes){ b<-bins[match(intersect(genes,univ),univ)]; tb<-table(b)
    unlist(lapply(names(tb), function(bi){ pool<-by_bin[[bi]]; n<-tb[[bi]]
      sample(pool, min(n,length(pool)), replace=length(pool)<n) })) }
  for(nm in names(sigs)){ s<-sigs[[nm]]; sc<-score_comps(M,s$up,s$down); if(is.null(sc)) next
    L<-(sc[osteo]-mean(sc))/sd(sc); rk<-rank(-sc)[osteo]
    Ln<-numeric(B)
    for(b in seq_len(B)){ ub<-draw_matched(s$up); x<-colMeans(M[ub,,drop=FALSE])
      if(length(intersect(s$down,univ))>=3){ db<-draw_matched(s$down); x<-x-colMeans(M[db,,drop=FALSE]) }
      Ln[b]<-(x[osteo]-mean(x))/sd(x) }
    res[[paste(tag,nm)]]<-data.table(atlas=tag, signature=nm, source_type=s$source_type,
      osteo_rank=unname(rk), transferability_z_coexpr=round((L-mean(Ln))/sd(Ln),3),
      p_coexpr=signif((1+sum(Ln>=L))/(B+1),3),
      localises_coexpr = unname(rk)==1 & (1+sum(Ln>=L))/(B+1) < 0.05) }
  rm(obj,dat,M); gc() }
res<-rbindlist(res)
save_tsv(as.data.frame(res), file.path(OUT,"invivo_localization_coexprnull.tsv"))

# provenance summary under the co-expression null
lz<-dcast(res, signature+source_type~atlas, value.var="localises_coexpr")
lz$n_loc<-rowSums(lz[,-(1:2)],na.rm=TRUE)
pk<-lz[source_type=="literature_GO"]; dd<-lz[source_type %in% c("own_meta","single_source_DE","excluded_perturbation_DE")]
cat("\n== co-expression-matched null: localises in all 3 atlases ==\n")
cat(sprintf("prior-knowledge GO: %d/%d ; data-driven: %d/%d ; positive control n_loc=%s\n",
  sum(pk$n_loc==3), nrow(pk), sum(dd$n_loc==3), nrow(dd), lz[source_type=="positive_control"]$n_loc))
print(res[order(source_type,-transferability_z_coexpr)][,.(signature,atlas,osteo_rank,transferability_z_coexpr,p_coexpr)])
cat("[wrote invivo_localization_coexprnull.tsv]\n")
