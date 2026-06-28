# R/real_run_step10_integration.R — REAL multi-evidence priority table over the
# 52-gene core program, combining the actually-computed outputs:
#   rra_strength (core RRA) | failure_specific (scRNA pseudobulk) | scrna_support
#   (per-celltype DE) | mr_causal (eQTLGen->eBMD MR) | druggability (DGIdb) |
#   jaw_context (GSE58474) | cmap (Enrichr/L1000 drug for the gene). Transparent weights.
.libPaths(c("F:/Rlib", .libPaths()))
readif <- function(p) if(file.exists(p)) read.delim(p, stringsAsFactors=FALSE) else NULL
mm <- function(x){ x<-suppressWarnings(as.numeric(x)); if(all(is.na(x))) return(rep(0,length(x)))
  r<-range(x,na.rm=TRUE); if(diff(r)==0) return(rep(0,length(x))); y<-(x-r[1])/diff(r); y[is.na(y)]<-0; y }

core <- read.delim("04_modules_wgcna/core_program.tsv")
tgt <- data.frame(gene=core$gene, direction=core$direction, stringsAsFactors=FALSE)
tgt$rra_strength <- mm(-core$rra_score[match(tgt$gene, core$gene)])

fs <- readif("05_projection/failure_specific_program.tsv")
tgt$failure_specific <- if(!is.null(fs)) mm(ifelse(fs$failure_specific[match(tgt$gene,fs$gene)], fs$abs_effect[match(tgt$gene,fs$gene)], 0)) else 0

sc <- readif("06_scrna/scrna_dysregulated_genes.tsv")
if(!is.null(sc)){ agg<-tapply(sc$abs_avg_log2FC[sc$is_dysregulated %in% c(TRUE,"TRUE")], sc$gene[sc$is_dysregulated %in% c(TRUE,"TRUE")], max)
  tgt$scrna_support <- mm(as.numeric(agg[tgt$gene])) } else tgt$scrna_support <- 0

# CORRECTED MR: lead-cis-SNP Wald (SMR-style), immune to the LD-clumping bug that polluted
# the old mr_results_ldclumped.tsv (now archived). A gene counts as genetic-causal only if
# (i) lead-SNP FDR<0.05 AND (ii) it is NOT excluded by the full-panel IVW + MR-Egger
# sensitivity (i.e. not IVW-null and not Egger-pleiotropic). GREM2 is excluded by this gate
# (IVW p=0.29; Egger intercept p=5.8e-11).
mr  <- readif("07_mr_coloc/mr_results_leadSNP.tsv")
sens<- readif("07_mr_coloc/mr_sensitivity_fullpanel_IVW.tsv")
pleio_excluded <- if(!is.null(sens)) sens$gene[ (!is.na(sens$egger_int_p) & sens$egger_int_p<0.05) |
                                                (!is.na(sens$pval) & sens$pval>=0.05) ] else character(0)
tgt$mr_p <- NA; tgt$mr_beta <- NA
if(!is.null(mr)){ mi<-match(tgt$gene, mr$gene); tgt$mr_p<-mr$padj[mi]; tgt$mr_beta<-mr$b[mi]
  robust <- !is.na(mr$padj[mi]) & mr$padj[mi] < 0.05 & !(tgt$gene %in% pleio_excluded)
  causal <- ifelse(robust, -log10(mr$padj[mi]), 0); tgt$mr_causal<-mm(causal)
} else tgt$mr_causal <- 0

dg <- readif("08_cmap/dgidb_core_program_drugs.tsv")
dgcol <- if(!is.null(dg)) intersect(c("gene_name","Gene"),colnames(dg))[1] else NA
drugcol <- if(!is.null(dg)) intersect(c("drug_name","Drug"),colnames(dg))[1] else NA
tgt$n_drugs <- if(!is.null(dg)) as.integer(table(toupper(dg[[dgcol]]))[toupper(tgt$gene)]) else 0
tgt$n_drugs[is.na(tgt$n_drugs)] <- 0
tgt$druggability <- mm(tgt$n_drugs)
tgt$example_drugs <- if(!is.null(dg)) sapply(tgt$gene, function(g){ d<-unique(dg[[drugcol]][toupper(dg[[dgcol]])==toupper(g)]); paste(head(d,4),collapse=";") }) else ""

jc <- readif("05_projection/jaw_position_context.tsv")
tgt$jaw_context <- if(!is.null(jc)){ eff<-jc$abs_jaw_effect[match(tgt$gene,jc$gene)]; cons<-jc$jaw_weak_direction_consistent[match(tgt$gene,jc$gene)]
  eff[!(cons %in% c(TRUE,"TRUE"))]<-0; mm(eff) } else 0

tgt$cmap_connectivity <- ifelse(tgt$n_drugs>0, 0.5, 0) + tgt$druggability*0.5   # drug-backed proxy

W <- c(rra_strength=0.20, failure_specific=0.20, scrna_support=0.10, mr_causal=0.25,
       druggability=0.10, jaw_context=0.05, cmap_connectivity=0.10)
for(d in names(W)) tgt[[paste0("c_",d)]] <- W[[d]]*tgt[[d]]
tgt$priority_score <- rowSums(tgt[,paste0("c_",names(W))])
# Tier names deliberately avoid "causal": the MR is a systemic-eBMD sanity check whose
# sign is not action-guiding and is NOT jaw/MRONJ-specific (see Discussion). "eBMD-MR-
# associated" denotes a robust lead-SNP+IVW/Egger association with heel-eBMD only.
tgt$evidence_tier <- ifelse(tgt$mr_causal>0 & tgt$n_drugs>0, "eBMD-MR-assoc + druggable (non-robust)",
                     ifelse(tgt$mr_causal>0,"eBMD-MR-associated", ifelse(tgt$failure_specific>0|tgt$scrna_support>0,"transcriptomic_hyp","program_member")))
tgt <- tgt[order(-tgt$priority_score),]
cols <- c("gene","direction","priority_score","evidence_tier",names(W),"mr_p","mr_beta","n_drugs","example_drugs")
dir.create("09_integration", showWarnings=FALSE)
write.table(tgt[,cols], "09_integration/priority_targets.tsv", sep="\t", quote=FALSE, row.names=FALSE)
cat(sprintf("Step10: %d core targets scored. tiers: %s\n", nrow(tgt), paste(names(table(tgt$evidence_tier)),table(tgt$evidence_tier),collapse=" ")))
cat("Top 15 priority targets:\n"); print(head(tgt[,c("gene","direction","priority_score","evidence_tier","n_drugs","example_drugs")],15))
