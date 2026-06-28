# R/real_run_cmap.R — REAL CMap coda via open L1000CDS2 API with the actual
# regeneration-competent signature. mimic mode = perturbagens that INDUCE the
# osteogenic program (pro-osteogenic repurposing candidates). Hypothesis-generating
# coda only (cf. Brey 2015 PNAS PMID 26420877); no efficacy claimed.
.libPaths(c("F:/Rlib", .libPaths()))
suppressMessages({ library(httr); library(jsonlite) })
sig <- read.delim("03_rra_meta/regeneration_competent_signature.tsv")
sig <- sig[order(sig$rra_score), ]
up <- head(sig$gene[sig$direction=="up"], 100)
dn <- head(sig$gene[sig$direction=="down"], 100)
cat("query:", length(up), "up /", length(dn), "down\n")

query_cds2 <- function(up, dn, aggravate) {
  body <- list(data=list(upGenes=as.list(up), dnGenes=as.list(dn)),
               config=list(aggravate=aggravate, searchMethod="geneSet", share=FALSE,
                           combination=FALSE, `db-version`="latest"), metadata=list())
  r <- tryCatch(POST("https://maayanlab.cloud/L1000CDS2/query",
                     body=toJSON(body, auto_unbox=TRUE), encode="raw",
                     content_type("application/json"), timeout(120)),
                error=function(e){ cat("POST err:", conditionMessage(e),"\n"); NULL })
  if (is.null(r) || http_error(r)) { cat("no response (aggravate=",aggravate,")\n"); return(NULL) }
  js <- fromJSON(content(r,"text",encoding="UTF-8"), simplifyVector=FALSE)
  tm <- js$topMeta; if (is.null(tm) || !length(tm)) return(NULL)
  do.call(rbind, lapply(tm, function(h) data.frame(
    drug=(h$pert_desc %||% h$pert_id) %||% NA, pert_id=h$pert_id %||% NA,
    cell=h$cell_id %||% NA, dose=h$pert_dose %||% NA, score=as.numeric(h$score %||% NA),
    mode=if(aggravate)"mimic" else "reverse", stringsAsFactors=FALSE)))
}
`%||%` <- function(a,b) if(is.null(a)||length(a)==0) b else a

mim <- query_cds2(up, dn, TRUE)
dir.create("08_cmap", showWarnings=FALSE)
if (!is.null(mim)) {
  mim <- mim[!is.na(mim$drug) & mim$drug!="-", ]
  write.table(mim, "08_cmap/cmap_L1000CDS2_mimic.tsv", sep="\t", quote=FALSE, row.names=FALSE)
  # aggregate by drug: how many cell-line signatures, mean score
  agg <- aggregate(score~drug, mim, function(x) c(n=length(x), mean=mean(x)))
  agg <- data.frame(drug=agg$drug, n_sig=agg$score[,"n"], mean_score=round(agg$score[,"mean"],3))
  agg <- agg[order(-agg$n_sig, agg$mean_score), ]
  write.table(agg, "08_cmap/cmap_mimic_drug_ranked.tsv", sep="\t", quote=FALSE, row.names=FALSE)
  cat("\n=== top pro-osteogenic (mimic) repurposing candidates (real L1000CDS2) ===\n")
  print(head(agg, 20))
  cat("\ntotal perturbagen signatures returned:", nrow(mim), "| unique drugs:", nrow(agg), "\n")
} else cat("L1000CDS2 returned nothing (API/offline).\n")
