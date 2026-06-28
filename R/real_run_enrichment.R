# R/real_run_enrichment.R — REAL functional enrichment of the regeneration-competent
# core program via Enrichr (open API; clusterProfiler's source deps unavailable here).
.libPaths(c("F:/Rlib", .libPaths())); suppressMessages(library(enrichR))
options(enrichR.base.address="https://maayanlab.cloud/Enrichr/")
core <- read.delim("04_modules_wgcna/core_program.tsv")
sig  <- read.delim("03_rra_meta/regeneration_competent_signature.tsv")
dbs <- c("GO_Biological_Process_2023","KEGG_2021_Human","WikiPathways_2024_Human","MSigDB_Hallmark_2020")
dir.create("04_modules_wgcna", showWarnings=FALSE)
run <- function(genes, tag){
  r <- tryCatch(enrichr(genes, dbs), error=function(e){cat("err",tag,conditionMessage(e),"\n");NULL})
  if(is.null(r)) return(invisible())
  for(nm in names(r)){ d<-r[[nm]]; if(is.null(d)||!nrow(d)) next
    d<-d[order(d$Adjusted.P.value), c("Term","Overlap","P.value","Adjusted.P.value","Genes")]
    write.table(head(d,40), sprintf("04_modules_wgcna/enrichment_%s_%s.tsv",tag,nm), sep="\t", quote=FALSE, row.names=FALSE)
    if(nm=="GO_Biological_Process_2023") { cat("\n===",tag,"GO BP top10 ===\n"); print(head(d[,c("Term","Adjusted.P.value")],10)) }
  }
}
run(core$gene, "core")
run(sig$gene, "signature")
cat("\nDONE. enrichment_*.tsv written.\n")
