# R/real_run_signaturesearch.R — offline LINCS CMap (signatureSearch), 4th repurposing
# cross-check. Maps the regeneration-competent signature to Entrez, runs LINCS GESS.
.libPaths(c("F:/Rlib", .libPaths()))
suppressMessages({ library(signatureSearch); library(org.Hs.eg.db); library(ExperimentHub); library(AnnotationDbi) })
dir.create("08_cmap", showWarnings=FALSE)
sig <- read.delim("03_rra_meta/regeneration_competent_signature.tsv")
toentrez <- function(s) na.omit(unname(mapIds(org.Hs.eg.db, s, "ENTREZID", "SYMBOL")))
up <- toentrez(head(sig$gene[order(sig$rra_score)][sig$direction=="up"], 150))
dn <- toentrez(head(sig$gene[order(sig$rra_score)][sig$direction=="down"], 150))
cat("entrez up/dn:", length(up), length(dn), "\n")

ok <- tryCatch({
  eh <- ExperimentHub()
  q <- AnnotationHub::query(eh, c("signatureSearchData"))
  cat("signatureSearchData resources:\n"); print(q$title)
  # prefer a LINCS reference h5 (lincs / lincs2)
  id <- names(q)[grep("lincs", tolower(q$title))][1]
  cat("using refdb resource:", id, q$title[match(id,names(q))], "\n")
  refdb <- eh[[id]]
  qs <- qSig(query=list(upset=as.character(up), downset=as.character(dn)), gess_method="LINCS", refdb=refdb)
  res <- gess_lincs(qs, sortby="NCS", tau=TRUE, cmp_annot_tb=NULL)
  rt <- result(res)
  write.table(head(rt, 200), "08_cmap/signaturesearch_lincs.tsv", sep="\t", quote=FALSE, row.names=FALSE)
  cat("\n=== signatureSearch LINCS top (pert / NCS / cell) ===\n")
  print(head(rt[, intersect(c("pert","NCS","WTCS","trend","cell"), colnames(rt))], 20))
  TRUE
}, error=function(e){ cat("signatureSearch FAILED:", conditionMessage(e), "\n"); FALSE })
cat("DONE ok=", ok, "\n")
