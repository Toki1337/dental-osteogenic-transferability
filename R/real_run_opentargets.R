# R/real_run_opentargets.R — annotate priority targets with Open Targets clinical
# phase + tractability (closes the candidate table's clinical_pipeline column).
.libPaths(c("F:/Rlib", .libPaths())); suppressMessages({ library(httr); library(jsonlite) })
API <- "https://api.platform.opentargets.org/api/v4/graphql"
gql <- function(q, v=NULL){ r <- tryCatch(POST(API, body=toJSON(list(query=q, variables=v), auto_unbox=TRUE),
   content_type_json(), timeout(60)), error=function(e) NULL)
  if(is.null(r)||http_error(r)) return(NULL); fromJSON(content(r,"text",encoding="UTF-8"), simplifyVector=FALSE) }

ensembl_of <- function(sym){ d <- gql('query($s:String!){search(queryString:$s,entityNames:["target"]){hits{id name entity}}}', list(s=sym))
  h <- d$data$search$hits; if(is.null(h)||!length(h)) return(NA); for(x in h) if(toupper(x$name)==toupper(sym)) return(x$id); h[[1]]$id }
annot <- function(id){ d <- gql('query($i:String!){target(ensemblId:$i){approvedSymbol tractability{modality value label} knownDrugs{count rows{phase drug{name} disease{name}}}}}', list(i=id))
  t <- d$data$target; if(is.null(t)) return(NULL)
  kd <- t$knownDrugs; maxph <- if(!is.null(kd) && length(kd$rows)) max(vapply(kd$rows, function(r) r$phase %||% 0, numeric(1))) else NA
  drugs <- if(!is.null(kd) && length(kd$rows)) paste(unique(vapply(kd$rows, function(r) r$drug$name %||% "", character(1)))[1:min(5,length(kd$rows))], collapse=";") else ""
  smtract <- if(!is.null(t$tractability)) any(vapply(t$tractability, function(x) (x$modality %||% "")=="SM" && isTRUE(x$value), logical(1))) else NA
  list(symbol=t$approvedSymbol, max_clinical_phase=maxph, n_known_drugs=if(!is.null(kd)) kd$count else 0, top_drugs=drugs, smallmolecule_tractable=smtract) }
`%||%` <- function(a,b) if(is.null(a)) b else a

pr <- read.delim("09_integration/priority_targets.tsv")
top <- head(pr$gene, 25)
rows <- list()
for(g in top){ id <- tryCatch(ensembl_of(g), error=function(e) NA); Sys.sleep(0.2)
  if(is.na(id)){ rows[[g]] <- data.frame(gene=g, ensembl=NA, max_clinical_phase=NA, n_known_drugs=NA, top_drugs="", smallmolecule_tractable=NA); next }
  a <- tryCatch(annot(id), error=function(e) NULL); Sys.sleep(0.2)
  rows[[g]] <- if(is.null(a)) data.frame(gene=g, ensembl=id, max_clinical_phase=NA, n_known_drugs=NA, top_drugs="", smallmolecule_tractable=NA) else
    data.frame(gene=g, ensembl=id, max_clinical_phase=a$max_clinical_phase, n_known_drugs=a$n_known_drugs, top_drugs=a$top_drugs, smallmolecule_tractable=a$smallmolecule_tractable)
  cat(sprintf("%-10s phase=%s drugs=%s SMtract=%s\n", g, rows[[g]]$max_clinical_phase, rows[[g]]$n_known_drugs, rows[[g]]$smallmolecule_tractable)) }
ot <- do.call(rbind, rows)
write.table(ot, "09_integration/opentargets_annotation.tsv", sep="\t", quote=FALSE, row.names=FALSE)
# merge into final candidate table
fin <- merge(pr, ot, by="gene", all.x=TRUE); fin <- fin[order(-fin$priority_score),]
write.table(fin, "09_integration/priority_targets_annotated.tsv", sep="\t", quote=FALSE, row.names=FALSE)
cat("\nDONE. opentargets_annotation.tsv + priority_targets_annotated.tsv written.\n")
