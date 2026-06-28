# R/real_run_repurpose_alt.R — CMap ALTERNATIVES (open, no clue.io account):
#  (A) Enrichr signature-connectivity: regeneration-competent up/down signature vs
#      LINCS_L1000_Chem_Pert (up/down) + DSigDB + Drug_Perturbations_from_GEO.
#      A drug that MIMICS the osteogenic program upregulates our up-genes AND
#      downregulates our down-genes => intersect the two enrichments.
#  (B) DGIdb target-based: which approved drugs target the core-program / hub genes.
# Hypothesis-generating; complements (not replaces) the L1000CDS2 query.
.libPaths(c("F:/Rlib", .libPaths()))
suppressMessages({ library(enrichR); library(data.table) })
options(enrichR.base.address = "https://maayanlab.cloud/Enrichr/")
dir.create("08_cmap", showWarnings = FALSE)

sig  <- read.delim("03_rra_meta/regeneration_competent_signature.tsv")
up   <- sig$gene[sig$direction == "up"]; dn <- sig$gene[sig$direction == "down"]
core <- read.delim("04_modules_wgcna/core_program.tsv")
hubs <- read.delim("04_modules_wgcna/ppi_hubs.tsv")$gene

drug_token <- function(term) toupper(trimws(sub("[ _-].*$", "", term)))  # leading drug name

## ---- (A) Enrichr ----
enr <- function(genes, dbs) tryCatch(enrichr(genes, dbs), error = function(e) { cat("enrichr err:", conditionMessage(e), "\n"); NULL })
up_res <- enr(up, c("LINCS_L1000_Chem_Pert_up", "DSigDB", "Drug_Perturbations_from_GEO_up"))
dn_res <- enr(dn, c("LINCS_L1000_Chem_Pert_down", "Drug_Perturbations_from_GEO_down"))

top_terms <- function(df, n = 40) if (is.null(df) || !nrow(df)) data.frame() else head(df[order(df$Adjusted.P.value), c("Term","Overlap","P.value","Adjusted.P.value")], n)

if (!is.null(up_res)) {
  for (nm in names(up_res)) {
    t <- top_terms(up_res[[nm]]); if (nrow(t)) { t$library <- nm
      write.table(t, file.path("08_cmap", paste0("enrichr_", nm, ".tsv")), sep="\t", quote=FALSE, row.names=FALSE) }
  }
}
if (!is.null(dn_res)) for (nm in names(dn_res)) { t <- top_terms(dn_res[[nm]]); if (nrow(t)) { t$library <- nm
  write.table(t, file.path("08_cmap", paste0("enrichr_", nm, ".tsv")), sep="\t", quote=FALSE, row.names=FALSE) }}

# mimic consensus: drug enriched in BOTH up-vs-L1000up AND down-vs-L1000down
mimic <- NULL
if (!is.null(up_res) && !is.null(dn_res) &&
    !is.null(up_res$LINCS_L1000_Chem_Pert_up) && !is.null(dn_res$LINCS_L1000_Chem_Pert_down)) {
  u <- up_res$LINCS_L1000_Chem_Pert_up; d <- dn_res$LINCS_L1000_Chem_Pert_down
  u <- u[u$P.value < 0.05, ]; d <- d[d$P.value < 0.05, ]
  ud <- drug_token(u$Term); dd <- drug_token(d$Term)
  both <- intersect(ud, dd)
  mimic <- data.frame(drug = both, n_up_sig = as.integer(table(ud)[both]), n_down_sig = as.integer(table(dd)[both]))
  mimic <- mimic[order(-(mimic$n_up_sig + mimic$n_down_sig)), ]
  write.table(mimic, "08_cmap/enrichr_LINCS_mimic_consensus.tsv", sep="\t", quote=FALSE, row.names=FALSE)
}
cat("\n=== (A) Enrichr DSigDB top drugs for UP signature ===\n")
if (!is.null(up_res) && !is.null(up_res$DSigDB)) print(head(up_res$DSigDB[order(up_res$DSigDB$Adjusted.P.value), c("Term","Overlap","Adjusted.P.value")], 15))
cat("\n=== (A) Enrichr LINCS mimic consensus (up&down) ===\n"); if (!is.null(mimic)) print(head(mimic, 20))

## ---- (B) DGIdb target-based ----
dgi_file <- "00_data/drug/dgidb_interactions.tsv"; dir.create("00_data/drug", showWarnings = FALSE, recursive = TRUE)
ok <- TRUE
if (!file.exists(dgi_file)) ok <- tryCatch({ download.file("https://www.dgidb.org/data/latest/interactions.tsv", dgi_file, quiet = TRUE, mode = "wb"); TRUE },
                                          error = function(e) tryCatch({ download.file("https://dgidb.org/data/latest/interactions.tsv", dgi_file, quiet=TRUE, mode="wb"); TRUE }, error = function(e2) { cat("DGIdb download failed:", conditionMessage(e2), "\n"); FALSE }))
if (ok && file.exists(dgi_file)) {
  dg <- tryCatch(fread(dgi_file), error = function(e) NULL)
  if (!is.null(dg)) {
    gcol <- intersect(c("gene_name","gene_claim_name","Gene"), colnames(dg))[1]
    dcol <- intersect(c("drug_name","drug_claim_name","Drug"), colnames(dg))[1]
    dg$.g <- toupper(dg[[gcol]])
    hit_core <- dg[dg$.g %in% toupper(core$gene), ]
    hit_hub  <- dg[dg$.g %in% toupper(hubs), ]
    sel <- unique(c("interaction_types","interaction_source_db_name","interaction_score"))
    sel <- intersect(sel, colnames(dg))
    outc <- unique(hit_core[, c(gcol, dcol, sel), with = FALSE])
    write.table(outc, "08_cmap/dgidb_core_program_drugs.tsv", sep="\t", quote=FALSE, row.names=FALSE)
    cat(sprintf("\n=== (B) DGIdb: %d drug-gene interactions hitting core program (%d genes); %d hitting hubs ===\n",
                nrow(hit_core), length(unique(hit_core$.g)), nrow(hit_hub)))
    print(head(unique(hit_hub[, c(gcol, dcol), with = FALSE]), 25))
  }
}
cat("\nDONE. CMap-alternative outputs in 08_cmap/ (enrichr_*.tsv, dgidb_*.tsv).\n")
