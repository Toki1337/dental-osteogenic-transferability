# R/tier1_build_signatures.R
# ============================================================================
# Tier 1-A: assemble a PANEL of independent in-vitro osteogenic / MSC-lineage
# signatures, so the in-vivo transferability audit can test whether
# "no in-vivo osteoprogenitor home" is SYSTEMATIC across signatures rather than
# a peculiarity of our own RRA meta-signature.
#
# Panel composition (all human symbols):
#   (i)   our RRA-744 pan-dental meta-signature                 [own meta]
#   (ii)  five single-source clean osteogenic-induction DE      [single-source]
#         signatures (PDLSC/DPSC/GSC/SHED/DFC)
#   (iii) curated prior-knowledge osteogenic gene sets from GO  [literature]
#         (osteoblast differentiation / ossification / biomineralisation /
#          regulation of bone mineralisation)
#   (iv)  re-mined signatures from EXCLUDED perturbation-osteogenic series
#         are appended by R/tier1_mine_excluded.R (kept separate so the
#         clean+literature core stands on its own).
#
# Outputs (10_transferability/):
#   signatures.rds          - named list; each = list(name, source_type,
#                             organism, up, down, citation, notes)
#   signature_registry.tsv  - one row per signature (tracked)
#   signatures_long.tsv     - (signature, direction, gene) portable form (tracked)
# ============================================================================
.libPaths(c("F:/Rlib", .libPaths()))
suppressWarnings(suppressMessages({
  library(data.table); library(AnnotationDbi); library(org.Hs.eg.db)
}))
setDTthreads(1)
source("R/utils.R")
OUT <- "10_transferability"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

CAP <- 250L            # cap genes per direction so sizes are comparable
MINDIR <- 20L          # if a direction has < MINDIR significant genes, fall back to top-|stat|

sigs <- list()

## ---- (i) our RRA-744 pan-dental meta-signature ----------------------------
rra <- fread("03_rra_meta/regeneration_competent_signature_dedup.tsv")
sigs[["RRA744_pan_dental"]] <- list(
  name = "RRA744_pan_dental", source_type = "own_meta", organism = "human",
  up   = rra[direction == "up", gene],
  down = rra[direction == "down", gene],
  citation = "this study (RRA across 5 dental sources, LOO-stable)",
  notes = "primary meta-signature under audit")

## ---- (ii) five single-source clean osteogenic-induction DE signatures ------
mk_de_sig <- function(f, name, src) {
  d <- fread(f)
  up <- d[padj < 0.05 & log2FC >  1][order(-stat), gene]
  dn <- d[padj < 0.05 & log2FC < -1][order( stat), gene]
  if (length(up) < MINDIR) up <- d[order(-stat), gene][seq_len(min(MINDIR, nrow(d)))]
  if (length(dn) < MINDIR) dn <- d[order( stat), gene][seq_len(min(MINDIR, nrow(d)))]
  up <- head(unique(up), CAP); dn <- head(unique(dn), CAP)
  list(name = name, source_type = "single_source_DE", organism = "human",
       up = up, down = dn, citation = src,
       notes = sprintf("top<=%d/dir by stat (sig padj<0.05,|lfc|>1 else top-|stat|)", CAP))
}
ss <- list(
  c("02_per_dataset_DE/GSE99958_ranked.tsv",      "DE_PDLSC_GSE99958",   "GSE99958 PDLSC osteo-vs-undiff"),
  c("02_per_dataset_DE/GSE226347_DPSC_ranked.tsv","DE_DPSC_GSE226347",   "GSE226347 DPSC growth-vs-osteo"),
  c("02_per_dataset_DE/GSE226347_GSC_ranked.tsv", "DE_GSC_GSE226347",    "GSE226347 GSC growth-vs-osteo (near-null)"),
  c("02_per_dataset_DE/GSE266257_ranked.tsv",     "DE_SHED_GSE266257",   "GSE266257 SHED ctrl-vs-Pi"),
  c("02_per_dataset_DE/GSE49007_ranked.tsv",      "DE_DFC_GSE49007",     "GSE49007 DFC ctrl-vs-diff"))
for (x in ss) sigs[[x[2]]] <- mk_de_sig(x[1], x[2], x[3])

## ---- (iii) curated prior-knowledge osteogenic gene sets (GO, human) --------
# GOALL maps each GO id PLUS all descendant terms to genes -> a complete program.
go_genes <- function(goid) {
  g <- tryCatch(AnnotationDbi::select(org.Hs.eg.db, keys = goid,
                  columns = "SYMBOL", keytype = "GOALL")$SYMBOL,
                error = function(e) character(0))
  sort(unique(g[!is.na(g) & nzchar(g)]))
}
go_sets <- list(
  c("GO_osteoblast_differentiation", "GO:0001649", "GO:0001649 osteoblast differentiation"),
  c("GO_ossification",               "GO:0001503", "GO:0001503 ossification"),
  c("GO_biomineral_tissue_dev",      "GO:0031214", "GO:0031214 biomineral tissue development"),
  c("GO_reg_bone_mineralization",    "GO:0030500", "GO:0030500 regulation of bone mineralization"))
for (x in go_sets) {
  gg <- go_genes(x[2])
  sigs[[x[1]]] <- list(name = x[1], source_type = "literature_GO", organism = "human",
    up = head(gg, 1000L), down = character(0),
    citation = x[3], notes = "prior-knowledge gene set (unsigned; treated as osteogenic up-program)")
}

## ---- write outputs ---------------------------------------------------------
save_rds(sigs, file.path(OUT, "signatures.rds"))

reg <- rbindlist(lapply(sigs, function(s) data.table(
  signature = s$name, source_type = s$source_type, organism = s$organism,
  n_up = length(s$up), n_down = length(s$down), n_total = length(s$up) + length(s$down),
  citation = s$citation, notes = s$notes)))
save_tsv(as.data.frame(reg), file.path(OUT, "signature_registry.tsv"))

long <- rbindlist(lapply(sigs, function(s) rbindlist(list(
  if (length(s$up))   data.table(signature = s$name, direction = "up",   gene = s$up)   else NULL,
  if (length(s$down)) data.table(signature = s$name, direction = "down", gene = s$down) else NULL))))
save_tsv(as.data.frame(long), file.path(OUT, "signatures_long.tsv"))

.log("Tier1-A built ", length(sigs), " signatures")
print(reg[, .(signature, source_type, n_up, n_down)])
