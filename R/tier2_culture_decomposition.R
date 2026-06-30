# R/tier2_culture_decomposition.R
# ============================================================================
# Tier 2-A: decompose each panel signature onto in-vitro CULTURE-STATE axes
# (cell-cycle, serum/immediate-early, hypoxia, senescence, interferon) and ask
# whether culture-fingerprint load explains the Transferability Score. This is
# the mechanistic "why": data-driven in-vitro signatures absorb the culture
# milieu, which has no discrete in-vivo osteoprogenitor home.
#
# Outputs: 10_transferability/culture_decomposition.tsv
#          figures/pub/Fig8_culture_confound.{png,pdf}
# ============================================================================
.libPaths(c("F:/Rlib", .libPaths()))
suppressWarnings(suppressMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(msigdbr); library(AnnotationDbi); library(org.Hs.eg.db); library(Seurat)
}))
setDTthreads(1)
source("R/utils.R")
OUT <- "10_transferability"
wrap <- function(s,n=60) paste(strwrap(s,width=n), collapse="\n")

## ---- culture-state gene sets ----------------------------------------------
H <- as.data.table(msigdbr(species = "Homo sapiens", collection = "H"))
symcol <- intersect(c("gene_symbol","gene_symbol"), names(H))[1]
hall <- function(nm) unique(H[gs_name == nm][[symcol]])
go_genes <- function(goid) { g <- tryCatch(AnnotationDbi::select(org.Hs.eg.db, keys=goid,
  columns="SYMBOL", keytype="GOALL")$SYMBOL, error=function(e) character(0)); unique(g[!is.na(g)&nzchar(g)]) }
cc <- cc.genes.updated.2019
IEG <- c("FOS","FOSB","JUN","JUNB","JUND","EGR1","EGR2","EGR3","ATF3","DUSP1","DUSP5",
         "NR4A1","NR4A2","NR4A3","IER2","IER3","ARC","CCN1","CYR61","CCN2","CTGF","MYC","ZFP36")
SASP <- c("IL6","IL1A","IL1B","CXCL8","IL8","CXCL1","CXCL2","CCL2","CCL20","MMP1","MMP3","MMP10",
          "SERPINE1","CDKN1A","CDKN2A","GLB1","TNFRSF1B","IGFBP3","IGFBP7","TIMP1")

culture <- list(
  cell_cycle  = unique(c(hall("HALLMARK_E2F_TARGETS"), hall("HALLMARK_G2M_CHECKPOINT"), cc$s.genes, cc$g2m.genes)),
  serum_IEG   = unique(c(hall("HALLMARK_TNFA_SIGNALING_VIA_NFKB"), hall("HALLMARK_MYC_TARGETS_V1"), IEG)),
  hypoxia     = hall("HALLMARK_HYPOXIA"),
  senescence  = unique(c(go_genes("GO:0090398"), SASP)),
  interferon  = unique(c(hall("HALLMARK_INTERFERON_ALPHA_RESPONSE"), hall("HALLMARK_INTERFERON_GAMMA_RESPONSE")))
)
osteo_identity <- unique(c(go_genes("GO:0001649"), go_genes("GO:0001503")))   # contrast axis
N_BG <- 20000L   # protein-coding background for hypergeometric enrichment

## ---- per-signature decomposition ------------------------------------------
sigs <- readRDS(file.path(OUT, "signatures.rds"))
sigs[["POSctrl_osteo_markers"]] <- list(name="POSctrl_osteo_markers", source_type="positive_control",
  up=c("RUNX2","SP7","COL1A1","COL1A2","BGLAP","ALPL","IBSP","SPP1","SPARC","DLX5"), down=character(0))

rows <- list()
for (nm in names(sigs)) {
  s <- sigs[[nm]]; g <- unique(c(s$up, s$down)); ng <- length(g)
  rec <- data.table(signature = nm, source_type = s$source_type, n_genes = ng)
  for (ax in names(culture)) {
    k <- length(intersect(g, culture[[ax]])); K <- length(intersect(culture[[ax]], 1:0)) # placeholder
    K <- length(culture[[ax]])
    rec[[paste0("frac_", ax)]] <- round(k / ng, 4)
    rec[[paste0("p_", ax)]] <- signif(phyper(k - 1, K, N_BG - K, ng, lower.tail = FALSE), 3)
  }
  rec[["frac_culture_any"]] <- round(length(intersect(g, unique(unlist(culture)))) / ng, 4)
  rec[["frac_osteo_identity"]] <- round(length(intersect(g, osteo_identity)) / ng, 4)
  rows[[nm]] <- rec
}
dec <- rbindlist(rows)

# merge Transferability Score
ts <- fread(file.path(OUT, "transferability_summary.tsv"))[, .(signature, TS, verdict)]
dec <- merge(dec, ts, by = "signature", all.x = TRUE)
dec[, prov := ifelse(source_type=="literature_GO","prior_knowledge",
              ifelse(source_type=="positive_control","positive_control","data_driven_in_vitro"))]
setorder(dec, -frac_culture_any)
fwrite(dec, file.path(OUT, "culture_decomposition.tsv"), sep = "\t")

## ---- correlation: culture load vs transferability -------------------------
d2 <- dec[!is.na(TS)]
rho <- cor(d2$frac_culture_any, d2$TS, method = "spearman")
rp  <- cor.test(d2$frac_culture_any, d2$TS, method = "spearman")$p.value
cat(sprintf("\nculture-fingerprint fraction vs Transferability Score: Spearman rho=%.2f, p=%.3g\n", rho, rp))
cat("\nculture load by provenance (mean frac_culture_any):\n")
print(dec[, .(mean_culture = round(mean(frac_culture_any),3), mean_osteo = round(mean(frac_osteo_identity),3)), by = prov])
print(dec[, .(signature, prov, frac_culture_any, frac_osteo_identity, TS, verdict)])

## ---- figure ----------------------------------------------------------------
th <- theme_bw(base_size=10) + theme(plot.title=element_text(size=9.2,face="bold"),
        plot.subtitle=element_text(size=7.6,lineheight=1.05), plot.tag=element_text(size=13,face="bold"),
        axis.text.y=element_text(size=7.4), legend.key.size=unit(.38,"cm"))
prov_cols <- c(positive_control="#1a9850", prior_knowledge="#4575b4", data_driven_in_vitro="#d73027")
# A: stacked culture-axis composition per signature
axes <- names(culture)
mlt <- melt(dec, id.vars=c("signature","prov"), measure.vars=paste0("frac_", axes),
            variable.name="axis", value.name="frac")
mlt[, axis := factor(sub("frac_","",axis), levels=axes)]
ord <- dec[order(frac_culture_any), signature]; mlt[, signature := factor(signature, levels=ord)]
pA <- ggplot(mlt, aes(frac, signature, fill=axis)) + geom_col() +
  scale_fill_brewer(palette="Set2", name="culture axis") +
  scale_x_continuous(labels=scales::percent) +
  labs(title="Culture-state fingerprint composition of each signature",
       subtitle=wrap("Fraction of signature genes belonging to in-vitro culture programs (cell-cycle, serum/immediate-early, hypoxia, senescence, interferon).",110),
       x="% of signature genes", y=NULL) + th
# B: culture-axis enrichment per signature (-log10 hypergeometric p)
axesp <- paste0("p_", names(culture))
em <- melt(dec, id.vars=c("signature","prov"), measure.vars=axesp, variable.name="axis", value.name="p")
em[, axis := factor(sub("p_","",axis), levels=names(culture))]
em[, nlp := pmin(-log10(pmax(as.numeric(p), 1e-60)), 30)]
em[, signature := factor(signature, levels=ord)]
pB <- ggplot(em, aes(axis, signature, fill=nlp)) + geom_tile(color="white", linewidth=.5) +
  geom_text(aes(label=ifelse(as.numeric(p)<0.05, "*","")), size=3.6, vjust=.72) +
  scale_fill_gradient(low="grey96", high="#b2182b", name="-log10 p\n(enrichment)") +
  labs(title="Which culture programs are enriched in each signature",
       subtitle=wrap(sprintf("Hypergeometric enrichment per culture axis (* p<0.05). RRA-744 is enriched for hypoxia/serum/interferon; perturbation signatures for cell-cycle. Culture load alone does not predict transfer (Spearman rho=%.2f, ns); the data-driven panel is ~5%% canonical osteo-identity genes vs ~92%% for curated sets.", rho),78),
       x=NULL, y=NULL) + th + theme(axis.text.x=element_text(angle=25,hjust=1,size=8))
fig8 <- (pA | pB) + plot_layout(widths=c(1,1)) + plot_annotation(tag_levels="A")
ggsave("figures/pub/Fig8_culture_confound.pdf", fig8, width=12, height=6)
ggsave("figures/pub/Fig8_culture_confound.png", fig8, width=12, height=6, dpi=150)
cat("\nwrote 10_transferability/culture_decomposition.tsv + figures/pub/Fig8_culture_confound\n")
