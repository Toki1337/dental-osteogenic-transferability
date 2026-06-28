# R/real_run_figures.R — real figures from computed results (Fig4 failure niche,
# Fig6 MR causal, Fig8 integrated priority).
.libPaths(c("F:/Rlib", .libPaths())); suppressMessages({ library(ggplot2); library(data.table) })
dir.create("figures", showWarnings=FALSE)

## Fig4a compartment shift
pt <- fread("06_scrna/sc_celltype_props.tsv")
m <- melt(pt, id.vars="cell_type", variable.name="condition", value.name="prop")
p4 <- ggplot(m, aes(reorder(cell_type,-prop), prop, fill=condition)) +
  geom_col(position="dodge") + scale_fill_manual(values=c(Control="#4575b4",BRONJ="#d73027")) +
  labs(title="Mandibular marrow compartment shift in BRONJ (failure niche)",
       subtitle="GSE295106 scRNA, 22,991 cells; osteoclasts & erythroid up, lymphoid down",
       x=NULL, y="proportion") + theme_bw(base_size=11) + theme(axis.text.x=element_text(angle=35,hjust=1))
ggsave("figures/Fig4a_compartment_shift.png", p4, width=7, height=4, dpi=130)

## Fig4b CellChat focus pathways
cc <- fread("06_scrna/cellchat_pathway_diff.tsv")
foc <- cc[is_focus==TRUE | abs(delta_BRONJ_minus_Control) > 0.2]
foc <- foc[order(delta_BRONJ_minus_Control)]
foc$pathway <- factor(foc$pathway, levels=foc$pathway)
p4b <- ggplot(foc, aes(pathway, delta_BRONJ_minus_Control, fill=delta_BRONJ_minus_Control>0)) +
  geom_col() + coord_flip() + scale_fill_manual(values=c(`TRUE`="#d73027",`FALSE`="#4575b4"), guide="none") +
  labs(title="Cell-cell communication rewiring in BRONJ (CellChat)",
       subtitle="up: SPP1/TNF (inflammatory); down: TGFb/BMP/PTN/PDGF (regenerative)",
       x=NULL, y="signaling strength delta (BRONJ - Control)") + theme_bw(base_size=11)
ggsave("figures/Fig4b_cellchat_rewiring.png", p4b, width=7, height=4.5, dpi=130)

## Fig6 MR top causal genes
mr <- fread("07_mr_coloc/mr_results.tsv")
mr <- mr[!is.na(b) & !is.na(padj)][order(padj)][1:25]
mr$lab <- factor(mr$gene, levels=rev(mr$gene))
p6 <- ggplot(mr, aes(b, lab, color=direction_signature)) +
  geom_vline(xintercept=0, linetype=2, color="grey50") +
  geom_point(aes(size=-log10(padj))) +
  geom_errorbarh(aes(xmin=b-1.96*se, xmax=b+1.96*se), height=0, na.rm=TRUE) +
  scale_color_manual(values=c(up="#d73027",down="#4575b4"), name="signature dir") +
  labs(title="MR: regeneration-competent genes causally affect bone density (eBMD)",
       subtitle="eQTLGen cis-instruments -> eBMD (Morris 2019); top by FDR. WNT16/MGP/SMAD3 etc.",
       x="MR causal effect on eBMD (beta)", y=NULL, size="-log10 FDR") + theme_bw(base_size=10)
ggsave("figures/Fig6_MR_causal_eBMD.png", p6, width=8, height=6, dpi=130)

## Fig8 integrated priority
pr <- fread("09_integration/priority_targets.tsv")[order(-priority_score)][1:20]
pr$gene <- factor(pr$gene, levels=rev(pr$gene))
p8 <- ggplot(pr, aes(priority_score, gene, fill=evidence_tier)) + geom_col() +
  labs(title="Integrated multi-evidence priority targets (top 20 of 52 core program)",
       subtitle="7 evidence dims (RRA/failure/scRNA/MR/druggable/jaw/CMap); double-evidence: GREM2/TLR4/LOXL2/EDNRB...",
       x="weighted priority score", y=NULL, fill="evidence tier") + theme_bw(base_size=10)
ggsave("figures/Fig8_priority_targets.png", p8, width=9, height=6, dpi=130)
cat("figures written: Fig4a, Fig4b, Fig6, Fig8\n")
