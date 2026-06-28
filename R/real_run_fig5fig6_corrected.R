# Regenerate Fig5 (MR) + Fig6 (priority) with the CORRECTED lead-SNP MR and honest framing.
.libPaths(c("F:/Rlib", .libPaths()))
suppressMessages({ library(ggplot2); library(data.table); library(patchwork) })
th <- theme_bw(base_size=10) + theme(plot.title=element_text(size=8.6,face="bold",lineheight=1.05),
        plot.subtitle=element_text(size=7.2,lineheight=1.05), plot.tag=element_text(size=13,face="bold"),
        legend.key.size=unit(.4,"cm"), plot.margin=margin(6,8,4,4))
wrap <- function(s,n=58) paste(strwrap(s,width=n), collapse="\n")
sv <- function(p, f, w, h){ ggsave(file.path("figures/pub",f), p, width=w, height=h, dpi=300)
  ggsave(sub("\\.pdf$",".png",file.path("figures/pub",f)), p, width=w, height=h, dpi=150) }

## ---- Fig5: MR sanity check (lead-cis-SNP Wald), known bone genes as positive controls ----
mr <- fread("07_mr_coloc/mr_results_leadSNP.tsv"); mr <- mr[!is.na(b)&!is.na(padj)][order(padj)]
known <- c("WNT16","KREMEN1","MGP","SMAD3","HSPG2","KLF12")
drug_de <- c("GREM2","SRGN")
top <- mr[1:22]
top$role <- ifelse(top$gene %in% known, "known bone (pos. control)",
            ifelse(top$gene %in% drug_de, "druggable core", "other signature gene"))
top$lab <- factor(top$gene, levels=rev(top$gene))
p5 <- ggplot(top, aes(b, lab, color=role)) +
  geom_vline(xintercept=0, linetype=2, color="grey60") +
  geom_errorbarh(aes(xmin=b-1.96*se, xmax=b+1.96*se), height=0, na.rm=TRUE, color="grey55") +
  geom_point(aes(size=-log10(padj))) +
  scale_color_manual(values=c(`known bone (pos. control)`="#1a9850",
        `druggable core`="#d73027", `other signature gene`="#4575b4"), name=NULL) +
  labs(title=wrap("Exploratory MR sanity check on heel-eBMD (lead cis-SNP Wald; systemic bone, not jaw/MRONJ)",62),
       subtitle=wrap("Recovers established bone genes (WNT16, KREMEN1, MGP, SMAD3, HSPG2) as positive controls; no protective direction implied",78),
       x="MR effect of higher expression on eBMD (beta)", y=NULL, size="-log10 FDR") + th +
  theme(legend.position="bottom")
sv(p5 + plot_annotation(tag_levels="A"), "Fig5_MR.pdf", 8.5, 6.2)
cat("Fig5 regenerated (lead-SNP MR; positive-control framing)\n")

## ---- Fig6: integrated priority (honest double-evidence note) ----
pr <- fread("09_integration/priority_targets.tsv")[order(-priority_score)][1:20]
pr$gene <- factor(pr$gene, levels=rev(pr$gene))
p6 <- ggplot(pr, aes(priority_score, gene, fill=evidence_tier)) + geom_col() +
  scale_fill_brewer(palette="Set1") +
  labs(title=wrap("Integrated priority across 52-gene core (transcriptomic-led)",60),
       subtitle=wrap("IGFBP5 leads on transcriptomic evidence. No robust genetic double-evidence: GREM2 lead-SNP signal is a pleiotropy artifact (Egger p=6e-11); SRGN weak/discordant (non-robust).",92),
       x="weighted priority score", y=NULL, fill="evidence tier") + th
sv(p6 + plot_annotation(tag_levels="A"), "Fig6_priority.pdf", 9, 6)
cat("Fig6 regenerated (honest double-evidence note)\n")
