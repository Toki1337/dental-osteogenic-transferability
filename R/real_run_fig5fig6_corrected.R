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
sv(p5 + plot_annotation(tag_levels="A"), "Fig7_MR.pdf", 8.5, 6.2)   # numbered Fig7 (cited in MR §3.6)
cat("Fig7_MR regenerated (lead-SNP MR; positive-control framing)\n")

## ---- Fig6a: reverse-connectivity convergence of four open CMap methods ----
# Computed directly from the released 08_cmap outputs (no hand-entered values):
# for each method, count how many distinct agents of each therapeutic class were
# recovered. Dexamethasone is the osteogenic-induction-medium agent = positive ctrl.
suppressWarnings(setDTthreads(1))
rd <- function(f) tryCatch(fread(file.path("08_cmap", f)), error=function(e) data.table())
m_l1000 <- rd("cmap_L1000CDS2_mimic.tsv")          # col: drug
m_dsig  <- rd("enrichr_DSigDB.tsv")                # col: Term
m_geo   <- rbindlist(list(rd("enrichr_Drug_Perturbations_from_GEO_up.tsv"),
                          rd("enrichr_Drug_Perturbations_from_GEO_down.tsv")), fill=TRUE)
m_ss    <- rd("signaturesearch_lincs.tsv")         # cols: pert + MOAss (mechanism)
cmap_cls <- list(
 "Dexamethasone\n(induction +ctrl)" = c("dexamethasone"),
 "Glucocorticoid\n(other)"          = c("triamcinolone","betamethasone","budesonide","fluticasone",
                                        "clobetasol","halometasone","fluocinolone","mometasone",
                                        "prednisolone","hydrocortisone","flunisolide"),
 "HDAC inhibitor"                   = c("vorinostat","trichostatin","entinostat","valproic","panobinostat",
                                        "belinostat","romidepsin","mocetinostat","SAHA","MS-275","\\bHDAC\\b"),
 "Steroid hormone\n(estradiol/prog.)"= c("estradiol","progesterone","estrone"))
cmap_strs <- list("L1000CDS2"=as.character(m_l1000$drug),
                  "Enrichr/\nDSigDB"=as.character(m_dsig$Term),
                  "Enrichr/\nGEO"=as.character(m_geo$Term),
                  "sigSearch-\nLINCS"=as.character(c(m_ss$pert, m_ss$MOAss)))
cnt_cls <- function(s, kws){ s <- s[!is.na(s)]; if(!length(s)) return(0L)
  sum(sapply(kws, function(k) any(grepl(k, s, ignore.case=TRUE)))) }
conv <- rbindlist(lapply(names(cmap_strs), function(meth)
  rbindlist(lapply(names(cmap_cls), function(cn)
    data.table(method=meth, class=cn, n=cnt_cls(cmap_strs[[meth]], cmap_cls[[cn]]))))))
conv$method <- factor(conv$method, levels=names(cmap_strs))
conv$class  <- factor(conv$class,  levels=rev(names(cmap_cls)))
pA <- ggplot(conv, aes(method, class, fill=n)) +
  geom_tile(color="white", linewidth=.7) +
  geom_text(aes(label=ifelse(n>0, n, "·"), color=n>0), size=3.1, show.legend=FALSE) +
  scale_color_manual(values=c(`TRUE`="white", `FALSE`="grey65")) +
  scale_fill_gradient(low="grey93", high="#2166ac", name="agents\nrecovered") +
  labs(title=wrap("Reverse-connectivity convergence of four open CMap methods",60),
       subtitle=wrap("Dexamethasone (induction-medium agent) + glucocorticoids recovered by all four methods; HDAC inhibitors by three of four — a methodological positive-control/sanity check, not a repurposing discovery",100),
       x=NULL, y=NULL) + th +
  theme(axis.text.x=element_text(size=7.4), panel.grid=element_blank(), legend.position="right")

## ---- Fig6b: integrated priority (honest double-evidence note) ----
pr <- fread("09_integration/priority_targets.tsv")[order(-priority_score)][1:20]
pr$gene <- factor(pr$gene, levels=rev(pr$gene))
p6 <- ggplot(pr, aes(priority_score, gene, fill=evidence_tier)) + geom_col() +
  scale_fill_brewer(palette="Set1") +
  labs(title=wrap("Integrated priority across 52-gene core (transcriptomic-led)",60),
       subtitle=wrap("IGFBP5 leads on transcriptomic evidence. No robust genetic double-evidence: GREM2 lead-SNP signal is a pleiotropy artifact (Egger p=6e-11); SRGN weak/discordant (non-robust).",92),
       x="weighted priority score", y=NULL, fill="evidence tier") + th
fig6 <- (pA / p6) + plot_layout(heights=c(0.85, 1.7)) + plot_annotation(tag_levels="A")
sv(fig6, "Fig8_priority.pdf", 9, 10.5)   # numbered Fig8 (cited in CMap/priority §3.7)
cat("Fig8_priority regenerated (a: CMap 4-method convergence; b: integrated priority)\n")
