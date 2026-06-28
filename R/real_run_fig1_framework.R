# R/real_run_fig1_framework.R — redraw the conceptual framework schematic (Fig 1)
# to match the current paper: in-vitro program -> human single-cell transferability
# audit -> exploratory genetic/repurposing, with the honest findings on the boxes.
.libPaths(c("F:/Rlib", .libPaths()))
suppressMessages({ library(ggplot2) })
box <- function(x,y,w,h,fill) annotate("rect", xmin=x-w/2, xmax=x+w/2, ymin=y-h/2, ymax=y+h/2,
                                        fill=fill, color="grey30", linewidth=.4, alpha=.92)
txt <- function(x,y,label,size=3,fontface="plain",col="black",hj=0.5)
  annotate("text", x=x, y=y, label=label, size=size, fontface=fontface, color=col, hjust=hj, lineheight=.95)
arr <- function(x,y1,y2) annotate("segment", x=x, xend=x, y=y1, yend=y2,
                                  arrow=grid::arrow(length=unit(.18,"cm"), type="closed"), linewidth=.5, color="grey25")

p <- ggplot() + xlim(0,10) + ylim(0,10) +
  txt(5,9.7,"A pan–dental-source in-vitro osteogenic program and a human single-cell audit of its in-vivo transferability",3.4,"bold") +

  # Box 1 — validated in-vitro program
  box(5,8.4,9.2,1.5,"#dbe7f3") +
  txt(5,8.85,"1 · IN-VITRO OSTEOGENIC PROGRAM  (validated deliverable)",3.1,"bold","#1f3864") +
  txt(5,8.25,"RRA meta-analysis across 5 clean dental-MSC sources (PDLSC·DPSC·GSC·SHED·DFC) + leave-one-out stability + WGCNA core",2.5) +
  txt(5,7.85,"744 LOO-stable genes (350↑/394↓) · 52-gene core · separable in vitro (GSE316449, AUC=1.0 but underpowered)",2.5,"italic","#333333") +

  arr(5,7.55,7.0) +

  # Box 2 — human transferability audit (the new pillar)
  box(5,5.7,9.2,2.4,"#f3dede") +
  txt(5,6.65,"2 · HUMAN SINGLE-CELL TRANSFERABILITY AUDIT  (the load-bearing result)",3.1,"bold","#7a1f1f") +
  txt(5,6.15,"GSE303003 human jaw: 4 MRONJ vs 4 cyst (65,454 cells)   +   GSE104473 mouse distraction osteogenesis",2.5) +
  txt(2.6,5.55,"✓ data validated:\nosteoclast/myeloid\nexpansion (p=0.029)",2.3,"plain","#1a7a1a",0.5) +
  txt(5,5.55,"✗ NOT suppressed in MRONJ\nmesenchyme (sign test 43.9%,\nbelow chance)",2.3,"plain","#7a1f1f",0.5) +
  txt(7.6,5.55,"✗ NO osteoprogenitor home:\nprogram highest in myeloid,\nlowest in mesenchymal\n(osteo-markers positive ctrl peak there)",2.3,"plain","#7a1f1f",0.5) +
  txt(5,4.75,"⇒ a reproducible in-vitro program need not correspond to any discrete in-vivo osteoprogenitor state (mouse + human converge)",2.5,"italic","#7a1f1f") +

  arr(5,4.35,3.8) +

  # Box 3 — exploratory genetic + repurposing
  box(5,2.9,9.2,1.7,"#e6e6e6") +
  txt(5,3.5,"3 · EXPLORATORY GENETIC + REPURPOSING  (no robust nomination)",3.1,"bold","#444444") +
  txt(5,3.0,"Druggable-genome MR (lead-cis-SNP Wald + full-panel IVW/MR-Egger) vs heel-eBMD; FinnGen osteoporosis replication",2.5) +
  txt(5,2.55,"GREM2 = directional-pleiotropy artifact (Egger p=6e-11) · only weak SRGN survives ⇒ no robust druggable target",2.4,"italic","#333333") +
  txt(5,2.15,"Connectivity map (4 methods): dexamethasone positive control + HDAC inhibitors — methodological sanity check",2.4,"italic","#333333") +

  # honest footer banner
  box(5,0.85,9.2,0.9,"#fff3cd") +
  txt(5,0.85,"Honest scope: reproducible in-vitro resource + human-anchored lesson that in-vitro transcriptomic reproducibility does not imply in-vivo cellular correspondence.\nISG culture signal (incl. core hub EPSTI1) persists through RRA+LOO+WGCNA.  No therapeutic or mechanistic claim.",2.3,"plain","#665c00") +

  theme_void()
ggsave("figures/pub/Fig1_framework.png", p, width=10, height=8, dpi=150)
ggsave("figures/pub/Fig1_framework.pdf", p, width=10, height=8)
cat("wrote Fig1_framework (redrawn for human-audit framing)\n")
