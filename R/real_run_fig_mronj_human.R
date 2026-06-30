# R/real_run_fig_mronj_human.R — figure for the human MRONJ single-cell transferability audit.
.libPaths(c("F:/Rlib", .libPaths()))
suppressMessages({ library(ggplot2); library(data.table); library(patchwork); library(Seurat) })
th <- theme_bw(base_size=10) + theme(plot.title=element_text(size=9,face="bold"),
        plot.tag=element_text(size=13,face="bold"), legend.key.size=unit(.32,"cm"),
        legend.text=element_text(size=7), legend.title=element_text(size=8))
sv <- function(p,f,w,h){ ggsave(file.path("figures/pub",f),p,width=w,height=h,dpi=300)
  ggsave(sub("\\.pdf$",".png",file.path("figures/pub",f)),p,width=w,height=h,dpi=150) }
OUT <- "06_scrna"
sc <- readRDS(file.path(OUT,"GSE303003_seurat.rds"))

# A) compartment UMAP
um <- as.data.frame(Embeddings(sc,"umap")); um$compartment <- sc$compartment
set.seed(1); um <- um[sample(nrow(um)),]
pA <- ggplot(um, aes(umap_1,umap_2,color=compartment))+geom_point(size=.04,alpha=.5)+
  guides(color=guide_legend(override.aes=list(size=2,alpha=1)))+
  labs(title="Human jaw (GSE303003): 65,454 cells, 4 MRONJ vs 4 cyst, 7 compartments",x="UMAP1",y="UMAP2")+th+
  theme(axis.text=element_blank(), axis.ticks=element_blank())

# B) compartment proportions MRONJ vs Cyst (replicated osteoclast/myeloid expansion)
pt <- fread(file.path(OUT,"GSE303003_compartment_proptest.tsv"))
ptm <- melt(pt, id.vars=c("compartment","wilcox_p"), measure.vars=c("Cyst_mean","MRONJ_mean"),
            variable.name="condition", value.name="prop")
ptm$condition <- factor(ifelse(ptm$condition=="Cyst_mean","Cyst","MRONJ"), levels=c("Cyst","MRONJ"))
ptm$lab <- ifelse(ptm$wilcox_p<0.05 & ptm$condition=="MRONJ", sprintf("p=%.3g", ptm$wilcox_p), "")
pB <- ggplot(ptm, aes(reorder(compartment,-prop), 100*prop, fill=condition))+
  geom_col(position=position_dodge(width=.8), width=.7)+
  geom_text(aes(label=lab), position=position_dodge(width=.8), vjust=-0.3, size=2.6)+
  scale_fill_manual(values=c(Cyst="#4575b4", MRONJ="#d73027"), name=NULL)+
  labs(title="Compartment proportions: osteoclast &\nmyeloid expansion in MRONJ (4v4, p=0.029)",
       x=NULL, y="% of cells")+th+theme(axis.text.x=element_text(angle=30,hjust=1,size=7),
       legend.position=c(.85,.8), plot.title=element_text(size=8.6,face="bold",lineheight=1.05),
       plot.margin=margin(6,10,4,4))

# C) program localization vs osteo-marker positive control (the key robustness panel)
loc <- fread(file.path(OUT,"GSE303003_localization_robust.tsv"))
# scale each series to [-1,1]-ish for visual comparison via z within series
zc <- function(x) (x-mean(x))/sd(x)
locm <- rbind(
  data.frame(compartment=loc$compartment, series="L1 program (AddModuleScore)", z=zc(loc$prog_AMS)),
  data.frame(compartment=loc$compartment, series="L1 program (UCell, rank)",    z=zc(loc$prog_UC)),
  data.frame(compartment=loc$compartment, series="osteo markers (RUNX2/SP7/COL1A1) — positive ctrl", z=zc(loc$osteo_AMS)))
locm$compartment <- factor(locm$compartment, levels=loc$compartment[order(loc$prog_AMS)])
pC <- ggplot(locm, aes(z, compartment, fill=series))+
  geom_col(position=position_dodge(width=.75), width=.7)+geom_vline(xintercept=0,linetype=2,color="grey50")+
  scale_fill_manual(values=c("#762a83","#9970ab","#1a9850"), name=NULL)+
  labs(title="Positive control: osteo markers peak in mesenchymal, but the L1 program does NOT",
       subtitle="two scorings agree (Spearman rho=1); program highest in myeloid, lowest in mesenchymal/endothelial",
       x="per-compartment score (z within series)", y=NULL)+th+
  theme(plot.subtitle=element_text(size=7.5), legend.position="bottom", legend.direction="vertical")

# D) per-sample mesenchymal program MRONJ vs Cyst (NOT suppressed) + sign test (deduped)
ps <- fread(file.path(OUT,"GSE303003_mes_program_persample_dedup.tsv"))
st <- fread(file.path(OUT,"GSE303003_program_signtest_dedup.tsv"))
pw <- wilcox.test(prog_AMS~condition, ps)$p.value
sgn <- st[set=="signature_744_dedup"]
pD <- ggplot(ps, aes(condition, prog_AMS, fill=condition))+
  geom_boxplot(outlier.shape=NA, width=.5, alpha=.6)+geom_jitter(width=.12, size=2)+
  scale_fill_manual(values=c(Cyst="#4575b4", MRONJ="#d73027"), guide="none")+
  labs(title="Program NOT suppressed in MRONJ mesenchyme",
       subtitle=sprintf("whole-program sign test: %d/%d (%.0f%%) suppressed-direction, NOT > 50%% (p=%.2g); Wilcoxon p=%.2g",
                        sgn$k, sgn$n, 100*sgn$frac, sgn$binom_p_greater, pw),
       x=NULL, y="per-sample mean program score")+th+theme(plot.subtitle=element_text(size=7))

fig <- (pA|pB)/(pC|pD) + plot_annotation(tag_levels="A")
sv(fig, "Fig4_failure_human.pdf", 12, 8.5)
cat("wrote Fig4_failure_human (human MRONJ single-cell transferability audit)\n")
