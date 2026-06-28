# R/real_run_pubfigures.R — publication-grade multi-panel figures (patchwork) from
# the real result tables. Panel letters via plot_annotation(tag_levels="A").
.libPaths(c("F:/Rlib", .libPaths()))
suppressMessages({ library(ggplot2); library(data.table); library(patchwork) })
dir.create("figures/pub", showWarnings=FALSE)
th <- theme_bw(base_size=10) + theme(plot.title=element_text(size=10,face="bold"),
        plot.tag=element_text(size=13,face="bold"), legend.key.size=unit(.4,"cm"))
sv <- function(p, f, w, h) { ggsave(file.path("figures/pub",f), p, width=w, height=h, dpi=300); ggsave(sub("\\.pdf$",".png",file.path("figures/pub",f)), p, width=w, height=h, dpi=150) }

## ---------- Fig2: regeneration-competent program ----------
sig <- fread("03_rra_meta/regeneration_competent_signature.tsv")
ranked <- list.files("02_per_dataset_DE","_ranked.tsv$",full.names=TRUE)
topg <- head(sig$gene[order(sig$rra_score)],30)
hm <- sapply(ranked, function(f){ d<-fread(f); setNames(d$stat,d$gene)[topg] }); rownames(hm)<-topg
colnames(hm)<-sub("_ranked.tsv","",basename(ranked)); hm[is.na(hm)]<-0
hmm <- melt(as.data.table(hm,keep.rownames="gene"),id.vars="gene")
hmm$gene <- factor(hmm$gene, levels=rev(topg))
p2a <- ggplot(hmm, aes(variable,gene,fill=pmax(pmin(value,8),-8))) + geom_tile() +
  scale_fill_gradient2(low="#2166ac",mid="white",high="#b2182b",name="DE stat") +
  labs(title="Top-30 in-vitro osteogenic program genes across 5 dental sources",x=NULL,y=NULL) + th +
  theme(axis.text.x=element_text(angle=35,hjust=1,size=7), axis.text.y=element_text(size=6))
loo <- fread("03_rra_meta/loo_stability.tsv")
p2b <- ggplot(loo, aes(loo_fraction, fill=direction)) + geom_histogram(bins=12,position="identity",alpha=.6) +
  scale_fill_manual(values=c(up="#b2182b",down="#2166ac")) +
  labs(title="Leave-one-dataset-out stability (744 stable genes)",x="LOO retention fraction",y="genes") + th
gof <- "04_modules_wgcna/enrichment_signature_GO_Biological_Process_2023.tsv"
p2c <- if(file.exists(gof)){ go<-fread(gof)[1:8]; go$Term<-factor(sub(" \\(GO.*","",go$Term),levels=rev(sub(" \\(GO.*","",go$Term)))
  ggplot(go,aes(-log10(Adjusted.P.value),Term))+geom_col(fill="#4d9221")+labs(title="GO-BP enrichment (signature)",x="-log10 FDR",y=NULL)+th+theme(axis.text.y=element_text(size=7))
} else plot_spacer()
sv((p2a | (p2b/p2c)) + plot_annotation(tag_levels="A"), "Fig2_program.pdf", 12, 6)

## ---------- Fig3: validation (in-sample / external / success-null) ----------
ins <- fread("05_projection/osteogenic_score_insample.tsv")
p3a <- ggplot(ins,aes(group,osteo_score,fill=group))+geom_boxplot(alpha=.6,outlier.shape=NA)+geom_jitter(width=.15,size=1)+
  facet_wrap(~dataset,nrow=1,scales="free_y")+scale_fill_manual(values=c(control="#4575b4",osteo="#d73027"),guide="none")+
  labs(title="In-sample consistency (pooled AUC=1.0; not independent)",x=NULL,y="osteogenic score")+th+theme(strip.text=element_text(size=6))
ext <- fread("05_projection/external_validation_GSE316449.tsv")
p3b <- ggplot(ext,aes(group,osteo_score,fill=group))+geom_boxplot(alpha=.6,outlier.shape=NA)+geom_jitter(width=.12,size=2)+
  scale_fill_manual(values=c(control="#4575b4",osteo="#d73027"),guide="none")+
  labs(title="External validation GSE316449 (AUC=1.0, p=0.071)",x=NULL,y="osteogenic score")+th
suc <- fread("05_projection/success_GSE104473_scores.tsv")
p3c <- ggplot(suc,aes(celltype,score,fill=celltype))+geom_boxplot(alpha=.6,outlier.shape=NA)+geom_jitter(aes(shape=condition),width=.15,size=1.6)+
  scale_fill_brewer(palette="Set2",guide="none")+
  labs(title="Does NOT transfer in vivo (GSE104473; rho=-0.14, p=0.56)",x="sorted skeletal cell",y="in-vitro osteogenic score")+th
sv((p3a/(p3b|p3c)) + plot_annotation(tag_levels="A"), "Fig3_validation.pdf", 11, 7)

## ---------- Fig4: failure niche ----------
pt <- melt(fread("06_scrna/sc_celltype_props.tsv"),id.vars="cell_type",variable.name="condition",value.name="prop")
p4a <- ggplot(pt,aes(reorder(cell_type,-prop),prop,fill=condition))+geom_col(position="dodge")+
  scale_fill_manual(values=c(Control="#4575b4",BRONJ="#d73027"))+
  labs(title="Marrow compartment shift in BRONJ (22,991 cells)",x=NULL,y="proportion")+th+theme(axis.text.x=element_text(angle=35,hjust=1,size=8))
cc <- fread("06_scrna/cellchat_pathway_diff.tsv"); cc <- cc[is_focus==TRUE | abs(delta_BRONJ_minus_Control)>0.2][order(delta_BRONJ_minus_Control)]
cc$pathway <- factor(cc$pathway,levels=cc$pathway)
p4b <- ggplot(cc,aes(pathway,delta_BRONJ_minus_Control,fill=delta_BRONJ_minus_Control>0))+geom_col()+coord_flip()+
  scale_fill_manual(values=c(`TRUE`="#d73027",`FALSE`="#4575b4"),guide="none")+
  labs(title="CellChat rewiring: SPP1/TNF up, TGFb/BMP down",x=NULL,y="signaling delta (BRONJ-Control)")+th
sv((p4a|p4b) + plot_annotation(tag_levels="A"), "Fig4_failure_niche.pdf", 11, 4.5)

## ---------- Fig5: MR forest ----------
# NOTE: see R/real_run_fig5fig6_corrected.R for the authoritative Fig5/Fig6 used in the
# manuscript (lead-cis-SNP Wald MR, positive-control framing). Block below is legacy.
mr <- fread("07_mr_coloc/mr_results_leadSNP.tsv"); mr <- mr[!is.na(b)&!is.na(padj)][order(padj)][1:22]
mr$lab <- factor(mr$gene,levels=rev(mr$gene))
p5 <- ggplot(mr,aes(b,lab,color=direction_signature))+geom_vline(xintercept=0,linetype=2,color="grey60")+
  geom_point(aes(size=-log10(padj)))+geom_errorbarh(aes(xmin=b-1.96*se,xmax=b+1.96*se),height=0,na.rm=TRUE)+
  scale_color_manual(values=c(up="#d73027",down="#4575b4"),name="signature dir")+
  labs(title="Exploratory MR sanity check on eBMD (lead cis-SNP Wald): recovers WNT16, KREMEN1, MGP, SMAD3",
       x="MR effect of higher expression on eBMD (beta)",y=NULL,size="-log10 FDR")+th
sv(p5 + plot_annotation(tag_levels="A"), "Fig5_MR.pdf", 8, 6)

## ---------- Fig6: integrated priority ----------
pr <- fread("09_integration/priority_targets.tsv")[order(-priority_score)][1:20]; pr$gene<-factor(pr$gene,levels=rev(pr$gene))
p6 <- ggplot(pr,aes(priority_score,gene,fill=evidence_tier))+geom_col()+
  scale_fill_brewer(palette="Set1")+
  labs(title="Integrated priority targets (top 20/52; no robust genetic double-evidence — see text)",
       x="weighted priority score",y=NULL,fill="evidence tier")+th
sv(p6 + plot_annotation(tag_levels="A"), "Fig6_priority.pdf", 9, 6)

cat("publication figures written to figures/pub/\n"); print(list.files("figures/pub"))
