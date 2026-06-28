# R/real_run_fig2fig3_corrected.R — regenerate ONLY Fig2 + Fig3 with the deduped
# 744-gene signature and corrected labels/numbers, writing to figures/pub/ WITHOUT
# clobbering the human Fig4 or the corrected Fig5/Fig6.
.libPaths(c("F:/Rlib", .libPaths()))
suppressMessages({ library(ggplot2); library(data.table); library(patchwork) })
th <- theme_bw(base_size=10) + theme(plot.title=element_text(size=8.4,face="bold",lineheight=1.05),
        plot.subtitle=element_text(size=7.2,lineheight=1.05),
        plot.tag=element_text(size=13,face="bold"), legend.key.size=unit(.4,"cm"),
        plot.margin=margin(6,8,4,4))
wrap <- function(s,n=42) paste(strwrap(s,width=n), collapse="\n")
sv <- function(p,f,w,h){ ggsave(file.path("figures/pub",f),p,width=w,height=h,dpi=300)
  ggsave(sub("\\.pdf$",".png",file.path("figures/pub",f)),p,width=w,height=h,dpi=150) }

## ---- Fig2: in-vitro osteogenic program (deduped 744) ----
sig <- fread("03_rra_meta/regeneration_competent_signature.tsv")  # canonical = 744 unique
ranked <- list.files("02_per_dataset_DE","_ranked.tsv$",full.names=TRUE)
topg <- head(sig$gene[order(sig$rra_score)],30)
hm <- sapply(ranked, function(f){ d<-fread(f); setNames(d$stat,d$gene)[topg] }); rownames(hm)<-topg
colnames(hm)<-sub("_ranked.tsv","",basename(ranked)); hm[is.na(hm)]<-0
hmm <- melt(as.data.table(hm,keep.rownames="gene"),id.vars="gene"); hmm$gene <- factor(hmm$gene, levels=rev(topg))
p2a <- ggplot(hmm, aes(variable,gene,fill=pmax(pmin(value,8),-8))) + geom_tile() +
  scale_fill_gradient2(low="#2166ac",mid="white",high="#b2182b",name="DE stat") +
  labs(title=wrap("Top-30 in-vitro osteogenic program genes across 5 dental sources",48),x=NULL,y=NULL) + th +
  theme(axis.text.x=element_text(angle=35,hjust=1,size=7), axis.text.y=element_text(size=6))
loo <- fread("03_rra_meta/loo_stability.tsv")
p2b <- ggplot(loo, aes(loo_fraction, fill=direction)) + geom_histogram(bins=12,position="identity",alpha=.6) +
  scale_fill_manual(values=c(up="#b2182b",down="#2166ac")) +
  labs(title=wrap("Leave-one-dataset-out stability (744 unique genes; 350 up / 394 down)",46),x="LOO retention fraction",y="genes") + th
gof <- "04_modules_wgcna/enrichment_signature_GO_Biological_Process_2023.tsv"
p2c <- if(file.exists(gof)){ go<-fread(gof)[1:8]; go$Term<-factor(sub(" \\(GO.*","",go$Term),levels=rev(sub(" \\(GO.*","",go$Term)))
  ggplot(go,aes(-log10(Adjusted.P.value),Term))+geom_col(fill="#4d9221")+labs(title=wrap("GO-BP enrichment (signature; note interferon/ISG cluster)",40),x="-log10 FDR",y=NULL)+th+theme(axis.text.y=element_text(size=7))
} else plot_spacer()
sv((p2a | (p2b/p2c)) + plot_annotation(tag_levels="A"), "Fig2_program.pdf", 12, 6)
cat("Fig2 regenerated (744 deduped, relabelled)\n")

## ---- Fig3: validation (deduped numbers) ----
ins <- fread("05_projection/osteogenic_score_insample.tsv")
p3a <- ggplot(ins,aes(group,osteo_score,fill=group))+geom_boxplot(alpha=.6,outlier.shape=NA)+geom_jitter(width=.15,size=1)+
  facet_wrap(~dataset,nrow=1,scales="free_y")+scale_fill_manual(values=c(control="#4575b4",osteo="#d73027"),guide="none")+
  labs(title=wrap("In-sample consistency (pooled AUC=1.0; not independent)",54),x=NULL,y="osteogenic score")+th+theme(strip.text=element_text(size=6))
ext <- fread("05_projection/external_validation_GSE316449.tsv")
p3b <- ggplot(ext,aes(group,osteo_score,fill=group))+geom_boxplot(alpha=.6,outlier.shape=NA)+geom_jitter(width=.12,size=2)+
  scale_fill_manual(values=c(control="#4575b4",osteo="#d73027"),guide="none")+
  labs(title=wrap("Out-of-sample in-vitro separability, GSE316449 (AUC=1.0, p=0.071; underpowered)",34),x=NULL,y="osteogenic score")+th
suc <- fread("05_projection/success_GSE104473_scores.tsv")
p3c <- ggplot(suc,aes(celltype,score,fill=celltype))+geom_boxplot(alpha=.6,outlier.shape=NA)+geom_jitter(aes(shape=condition),width=.15,size=1.6)+
  scale_fill_brewer(palette="Set2",guide="none")+
  labs(title=wrap("Does NOT transfer in vivo (GSE104473; rho=-0.14, p=0.56)",34),x="sorted skeletal cell",y="in-vitro osteogenic score")+th
sv((p3a/(p3b|p3c)) + plot_annotation(tag_levels="A"), "Fig3_validation.pdf", 11, 7)
cat("Fig3 regenerated (rho=-0.14/p=0.56, deduped external scores)\n")
