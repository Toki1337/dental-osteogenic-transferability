# R/real_run_mronj_robust.R — referee R5 M3+M4: re-do the human MRONJ projection on
# the DEDUPLICATED 744-gene signature, and harden the "no in-vivo osteoprogenitor
# home" claim with (1) UCell rank-based scoring as a second method, (2) a positive
# control that canonical osteo markers DO localize to the mesenchymal compartment
# (ruling out a depth/library-size artifact), (3) per-compartment Wilcoxon tests on
# per-sample means rather than bare global mean ranking, and (4) a RUNX2+/SP7+
# osteoprogenitor SUBcluster test within the mesenchymal compartment.
.libPaths(c("F:/Rlib", .libPaths()))
suppressMessages({ library(Seurat); library(data.table); library(UCell) })
set.seed(1)
OUT <- "06_scrna"
fl <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(),"%H:%M:%S"), paste0(...)))
sc <- readRDS(file.path(OUT,"GSE303003_seurat.rds"))
present <- function(g) g[g %in% rownames(sc)]

sig <- fread("03_rra_meta/regeneration_competent_signature_dedup.tsv")   # 744 unique
up <- present(sig$gene[sig$direction=="up"]); dn <- present(sig$gene[sig$direction=="down"])
fl("deduped signature mapped: ", length(up), " up, ", length(dn), " down")

# ---- dual scoring: AddModuleScore (already-style) + UCell (rank-based) ----
sc <- AddModuleScore(sc, features=list(up), name="dUp_", seed=1)
sc <- AddModuleScore(sc, features=list(dn), name="dDn_", seed=1)
sc$prog_AMS <- sc$dUp_1 - sc$dDn_1
sc <- AddModuleScore_UCell(sc, features=list(progUp=up, progDn=dn), name="_UC")
sc$prog_UC <- sc$progUp_UC - sc$progDn_UC
# positive-control osteo markers (should localize to mesenchymal if scoring is real)
osteo <- present(c("RUNX2","SP7","COL1A1","BGLAP","ALPL","IBSP","SPP1","COL1A2"))
sc <- AddModuleScore(sc, features=list(osteo), name="osteoMk_", seed=1)
sc$osteo_AMS <- sc$osteoMk_1
sc <- AddModuleScore_UCell(sc, features=list(osteo=osteo), name="_UC")   # creates column 'osteo_UC'

# ---- per-compartment localization (both methods) + positive control ----
agg <- function(v) tapply(v, sc$compartment, mean)
loc <- data.frame(compartment=names(agg(sc$prog_AMS)),
                  prog_AMS=as.numeric(agg(sc$prog_AMS)),
                  prog_UC =as.numeric(agg(sc$prog_UC)),
                  osteo_AMS=as.numeric(agg(sc$osteo_AMS)),
                  osteo_UC =as.numeric(agg(sc$osteo_UC)))
loc <- loc[order(-loc$prog_AMS),]
fwrite(loc, file.path(OUT,"GSE303003_localization_robust.tsv"), sep="\t")
fl("=== per-compartment scores: program (AMS,UCell) vs osteo-marker positive control ===")
print(loc)
# Spearman agreement of the two program scorings across compartments
fl("program AMS vs UCell per-compartment Spearman rho = ",
   round(cor(loc$prog_AMS, loc$prog_UC, method="spearman"),3))
fl("Does the osteo-marker positive control peak in Mesenchymal_osteo? ",
   loc$compartment[which.max(loc$osteo_AMS)], " (AMS), ", loc$compartment[which.max(loc$osteo_UC)], " (UCell)")

# ---- per-compartment proper test: per-sample mean, Wilcoxon MRONJ vs Cyst, both scores ----
cells <- data.frame(sample=sc$sample, condition=sc$condition, compartment=sc$compartment,
                    prog_AMS=sc$prog_AMS, prog_UC=sc$prog_UC)
ptab <- do.call(rbind, lapply(unique(cells$compartment), function(cp){
  d <- cells[cells$compartment==cp,]
  ps <- aggregate(cbind(prog_AMS,prog_UC) ~ sample+condition, d, mean)
  if(length(unique(ps$condition))<2) return(NULL)
  data.frame(compartment=cp, n_cells=nrow(d),
             AMS_Cyst=mean(ps$prog_AMS[ps$condition=="Cyst"]), AMS_MRONJ=mean(ps$prog_AMS[ps$condition=="MRONJ"]),
             AMS_p=tryCatch(wilcox.test(prog_AMS~condition,ps)$p.value,error=function(e)NA),
             UC_p =tryCatch(wilcox.test(prog_UC ~condition,ps)$p.value,error=function(e)NA))
}))
fwrite(ptab, file.path(OUT,"GSE303003_percompartment_MRONJ_test.tsv"), sep="\t")
fl("=== per-compartment MRONJ vs Cyst program test (per-sample Wilcoxon) ==="); print(ptab)

# ---- RUNX2+/SP7+/ALPL+ osteoprogenitor SUBcluster within mesenchymal ----
mes <- WhichCells(sc, expression = compartment=="Mesenchymal_osteo")
opmk <- present(c("RUNX2","SP7","ALPL","IBSP","BGLAP"))
expr <- GetAssayData(sc, layer="data")[opmk, mes, drop=FALSE]
is_op <- Matrix::colSums(expr>0) >= 2                      # >=2 osteoprogenitor markers detected
fl("mesenchymal cells: ", length(mes), "; RUNX2/SP7/ALPL/IBSP/BGLAP osteoprogenitor (>=2 markers+): ",
   sum(is_op), " (", round(100*mean(is_op),1), "%)")
opdf <- data.frame(sample=sc$sample[mes], condition=sc$condition[mes],
                   prog_AMS=sc$prog_AMS[mes], prog_UC=sc$prog_UC[mes], is_op=is_op)
if(sum(is_op)>=20){
  ops <- aggregate(cbind(prog_AMS,prog_UC) ~ sample+condition, opdf[opdf$is_op,], mean)
  fwrite(ops, file.path(OUT,"GSE303003_osteoprogenitor_subcluster.tsv"), sep="\t")
  fl("=== program in RUNX2+/SP7+ osteoprogenitor subcluster, per sample ==="); print(ops)
  fl("osteoprogenitor subcluster MRONJ vs Cyst: AMS Wilcoxon p=",
     signif(tryCatch(wilcox.test(prog_AMS~condition,ops)$p.value,error=function(e)NA),3),
     "  UCell p=", signif(tryCatch(wilcox.test(prog_UC~condition,ops)$p.value,error=function(e)NA),3))
  # does the program score higher in osteoprogenitors than in non-OP mesenchymal?
  fl("program AMS: osteoprogenitor mean=", round(mean(opdf$prog_AMS[opdf$is_op]),4),
     " vs non-OP mesenchymal mean=", round(mean(opdf$prog_AMS[!opdf$is_op]),4))
}

# ---- deduped sign test (save to TSV; M3) ----
suppressMessages(library(DESeq2))
mescells <- colnames(sc)[sc$compartment=="Mesenchymal_osteo"]
cnt <- GetAssayData(sc, layer="counts")[, mescells, drop=FALSE]
grp <- factor(sc$sample[sc$compartment=="Mesenchymal_osteo"])
pb <- t(rowsum(t(as.matrix(cnt)), grp)); cond <- ifelse(grepl("MRONJ",colnames(pb)),"MRONJ","Cyst")
pb <- pb[rowSums(pb)>=10,]
dds <- DESeq(DESeqDataSetFromMatrix(round(pb), data.frame(row.names=colnames(pb),
        condition=factor(cond,levels=c("Cyst","MRONJ"))), ~condition), quiet=TRUE)
res <- as.data.frame(results(dds, contrast=c("condition","MRONJ","Cyst"))); res$gene <- rownames(res)
st <- function(gu,gd,lab){ ru<-res[res$gene%in%gu&!is.na(res$log2FoldChange),]; rd<-res[res$gene%in%gd&!is.na(res$log2FoldChange),]
  s<-c(ru$log2FoldChange<0, rd$log2FoldChange>0); k<-sum(s); n<-length(s)
  data.frame(set=lab,k=k,n=n,frac=round(k/n,3),binom_p_greater=signif(binom.test(k,n,.5,alt="greater")$p.value,3)) }
core <- fread("04_modules_wgcna/core_program.tsv")
signtab <- rbind(st(up,dn,"signature_744_dedup"),
                 st(present(core$gene[core$direction=="up"]), present(core$gene[core$direction=="down"]), "core_52"))
fwrite(signtab, file.path(OUT,"GSE303003_program_signtest_dedup.tsv"), sep="\t")
fl("=== deduped whole-program sign test ==="); print(signtab)

# also per-sample mesenchymal program (deduped) + ISG-removed, save TSV
ISG <- c("MX1","MX2","ISG15","ISG20","IFIT1","IFIT2","IFIT3","OAS1","OAS2","OAS3","OASL","IFI44","IFI44L","IFI6","IFI27","RSAD2","BST2","IRF7","STAT1","EPSTI1","HERC5","USP18","XAF1","SAMD9","SAMD9L","DDX58","IFIH1","CMPK2","LY6E","PLSCR1")
sc <- AddModuleScore(sc, features=list(setdiff(up,ISG)), name="nUp_", seed=1)
sc <- AddModuleScore(sc, features=list(setdiff(dn,ISG)), name="nDn_", seed=1)
sc$prog_noISG <- sc$nUp_1 - sc$nDn_1
mm <- data.frame(sample=sc$sample[sc$compartment=="Mesenchymal_osteo"], condition=sc$condition[sc$compartment=="Mesenchymal_osteo"],
                 prog_AMS=sc$prog_AMS[sc$compartment=="Mesenchymal_osteo"], prog_UC=sc$prog_UC[sc$compartment=="Mesenchymal_osteo"],
                 prog_noISG=sc$prog_noISG[sc$compartment=="Mesenchymal_osteo"])
mps <- aggregate(cbind(prog_AMS,prog_UC,prog_noISG)~sample+condition, mm, mean)
fwrite(mps, file.path(OUT,"GSE303003_mes_program_persample_dedup.tsv"), sep="\t")
fl("mesenchymal per-sample (deduped): AMS p=", signif(wilcox.test(prog_AMS~condition,mps)$p.value,3),
   " UCell p=", signif(wilcox.test(prog_UC~condition,mps)$p.value,3),
   " ISG-removed p=", signif(wilcox.test(prog_noISG~condition,mps)$p.value,3))
fl("DONE robust re-analysis.")
