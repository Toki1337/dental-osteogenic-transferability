# R/real_run_mronj_isg_check.R — does the in-vitro program's myeloid-high / mesenchymal-low
# in-vivo localization come from its interferon-stimulated-gene (ISG) content? Re-score the
# 768 program with vs without ISGs and recompute per-compartment localization.
.libPaths(c("F:/Rlib", .libPaths()))
suppressMessages({ library(Seurat); library(data.table) })
OUT <- "06_scrna"
sc <- readRDS(file.path(OUT,"GSE303003_seurat.rds"))
sig <- fread("03_rra_meta/regeneration_competent_signature.tsv")
present <- function(g) g[g %in% rownames(sc)]

# ISG set (interferon/antiviral) — from the signature's GO 'defense response' leading edge
ISG <- c("MX1","MX2","ISG15","ISG20","IFIT1","IFIT2","IFIT3","OAS1","OAS2","OAS3","OASL",
         "IFI44","IFI44L","IFI6","IFI27","RSAD2","BST2","IRF7","STAT1","EPSTI1","HERC5",
         "USP18","XAF1","SAMD9","SAMD9L","DDX58","IFIH1","CMPK2","LY6E","PLSCR1")
up <- present(sig$gene[sig$direction=="up"]); dn <- present(sig$gene[sig$direction=="down"])
isg_in_prog <- intersect(union(up,dn), ISG)
cat("ISG genes present in the 768 program:", length(isg_in_prog), "->", paste(isg_in_prog, collapse=", "), "\n")

up_noisg <- setdiff(up, ISG); dn_noisg <- setdiff(dn, ISG)
sc <- AddModuleScore(sc, features=list(up_noisg), name="upN_", seed=1)
sc <- AddModuleScore(sc, features=list(dn_noisg), name="dnN_", seed=1)
sc$prog_net_noisg <- sc$upN_1 - sc$dnN_1

loc_full <- aggregate(prog_net      ~ compartment, data.frame(prog_net=sc$prog_net,       compartment=sc$compartment), mean)
loc_noisg<- aggregate(prog_net_noisg~ compartment, data.frame(prog_net_noisg=sc$prog_net_noisg, compartment=sc$compartment), mean)
loc <- merge(loc_full, loc_noisg, by="compartment")
loc <- loc[order(-loc$prog_net),]
fwrite(loc, file.path(OUT,"GSE303003_localization_ISG_check.tsv"), sep="\t")
cat("\n=== program score by compartment: full vs ISG-removed ===\n"); print(loc)

# correlation of per-compartment ranking, and the myeloid-minus-mesenchymal gap
gap_full  <- loc_full$prog_net[loc_full$compartment=="Myeloid"]       - loc_full$prog_net[loc_full$compartment=="Mesenchymal_osteo"]
gap_noisg <- loc_noisg$prog_net_noisg[loc_noisg$compartment=="Myeloid"] - loc_noisg$prog_net_noisg[loc_noisg$compartment=="Mesenchymal_osteo"]
cat(sprintf("\nMyeloid-minus-Mesenchymal program-score gap: full=%.4f  ISG-removed=%.4f  (%.0f%% reduction)\n",
            gap_full, gap_noisg, 100*(1-gap_noisg/gap_full)))

# re-test mesenchymal suppression with ISG-removed program (per-sample 4v4)
mes <- sc$compartment=="Mesenchymal_osteo"
md <- data.frame(sample=sc$sample[mes], condition=sc$condition[mes], s=sc$prog_net_noisg[mes])
ps <- aggregate(s ~ sample + condition, md, mean)
cat("\n=== ISG-removed program in mesenchymal compartment, per sample ===\n"); print(ps)
cat(sprintf("Wilcoxon MRONJ vs Cyst (ISG-removed, mesenchymal): p=%.3g; MRONJ_lower=%s\n",
            wilcox.test(s~condition, ps)$p.value,
            mean(ps$s[ps$condition=="MRONJ"]) < mean(ps$s[ps$condition=="Cyst"])))
cat("DONE ISG localization check.\n")
