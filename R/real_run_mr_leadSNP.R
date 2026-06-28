# R/real_run_mr_leadSNP.R
# CORRECTED cis-MR after referee R3 (M1/M2/M5): the previous "LD-clumped" run
# (real_run_harden.R) clumped against a reduced candidate-only panel and, on
# clump failure, passed ALL correlated SNPs into IVW -> impossible instrument
# counts (A2M nsnp=410) and spuriously inflated significance. We replace it with
# a single-lead-cis-SNP Wald-ratio (SMR-style) design that needs NO LD reference
# and is immune to that bug class. Palindromic (strand-ambiguous) SNPs are dropped
# (eQTLGen carries no allele frequency to resolve strand). Conservative by design.
.libPaths(c("F:/Rlib", .libPaths()))
suppressMessages({ library(data.table) })
setDTthreads(0)
G <- "00_data/genetic"
fl <- function(...) { cat(sprintf("[%s] %s\n", format(Sys.time(),"%H:%M:%S"), paste0(...))); flush.console() }

# ---- inputs ----
sig <- fread("03_rra_meta/regeneration_competent_signature.tsv"); cand <- unique(sig$gene)
fl("reading eQTLGen ..."); eq <- fread(file.path(G,"eqtlgen_cis_sig.txt.gz"))
eq <- eq[GeneSymbol %in% cand & Pvalue < 5e-8]
fl("instruments: ", nrow(eq), " across ", length(unique(eq$GeneSymbol)), " genes")
fl("reading eBMD ..."); eb <- fread(file.path(G,"ebmd_GCST006979.h.tsv.gz"),
  select=c("hm_rsid","hm_effect_allele","hm_other_allele","hm_beta","standard_error","hm_effect_allele_frequency","p_value","n"))
setnames(eb, c("rsid","ea_o","oa_o","beta_o","se_o","eaf_o","p_o","n_o"))
eb <- eb[!is.na(rsid) & !is.na(beta_o) & !is.na(eaf_o)]; setkey(eb, rsid)

# ---- Zhu et al. 2016 Z -> beta/se (using outcome EAF as MAF proxy) ----
zhu  <- function(z,p,N){ d<-2*p*(1-p)*(N+z^2); ifelse(d>0, z/sqrt(d), NA) }
zhse <- function(p,N,z){ d<-2*p*(1-p)*(N+z^2); ifelse(d>0, 1/sqrt(d), NA) }
is_palindromic <- function(a1,a2){ a<-paste0(toupper(a1),toupper(a2)); a %in% c("AT","TA","CG","GC") }

run_gene <- function(g, outcome=eb, label="eBMD"){
  e <- eq[GeneSymbol==g]; m <- merge(e, outcome, by.x="SNP", by.y="rsid")
  if(!nrow(m)) return(NULL)
  # allele match (same or flipped)
  match_same <- toupper(m$AssessedAllele)==toupper(m$ea_o) & toupper(m$OtherAllele)==toupper(m$oa_o)
  match_flip <- toupper(m$AssessedAllele)==toupper(m$oa_o) & toupper(m$OtherAllele)==toupper(m$ea_o)
  ok <- match_same | match_flip
  # DROP palindromic / strand-ambiguous SNPs (no eQTLGen EAF to resolve strand)
  palin <- is_palindromic(m$AssessedAllele, m$OtherAllele)
  m <- m[ok & !palin]; if(!nrow(m)) return(NULL)
  flip <- toupper(m$AssessedAllele)!=toupper(m$ea_o)
  eaf  <- ifelse(flip, 1-m$eaf_o, m$eaf_o)
  d <- data.frame(rsid=m$SNP, pval=m$Pvalue,
                  bx =zhu (m$Zscore, eaf, m$NrSamples),
                  bxse=zhse(eaf, m$NrSamples, m$Zscore),
                  by =ifelse(flip, -m$beta_o, m$beta_o), byse=m$se_o)
  d <- d[is.finite(d$bx)&is.finite(d$bxse)&is.finite(d$by)&is.finite(d$byse)&d$bxse>0&d$byse>0,]
  if(!nrow(d)) return(NULL)
  # LEAD cis-SNP = strongest instrument (smallest eQTLGen p)
  d <- d[order(d$pval),][1,,drop=FALSE]
  b  <- d$by/d$bx; se <- abs(d$byse/d$bx); p <- 2*pnorm(-abs(b/se))
  data.frame(gene=g, lead_snp=d$rsid, lead_eqtl_p=d$pval, method="Wald_leadSNP",
             b=b, se=se, pval=p, label=label, stringsAsFactors=FALSE)
}

genes <- unique(eq$GeneSymbol); fl("lead-cis-SNP Wald MR over ", length(genes), " genes ...")
res <- do.call(rbind, lapply(genes, function(g) tryCatch(run_gene(g), error=function(e) NULL)))
res <- res[!is.na(res$b),]
res$padj <- p.adjust(res$pval,"BH")
res$direction_signature <- sig$direction[match(res$gene, sig$gene)]
res <- res[order(res$pval),]
fwrite(res, "07_mr_coloc/mr_results_leadSNP.tsv", sep="\t")
fl("lead-SNP MR: ", nrow(res), " genes tested; FDR<0.05: ", sum(res$padj<0.05,na.rm=TRUE))
fl("top 12: ", paste(head(res$gene,12), collapse=", "))

# ---- known-bone-gene positive controls present among hits ----
known <- c("WNT16","MGP","SMAD3","HSPG2","KREMEN1","KLF12","SOST","LRP5","TNFRSF11B","SP7","RUNX2","SPP1","COL1A1")
kk <- res[res$gene %in% known, c("gene","b","pval","padj","direction_signature")]
fl("known bone genes recovered (any p): ", paste(kk$gene, collapse=", "))
print(kk)

# ---- re-derive druggable double-evidence under the corrected MR ----
core <- fread("04_modules_wgcna/core_program.tsv")
dg   <- fread("08_cmap/dgidb_core_program_drugs.tsv")
druggable_core <- unique(dg$gene_name[!is.na(dg$drug_name) & dg$drug_name!="" & dg$drug_name!="NULL"])
mr_sig <- res$gene[res$padj < 0.05]
de <- intersect(intersect(mr_sig, druggable_core), core$gene)
fl("FDR-sig MR genes: ", length(mr_sig))
fl("druggable core genes (DGIdb): ", length(intersect(druggable_core, core$gene)))
fl("DOUBLE-EVIDENCE (MR-FDR-sig AND druggable AND core): ",
   if(length(de)) paste(de, collapse=", ") else "NONE")
de_tab <- res[res$gene %in% de, ]
if(nrow(de_tab)) { de_tab$druggable<-TRUE; de_tab$in_core<-TRUE
  fwrite(de_tab, "09_integration/double_evidence_leadSNP.tsv", sep="\t") } else {
  fwrite(data.frame(note="No gene is simultaneously MR-FDR-significant (lead-SNP Wald), DGIdb-druggable, and in the 52-gene core program."),
         "09_integration/double_evidence_leadSNP.tsv", sep="\t") }

# ---- where did the provisional R2 genes land now? ----
prov <- c("GREM2","TLR4","LOXL2","EDNRB","PDGFD","IRS1","SRGN")
pj <- res[res$gene %in% prov, c("gene","lead_snp","b","pval","padj","direction_signature")]
fl("provisional R2 double-evidence genes under corrected MR:")
print(pj)
fwrite(res[res$gene %in% prov,], "07_mr_coloc/provisional_genes_leadSNP.tsv", sep="\t")
fl("DONE corrected lead-SNP MR.")
