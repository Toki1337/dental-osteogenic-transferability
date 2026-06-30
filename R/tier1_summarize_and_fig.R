# R/tier1_summarize_and_fig.R
# ============================================================================
# Tier 1 synthesis: combine the scRNA-atlas localisation (now THREE jaw-disease
# atlases: GSE303003 human MRONJ, GSE295106 mouse BRONJ, GSE269255 mouse ORNJ)
# and the bulk lineage (GSE104473) into a per-signature TRANSFERABILITY SCORE +
# verdict, and emit the benchmark figure.
#
# An atlas is treated as positive-control-validated when the canonical osteoblast
# marker set localises to its osteoprogenitor compartment (rank 1 & p<0.05).
# Transferability Score (TS) = mean osteoprogenitor-localisation transferability_z
# across the validated scRNA atlases. Verdict "transfers" = osteoprogenitor rank 1
# AND p_localizes_more<0.05 in ALL validated scRNA atlases (conservative, replicated).
#
# Outputs: 10_transferability/transferability_summary.tsv
#          figures/pub/Fig7_transferability.{png,pdf}
# ============================================================================
.libPaths(c("F:/Rlib", .libPaths()))
suppressWarnings(suppressMessages({ library(data.table); library(ggplot2); library(patchwork) }))
setDTthreads(1)
OUT <- "10_transferability"
wrap <- function(s,n=60) paste(strwrap(s,width=n), collapse="\n")

loc <- fread(file.path(OUT, "invivo_localization.tsv"))
lin <- fread(file.path(OUT, "invivo_lineage_GSE104473.tsv"))
loc[, localizes := osteo_rank == 1 & p_localizes_more < 0.05]

# atlases, in display order, and which are positive-control-validated
atlas_order <- c("GSE303003_human_MRONJ","GSE295106_mouse_BRONJ","GSE269255_mouse_ORNJ")
atlas_order <- intersect(atlas_order, unique(loc$atlas))
pos <- loc[signature == "POSctrl_osteo_markers"]
validated <- pos[localizes == TRUE, atlas]
cat("positive-control-validated scRNA atlases:", paste(validated, collapse=", "), "\n")

sc <- dcast(loc, signature + source_type ~ atlas,
            value.var = c("transferability_z","osteo_rank","p_localizes_more","localizes"))
tzc  <- paste0("transferability_z_", validated)
locc <- paste0("localizes_", validated)
sc[, TS := rowMeans(as.matrix(.SD), na.rm = TRUE), .SDcols = tzc]
sc[, n_loc := rowSums(as.matrix(.SD) == TRUE, na.rm = TRUE), .SDcols = locc]
sc[, verdict := ifelse(n_loc == length(validated), "transfers",
                ifelse(n_loc > 0, "partial", "does_not_transfer"))]
sc <- merge(sc, lin[, .(signature, lineage_rho, lineage_TS = transferability_z)], by = "signature", all.x = TRUE)
setorder(sc, -TS)
fwrite(sc, file.path(OUT, "transferability_summary.tsv"), sep = "\t")

## ---- provenance + systematic-transfer statistic ---------------------------
sc[, prov := ifelse(source_type=="literature_GO","prior_knowledge",
             ifelse(source_type=="positive_control","positive_control","data_driven_in_vitro"))]
dd <- sc[prov=="data_driven_in_vitro"]; pk <- sc[prov=="prior_knowledge"]
ft <- fisher.test(matrix(c(sum(pk$verdict=="transfers"), sum(pk$verdict!="transfers"),
                           sum(dd$verdict=="transfers"), sum(dd$verdict!="transfers")), nrow=2))
cat(sprintf("\n%d validated scRNA atlases. Prior-knowledge transferring: %d/%d ; data-driven: %d/%d ; Fisher p=%.4f\n",
            length(validated), sum(pk$verdict=="transfers"), nrow(pk),
            sum(dd$verdict=="transfers"), nrow(dd), ft$p.value))
cat("\n==== Transferability summary ====\n")
print(sc[, .(signature, source_type, TS=round(TS,2), n_loc, n_val=length(validated), verdict)])

## ---- Figure 7 --------------------------------------------------------------
th <- theme_bw(base_size=10) + theme(plot.title=element_text(size=9.4,face="bold"),
        plot.subtitle=element_text(size=7.4,lineheight=1.05), plot.tag=element_text(size=13,face="bold"),
        axis.text.y=element_text(size=7.6), legend.key.size=unit(.4,"cm"))
prov_cols <- c(positive_control="#1a9850", prior_knowledge="#4575b4", data_driven_in_vitro="#d73027")
atlas_lab <- c(GSE303003_human_MRONJ="Human\nMRONJ", GSE295106_mouse_BRONJ="Mouse\nBRONJ",
               GSE269255_mouse_ORNJ="Mouse\nORNJ", lineage_TS="Mouse lineage\n(GSE104473)")

# panel A: per-atlas transferability_z heat (3 scRNA atlases + lineage)
ha <- melt(sc, id.vars=c("signature","prov"),
           measure.vars=c(paste0("transferability_z_", atlas_order), "lineage_TS"),
           variable.name="atlas", value.name="tz")
ha[, atlas_key := sub("transferability_z_", "", as.character(atlas))]
ha[, atlas := factor(atlas_lab[atlas_key], levels=atlas_lab[c(atlas_order,"lineage_TS")])]
ord <- sc[order(TS), signature]; ha[, signature := factor(signature, levels=ord)]
pA <- ggplot(ha, aes(atlas, signature, fill=tz)) + geom_tile(color="white", linewidth=.5) +
  geom_text(aes(label=sprintf("%.1f", tz)), size=2.6) +
  scale_fill_gradient2(low="#762a83", mid="grey95", high="#1b7837", midpoint=0,
        name="osteoprog.\nlocalisation z\n(vs random)") +
  labs(title="Osteoprogenitor localisation of each signature across in-vivo contexts",
       subtitle=wrap("Calibrated z vs size-matched permutation null (>0 = localises to the osteoprogenitor/mesenchymal compartment more than random). The three scRNA jaw-disease atlases are positive-control-validated; the bulk lineage is not (even canonical markers do not rise SSC->BCSP).",118),
       x=NULL, y=NULL) + th + theme(axis.text.x=element_text(size=7.4))
# panel B: combined Transferability Score, coloured by provenance, verdict
sc[, signature := factor(signature, levels=ord)]
pB <- ggplot(sc, aes(TS, signature, fill=prov)) + geom_col() +
  geom_vline(xintercept=0, linetype=2, color="grey50") +
  geom_text(aes(label=ifelse(verdict=="transfers","localises", ifelse(verdict=="partial","partial","")),
                x=pmax(TS,0)+0.1), hjust=0, size=2.5, color="grey30") +
  scale_fill_manual(values=prov_cols, name="signature provenance") +
  scale_x_continuous(expand=expansion(mult=c(0.03, 0.22))) +
  labs(title=sprintf("Transferability Score (mean osteoprogenitor-localisation z across %d validated scRNA atlases)", length(validated)),
       subtitle=wrap(sprintf("Prior-knowledge osteogenic gene sets and canonical markers localise to in-vivo osteoprogenitors in ALL validated atlases; data-driven in-vitro signatures (incl. the robust RRA-744 meta) do not (transfer by provenance: prior-knowledge %d/%d vs data-driven %d/%d, Fisher p=%.3f).",
                     sum(pk$verdict=="transfers"), nrow(pk), sum(dd$verdict=="transfers"), nrow(dd), ft$p.value), 118),
       x="Transferability Score (z vs random null)", y=NULL) + th
fig7 <- (pA / pB) + plot_layout(heights=c(1,1)) + plot_annotation(tag_levels="A")
ggsave("figures/pub/Fig7_transferability.pdf", fig7, width=9.5, height=11)
ggsave("figures/pub/Fig7_transferability.png", fig7, width=9.5, height=11, dpi=150)
cat("\nwrote figures/pub/Fig7_transferability + 10_transferability/transferability_summary.tsv\n")
