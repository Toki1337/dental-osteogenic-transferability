# Single-cell RNA-seq reporting checklist

Reporting elements for the human MRONJ single-cell transferability audit
(GSE303003) and the supporting mouse atlas (GSE295106), following community
best practice (e.g., Füllgrabe et al. *Nat Biotechnol* 2020 MINSEQE-aligned
scRNA reporting; scverse/Seurat conventions).

| Element | GSE303003 (human MRONJ, primary) | GSE295106 (mouse BRONJ, secondary) |
|---|---|---|
| **Source / accession** | GSE303003, public 10x scRNA-seq | GSE295106, public 10x scRNA-seq |
| **Design / groups** | 4 MRONJ vs 4 odontogenic-cyst controls (biological replicates) | 1 BRONJ vs 1 Control animal (single-animal; demoted to directional corroboration only) |
| **Cells passing QC** | 65,454 | 22,991 |
| **QC thresholds** | >200 & <7,000 features/cell; <20% mitochondrial reads | same pipeline (`config/params.R`: sc_min_genes=200, sc_max_mt_pct=20) |
| **Normalization** | log-normalization (Seurat) | log-normalization |
| **Feature selection** | 2,000 highly variable features | 2,000 HVF |
| **Dimensionality reduction** | 30 PCs | 30 PCs |
| **Batch integration** | Harmony, integrated by sample | per-object processing |
| **Clustering** | Louvain, resolution 0.5 | resolution 0.5 |
| **Annotation** | 7 compartments by canonical human markers (COL1A1/PDGFRA/RUNX2/SP7/ALPL mesenchymal-osteo; CTSK/ACP5/MMP9 osteoclast; etc.) | 8 cell types incl. Osteo_MSC; mouse markers; mouse→human via babelgene 1:1 (`to_human_symbols`, human=FALSE) for signature scoring |
| **Compartment sizes** | Mesenchymal_osteo 13,991; Myeloid 6,719; Osteoclast 412; Tcell 8,339; Bcell_plasma 25,301; Endothelial 5,691; Epithelial 5,001 | Osteo_MSC 774; Macrophage 2,310; Neutrophil 10,069; Osteoclast 2,454; etc. |
| **Positive control (data captures disease)** | osteoclast fraction +~45× (p=0.029), myeloid >2× (p=0.029), Wilcoxon 4v4 | osteoclast expansion (single-animal, directional only) |
| **Signature scoring** | Seurat AddModuleScore + rank-based UCell (concordant, Spearman ρ=1); fast mean-z (genes×compartments) for the panel + permutation null | same fast mean-z + permutation null |
| **Pre-specified tests** | (i) compartment proportions; (ii) mesenchymal program score MRONJ vs ctrl; (iii) whole-program (744-gene) sign test | localization only (no inference, single animal) |
| **Multi-signature transferability** | 13-signature panel + positive control; osteoprogenitor-localization z vs size-matched permutation null (B=2,000); `10_transferability/` | same panel + null |
| **Statistical method (DE)** | mesenchymal pseudobulk DESeq2 (4v4) | not used for inference |
| **Replication** | mouse BRONJ atlas (GSE295106) reproduces the localization pattern (positive control + literature sets localize; data-driven signatures do not) | — |
| **Software / versions** | Seurat 5.5.1, Harmony 2.0.5, UCell 2.16.0, DESeq2 1.52.0 (`supp/sessionInfo.txt`) | same |
| **Code** | `R/real_run_mronj_human_projection.R`, `R/real_run_fig_mronj_human.R`, `R/tier1_project_invivo.R` | `R/real_run_step56_scrna.R`, `R/tier1_project_invivo.R` |
| **Limitations** | control = inflammatory cyst granulation tissue (not healthy bone); "not suppressed" stated relative to that baseline; MRONJ subtypes beyond anti-resorptive untested | single animal (n=1 vs 1); no statistical inference; CellChat used human DB on uppercased mouse symbols (illustrative) |
