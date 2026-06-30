# Figures — provenance (current Fig 1–8 system)

The submission figure set is **`figures/pub/Fig1`–`Fig8`** (PNG + PDF). Everything in
`figures/legacy/` is superseded and not referenced by the manuscript.

| Figure | Content | Generating script |
|---|---|---|
| **Fig1_framework** | Conceptual schematic: in-vitro program → human single-cell transferability audit → exploratory genetic/repurposing, with honest-scope footer | `R/real_run_fig1_framework.R` (programmatic, not hand-drawn) |
| **Fig2_program** | Pan-dental in-vitro osteogenic program: cross-dataset DE of top RRA genes, LOO stability (744 unique genes; 350 up / 394 down), GO-BP enrichment (incl. interferon/ISG cluster) | `R/real_run_fig2fig3_corrected.R` |
| **Fig3_validation** | In-sample consistency (AUC=1.0, not independent); out-of-sample in-vitro separability on GSE316449 (underpowered, p=0.071); negative in-vivo projection on GSE104473 (rho=−0.14, p=0.56) | `R/real_run_fig2fig3_corrected.R` |
| **Fig4_failure_human** | Human single-cell transferability audit in MRONJ (GSE303003, 4 vs 4; 65,454 cells): compartment UMAP; replicated osteoclast/myeloid expansion (p=0.029); program localization vs osteo-marker positive control (UCell+AddModuleScore, Spearman ρ=1); program not suppressed (sign test 43.9%) | `R/real_run_fig_mronj_human.R` |
| **Fig5_MR** | Exploratory lead-cis-SNP Wald MR on heel-eBMD; recovers known bone genes (WNT16/KREMEN1/MGP/SMAD3/HSPG2) as positive controls; not jaw/MRONJ-specific | `R/real_run_fig5fig6_corrected.R` |
| **Fig6_priority** | (a) CMap 4-method convergence (dexamethasone + HDAC inhibitors); (b) integrated priority; no robust genetic double-evidence (GREM2 pleiotropy artifact; only weak SRGN); tiers labelled "eBMD-MR-associated" (not "causal") | `R/real_run_fig5fig6_corrected.R` |
| **Fig7_transferability** | Multi-signature transferability benchmark: (a) osteoprogenitor-localisation z across 3 in-vivo contexts; (b) Transferability Score by provenance (prior-knowledge 4/4 vs data-driven 1/9 localise; Fisher p=0.007) | `R/tier1_summarize_and_fig.R` (panel from `R/tier1_build_signatures.R`, `R/tier1_mine_excluded.R`, `R/tier1_project_invivo.R`, `R/tier1_project_lineage.R`) |
| **Fig8_culture_confound** | Culture-state decomposition: (a) culture-program composition per signature; (b) hypergeometric enrichment (RRA-744 → hypoxia/serum/IFN; perturbation → cell-cycle); data-driven ≈5% osteo-identity vs ≈92% curated | `R/tier2_culture_decomposition.R` |

**Honesty note.** Fig 1 is a conceptual diagram (drawn programmatically, no invented data points).
Fig 2–6 are emitted by the scripts above from the public data; none are mocked. Earlier
numbering (Fig7/Fig8, mouse-BRONJ Fig4/Fig5, distance-pruned MR forest) is retained only under
`figures/legacy/` and is **not** part of the manuscript.
