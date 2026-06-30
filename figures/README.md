# Figures — provenance (current Fig 1–8 system)

The submission figure set is **`figures/pub/Fig1`–`Fig8`** (PNG + PDF). Everything in
`figures/legacy/` is superseded and not referenced by the manuscript.

| Figure | Content | Generating script |
|---|---|---|
| **Fig1_framework** | Conceptual schematic: in-vitro program → human single-cell transferability audit → exploratory genetic/repurposing, with honest-scope footer | `R/real_run_fig1_framework.R` (programmatic, not hand-drawn) |
| **Fig2_program** | Pan-dental in-vitro osteogenic program: cross-dataset DE of top RRA genes, LOO stability (744 unique genes; 350 up / 394 down), GO-BP enrichment (incl. interferon/ISG cluster) | `R/real_run_fig2fig3_corrected.R` |
| **Fig3_validation** | In-sample consistency (AUC=1.0, not independent); out-of-sample in-vitro separability on GSE316449 (underpowered, p=0.071); negative in-vivo projection on GSE104473 (rho=−0.14, p=0.56) | `R/real_run_fig2fig3_corrected.R` |
| **Fig4_failure_human** | Human single-cell transferability audit in MRONJ (GSE303003, 4 vs 4; 65,454 cells): compartment UMAP; replicated osteoclast/myeloid expansion (p=0.029); program localization vs osteo-marker positive control (UCell+AddModuleScore, Spearman ρ=1); program not suppressed (sign test 43.9%) | `R/real_run_fig_mronj_human.R` |
| **Fig5_transferability** | Multi-signature transferability benchmark: (a) osteoprogenitor-localisation z across 4 in-vivo contexts (3 validated scRNA atlases + lineage); (b) Transferability Score by provenance (prior-knowledge 4/4 vs data-driven 1/9 localise; Fisher p=0.007). *First cited in §3.5.* | `R/tier1_summarize_and_fig.R` (panel from `R/tier1_build_signatures.R`, `R/tier1_mine_excluded.R`, `R/tier1_atlas_GSE269255.R`, `R/tier1_project_invivo.R`, `R/tier1_project_lineage.R`) |
| **Fig6_culture_confound** | Culture-state decomposition: (a) culture-program composition per signature; (b) hypergeometric enrichment (RRA-744 → hypoxia/serum/IFN; perturbation → cell-cycle); data-driven ≈5% osteo-identity vs ≈92% curated. *Cited in §3.5.* | `R/tier2_culture_decomposition.R` |
| **Fig7_MR** | Exploratory lead-cis-SNP Wald MR on heel-eBMD; recovers known bone genes (WNT16/KREMEN1/MGP/SMAD3/HSPG2) as positive controls; not jaw/MRONJ-specific. *Cited in §3.6.* | `R/real_run_fig5fig6_corrected.R` |
| **Fig8_priority** | (a) CMap 4-method convergence (dexamethasone + HDAC inhibitors); (b) integrated priority; no robust genetic double-evidence (GREM2 pleiotropy artifact; only weak SRGN). *Cited in §3.7.* | `R/real_run_fig5fig6_corrected.R` |

**Numbering note.** Figures are numbered in **citation order**. The multi-signature benchmark
(Fig5) and culture decomposition (Fig6) are cited in §3.5, which precedes the exploratory MR
(Fig7, §3.6) and CMap/priority (Fig8, §3.7) — so the benchmark/culture figures come *before* the
MR/CMap figures even though they were added later. Filenames encode the content (`Fig5_transferability`,
`Fig7_MR`, …) so a rename is needed if the section order ever changes.

**Honesty note.** Fig 1 is a conceptual diagram (drawn programmatically, no invented data points).
Fig 2–8 are emitted by the scripts above from the public data; none are mocked. A superseded
earlier figure scheme (mouse-BRONJ panels, distance-pruned MR forest) is retained only under
`figures/legacy/` and is **not** part of the manuscript.
