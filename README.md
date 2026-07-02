# Dental osteogenic program — in-vivo transferability audit

[![DOI](https://zenodo.org/badge/1286692569.svg)](https://doi.org/10.5281/zenodo.21126641)

Reproducible, public-data analysis package for the study:

> **A reproducible pan–dental-source in-vitro osteogenic program and a human single-cell
> audit of its in-vivo transferability in jaw regeneration and MRONJ.**

This repository contains the full analysis code, curation provenance, derived result
tables, and publication figures. It uses **only publicly available data** (GEO, GWAS
summary statistics, eQTLGen, open drug/connectivity resources); no new data were
generated and no wet-lab experiments were performed.

---

## Summary

1. **In-vitro osteogenic program.** A pan–dental-source osteogenic-induction signature is
   distilled by robust rank aggregation (RRA) across five cleanly curated dental-MSC
   datasets (PDLSC, DPSC, gingival GSC, SHED, DFC), with leave-one-dataset-out (LOO)
   stability and a WGCNA core. Result: **744 unique LOO-stable genes (350 up / 394 down)**
   and a **52-gene core**, recovering canonical osteogenic regulators (ALPL, BMP2, DLX5,
   MSX2 up; GREM2 down). It is separable out-of-sample *in vitro* (GSE316449; AUC = 1.0,
   but underpowered at n = 2 vs 6).

2. **Human single-cell transferability audit (the load-bearing result).** Projected onto
   human jaw single-cell data — **4 MRONJ vs 4 odontogenic-cyst controls (GSE303003,
   65,454 cells)** — the data reproduce MRONJ pathology (osteoclast/myeloid expansion,
   Wilcoxon p = 0.029), yet the program is **not suppressed** in the MRONJ mesenchyme
   (pre-specified whole-program sign test 43.9%, below chance; two scoring methods
   concordant, Spearman ρ = 1) and has **no in-vivo osteoprogenitor home** (it scores
   highest in myeloid/immune cells; a positive control confirms canonical osteo markers
   RUNX2/SP7/COL1A1 *do* peak in mesenchymal cells). This converges with a negative
   projection onto mouse mandibular distraction osteogenesis (GSE104473, ρ = −0.14).

3. **Exploratory genetics + repurposing.** Druggable-genome Mendelian randomization
   (lead-cis-SNP Wald + full-panel-clumped IVW/MR-Egger, eQTLGen → heel-eBMD; FinnGen
   osteoporosis replication) recovers known bone genes as positive controls but, after a
   pleiotropy-aware sensitivity analysis (GREM2 = directional-pleiotropy artifact),
   nominates **no robust druggable target**. Connectivity-map repurposing is a
   methodological sanity check (dexamethasone positive control + HDAC inhibitors).

### Scope and claim boundaries (please read before reuse)
- This is a **hypothesis-generating computational resource**, not a claim of therapeutic
  efficacy or mechanism.
- The take-home is a **cautionary, human-anchored lesson**: *in-vitro transcriptomic
  reproducibility does not imply in-vivo cellular correspondence.* We deliberately do
  **not** call the program "regeneration-competent," because the in-vivo audit does not
  support that.
- A reproducible interferon-stimulated-gene (ISG) culture signal persists through
  RRA + LOO + WGCNA (including the core hub *EPSTI1*) — disclosed, not hidden.
- The MR layer is anchored to **systemic** heel-eBMD via **blood** eQTL; it is explicitly
  **not** jaw- or MRONJ-specific, and the MR sign is not action-guiding.

---

## Repository layout

```
R/                 analysis scripts (RRA, WGCNA, scoring, single-cell, MR, CMap, integration)
config/            dataset registry (datasets.tsv), params, resources, druggable genome, drug-target map
01_curation/       human comparison-arm curation (sample_arms_final.tsv) + usable-dataset registry
02_per_dataset_DE/ … 09_integration/   derived result tables (small; reproduce the figures)
figures/pub/       publication figures Fig1–Fig6 (PNG + PDF); see figures/README.md
audit/             data-provenance audit (DATA_AUDIT.md) + GEO E-utilities verification
supp/              limitations & robustness audit; sessionInfo.txt (package versions)
env/               install_packages.R
manuscript/        data_availability.md (data/code availability statement)
run_all.R, Makefile   pipeline entry points
```

Raw downloaded data (`00_data/`), serialized objects (`*.rds`), and large GWAS/eQTL files
are **not** tracked; all are re-fetchable from the accessions in `config/datasets.tsv` and
`manuscript/data_availability.md`.

## Reproducing the analysis

1. **Environment.** R ≥ 4.6 with Bioconductor. Install packages via
   `Rscript env/install_packages.R` (exact versions used are recorded in
   `supp/sessionInfo.txt`).
2. **Curation gate.** `01_curation/sample_arms_final.tsv` (human-confirmed clean
   osteogenic-vs-control arms) and `rra_usable_datasets.tsv` are provided; `run_all.R`
   stops before the quantitative stages until these exist.
3. **Run.** `Rscript run_all.R` (or `make`) drives the pipeline. Individual stages can be
   run from `R/` (script names are self-describing).

> Note: scripts reference a local R library via `.libPaths("F:/Rlib")` for convenience on
> the original machine; adjust or remove for your environment.

## Data sources
All accessions and resources are listed in `config/datasets.tsv`, verified in
`audit/DATA_AUDIT.md`, and stated in `manuscript/data_availability.md`. Key datasets:
GSE99958, GSE226347, GSE266257, GSE49007 (training); GSE316449 (out-of-sample in-vitro);
GSE303003 (human MRONJ single-cell); GSE104473 (mouse distraction osteogenesis); GSE58474
(jaw position); heel-eBMD GWAS (GCST006979); eQTLGen; FinnGen; DGIdb; L1000CDS2; Enrichr.

## Citation
See [`CITATION.cff`](CITATION.cff). A permanent archived snapshot with a citable DOI will
be minted via Zenodo on release.

## Contact
Wendi Dou — Institute of Dentistry, I.M. Sechenov First Moscow State Medical University
(Sechenov University), Moscow, Russia.
ORCID: [0009-0008-2815-3726](https://orcid.org/0009-0008-2815-3726).

## License
Code is released under the MIT License (see [`LICENSE`](LICENSE)). Derived data tables and
figures may be reused under CC-BY 4.0 with attribution.
