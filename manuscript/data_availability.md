# Data and Code Availability

This study analysed **only publicly available data**; no new data were generated and no human/animal experiments were performed.

## Transcriptomic datasets (GEO; verified 2026-06, see `audit/DATA_AUDIT.md`)
Success side — dental-MSC osteogenic differentiation (candidate training datasets; final clean-arm count reported from `01_curation/rra_usable_datasets.tsv`, with the audit predicting ~11–12 usable arms after per-arm curation): GSE99958, GSE163354, GSE226347, GSE296018, GSE271641, GSE299041, GSE286540, GSE236009, GSE266150, GSE159507, GSE266257, GSE49007 (+ GSE159508 miRNA as side evidence; GSE160451 as negative control; GSE105145 as cross-source anchor). External osteogenic validation: GSE316449/GSE316447 (α-KG, RNA-seq + ATAC).
Success side — jaw regeneration: GSE104473 (mandibular distraction osteogenesis, RNA-seq + ATAC), GSE223778 (socket healing).
Failure side — MRONJ/ORNJ: **GSE303003 (human jaw scRNA, 4 MRONJ vs 4 odontogenic-cyst control — PRIMARY failure-side transferability audit; 65,454 cells)**; GSE295106 (mouse mandibular marrow, control vs BRONJ scRNA — secondary, single-animal directional corroboration only); GSE269255 (mouse mandibular irradiation/ORNJ time-course, Control/Day1/Day7 — processed into a compartment-annotated atlas and used as the third jaw-disease single-cell context in the transferability benchmark; `R/tier1_atlas_GSE269255.R`). Supportive only: GSE7116 (multiple-myeloma ONJ peripheral blood; confound acknowledged).
Jaw position-specific background: GSE58474 (mandibular vs iliac osteoblasts), GSE30167 (jaw/alveolar vs long bone).

## Raw-sequence submission (DDBJ/NCBI-SRA) — considered but NOT analysed
DRA006607 — maxilla/mandible/iliac MSC position-memory (NCBI SRA uid 10220584); raw FASTQ requiring alignment/quantification. This record was identified during dataset scoping but was **not** analysed in the present pipeline (no alignment/quantification was performed); it is listed only for completeness and is not the basis of any reported result.

## Genetic, perturbation, and annotation resources
- GWAS: heel-eBMD (Morris & Kemp 2019; GWAS Catalog GCST006979); GEFOS DXA-BMD; FinnGen (fracture/osteoporosis endpoints).
- eQTL: eQTLGen (blood cis-eQTL, n=31,684); GTEx v8.
- Druggability/targets: Finan et al. 2017 druggable genome; Open Targets; DGIdb; ChEMBL; DrugBank.
- Perturbation: CMap/LINCS L1000 (clue.io; account required); L1000CDS2; L1000FWD (open APIs).
- Gene sets/networks: MSigDB; WikiPathways; STRING v12.

## Code
All analysis code (`R/step00`–`step11`), central parameters (`config/params.R`), the dataset registry (`config/datasets.tsv`), resource map (`config/resources.tsv`), and the data-provenance audit (`audit/`) are available at https://github.com/Toki1337/dental-osteogenic-transferability (a citable Zenodo DOI will be minted on release). The pipeline is driven by `run_all.R` / `Makefile`; dependencies install via `env/install_packages.R`; the exact package versions are recorded in `supp/sessionInfo.txt`.

The multi-signature transferability benchmark is implemented in `R/tier1_build_signatures.R`, `R/tier1_mine_excluded.R`, `R/tier1_project_invivo.R`, `R/tier1_project_lineage.R`, `R/tier1_summarize_and_fig.R` and `R/tier2_culture_decomposition.R`, writing to `10_transferability/`; its calibrated Transferability Score is released as a standalone, installable R package with a vignette (`tools/transferaudit/`). The MR sensitivity suite (`R/tier2_mr_sensitivity_suite.R` → `07_mr_coloc/mr_sensitivity_suite.tsv`), a STROBE-MR checklist, a single-cell reporting checklist, and a reusable in-vivo-correspondence checklist are in `supp/`.

**Honesty note.** The released code computes all reported quantities on the public data above. No result, figure, or statistic in this project was fabricated or hand-entered; every numeric output is produced by running the released scripts.
