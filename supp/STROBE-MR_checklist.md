# STROBE-MR reporting checklist

Strengthening the Reporting of Observational Studies in Epidemiology — Mendelian
Randomization (STROBE-MR, Skrivankova et al. *BMJ* 2021). Item-by-item mapping to
where each is addressed. The MR layer of this study is explicitly an **exploratory
genetic-prioritization sanity check** (blood-eQTL exposure → systemic heel-eBMD
outcome; not jaw- or MRONJ-specific), not a confirmatory causal analysis.

| # | STROBE-MR item | Addressed |
|---|----------------|-----------|
| 1 | Title/abstract — MR design stated | Abstract: "druggable-genome Mendelian randomization (MR; blood cis-eQTL → heel-eBMD)"; framed exploratory. |
| 2 | Background/rationale | Methods §"Druggable-genome MR"; Intro notes the MR layer is *not* an MRONJ/jaw-specific causal analysis. |
| 3 | Objectives / hypotheses | Prioritization sanity check: do signature genes recover known bone-density biology and nominate druggable nodes. |
| 4 | Study design + MR assumptions (relevance, independence, exclusion-restriction) | Methods: cis-eQTL instruments (relevance); LD-clumping for independence; MR-Egger/weighted-median/MR-Lasso + Steiger for exclusion-restriction/pleiotropy (`07_mr_coloc/mr_sensitivity_suite.tsv`). |
| 5 | Data sources — exposure & outcome GWAS | Exposure: eQTLGen whole-blood cis-eQTL (n=31,684). Outcome: heel-eBMD GWAS (Morris & Kemp 2019, GCST006979). Replication: FinnGen M13_OSTEOPOROSIS. Table 1 + Data Availability. |
| 6 | Populations / ancestry overlap | Both predominantly European-ancestry; LD reference = 1000G EUR. Sample-overlap minimal (blood-eQTL vs UK-Biobank eBMD); noted as a residual limitation. |
| 7 | Exposure/outcome harmonization | Methods: allele harmonization; strand-ambiguous palindromic SNPs (A/T, C/G) dropped (eQTLGen lacks EAF); eQTLGen Z→β via Zhu et al. 2016 using eBMD EAF & eQTLGen N. |
| 8 | Instrument selection | cis-eQTL p<5×10⁻⁸; **primary** = single lead-cis-SNP Wald (SMR-style, no LD ref); **sensitivity** = full-1000G-EUR LD-clumped (r²<0.001, 10 Mb) IVW. |
| 9 | Assessment of assumptions / diagnostics | F-statistic implied by genome-wide cis instruments; MR-Egger intercept (directional pleiotropy), weighted-median, MR-Lasso (outlier-robust; MR-PRESSO analogue — MRPRESSO not installable), Steiger directionality. `07_mr_coloc/mr_sensitivity_suite.tsv`. |
| 10 | Statistical methods | Wald ratio (1 instrument), IVW (≥2), MR-Egger/weighted-median/MR-Lasso/mode-based (≥3); per-gene Benjamini–Hochberg FDR. |
| 11 | Sensitivity / robustness | Full suite run for all discussed genes (item 9). Colocalization NOT run (eQTLGen significant-only release lacks full cis regions) — disclosed. SMR-HEIDI not run — disclosed. |
| 12 | Effect measures | β = effect of higher blood expression on SD-units eBMD; explicitly **not** action-guiding (WNT16 β negative despite pro-bone role). |
| 13 | Descriptive: n instruments per gene | 490/744 genes with a usable lead-SNP estimate; per-gene post-clump nsnp reported (`mr_sensitivity_suite.tsv`, `mr_results_leadSNP.tsv`). |
| 14 | Main results | WNT16/KREMEN1/MGP/SMAD3/HSPG2 recovered (known bone genes); 105/490 FDR<0.05. Results §3.5. |
| 15 | Sensitivity results | GREM2 = directional-pleiotropy artifact (Egger intercept p=5.8×10⁻¹¹; IVW p=0.29, weighted-median p=0.24); IRS1 pleiotropic (Egger p=0.01); HSPG2 robust (IVW/wmedian/lasso p=8×10⁻⁹/6×10⁻⁸/3×10⁻¹²); TLR4/LOXL2/PDGFD null. Steiger: correct direction for all. |
| 16 | Replication | FinnGen clinical osteoporosis: 0 genes FDR<0.05, but eBMD-significant set directionally concordant above chance (62.9%, binomial p=0.005). |
| 17 | Limitations | Tissue mismatch (blood eQTL), phenotype mismatch (systemic eBMD, not jaw/MRONJ), single-instrument estimates lack per-gene pleiotropy test, no colocalization, underpowered FinnGen. Discussion + Table 1 footnote. |
| 18 | Interpretation | MR functions as a sanity check recovering known bone biology; **no robust druggable target nominated**. |
| 19 | Generalizability | Bounded to systemic bone density; jaw/MRONJ-specific GWAS would be required for context-matched inference. |
| 20 | Data/code availability | `07_mr_coloc/` outputs + `R/real_run_mr_*.R`, `R/tier2_mr_sensitivity_suite.R`; public GWAS/eQTL sources in Data Availability. |
