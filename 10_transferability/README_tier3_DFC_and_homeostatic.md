# Tier-3 extensions — DFC deep-dive (A) and homeostatic reference (B)

Reviewer-requested additions. **A1 is computed** (below); **A2/A3 and B are scripted and
ready to run** but require the R + Seurat + `transferaudit` environment and the atlas
`.rds` objects — they were NOT re-run in the drafting session (the original R package
library lived on a drive that is no longer mounted). No un-run result is reported in the
manuscript; the manuscript states only the A1 numbers and points to the released scripts.

## A1 — DFC vs RRA-744 gene-content contrast (computed, in `dfc_vs_rra_contrast.txt`)
Source: `signatures_long.tsv` + `culture_decomposition.tsv` (no atlas needed).
- DFC = 394 genes (250 up / 144 down); shares only **193/394 (49%)** with RRA-744; **201 DFC-unique**; direction-concordant overlap 178, discordant 15 (Jaccard 0.20).
- **Osteoblast-identity content is essentially equal** — DFC 5.3% vs RRA-744 5.8% — so identity content is **not** the discriminator.
- DFC does **not** overlap the curated GO osteo sets more than RRA-744 (Jaccard to GO_osteoblast_diff 0.015 vs 0.028; to GO_ossification 0.029 vs 0.041).
- DFC carries **less culture contamination**: interferon 2.8% vs 5.9%, cell-cycle 0.25% vs 1.6%, senescence 1.5% vs 2.2% (hypoxia/serum similar); total culture load 11.2% vs 16.0%.
- But across the data-driven panel culture load does **not** predict TS (e.g. simvastatin transfers partially with 24% culture load; CBD has 0% culture load and does not transfer), consistent with the manuscript's ρ = 0.16 (ns).
- **Conclusion:** summary gene-content statistics do not explain why DFC alone transfers; the answer must be *which specific genes* it contains and whether they are expressed in the in-vivo osteoprogenitor compartment → resolved by the per-gene leading-edge step (A2/A3).
- Gene lists written: `dfc_genes_shared_with_rra.txt` (193), `dfc_genes_unique.txt` (201).

## A2 / A3 — leading-edge + correspondent subprogram — **COMPUTED** (`R/tier3_dfc_deepdive.R`)
Ran on the three atlases (R 4.6.1 + SeuratObject). Results:
- Leading-edge (`dfc_rra_leading_edge.tsv`): the top in-vivo-osteoprogenitor-expressed genes are
  **shared** between DFC and RRA-744 (COL3A1, COL11A1, GPX3, IGFBP5, SLIT3, PRELP, SERPING1, BMP2, OMD…).
- **In-vivo-correspondent subprogram** (`rra744_correspondent_subprogram_genes.txt`, 27 genes):
  ADAMTS1, AEBP1, **ALPL, BMP2, MGP, IGFBP5**, ASPN, C11orf96, CFH, COL11A1, COL3A1, CYP26B1, ENPP2,
  GPX3, IRX3, MFAP4, MFGE8, MME, MT2A, OMD, PMP22, PRELP, SCARA3, SERPING1, SLIT3, SRPX, TIMP1.
  = only **27 / 350 (~8%)** of RRA-744 up-genes reach canonical-marker-level osteoprogenitor
  expression (mean across the 3 atlases).
- **Re-projection** (`rra_correspondent_subprogram_projection.tsv`): the 27-gene subprogram localises
  to the osteoprogenitor compartment in **all three atlases (rank 1; Transferability z = 3.2 / 3.1 / 2.1;
  p = 5e-4 / 5e-4 / 1e-3)**, whereas the full RRA-744 does **not** (rank 4; p = 0.64 / 0.43 / 0.38);
  DFC_full and the positive control also localise (rank 1). So the meta-signature *contains* an
  in-vivo-correspondent core diluted ~13-fold by induction-context genes.
- **Caveat (stated in the manuscript):** the subprogram was defined on these atlases, so its
  re-localisation is expected by construction; the load-bearing, non-circular quantity is the DILUTION
  (only ~8% of a highly reproducible signature is osteoprogenitor-expressed in vivo).
Run log: `tier3_dfc_run.log`.

## B — homeostatic orofacial-bone reference — **COMPUTED** (`R/tier3_add_homeostatic_GSE316924.R`)
Downloaded **GSE316924** (GSM9462553, mouse mandible, 4,287 cells) to
`00_data/homeostatic/GSE316924/`, built a Seurat object, annotated compartments, and scored
the panel. Output: `invivo_localization_GSE316924_homeostatic.tsv`; run log
`tier3_homeostatic_run.log`. Result:
- **Gate PASSED** — the canonical-marker positive control localises to the osteolineage
  compartment (osteo_rank 1, Transferability z = 3.77, p = 0.0095), so the atlas is method-validated.
- **Provenance split reproduced in healthy bone:** all four prior-knowledge GO sets are
  significantly enriched toward the osteoprogenitor compartment (z = 1.8 / 4.8 / 4.9 / 7.1;
  p = 0.027 / 5e-4 / 5e-4 / 5e-4), while **none of the nine data-driven signatures is**
  (RRA-744 z = 1.07, p = 0.16; all data-driven p ≥ 0.16, including DFC here). So the
  non-correspondence is **not** an artifact of diseased/inflammatory tissue.
- Reported as an *additive reference-localisation panel* — single marrow-dominated sample; no
  differential stats; **not** merged into the 4/4-vs-1/9 disease-atlas Fisher test. Because no
  signature but the positive control ranked the osteoprogenitor compartment strictly first
  in this atlas, the calibrated null-based Transferability z (not the rank-1 criterion) is the readout.

## C — co-expression-preserving null (Major 3) — **COMPUTED** (`R/tier3_coexpr_null.R`)
Recomputed every Transferability Score against an **expression-matched null** (null gene sets drawn to
match each signature's per-gene mean-expression profile in 20 quantile bins → similar co-expression /
score variance as the real signature), addressing the concern that the random-gene-set null under-estimates
null variance and inflates |z|. Output: `invivo_localization_coexprnull.tsv`. Result:
- The mild p-inflation of the random null is confirmed (e.g. positive-control human-MRONJ p 0.019 → 0.025;
  GO_ossification ORNJ 5e-4 → 0.0015; DFC BRONJ 0.037 → 0.209) — but GO sets stay significant (p ≤ 0.036)
  and data-driven stay non-significant.
- **Split preserved and sharpened**: prior-knowledge **4/4**, data-driven **0/9** localising in all three
  atlases (the lone random-null exception, DFC, drops out — it fails BRONJ under the matched null);
  positive control still 3/3. RRA-744 still does not localise (z = −0.4 / −0.2 / +0.7).
- ⇒ the qualitative conclusion does not depend on the random-null calibration.

## D — AUCell cross-tool consistency (label-transfer item) — `R/tier3_aucell_check.R`
Cross-checks the transferaudit per-compartment score against an independent scorer (AUCell rank-AUC),
correlating per-compartment AUC with the M-based per-compartment score for every signature/atlas.
Output: `aucell_consistency.tsv` (complements the in-paper AddModuleScore-vs-UCell ρ = 1.0 with a third method).

## Not done (deliberately — see triage)
- GSE223778 (Mertk extraction-socket, bulk) — reintroduces compartment mixing; future-work mention only.
- DRA006607 (maxilla/mandible/ilium positional) — expands a deliberately-minimised negative side-thread; content unverified.
- Full ATAC/regulatory layer — separate follow-up. Note the in-vivo ATAC is already in hand (GSE104473); the proposed in-vitro "ATAC" GSE316447 is actually **H3K27me3 occupancy**, not ATAC.
