# In-vivo correspondence checklist for in-vitro transcriptomic signatures

A reusable, actionable test to run **before** claiming that an in-vitro–derived
transcriptomic signature reflects an in-vivo cell state or has regenerative /
disease relevance. Cross-dataset reproducibility (RRA, leave-one-out, WGCNA, hub
connectivity) certifies *stability across in-vitro datasets* — it does **not**
certify *correspondence to an in-vivo cell identity*. These are separate claims
needing separate evidence.

Motivation: in this study a 744-gene, leave-one-dataset-out-stable, WGCNA-refined
pan-dental osteogenic signature passed every reproducibility filter yet did **not**
localise to in-vivo osteoprogenitors in either human MRONJ or mouse jaw single-cell
data — whereas prior-knowledge osteogenic gene sets and canonical markers did
(4/4 vs 1/9 data-driven; Fisher p = 0.007). The gap between "reproducible" and
"transferable" is the rule the checklist guards against.

## The checklist

1. **State the claim precisely.** Distinguish "stable across in-vitro datasets"
   (reproducibility) from "marks an in-vivo cell state" (correspondence) from
   "causal/therapeutic" (mechanism). Each needs its own evidence.

2. **Project onto ≥1 annotated in-vivo single-cell atlas** of the relevant tissue.
   Score the signature per cell (e.g., AddModuleScore + a rank-based method such as
   UCell; confirm the two agree) and per compartment.

3. **Include a positive control.** Score canonical identity markers of the target
   cell type the same way. If the positive control does **not** localise to the
   expected compartment in your atlas, the atlas/axis is non-discriminating — do
   not draw conclusions from it (e.g., an early SSC→BCSP bulk lineage where even
   canonical markers do not rise).

4. **Calibrate against a permutation null.** Compare the signature's
   target-compartment localisation to size-matched random gene sets
   (a *Transferability Score* = z vs null + empirical p). "Highest in compartment X"
   is meaningless without the null — many random sets rank some compartment first.
   (See the `transferaudit` R package.)

5. **Test ≥2 contexts / species where possible.** Replication across atlases
   (e.g., human disease + mouse model) guards against atlas-specific artifacts.

6. **Decompose confounds.** Quantify how much of the signature is culture-state
   programme (cell-cycle, serum/immediate-early, hypoxia, senescence, interferon)
   versus canonical identity genes. Data-driven induction signatures often capture
   induction-context transcription rather than cell identity.

7. **Report negatives as load-bearing.** A reproducible signature that does not
   localise in vivo is an honest, informative result — report it, scope the claim
   to where it holds, and do not infer in-vivo or therapeutic relevance from
   in-vitro reproducibility alone.

## One-line rule
> Reproducible ≠ transferable. Before attaching an in-vivo or regenerative claim to
> an in-vitro signature, show — with a positive control and a permutation null —
> that it localises to the expected in-vivo cell compartment.
