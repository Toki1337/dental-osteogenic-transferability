# transferaudit

**In-vivo transferability audit of transcriptomic signatures.**

A reproducible in-vitro transcriptomic signature need not correspond to any
discrete in-vivo cell state. `transferaudit` makes that question quantitative and
calibrated: given a signature and an annotated single-cell reference atlas, it asks
whether the signature localises to its expected compartment **more than a
size-matched random gene set would**, and returns a **Transferability Score**
(localisation z vs a permutation null).

It was built for, and is the analysis engine of, a dental mesenchymal-stromal-cell
osteogenic program audited across human MRONJ (GSE303003) and mouse jaw atlases,
but it is generic: any signature × compartment × atlas.

## Why
Cross-dataset reproducibility filters (RRA, leave-one-out, co-expression modules)
certify that a signature is *stable*, not that it marks an *in-vivo cell identity*.
In a 13-signature benchmark, prior-knowledge osteogenic gene sets and canonical
markers localised to in-vivo osteoprogenitors (4/4), while data-driven in-vitro
signatures — including a robust leave-one-out-stable meta-signature — systematically
did not (1/9; Fisher p = 0.007). `transferaudit` is the calibrated test that
distinguishes the two.

## Install
```r
# from a local clone of the repository
install.packages("Matrix")
# install.packages("remotes"); remotes::install_local("tools/transferaudit")
devtools::load_all("tools/transferaudit")   # or build & INSTALL
```
Depends only on `Matrix` and `stats`.

## Quick start
```r
library(transferaudit)

# expr: genes x cells log-normalized matrix (dense or Matrix::dgCMatrix)
# comp: per-cell compartment labels (e.g. from a Seurat object's metadata)
M <- build_compartment_matrix(expr, comp)

# one signature against the osteoprogenitor compartment
res <- transferability_score(M, up = my_up_genes, down = my_down_genes,
                             target = "Mesenchymal_osteo", n_perm = 2000)
res$transferability_z      # z vs random null  (>0 = localises more than random)
res$p_localizes_more       # empirical one-sided p

# a whole panel at once
panel <- list(
  my_signature = list(up = my_up_genes, down = my_down_genes),
  positive_ctrl = list(up = c("RUNX2","SP7","COL1A1","BGLAP","ALPL"))
)
audit_panel(M, panel, target = "Mesenchymal_osteo")

# bulk ordinal lineage (e.g. sorted SSC -> BCSP)
lineage_transferability(bulk_expr, lineage_rank = c(1,1,2,2), up = my_up_genes)
```

See `vignette("transferaudit")` for a fully self-contained worked example.

## The Transferability Score
For a signature with `n_up`/`n_down` genes and an annotated atlas:
1. score each compartment (mean z of up-genes minus down-genes, from `M`);
2. take the target compartment's z among compartments as the localisation `L`;
3. draw `n_perm` size-matched random gene sets, recompute `L`;
4. **Transferability Score** = `(L - mean(L_null)) / sd(L_null)`; empirical p from the null.

A signature "transfers" when it ranks the target compartment first **and**
`p_localizes_more < 0.05`.

## Citation
See `inst/CITATION` (`citation("transferaudit")`). License: MIT.
