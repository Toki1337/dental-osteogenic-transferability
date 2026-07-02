# SUPERSEDED legacy scripts (archived)

These two scripts predate the corrected mouse->human ortholog mapping and retain the buggy
`babelgene::orthologs(..., human = TRUE)` call (plus an uppercase fallback). They are NOT part of
the released pipeline: `run_all.R` sources `R/step04_success_projection.R` and
`R/step08_position_context.R`, which use the corrected `utils::to_human_symbols(from='mouse')`
(deterministic 1:1, `human = FALSE`). Kept only for provenance; do not run.
