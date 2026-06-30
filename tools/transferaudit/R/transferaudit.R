# transferaudit: in-vivo transferability audit of transcriptomic signatures
# ----------------------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

#' Build a genes x compartments standardized-mean expression matrix
#'
#' The per-compartment mean of a per-cell z-scored signature score equals the
#' signature-averaged per-compartment mean of each gene's z-scored expression.
#' This function precomputes that small genes x compartments matrix `M` from a
#' (sparse) log-normalized expression matrix, so signature scoring and the
#' permutation null reduce to fast row-subset column means.
#'
#' @param expr genes x cells matrix (base matrix or Matrix::dgCMatrix), log-normalized.
#' @param compartment character/factor vector of length ncol(expr): cell compartment labels.
#' @return matrix M (genes x compartments) of standardized mean expression; genes
#'   with zero variance are dropped.
#' @export
build_compartment_matrix <- function(expr, compartment) {
  compartment <- as.character(compartment)
  if (length(compartment) != ncol(expr)) stop("length(compartment) must equal ncol(expr)")
  mu  <- Matrix::rowMeans(expr)
  m2  <- Matrix::rowMeans(expr * expr)
  sdg <- sqrt(pmax(m2 - mu^2, 0))
  cl  <- sort(unique(compartment))
  cmean <- vapply(cl, function(c) as.numeric(Matrix::rowMeans(expr[, compartment == c, drop = FALSE])),
                  numeric(nrow(expr)))
  rownames(cmean) <- rownames(expr); colnames(cmean) <- cl
  M <- (cmean - mu) / sdg
  M[is.finite(rowSums(M)) & sdg > 0, , drop = FALSE]
}

# internal: per-compartment raw signature score (up minus down)
.score_comps <- function(M, up, down) {
  up <- intersect(up, rownames(M)); down <- intersect(down, rownames(M))
  if (length(up) < 3) return(NULL)
  s <- colMeans(M[up, , drop = FALSE])
  if (length(down) >= 3) s <- s - colMeans(M[down, , drop = FALSE])
  s
}

#' Transferability Score of one signature against a target compartment
#'
#' Scores each compartment for the signature, takes the target compartment's
#' z among compartments as the localisation statistic L, and calibrates it
#' against a size-matched random-gene-set permutation null. A signature that
#' genuinely marks the target compartment sits in the right tail (transferability
#' z >> 0, small `p_localizes_more`); a non-transferring signature does not.
#'
#' @param M genes x compartments matrix from \code{build_compartment_matrix}.
#' @param up,down character vectors of up- / down-regulated genes (down optional).
#' @param target name of the target compartment (a column of M).
#' @param n_perm number of permutations (default 2000).
#' @param seed RNG seed for reproducibility.
#' @param universe gene pool to sample the null from (default all rows of M).
#' @return a list with localization_z, target_rank, top_compartment,
#'   transferability_z, p_localizes_more, p_localizes_less, compartment_scores.
#' @export
transferability_score <- function(M, up, down = character(0), target,
                                   n_perm = 2000L, seed = 1L, universe = rownames(M)) {
  if (!target %in% colnames(M)) stop("target '", target, "' is not a compartment of M")
  s <- .score_comps(M, up, down)
  if (is.null(s)) stop("fewer than 3 up-genes map to M")
  z <- (s - mean(s)) / stats::sd(s)
  L <- z[[target]]; rk <- rank(-s)[[target]]
  nu <- length(intersect(up, universe)); nd <- length(intersect(down, universe))
  set.seed(seed); Ln <- numeric(n_perm)
  for (b in seq_len(n_perm)) {
    sb <- colMeans(M[sample(universe, nu), , drop = FALSE])
    if (nd >= 3) sb <- sb - colMeans(M[sample(universe, nd), , drop = FALSE])
    Ln[b] <- (sb[[target]] - mean(sb)) / stats::sd(sb)
  }
  list(localization_z = unname(L), target_rank = unname(rk), n_compartments = length(s),
       top_compartment = names(s)[which.max(s)],
       transferability_z = unname((L - mean(Ln)) / stats::sd(Ln)),
       p_localizes_more = (1 + sum(Ln >= L)) / (n_perm + 1),
       p_localizes_less = (1 + sum(Ln <= L)) / (n_perm + 1),
       compartment_scores = s)
}

#' Audit a panel of signatures against one target compartment
#'
#' @param M genes x compartments matrix.
#' @param signatures named list; each element a list with `up` and (optional) `down`.
#' @param target target compartment name.
#' @param n_perm,seed passed to \code{transferability_score}.
#' @return a data.frame, one row per signature.
#' @export
audit_panel <- function(M, signatures, target, n_perm = 2000L, seed = 1L) {
  rows <- lapply(names(signatures), function(nm) {
    s <- signatures[[nm]]
    r <- tryCatch(transferability_score(M, s$up, s$down %||% character(0), target, n_perm, seed),
                  error = function(e) NULL)
    if (is.null(r)) return(NULL)
    data.frame(signature = nm, target_rank = r$target_rank, n_compartments = r$n_compartments,
               top_compartment = r$top_compartment, localization_z = round(r$localization_z, 3),
               transferability_z = round(r$transferability_z, 3),
               p_localizes_more = signif(r$p_localizes_more, 3),
               transfers = r$target_rank == 1 & r$p_localizes_more < 0.05,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' Lineage transferability for a bulk ordinal lineage (e.g. sorted progenitors)
#'
#' Tests whether a signature score rises along an ordinal lineage, calibrated by
#' a size-matched permutation null on Spearman's rho.
#'
#' @param expr genes x samples log-expression matrix.
#' @param lineage_rank integer/numeric ordinal rank per sample (length ncol(expr)).
#' @param up,down signature gene vectors.
#' @param n_perm,seed permutation settings.
#' @return list(rho, transferability_z, p_rises_more).
#' @export
lineage_transferability <- function(expr, lineage_rank, up, down = character(0),
                                    n_perm = 2000L, seed = 1L) {
  z <- t(scale(t(expr))); z[is.na(z)] <- 0
  univ <- rownames(z)
  scr <- function(u, d) { u <- intersect(u, univ); d <- intersect(d, univ)
    if (length(u) < 3) stop("fewer than 3 up-genes map")
    s <- colMeans(z[u, , drop = FALSE]); if (length(d) >= 3) s <- s - colMeans(z[d, , drop = FALSE]); s }
  rho <- stats::cor(lineage_rank, scr(up, down), method = "spearman")
  nu <- length(intersect(up, univ)); nd <- length(intersect(down, univ))
  set.seed(seed); rn <- numeric(n_perm)
  for (b in seq_len(n_perm)) {
    sb <- colMeans(z[sample(univ, nu), , drop = FALSE])
    if (nd >= 3) sb <- sb - colMeans(z[sample(univ, nd), , drop = FALSE])
    rn[b] <- stats::cor(lineage_rank, sb, method = "spearman")
  }
  list(rho = unname(rho), transferability_z = unname((rho - mean(rn)) / stats::sd(rn)),
       p_rises_more = (1 + sum(rn >= rho)) / (n_perm + 1))
}
