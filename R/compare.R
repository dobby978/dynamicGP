# Accuracy metrics between two 3-D state arrays.

# Internal RV coefficient (Robert & Escoufier 1976).
# X, Y: n x p and n x q matrices (rows = observations, columns = variables).
.rv_coef <- function(X, Y) {
  Xc    <- scale(X, center = TRUE, scale = FALSE)
  Yc    <- scale(Y, center = TRUE, scale = FALSE)
  VX    <- tcrossprod(Xc)
  VY    <- tcrossprod(Yc)
  denom <- sqrt(sum(VX * VX) * sum(VY * VY))
  if (!is.finite(denom) || denom == 0) return(NA_real_)
  sum(VX * VY) / denom
}

#' Compare Two 3-D State Arrays
#'
#' Computes prediction accuracy between an observed array (e.g. `X_arr`) and
#' a predicted array (e.g. `X_pred_arr` from
#' [reconstruct_dmd()]) along three complementary axes.
#'
#' @section Accuracy types:
#' \describe{
#'   \item{**Snapshot accuracy** (`$snapshot`)}{For each (trait, time-point)
#'     pair: the similarity *across genotypes* between observed and predicted
#'     values.  Metrics: Pearson correlation (PCC) with p-value and MSE.
#'     At least three genotypes are required for a non-trivial PCC.}
#'   \item{**Longitudinal accuracy** (`$longitudinal`)}{For each (trait,
#'     genotype) pair: the similarity *over time* between the observed and
#'     predicted trajectory.  Metrics: PCC, Spearman correlation (both with
#'     p-values), MSE, and optionally DTW distance.}
#'   \item{**RV by trait** (`$rv_by_trait`, optional)}{For each trait: the RV
#'     coefficient (Robert & Escoufier 1976) between the observed and predicted
#'     \eqn{[G \times N]} matrices -- a multivariate similarity measure.}
#'   \item{**RV by genotype** (`$rv_by_geno`, optional)}{For each genotype:
#'     the RV coefficient between the observed and predicted
#'     \eqn{[N \times p]} matrices.}
#' }
#'
#' @param X_arr Observed array \eqn{[p \times N \times G]}, or a
#'   \eqn{[p \times N]} matrix (single genotype, `G = 1`).
#' @param X_pred_arr Predicted array with the same `p` and `G` dimensions.
#'   If the `N` dimension differs, both arrays are truncated to the shorter
#'   length.
#' @param rv_by_trait Logical.  Compute per-trait RV coefficients.
#'   Default `FALSE`.
#' @param rv_by_geno Logical.  Compute per-genotype RV coefficients.
#'   Default `FALSE`.
#' @param use_dtw Logical.  Compute DTW distances using the \pkg{dtw} package
#'   (listed in `Suggests`).  Silently skipped with a message if the package
#'   is not installed.  Default `TRUE`.
#'
#' @return A named list with two mandatory elements (`snapshot`,
#'   `longitudinal`) and up to two optional elements:
#' \describe{
#'   \item{`snapshot`}{Data frame: `trait`, `time_step`, `n_obs`, `pcc`,
#'     `pcc_pval`, `mse`.}
#'   \item{`longitudinal`}{Data frame: `trait`, `genotype`, `n_obs`, `pcc`,
#'     `pcc_pval`, `spearman`, `spearman_pval`, `mse`, `dtw_dist` (normalised
#'     DTW distance in \eqn{[0, 1]}, divided by warping path length).}
#'   \item{`rv_by_trait`}{Named numeric vector (present when
#'     `rv_by_trait = TRUE`).}
#'   \item{`rv_by_geno`}{Named numeric vector (present when
#'     `rv_by_geno = TRUE`).}
#' }
#'
#' @references
#' Robert, P. & Escoufier, Y. (1976). A unifying tool for linear multivariate
#' statistical methods: the RV coefficient. *Applied Statistics*, 25(3),
#' 257--265. \doi{10.2307/2347233}
#'
#' Sakoe, H. & Chiba, S. (1978). Dynamic programming algorithm optimization
#' for spoken word recognition. *IEEE Transactions on Acoustics, Speech, and
#' Signal Processing*, 26(1), 43--49.
#'
#' @seealso [reconstruct_dmd()], [plot_trajectories()]
#' @export
compare_arrays <- function(X_arr,
                            X_pred_arr,
                            rv_by_trait = FALSE,
                            rv_by_geno  = FALSE,
                            use_dtw     = TRUE) {
  #  Promote 2-D inputs 
  X_arr      <- .to_3d_arr(X_arr)
  X_pred_arr <- .to_3d_arr(X_pred_arr)

  p <- dim(X_arr)[1L]
  G <- dim(X_arr)[3L]
  N <- min(dim(X_arr)[2L], dim(X_pred_arr)[2L])

  if (dim(X_pred_arr)[1L] != p)
    stop("compare_arrays: X_arr and X_pred_arr have different numbers of traits.")
  if (dim(X_pred_arr)[3L] != G)
    stop("compare_arrays: X_arr and X_pred_arr have different numbers of genotypes.")

  t_names <- dimnames(X_arr)[[1L]]
  if (is.null(t_names)) t_names <- paste0("T", seq_len(p))
  g_names <- dimnames(X_arr)[[3L]]
  if (is.null(g_names)) g_names <- paste0("G", seq_len(G))

  X_obs  <- X_arr[,      seq_len(N), , drop = FALSE]
  X_pred <- X_pred_arr[, seq_len(N), , drop = FALSE]

  #  DTW availability 
  dtw_ok <- use_dtw && requireNamespace("dtw", quietly = TRUE)
  if (use_dtw && !dtw_ok)
    message("compare_arrays: 'dtw' package not available; DTW distances set to NA.")
  if (G < 3L)
    message("compare_arrays: fewer than 3 genotypes -- snapshot PCC may be unreliable.")

  #  Snapshot accuracy [p x N rows] 
  snap_rows <- vector("list", p * N)
  k <- 0L
  for (i in seq_len(p)) {
    for (t in seq_len(N)) {
      k <- k + 1L
      obs_it  <- X_obs[i, t, ]
      pred_it <- X_pred[i, t, ]
      fin     <- is.finite(obs_it) & is.finite(pred_it)
      n_obs   <- sum(fin)
      pcc <- NA_real_; pval <- NA_real_
      if (n_obs >= 3L) {
        ct   <- stats::cor.test(obs_it[fin], pred_it[fin], method = "pearson")
        pcc  <- as.numeric(ct$estimate)
        pval <- ct$p.value
      }
      mse <- if (n_obs >= 1L) mean((obs_it[fin] - pred_it[fin])^2) else NA_real_
      snap_rows[[k]] <- data.frame(
        trait     = t_names[i],
        time_step = t,
        n_obs     = n_obs,
        pcc       = pcc,
        pcc_pval  = pval,
        mse       = mse,
        stringsAsFactors = FALSE
      )
    }
  }
  snapshot <- do.call(rbind, snap_rows)
  rownames(snapshot) <- NULL

  #  Longitudinal accuracy [p x G rows] 
  long_rows <- vector("list", p * G)
  k <- 0L
  for (i in seq_len(p)) {
    for (gi in seq_len(G)) {
      k <- k + 1L
      obs_ig  <- X_obs[i, , gi]
      pred_ig <- X_pred[i, , gi]
      fin     <- is.finite(obs_ig) & is.finite(pred_ig)
      n_obs   <- sum(fin)
      pcc  <- NA_real_; pcc_p  <- NA_real_
      spr  <- NA_real_; spr_p  <- NA_real_
      mse  <- NA_real_; dtw_d  <- NA_real_
      if (n_obs >= 3L) {
        ct_p  <- stats::cor.test(obs_ig[fin], pred_ig[fin], method = "pearson")
        pcc   <- as.numeric(ct_p$estimate); pcc_p <- ct_p$p.value
        ct_s  <- stats::cor.test(obs_ig[fin], pred_ig[fin], method = "spearman",
                                  exact = FALSE)
        spr   <- as.numeric(ct_s$estimate); spr_p <- ct_s$p.value
        mse   <- mean((obs_ig[fin] - pred_ig[fin])^2)
        if (dtw_ok) {
          dtw_res <- tryCatch(
            dtw::dtw(obs_ig[fin], pred_ig[fin], distance.only = TRUE),
            error = function(e) NULL
          )
          if (!is.null(dtw_res)) {
            path_len <- dtw_res$N + dtw_res$M - 1L
            dtw_d <- dtw_res$distance / max(path_len, 1L)
          }
        }
      }
      long_rows[[k]] <- data.frame(
        trait         = t_names[i],
        genotype      = g_names[gi],
        n_obs         = n_obs,
        pcc           = pcc,
        pcc_pval      = pcc_p,
        spearman      = spr,
        spearman_pval = spr_p,
        mse           = mse,
        dtw_dist      = dtw_d,
        stringsAsFactors = FALSE
      )
    }
  }
  longitudinal <- do.call(rbind, long_rows)
  rownames(longitudinal) <- NULL

  out <- list(snapshot = snapshot, longitudinal = longitudinal)

  #  RV by trait 
  if (rv_by_trait) {
    rv_t <- stats::setNames(
      vapply(seq_len(p), function(i) {
        X_mat <- t(X_obs[i, , ])   # [G x N]
        P_mat <- t(X_pred[i, , ])  # [G x N]
        ok    <- apply(X_mat, 1L, function(r) all(is.finite(r))) &
                 apply(P_mat, 1L, function(r) all(is.finite(r)))
        if (sum(ok) < 3L) return(NA_real_)
        .rv_coef(X_mat[ok, , drop = FALSE], P_mat[ok, , drop = FALSE])
      }, numeric(1L)),
      t_names
    )
    out$rv_by_trait <- rv_t
  }

  #  RV by genotype 
  if (rv_by_geno) {
    rv_g <- stats::setNames(
      vapply(seq_len(G), function(gi) {
        X_mat <- t(X_obs[, , gi])   # [N x p]
        P_mat <- t(X_pred[, , gi])  # [N x p]
        ok    <- apply(X_mat, 1L, function(r) all(is.finite(r))) &
                 apply(P_mat, 1L, function(r) all(is.finite(r)))
        if (sum(ok) < 3L) return(NA_real_)
        .rv_coef(X_mat[ok, , drop = FALSE], P_mat[ok, , drop = FALSE])
      }, numeric(1L)),
      g_names
    )
    out$rv_by_geno <- rv_g
  }

  out
}
