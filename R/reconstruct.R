# Multi-step state reconstruction from DMD operators.

# Suppress R CMD CHECK NOTEs for ggplot2 aes() column references.
utils::globalVariables(c("time_step", "bu_val", "trait_fac", "factor_name"))

#' Multi-Step DMD State Reconstruction
#'
#' Propagates state trajectories forward in time using the fitted DMD operator.
#' Two modes are available:
#'
#' * **Multi-step reconstruction** from an initial state \eqn{x_0}:
#'   \deqn{\hat{x}(1) = x_0, \quad \hat{x}(t+1) = A\,\hat{x}(t)}
#'   (discrete) or
#'   \deqn{\hat{x}(t+\Delta t) = e^{A_c \Delta t}\hat{x}(t)}
#'   (continuous).
#'   This produces long-range forecasts and is always returned as
#'   `X_pred_arr`.
#'
#' * **Single-step (one-ahead) reconstruction** from the *observed* state:
#'   \deqn{\hat{x}(t+1) = A\,x_{\text{obs}}(t)}
#'   (or its continuous-time analogue).
#'   This measures the one-step model residual and is returned as
#'   `X_pred_1step_arr` when `return_single_step = TRUE` and `X_obs` is
#'   provided.
#'
#' @param A_arr State-transition operator.  A 3-D array
#'   \eqn{[p \times p \times G]} or a single \eqn{[p \times p]} matrix
#'   (single-genotype convenience form).
#' @param x0 Initial state vector or matrix.  Accepted forms:
#'   \itemize{
#'     \item A numeric vector of length \eqn{p} -- same starting state for
#'       every genotype.
#'     \item A \eqn{[p \times G]} matrix -- per-genotype starting states.
#'   }
#'   When `NULL` (default), the first time-column of `X_obs` is used.
#' @param X_obs Optional observed state array \eqn{[p \times N \times G]}
#'   (or \eqn{[p \times N]} matrix for a single genotype).  Used to extract
#'   `x0` when `x0 = NULL`, and as the driving states for single-step
#'   reconstruction.
#' @param n_steps Integer number of prediction steps forward from `x0`.
#'   Defaults to `ncol(X_obs) - 1` when `X_obs` is provided.
#' @param return_single_step Logical.  If `TRUE` and `X_obs` is supplied,
#'   compute and return `X_pred_1step_arr`.  Default `FALSE`.
#' @param operator_type Character; `"discrete"` (default) or `"continuous"`.
#'   Set to `"continuous"` when `A_arr` is the continuous-time
#'   operator. The continuous-time propagation uses the matrix
#'   exponential: \deqn{\hat{x}(t+\Delta t) = e^{A_c \Delta t}\hat{x}(t).}
#'   Requires the \pkg{expm} package.
#' @param ts_days Numeric vector of time points (in days) of length
#'   \eqn{n\_steps + 1} for multi-step reconstruction or \eqn{N} for
#'   single-step reconstruction.  Used to compute per-step \eqn{\Delta t_k =
#'   ts\_days_{k+1} - ts\_days_k}.  Provide this **or** `dt` when
#'   `operator_type = "continuous"`.  Ignored for `"discrete"`.
#' @param dt Scalar uniform time step in days.  Used when `ts_days` is
#'   `NULL` and `operator_type = "continuous"`.  Ignored for `"discrete"`.
#'
#' @return A named list containing one or more of:
#' \describe{
#'   \item{`X_pred_arr`}{\eqn{[p \times (n\_steps+1) \times G]} multi-step
#'     reconstruction.  Column 1 equals `x0`.}
#'   \item{`X_pred_1step_arr`}{\eqn{[p \times (N-1) \times G]} single-step
#'     predictions (present only when `return_single_step = TRUE` and
#'     `X_obs` is supplied).}
#' }
#'
#' @seealso [compare_arrays()], [plot_trajectories()]
#' @export
reconstruct_dmd <- function(
  A_arr,
  x0                = NULL,
  X_obs             = NULL,
  n_steps           = NULL,
  return_single_step = FALSE,
  operator_type      = "discrete",
  ts_days            = NULL,
  dt                 = NULL
) {
  #  Promote 2-D / vector inputs to 3-D 
  A_arr <- .to_3d_arr(A_arr)
  X_obs <- .to_3d_arr(X_obs)


  p <- dim(A_arr)[1L]
  G <- dim(A_arr)[3L]

  geno_names <- dimnames(A_arr)[[3L]]
  if (is.null(geno_names)) geno_names <- paste0("G", seq_len(G))
  t_names <- dimnames(A_arr)[[1L]]
  if (is.null(t_names)) t_names <- paste0("T", seq_len(p))

  #  Build x0 matrix [p x G] 
  if (!is.null(x0)) {
    if (is.null(dim(x0))) {
      # 1-D vector: replicate across genotypes
      x0_mat <- matrix(rep(as.numeric(x0), G), nrow = p, ncol = G)
    } else {
      x0_mat <- as.matrix(x0)   # [p x G] or [p x 1]
      if (ncol(x0_mat) == 1L && G > 1L)
        x0_mat <- matrix(rep(x0_mat[, 1L], G), nrow = p, ncol = G)
    }
  } else if (!is.null(X_obs)) {
    x0_mat <- matrix(X_obs[, 1L, ], nrow = p, ncol = G)
  } else {
    stop("reconstruct_dmd: either x0 or X_obs must be provided.")
  }
  colnames(x0_mat) <- geno_names

  #  Determine n_steps 
  if (is.null(n_steps)) {
    n_steps <- if (!is.null(X_obs)) {
      dim(X_obs)[2L] - 1L
    } else if (has_ctrl) {
      dim(U_arr)[2L]
    } else {
      stop("reconstruct_dmd: n_steps must be specified when neither X_obs nor U_arr is provided.")
    }
  }
  n_steps <- as.integer(n_steps)

  # ── Continuous-time setup ─────────────────────────────────────────────────
  is_continuous <- identical(operator_type, "continuous")
  if (is_continuous) {
    if (!requireNamespace("expm", quietly = TRUE))
      stop("reconstruct_dmd: operator_type='continuous' requires the 'expm' package. ",
           "Install it with install.packages('expm').")
    if (is.null(ts_days) && is.null(dt))
      stop("reconstruct_dmd: operator_type='continuous' requires either 'ts_days' or 'dt'.")
    if (!is.null(ts_days)) {
      ts_days <- as.numeric(ts_days)
      # ts_days must cover (n_steps + 1) points
      if (length(ts_days) != n_steps + 1L)
        stop(sprintf(
          "reconstruct_dmd: ts_days must have length n_steps + 1 = %d, got %d.",
          n_steps + 1L, length(ts_days)))
      dt_vec <- diff(ts_days)   # length n_steps
    } else {
      dt_vec <- rep(as.numeric(dt), n_steps)
    }
    if (any(dt_vec <= 0))
      stop("reconstruct_dmd: all time steps (dt_vec) must be positive.")

    # Single-step ts_days (length N_obs) for 1-step reconstruction
    if (return_single_step && !is.null(X_obs)) {
      N_obs_ts <- dim(X_obs)[2L]
      if (!is.null(ts_days) && length(ts_days) == N_obs_ts) {
        dt_vec_1step <- diff(ts_days[seq_len(N_obs_ts)])
      } else if (!is.null(ts_days)) {
        # ts_days was for multi-step; fall back to scalar dt if available
        if (!is.null(dt)) {
          dt_vec_1step <- rep(as.numeric(dt), N_obs_ts - 1L)
        } else {
          warning("reconstruct_dmd: ts_days length doesn't match X_obs; ",
                  "using mean dt for single-step reconstruction.")
          dt_vec_1step <- rep(mean(dt_vec), N_obs_ts - 1L)
        }
      } else {
        dt_vec_1step <- rep(as.numeric(dt), N_obs_ts - 1L)
      }
    }
  }

  # Helper: compute the matrix exponential propagator.
  # Returns eAdT = expm(A*dT)
  .ct_step <- function(A_g, dT) {
    eAdT <- expm::expm(A_g * dT)
    list(eAdT = eAdT)
  }

  #  Multi-step reconstruction 
  X_pred <- array(NA_real_,
                  dim      = c(p, n_steps + 1L, G),
                  dimnames = list(t_names, NULL, geno_names))

  for (gi in seq_len(G)) {
    A_g <- matrix(A_arr[, , gi], nrow = p, ncol = p)
    X_pred[, 1L, gi] <- x0_mat[, gi]
    for (t in seq_len(n_steps)) {
      if (is_continuous) {
        cs      <- .ct_step(A_g, dt_vec[t])
        x_next  <- cs$eAdT %*% X_pred[, t, gi]
      } else {
        x_next <- A_g %*% X_pred[, t, gi]
      }
      X_pred[, t + 1L, gi] <- Re(as.vector(x_next))
    }
  }

  #  Single-step reconstruction (requires X_obs) 
  X_pred_1step <- NULL
  if (return_single_step && !is.null(X_obs)) {
    N_obs <- dim(X_obs)[2L]
    X_pred_1step <- array(NA_real_,
                          dim      = c(p, N_obs, G),
                          dimnames = list(t_names, NULL, geno_names))
    for (gi in seq_len(G)) {
      X_pred_1step[, 1L, gi] <- X_obs[, 1L, gi]
      A_g <- matrix(A_arr[, , gi], nrow = p, ncol = p)
      for (t in seq_len(N_obs - 1L) + 1) {
        if (is_continuous) {
          cs      <- .ct_step(A_g, dt_vec_1step[t])
          x_next  <- cs$eAdT %*% as.numeric(X_obs[, t - 1, gi])
        } else {
          x_next <- A_g %*% as.numeric(X_obs[, t - 1, gi])
        }
        X_pred_1step[, t, gi] <- Re(as.vector(x_next))
      }
    }
  }

  # Return only non-NULL elements
  out <- list(
    X_pred_arr       = X_pred,
    X_pred_1step_arr = X_pred_1step
  )
  Filter(Negate(is.null), out)
}
