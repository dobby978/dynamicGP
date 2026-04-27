# Multi-step state reconstruction from DMD operators.

# Suppress R CMD CHECK NOTEs for ggplot2 aes() column references.
utils::globalVariables(c("time_step", "bu_val", "trait_fac", "factor_name"))

#' Multi-Step DMD State Reconstruction
#'
#' Propagates state trajectories forward in time using fitted DMD operators.
#' Two modes are available:
#'
#' * **Multi-step reconstruction** from an initial state \eqn{x_0}:
#'   \deqn{\hat{x}(1) = x_0, \quad \hat{x}(t+1) = A\,\hat{x}(t) + B\,u(t)}
#'   (discrete) or
#'   \deqn{\hat{x}(t+\Delta t) = e^{A_c \Delta t}\hat{x}(t) +
#'         A_c^{-1}(e^{A_c \Delta t}-I)\,B_c\,u(t)}
#'   (continuous, zero-order-hold for \eqn{B}).
#'   This produces long-range forecasts and is always returned as
#'   `X_pred_arr`.
#'
#' * **Single-step (one-ahead) reconstruction** from the *observed* state:
#'   \deqn{\hat{x}(t+1) = A\,x_{\text{obs}}(t) + B\,u(t)}
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
#' @param B_arr Optional control-input operator \eqn{[p \times q \times G]},
#'   or a single \eqn{[p \times q]} matrix that is broadcast to all `G`
#'   genotypes.
#' @param U_arr Optional control-input array \eqn{[q \times N_u \times G]},
#'   or a single \eqn{[q \times N_u]} matrix that is broadcast to all `G`
#'   genotypes.  Must be provided together with `B_arr`.
#' @param n_steps Integer number of prediction steps forward from `x0`.
#'   Defaults to `ncol(X_obs) - 1` when `X_obs` is provided, or
#'   `ncol(U_arr)` when only `U_arr` is provided.
#' @param return_single_step Logical.  If `TRUE` and `X_obs` is supplied,
#'   compute and return `X_pred_1step_arr`.  Default `FALSE`.
#' @param plot_control_influence Logical.  If `TRUE` (default) and both
#'   `B_arr` and `U_arr` are supplied, include heatmap(s) of the per-factor
#'   signed control influence \eqn{B_{:,j}\,u_j(t)} per trait and time step
#'   as `ctrl_plot`.
#' @param ctrl_factors Optional character vector of control-factor names or
#'   integer indices selecting which columns of \code{U}/\code{B} to plot.
#'   \code{NULL} (default) plots all \code{q} factors.  Ignored when
#'   \code{plot_control_influence = FALSE}.
#' @param facet_by_factor Logical.  When multiple factors are selected,
#'   combine all factor heatmaps into a single \pkg{ggplot2} plot faceted by
#'   factor name (\code{TRUE}).  The default (\code{FALSE}) returns a named
#'   list of individual heatmaps, one per factor.  Ignored when
#'   \code{plot_control_influence = FALSE} or only one factor is selected.
#' @param operator_type Character; `"discrete"` (default) or `"continuous"`.
#'   Set to `"continuous"` when `A_arr` and `B_arr` are the continuous-time
#'   operators returned by `method = "continuous"` in
#'   [run_dmd()].  The continuous-time propagation uses the matrix
#'   exponential with a zero-order-hold (ZOH) approximation for `B`:
#'   \deqn{\hat{x}(t+\Delta t) = e^{A_c \Delta t}\hat{x}(t) +
#'         A_c^{-1}(e^{A_c \Delta t}-I)\,B_c\,u(t).}
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
#'   \item{`ctrl_plot`}{Control-influence heatmap(s) for the selected
#'     factor(s).  When a single factor is selected or
#'     \code{facet_by_factor = TRUE}, a single \pkg{ggplot2} object.
#'     Otherwise a named list of \pkg{ggplot2} objects, one per selected
#'     factor.  Each heatmap shows the signed per-factor contribution
#'     \eqn{B_{:,j}\,u_j(t)} to each trait at each time step, averaged
#'     across genotypes.  Present only when \code{B_arr}, \code{U_arr},
#'     and \code{plot_control_influence = TRUE}.}
#' }
#'
#' @seealso [run_dmd()], [compare_arrays()], [plot_trajectories()]
#' @export
reconstruct_dmd <- function(
    A_arr,
    x0                    = NULL,
    X_obs                 = NULL,
    B_arr                 = NULL,
    U_arr                 = NULL,
    n_steps               = NULL,
    return_single_step     = FALSE,
    plot_control_influence = TRUE,
    ctrl_factors           = NULL,
    facet_by_factor        = FALSE,
    operator_type          = "discrete",
    ts_days                = NULL,
    dt                     = NULL
) {
  #  Promote 2-D / vector inputs to 3-D 
  A_arr <- .to_3d_arr(A_arr)
  B_arr <- .to_3d_arr(B_arr)
  U_arr <- .to_3d_arr(U_arr)
  X_obs <- .to_3d_arr(X_obs)

  if (!is.null(B_arr) && is.null(U_arr))
    stop("reconstruct_dmd: B_arr supplied but U_arr is NULL -- both are required.")

  has_ctrl <- !is.null(B_arr) && !is.null(U_arr)

  p <- dim(A_arr)[1L]
  G <- dim(A_arr)[3L]

  # Broadcast a single-genotype (2-D promoted) U_arr or B_arr to all G genotypes.
  if (has_ctrl && dim(U_arr)[3L] == 1L && G > 1L) {
    U_arr <- array(rep(U_arr, G),
                   dim      = c(dim(U_arr)[1L], dim(U_arr)[2L], G),
                   dimnames = list(dimnames(U_arr)[[1L]], dimnames(U_arr)[[2L]], NULL))
  }
  if (has_ctrl && dim(B_arr)[3L] == 1L && G > 1L) {
    B_arr <- array(rep(B_arr, G),
                   dim      = c(dim(B_arr)[1L], dim(B_arr)[2L], G),
                   dimnames = list(dimnames(B_arr)[[1L]], dimnames(B_arr)[[2L]], NULL))
  }

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
  q       <- if (has_ctrl) dim(B_arr)[2L] else 0L
  N_u     <- if (has_ctrl) dim(U_arr)[2L] else 0L

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

  # Helper: compute the matrix exponential propagator and ZOH B contribution.
  # Returns list(eAdT = expm(A*dT), B_eff = A^{-1}*(eAdT - I)*B  or B*dT fallback)
  .ct_step <- function(A_g, B_g, dT) {
    eAdT <- expm::expm(A_g * dT)
    if (!is.null(B_g)) {
      # ZOH: B_eff = A_c^{-1} (expm(A_c*dt) - I) B_c
      IeAdT <- eAdT - diag(nrow(A_g))
      B_eff <- tryCatch(
        solve(A_g, IeAdT) %*% B_g,
        error = function(e) B_g * dT   # singular A fallback
      )
    } else {
      B_eff <- NULL
    }
    list(eAdT = eAdT, B_eff = B_eff)
  }

  #  Multi-step reconstruction 
  X_pred <- array(NA_real_,
                  dim      = c(p, n_steps + 1L, G),
                  dimnames = list(t_names, NULL, geno_names))

  for (gi in seq_len(G)) {
    A_g <- matrix(A_arr[, , gi], nrow = p, ncol = p)
    B_g <- if (has_ctrl) matrix(B_arr[, , gi], nrow = p, ncol = q) else NULL
    X_pred[, 1L, gi] <- x0_mat[, gi]
    for (t in seq_len(n_steps)) {
      if (is_continuous) {
        cs      <- .ct_step(A_g, B_g, dt_vec[t])
        x_next  <- cs$eAdT %*% X_pred[, t, gi]
        if (has_ctrl && t <= N_u)
          x_next <- x_next + cs$B_eff %*% as.numeric(U_arr[, t, gi])
      } else {
        x_next <- A_g %*% X_pred[, t, gi]
        if (has_ctrl && t <= N_u)
          x_next <- x_next + B_g %*% as.numeric(U_arr[, t, gi])
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
      B_g <- if (has_ctrl) matrix(B_arr[, , gi], nrow = p, ncol = q) else NULL
      for (t in seq_len(N_obs - 1L) + 1) {
        if (is_continuous) {
          cs      <- .ct_step(A_g, B_g, dt_vec_1step[t])
          x_next  <- cs$eAdT %*% as.numeric(X_obs[, t - 1, gi])
          if (has_ctrl && t <= N_u)
            x_next <- x_next + cs$B_eff %*% as.numeric(U_arr[, t - 1, gi])
        } else {
          x_next <- A_g %*% as.numeric(X_obs[, t - 1, gi])
          if (has_ctrl && t <= N_u)
            x_next <- x_next + B_g %*% as.numeric(U_arr[, t - 1, gi])
        }
        X_pred_1step[, t, gi] <- Re(as.vector(x_next))
      }
    }
  }

  #  Control influence heatmap (per-factor)
  ctrl_plot <- NULL
  if (has_ctrl && plot_control_influence) {
    N_plot <- min(N_u, n_steps)
    u_nms  <- dimnames(U_arr)[[1L]]
    if (is.null(u_nms)) u_nms <- paste0("U", seq_len(q))

    # Resolve ctrl_factors to a character vector of factor names.
    if (!is.null(ctrl_factors)) {
      if (is.numeric(ctrl_factors))
        ctrl_factors <- u_nms[as.integer(ctrl_factors)]
      bad <- setdiff(ctrl_factors, u_nms)
      if (length(bad))
        warning(sprintf(
          "reconstruct_dmd: ctrl_factors not found in U_arr dimension names: %s",
          paste(bad, collapse = ", ")
        ))
      ctrl_factors <- intersect(ctrl_factors, u_nms)
    } else {
      ctrl_factors <- u_nms
    }
    fac_idx <- match(ctrl_factors, u_nms)
    n_fac   <- length(fac_idx)

    if (n_fac > 0L) {
      # Bu_j[i, t] = (1/G) * sum_g  B_g[i, j] * U_arr[j, t, g]
      bu_rows <- lapply(seq_along(fac_idx), function(fi) {
        j   <- fac_idx[fi]
        mat <- matrix(0, nrow = p, ncol = N_plot)
        for (gi in seq_len(G)) {
          B_g <- matrix(B_arr[, , gi], nrow = p, ncol = q)
          b_j <- B_g[, j]
          for (t in seq_len(N_plot))
            mat[, t] <- mat[, t] + b_j * as.numeric(U_arr[j, t, gi])
        }
        mat <- mat / G
        data.frame(
          factor_name = ctrl_factors[fi],
          trait_fac   = rep(factor(t_names, levels = rev(t_names)), times = N_plot),
          time_step   = rep(seq_len(N_plot), each = p),
          bu_val      = as.vector(mat),
          stringsAsFactors = FALSE
        )
      })
      bu_df <- do.call(rbind, bu_rows)
      bu_df$factor_name <- factor(bu_df$factor_name, levels = ctrl_factors)

      show_y  <- (p <= 40L)
      abs_max <- max(abs(bu_df$bu_val[is.finite(bu_df$bu_val)]), na.rm = TRUE)
      if (!is.finite(abs_max) || abs_max == 0) abs_max <- 1

      .ctrl_hm <- function(df_sub, ttl) {
        ggplot2::ggplot(df_sub,
                        ggplot2::aes(x = time_step, y = trait_fac, fill = bu_val)) +
          ggplot2::geom_tile() +
          ggplot2::scale_fill_gradient2(
            low      = "#2166AC",
            mid      = "white",
            high     = "#D6604D",
            midpoint = 0,
            limits   = c(-abs_max, abs_max),
            name     = "B[,j]\u00d7u[j,t]"
          ) +
          ggplot2::labs(
            title    = ttl,
            subtitle = sprintf("Signed mean across %d genotype(s)", G),
            x        = "Time step",
            y        = "State trait"
          ) +
          ggplot2::theme_bw(base_size = 10) +
          ggplot2::theme(
            axis.text.y  = if (show_y) ggplot2::element_text(size = 7)
                           else        ggplot2::element_blank(),
            axis.ticks.y = if (show_y) ggplot2::element_line()
                           else        ggplot2::element_blank(),
            panel.grid   = ggplot2::element_blank()
          )
      }

      if (n_fac == 1L || facet_by_factor) {
        plt <- .ctrl_hm(bu_df, "Control influence  B[,j] u[j,t]  per trait")
        if (n_fac > 1L)
          plt <- plt + ggplot2::facet_wrap(~ factor_name)
        ctrl_plot <- plt
      } else {
        ctrl_plot <- stats::setNames(
          lapply(ctrl_factors, function(fn)
            .ctrl_hm(
              bu_df[bu_df$factor_name == fn, , drop = FALSE],
              sprintf("Control influence  [%s]  per trait", fn)
            )
          ),
          ctrl_factors
        )
      }
    }
  }

  # Return only non-NULL elements
  out <- list(
    X_pred_arr       = X_pred,
    X_pred_1step_arr = X_pred_1step,
    ctrl_plot        = ctrl_plot
  )
  Filter(Negate(is.null), out)
}
