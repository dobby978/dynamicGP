#' Fit DMD Models Across Multiple Genotypes (or Samples)
#'
#' The DMD-only counterpart of [run_dmdc_all()].  Accepts 2-D (single
#' genotype) or 3-D (multiple genotypes) state arrays, normalises the data,
#' and fits a Dynamic Mode Decomposition (DMD) model for each genotype using
#' the chosen method.
#'
#' @details
#' **Input shapes:**
#' \itemize{
#'   \item `X_array`: `[p × N]` (single genotype, promoted internally to
#'     `[p × N × 1]`) or `[p × N × G]` (G genotypes).
#'   \item `ts_list`: a `POSIXct` vector (single genotype) or named list of
#'     `POSIXct` vectors (multiple genotypes).
#' }
#'
#' **Normalisation** is applied globally across *all* genotypes before
#' fitting so that operator magnitudes are directly comparable:
#' \itemize{
#'   \item `"minmax"` (default): map to the range \[0, 1\] using global min and
#'     range.  Sensitive to outliers; prefer `"zscore"` when cross-genotype
#'     magnitude comparability is less important than robustness.
#'   \item `"zscore"`: subtract global mean, divide by global SD.
#'     Zero-SD traits are set to 0.
#'   \item `"none"`: no normalisation.
#' }
#'
#' **Methods** (selected via the `method` argument):
#' \describe{
#'   \item{`"standard"` (default)}{SVD-truncated exact DMD.
#'     Fits \eqn{x(t+1) = Ax(t)}.}
#'   \item{`"schur-dmd"`}{Schur-decomposition DMD (requires \pkg{Matrix}).
#'     Falls back to `"standard"` if \pkg{Matrix} is unavailable.}
#'   \item{`"direct"`}{Full pseudoinverse solution (requires \pkg{pracma}).}
#' }
#'
#' @param X_array Numeric array `[p × N]` or `[p × N × G]` of state
#'   (trait) snapshots.  Traits are rows, time is columns, genotypes are the
#'   third dimension.  `dimnames(X_array)[[1]]` (trait names) and
#'   `dimnames(X_array)[[3]]` (genotype names) are used when present.
#' @param ts_list `POSIXct` vector or named list of vectors giving the
#'   timestamp for each column of `X_array`.  `NULL` falls back to integer
#'   indices.
#' @param trait_names Optional character vector of trait display names.
#'   Falls back to `dimnames(X_array)[[1]]`.
#' @param var_thresh Cumulative-variance threshold for automatic SVD rank
#'   selection (default `0.95`).
#' @param r Integer output rank override.  `NULL` (default) → auto.
#'   Ignored when `method = "direct"`.
#' @param return_intermediates Logical.  If `TRUE` (default) also return
#'   reduced-space intermediate matrices (Phi, A_til, U_hat, V_r, sig_inv)
#'   as 3-D arrays in `$inter`.
#' @param break_after Gap specification.  Can be:
#'   \describe{
#'     \item{`NULL`}{No exclusions (default).}
#'     \item{integer vector}{Same indices applied to all genotypes.}
#'     \item{named list}{Per-genotype integer vectors,
#'       e.g. `list(D04b = c(480L, 960L), I16b = c(480L))`.}
#'   }
#'   1-based indices in the per-genotype filtered (complete-cases) time
#'   series.  Snapshot pair `(X[,i], X[,i+1])` is excluded.
#' @param norm_method Normalisation for state `X`: `"minmax"` (default),
#'   `"zscore"`, or `"none"`.
#' @param norm_stats Pre-computed normalisation statistics (list with fields
#'   `state_mean`, `state_sd`, `state_min`, `state_range`) injected by an
#'   external wrapper when data should be scaled relative to a pooled
#'   reference.  `NULL` (default) computes statistics from `X_array`.
#' @param method Fitting method: `"standard"` (default), `"schur-dmd"`, or
#'   `"direct"`.
#' @param irregular_handling Action when irregular time spacing is detected:
#'   `"error"` (default, stop with an informative message) or `"warn"`
#'   (emit a warning and continue with the chosen `method`).
#'
#' @return A named list with five elements:
#' \describe{
#'   \item{`dat`}{Named list of per-genotype result lists, each containing:
#'     timestamps (`ts`, `ts_days`), normalisation parameters
#'     (`norm_method`, `state_center`, `state_scale`), spectral results
#'     (`r`, `rho_A`, `resid_F`, `lam`, `omega_c`, `freq_cpd`,
#'     `growth_day`, `amp`, `Phi`), and metadata (`method`, `break_after`).}
#'   \item{`A_arr`}{3-D array `[p × p × G]` of state-transition matrices.}
#'   \item{`X_arr`}{3-D array `[p × N_max × G]` of normalised state data.}
#'   \item{`inter`}{List of reduced-space 3-D arrays (Phi, A_til, U_hat,
#'     V_r, sig_inv), or `NULL` when `return_intermediates = FALSE`.}
#' }
#'
#' @seealso [run_dmd()], [reconstruct_dmd()], [plot_trajectories()]
#' @export
run_dmd_all <- function(
    X_array              = NULL,
    ts_list              = NULL,
    trait_names          = NULL,
    var_thresh           = 0.95,
    r                    = NULL,
    return_intermediates = TRUE,
    break_after          = NULL,
    norm_method          = "minmax",
    norm_stats           = NULL,
    method               = "standard",
    irregular_handling   = "error"
) {
  # Promote 2-D matrices to 3-D arrays (single-genotype convenience form)
  if (is.matrix(X_array) || length(dim(X_array)) == 2L) {
    gname   <- if (!is.null(names(ts_list))) names(ts_list)[[1L]] else "G1"
    X_array <- array(X_array, dim = c(dim(X_array), 1L),
                     dimnames = list(dimnames(X_array)[[1L]], NULL, gname))
  }
  if (!is.null(ts_list) && !is.list(ts_list))
    ts_list <- list(ts_list)

  stopifnot(is.array(X_array), length(dim(X_array)) == 3L)

  norm_method        <- match.arg(norm_method, c("zscore", "minmax", "none"))
  method             <- match.arg(method, c("standard", "schur-dmd", "direct"))
  irregular_handling <- match.arg(irregular_handling, c("error", "warn"))

  # ── Global normalisation statistics ────────────────────────────────────────
  if (!is.null(norm_stats)) {
    if (norm_method == "zscore") {
      state_global_mean <- norm_stats$state_mean
      state_global_sd   <- norm_stats$state_sd
    } else if (norm_method == "minmax") {
      state_global_min   <- norm_stats$state_min
      state_global_range <- norm_stats$state_range
    }
  } else {
    if (norm_method == "zscore") {
      state_global_mean <- apply(X_array, 1L, function(v) mean(v[is.finite(v)], na.rm = TRUE))
      state_global_sd   <- apply(X_array, 1L, function(v) sd(v[is.finite(v)],   na.rm = TRUE))
      state_global_sd[state_global_sd == 0 | !is.finite(state_global_sd)] <- NA
    } else if (norm_method == "minmax") {
      state_global_min   <- apply(X_array, 1L, function(v) min(v[is.finite(v)], na.rm = TRUE))
      state_global_max   <- apply(X_array, 1L, function(v) max(v[is.finite(v)], na.rm = TRUE))
      state_global_range <- state_global_max - state_global_min
      state_global_range[state_global_range == 0 | !is.finite(state_global_range)] <- NA
    }
  }

  G          <- dim(X_array)[3L]
  geno_names <- if (!is.null(dimnames(X_array)[[3L]])) dimnames(X_array)[[3L]]
                else if (!is.null(names(ts_list)))      names(ts_list)
                else paste0("G", seq_len(G))

  t_names <- if (!is.null(trait_names))                  trait_names
             else if (!is.null(dimnames(X_array)[[1L]])) dimnames(X_array)[[1L]]
             else paste0("T", seq_len(dim(X_array)[1L]))

  input_data <- setNames(lapply(seq_len(G), function(gi) {
    state_mat_raw           <- t(X_array[, , gi])
    colnames(state_mat_raw) <- t_names
    ts_raw <- if (!is.null(ts_list)) {
      ts_list[[gi]]
    } else {
      col_nms <- dimnames(X_array)[[2L]]
      if (!is.null(col_nms)) {
        num_nms <- suppressWarnings(as.numeric(col_nms))
        if (all(is.finite(num_nms))) num_nms else seq_len(dim(X_array)[2L])
      } else {
        seq_len(dim(X_array)[2L])
      }
    }
    list(state_mat = state_mat_raw, ts = ts_raw,
         n_total = dim(X_array)[2L], g = geno_names[gi])
  }), geno_names)

  cat("\n=============================================================\n")
  cat("Running DMD per genotype\n")
  cat("=============================================================\n")

  dmd_dat     <- list()
  X_list      <- list()
  inter_stash <- list()

  p_dim <- length(t_names)
  G_dim <- length(geno_names)

  A_arr <- array(NA_real_, dim = c(p_dim, p_dim, G_dim),
                 dimnames = list(t_names, t_names, geno_names))

  for (item in input_data) {
    g         <- item$g
    state_mat <- item$state_mat
    ts_all    <- item$ts
    n_total   <- item$n_total
    cat(sprintf("\n-- %s --------------------------------------------\n", g))

    state_mat[!is.finite(state_mat)] <- NA
    valid    <- complete.cases(state_mat)
    state_ok <- state_mat[valid, , drop = FALSE]
    ts_ok    <- ts_all[valid]
    N        <- nrow(state_ok)

    cat(sprintf("  Complete timesteps : %d / %d  (%.1f%%)\n",
                N, n_total, 100 * N / n_total))

    # ── Detect irregular time spacing ──────────────────────────────────────
    ts_days_g <- if (inherits(ts_ok, "POSIXct"))
                   as.numeric(ts_ok, units = "days") / 86400
                 else
                   as.numeric(ts_ok)
    dts_g    <- diff(ts_days_g)
    cv_g     <- if (length(dts_g) >= 2L && mean(dts_g) > 0)
                  sd(dts_g) / mean(dts_g) else 0
    is_irreg <- cv_g > 0.01

    if (is_irreg) {
      msg <- sprintf(
        paste0("run_dmd_all: genotype '%s' has irregular time spacing ",
               "(CV=%.1f%%, mean dt=%.3f days).  Use run_dmdc_all() with ",
               "method='cOptDMDc' for irregular series, or set ",
               "irregular_handling='warn' to continue anyway."),
        g, 100 * cv_g, mean(dts_g)
      )
      if (irregular_handling == "error") stop(msg)
      else warning(msg)
    }

    cat(sprintf("  Time spacing       : %s (CV=%.2f%%), method = %s\n",
                if (is_irreg) "irregular" else "regular", 100 * cv_g, method))

    dt_mean_g <- if (length(dts_g) > 0L) mean(dts_g) else 1

    # ── Normalise state X ───────────────────────────────────────────────────
    if (norm_method == "minmax") {
      cat("  Normalisation      : min-max scaling X to [0, 1]\n")
      state_sc <- sweep(sweep(state_ok, 2L, state_global_min, "-"),
                        2L, state_global_range, "/")
      state_sc[!is.finite(state_sc)] <- 0
      const_state <- names(which(is.na(state_global_range)))
      if (length(const_state) > 0)
        cat(sprintf("  WARNING: zero-range trait(s) set to 0: %s\n",
                    paste(const_state, collapse = ", ")))
      state_center <- state_global_min
      state_scale  <- state_global_range
    } else if (norm_method == "none") {
      cat("  Normalisation      : none\n")
      state_sc     <- state_ok
      state_center <- setNames(rep(0, ncol(state_ok)), colnames(state_ok))
      state_scale  <- setNames(rep(1, ncol(state_ok)), colnames(state_ok))
    } else {
      cat("  Normalisation      : z-score scaling X\n")
      state_sc <- sweep(sweep(state_ok, 2L, state_global_mean, "-"),
                        2L, state_global_sd, "/")
      state_sc[!is.finite(state_sc)] <- 0
      const_state <- names(which(is.na(state_global_sd)))
      if (length(const_state) > 0)
        cat(sprintf("  WARNING: zero-sd trait(s) set to 0: %s\n",
                    paste(const_state, collapse = ", ")))
      state_center <- state_global_mean
      state_scale  <- state_global_sd
    }

    X  <- t(state_sc)
    X1 <- X[, -N]
    X2 <- X[, -1L]

    # ── Drop snapshot pairs spanning time-series gaps ──────────────────────
    ba <- if (is.list(break_after)) break_after[[g]] else break_after
    if (!is.null(ba) && length(ba) > 0L) {
      ba_valid <- ba[ba >= 1L & ba <= ncol(X1)]
      if (length(ba_valid) < length(ba))
        warning(sprintf("%s: %d break_after index(es) out of range [1, %d] - ignored.",
                        g, length(ba) - length(ba_valid), ncol(X1)))
      if (length(ba_valid) > 0L) {
        keep <- setdiff(seq_len(ncol(X1)), ba_valid)
        X1   <- X1[, keep, drop = FALSE]
        X2   <- X2[, keep, drop = FALSE]
        cat(sprintf("  Gap exclusions     : %d pair(s) removed at index(es): %s\n",
                    length(ba_valid), paste(sort(ba_valid), collapse = ", ")))
      }
    }

    ba_g <- if (!is.null(ba) && length(ba) > 0L)
              sort(ba[ba >= 1L & ba <= (N - 1L)]) else integer(0L)

    res <- run_dmd(X1, X2, r = r, var_thresh = var_thresh,
                   method = method, dt_days = dt_mean_g)

    rownames(res$A_full) <- t_names
    colnames(res$A_full) <- t_names

    if (method != "direct") {
      cat(sprintf("  Output rank r      : %d / %d max\n", res$r, p_dim))
    } else {
      cat("  Output rank r      : full (direct pseudoinverse)\n")
    }
    cat(sprintf("  Spectral radius rho: %.6f\n", res$rho_A))
    cat(sprintf("  Relative residual  : %.4e  (||X2 - AX1||_F / ||X2||_F)\n",
                res$resid_F))

    A_arr[,, g] <- res$A_full
    X_list[[g]] <- X
    if (return_intermediates)
      inter_stash[[g]] <- list(A_til   = res$A_til,
                               U_hat   = res$U_hat,
                               V_r     = res$V_r,
                               sig_inv = res$sig_inv)

    dmd_dat[[g]] <- list(
      ts           = ts_ok,
      ts_days      = ts_days_g,
      break_after  = ba_g,
      norm_method  = norm_method,
      state_center = state_center,
      state_scale  = state_scale,
      method       = method,
      r            = res$r,
      eps          = res$eps,
      rho_A        = res$rho_A,
      resid_F      = res$resid_F,
      pv_om        = res$pv_om,
      pv_x2        = res$pv_x2,
      lam          = res$lam,
      omega_c      = res$omega_c,
      freq_cpd     = res$freq_cpd,
      growth_day   = res$growth_day,
      amp          = res$amp,
      Phi          = res$Phi
    )
  }

  # ── Assemble state array [variable x time x genotype] ─────────────────────
  N_vec <- sapply(geno_names, function(g) length(dmd_dat[[g]]$ts))
  N_max <- max(N_vec)

  X_arr <- array(NA_real_, dim = c(p_dim, N_max, G_dim),
                 dimnames = list(t_names, NULL, geno_names))
  for (g in geno_names)
    X_arr[, seq_len(N_vec[[g]]), g] <- X_list[[g]]
  rm(X_list)

  # ── Optionally assemble intermediate reduced-space 3-D arrays ─────────────
  inter <- NULL
  if (return_intermediates) {

    r_vec <- sapply(geno_names, function(g) dmd_dat[[g]]$r)
    r_max <- max(r_vec)

    Phi_arr <- array(NA_complex_, dim = c(p_dim, r_max, G_dim),
                     dimnames = list(t_names, NULL, geno_names))
    for (g in geno_names)
      Phi_arr[, seq_len(r_vec[[g]]), g] <- dmd_dat[[g]]$Phi

    if (method %in% c("standard", "schur-dmd")) {
      A_til_arr <- array(NA_real_, dim = c(r_max, r_max, G_dim),
                         dimnames = list(NULL, NULL, geno_names))
      U_hat_arr <- array(NA_real_, dim = c(p_dim, r_max, G_dim),
                         dimnames = list(t_names, NULL, geno_names))
      for (g in geno_names) {
        r_g <- r_vec[[g]]
        if (!is.null(inter_stash[[g]]$A_til))
          A_til_arr[seq_len(r_g), seq_len(r_g), g] <- inter_stash[[g]]$A_til
        if (!is.null(inter_stash[[g]]$U_hat))
          U_hat_arr[, seq_len(r_g), g] <- inter_stash[[g]]$U_hat
      }
    } else {
      A_til_arr <- NULL
      U_hat_arr <- NULL
    }

    if (method %in% c("standard", "schur-dmd")) {
      vr_ncol_vec <- sapply(geno_names, function(g) {
        vr <- inter_stash[[g]]$V_r; if (is.null(vr)) NA_integer_ else ncol(vr)
      })
      si_nrow_vec <- sapply(geno_names, function(g) {
        si <- inter_stash[[g]]$sig_inv; if (is.null(si)) NA_integer_ else nrow(si)
      })
      vr_max <- if (all(is.na(vr_ncol_vec))) 0L else max(vr_ncol_vec, na.rm = TRUE)
      si_max <- if (all(is.na(si_nrow_vec))) 0L else max(si_nrow_vec, na.rm = TRUE)
      vr_nrow_max <- max(sapply(geno_names, function(g) {
        vr <- inter_stash[[g]]$V_r; if (is.null(vr)) 0L else nrow(vr)
      }))

      V_r_arr <- if (vr_max > 0L) {
        arr <- array(NA_real_, dim = c(vr_nrow_max, vr_max, G_dim),
                     dimnames = list(NULL, NULL, geno_names))
        for (g in geno_names) {
          vr <- inter_stash[[g]]$V_r
          if (!is.null(vr))
            arr[seq_len(nrow(vr)), seq_len(ncol(vr)), g] <- vr
        }
        arr
      } else NULL

      sig_inv_arr <- if (si_max > 0L) {
        arr <- array(NA_real_, dim = c(si_max, si_max, G_dim),
                     dimnames = list(NULL, NULL, geno_names))
        for (g in geno_names) {
          si <- inter_stash[[g]]$sig_inv
          if (!is.null(si)) { d <- nrow(si); arr[seq_len(d), seq_len(d), g] <- si }
        }
        arr
      } else NULL
    } else {
      V_r_arr     <- NULL
      sig_inv_arr <- NULL
    }

    inter <- list(
      Phi     = Phi_arr,
      A_til   = A_til_arr,
      U_hat   = U_hat_arr,
      V_r     = V_r_arr,
      sig_inv = sig_inv_arr
    )
    rm(inter_stash)
  }

  list(
    dat   = dmd_dat,
    A_arr = A_arr,
    X_arr = X_arr,
    inter = inter
  )
}
