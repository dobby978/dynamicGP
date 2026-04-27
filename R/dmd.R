#' Dynamic Mode Decomposition (DMD)
#'
#' Fits the linear model \eqn{\mathbf{x}(t+1) = A \cdot \mathbf{x}(t)} by
#' SVD-truncated regression or full pseudoinverse.
#'
#' Two decomposition methods are available:
#' \describe{
#'   \item{`"direct"`}{Full pseudoinverse.  Computes \eqn{A = X_2 X_1^+} via
#'     `pracma::pinv()` with no rank truncation.  The rank of the result equals
#'     the numerical rank of `X1`.}
#'   \item{`"standard"`}{Standard exact DMD.  Computes the reduced operator
#'     \eqn{\tilde{A} = \hat{U}^T X_2 \hat{V} \hat{\Sigma}^{-1}} and extracts
#'     eigenvectors \eqn{W} and eigenvalues \eqn{\lambda} via `eigen()`.
#'     Modes are lifted back to full space:
#'     \eqn{\Phi = X_2 \hat{V} \hat{\Sigma}^{-1} W / \lambda}.}
#'   \item{`"schur-dmd"`}{Same \eqn{\tilde{A}}, but decomposed via the real Schur
#'     form \eqn{\tilde{A} = Q T Q^T} using [Matrix::Schur()].  More stable
#'     when eigenvalues are nearly degenerate.  Falls back to `eigen()` with a
#'     warning if the **Matrix** package is unavailable.}
#' }
#'
#' @param X1 Numeric matrix of shape `p × m`: state snapshots at times
#'   `t = 1, ..., m`.
#' @param X2 Numeric matrix of shape `p × m`: state snapshots at times
#'   `t = 2, ..., m+1` (one step ahead of `X1`).
#' @param r Integer DMD rank.  `NULL` (default) triggers automatic selection
#'   based on `var_thresh`.  Ignored when `method = "direct"`.
#' @param var_thresh Cumulative-variance threshold for automatic rank selection
#'   (default `0.95`).  The smallest `r` such that the leading `r` singular
#'   values of `X1` explain at least `var_thresh` of total variance is used.
#'   Ignored when `method = "direct"`.
#' @param use_pinv Deprecated.  Use `method = "direct"` instead.
#' @param break_after Integer vector of 1-based column indices of `X1` after
#'   which a time-series gap occurs.  The snapshot pair `(X1[,k], X2[,k])` is
#'   excluded so that cross-gap transitions do not enter the regression.
#'   Default `integer(0)`.
#' @param method Character string: `"direct"` (full pseudoinverse),
#'   `"standard"` (default), or `"schur-dmd"` (Schur decomposition).
#' @param dt_days Timestep in days, used only to convert discrete eigenvalues
#'   \eqn{\lambda} to continuous-time frequencies and growth rates via
#'   \eqn{\omega = \log(\lambda) / \Delta t}.  Default `1` (eigenvalues
#'   reported in per-timestep units).
#'
#' @return A named list with components:
#' \describe{
#'   \item{`A_full`}{Full p \eqn{\times} p state-transition matrix.}
#'   \item{`A_til`}{Reduced r \eqn{\times} r operator (NULL when `use_pinv = TRUE`).}
#'   \item{`B_full`, `B_til`}{`NULL` (no control inputs for plain DMD).}
#'   \item{`lam`}{Complex vector of `r` DMD eigenvalues.}
#'   \item{`Phi`}{Complex `p × r` matrix of DMD modes.}
#'   \item{`omega_c`}{Continuous-time eigenvalues \eqn{\log(\lambda)/\Delta t}.}
#'   \item{`freq_cpd`}{Frequencies in cycles per day.}
#'   \item{`growth_day`}{Growth/decay rates per day.}
#'   \item{`amp`}{Mode amplitudes.}
#'   \item{`rho_A`}{Spectral radius of `A_full`.}
#'   \item{`resid_F`}{Relative Frobenius residual
#'     \eqn{\|X_2 - A X_1\|_F / \|X_2\|_F}.}
#'   \item{`r`, `eps`}{Chosen output rank (eps = r for DMD).}
#'   \item{`pv_om`, `pv_x2`}{Cumulative variance explained by the SVD of
#'     `X1` and `X2` respectively.}
#'   \item{`method`}{Character tag: `"direct"`, `"standard"`, or `"schur-dmd"`.}
#' }
#'
#' @references
#' Schmid, P.J. (2010). Dynamic mode decomposition of numerical and
#' experimental data. *Journal of Fluid Mechanics*, 656, 5–28.
#' \doi{10.1017/S0022112010001217}
#'
#' Tu, J.H. *et al.* (2014). On dynamic mode decomposition: theory and
#' applications. *Journal of Computational Dynamics*, 1(2), 391–421.
#' \doi{10.3934/jcd.2014.1.391}
#'
#' @seealso [reconstruct_dmd()]
#' @export
run_dmd <- function(X1, X2,
                    r           = NULL,
                    var_thresh  = 0.95,
                    use_pinv    = FALSE,
                    break_after = integer(0L),
                    method      = c("direct", "standard", "schur-dmd"),
                    dt_days     = 1) {
  method <- match.arg(method)
  if (use_pinv) {
    warning("run_dmd: `use_pinv` is deprecated; using method = \"direct\" instead.")
    method <- "direct"
  }

  # -- Gap-exclusion ----------------------------------------------------------
  if (length(break_after) > 0L) {
    ba_v <- as.integer(break_after)
    ba_v <- ba_v[ba_v >= 1L & ba_v <= ncol(X1)]
    if (length(ba_v) > 0L) {
      keep <- setdiff(seq_len(ncol(X1)), ba_v)
      X1   <- X1[, keep, drop = FALSE]
      X2   <- X2[, keep, drop = FALSE]
    }
  }

  p <- nrow(X1)
  m <- ncol(X1)
  if (!all(is.finite(X1))) stop("run_dmd: X1 contains non-finite values (NA/Inf/NaN).")
  if (!all(is.finite(X2))) stop("run_dmd: X2 contains non-finite values (NA/Inf/NaN).")

  sv_x1 <- svd(X1)
  pv_x1 <- cumsum(sv_x1$d^2) / sum(sv_x1$d^2)
  sv_x2 <- svd(X2, nu = p, nv = min(p, m))
  pv_x2 <- cumsum(sv_x2$d^2) / sum(sv_x2$d^2)

  if (method == "direct") {
    # -- Pseudoinverse path: A = X2 * pinv(X1) -----------------------------
    if (!requireNamespace("pracma", quietly = TRUE))
      stop("run_dmd: method = 'direct' requires the 'pracma' package.\n",
           "  Install it with: install.packages('pracma')")
    tol_pinv  <- max(dim(X1)) * .Machine$double.eps * sv_x1$d[1L]
    rank_pinv <- sum(sv_x1$d > tol_pinv)
    rk        <- seq_len(rank_pinv)
    pinv_X1   <- pracma::pinv(X1)
    A_full    <- X2 %*% pinv_X1
    eig_obj   <- eigen(A_full)
    lam       <- eig_obj$values
    Phi       <- eig_obj$vectors
    r         <- p
    eps_out   <- rank_pinv
    A_til     <- NULL;  U_hat <- NULL
  } else {
    if (is.null(r)) {
      r <- which(pv_x1 >= var_thresh)[1L]
      if (is.na(r)) r <- length(sv_x1$d)
    }
    r       <- min(as.integer(r), length(sv_x1$d), p)
    U_hat   <- sv_x1$u[, seq_len(r), drop = FALSE]   # p x r
    V_r     <- sv_x1$v[, seq_len(r), drop = FALSE]   # m x r
    sig_inv <- diag(1 / sv_x1$d[seq_len(r)], nrow = r)

    # Reduced operator  A_til = U_hat^T X2 V_r Sig_inv   (r x r)
    A_til <- t(U_hat) %*% X2 %*% V_r %*% sig_inv

    # Downgrade schur-dmd to standard if Matrix is unavailable
    if (method == "schur-dmd" && !requireNamespace("Matrix", quietly = TRUE)) {
      warning("run_dmd: 'Matrix' package not available; falling back to eigen() for schur-dmd.")
      method <- "standard"
    }

    if (method == "schur-dmd") {
      # -- Real Schur decomposition: A_til = Q T Q^T ---------------------
      if (!requireNamespace("pracma", quietly = TRUE))
        stop("run_dmd: method = 'schur-dmd' requires the 'pracma' package.\n",
             "  Install it with: install.packages('pracma')")
      sch_obj    <- Matrix::Schur(A_til)
      lam        <- sch_obj$EValues
      Q_sch      <- sch_obj$Q
      R_shr      <- sch_obj$T
      Phi        <- X2 %*% V_r %*% sig_inv %*% Q_sch    # p x r
      A_full     <- Re(Phi %*% R_shr %*% pracma::pinv(Phi))
      A_full_til <- U_hat %*% A_til %*% t(U_hat)
    } else {
      # -- Standard exact DMD (also the Matrix-unavailable fallback) -----
      eig_r <- eigen(A_til)
      lam   <- eig_r$values
      W     <- eig_r$vectors
      Phi   <- X2 %*% V_r %*% sig_inv %*% W
      for (j in seq_len(r)) Phi[, j] <- Phi[, j] / lam[j]
      A_full     <- Phi %*% W %*% Conj(t(Phi))
      A_full_til <- U_hat %*% A_til %*% t(U_hat)
    }

    eps_out <- r
  }

  # -- Fit diagnostics --------------------------------------------------------
  pred_X2    <- Re(A_full %*% X1)
  resid_F    <- norm(X2 - pred_X2, "F") / norm(X2, "F")
  rho_A      <- max(Mod(lam))
  omega_c    <- log(as.complex(lam)) / dt_days
  freq_cpd   <- Im(omega_c) / (2 * pi)
  growth_day <- Re(omega_c)
  amp        <- apply(Phi, 2L, function(v) sqrt(sum(Mod(v)^2)))

  if (!exists("A_full_til")) A_full_til <- NULL

  list(
    A_til = A_til, B_til = NULL,
    A_full = A_full, A_full_til = A_full_til, B_full = NULL,
    eig   = if (method == "direct") eig_obj else list(values = lam),
    lam = lam, Phi = Phi,
    omega_c = omega_c, freq_cpd = freq_cpd, growth_day = growth_day, amp = amp,
    U_hat = U_hat, Ut = NULL,
    V_r     = if (method == "direct") NULL else V_r,
    sig_inv = if (method == "direct") NULL else sig_inv,
    sv_om = sv_x1, sv_x2 = sv_x2,
    pv_om = pv_x1, pv_x2 = pv_x2,
    r = r, eps = eps_out,
    resid_F = resid_F, rho_A = rho_A,
    method  = method
  )
}
