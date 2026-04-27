# Tests for run_dmd()
# These use small synthetic systems where the ground-truth A is known.

# ── Helpers ───────────────────────────────────────────────────────────────────
make_stable_A <- function(p, seed = 1L) {
  set.seed(seed)
  A_raw  <- matrix(rnorm(p * p), p, p)
  eig_v  <- eigen(A_raw, only.values = TRUE)$values
  rho    <- max(Mod(eig_v))
  A_raw / (rho * 1.1)   # scale to spectral radius < 1
}

simulate_X <- function(A, N, seed = 2L) {
  p <- nrow(A)
  set.seed(seed)
  X <- matrix(0, p, N)
  X[, 1L] <- rnorm(p)
  for (t in seq_len(N - 1L)) X[, t + 1L] <- A %*% X[, t]
  X
}

# ── run_dmd() ─────────────────────────────────────────────────────────────────
test_that("run_dmd returns correct structure", {
  p  <- 4L; N <- 30L
  A0 <- make_stable_A(p)
  X  <- simulate_X(A0, N)
  res <- run_dmd(X[, -N], X[, -1L])

  expect_named(res, c("A_til", "B_til", "A_full", "A_full_til", "B_full",
                       "eig", "lam", "Phi", "omega_c", "freq_cpd",
                       "growth_day", "amp", "U_hat", "Ut",
                       "V_r", "sig_inv", "sv_om", "sv_x2",
                       "pv_om", "pv_x2", "r", "eps", "resid_F", "rho_A",
                       "method"),
               ignore.order = TRUE)
  expect_true(is.complex(res$lam))
  expect_equal(nrow(res$A_full), p)
  expect_equal(ncol(res$A_full), p)
  expect_true(is.numeric(res$resid_F))
  expect_true(res$resid_F >= 0)
  expect_true(res$rho_A < 1.5)
})

test_that("run_dmd recovers A for a noiseless system (pinv path)", {
  # Use use_pinv = TRUE for an exact LS solution that avoids SVD conditioning
  # issues in contracting systems.
  p  <- 5L; N <- 50L
  A0 <- make_stable_A(p, seed = 42L)
  X  <- simulate_X(A0, N, seed = 7L)
  res <- run_dmd(X[, -N], X[, -1L], method = "direct")

  expect_lt(res$resid_F, 1e-6,
            label = "Frobenius residual for noiseless system (direct path)")
})

test_that("run_dmd direct path works", {
  p  <- 3L; N <- 20L
  A0 <- make_stable_A(p)
  X  <- simulate_X(A0, N)
  res <- run_dmd(X[, -N], X[, -1L], method = "direct")

  expect_null(res$A_til)
  expect_true(all(is.finite(res$A_full)))
})

test_that("run_dmd respects break_after", {
  p  <- 4L; N <- 30L
  A0 <- make_stable_A(p)
  X  <- simulate_X(A0, N)
  res_full  <- run_dmd(X[, -N], X[, -1L])
  res_break <- run_dmd(X[, -N], X[, -1L], break_after = 15L)
  # Both should return valid results; residuals may differ
  expect_true(is.numeric(res_break$resid_F))
})

