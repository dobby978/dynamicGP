# Internal utilities shared across dynamicGP functions.
# None of these are exported.

# Package-level importFrom declarations (avoids R CMD CHECK NOTEs for
# frequently-used stats/base functions that don't receive :: qualification).
#' @importFrom stats complete.cases sd setNames cor median cor.test
#' @importFrom utils read.table write.csv write.table head capture.output
NULL

# Promote a vector, matrix, or 3-D array to a canonical [d1 × d2 × G] array.
# • NULL     -> NULL (pass-through)
# • vector   -> [p × 1 × 1]   (single trait, single time-point, single genotype)
# • matrix   -> [p × N × 1]   (single genotype)
# • 3-D array -> unchanged
.to_3d_arr <- function(arr, g_label = "G1") {
  if (is.null(arr)) return(NULL)
  if (is.numeric(arr) && is.null(dim(arr))) {
    return(array(arr, dim = c(length(arr), 1L, 1L),
                 dimnames = list(names(arr), NULL, g_label)))
  }
  if (is.matrix(arr) || (is.array(arr) && length(dim(arr)) == 2L)) {
    dn <- dimnames(arr)
    return(array(arr,
                 dim      = c(dim(arr)[1L], dim(arr)[2L], 1L),
                 dimnames = list(dn[[1L]], dn[[2L]], g_label)))
  }
  arr   # already 3-D (or higher)
}

# Parse timestamp strings to POSIXct using common date-time formats.
parse_ts <- function(x) {
  lubridate::parse_date_time(
    x,
    orders = c("dmy HM", "dmy HMS", "ymd HM", "ymd HMS"),
    tz     = "UTC"
  )
}

# Diverging heatmap theme used for A and B operator plots.
theme_heatmap <- function(base_size = 10) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 0, size = 8),
      axis.text.y      = ggplot2::element_text(size = 8),
      panel.grid       = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey92")
    )
}

#' Enforce Spectral Stability via Eigenvalue Clipping
#'
#' Clips the eigenvalues of one or more state-transition operator matrices
#' to ensure spectral stability (i.e., all eigenvalues have magnitude strictly
#' less than 1 for discrete-time systems). This is useful for enforcing
#' physically-constrained dynamics when DMD operators are estimated
#' without explicit stability constraints.
#'
#' For each eigenvalue \\eqn{\\lambda_i}, if \\eqn{|\\lambda_i| > 1 - \\epsilon},
#' it is rescaled to \\eqn{(1 - \\epsilon) \\cdot \\operatorname{sign}(\\lambda_i)}
#' (preserving argument/sign). Eigenvalues with \\eqn{|\\lambda_i| \\le 1 - \\epsilon}
#' are retained unchanged.
#'
#' @param A_all Numeric matrix `[p × p]` (single genotype) or 3-D array
#'   `[p × p × G]` (multiple genotypes) of state-transition operator matrices.
#' @param epsilon Numeric scalar or vector of length `G`.  Stability margin:
#'   clipping threshold is \\eqn{1 - \\epsilon}. Default `0.00001` (very tight
#'   margin from unit circle).  If a vector, element `i` applies to genotype `i`.
#'
#' @return A named list:
#' \\describe{
#'   \\item{`A_stable_all`}{3-D array `[p × p × G]` (or `[p × p × 1]` for
#'     matrix input) of stabilised state-transition operators, with
#'     all eigenvalues guaranteed to satisfy \\eqn{|\\lambda_i| \\le 1 - \\epsilon}.}
#'   \\item{`eigenvalues`}{Data frame with columns:
#'     \\itemize{
#'       \\item `Lambda`: original eigenvalue magnitude.
#'       \\item `evs_clipped`: clipped eigenvalue magnitude.
#'       \\item Column 3: genotype identifier.
#'     }
#'     Useful for diagnosing which eigenvalues were modified.}
#' }
#'
#' @seealso [run_dmd()]
#' @export
spectral_clipper <- function(A_all, epsilon = 0.00001) {

  # Promote 2-D matrix to 3-D array for uniform processing
  if (length(dim(A_all)) == 2) {
    A_all_3D <- array(A_all,
                      dim = c(nrow(A_all), ncol(A_all), 1),
                      dimnames = list(rownames(A_all), colnames(A_all), "Line_1"))
    A_all <- A_all_3D
  }

  # Extract dimensions and genotype labels
  n_genotypes <- dim(A_all)[3]
  genotypes <- dimnames(A_all)[[3]]

  # Initialise output arrays and data frame
  ev_df <- data.frame()
  A_stable_all <- A_all
  A_stable_all[] <- NA  # Clear to ensure we only write computed values

  # Process each genotype independently
  for (genotype in seq_len(n_genotypes)) {
    # Extract A matrix for this genotype
    A <- A_all[, , genotype]

    # Compute eigendecomposition A = V Λ V⁻¹
    ED <- eigen(A)
    Lambda <- ED$values
    V <- ED$vectors

    # Initialise clipped eigenvalues
    evs_clipped <- Lambda
    evs_clipped[] <- NA

    # Determine stability threshold for this genotype
    if (length(epsilon) == 1) {
      one_minus_epsilon <- 1 - epsilon
    } else {
      one_minus_epsilon <- 1 - epsilon[genotype]
    }

    # Clip each eigenvalue
    for (i in seq_len(length(Lambda))) {
      ev <- Lambda[i]
      # If eigenvalue is already stable, keep it unchanged
      if (abs(ev) <= one_minus_epsilon) {
        evs_clipped[i] <- ev
      } else {
        # Otherwise, rescale to stability threshold, preserving sign/argument
        evs_clipped[i] <- one_minus_epsilon * (ev / abs(ev))
      }
    }

    # Reconstruct A from clipped eigenvalues: A_stable = V Λ_clipped V⁻¹
    Lambda_clipped <- diag(evs_clipped)
    A_stable <- V %*% Lambda_clipped %*% solve(V)
    dimnames(A_stable) <- dimnames(A)
    A_stable_all[, , genotype] <- A_stable

    # Log original vs clipped eigenvalues for diagnostics
    evs <- cbind.data.frame(Lambda, evs_clipped, genotypes[genotype])
    ev_df <- rbind.data.frame(ev_df, evs)
  }

  # Return stabilised operators and eigenvalue comparison
  list(A_stable_all = A_stable_all, eigenvalues = ev_df)
}
