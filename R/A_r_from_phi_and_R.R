# Reconstruct a state-transition operator from DMD modes and eigenvalues.

#' Reconstruct a State-Transition Operator from DMD Modes and Eigenvalues
#'
#' Computes the state-transition matrix \eqn{A = \Phi R \Phi^+} from a matrix
#' of DMD modes \eqn{\Phi} and the corresponding eigenvalue (or recurrence)
#' matrix \eqn{R}, where \eqn{\Phi^+} is the Moore-Penrose pseudoinverse.
#' Accepts either a single pair of 2-D matrices or two 3-D tensors (one mode /
#' eigenvalue matrix per genotype along the third dimension).
#'
#' @section Formula:
#' \deqn{A = \Phi \, R \, \Phi^+}
#' In DMD, \eqn{\Phi} is the matrix of dynamic modes and \eqn{R} is the
#' diagonal matrix of eigenvalues, so this inverts the eigendecomposition to
#' recover the linear operator \eqn{A}.  The pseudoinverse is computed via
#' \code{\link[pracma]{pinv}}.
#'
#' @param phi_all Either a \eqn{[p \times p]} numeric matrix of DMD modes, or
#'   a 3-D array \eqn{[p \times p \times G]} with one mode matrix per
#'   genotype/sample in the third dimension.  Row names (and, for arrays,
#'   third-dimension names) should be set.
#' @param R_all Either a \eqn{[p \times p]} numeric eigenvalue matrix
#'   (typically diagonal), or a 3-D array \eqn{[p \times p \times G]} whose
#'   third-dimension names align with those of `phi_all`.  Only genotypes
#'   present in both tensors are processed.
#'
#' @return
#' * **2-D inputs**: a single \eqn{[p \times p]} matrix with row and column
#'   names from `rownames(phi_all)`.
#' * **3-D inputs**: a 3-D array \eqn{[p \times p \times G]} where \eqn{G}
#'   is the number of genotypes common to both tensors.  The first two
#'   dimension names are taken from `phi_all`; the third from the intersection
#'   of the two tensors' third-dimension names.  Genotypes that could not be
#'   processed retain their initial `NA` values.
#'
#' @examples
#' \dontrun{
#' # Single pair (2-D) ---------------------------------------------------
#' set.seed(1)
#' phi <- matrix(rnorm(16), 4, 4,
#'               dimnames = list(paste0("T", 1:4), paste0("T", 1:4)))
#' R   <- diag(c(0.9, 0.8, 0.7, 0.6))
#' A_rec <- operator_from_modes(phi, R)
#'
#' # Tensor (multi-genotype, 3-D) ----------------------------------------
#' n_geno  <- 10
#' phi_arr <- array(rnorm(4 * 4 * n_geno), dim = c(4, 4, n_geno),
#'                  dimnames = list(paste0("T", 1:4), paste0("T", 1:4),
#'                                  paste0("G", 1:n_geno)))
#' R_arr   <- array(0, dim = dim(phi_arr), dimnames = dimnames(phi_arr))
#' for (g in seq_len(n_geno)) diag(R_arr[, , g]) <- runif(4, 0.5, 1)
#' A_arr <- operator_from_modes(phi_arr, R_arr)
#' dim(A_arr)   # 4 x 4 x 10
#' }
#'
#' @seealso [reconstruct_dmd()]
#' @export
operator_from_modes <- function(phi_all, R_all) {

  if (!requireNamespace("pracma", quietly = TRUE)) {
    stop("operator_from_modes: the 'pracma' package is required.\n",
         "  Install it with: install.packages('pracma')")
  }

  if (length(dim(phi_all)) > 2L) {
    # ── Tensor path: reconstruct one A matrix per genotype ─────────────────
    lines    <- intersect(dimnames(phi_all)[[3L]], dimnames(R_all)[[3L]])
    n_traits <- dim(phi_all)[[1L]]
    trait_dn <- dimnames(phi_all)[[1L]]

    A_all <- array(
      NA_real_,
      dim      = c(n_traits, n_traits, length(lines)),
      dimnames = list(trait_dn, trait_dn, lines)
    )

    for (line in lines) {
      phi <- phi_all[, , line]
      R   <- R_all[, , line]
      A   <- phi %*% R %*% pracma::pinv(phi)
      rownames(A) <- colnames(A) <- rownames(phi)
      A_all[, , line] <- A
    }

    return(A_all)

  } else {
    # ── Single-matrix path ─────────────────────────────────────────────────
    A <- phi_all %*% R_all %*% pracma::pinv(phi_all)
    rownames(A) <- colnames(A) <- rownames(phi_all)
    return(A)
  }
}

