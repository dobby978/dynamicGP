# Broadcast a matrix to a 3-D tensor with identical slices.

#' Broadcast a Matrix to a 3-D Tensor with Identical Slices
#'
#' Creates a 3-D array \eqn{[p \times q \times n]} where every slice along
#' the third dimension is a copy of `avg_matrix`.  This is the inverse
#' operation to [tensor_to_mean_matrix()]: given a mean or baseline matrix,
#' it produces a constant tensor suitable for element-wise comparison or
#' combination with a per-genotype tensor.
#'
#' @param avg_matrix A numeric matrix \eqn{[p \times q]} to broadcast.  Row
#'   and column names, if present, are propagated to the first two dimensions
#'   of the output tensor.
#' @param n Positive integer.  Number of slices (genotypes/samples) in the
#'   output tensor.
#' @param line_names Character vector of length `n` providing names for the
#'   third dimension (e.g. genotype IDs).  `NULL` (default) leaves the third
#'   dimension unnamed.
#'
#' @return A 3-D numeric array of dimensions \eqn{[p \times q \times n]}
#'   where every third-dimension slice equals `avg_matrix`.  `dimnames` on
#'   the first two dimensions are taken from `avg_matrix`; the third dimension
#'   is named by `line_names`.
#'
#' @examples
#' set.seed(1)
#' A_mean   <- matrix(rnorm(16), 4, 4,
#'                    dimnames = list(paste0("T", 1:4), paste0("T", 1:4)))
#' A_tensor <- mean_matrix_to_tensor(A_mean, n = 20,
#'                                   line_names = paste0("G", 1:20))
#' dim(A_tensor)              # 4 x 4 x 20
#' all(A_tensor[, , 1] == A_mean)   # TRUE
#'
#' @seealso [tensor_to_mean_matrix()]
#' @export
mean_matrix_to_tensor <- function(avg_matrix, n, line_names = NULL) {

  if (!is.matrix(avg_matrix)) {
    stop("mean_matrix_to_tensor: 'avg_matrix' must be a matrix.")
  }
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 1L) {
    stop("mean_matrix_to_tensor: 'n' must be a single positive integer.")
  }
  if (!is.null(line_names) && length(line_names) != n) {
    stop("mean_matrix_to_tensor: 'line_names' must have length n (", n, ").")
  }

  # Recycle the column-major data vector of avg_matrix n times to fill the
  # third dimension without an explicit loop.
  tensor <- array(
    rep(avg_matrix, times = n),
    dim      = c(nrow(avg_matrix), ncol(avg_matrix), n),
    dimnames = list(rownames(avg_matrix), colnames(avg_matrix), line_names)
  )

  tensor
}