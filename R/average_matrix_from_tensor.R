# Compute the element-wise mean matrix across the third dimension of a 3-D array.

#' Compute the Element-Wise Mean Matrix from a 3-D Array
#'
#' Reduces a 3-D array \eqn{[p \times q \times G]} to a single
#' \eqn{[p \times q]} matrix by taking the mean across the third
#' (genotype/sample) dimension at every \eqn{(i, j)} position.
#'
#' @param tensor A 3-D numeric array \eqn{[p \times q \times G]}.  The first
#'   two dimensions define the matrix structure; the third dimension
#'   (genotypes/samples) is averaged over.  `dimnames` on the first two
#'   dimensions are preserved in the output.
#' @param na.rm Logical.  If `TRUE` (default), `NA` values are excluded when
#'   computing each element mean.  Set to `FALSE` to propagate `NA` whenever
#'   any slice has a missing entry at that position.
#'
#' @return A numeric matrix of dimensions \eqn{[p \times q]} whose
#'   \eqn{(i,j)} entry is the mean of `tensor[i, j, ]`.  Row and column
#'   names are taken from `dimnames(tensor)[[1]]` and `dimnames(tensor)[[2]]`.
#'
#' @examples
#' set.seed(1)
#' A_arr <- array(rnorm(4 * 4 * 20), dim = c(4, 4, 20),
#'                dimnames = list(paste0("T", 1:4), paste0("T", 1:4),
#'                                paste0("G", 1:20)))
#' A_mean <- tensor_to_mean_matrix(A_arr)
#' dim(A_mean)   # 4 x 4
#'
#' @seealso [mean_matrix_to_tensor()]
#' @export
tensor_to_mean_matrix <- function(tensor, na.rm = TRUE) {

  if (length(dim(tensor)) != 3L) {
    stop("tensor_to_mean_matrix: 'tensor' must be a 3-D array.")
  }

  # apply() over the first two margins computes the mean across slices.
  avg_matrix <- apply(tensor, c(1L, 2L), mean, na.rm = na.rm)

  # Explicitly restore dimnames in case apply() drops them on degenerate arrays.
  rownames(avg_matrix) <- dimnames(tensor)[[1L]]
  colnames(avg_matrix) <- dimnames(tensor)[[2L]]

  avg_matrix
}