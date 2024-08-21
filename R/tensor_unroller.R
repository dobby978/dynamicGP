#' "Unroll" tensors
#'
#' This function reshapes the tensor objects output by find_DMD_components() into matrices with the genotypes as rows and the individual elements of the component matrices as columns, which can then be used for genomic prediction and heritability analysis.

tensor_unroller <- function (tensor_all) {
  unrolled <- as.data.frame(dimnames(tensor_all)[3])
  cols <- colnames(tensor_all[, , 1])
  rows <- rownames(tensor_all[, , 1])

  count <- 1

  for (row in rows) {
    for (col in cols) {
      unrolled <- cbind.data.frame(unrolled, tensor_all[row, col, ])
      colnames(unrolled)[ncol(unrolled)] <- paste0(row, "_-_-_", col)

      print(paste0(count, " of ", dim(tensor_all)[1] * dim(tensor_all)[2]))
      count <- count + 1
    }
  }
  rownames(unrolled) <- unrolled[, 1]
  unrolled <- unrolled[, -1]
  return (unrolled)
}
