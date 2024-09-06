#' Calculate A_r from phi and R
#'
#' This function computes A_r matrices for a set of lines given a set of phi and R matrices. Matrices can be given individually or collected into a tensor.

A_r_from_phi_and_R <- function (phi_all, R_all) {

  if (length(dim(phi_all)) > 2) {
    lines <- intersect(dimnames(phi_all)[[3]], dimnames(R_all)[[3]])
    n_traits <- dim(phi_all)[1]
    A_reconst_all <- array(NA, dim = c(n_traits, n_traits, length(lines)), dimnames = list(dimnames(phi_all)[[1]], dimnames(phi_all)[[1]], lines))

    for (line in lines) {
      print(line)

      phi <- phi_all[, , line]
      R <- R_all[, , line]

      A <- phi %*% R %*% pracma::pinv(phi)
      colnames(A) <- rownames(A) <- rownames(phi)

      A_reconst_all[, , line] <- A
    }
    return(A_reconst_all)
  } else {
    A <- phi_all %*% R_all %*% pracma::pinv(phi_all)
    colnames(A) <- rownames(A) <- rownames(phi_all)
    return(A)
  }
}

