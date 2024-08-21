#' Find DMD componets
#'
#' This function takes as input a data frame in which the first two columns represent two different naming standards for the genotype - "bio_ID" is used and the second column is ignored.
#' The third column represents the time point, notated here as "DAS".
#' All subsequent columns are the traits of interest which are then organized into X matrices and used to find A directly or through the Schur-DMD algorithm. These two options are selected by setting method to "calculated.A" or "schur.DMD", respectively.
#' If "calculated.A" is used the function returns a list of tensors corresponding to the A and X matrices for each line. This corresponds to Algorithm 1 from Paper.
#' If "Schur.DMD" is used, the function returns a list of tensors corresponding to the components of the Schur-DMD algorithm, this corresponds to Algorithm 2 from Paper.
#' r indicates the number of singular values and vectors included in the truncated SVD step of Schur-DMD.
#' scale = "min.max" min-max scales the data for each trait and accession.

find_DMD_components <- function (data, r = 25, method = "calculated.A", scale = "min.max") {

  accessions <- unique(data$bio_ID)

  n_traits <- ncol(data) - 3
  n_timepoints <- length(unique(data$DAS))
  A_all <- array(NA, dim = c(n_traits, n_traits, length(accessions)), dimnames = list(colnames(data)[-c(1, 2, 3)], colnames(data)[-c(1, 2, 3)], accessions))
  X_all <- array(NA, dim = c(n_traits, n_timepoints, length(accessions)), dimnames = list(colnames(data)[-c(1, 2, 3)], paste0("t", unique(data$DAS)), accessions))

  if (method == "schur.DMD") {
    U_r_transpose_all <- array(NA, dim = c(r, n_traits, length(accessions)), dimnames = list(paste0("SV", 1:r), colnames(data)[-c(1, 2, 3)], accessions))

    # use this when not removing time points before gaps:
    V_r_all           <- array(NA, dim = c(n_timepoints - 1, r, length(accessions)), dimnames = list(paste0("t", 1:(n_timepoints - 1)), paste0("SV", 1:r), accessions))

    # use this when not removing time points before gaps:
    # V_r_all           <- array(NA, dim = c(n_timepoints - 1 - 4, r, length(accessions)), dimnames = list(paste0("t", 1:(n_timepoints - 1 - 4)), paste0("SV", 1:r), accessions))

    sigma_r_all       <- array(NA, dim = c(r, r, length(accessions)), dimnames = list(paste0("SV", 1:r), paste0("SV", 1:r), accessions))
    A_tilde_all       <- array(NA, dim = c(r, r, length(accessions)), dimnames = list(paste0("SV", 1:r), paste0("SV", 1:r), accessions))
    Q_all             <- array(NA, dim = c(r, r, length(accessions)), dimnames = list(paste0("SV", 1:r), paste0("SV", 1:r), accessions))
    R_all             <- array(NA, dim = c(r, r, length(accessions)), dimnames = list(paste0("SV", 1:r), paste0("SV", 1:r), accessions))
    phi_all           <- array(NA, dim = c(n_traits, r, length(accessions)), dimnames = list(colnames(data)[-c(1, 2, 3)], paste0("SV", 1:r), accessions))
  }

  if (scale == "min.max") {
    print("min.max scaling...")
    norm_minmax <- function(x) {
        (x - min(x)) / (max(x) - min(x))
    }
    data[, -c(1:3)] <- as.data.frame(lapply(data[, -c(1:3)] , norm_minmax))
  }

  for (accession in accessions) {

    X <- t(data[which(data$bio_ID == accession), -c(1, 2, 3)])

    X_1 <- X[, -(ncol(X))]
    X_2 <- X[, -1]

    ## # for removing the step before and step after the two day gap from X_1 and X_2, respectively.
    ## X_1 <- X_1[, -seq(5, 24, 5)]
    ## X_2 <- X_2[, -seq(5, 24, 5)]

    if (method == "schur.DMD") {
      # Thitsa et al. 2021 III. MAIN RESULTS A. Schur-Based DMD

      x <- svd(X_1)
      sigma <- matrix(0, ncol(x$u), ncol(x$u))
      diag(sigma) <- x$d

      U <- x$u
      V <- t(x$v)

      U_r_transpose <- Conj(t(U[, 1:r]))
      V_r <- t(V[1:r, ])
      sigma_r <- sigma[1:r, 1:r]
      diag(sigma_r) <- 1 / diag(sigma_r)
      diag(sigma_r)[-c(1:r)] <- 0

      A_tilde <- U_r_transpose %*% X_2 %*% V_r %*% sigma_r

      schur <- Matrix::Schur(A_tilde, vectors = TRUE)

      Q <- schur$Q
      R <- schur$T

      phi <- X_2 %*% V_r %*% sigma_r %*% Q
      A <- phi %*% R %*% pracma::pinv(phi)

      U_r_transpose_all[, , accession] <- U_r_transpose
      V_r_all[, , accession] <- V_r
      sigma_r_all[, , accession] <- sigma_r
      A_tilde_all[, , accession] <- A_tilde

      Q_all[, , accession] <- Q
      R_all[, , accession] <- R

      phi_all[, , accession] <- phi

    } else if (method == "calculated.A") {
      # discussion of different R packages for performing pseudoinversion: https://stackoverflow.com/questions/48097973/pseudoinverse-different-results-in-r-c-and-python
      # A_corpcor <- psi %*% phi %*% corpcor::pseudoinverse(psi, tol = (1e-10)*10)
      # A_MASS <- psi %*% phi %*% MASS::ginv(psi)

      X_1_pseudo <- pracma::pinv(X_1)
      A <- X_2 %*% X_1_pseudo
    }

    A_all[, , accession] <- A
    X_all[, , accession] <- X
  }

  if (method == "schur.DMD") {
    return (list(A_all = A_all, X_all = X_all, U_r_transpose_all = U_r_transpose_all, V_r_all = V_r_all, sigma_r_all = sigma_r_all, A_tilde_all = A_tilde_all, Q_all = Q_all, R_all = R_all, phi_all = phi_all))
  } else {
    return (list(A_all = A_all, X_all = X_all))
  }
}
