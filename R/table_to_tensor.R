#' Transform table into tensor
#'
#' This function takes a table with rows as genotypes and columns as individual matrix entries with column names corresponding to the row and column name in the component matrix separated by "_._._" and arranges it into a tensor.

table_to_tensor <- function (table, component = "phi", r = 2, n_timepoints = NULL) {
  lines <- rownames(table)

  if (component == "A") {
    n_traits <- sqrt(ncol(table))
    A_all <- array(unlist(t(table)), dim = c(n_traits, n_traits, length(lines)), dimnames = list(unique(matrix(unlist(strsplit(colnames(table), "_._._", )), length(table), 2, byrow = TRUE)[, 1]), unique(matrix(unlist(strsplit(colnames(table), "_._._", )), length(table), 2, byrow = TRUE)[, 1]), lines))
    A_all <- aperm(A_all, c(2, 1, 3))
    return(A_all)

  } else if (component == "U") {

    n_traits <- ncol(table) / r
    U_r_all <- array(unlist(t(table)), dim = c(n_traits, r, length(lines)), dimnames = list(substr(colnames(table)[1:n_traits], 9, stringr::str_length(colnames(table)[1:n_traits])), paste0("SV", 1:r), lines))
    U_r_all <- aperm(U_r_all, c(2, 1, 3))
    return(U_r_all)

  } else if (component == "V") {

    if (!is.numeric(n_timepoints)) {

      print("n_timepoints is required")
      break

    } else {

      n_timepoints <- n_timepoints - 1
      V_r_all <- array(unlist(t(table)), dim = c(r, n_timepoints, length(lines)), dimnames = list(paste0("SV", 1:r), substr(colnames(table)[seq(1, r * n_timepoints, r)], 1, stringr::str_length(colnames(table)[seq(1, r * n_timepoints, r)]) - 8), lines))
      V_r_all <- aperm(V_r_all, c(2, 1, 3))
      return(V_r_all)

    }
  } else if (component == "sigma") {

    sigma_all <- array(unlist(t(table)), dim = c(r, r, length(lines)), dimnames = list(paste0("SV", 1:r), paste0("SV", 1:r), lines))
    sigma_all <- aperm(sigma_all, c(2, 1, 3))
    return(sigma_all)

  } else if (component == "A_tilde") {

    A_tilde_all <- array(unlist(t(table)), dim = c(r, r, length(lines)), dimnames = list(paste0("SV", 1:r), paste0("SV", 1:r), lines))
    A_tilde_all <- aperm(A_tilde_all, c(2, 1, 3))
    return(A_tilde_all)

  } else if (component == "Q") {

    Q_all <- array(unlist(t(table)), dim = c(r, r, length(lines)), dimnames = list(paste0("SV", 1:r), paste0("SV", 1:r), lines))
    Q_all <- aperm(Q_all, c(2, 1, 3))
    return(Q_all)

  } else if (component == "R") {

    R_all <- array(unlist(t(table)), dim = c(r, r, length(lines)), dimnames = list(paste0("SV", 1:r), paste0("SV", 1:r), lines))
    R_all <- aperm(R_all, c(2, 1, 3))
    return(R_all)

  } else if (component == "phi") {

    n_traits <- ncol(table) / r
    phi_all <- array(unlist(t(table)), dim = c(r, n_traits, length(lines)), dimnames = list(paste0("SV", 1:r), unique(substr(colnames(table), 1, stringr::str_length(colnames(table)) - 8)), lines))
    phi_all <- aperm(phi_all, c(2, 1, 3))
    return(phi_all)

  }
}


# length(which(table_to_tensor(tensor_unroller(DMD_components$U_r_transpose_all), component = "U",       r = 6) != DMD_components$U_r_transpose_all))
# length(which(table_to_tensor(tensor_unroller(DMD_components$V_r_all),           component = "V",       r = 6, n_timepoints = 25) != DMD_components$V_r_all))
# length(which(table_to_tensor(tensor_unroller(DMD_components$sigma_r_all),       component = "sigma",   r = 6) != DMD_components$sigma))
# length(which(table_to_tensor(tensor_unroller(DMD_components$A_tilde_all),       component = "A_tilde", r = 6) != DMD_components$A_tilde_all))
# length(which(table_to_tensor(tensor_unroller(DMD_components$Q_all),             component = "Q",       r = 6) != DMD_components$Q_all))
# length(which(table_to_tensor(tensor_unroller(DMD_components$R_all),             component = "R",       r = 6) != DMD_components$R_all))
# length(which(table_to_tensor(tensor_unroller(DMD_components$phi_all),           component = "phi",     r = 6) != DMD_components$phi_all))