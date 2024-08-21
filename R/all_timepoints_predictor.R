#' Run full time series predictions
#'
#' Function takes two inputs, the A_all and X_all objects output by find_DMD_compoents(), and outputs iterative and recursive predictions for each trait and each time point.

all_timepoints_predictor <- function (A_all, X_all) {

  lines <- unlist(dimnames(A_all)[3])
  traits <- unlist(dimnames(A_all)[1])
  timepoints <- dim(X_all)[2]
  res_all <- res_all_rec <- array(NA, dim = c(length(traits), timepoints, length(lines)), dimnames = list(colnames(data)[-c(1, 2, 3)], 1:timepoints, lines))

  for (line in lines) {

    X <- X_all[, , line]
    A <- A_all[, , line]

    mat <- mat_rec <- matrix(0, nrow(X), timepoints) # + 1)
    mat[, 1] <- mat_rec[, 1] <- X[, 1]

    for (t in 1:(timepoints - 1)) {

          mat[, t + 1] <-     Re(A %*%       X[, t])
      mat_rec[, t + 1] <-     Re(A %*% mat_rec[, t])
    }

        res_all[, , line] <- mat
    res_all_rec[, , line] <- mat_rec
  }

  res_trait_by_timepoint <- res_trait_by_timepoint_rec <- matrix(0, length(traits), timepoints)
  rownames(res_trait_by_timepoint) <- rownames(res_trait_by_timepoint_rec) <- traits
  colnames(res_trait_by_timepoint) <- colnames(res_trait_by_timepoint_rec) <- 1:timepoints

  for (trait in traits) {
    for (t in 1:timepoints) {

      print(paste0("Correlation for ", trait, " at time point ", t))
      c     <- cor(res_all[trait, t, ],     X_all[trait, t, ], use = "complete.obs")
      c_rec <- cor(res_all_rec[trait, t, ], X_all[trait, t, ], use = "complete.obs")

          res_trait_by_timepoint[trait, t] <- c
      res_trait_by_timepoint_rec[trait, t] <- c_rec
    }
  }

  res_by_trait     <- rowMeans(res_trait_by_timepoint[, -1])
  res_by_trait_rec <- rowMeans(res_trait_by_timepoint_rec[, -1])
  res_by_timepont     <- colMeans(res_trait_by_timepoint[, -1])
  res_by_timepont_rec <- colMeans(res_trait_by_timepoint_rec[, -1])

  colnames(res_trait_by_timepoint) <- colnames(res_trait_by_timepoint_rec) <- paste0("t", 1:timepoints)
  return (list(res_trait_by_timepoint = res_trait_by_timepoint,
               res_trait_by_timepoint_rec = res_trait_by_timepoint_rec,
               res_by_trait = res_by_trait,
               res_by_trait_rec = res_by_trait_rec,
               all_trait_mean = mean(res_by_trait, na.rm = TRUE),
               all_trait_sd = sd(res_by_trait, na.rm = TRUE),
               all_trait_mean_rec = mean(res_by_trait_rec, na.rm = TRUE),
               all_trait_sd_rec = sd(res_by_trait_rec, na.rm = TRUE)))
}
