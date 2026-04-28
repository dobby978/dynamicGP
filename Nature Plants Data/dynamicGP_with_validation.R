devtools::load_all("dynamicGP - private")

########################################################################################################################
# Heritability Analysis
########################################################################################################################

data <- read.csv("Dynamic_Mode_Decomposition_Project/Arabidopsis_Data/BLUE_Experiment_3_IAP_NA_imputed_selected_traits.csv")
data <- data[, -1]
# data <- cbind.data.frame(data[, 1], data[, 1], data[, -1])
# colnames(data)[c(1:3)] <- c("bio_ID", "accession", "DAS")

cv <- read.csv("Dynamic_Mode_Decomposition_Project/Arabidopsis_Data/cross_validation_folds_w_validation_set_arabidopsis.csv")

common_lines <- intersect(cv$X, unique(data$bio_ID))

data <- data[data$bio_ID %in% common_lines, ]

iterations <- 10

validation_heritability <- training_heritability <- data.frame()

dir.create("Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Heritability")

for (k in 1:10) {
  cv_loop <- cv[, c(1, k + 1)]
  rownames(cv_loop) <- cv_loop[, 1]

  validation_data <- data[data$bio_ID %in% cv_loop$X[which(cv_loop[, 2]  == 6)], ]
  validation_set_components <- find_DMD_components(validation_data, r = 2, method = "schur.DMD")

  validation_U <- tensor_unroller(validation_set_components$U_r_transpose_all)
  validation_U_h <- cbind.data.frame("U", "Validation", 2, k, GCTA(validation_U, kinship_filename = "GCTA_GRM/myGD_w_h_GRM_maf0.05_filtered", wd = "/home/david/PycharmProjects/PhD_Projekt/Dynamic_Mode_Decomposition_Project/Arabidopsis_Data"))
  validation_U_h$Element <- rownames(validation_U_h)
  rownames(validation_U_h) <- NULL
  colnames(validation_U_h) <- c("Component", "TV", "r", "Iteration", "Heritability", "Element")

  validation_V <- tensor_unroller(validation_set_components$V_r_all)
  validation_V_h <- cbind.data.frame("V", "Validation", 2, k, GCTA(validation_V, kinship_filename = "GCTA_GRM/myGD_w_h_GRM_maf0.05_filtered", wd = "/home/david/PycharmProjects/PhD_Projekt/Dynamic_Mode_Decomposition_Project/Arabidopsis_Data"))
  validation_V_h$Element <- rownames(validation_V_h)
  rownames(validation_V_h) <- NULL
  colnames(validation_V_h) <- c("Component", "TV", "r", "Iteration", "Heritability", "Element")

  validation_sigma <- tensor_unroller(validation_set_components$sigma_r_all)
  validation_sigma_h <- cbind.data.frame("sigma", "Validation", 2, k, GCTA(validation_sigma, kinship_filename = "GCTA_GRM/myGD_w_h_GRM_maf0.05_filtered", wd = "/home/david/PycharmProjects/PhD_Projekt/Dynamic_Mode_Decomposition_Project/Arabidopsis_Data"))
  validation_sigma_h$Element <- rownames(validation_sigma_h)
  rownames(validation_sigma_h) <- NULL
  colnames(validation_sigma_h) <- c("Component", "TV", "r", "Iteration", "Heritability", "Element")

  validation_A_tilde <- tensor_unroller(validation_set_components$A_tilde_all)
  validation_A_tilde_h <- cbind.data.frame("A_tilde", "Validation", 2, k, GCTA(validation_A_tilde, kinship_filename = "GCTA_GRM/myGD_w_h_GRM_maf0.05_filtered", wd = "/home/david/PycharmProjects/PhD_Projekt/Dynamic_Mode_Decomposition_Project/Arabidopsis_Data"))
  validation_A_tilde_h$Element <- rownames(validation_A_tilde_h)
  rownames(validation_A_tilde_h) <- NULL
  colnames(validation_A_tilde_h) <- c("Component", "TV", "r", "Iteration", "Heritability", "Element")

  validation_Q <- tensor_unroller(validation_set_components$Q_all)
  validation_Q_h <- cbind.data.frame("Q", "Validation", 2, k, GCTA(validation_Q, kinship_filename = "GCTA_GRM/myGD_w_h_GRM_maf0.05_filtered", wd = "/home/david/PycharmProjects/PhD_Projekt/Dynamic_Mode_Decomposition_Project/Arabidopsis_Data"))
  validation_Q_h$Element <- rownames(validation_Q_h)
  rownames(validation_Q_h) <- NULL
  colnames(validation_Q_h) <- c("Component", "TV", "r", "Iteration", "Heritability", "Element")

  validation_R <- tensor_unroller(validation_set_components$R_all)
  validation_R_h <- cbind.data.frame("R", "Validation", 2, k, GCTA(validation_R, kinship_filename = "GCTA_GRM/myGD_w_h_GRM_maf0.05_filtered", wd = "/home/david/PycharmProjects/PhD_Projekt/Dynamic_Mode_Decomposition_Project/Arabidopsis_Data"))
  validation_R_h$Element <- rownames(validation_R_h)
  rownames(validation_R_h) <- NULL
  colnames(validation_R_h) <- c("Component", "TV", "r", "Iteration", "Heritability", "Element")

  validation_phi <- tensor_unroller(validation_set_components$phi_all)
  validation_phi_h <- cbind.data.frame("phi", "Validation", 2, k, GCTA(validation_phi, kinship_filename = "GCTA_GRM/myGD_w_h_GRM_maf0.05_filtered", wd = "/home/david/PycharmProjects/PhD_Projekt/Dynamic_Mode_Decomposition_Project/Arabidopsis_Data"))
  validation_phi_h$Element <- rownames(validation_phi_h)
  rownames(validation_phi_h) <- NULL
  colnames(validation_phi_h) <- c("Component", "TV", "r", "Iteration", "Heritability", "Element")

  validation_heritability <- rbind.data.frame(validation_heritability, validation_U_h, validation_V_h, validation_sigma_h, validation_A_tilde_h, validation_Q_h, validation_R_h, validation_phi_h)

  training_data <- data[data$bio_ID %in% cv_loop$X[which(cv_loop[, 2] != 6)], ]

  for (i in 2:13) {
    training_components <- find_DMD_components(training_data, r = i, method = "schur.DMD")

    if (i == 13) {
      training_U <- tensor_unroller(training_components$U_r_transpose_all)
      training_U_h <- cbind.data.frame("U", "Training", i, k, GCTA(training_U, kinship_filename = "GCTA_GRM/myGD_w_h_GRM_maf0.05_filtered", wd = "/home/david/PycharmProjects/PhD_Projekt/Dynamic_Mode_Decomposition_Project/Arabidopsis_Data"))
      training_U_h$Element <- rownames(training_U_h)
      rownames(training_U_h) <- NULL
      colnames(training_U_h) <- c("Component", "TV", "r", "Iteration", "Heritability", "Element")

      training_V <- tensor_unroller(training_components$V_r_all)
      training_V_h <- cbind.data.frame("V", "Training", i, k, GCTA(training_V, kinship_filename = "GCTA_GRM/myGD_w_h_GRM_maf0.05_filtered", wd = "/home/david/PycharmProjects/PhD_Projekt/Dynamic_Mode_Decomposition_Project/Arabidopsis_Data"))
      training_V_h$Element <- rownames(training_V_h)
      rownames(training_V_h) <- NULL
      colnames(training_V_h) <- c("Component", "TV", "r", "Iteration", "Heritability", "Element")

      training_sigma <- tensor_unroller(training_components$sigma_r_all)
      training_sigma_h <- cbind.data.frame("sigma", "Training", i, k, GCTA(training_sigma, kinship_filename = "GCTA_GRM/myGD_w_h_GRM_maf0.05_filtered", wd = "/home/david/PycharmProjects/PhD_Projekt/Dynamic_Mode_Decomposition_Project/Arabidopsis_Data"))
      training_sigma_h$Element <- rownames(training_sigma_h)
      rownames(training_sigma_h) <- NULL
      colnames(training_sigma_h) <- c("Component", "TV", "r", "Iteration", "Heritability", "Element")

      training_A_tilde <- tensor_unroller(training_components$A_tilde_all)
      training_A_tilde_h <- cbind.data.frame("A_tilde", "Training", i, k, GCTA(training_A_tilde, kinship_filename = "GCTA_GRM/myGD_w_h_GRM_maf0.05_filtered", wd = "/home/david/PycharmProjects/PhD_Projekt/Dynamic_Mode_Decomposition_Project/Arabidopsis_Data"))
      training_A_tilde_h$Element <- rownames(training_A_tilde_h)
      rownames(training_A_tilde_h) <- NULL
      colnames(training_A_tilde_h) <- c("Component", "TV", "r", "Iteration", "Heritability", "Element")
    }

    training_Q <- tensor_unroller(training_components$Q_all)
    training_Q_h <- cbind.data.frame("Q", "Training", i, k, GCTA(training_Q, kinship_filename = "GCTA_GRM/myGD_w_h_GRM_maf0.05_filtered", wd = "/home/david/PycharmProjects/PhD_Projekt/Dynamic_Mode_Decomposition_Project/Arabidopsis_Data"))
    training_Q_h$Element <- rownames(training_Q_h)
    rownames(training_Q_h) <- NULL
    colnames(training_Q_h) <- c("Component", "TV", "r", "Iteration", "Heritability", "Element")

    training_R <- tensor_unroller(training_components$R_all)
    training_R_h <- cbind.data.frame("R", "Training", i, k, GCTA(training_R, kinship_filename = "GCTA_GRM/myGD_w_h_GRM_maf0.05_filtered", wd = "/home/david/PycharmProjects/PhD_Projekt/Dynamic_Mode_Decomposition_Project/Arabidopsis_Data"))
    training_R_h$Element <- rownames(training_R_h)
    rownames(training_R_h) <- NULL
    colnames(training_R_h) <- c("Component", "TV", "r", "Iteration", "Heritability", "Element")

    training_phi <- tensor_unroller(training_components$phi_all)
    training_phi_h <- cbind.data.frame("phi", "Training", i, k, GCTA(training_phi, kinship_filename = "GCTA_GRM/myGD_w_h_GRM_maf0.05_filtered", wd = "/home/david/PycharmProjects/PhD_Projekt/Dynamic_Mode_Decomposition_Project/Arabidopsis_Data"))
    training_phi_h$Element <- rownames(training_phi_h)
    rownames(training_phi_h) <- NULL
    colnames(training_phi_h) <- c("Component", "TV", "r", "Iteration", "Heritability", "Element")

    if (i == 13) {
      training_heritability <- rbind.data.frame(training_heritability, training_U_h, training_V_h, training_sigma_h, training_A_tilde_h, training_Q_h, training_R_h, training_phi_h)
    } else {
      training_heritability <- rbind.data.frame(training_heritability, training_Q_h, training_R_h, training_phi_h)
    }
  }
}

heritability <- rbind.data.frame(validation_heritability, training_heritability)
# write.csv(heritability, "Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Heritability/Heritability_all_folds_10iterations_arabidopsis.csv")


heritability <- read.csv("Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Heritability/Heritability_all_folds_10iterations_arabidopsis.csv")
heritability <- heritability[which(heritability$TV != "Validation"), -1]
heritability <- heritability[-which(is.na(heritability$Heritability)), ]


# heritability <- heritability[which(heritability$TV != "Validation"), ]

heritability$Method <- "Schur"
heritability$Method[which(heritability$Component == "U")] <- "SVD"
heritability$Method[which(heritability$Component == "V")] <- "SVD"
heritability$Method[which(heritability$Component == "sigma")] <- "SVD"
heritability$Method[which(heritability$Component == "A_tilde")] <- "SVD"

last <- substr(heritability$Element, stringr::str_length(heritability$Element) - 2, stringr::str_length(heritability$Element) - 1)
first <- substr(heritability$Element, 1, 2)

heritability$R <- NA
for (i in 1:nrow(heritability)) {
  if (heritability$Method[i] == "SVD") {
    if (first[i] == "SV" && last[i] == "SV") {
      heritability$R[i] <- as.numeric(max(substr(heritability$Element[i], 3, 3), substr(heritability$Element[i], stringr::str_length(heritability$Element[i]), stringr::str_length(heritability$Element[i]))))
    } else if (first[i] == "SV" && last[i] != "SV") {
      heritability$R[i] <- as.numeric(substr(heritability$Element[i], 3, 3))
    } else if (first[i] != "SV" && last[i] == "SV") {
      heritability$R[i] <- as.numeric(substr(heritability$Element[i], stringr::str_length(heritability$Element[i]), stringr::str_length(heritability$Element[i])))
    }
  } else {
    heritability$R[i] <- heritability$r[i]
  }
}

heritability$r <- heritability$R
heritability <- heritability[, -8]
# colnames(heritability)[5] <- "Component"
colnames(heritability)[5] <- "Heritability"
# heritability <- heritability[-which(is.na(heritability$Heritability) == TRUE), ]



title_size <- 7
axis_text_size <- 7

# heritability$r <- factor(heritability$r, levels = sort(unique(heritability$r)))
heritability$Component[which(heritability$Component == "A_tilde")] <- "\U00C3"
heritability$Component[which(heritability$Component == "sigma")] <- "\U03A3"
heritability$Component[which(heritability$Component == "phi")] <- "\U03D5"


heritability <- heritability[-which(as.numeric(heritability$r) > 9), ]
heritability <- heritability[-which(is.na(heritability$r)), ]

heritability$r <- factor(heritability$r, levels = 1:9)
heritability$Component <- factor(heritability$Component, levels = c("\U03D5", "Q", "R", "\U00C3", "U", "V", "\U03A3"))

library("ggplot2")

plot <- ggplot(data = heritability[which(heritability$Method == "SVD"), ], aes(y = Heritability, x = r, fill = Component)) +
  # geom_violin()
  geom_boxplot(position = position_dodge(), lwd = 0.25, outlier.size = 0.05) +
  scale_fill_manual(values = c('#00120B', '#6B818C', '#709BFF', "#17DE6D")) + # '#00120B',
  xlab("Corresponding singular vector") +
  #xlim(c(0, 6)) +
  theme_classic() +
  theme(legend.position = c(0.85, 0.85),
        # legend.box.background = element_rect(colour = "black"),
        legend.key.width = unit(0.35, 'cm'),
        legend.key.height = unit(0.35, "cm"),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        axis.title = element_text(size = title_size, face = "plain"),
        axis.text = element_text(size = axis_text_size)
  )

plot

ggsave(filename = "Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Plots/Figure S11 - Heritability of SVD components - after validation.png",
       plot = plot,
       width = 88,
       height = 88,
       units = "mm",
       dpi = 1000
)


plot <- ggplot(data = heritability[which(heritability$Method == "Schur"), ], aes(y = Heritability, x = r, fill = Component)) +
  # geom_violin()
  geom_boxplot(position = position_dodge(), lwd = 0.25, outlier.size = 0.05) +
  scale_fill_manual(values = c('#00120B', '#6B818C', '#709BFF', "#17DE6D")) + # '#00120B',
  xlab("r") +
  theme_classic() +
  theme(legend.position = c(0.85, 0.85),
        # legend.box.background = element_rect(colour = "black"),
        legend.key.width = unit(0.35, 'cm'),
        legend.key.height = unit(0.35, "cm"),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        axis.title = element_text(size = title_size, face = "plain"),
        axis.text = element_text(size = axis_text_size)
  )

plot

ggsave(filename = "Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Plots/Figure S12 - Heritability of Schur components - after validation.svg",
       plot = plot,
       width = 88,
       height = 88,
       units = "mm",
       dpi = 1000
)


########################################################################################################################
# Prediction of R and Phi
########################################################################################################################

CV_looper <- function (pheno_data, snps, cv, component) {
  trait_result_list <- vector(mode = "list", length = 4)
  # trait <- colnames(training_R)[1]

  # predicted_components <- data.frame()

  for (i in 1:ncol(pheno_data)) {

    print(paste0(i, " of ", ncol(pheno_data)))

    trait <- colnames(pheno_data)[i]
    y <- pheno_data[, i]
    names(y) <- rownames(pheno_data)

    fold_result_list <- vector(mode = "list", length = 5)

    for (j in 1:5) {
      training_idx <- rownames(cv)[which(cv[, 2] != j)]
      testing_idx  <- rownames(cv)[which(cv[, 2] == j)]

      training_y <- as.matrix(   y[training_idx])
      training_z <- as.matrix(snps[training_idx, ])

      testing_y <- as.matrix(   y[testing_idx])
      testing_z <- as.matrix(snps[testing_idx, ])

      model <- rrBLUP::mixed.solve(y = training_y, Z = training_z)

      marker_effects <- t(as.matrix(model$u))
      BLUE <- model$beta

      predicted_train <- as.matrix(training_z) %*% as.vector(marker_effects)
      predicted_train_result <- predicted_train + as.vector(BLUE)

      predicted_test <- as.matrix(testing_z) %*% as.vector(marker_effects)
      predicted_test_result <- predicted_test + as.vector(BLUE)

      train_cor <- cor(predicted_train_result, training_y, use = 'complete.obs')
      test_cor  <- cor(predicted_test_result,  testing_y,  use = 'complete.obs')

      print(test_cor)

      fold_result_list[[j]] <- list(Component = component, Trait = trait, Fold = j, Training_Accuracy = train_cor, Testing_Accuracy = test_cor, BLUE = BLUE, Marker_Effects = marker_effects)
    }

    trait_result_list[[i]] <- fold_result_list

  }

  return(trait_result_list)
}

validator <- function (pheno_data, snps, predictors, iteration, best_fold, component) {

  predictors <- predictors[which(predictors[, 3] == best_fold), ]
  trait_result_list <- vector(mode = "list", length = 4)
  # trait <- colnames(training_R)[1]

  predicted_component <- data.frame()

  for (i in 1:ncol(pheno_data)) {

    trait <- colnames(pheno_data)[i]
    y <- pheno_data[, i]
    names(y) <- rownames(pheno_data)

    # fold_result_list <- vector(mode = "list", length = 5)

    BLUE <- as.numeric(predictors[which(predictors[, 2] == trait), 4])

    marker_effects <- as.numeric(predictors[intersect(which(predictors[, 2] == trait), which(predictors[, 3] == best_fold)), -c(1:4)])

    predicted_test <- as.matrix(snps) %*% as.vector(marker_effects)
    predicted_test_result <- predicted_test + as.vector(BLUE)

    validation_cor <- cor(y, predicted_test_result)
    print(validation_cor)

    predicted_component <- rbind.data.frame(predicted_component, t(predicted_test_result))

          # fold_result_list[[j]] <- list(Component = component, Trait = trait, Fold = j, Training_Accuracy = train_cor, Testing_Accuracy = test_cor, BLUE = BLUE, Marker_Effects = marker_effects)

    trait_result_list[[i]] <- list(Component = component, Iteration = iteration, Trait = trait, Validation_Accuracy = validation_cor)

  }

  predicted_component <- t(predicted_component)
  colnames(predicted_component) <- colnames(pheno_data)
  write.csv(predicted_component, paste0("Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Training_Validation_Results/Validation_Results/Iteration_", iteration, "/", component, "_predicted_component_validated_iteration_", iteration, ".csv"))

  return(trait_result_list)
}

iterations <- 10

cv <- read.csv("Dynamic_Mode_Decomposition_Project/Arabidopsis_Data/cross_validation_folds_w_validation_set.csv")

snps <- read.table("Dynamic_Mode_Decomposition_Project/Arabidopsis_Data/myGD_w_h_maf0.05_filtered.txt")
# snps <- read.csv("Dynamic_Mode_Decomposition_Project/Arabidopsis_Data/myGD_w_h.csv")
# lines <- snps[, 1]
# # rownames(snps)
# snps <- data.matrix(snps, rownames.force = NA)
# rownames(snps) <- lines
# snps <- snps - 1

dir.create("Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Training_Validation_Results")
dir.create("Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Training_Validation_Results/Training_Prediction_Accuracy")
dir.create("Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Training_Validation_Results/Validation_Prediction_Accuracy")
dir.create("Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Training_Validation_Results/Training_Results")
dir.create("Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Training_Validation_Results/Validation_Results")
dir.create("Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Training_Validation_Results/Training_Marker_Effects")

for (k in 1:iterations) {
  dir.create(paste0("Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Training_Validation_Results/Training_Prediction_Accuracy/Iteration_", k))
  dir.create(paste0("Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Training_Validation_Results/Training_Results/Iteration_", k))
  dir.create(paste0("Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Training_Validation_Results/Training_Marker_Effects/Iteration_", k))

  cv_loop <- cv[, c(1, k + 1)]
  rownames(cv_loop) <- cv_loop[, 1]
  training_data <- data[data$bio_ID %in% cv_loop$X[which(cv_loop[, 2] != 6)], ]

  training_components <- find_DMD_components(training_data, r = 2, method = "schur.DMD")
  training_R       <- tensor_unroller(training_components$R_all)
  training_phi     <- tensor_unroller(training_components$phi_all)

  snps_training   <- snps[rownames(cv_loop)[which(cv_loop[, 2] != 6)], ]

  R_PAs <- CV_looper(training_R, snps, cv_loop[-which(cv_loop[, 2] == 6), ], "R")
  R_PAs <- matrix(unlist(R_PAs), ncol = (6 + ncol(snps)), byrow = TRUE)
  R_MEs <- R_PAs[, -c(4:5)]
  R_PAs <- as.data.frame(R_PAs[, 1:5])
  colnames(R_PAs) <- c("Component", "Element", "Fold", "Training_Accuracy", "Testing_Accuracy")
  R_PAs$Testing_Accuracy <- as.numeric(R_PAs$Testing_Accuracy)

  write.csv(R_PAs, paste0("Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Training_Validation_Results/Training_Prediction_Accuracy/Iteration_", k, "/R_Prediction_Accuracy_iteration_", k, ".csv"))
  write.csv(R_MEs, paste0("Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Training_Validation_Results/Training_Marker_Effects/Iteration_",      k, "/R_Marker_Effects_iteration_",      k, ".csv"))

  phi_PAs <- CV_looper(training_phi, snps, cv_loop[-which(cv_loop[, 2] == 6), ], "phi")
  phi_PAs <- matrix(unlist(phi_PAs), ncol = (6 + ncol(snps)), byrow = TRUE)
  phi_MEs <- phi_PAs[, -c(4:5)]
  phi_PAs <- as.data.frame(phi_PAs[, 1:5])
  colnames(phi_PAs) <- c("Component", "Element", "Fold", "Training_Accuracy", "Testing_Accuracy")
  phi_PAs$Testing_Accuracy <- as.numeric(phi_PAs$Testing_Accuracy)

  write.csv(phi_PAs, paste0("Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Training_Validation_Results/Training_Prediction_Accuracy/Iteration_", k, "/phi_Prediction_Accuracy_iteration_", k, ".csv"))
  write.csv(phi_MEs, paste0("Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Training_Validation_Results/Training_Marker_Effects/Iteration_",      k, "/phi_Marker_Effects_iteration_",      k, ".csv"))

  PAs <- rbind.data.frame(R_PAs, phi_PAs)
  PAs_agg <- aggregate(PAs$Testing_Accuracy, by = list(PAs$Fold), FUN = mean)
  colnames(PAs_agg) <- c("Fold", "Accuracy")

  write.csv(PAs_agg, paste0("Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Training_Validation_Results/Training_Prediction_Accuracy/Iteration_", k, "/Aggregated_Prediction_Accuracy_iteration_", k, ".csv"))

  rm(snps_training, R_PAs, phi_PAs)
  gc()

  best_fold <- which.max(PAs_agg$Accuracy)

  dir.create(paste0("Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Training_Validation_Results/Validation_Prediction_Accuracy/Iteration_", k))
  dir.create(paste0("Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Training_Validation_Results/Validation_Results/Iteration_", k))

  snps_validation <- snps[rownames(cv_loop)[which(cv_loop[, 2] == 6)], ]

  validation_data <- data[data$bio_ID %in% cv_loop$X[which(cv_loop[, 2] == 6)], ]
  validation_set_components <- find_DMD_components(validation_data, r = 2, method = "schur.DMD")

  validation_R <- tensor_unroller(validation_set_components$R_all)
  validation_phi <- tensor_unroller(validation_set_components$phi_all)

  # if the script crashes this loads for the saved marker effects
  best_fold <- read.csv(paste0("Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Training_Validation_Results/Training_Prediction_Accuracy/Iteration_", k, "/Aggregated_Prediction_Accuracy_iteration_", k, ".csv"))
  best_fold <- which.max(best_fold$Accuracy)
  R_MEs <- read.csv(paste0("Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Training_Validation_Results/Training_Marker_Effects/Iteration_", k, "/R_Marker_Effects_iteration_", k, ".csv"))
  R_MEs <- R_MEs[, -1]
  phi_MEs <- read.csv(paste0("Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Training_Validation_Results/Training_Marker_Effects/Iteration_", k, "/phi_Marker_Effects_iteration_", k, ".csv"))
  phi_MEs <- phi_MEs[, -1]

  R_validated <- validator(validation_R, snps_validation, R_MEs, iteration = k, best_fold, "R")
  phi_validated <- validator(validation_phi, snps_validation, phi_MEs, iteration = k, best_fold, "phi")

  R_validation_PAs   <- matrix(unlist(R_validated),   ncol = 4, byrow = TRUE)
  phi_validation_PAs <- matrix(unlist(phi_validated), ncol = 4, byrow = TRUE)

  rm(R_MEs, phi_MEs, snps_validation, validation_data, validation_set_components, validation_R, validation_phi, R_validated, phi_validated)

  write.csv(R_validation_PAs,   paste0("Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Training_Validation_Results/Validation_Prediction_Accuracy/Iteration_", k, "/R_Prediction_Accuracy_iteration_",   k, ".csv"))
  write.csv(phi_validation_PAs, paste0("Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Training_Validation_Results/Validation_Prediction_Accuracy/Iteration_", k, "/phi_Prediction_Accuracy_iteration_", k, ".csv"))

  rm(R_validation_PAs, phi_validation_PAs)
  gc()
}


########################################################################################################################
# rrBLUP benchmark predictions
########################################################################################################################


rrBLUP_validator <- function (pheno_data, snps, predictors, iteration, best_fold, component, n_timepoints, n_traits) {

  predictors <- predictors[which(predictors[, 3] == best_fold), ]
  trait_result_list <- vector(mode = "list", length = ncol(pheno_data))
  # trait <- colnames(training_R)[1]

  predicted_component <- data.frame()
  counter <- 1

  for (i in seq(1, n_timepoints * n_traits, n_timepoints)) {
  # for (i in 1:ncol(pheno_data[seq(1, n_timepoints * n_traits, n_timepoints)])) {

    trait_length <- stringr::str_length(colnames(pheno_data)[i])
    trait <- substr(colnames(pheno_data)[i], 1, trait_length - 8)

    single_trait_time_series <- pheno_data[, which(substr(colnames(pheno_data), 1, trait_length - 8) == trait)]

    BLUE <- as.numeric(predictors[which(substr(predictors[, 2], 1, trait_length - 8) == trait), 4])

    marker_effects <- as.numeric(predictors[intersect(which(substr(predictors[, 2], 1, trait_length - 8) == trait), which(predictors[, 3] == best_fold)), -c(1:4)])

    predicted_test <- as.matrix(snps) %*% as.vector(marker_effects)
    predicted_test_result <- predicted_test + as.vector(BLUE)

    for (j in 1:ncol(single_trait_time_series)) {

      print(paste0(counter, " of ", n_timepoints * n_traits))


      y <- single_trait_time_series[, j]
      names(y) <- rownames(single_trait_time_series)

      validation_cor <- cor(y, predicted_test_result)
      print(validation_cor)

      predicted_component <- rbind.data.frame(predicted_component, t(predicted_test_result))

      trait_result_list[[counter]] <- list(Component = component, Iteration = iteration, Trait = colnames(single_trait_time_series)[j], Validation_Accuracy = validation_cor)

      counter <- counter + 1
    }
  }

  predicted_component <- t(predicted_component)
  colnames(predicted_component) <- colnames(pheno_data)
  write.csv(predicted_component, paste0("Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Training_Validation_Results/rrBLUP_benchmark_Validation_Results/Iteration_", iteration, "/rrBLUP_benchmark_", component, "_predicted_component_validated_iteration_", iteration, ".csv"))

  return(trait_result_list)
}

# data <- read.csv("~/PycharmProjects/PhD_Projekt/Dynamic_Mode_Decomposition_Project/Maize_Data/trait/maize_trait_BLUP_50_NArmv_traits_w0SDrmv_afterClustering.csv")[, -1]
# days <- unique(data$DAS)

cv <- read.csv("Dynamic_Mode_Decomposition_Project/Arabidopsis_Data/cross_validation_folds_w_validation_set.csv")

snps <- read.csv("Dynamic_Mode_Decomposition_Project/Arabidopsis_Data/myGD_w_h.csv")
lines <- snps[, 1]
# rownames(snps)
snps <- data.matrix(snps, rownames.force = NA)
rownames(snps) <- lines
snps <- snps - 1

iterations <- 10

dir.create("Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Training_Validation_Results")
dir.create("Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Training_Validation_Results/rrBLUP_benchmark_Training_Prediction_Accuracy")
dir.create("Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Training_Validation_Results/rrBLUP_benchmark_Validation_Prediction_Accuracy")
dir.create("Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Training_Validation_Results/rrBLUP_benchmark_Training_Results")
dir.create("Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Training_Validation_Results/rrBLUP_benchmark_Validation_Results")
dir.create("Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Training_Validation_Results/rrBLUP_benchmark_Training_Marker_Effects")

for (k in 1:iterations) {
  dir.create(paste0("Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Training_Validation_Results/rrBLUP_benchmark_Training_Prediction_Accuracy/Iteration_", k))
  dir.create(paste0("Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Training_Validation_Results/rrBLUP_benchmark_Training_Results/Iteration_", k))
  dir.create(paste0("Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Training_Validation_Results/rrBLUP_benchmark_Training_Marker_Effects/Iteration_", k))

  cv_loop <- cv[, c(1, k + 1)]
  rownames(cv_loop) <- cv_loop[, 1]
  training_data <- data[data$bio_ID %in% cv_loop$X[which(cv_loop[, 2] != 6)], ]
  training_components <- find_DMD_components(training_data, r = 2, method = "schur.DMD")
  training_X <- tensor_unroller(training_components$X_all)

  snps_training   <- snps[rownames(cv_loop)[which(cv_loop[, 2] != 6)], ]

  X_PAs <- CV_looper(training_X[seq(1, 450, 15)], snps, cv_loop, "RR-BLUP benchmark")
  X_PAs <- matrix(unlist(X_PAs), ncol = (6 + ncol(snps)), byrow = TRUE)
  X_MEs <- X_PAs[, -c(4:5)]
  X_PAs <- as.data.frame(X_PAs[, 1:5])
  colnames(X_PAs) <- c("Component", "Element", "Fold", "Training_Accuracy", "Testing_Accuracy")
  X_PAs$Testing_Accuracy <- as.numeric(X_PAs$Testing_Accuracy)

  write.csv(X_PAs, paste0("Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Training_Validation_Results/rrBLUP_benchmark_Training_Prediction_Accuracy/Iteration_", k, "/rrBLUP_benchmark_Prediction_Accuracy_iteration_", k, ".csv"))
  write.csv(X_MEs, paste0("Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Training_Validation_Results/rrBLUP_benchmark_Training_Marker_Effects/Iteration_",      k, "/rrBLUP_benchmark_Marker_Effects_iteration_",      k, ".csv"))

  X_PAs_agg <- aggregate(X_PAs$Testing_Accuracy, by = list(X_PAs$Fold), FUN = mean)
  colnames(X_PAs_agg) <- c("Fold", "Accuracy")

  best_fold <- which.max(X_PAs_agg$Accuracy)

  write.csv(X_PAs_agg, paste0("Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Training_Validation_Results/rrBLUP_benchmark_Training_Prediction_Accuracy/Iteration_", k, "/rrBLUP_benchmark_Aggregated_Prediction_Accuracy_iteration_", k, ".csv"))

  rm(X_PAs, X_PAs_agg, snps_training, training_data, training_components, training_X)

  dir.create(paste0("Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Training_Validation_Results/rrBLUP_benchmark_Validation_Prediction_Accuracy/Iteration_", k))
  dir.create(paste0("Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Training_Validation_Results/rrBLUP_benchmark_Validation_Results/Iteration_", k))

  snps_validation <- snps[rownames(cv_loop)[which(cv_loop[, 2] == 6)], ]

  validation_data <- data[data$bio_ID %in% cv_loop$X[which(cv_loop[, 2] == 6)], ]
  validation_set_components <- find_DMD_components(validation_data, r = 2, method = "schur.DMD")

  validation_X <- tensor_unroller(validation_set_components$X_all)

  X_validated <- rrBLUP_validator(validation_X, snps_validation, X_MEs, iteration = k, best_fold, "rrBLUP_benchmark", n_timepoints = 15, n_traits = 30)

  X_validation_PAs   <- matrix(unlist(X_validated),   ncol = 4, byrow = TRUE)

  # rm(R_MEs, phi_MEs, snps_validation, validation_data, validation_set_components, validation_R, validation_phi, R_validated, phi_validated)

  write.csv(X_validation_PAs,   paste0("Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Training_Validation_Results/rrBLUP_benchmark_Validation_Prediction_Accuracy/Iteration_", k, "/rrBLUP_benchmark_Aggregated_Prediction_Accuracy_iteration_", k, ".csv"))

  rm(snps_validation, validation_data, validation_set_components, validation_X, X_validated, X_validation_PAs)
}


########################################################################################################################
# dynamicGP predictions
########################################################################################################################

devtools::load_all("dynamicGP - private")

cv <- read.csv("Dynamic_Mode_Decomposition_Project/Arabidopsis_Data/cross_validation_folds_w_validation_set.csv")

iterations <- 10

res_all <- res_rec_all <- MSE_all <- MSE_rec_all <- data.frame()

for (k in 1:iterations) {

  cv_loop <- cv[, c(1, k + 1)]
  rownames(cv_loop) <- cv_loop[, 1]

  R_all <- read.csv(paste0("Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Run_1_45_selected_traits/Training_Validation_Results/Validation_Results/Iteration_", k, "/R_predicted_component_validated_iteration_", k, ".csv"))
  rownames(R_all) <- R_all[, 1]
  R_all <- R_all[, -1]

  phi_all <- read.csv(paste0("Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Run_1_45_selected_traits/Training_Validation_Results/Validation_Results/Iteration_", k, "/phi_predicted_component_validated_iteration_", k, ".csv"))
  rownames(phi_all) <- phi_all[, 1]
  phi_all <- phi_all[, -1]

  R_all <- table_to_tensor(R_all, component = "R")
  phi_all <- table_to_tensor(phi_all, component = "phi")

  A_r_all <- A_r_from_phi_and_R(phi_all, R_all)

  validation_data <- data[data$bio_ID %in% cv_loop$X[which(cv_loop[, 2] == 6)], ]
  validation_set_components <- find_DMD_components(validation_data, r = 2, method = "schur.DMD")

  X_all_validation <- validation_set_components$X_all

  res <- all_timepoints_predictor(A_r_all, X_all_validation)
  res_melt     <- reshape2::melt(res$res_trait_by_timepoint)
  res_melt_rec <- reshape2::melt(res$res_trait_by_timepoint_rec)

  res_all <- rbind.data.frame(res_all, res_melt)
  res_rec_all <- rbind.data.frame(res_rec_all, res_melt_rec)

  MSE_melt     <- reshape2::melt(res$MSE_trait_by_timepoint)
  MSE_melt_rec <- reshape2::melt(res$MSE_trait_by_timepoint_rec)

  MSE_all <- rbind.data.frame(MSE_all, MSE_melt)
  MSE_rec_all <- rbind.data.frame(MSE_rec_all, MSE_melt_rec)
}

colnames(res_all) <- colnames(res_rec_all) <- c("Trait", "Time", "Accuracy")
res_all$Method <- "iterative"
res_rec_all$Method <- "recursive"

colnames(MSE_all) <- colnames(MSE_rec_all) <- c("Trait", "Time", "MSE")
MSE_all$Method <- "iterative"
MSE_rec_all$Method <- "recursive"

res_benchmark <- data.frame()

for (k in 1:iterations) {
  res <- read.csv(paste0("Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Run_1_45_selected_traits/Training_Validation_Results/rrBLUP_benchmark_Validation_Prediction_Accuracy/Iteration_", k, "/rrBLUP_benchmark_Aggregated_Prediction_Accuracy_iteration_", k, ".csv"))
  res_benchmark <- rbind.data.frame(res_benchmark, res)
}

res_benchmark <- cbind.data.frame(res_benchmark[, 2:3], stringr::str_split_fixed(res_benchmark[, 4], "_-_-_", 2), res_benchmark[, 5])
res_benchmark[, 4] <- substr(res_benchmark[, 4], 2, 3)

colnames(res_benchmark) <- c("Method", "Iteration", "Trait", "Time", "Accuracy")

res_benchmark_agg <- aggregate(res_benchmark$Accuracy, by = list(res_benchmark$Time, res_benchmark$Method), FUN = mean)
res_benchmark_agg$SD <- aggregate(res_benchmark$Accuracy, by = list(res_benchmark$Time, res_benchmark$Method), FUN = sd)[, 3]

colnames(res_benchmark_agg) <- c("Time", "Method", "Mean", "SD")


prediction_accuracies <- rbind.data.frame(res_all, res_rec_all)
# write.csv(prediction_accuracies, "Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Run_1_45_selected_traits/Training_Validation_Results/dynamicGP_prediction_accuracies_arabidopsis.csv")

  accuracy_plot <- ggplot2::ggplot(data = prediction_accuracies, ggplot2::aes(x = Time, y = Accuracy, fill = Method)) +
    ggplot2::geom_boxplot(position = ggplot2::position_dodge(), lwd = 0.25) +
    # ggplot2::ylim(-0.5, 1) +
    ggplot2::theme_classic() +
    ggplot2::ylab("Prediction accuracy") +
    ggplot2::scale_fill_manual(values = c('#709BFF', "#17DE6D")) + # '#00120B',
    ggplot2::theme(#legend.position = c(1-0.075, 0.15), #legend.position = c(0.9, 0.15),
          legend.title = ggplot2::element_text(size = 7),
          legend.text = ggplot2::element_text(size = 6),
          axis.title = ggplot2::element_text(size = 7, face = "plain"),
          axis.text = ggplot2::element_text(size = 6),
          legend.box.background = ggplot2::element_rect(colour = "black")
         )

prediction_accuracies_agg <- aggregate(prediction_accuracies$Accuracy, by = list(prediction_accuracies$Time, prediction_accuracies$Method), FUN = mean)
prediction_accuracies_agg$SD <- aggregate(prediction_accuracies$Accuracy, by = list(prediction_accuracies$Time, prediction_accuracies$Method), FUN = sd)[, 3]

colnames(prediction_accuracies_agg) <- c("Time", "Method", "Mean", "SD")

prediction_accuracies_agg <- microeco::dropallfactors(prediction_accuracies_agg)

prediction_accuracies_agg$Time[which(prediction_accuracies_agg$Time == "t1")] <- 8
prediction_accuracies_agg$Time[which(prediction_accuracies_agg$Time == "t2")] <- 9
prediction_accuracies_agg$Time[which(prediction_accuracies_agg$Time == "t3")] <- 10
prediction_accuracies_agg$Time[which(prediction_accuracies_agg$Time == "t4")] <- 11
prediction_accuracies_agg$Time[which(prediction_accuracies_agg$Time == "t5")] <- 12
prediction_accuracies_agg$Time[which(prediction_accuracies_agg$Time == "t6")] <- 13
prediction_accuracies_agg$Time[which(prediction_accuracies_agg$Time == "t7")] <- 14
prediction_accuracies_agg$Time[which(prediction_accuracies_agg$Time == "t8")] <- 15
prediction_accuracies_agg$Time[which(prediction_accuracies_agg$Time == "t9")] <- 16
prediction_accuracies_agg$Time[which(prediction_accuracies_agg$Time == "t10")] <-17
prediction_accuracies_agg$Time[which(prediction_accuracies_agg$Time == "t11")] <-18
prediction_accuracies_agg$Time[which(prediction_accuracies_agg$Time == "t12")] <-19
prediction_accuracies_agg$Time[which(prediction_accuracies_agg$Time == "t13")] <-20
prediction_accuracies_agg$Time[which(prediction_accuracies_agg$Time == "t14")] <-21
# prediction_accuracies_agg$Time[which(prediction_accuracies_agg$Time == "t15")] <-21


prediction_accuracies_agg$Linetype <- "dynamicGP"
res_benchmark_agg$Linetype <- "rrBLUP"

prediction_accuracies_agg <- rbind.data.frame(prediction_accuracies_agg, res_benchmark_agg)


prediction_accuracies_agg <- prediction_accuracies_agg[-which(prediction_accuracies_agg$Time == 8), ]

prediction_accuracies_agg$Time <-factor(prediction_accuracies_agg$Time, levels = unique(prediction_accuracies_agg$Time))

shades <- data.frame(xmin=seq(1, length(unique(prediction_accuracies_agg$Time)),     2) - 0.5,
                     xmax=seq(1, length(unique(prediction_accuracies_agg$Time)) , 2) + 0.5,
                     ymin=-Inf, ymax=Inf)

plot <- ggplot(data = prediction_accuracies_agg, aes(x = Time, y = Mean, colour = Method)) + # , group = Method
  # geom_boxplot(stat = "identity") +
  # ggstats::geom_stripped_cols() +
  # ggforestplot::geom_stripes(stroke = 2) +
  geom_rect(inherit.aes = F, data = shades, mapping = aes(xmin=xmin, xmax=xmax, ymin = ymin, ymax = ymax), alpha = 0.2) +
  geom_line(stat = "identity", aes(linetype = Linetype, group = interaction(Method)), size = 0.5, position = position_dodge(.9)) + # interaction(Method, Trait))
  geom_point(stat = "identity", position = position_dodge(.9), stroke = 0.25) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 1,
                position = position_dodge(.9)) +
  scale_color_manual(values = c( '#709BFF', "#17DE6D", '#00120B'), labels = c("iterative", "recursive", "baseline")) +
  ylim(c(-0.2, 0.6)) +
  xlab("DAS") +
  theme_classic() +
  ylab("Prediction accuracy") +
  theme(legend.position = c(0.2, 0.08),
        legend.box = "horizontal",
        legend.key.height = unit(0.2, "cm"),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 7),
        legend.title = element_text(size = 5, face = "bold"),
        #title = element_text(size = title_size),
        legend.text = element_text(size = 5)) +
        #legend.box.background = element_rect(colour = "black")) +
  guides(colour = guide_legend(title = "Method", ncol = 3),
         linetype = FALSE)

plot

ggsave(filename = "Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Run_1_45_selected_traits/Training_Validation_Results/Figure S13 - Longitudinal Prediction Accuracies.svg",
       plot = plot,
       width = 180,
       height = 90,
       units = "mm",
       dpi = 1000
)

ggsave(filename = "Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Run_1_45_selected_traits/Training_Validation_Results/Figure S13 - Longitudinal Prediction Accuracies.png",
       plot = plot,
       width = 180,
       height = 90,
       units = "mm",
       dpi = 1000
)











iterations <- 10

MSE <- matrix(NA, iterations, ncol(validation_X))

for (k in 1:iterations) {

  cv_loop <- cv[, c(1, k + 1)]
  rownames(cv_loop) <- cv_loop[, 1]

  pred <- read.csv(paste0("Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Run_1_45_selected_traits/Training_Validation_Results/rrBLUP_benchmark_Validation_Results/Iteration_", k, "/rrBLUP_benchmark_rrBLUP_benchmark_predicted_component_validated_iteration_", k, ".csv"))
  rownames(pred) <- pred[, 1]
  pred <- pred[, -1]

  validation_data <- data[data$bio_ID %in% cv_loop$X[which(cv_loop[, 2] == 6)], ]
  validation_set_components <- find_DMD_components(validation_data, r = 2, method = "schur.DMD")

  validation_X <- tensor_unroller(validation_set_components$X_all)

  colnames(MSE) <- colnames(validation_X)
  # cor(pred[, i], validation_X[, i])
  for (i in 1:ncol(pred)) {
    MSE[k, i] <- mean((validation_X[, i] - pred[, i]) ^ 2)
  }
}

MSE_benchmark <- reshape2::melt(MSE)
MSE_benchmark <- cbind.data.frame(MSE_benchmark[, 1], stringr::str_split_fixed(MSE_benchmark[, 2], "_-_-_", 2), MSE_benchmark[, 3])
colnames(MSE_benchmark) <- c("Iteration", "Trait", "Time", "MSE")



MSEs <- rbind.data.frame(MSE_all, MSE_rec_all)
MSEs <- microeco::dropallfactors(MSEs)

MSE_benchmark$Iteration <- "benchmark"
colnames(MSE_benchmark)[1] <- "Method"

MSE_benchmark$Time[which(MSE_benchmark$Time == "t8")]  <- 8
MSE_benchmark$Time[which(MSE_benchmark$Time == "t9")]  <- 9
MSE_benchmark$Time[which(MSE_benchmark$Time == "t10")] <- 10
MSE_benchmark$Time[which(MSE_benchmark$Time == "t11")] <- 11
MSE_benchmark$Time[which(MSE_benchmark$Time == "t12")] <- 12
MSE_benchmark$Time[which(MSE_benchmark$Time == "t13")] <- 13
MSE_benchmark$Time[which(MSE_benchmark$Time == "t14")] <- 14
MSE_benchmark$Time[which(MSE_benchmark$Time == "t15")] <- 15
MSE_benchmark$Time[which(MSE_benchmark$Time == "t16")] <- 16
MSE_benchmark$Time[which(MSE_benchmark$Time == "t17")] <- 17
MSE_benchmark$Time[which(MSE_benchmark$Time == "t18")] <- 18
MSE_benchmark$Time[which(MSE_benchmark$Time == "t19")] <- 19
MSE_benchmark$Time[which(MSE_benchmark$Time == "t20")] <- 20
MSE_benchmark$Time[which(MSE_benchmark$Time == "t21")] <- 21

MSEs$Time[which(MSEs$Time == "t1")] <- 8
MSEs$Time[which(MSEs$Time == "t2")] <- 9
MSEs$Time[which(MSEs$Time == "t3")] <- 10
MSEs$Time[which(MSEs$Time == "t4")] <- 11
MSEs$Time[which(MSEs$Time == "t5")] <- 12
MSEs$Time[which(MSEs$Time == "t6")] <- 13
MSEs$Time[which(MSEs$Time == "t7")] <- 14
MSEs$Time[which(MSEs$Time == "t8")] <- 15
MSEs$Time[which(MSEs$Time == "t9")] <- 16
MSEs$Time[which(MSEs$Time == "t10")] <- 17
MSEs$Time[which(MSEs$Time == "t11")] <- 18
MSEs$Time[which(MSEs$Time == "t12")] <- 19
MSEs$Time[which(MSEs$Time == "t13")] <- 20
MSEs$Time[which(MSEs$Time == "t14")] <- 21

MSEs <- rbind.data.frame(MSEs, MSE_benchmark)


MSEs <- MSEs[which(MSEs$Time != 8), ]

MSEs$Time <- factor(MSEs$Time, levels = 9:21)

shades <- data.frame(xmin=seq(1, length(unique(MSEs$Time)),     2) - 0.5,
                     xmax=seq(1, length(unique(MSEs$Time)), 2) + 0.5,
                     ymin=-Inf, ymax=Inf)

plot <- ggplot(data = MSEs, aes(x = Time, y = log(MSE), fill = Method)) + # , group = Method
  geom_boxplot(lwd = 0.25, outlier.size = 0.05) +
  geom_rect(inherit.aes = F, data = shades, mapping = aes(xmin=xmin, xmax=xmax, ymin = ymin, ymax = ymax), alpha = 0.2) +
  scale_fill_manual(values = c('#696b77', '#709BFF', "#17DE6D"), labels = c("baseline", "iterative", "recursive")) +
  # ylim(c(-0.25, 0.75)) +
  theme_classic() +
  ylab("log mean squared error") +
  xlab("DAS") +
  theme(legend.position = c(0.8, 0.08),
        legend.box = "horizontal",
        legend.key.height = unit(0.2, "cm"),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 7),
        legend.title = element_text(size = 5, face = "bold"),
        #title = element_text(size = title_size),
        legend.text = element_text(size = 5)) +
        #legend.box.background = element_rect(colour = "black")) +
  guides(fill = guide_legend(title = "Method", ncol = 3),
         linetype = FALSE)

plot

ggsave(filename = "Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Plots/Figure S18 - Longitudinal MSE.svg",
       plot = plot,
       width = 180,
       height = 90,
       units = "mm",
       dpi = 1000
)

ggsave(filename = "Dynamic_Mode_Decomposition_Project/Maize_Data/trait/DMD_4/Arabidopsis_Data/Plots/Figure S18 - Longitudinal MSE.png",
       plot = plot,
       width = 180,
       height = 90,
       units = "mm",
       dpi = 1000
)