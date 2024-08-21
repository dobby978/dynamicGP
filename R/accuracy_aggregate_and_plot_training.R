#' Prediction accuracies aggregator
#'
#' This function takes as input the A_all and X_all tensor objects output by find_DMD_components() and aggregates the iterative and recursive prediction accuracies into a single data frame as well as a longitudinal plot of prediction accuracies across all traits.

accuracy_aggregate_and_plot_training <- function (A_all, X_all) {

  accuracy_all_timepoints <- all_timepoints_predictor(A_all, X_all)

  accuracy_all_timepoints_melt <- reshape2::melt(accuracy_all_timepoints$res_trait_by_timepoint[complete.cases(accuracy_all_timepoints$res_trait_by_timepoint), ])
  accuracy_all_timepoints_melt$Method <- "iterative"

  accuracy_all_timepoints_rec_melt <- reshape2::melt(accuracy_all_timepoints$res_trait_by_timepoint_rec[complete.cases(accuracy_all_timepoints$res_trait_by_timepoint_rec), ])
  accuracy_all_timepoints_rec_melt$Method <- "recursive"

  accuracies_true <- rbind.data.frame(accuracy_all_timepoints_melt, accuracy_all_timepoints_rec_melt)
  colnames(accuracies_true)[1:3] <- c("Trait", "Time", "Accuracy")

  accuracies_true <- accuracies_true[-which(accuracies_true$Time == "t1"), ]

  accuracy_plot <- ggplot2::ggplot(data = accuracies_true, ggplot2::aes(x = Time, y = Accuracy, fill = Method)) +
    ggplot2::geom_boxplot(position = ggplot2::position_dodge(), lwd = 0.25) +
    ggplot2::ylim(-0.5, 1) +
    ggplot2::theme_classic() +
    ggplot2::ylab("Prediction accuracy") +
    ggplot2::scale_fill_manual(values = c('#709BFF', "#17DE6D")) + # '#00120B',
    ggplot2::theme(legend.position = c(1-0.075, 0.15), #legend.position = c(0.9, 0.15),
          legend.title = ggplot2::element_text(size = 7),
          legend.text = ggplot2::element_text(size = 6),
          axis.title = ggplot2::element_text(size = 7, face = "plain"),
          axis.text = ggplot2::element_text(size = 6),
          legend.box.background = ggplot2::element_rect(colour = "black")
         )

  return(list(accuracies_true, accuracy_plot))
}
