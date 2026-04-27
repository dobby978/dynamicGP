# Mean-trajectory plots with standard-deviation ribbons.

# Suppress R CMD CHECK NOTEs for ggplot2 aes() column references.
utils::globalVariables(c("time_step", "trait", "value", "genotype",
                          "mean_val", "sd_val", "source", "geno_src"))

#' Plot Trait Trajectories Over Time
#'
#' Creates mean-trajectory plots (line \eqn{\pm} 1 SD ribbon) for one or more
#' traits across genotypes.  Accepts a single 3-D array **or** a named list
#' of 3-D arrays (e.g. observed and predicted), which are then overlaid on
#' the same axes with a colour legend.
#'
#' @param X_list A 3-D numeric array \eqn{[p \times N \times G]}, a
#'   \eqn{[p \times N]} matrix (single genotype), or a **named list** of such
#'   arrays.  Unnamed list elements are labelled `"Source 1"`, `"Source 2"`,
#'   etc.
#' @param traits Optional character vector of trait names to plot.  Must match
#'   `dimnames(X_list)[[1]]` (or `dimnames(X_list[[1]])[[1]]` for a list).
#'   `NULL` (default) plots all traits.
#' @param time_index Optional numeric vector of length `N` used as the x-axis
#'   values.  Defaults to integer indices `1, 2, ..., N`.
#' @param show_individuals Logical.  Draw a thin, semi-transparent line for
#'   each individual genotype in addition to the mean.  Default `FALSE`.
#' @param facet_by_trait Logical.  If `TRUE` and the number of selected traits
#'   is fewer than 30, produce a single faceted plot rather than one plot per
#'   trait.  Default `FALSE`.
#' @param save_plots Logical.  Save all plots to `plots_dir`.  Default `FALSE`.
#' @param plots_dir Character.  Output directory (created recursively if
#'   absent).  Default `"."`.
#' @param plot_format Character file format: `"png"`, `"svg"`, or `"pdf"`.
#'   Default `"png"`.
#' @param width,height Plot dimensions in inches.  Defaults `10` and `6`.
#'
#' @return Invisibly returns a named list of \pkg{ggplot2} objects.  Names
#'   are trait names, or `"faceted"` for a faceted plot.
#'
#' @examples
#' \dontrun{
#' # Overlay observed and reconstructed trajectories
#' result <- run_dmd(X_array[, -ncol(X_array)], X_array[, -1L])
#' pred   <- reconstruct_dmd(result$A_arr, X_obs = result$X_arr)
#' plot_trajectories(
#'   list(Observed = result$X_arr, Predicted = pred$X_pred_arr),
#'   show_individuals = TRUE,
#'   facet_by_trait   = TRUE
#' )
#' }
#'
#' @seealso [reconstruct_dmd()], [compare_arrays()]
#' @export
plot_trajectories <- function(
    X_list,
    traits           = NULL,
    time_index       = NULL,
    show_individuals = FALSE,
    facet_by_trait   = FALSE,
    save_plots       = FALSE,
    plots_dir        = ".",
    plot_format      = "png",
    width            = 10,
    height           = 6
) {
  # ── Normalise to a named list of 3-D arrays ───────────────────────────────
  if (!is.list(X_list) || is.array(X_list)) {
    X_list <- list(Data = .to_3d_arr(X_list))
  } else {
    X_list <- lapply(X_list, .to_3d_arr)
    nms    <- names(X_list)
    if (is.null(nms) || any(nms == ""))
      names(X_list) <- paste0("Source ", seq_along(X_list))
  }
  n_src     <- length(X_list)
  src_names <- names(X_list)

  # ── Metadata from first source ────────────────────────────────────────────
  arr1 <- X_list[[1L]]
  p    <- dim(arr1)[1L]
  N    <- dim(arr1)[2L]
  G    <- dim(arr1)[3L]

  all_traits <- dimnames(arr1)[[1L]]
  if (is.null(all_traits)) all_traits <- paste0("T", seq_len(p))
  g_names    <- dimnames(arr1)[[3L]]
  if (is.null(g_names))    g_names    <- paste0("G", seq_len(G))

  # ── Subset traits ─────────────────────────────────────────────────────────
  if (!is.null(traits)) {
    bad <- setdiff(traits, all_traits)
    if (length(bad))
      warning(sprintf("plot_trajectories: trait(s) not found and skipped: %s",
                      paste(bad, collapse = ", ")))
    trait_sel <- intersect(traits, all_traits)
  } else {
    trait_sel <- all_traits
  }
  if (length(trait_sel) == 0L)
    stop("plot_trajectories: no valid traits selected.")

  trait_idx <- match(trait_sel, all_traits)
  p_sel     <- length(trait_sel)

  # ── Time axis ─────────────────────────────────────────────────────────────
  if (is.null(time_index)) {
    # Try to extract from dimnames of the first array
    dn2 <- dimnames(arr1)[[2L]]
    if (!is.null(dn2)) {
      time_index <- tryCatch(as.numeric(dn2), warning = function(w) NULL)
    } else {
      time_index <- NULL
    }
    # Fall back to integer sequence if extraction failed
    if (is.null(time_index)) time_index <- seq_len(N)
  }
  if (length(time_index) != N)
    stop(sprintf("plot_trajectories: length(time_index) = %d but N = %d.",
                 length(time_index), N))
  # Convert time_index to character labels for discrete x-axis
  time_labels <- as.character(time_index)

  # ── Build individual-level long data frame ────────────────────────────────
  ind_rows <- list()
  for (si in seq_len(n_src)) {
    arr <- X_list[[si]]
    N_s <- min(dim(arr)[2L], N)
    for (ii in seq_along(trait_sel)) {
      i  <- trait_idx[ii]
      tr <- trait_sel[ii]
      for (gi in seq_len(G)) {
        ind_rows[[length(ind_rows) + 1L]] <- data.frame(
          source    = src_names[si],
          trait     = tr,
          genotype  = g_names[gi],
          time_step = as.numeric(time_labels[seq_len(N_s)]), #factor(time_labels[seq_len(N_s)], levels = time_labels),
          value     = as.numeric(arr[i, seq_len(N_s), gi]),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  df_long          <- do.call(rbind, ind_rows)
  df_long$source   <- factor(df_long$source,   levels = src_names)
  df_long$trait    <- factor(df_long$trait,    levels = trait_sel)
  df_long$geno_src <- paste(df_long$genotype, df_long$source, sep = "|")

  # ── Summary: mean ± SD per source × trait × time_step ────────────────────
  sum_rows <- list()
  for (si in seq_len(n_src)) {
    arr <- X_list[[si]]
    N_s <- min(dim(arr)[2L], N)
    for (ii in seq_along(trait_sel)) {
      i  <- trait_idx[ii]
      tr <- trait_sel[ii]
      for (t_idx in seq_len(N_s)) {
        vals <- as.numeric(arr[i, t_idx, ])
        fin  <- vals[is.finite(vals)]
        sum_rows[[length(sum_rows) + 1L]] <- data.frame(
          source    = src_names[si],
          trait     = tr,
          time_step = as.numeric(time_labels[t_idx]), #factor(time_labels[t_idx], levels = time_labels),
          mean_val  = if (length(fin) >= 1L) mean(fin) else NA_real_,
          sd_val    = if (length(fin) >= 2L) stats::sd(fin) else 0,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  df_sum         <- do.call(rbind, sum_rows)
  df_sum$source  <- factor(df_sum$source, levels = src_names)
  df_sum$trait   <- factor(df_sum$trait,  levels = trait_sel)

  # ── Colour palette ────────────────────────────────────────────────────────
  src_colours <- stats::setNames(
    scales::hue_pal()(n_src), src_names
  )

  # ── Per-trait plot builder ────────────────────────────────────────────────
  .make_plot <- function(tr) {
    ds <- df_sum[df_sum$trait == tr, , drop = FALSE]
    di <- df_long[df_long$trait == tr, , drop = FALSE]

    plt <- ggplot2::ggplot(
      ds,
      ggplot2::aes(x = time_step, y = mean_val,
                   colour = source, fill = source) #, group = source
    ) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
        alpha = 0.15, colour = NA
      ) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::geom_point(stat = "identity") 

    if (show_individuals) {
      plt <- plt +
        ggplot2::geom_line(
          data = di,
          ggplot2::aes(x = time_step, y = value,
                       group = geno_src, colour = source),
          linewidth = 0.25, alpha = 0.35, inherit.aes = FALSE
        )
    }
    plt +
      ggplot2::scale_colour_manual(values = src_colours, name = NULL) +
      ggplot2::scale_fill_manual(  values = src_colours, guide = "none") +
      ggplot2::scale_x_continuous() +
      ggplot2::labs(
        title = as.character(tr),
        x     = "Time",
        y     = "Value"
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        legend.position  = if (n_src > 1L) "bottom" else "none"
      )
  }

  # ── Save helper ───────────────────────────────────────────────────────────
  if (save_plots && !dir.exists(plots_dir))
    dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

  .save <- function(plt, name) {
    fpath <- file.path(plots_dir,
                       sprintf("%s.%s", make.names(name), plot_format))
    ggplot2::ggsave(fpath, plot = plt, device = plot_format,
                    width = width, height = height, units = "in")
    invisible(plt)
  }

  # ── Produce and collect plots ─────────────────────────────────────────────
  out_plots <- list()

  if (facet_by_trait && p_sel < 30L) {
    # Single faceted plot
    plt_fac <- ggplot2::ggplot(
      df_sum,
      ggplot2::aes(x = time_step, y = mean_val,
                   colour = source, fill = source)
    ) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
        alpha = 0.15, colour = NA
      ) +
      ggplot2::geom_line(linewidth = 0.7)

    if (show_individuals) {
      plt_fac <- plt_fac +
        ggplot2::geom_line(
          data = df_long,
          ggplot2::aes(x = time_step, y = value,
                       group = geno_src, colour = source),
          linewidth = 0.2, alpha = 0.3, inherit.aes = FALSE
        )
    }
    plt_fac <- plt_fac +
      ggplot2::scale_colour_manual(values = src_colours, name = NULL) +
      ggplot2::scale_fill_manual(  values = src_colours, guide = "none") +
      ggplot2::scale_x_continuous() +
      ggplot2::facet_wrap(~ trait, scales = "free_y") +
      ggplot2::labs(title = "Trait trajectories", x = "Time", y = "Value") +
      ggplot2::theme_bw(base_size = 10) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        strip.background = ggplot2::element_rect(fill = "grey92"),
        legend.position  = if (n_src > 1L) "bottom" else "none"
      )

    out_plots[["faceted"]] <- plt_fac
    if (save_plots) .save(plt_fac, "trajectories_faceted")

  } else {
    for (tr in trait_sel) {
      plt              <- .make_plot(tr)
      out_plots[[tr]]  <- plt
      if (save_plots) .save(plt, paste0("trajectory_", tr))
    }
  }

  invisible(out_plots)
}
