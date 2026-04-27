# Trajectory property analysis for 3-D phenotype arrays.

# Suppress R CMD CHECK NOTEs for ggplot2 aes() column references.
utils::globalVariables(c("metric", "value", ".array", ".data"))

#' @importFrom stats lm coef var
NULL

# ── Internal metric computation ──────────────────────────────────────────────

# Compute trajectory properties for a single numeric time-series vector x.
# Returns a named list of scalar metrics. NA is returned for metrics that
# require more observations than are available after removing non-finite values.
.traj_metrics <- function(x) {
  fin <- which(is.finite(x))
  xf  <- x[fin]
  n   <- length(xf)

  na_metrics <- list(
    total_variation     = NA_real_,
    quadratic_variation = NA_real_,
    cv                  = NA_real_,
    rms_roughness       = NA_real_,
    kurtosis            = NA_real_,
    mean_second_diff    = NA_real_,
    quadratic_coef      = NA_real_,
    mean_chord_dev      = NA_real_,
    trend_slope         = NA_real_,
    lag1_acf            = NA_real_,
    skewness            = NA_real_,
    range               = NA_real_,
    n_obs               = 0L
  )

  if (n < 2L) {
    na_metrics$n_obs <- n
    return(na_metrics)
  }

  diffs1 <- diff(xf)
  mu     <- mean(xf)
  s      <- sd(xf)

  # a. Total variation: sum of absolute first differences
  tv <- sum(abs(diffs1))

  # b. Quadratic variation: variance of first differences
  #    (measures irregularity / variability of the increments)
  qv <- if (n >= 3L) var(diffs1) else NA_real_

  # c. Coefficient of variation: SD / |mean|
  cv_val <- if (mu != 0) s / abs(mu) else NA_real_

  # d. RMS roughness: root mean square of first differences
  #    (un-centred; differs from quadratic_variation when there is a trend)
  rms_r <- sqrt(mean(diffs1^2))

  # e. Kurtosis: 4th standardised central moment (Gaussian = 3)
  kurt <- if (n >= 4L && s > 0) mean((xf - mu)^4) / s^4 else NA_real_

  # f. Mean second difference: mean of the second-order differences
  #    (captures average acceleration / curvature)
  msd <- if (n >= 3L) mean(diff(diffs1)) else NA_real_

  # g. Quadratic coefficient: coefficient of t^2 in OLS fit x ~ t + t^2
  #    Uses the original (possibly irregular) time indices from `fin`.
  t_idx  <- fin
  quad_c <- if (n >= 3L) {
    tryCatch({
      fit <- lm(xf ~ t_idx + I(t_idx^2))
      unname(coef(fit)[3L])
    }, error   = function(e) NA_real_,
       warning = function(w) NA_real_)
  } else NA_real_

  # h. Mean chord deviation: mean |x[t] - chord(t)|
  #    where the chord connects (fin[1], xf[1]) to (fin[n], xf[n]).
  mcd <- if (fin[n] != fin[1L]) {
    chord <- xf[1L] + (xf[n] - xf[1L]) * (fin - fin[1L]) / (fin[n] - fin[1L])
    mean(abs(xf - chord))
  } else {
    NA_real_
  }

  # i. Trend slope: slope from OLS fit x ~ t  (linear trend per unit time)
  ts_slope <- tryCatch({
    fit <- lm(xf ~ t_idx)
    unname(coef(fit)[2L])
  }, error   = function(e) NA_real_,
     warning = function(w) NA_real_)

  # j. Lag-1 autocorrelation: cor(x[1..n-1], x[2..n])
  lag1 <- if (n >= 3L) {
    v <- cor(xf[-n], xf[-1L])
    if (is.finite(v)) v else NA_real_
  } else NA_real_

  # k. Skewness: 3rd standardised central moment
  skew <- if (n >= 3L && s > 0) mean((xf - mu)^3) / s^3 else NA_real_

  # l. Range: max - min
  rng <- max(xf) - min(xf)

  list(
    total_variation     = tv,
    quadratic_variation = qv,
    cv                  = cv_val,
    rms_roughness       = rms_r,
    kurtosis            = kurt,
    mean_second_diff    = msd,
    quadratic_coef      = quad_c,
    mean_chord_dev      = mcd,
    trend_slope         = ts_slope,
    lag1_acf            = lag1,
    skewness            = skew,
    range               = rng,
    n_obs               = n
  )
}

# ── Main exported function ───────────────────────────────────────────────────

#' Compute Trajectory Properties for 3-D Phenotype Arrays
#'
#' For every combination of array × trait × genotype, extracts the univariate
#' time series and computes a suite of trajectory characteristics.  The result
#' is a tidy data frame — one row per (array, trait, genotype) — that can be
#' used directly for downstream modelling or visualised with the optional
#' built-in plots.
#'
#' @section Metrics computed:
#' \describe{
#'   \item{`total_variation`}{Sum of absolute first differences
#'     \eqn{\sum_t |x_{t+1} - x_t|} — total distance travelled.}
#'   \item{`quadratic_variation`}{Variance of first differences
#'     \eqn{\mathrm{Var}(\Delta x)} — irregularity of increments.}
#'   \item{`cv`}{Coefficient of variation \eqn{\sigma / |\mu|} — relative
#'     spread.  `NA` when the trajectory mean is zero.}
#'   \item{`rms_roughness`}{\eqn{\sqrt{\mathrm{mean}(\Delta x^2)}} —
#'     un-centred RMS of increments; measures typical step magnitude.}
#'   \item{`kurtosis`}{4th standardised central moment
#'     \eqn{\mathrm{E}[(x - \mu)^4] / \sigma^4} (Gaussian = 3).}
#'   \item{`mean_second_diff`}{Mean of the second-order differences
#'     \eqn{\mathrm{mean}(\Delta^2 x)} — average acceleration/curvature.}
#'   \item{`quadratic_coef`}{Coefficient of \eqn{t^2} from OLS fit
#'     \eqn{x \sim 1 + t + t^2} — non-linearity of the trend.}
#'   \item{`mean_chord_dev`}{Mean absolute deviation from the straight-line
#'     chord connecting the first and last observed time points.}
#'   \item{`trend_slope`}{Slope from OLS fit \eqn{x \sim 1 + t} — linear
#'     trend per unit time step.}
#'   \item{`lag1_acf`}{Lag-1 autocorrelation \eqn{\rho(1) = \mathrm{cor}(x_t,
#'     x_{t+1})} — temporal persistence.}
#'   \item{`skewness`}{3rd standardised central moment
#'     \eqn{\mathrm{E}[(x-\mu)^3] / \sigma^3} — asymmetry.}
#'   \item{`range`}{\eqn{\max(x) - \min(x)} — peak-to-peak amplitude.}
#'   \item{`n_obs`}{Number of finite (non-`NA`) observations used.}
#' }
#'
#' @param X_list A 3-D numeric array \eqn{[p \times N \times G]} **or** a
#'   named list of such arrays (one per dataset / environment / period).  A
#'   single unnamed array is wrapped into a one-element list labelled
#'   `"arr1"`.  Each array element may also be a 2-D matrix, promoted to
#'   \eqn{[p \times N \times 1]} internally.
#' @param plot Logical.  When `TRUE`, a ggplot2 figure is produced and
#'   returned as an attribute (accessible via `attr(result, "plot")`).
#'   Default `FALSE`.
#' @param plot_type Character; one of `"boxplot"` (default), `"violin"`, or
#'   `"histogram"`.  Ignored when `plot = FALSE`.
#' @param metrics Character vector naming which metrics to include in the plot.
#'   Default `NULL` uses all metrics except `n_obs`.
#' @param facet_arrays Logical.  When `TRUE` **and** `X_list` has more than
#'   one element, the plot is faceted by both metric and array.  When `FALSE`
#'   (default) arrays are distinguished by fill colour.  Ignored when
#'   `plot = FALSE`.
#' @param save_plot Logical.  When `TRUE` the plot is saved to `plots_dir`.
#'   Default `FALSE`.
#' @param plots_dir Character path to the output directory for saved plots.
#'   Required when `save_plot = TRUE`.
#' @param width,height Plot dimensions in inches (default 10 × 8).
#'
#' @return A `data.frame` with columns `array`, `trait`, `genotype`, and one
#'   column per metric listed above.  When `plot = TRUE` the plot object is
#'   attached as `attr(result, "plot")`.
#'
#' @examples
#' set.seed(42)
#' p <- 3L; N <- 10L; G <- 5L
#' X <- array(rnorm(p * N * G), dim = c(p, N, G),
#'            dimnames = list(paste0("T", 1:p),
#'                            paste0("t", 1:N),
#'                            paste0("G", 1:G)))
#' df <- trajectory_properties(X)
#' head(df)
#'
#' # Multiple arrays (e.g. stress vs. stable environment)
#' X2 <- array(rnorm(p * N * G), dim = c(p, N, G),
#'             dimnames = dimnames(X))
#' df2 <- trajectory_properties(list(stress = X, stable = X2), plot = TRUE)
#' attr(df2, "plot")
#'
#' @seealso [concat_arrays()]
#' @export
trajectory_properties <- function(
    X_list,
    plot         = FALSE,
    plot_type    = c("boxplot", "violin", "histogram"),
    metrics      = NULL,
    facet_arrays = FALSE,
    save_plot    = FALSE,
    plots_dir    = NULL,
    width        = 10,
    height       = 8
) {
  plot_type <- match.arg(plot_type)

  # ── Normalise input to a list of 3-D arrays ──────────────────────────────
  if (is.array(X_list) || is.matrix(X_list)) {
    X_list <- list(arr1 = X_list)
  }
  if (!is.list(X_list))
    stop("trajectory_properties: X_list must be a 3-D array or a list of 3-D arrays.")

  X_list    <- lapply(X_list, .to_3d_arr)
  arr_names <- if (!is.null(names(X_list)) && !all(nzchar(names(X_list)) == 0L))
    names(X_list)
  else
    paste0("arr", seq_along(X_list))
  names(X_list) <- arr_names

  # ── Compute metrics for each (array, trait, genotype) combination ─────────
  all_rows <- lapply(seq_along(X_list), function(i) {
    a       <- X_list[[i]]
    p_dim   <- dim(a)[1L]
    G_dim   <- dim(a)[3L]
    t_nms   <- if (!is.null(dimnames(a)[[1L]])) dimnames(a)[[1L]]
               else paste0("T", seq_len(p_dim))
    g_nms   <- if (!is.null(dimnames(a)[[3L]])) dimnames(a)[[3L]]
               else paste0("G", seq_len(G_dim))

    rows_ij <- vector("list", p_dim * G_dim)
    k <- 0L
    for (j in seq_len(p_dim)) {
      for (gi in seq_len(G_dim)) {
        k       <- k + 1L
        x       <- a[j, , gi]
        m       <- .traj_metrics(x)
        rows_ij[[k]] <- data.frame(
          .array   = arr_names[i],
          trait    = t_nms[j],
          genotype = g_nms[gi],
          as.data.frame(m, stringsAsFactors = FALSE),
          stringsAsFactors = FALSE
        )
      }
    }
    do.call(rbind, rows_ij)
  })

  out           <- do.call(rbind, all_rows)
  rownames(out) <- NULL
  names(out)[names(out) == ".array"] <- "array"

  # ── Optional visualisation ────────────────────────────────────────────────
  if (plot) {
    all_metric_names <- c(
      "total_variation", "quadratic_variation", "cv", "rms_roughness",
      "kurtosis", "mean_second_diff", "quadratic_coef", "mean_chord_dev",
      "trend_slope", "lag1_acf", "skewness", "range"
    )
    plot_metrics <- if (!is.null(metrics)) metrics else all_metric_names
    # Only keep metrics that exist as columns
    plot_metrics <- intersect(plot_metrics, names(out))

    if (length(plot_metrics) == 0L) {
      warning("trajectory_properties: no valid metric columns found for plotting.")
    } else {
      # Reshape to long format
      df_long <- tidyr::pivot_longer(
        out[, c("array", "trait", "genotype", plot_metrics)],
        cols      = tidyr::all_of(plot_metrics),
        names_to  = "metric",
        values_to = "value"
      )
      df_long$metric <- factor(df_long$metric, levels = plot_metrics)

      multi_arr <- length(unique(out$array)) > 1L

      p_obj <- if (plot_type == "boxplot") {
        ggplot2::ggplot(df_long,
          if (multi_arr)
            ggplot2::aes(x = .data[["array"]], y = value, fill = .data[["array"]])
          else
            ggplot2::aes(x = trait, y = value, fill = trait)) +
          ggplot2::geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.4) +
          ggplot2::facet_wrap(
            if (multi_arr && facet_arrays)
              ggplot2::vars(.data[["metric"]], .data[["array"]])
            else
              ggplot2::vars(.data[["metric"]]),
            scales = "free_y"
          )

      } else if (plot_type == "violin") {
        ggplot2::ggplot(df_long,
          if (multi_arr)
            ggplot2::aes(x = .data[["array"]], y = value, fill = .data[["array"]])
          else
            ggplot2::aes(x = trait, y = value, fill = trait)) +
          ggplot2::geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), trim = TRUE) +
          ggplot2::facet_wrap(ggplot2::vars(.data[["metric"]]), scales = "free_y")

      } else {  # histogram
        ggplot2::ggplot(df_long,
          if (multi_arr)
            ggplot2::aes(x = value, fill = .data[["array"]])
          else
            ggplot2::aes(x = value, fill = trait)) +
          ggplot2::geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
          ggplot2::facet_wrap(ggplot2::vars(.data[["metric"]]), scales = "free")
      }

      p_obj <- p_obj +
        ggplot2::labs(
          title = "Trajectory properties",
          x     = NULL, y = NULL, fill = NULL
        ) +
        ggplot2::theme_bw(base_size = 10) +
        ggplot2::theme(
          strip.background = ggplot2::element_rect(fill = "grey92"),
          legend.position  = "bottom"
        )

      if (save_plot) {
        if (is.null(plots_dir))
          stop("trajectory_properties: plots_dir must be specified when save_plot = TRUE.")
        if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
        fpath <- file.path(plots_dir, "trajectory_properties.png")
        ggplot2::ggsave(fpath, plot = p_obj, width = width, height = height,
                        dpi = 150)
        message("trajectory_properties: plot saved to ", fpath)
      }

      attr(out, "plot") <- p_obj
    }
  }

  out
}
