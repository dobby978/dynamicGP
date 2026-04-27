# Concatenation of 3-D snapshot arrays along the time dimension.

#' Concatenate Multiple 3-D Snapshot Arrays Along the Time Dimension
#'
#' Joins a named list of 3-D phenotype arrays \eqn{[p \times N_i \times G]}
#' into a single array \eqn{[p \times N_{total} \times G]} by stacking along
#' the time (second) dimension.  An optional parallel list of environmental /
#' control arrays is processed identically.  The returned `break_after` vector
#' marks the concatenation seams so that cross-array snapshot pairs can be
# excluded when fitting DMD models.
#'
#' All input arrays must share the **same** trait dimension \eqn{p} and
#' genotype dimension \eqn{G}, with matching `dimnames` on those axes (when
#' present).  Time-point labels from each element are prefixed with the
#' source-array name and `sep` to avoid duplicate label collisions.
#'
#' @param X_list Named list of 3-D numeric arrays, each of shape
#'   \eqn{[p \times N_i \times G]}.  A 2-D matrix is accepted as a
#'   single-genotype convenience form and promoted internally.  The names of
#'   the list become source-array labels embedded in the output time dimnames
#'   and in `array_labels`.
#' @param U_list Optional named list of control / environmental arrays,
#'   parallel to `X_list` (must have the same length).  Each element may be
#'   either a matrix \eqn{[q \times N_i]} (same controls for every genotype,
#'   replicated internally along the genotype dimension) or a 3-D array
#'   \eqn{[q \times N_i \times G]}.  All elements must share the same number
#'   of control factors \eqn{q}.  `NULL` (default) — no control inputs.
#' @param sep Character string inserted between the source-array label and the
#'   original time-point label when constructing the combined time dimnames.
#'   Default `":"`.
#'
#' @return A named list with four elements:
#' \describe{
#'   \item{`X_arr`}{3-D numeric array \eqn{[p \times N_{total} \times G]} of
#'     concatenated state snapshots, with `dimnames` set on all three axes.}
#'   \item{`U_arr`}{3-D numeric array \eqn{[q \times N_{total} \times G]} of
#'     concatenated control inputs, or `NULL` when `U_list` is `NULL`.}
#'   \item{`break_after`}{Integer vector of 1-based column indices in the
#'     snapshot-pair matrix \eqn{X_1} (which has \eqn{N_{total} - 1} columns)
#'     after which a concatenation seam lies.  For a two-array concatenation
#'     this is just \eqn{N_1}.  Pass directly to the `break_after` argument
#'     of [run_dmd()] to exclude cross-seam pairs from model fitting.}
#'   \item{`array_labels`}{Character vector of length \eqn{N_{total}} mapping
#'     each time point back to its source-array name.}
#' }
#'
#' @examples
#' p <- 4L; G <- 10L
#' set.seed(1)
#' X_env1 <- array(rnorm(p * 5L * G), dim = c(p, 5L, G),
#'                 dimnames = list(paste0("T", 1:p),
#'                                 paste0("t", 1:5),
#'                                 paste0("G", 1:G)))
#' X_env2 <- array(rnorm(p * 7L * G), dim = c(p, 7L, G),
#'                 dimnames = list(paste0("T", 1:p),
#'                                 paste0("t", 1:7),
#'                                 paste0("G", 1:G)))
#' res <- concat_arrays(list(stress = X_env1, stable = X_env2))
#' dim(res$X_arr)      # 4 x 12 x 10
#' res$break_after     # 5  (seam after stress time point 5)
#'
#' @seealso [trajectory_properties()]
#' @export
concat_arrays <- function(X_list, U_list = NULL, sep = ":") {

  # ── Input validation ────────────────────────────────────────────────────────
  if (!is.list(X_list))
    stop("concat_arrays: X_list must be a list of 3-D arrays.")
  if (length(X_list) < 1L)
    stop("concat_arrays: X_list must contain at least one array.")

  X_list <- lapply(X_list, .to_3d_arr)

  d_list <- lapply(X_list, dim)
  p      <- d_list[[1L]][1L]
  G      <- d_list[[1L]][3L]

  for (i in seq_along(X_list)) {
    if (length(d_list[[i]]) != 3L)
      stop(sprintf("concat_arrays: X_list[[%d]] is not a 3-D array.", i))
    if (d_list[[i]][1L] != p)
      stop(sprintf(
        "concat_arrays: X_list[[%d]] has %d trait(s); expected %d (from X_list[[1]]).",
        i, d_list[[i]][1L], p))
    if (d_list[[i]][3L] != G)
      stop(sprintf(
        "concat_arrays: X_list[[%d]] has %d genotype(s); expected %d (from X_list[[1]]).",
        i, d_list[[i]][3L], G))
  }

  # Check that trait and genotype dimnames are consistent across arrays
  t_names <- dimnames(X_list[[1L]])[[1L]]
  g_names <- dimnames(X_list[[1L]])[[3L]]

  for (i in seq_along(X_list)[-1L]) {
    tn_i <- dimnames(X_list[[i]])[[1L]]
    gn_i <- dimnames(X_list[[i]])[[3L]]
    if (!is.null(t_names) && !is.null(tn_i) && !identical(t_names, tn_i))
      stop(sprintf(
        "concat_arrays: X_list[[%d]] has different trait names from X_list[[1]].", i))
    if (!is.null(g_names) && !is.null(gn_i) && !identical(g_names, gn_i))
      stop(sprintf(
        "concat_arrays: X_list[[%d]] has different genotype names from X_list[[1]].", i))
  }

  # ── Source labels ────────────────────────────────────────────────────────────
  arr_labels <- if (!is.null(names(X_list)) && !all(nzchar(names(X_list)) == 0L))
    names(X_list)
  else
    paste0("arr", seq_along(X_list))

  # ── Time-point labels for the concatenated array ────────────────────────────
  N_vec <- vapply(X_list, function(a) dim(a)[2L], integer(1L))

  tp_labels <- unlist(lapply(seq_along(X_list), function(i) {
    tp_i <- dimnames(X_list[[i]])[[2L]]
    if (is.null(tp_i)) tp_i <- as.character(seq_len(N_vec[i]))
    paste0(arr_labels[i], sep, tp_i)
  }), use.names = FALSE)

  source_labels <- rep(arr_labels, times = N_vec)

  # ── Assemble X_arr ───────────────────────────────────────────────────────────
  N_total <- sum(N_vec)
  X_arr   <- array(NA_real_,
                   dim      = c(p, N_total, G),
                   dimnames = list(t_names, tp_labels, g_names))

  col_off <- 0L
  for (i in seq_along(X_list)) {
    idx            <- col_off + seq_len(N_vec[i])
    X_arr[, idx, ] <- X_list[[i]]
    col_off        <- col_off + N_vec[i]
  }

  # ── break_after: seam positions in X1 (length N_total - 1) ─────────────────
  # X1 = X[, -N_total], X2 = X[, -1].  The pair spanning seam i is at
  # column cumN[i] of X1 (time point cumN[i] paired with cumN[i]+1).
  cumN        <- cumsum(N_vec)
  break_after <- as.integer(cumN[-length(cumN)])

  # ── U_list processing ────────────────────────────────────────────────────────
  U_arr <- NULL

  if (!is.null(U_list)) {
    if (!is.list(U_list))
      stop("concat_arrays: U_list must be a list (or NULL).")
    if (length(U_list) != length(X_list))
      stop(sprintf(
        "concat_arrays: U_list has %d element(s) but X_list has %d.",
        length(U_list), length(X_list)))

    # Promote each U element to 3-D [q x N_i x G]
    U3d <- lapply(seq_along(U_list), function(i) {
      u   <- U_list[[i]]
      N_i <- N_vec[i]

      if (is.null(u))
        stop(sprintf("concat_arrays: U_list[[%d]] is NULL; provide an array or remove U_list.", i))

      if (is.matrix(u) || (is.array(u) && length(dim(u)) == 2L)) {
        if (ncol(u) != N_i)
          stop(sprintf(
            "concat_arrays: U_list[[%d]] has %d time column(s); X_list[[%d]] has %d.",
            i, ncol(u), i, N_i))
        # Replicate [q x N] -> [q x N x G]
        array(rep(as.vector(u), times = G),
              dim      = c(nrow(u), N_i, G),
              dimnames = list(rownames(u), NULL, g_names))

      } else if (is.array(u) && length(dim(u)) == 3L) {
        if (dim(u)[2L] != N_i)
          stop(sprintf(
            "concat_arrays: U_list[[%d]] has %d time column(s); X_list[[%d]] has %d.",
            i, dim(u)[2L], i, N_i))
        if (dim(u)[3L] != G)
          stop(sprintf(
            "concat_arrays: U_list[[%d]] has %d genotype slice(s); expected %d.",
            i, dim(u)[3L], G))
        u

      } else {
        stop(sprintf(
          "concat_arrays: U_list[[%d]] must be a matrix [q x N] or 3-D array [q x N x G].", i))
      }
    })

    q_vals <- vapply(U3d, function(u) dim(u)[1L], integer(1L))
    if (length(unique(q_vals)) > 1L)
      stop("concat_arrays: U_list elements have different numbers of control factors (q).")

    q       <- q_vals[1L]
    u_names <- dimnames(U3d[[1L]])[[1L]]

    U_arr <- array(NA_real_,
                   dim      = c(q, N_total, G),
                   dimnames = list(u_names, tp_labels, g_names))
    col_off <- 0L
    for (i in seq_along(U3d)) {
      idx             <- col_off + seq_len(N_vec[i])
      U_arr[, idx, ]  <- U3d[[i]]
      col_off         <- col_off + N_vec[i]
    }
  }

  list(
    X_arr        = X_arr,
    U_arr        = U_arr,
    break_after  = break_after,
    array_labels = source_labels
  )
}
