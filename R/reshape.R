# Reshape helpers: 3-D arrays <-> 2-D genotype-by-entry tables.

#' Flatten a 3-D Array to a 2-D Genotype-by-Entry Data Frame
#'
#' Converts any 3-D array \eqn{[d_1 \times d_2 \times G]} — such as `A_arr`,
#' `B_arr`, `X_arr` — into a 2-D
#' data frame with one row per genotype and one column per
#' \eqn{(d_1, d_2)} entry pair.
#'
#' Column names are formed by pasting the \eqn{d_1} and \eqn{d_2} dimnames
#' together with `sep`.  When dimnames are absent, fallback labels
#' `"d1_<i>"` / `"d2_<j>"` are used.  Row names are the genotype names from
#' `dimnames(arr)[[3]]`.
#'
#' The vectorisation order matches R's column-major convention: for a
#' \eqn{d_1 \times d_2} slice, element \eqn{(i, j)} appears at position
#' \eqn{(j - 1) \cdot d_1 + i}.
#'
#' @param arr A 3-D numeric array \eqn{[d_1 \times d_2 \times G]}, or a 2-D
#'   matrix (treated as a single genotype and returned as a 1-row data frame).
#' @param sep Character separator between the two dimension labels in the
#'   column names.  Default `"_._._"`.
#'
#' @return A `data.frame` with `G` rows (genotypes) and \eqn{d_1 \cdot d_2}
#'   columns.  Row names equal the genotype names.
#'
#' @examples
#' # Flatten A_arr [4 x 4 x 3] -> data frame [3 x 16]
#' A <- array(rnorm(4 * 4 * 3), dim = c(4, 4, 3),
#'            dimnames = list(paste0("T", 1:4), paste0("T", 1:4),
#'                            paste0("G", 1:3)))
#' df <- arr_to_table(A)
#' dim(df)   # 3 x 16
#'
#' @seealso [table_to_arr()]
#' @export
arr_to_table <- function(arr, sep = "_._._") {
  arr <- .to_3d_arr(arr)
  stopifnot(is.array(arr), length(dim(arr)) == 3L)

  d  <- dim(arr)
  dn <- dimnames(arr)
  G  <- d[3L]

  n1 <- if (!is.null(dn[[1L]])) dn[[1L]] else paste0("d1_", seq_len(d[1L]))
  n2 <- if (!is.null(dn[[2L]])) dn[[2L]] else paste0("d2_", seq_len(d[2L]))
  col_nms <- as.vector(outer(n1, n2, paste, sep = sep))

  g_nms <- if (!is.null(dn[[3L]])) dn[[3L]] else paste0("G", seq_len(G))

  mat <- do.call(rbind, lapply(seq_len(G), function(gi) as.vector(arr[, , gi])))
  rownames(mat) <- g_nms
  colnames(mat) <- col_nms
  as.data.frame(mat, check.names = FALSE)
}


#' Reconstruct a 3-D Array from a 2-D Genotype-by-Entry Table
#'
#' Reverses [arr_to_table()]: given a data frame with one row per genotype
#' and \eqn{d_1 \cdot d_2} columns, returns a \eqn{[d_1 \times d_2 \times G]}
#' array.
#'
#' @param tbl A data frame or matrix with `G` rows and \eqn{d_1 \cdot d_2}
#'   columns.  Row names are used as genotype names when available.
#' @param d1 Integer: size of the first array dimension.  When `NULL`
#'   (default) it is inferred from the number of unique first-part labels
#'   obtained by splitting `colnames(tbl)` at `sep`.
#' @param d2 Integer: size of the second array dimension.  When `NULL`
#'   (default) it is inferred analogously from the unique second-part labels.
#' @param sep Character separator used to split column names into their
#'   \eqn{(d_1, d_2)} label components.  Must match the `sep` used when
#'   [arr_to_table()] produced the table.  Default `"_._._"`.
#' @param dim1_names Optional character vector of length `d1` for
#'   `dimnames(arr)[[1]]`.  When `NULL` and `sep`-splitting succeeds, the
#'   unique first-part labels are used automatically.
#' @param dim2_names Optional character vector of length `d2` for
#'   `dimnames(arr)[[2]]`.  When `NULL` and `sep`-splitting succeeds, the
#'   unique second-part labels are used automatically.
#' @param geno_names Optional character vector of length `G` for
#'   `dimnames(arr)[[3]]`.  Defaults to `rownames(tbl)`.
#'
#' @return A 3-D numeric array \eqn{[d_1 \times d_2 \times G]}.
#'
#' @examples
#' A <- array(rnorm(4 * 4 * 3), dim = c(4, 4, 3),
#'            dimnames = list(paste0("T", 1:4), paste0("T", 1:4),
#'                            paste0("G", 1:3)))
#' df <- arr_to_table(A)
#' # d1 / d2 inferred automatically from column names:
#' A2 <- table_to_arr(df)
#' all.equal(A, A2)  # TRUE
#'
#' @seealso [arr_to_table()]
#' @export
table_to_arr <- function(tbl, d1 = NULL, d2 = NULL,
                          sep        = "_._._",
                          dim1_names = NULL,
                          dim2_names = NULL,
                          geno_names = NULL) {
  mat  <- as.matrix(tbl)
  G    <- nrow(mat)
  cns  <- colnames(mat)

  # ── Auto-infer d1 / d2 from column names when not supplied ────────────────
  inferred <- FALSE
  if (is.null(d1) || is.null(d2)) {
    if (is.null(cns))
      stop("table_to_arr: d1/d2 not supplied and tbl has no column names to infer from.")
    parts <- strsplit(cns, sep, fixed = TRUE)
    n_parts <- lengths(parts)
    if (any(n_parts < 2L))
      stop(sprintf(
        "table_to_arr: could not split all column names using sep = %s; supply d1 and d2 explicitly.",
        dQuote(sep)))
    lhs <- vapply(parts, `[[`, character(1L), 1L)
    rhs <- vapply(parts, function(p) paste(p[-1L], collapse = sep), character(1L))

    ul1 <- unique(lhs)
    ul2 <- unique(rhs)

    if (is.null(d1)) d1 <- length(ul1)
    if (is.null(d2)) d2 <- length(ul2)

    if (is.null(dim1_names)) dim1_names <- ul1
    if (is.null(dim2_names)) dim2_names <- ul2
    inferred <- TRUE
  }

  d1 <- as.integer(d1)
  d2 <- as.integer(d2)
  if (ncol(mat) != d1 * d2)
    stop(sprintf(
      "table_to_arr: tbl has %d column(s) but d1 * d2 = %d * %d = %d.%s",
      ncol(mat), d1, d2, d1 * d2,
      if (inferred) " Check that sep matches the separator used in arr_to_table()." else ""))

  g_nms <- if (!is.null(geno_names)) {
    geno_names
  } else if (!is.null(rownames(mat))) {
    rownames(mat)
  } else {
    paste0("G", seq_len(G))
  }

  arr <- array(NA_real_,
               dim      = c(d1, d2, G),
               dimnames = list(dim1_names, dim2_names, g_nms))
  for (gi in seq_len(G))
    arr[, , gi] <- matrix(mat[gi, ], nrow = d1, ncol = d2)
  arr
}
