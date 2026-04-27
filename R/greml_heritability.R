# =============================================================================
# greml_heritability()  —  SNP-based heritability via GREML (REML estimation)

# Suppress R CMD CHECK NOTE for ggplot2 aes() column reference.
utils::globalVariables("h2")
#
# Two computation backends are supported:
#
#   backend = "r"    (default)
#     Uses the 'gaston' package to run Average Information REML (AI-REML)
#     entirely within R.  No external software is required.  Fast for typical
#     plant/animal breeding datasets (n < ~3 000 genotypes).
#
#   backend = "gcta"  (opt-in)
#     Shells out to the GCTA binary (gcta64 / gcta-1.94.1) for each trait.
#     Requires GCTA to be installed and the path provided via 'gcta_bin'.
#     Preferable only when n > ~5 000 genotypes and GCTA is already set up.
#
# Both backends return the same output data frame, so switching is seamless.
# =============================================================================

#' Compute GREML SNP-based Heritability
#'
#' Estimates narrow-sense SNP heritability (\eqn{h^2_{SNP}}) for each column
#' of `matrix_unrolled` using Genomic REML (GREML).  Two backends are
#' available: a pure-R implementation via the \pkg{gaston} package
#' (`backend = "r"`, default) and an optional GCTA binary backend
#' (`backend = "gcta"`).
#'
#' @section Model:
#' For each trait \eqn{y}, the univariate mixed model
#' \deqn{y = 1\mu + u + e, \quad u \sim \mathcal{N}(0, G\sigma^2_g),
#'        \quad e \sim \mathcal{N}(0, I\sigma^2_e)}
#' is fitted by REML.  Heritability is then
#' \eqn{h^2 = \sigma^2_g / (\sigma^2_g + \sigma^2_e)}.
#'
#' @section Backend choice:
#' \describe{
#'   \item{`"r"` (default)}{Runs AI-REML via \code{\link[gaston]{lmm.aireml}}.
#'     No external software needed.  Install with
#'     \code{install.packages("gaston")}.  Fastest for \eqn{n < 3\,000}.}
#'   \item{`"gcta"`}{Calls the GCTA binary (usually \code{gcta64}) via
#'     \code{system()}.  Requires GCTA to be installed on your system and the
#'     path supplied via \code{gcta_bin}.  Advantageous only when
#'     \eqn{n > 5\,000}.  See \url{https://yanglab.westlake.edu.cn/software/gcta/}.}
#' }
#'
#' @param matrix_unrolled A numeric data frame or matrix with \eqn{G} rows
#'   (genotypes) and \eqn{q} columns (traits).  Row names must be genotype
#'   identifiers matching the row/column names of `kinship_matrix`.  This is
#'   typically the output of [arr_to_table()].
#' @param kinship_matrix A square numeric matrix or data frame of dimension
#'   \eqn{G \times G} containing the genomic relationship matrix (GRM).  Row
#'   and column names should be the same genotype identifiers as
#'   `rownames(matrix_unrolled)`.  Required for `backend = "r"`.  For
#'   `backend = "gcta"` either this or `kinship_filename` must be provided.
#' @param kinship_filename Character.  Path prefix to a pre-existing GCTA
#'   binary GRM (\code{.grm.bin}, \code{.grm.N.bin}, \code{.grm.id}).  Only
#'   used when `backend = "gcta"`.  If `kinship_matrix` is provided instead,
#'   the binary files are written to `run_name` automatically using
#'   \pkg{plinkFile}.
#' @param backend Character; one of `"r"` (default) or `"gcta"`.  Controls
#'   which computation engine is used.  See the *Backend choice* section.
#' @param gcta_bin Character.  Name or full path of the GCTA executable.
#'   Only used when `backend = "gcta"`.  Defaults to `"gcta64"` (assumes GCTA
#'   is on the system \code{PATH}).  Set to the absolute path if needed, e.g.
#'   `gcta_bin = "/usr/local/bin/gcta64"`.
#' @param run_name Character.  Name of the temporary working directory created
#'   when `backend = "gcta"`.  Must be a plain directory name with no path
#'   separators (e.g. `"greml_temp"`, not `"../greml_temp"`).  Ignored for
#'   `backend = "r"`.  Default `"greml_temp"`.
#' @param save_run Logical.  When `TRUE` the temporary GCTA directory (and its
#'   `.hsq` files) is kept after the run.  Default `FALSE`.
#' @param save_filename Character or `NULL`.  If a file path is provided the
#'   results data frame is written there as a CSV.
#' @param plot Logical.  When `TRUE` a histogram of \eqn{\hat{h}^2} values is
#'   attached to the return value as `attr(result, "plot")`.  Default `TRUE`.
#' @param progress_bar Logical.  Show a trait-level progress bar.  Default
#'   `FALSE`.
#' @param verbose Logical.  When `TRUE`, GCTA console output is echoed and
#'   gaston convergence messages are shown.  Default `FALSE`.
#' @param tag Optional character label appended to progress bar and plot
#'   titles (e.g. an environment or time-period name).  Default `""`.
#'
#' @return A `data.frame` with one row per trait and columns:
#' \describe{
#'   \item{`trait`}{Trait name (from `colnames(matrix_unrolled)`).}
#'   \item{`h2`}{GREML estimate of \eqn{h^2_{SNP}}.}
#'   \item{`se_h2`}{Standard error of \eqn{\hat{h}^2}. For `backend = "r"`
#'     this is computed via the delta method from the AI matrix; for
#'     `backend = "gcta"` it is taken directly from the `.hsq` file.}
#'   \item{`V_G`}{Estimated genetic variance \eqn{\hat{\sigma}^2_g}.}
#'   \item{`V_e`}{Estimated residual variance \eqn{\hat{\sigma}^2_e}.}
#'   \item{`V_p`}{Total phenotypic variance \eqn{\hat{\sigma}^2_g +
#'     \hat{\sigma}^2_e}.}
#' }
#' When `plot = TRUE`, the ggplot2 histogram is attached as
#' `attr(result, "plot")`.
#'
#' @examples
#' \dontrun{
#' # ── R backend (default) ──────────────────────────────────────────────────
#' # 'unrolled' is a [G x q] matrix from arr_to_table()
#' # 'K' is a [G x G] genomic relationship matrix
#' h2_df <- greml_heritability(unrolled, kinship_matrix = K)
#' attr(h2_df, "plot")   # view the histogram
#'
#' # ── GCTA backend (opt-in) ────────────────────────────────────────────────
#' h2_df <- greml_heritability(
#'   unrolled,
#'   kinship_matrix = K,
#'   backend  = "gcta",
#'   gcta_bin = "/usr/local/bin/gcta64"
#' )
#' }
#'
#' @seealso [arr_to_table()]
#' @export
greml_heritability <- function(
    matrix_unrolled,
    kinship_matrix   = NULL,
    kinship_filename = NULL,
    backend          = c("r", "gcta"),
    gcta_bin         = "gcta64",
    run_name         = "greml_temp",
    save_run         = FALSE,
    save_filename    = NULL,
    plot             = TRUE,
    progress_bar     = FALSE,
    verbose          = FALSE,
    tag              = ""
) {

  # ── 1. Resolve and validate the backend choice ───────────────────────────
  backend <- match.arg(backend)

  # At least one kinship source must be provided regardless of backend.
  if (is.null(kinship_matrix) && is.null(kinship_filename)) {
    stop("greml_heritability: provide either 'kinship_matrix' (an [n x n] R ",
         "matrix) or 'kinship_filename' (a path prefix to GCTA binary GRM files).")
  }

  # Backend-specific pre-flight checks -------------------------------------------

  if (backend == "r") {
    # The 'gaston' package is required for the R backend but listed in Suggests
    # so that users who only want the GCTA backend need not install it.
    if (!requireNamespace("gaston", quietly = TRUE)) {
      stop("greml_heritability: backend = 'r' requires the 'gaston' package.\n",
           "  Install it with: install.packages('gaston')\n",
           "  Alternatively, use backend = 'gcta' if GCTA is installed.")
    }
    if (is.null(kinship_matrix)) {
      stop("greml_heritability: backend = 'r' requires 'kinship_matrix' to be ",
           "provided as an in-memory R matrix.  'kinship_filename' points to GCTA ",
           "binary files and cannot be read by the R backend.")
    }
  }

  if (backend == "gcta") {
    # The user must explicitly supply the GCTA binary path (or have it on PATH).
    gcta_found <- nchar(Sys.which(gcta_bin)) > 0L || file.exists(gcta_bin)
    if (!gcta_found) {
      stop("greml_heritability: GCTA binary not found.\n",
           "  Searched for: '", gcta_bin, "'\n",
           "  Either add GCTA to your system PATH, or pass the full path via ",
           "gcta_bin = '/path/to/gcta64'.\n",
           "  Download GCTA from: https://yanglab.westlake.edu.cn/software/gcta/")
    }
  }

  # ── 2. Validate genotype overlap between matrix_unrolled and kinship_matrix
  if (!is.null(kinship_matrix)) {
    pheno_genos <- rownames(matrix_unrolled)
    grm_genos   <- rownames(kinship_matrix)

    if (is.null(pheno_genos)) {
      stop("greml_heritability: 'matrix_unrolled' must have row names ",
           "(genotype identifiers).")
    }
    if (is.null(grm_genos)) {
      stop("greml_heritability: 'kinship_matrix' must have row and column names ",
           "(matching genotype identifiers from 'matrix_unrolled').")
    }

    # Find common genotypes
    common_genos <- intersect(pheno_genos, grm_genos)
    if (length(common_genos) == 0L) {
      stop("greml_heritability: no overlap between genotype IDs in ",
           "'matrix_unrolled' and 'kinship_matrix'.\n",
           "  matrix_unrolled has ", length(pheno_genos), " genotypes.\n",
           "  kinship_matrix has ", length(grm_genos), " genotypes.\n",
           "  Common genotypes: 0. Check that row/column names match.")
    }

    # Optionally warn about missing genotypes (will be handled by subsetting to common)
    missing_in_grm <- setdiff(pheno_genos, grm_genos)
    if (length(missing_in_grm) > 0L && verbose) {
      message("greml_heritability: ", length(missing_in_grm), " genotype(s) in ",
              "'matrix_unrolled' not found in 'kinship_matrix'. Proceeding with ",
              length(common_genos), " common genotype(s).")
    }

    # Optionally warn about extra genotypes in kinship matrix
    extra_in_grm <- setdiff(grm_genos, pheno_genos)
    if (length(extra_in_grm) > 0L && !all(is.na(extra_in_grm)) && verbose) {
      message("greml_heritability: ", length(extra_in_grm), " genotype(s) in ",
              "'kinship_matrix' not found in 'matrix_unrolled'. These will be ignored.")
    }

    # Subset both inputs to the common set of genotypes (preserving order from
    # matrix_unrolled, which is the phenotype table ordering).
    matrix_unrolled <- matrix_unrolled[common_genos, , drop = FALSE]
    kinship_matrix  <- kinship_matrix[common_genos, common_genos, drop = FALSE]
  }

  # ── 3. Coerce kinship_matrix to a plain numeric matrix ───────────────────
  K_mat <- if (!is.null(kinship_matrix)) as.matrix(kinship_matrix) else NULL

  # ── 4. Set up an optional progress bar ───────────────────────────────────
  n_traits <- ncol(matrix_unrolled)
  pb <- if (progress_bar) {
    progress::progress_bar$new(
      format = paste0(
        "  GREML", if (nzchar(tag)) paste0(" [", tag, "]") else "",
        " [:bar] :current/:total  ETA: :eta"
      ),
      total = n_traits, width = 100, clear = FALSE,
      force = TRUE   # always update in place; prevents re-printing on each tick
    )
  } else NULL

  # ── 5. Dispatch to the selected backend ──────────────────────────────────
  results <- if (backend == "r") {
    .greml_backend_r(matrix_unrolled, K_mat, pb, verbose)
  } else {
    .greml_backend_gcta(
      matrix_unrolled, K_mat, kinship_filename,
      gcta_bin, run_name, save_run, pb, verbose
    )
  }

  # ── 5. Attach trait names and tidy column order ──────────────────────────
  results$trait <- colnames(matrix_unrolled)
  results       <- results[, c("trait", "h2", "se_h2", "V_G", "V_e", "V_p")]
  rownames(results) <- NULL

  # ── 6. Optional CSV export ────────────────────────────────────────────────
  if (!is.null(save_filename) && is.character(save_filename)) {
    write.csv(results, save_filename, row.names = FALSE)
  }

  # ── 7. Optional heritability histogram ───────────────────────────────────
  if (plot) {
    p_obj <- ggplot2::ggplot(results, ggplot2::aes(x = h2)) +
      ggplot2::geom_histogram(binwidth = 0.05, fill = "steelblue",
                              colour = "white", na.rm = TRUE) +
      ggplot2::scale_x_continuous(limits = c(0, 1),
                                  name   = expression(hat(h)[SNP]^2)) +
      ggplot2::labs(
        title = paste0("SNP-based heritability",
                       if (nzchar(tag)) paste0("  \u2014  ", tag) else ""),
        y = "Number of traits"
      ) +
      ggplot2::theme_bw(base_size = 11)
    attr(results, "plot") <- p_obj
  }

  results
}


# =============================================================================
# Internal backend: R / gaston AI-REML
# =============================================================================
#
# For each trait column of matrix_unrolled, fits the model
#   y = mu*1 + u + e,   u ~ N(0, K * sigma_g^2),   e ~ N(0, I * sigma_e^2)
# using gaston::lmm.aireml().
#
# Standard errors for h^2 are obtained via the delta method applied to the
# variance-covariance matrix of the variance component estimates (the inverse
# of the Average Information matrix returned by gaston).
#
# Returns a data frame with columns: h2, se_h2, V_G, V_e, V_p.
# =============================================================================
.greml_backend_r <- function(matrix_unrolled, K_mat, pb, verbose) {

  n_traits <- ncol(matrix_unrolled)

  # Pre-align kinship matrix rows/cols to match the genotype order in the
  # phenotype table.  This guards against cases where the GRM was constructed
  # from a differently-ordered sample.
  geno_ids <- rownames(matrix_unrolled)
  if (!is.null(rownames(K_mat)) && !is.null(geno_ids)) {
    K_mat <- K_mat[geno_ids, geno_ids]
  }

  # Initialise output vectors (NA = REML failed or insufficient data)
  h2_vec <- se_vec <- Vg_vec <- Ve_vec <- Vp_vec <- rep(NA_real_, n_traits)

  for (i in seq_len(n_traits)) {

    y <- as.numeric(matrix_unrolled[, i])

    # Drop individuals with missing phenotype; adjust K accordingly.
    keep <- which(is.finite(y))
    y_k  <- y[keep]
    K_k  <- K_mat[keep, keep]

    # Fit the mixed model.  lmm.aireml() returns:
    #   $tau    — genetic variance component (sigma_g^2)
    #   $sigma2 — residual variance component (sigma_e^2)
    #   $P      — n x n precision matrix (only when get.P = TRUE)
    #   $Py     — n-vector P * y (always returned)
    #
    # gaston's C backend emits diagnostic text (including the non-fatal
    # "EM step failed to improve likelihood" notice) via Rprintf() directly
    # to stdout, which bypasses R's warning/message system.  We therefore
    # use three complementary suppression layers when verbose = FALSE:
    #   1. suppressWarnings()           — muffles R-level warnings
    #   2. suppressMessages()           — muffles R-level messages
    #   3. capture.output(type="output")— captures stdout from C/cat()
    # When verbose = TRUE everything is forwarded so the user can see it.
    fit <- local({
      fit_inner <- NULL
      tryCatch({
        if (!verbose) {
          suppressWarnings(suppressMessages(
            utils::capture.output(
              {fit_inner <- gaston::lmm.aireml(
                Y       = y_k,
                K       = list(K_k),
                verbose = FALSE,
                get.P   = TRUE
              )},
              type = "output"
            )
          ))
        } else {
          fit_inner <- gaston::lmm.aireml(
            Y       = y_k,
            K       = list(K_k),
            verbose = TRUE,
            get.P   = TRUE
          )
        }
        fit_inner
      }, error = function(e) {
        if (verbose) message("gaston failed: ", conditionMessage(e))
        NULL
      })
    })

    if (!is.null(fit)) {
      tau    <- fit$tau     # sigma_g^2 estimate
      sigma2 <- fit$sigma2  # sigma_e^2 estimate
      Vp     <- tau + sigma2

      # h^2 estimate
      h2_val <- tau / Vp

      # Delta-method SE for h^2 = tau / (tau + sigma2):
      #   grad = [ d(h2)/d(tau),   d(h2)/d(sigma2) ]
      #        = [ sigma2/Vp^2,   -tau/Vp^2 ]
      #   Var(h2) = grad^T %*% Cov([tau, sigma2]) %*% grad
      # Cov([tau, sigma2]) = solve(AI), where the 2x2 Average Information matrix is
      #   AI[i,j] = (1/2) y^T P V_i P V_j P y,  V_tau = K, V_sigma2 = I
      # gaston returns P via get.P = TRUE and Py as fit$Py.
      se_val <- tryCatch({
        Py_vec <- fit$Py            # n-vector: P * y
        KPy    <- K_k %*% Py_vec   # n-vector: K * P * y
        PKPy   <- fit$P %*% KPy    # n-vector: P * K * P * y
        PPy    <- fit$P %*% Py_vec # n-vector: P * P * y
        # AI matrix entries (each a scalar)
        AI_mat <- matrix(c(
          0.5 * as.numeric(crossprod(KPy,    PKPy)),   # AI[tau,   tau  ]
          0.5 * as.numeric(crossprod(KPy,    PPy)),    # AI[sigma2,tau  ] (symmetric)
          0.5 * as.numeric(crossprod(KPy,    PPy)),    # AI[tau,   sigma2]
          0.5 * as.numeric(crossprod(Py_vec, PPy))     # AI[sigma2,sigma2]
        ), nrow = 2)
        Cov_mat <- solve(AI_mat)   # 2x2 Cov([tau, sigma2])
        grad <- c(sigma2 / Vp^2, -tau / Vp^2)
        sqrt(as.numeric(t(grad) %*% Cov_mat %*% grad))
      }, error = function(e) NA_real_)

      h2_vec[i] <- h2_val
      se_vec[i] <- se_val
      Vg_vec[i] <- tau
      Ve_vec[i] <- sigma2
      Vp_vec[i] <- Vp
    }
    # If fit is NULL, the NA_real_ defaults from above are retained.

    if (!is.null(pb)) pb$tick()
  }

  data.frame(h2 = h2_vec, se_h2 = se_vec,
             V_G = Vg_vec, V_e = Ve_vec, V_p = Vp_vec,
             stringsAsFactors = FALSE)
}


# =============================================================================
# Internal backend: GCTA binary
# =============================================================================
#
# Writes a GCTA-format phenotype file once, then loops over traits calling
#   gcta64 --grm <prefix> --pheno <file> --mpheno <i> --reml --out <out>
# for each trait.
#
# The .hsq output is parsed by field name (not character position), making it
# robust to GCTA version differences.
#
# If kinship_matrix is provided (R matrix), it is converted to the GCTA
# binary GRM format using plinkFile::saveGRM().
#
# Returns a data frame with columns: h2, se_h2, V_G, V_e, V_p.
# =============================================================================
.greml_backend_gcta <- function(matrix_unrolled, K_mat, kinship_filename,
                                gcta_bin, run_name, save_run, pb, verbose) {

  # Sanitise run_name: strip any path separators so it can only name a
  # subdirectory of the current working directory, preventing path traversal.
  run_name_safe <- basename(run_name)
  if (!nzchar(run_name_safe) || run_name_safe %in% c(".", "..")) {
    stop("greml_heritability: 'run_name' must be a plain directory name ",
         "(no path separators, '.', or '..').")
  }

  # Create temp directory (disambiguate if it already exists)
  run_dir <- run_name_safe
  if (dir.exists(run_dir)) {
    run_dir <- paste0(run_dir, "_", format(Sys.time(), "%H%M%S"))
  }
  dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)

  # Register cleanup: remove the temp directory on exit unless save_run = TRUE
  on.exit(
    if (!save_run) unlink(run_dir, recursive = TRUE),
    add = TRUE
  )

  # ── Convert in-memory kinship matrix to GCTA binary GRM if needed ─────────
  # plinkFile::saveGRM() writes the three binary files (.grm.bin, .grm.N.bin,
  # .grm.id) that GCTA expects.
  if (!is.null(K_mat)) {
    if (!requireNamespace("plinkFile", quietly = TRUE)) {
      stop("greml_heritability (gcta backend): converting a kinship matrix to ",
           "GCTA binary GRM format requires the 'plinkFile' package.\n",
           "  Install it with: install.packages('plinkFile')\n",
           "  Alternatively, pre-convert the GRM with GCTA and supply the ",
           "path prefix via 'kinship_filename'.")
    }
    grm_prefix       <- file.path(run_dir, "kinship_matrix")
    plinkFile::saveGRM(grm_prefix, K_mat)
    kinship_filename <- grm_prefix
  }

  if (is.null(kinship_filename)) {
    stop("greml_heritability (gcta backend): 'kinship_filename' is required when ",
         "no 'kinship_matrix' is provided.")
  }

  # ── Build the GCTA phenotype file ─────────────────────────────────────────
  # GCTA format: two ID columns (FID, IID) followed by one column per trait.
  # Missing values are represented as "NA" (GCTA ignores them automatically).
  # Trait columns are renamed Y1, Y2, ... because GCTA selects them by index
  # via --mpheno rather than by name.

  # Split genotype IDs on "_" to obtain family (FID) and individual (IID) IDs.
  # If the ID has no "_", the same value is used for both FID and IID.
  raw_ids  <- rownames(matrix_unrolled)
  id_parts <- strsplit(raw_ids, "_", fixed = TRUE)
  fam_id   <- vapply(id_parts, `[[`, character(1L), 1L)
  indiv_id <- vapply(id_parts,
                     function(x) if (length(x) >= 2L) x[[2L]] else x[[1L]],
                     character(1L))

  pheno_tbl <- cbind.data.frame(
    FID = fam_id,
    IID = indiv_id,
    as.data.frame(matrix_unrolled, check.names = FALSE),
    stringsAsFactors = FALSE
  )
  # Rename trait columns to Y1, Y2, ... as required by GCTA --mpheno indexing
  colnames(pheno_tbl)[-(1:2)] <- paste0("Y", seq_len(ncol(matrix_unrolled)))

  pheno_file <- file.path(run_dir, "phenotypes.txt")
  write.table(pheno_tbl, file = pheno_file, sep = "\t",
              quote = FALSE, row.names = FALSE, na = "NA")

  # ── Loop over traits ───────────────────────────────────────────────────────
  n_traits <- ncol(matrix_unrolled)
  h2_vec <- se_vec <- Vg_vec <- Ve_vec <- Vp_vec <- rep(NA_real_, n_traits)

  for (i in seq_len(n_traits)) {

    out_prefix <- file.path(run_dir, paste0("trait_", i))

    # Build the GCTA command using shQuote() so paths with spaces are safe
    cmd <- paste(
      shQuote(gcta_bin),
      "--grm",    shQuote(kinship_filename),
      "--pheno",  shQuote(pheno_file),
      "--mpheno", i,
      "--reml",
      "--out",    shQuote(out_prefix)
    )

    if (verbose) {
      system(cmd)
    } else {
      # intern = TRUE captures stdout; ignore.stderr suppresses the console
      invisible(system(cmd, intern = TRUE, ignore.stderr = TRUE))
    }

    # ── Parse the .hsq output file ──────────────────────────────────────────
    # GCTA .hsq format (tab-separated, with header):
    #   Source    Variance    SE
    #   V(G)      <val>       <SE>
    #   V(e)      <val>       <SE>
    #   Vp        <val>       <SE>
    #   V(G)/Vp   <val>       <SE>   <- this is h^2
    #   ... (possibly more rows for n, logL, etc.)
    #
    # We look up rows by the "Source" field so the code is immune to GCTA
    # adding extra rows in future versions.
    hsq_file <- paste0(out_prefix, ".hsq")

    if (file.exists(hsq_file)) {
      hsq <- tryCatch(
        read.table(hsq_file, header = TRUE, sep = "\t",
                   fill = TRUE, stringsAsFactors = FALSE),
        error = function(e) NULL
      )

      if (!is.null(hsq) && "Source" %in% names(hsq)) {
        # Extract each row by name; nrow() guard handles missing rows gracefully
        get_field <- function(src, col) {
          row <- hsq[hsq$Source == src, ]
          if (nrow(row) == 1L) as.numeric(row[[col]]) else NA_real_
        }
        h2_vec[i] <- get_field("V(G)/Vp", "Variance")
        se_vec[i] <- get_field("V(G)/Vp", "SE")
        Vg_vec[i] <- get_field("V(G)",    "Variance")
        Ve_vec[i] <- get_field("V(e)",    "Variance")
        Vp_vec[i] <- get_field("Vp",      "Variance")
      }
    }
    # If .hsq is absent (REML failed / not converged) the NA defaults are kept.

    if (!is.null(pb)) pb$tick()
  }

  data.frame(h2 = h2_vec, se_h2 = se_vec,
             V_G = Vg_vec, V_e = Ve_vec, V_p = Vp_vec,
             stringsAsFactors = FALSE)
}


