# dynamicGP

**Dynamic Mode Decomposition for Plant Phenomics**

<!-- badges: start -->
<!-- badges: end -->

`dynamicGP` fits Dynamic Mode Decomposition (DMD) models to multivariate plant phenotyping time-series data. The package was developed to support the analysis described in:

> Hobby D., Tong H., Heuermann M., Mbebi A.J., Laitinen R.A.E., Dell'Acqua M., Altman T. & Nikoloski Z. (2025). Predicting plant trait dynamics from genetic markers. *Nature Plants*, **11**(5), 1018–1027. https://doi.org/10.1038/s41477-025-01986-y

---

## Overview

`run_dmd()` decomposes a multivariate time series into DMD modes, eigenvalues, growth rates, and oscillation frequencies, returning the reconstructed state-transition operator **A**.

| Output field | Contents |
|---|---|
| `$A_full` | State-transition operator **A** `[traits × traits]` |
| `$lam` | DMD eigenvalues |
| `$Phi` | DMD modes |
| `$freq_cpd`, `$growth_day` | Oscillation frequencies and growth rates |
| `$resid_F` | Relative Frobenius residual |

Two fitting methods are available:

| Method | Description |
|---|---|
| `"standard"` | SVD-truncated eigen-decomposition (default) |
| `"direct"` | Direct pseudoinverse solution (`pracma` required) |

---

## Installation

```r
# Install from a local clone:
# install.packages("remotes")
remotes::install_local("dobby978/dynamicGP")

# Or via devtools:
devtools::install("dobby978/dynamicGP")
```

Required packages: `ggplot2 (>= 3.4.0)`, `lubridate (>= 1.9.0)`, `scales (>= 1.2.0)`, `tidyr (>= 1.3.0)`.

Optional: `Matrix (>= 1.5.0)` and `pracma (>= 2.4.0)` (Schur/direct DMD paths); `gaston` (GREML R backend); `plinkFile` (GREML GCTA backend); `dtw` (DTW distance in `compare_arrays()`); `expm` (continuous-time reconstruction).

---

## Quick start

```r
library(dynamicGP)

# Fit DMD  (X: traits × time matrix)
result <- run_dmd(X[, -ncol(X)], X[, -1L], var_thresh = 0.95)

# Key outputs
result$A_full   # state-transition operator
result$rho_A    # spectral radius
result$resid_F  # relative Frobenius residual

# Multi-step reconstruction
rec <- reconstruct_dmd(result$A_full, x0 = X[, 1L], n_steps = 30L)
rec$X_pred_arr
```

---

## Function reference

| Function | Purpose |
|---|---|
| `run_dmd()` | Single-genotype DMD |
| `greml_heritability()` | GREML SNP-based heritability (R or GCTA backend) |
| `arr_to_table()` | Flatten a 3-D operator array to a 2-D table |
| `table_to_arr()` | Reconstruct a 3-D array from a 2-D table |
| `concat_arrays()` | Concatenate 3-D arrays along the time dimension |
| `trajectory_properties()` | Compute trajectory statistics per trait per genotype |
| `spectral_clipper()` | Enforce spectral stability by eigenvalue clipping |
| `reconstruct_dmd()` | Multi-step state reconstruction from fitted operators |
| `compare_arrays()` | Accuracy metrics between observed and predicted arrays |
| `plot_trajectories()` | Mean-trajectory plots with SD ribbons |

---

## Citation

If you use `dynamicGP` in published work, please cite the associated paper:

```r
citation("dynamicGP")
```
