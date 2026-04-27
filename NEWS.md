# dynamicGP News

## dynamicGP 0.2.0

*First public release.*

### New functions

* `run_dmd()` — fits DMD to a single multivariate time series, returning the
  state-transition operator **A**, DMD modes, eigenvalues, growth rates,
  oscillation frequencies, and fit diagnostics.
* `reconstruct_dmd()` — multi-step forward reconstruction from fitted DMD
  operators.
* `greml_heritability()` — GREML SNP-based heritability using a pure-R
  (`gaston`) or GCTA binary backend.
* `compare_arrays()`, `plot_trajectories()` — accuracy metrics and
  trajectory visualisation.
* `arr_to_table()`, `table_to_arr()`, `concat_arrays()`,
  `trajectory_properties()`, `spectral_clipper()` — utility helpers.

### Bug fixes

* **`run_dmd()` residual computation** — `pred_X2` is now taken as `Re()`
  before computing the Frobenius norm, avoiding spurious "imaginary parts
  discarded" warnings.
* **`run_dmd()` eigenvalue log** — `log(lam)` now coerces to complex before
  taking the logarithm, preventing NaN for real negative eigenvalues.
* **`greml_heritability()`** — `run_name` is now sanitised with `basename()`
  before use, preventing path-traversal values such as `"../../../tmp"` from
  writing files outside the working directory.  An explicit error is raised
  if the sanitised name is empty, `"."`, or `".."`.
