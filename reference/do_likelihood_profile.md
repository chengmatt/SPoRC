# Run Likelihood Profile

Run Likelihood Profile

## Usage

``` r
do_likelihood_profile(
  data,
  parameters,
  mapping,
  random = NULL,
  what,
  idx = NULL,
  min_val,
  max_val,
  inc = 0.05,
  do_par = FALSE,
  n_cores = NULL
)
```

## Arguments

- data:

  data list from model

- parameters:

  parameter list from model

- mapping:

  mapping list from model

- random:

  character vector of random effects to estimate

- what:

  parameter name we want to profile

- idx:

  Index for an parameter array, pointing to the value we want to map off
  (index is relative to a flattened array)

- min_val:

  minimum value of profile

- max_val:

  maximum value of profile

- inc:

  increment value between min and max value

- do_par:

  logical, whether to use parallel processing (default FALSE)

- n_cores:

  integer, number of cores to use for parallel processing (default is
  detectCores() - 1)

## Value

Returns a list of likelihood profiled values for each data component
with their respective dimensions (e.g., likelihood profiles by fleet,
region, year, etc.) as well likelihood profiles for each data component,
aggregated across all their respective dimensions.

## See also

Other Model Diagnostics:
[`do_jitter()`](https://chengmatt.github.io/SPoRC/reference/do_jitter.md),
[`do_retrospective()`](https://chengmatt.github.io/SPoRC/reference/do_retrospective.md),
[`do_runs_test()`](https://chengmatt.github.io/SPoRC/reference/do_runs_test.md),
[`get_catch_fits_plot()`](https://chengmatt.github.io/SPoRC/reference/get_catch_fits_plot.md),
[`get_comp_prop()`](https://chengmatt.github.io/SPoRC/reference/get_comp_prop.md),
[`get_idx_fits()`](https://chengmatt.github.io/SPoRC/reference/get_idx_fits.md),
[`get_idx_fits_plot()`](https://chengmatt.github.io/SPoRC/reference/get_idx_fits_plot.md),
[`get_model_rep_from_mcmc()`](https://chengmatt.github.io/SPoRC/reference/get_model_rep_from_mcmc.md),
[`get_nLL_plot()`](https://chengmatt.github.io/SPoRC/reference/get_nLL_plot.md),
[`get_osa()`](https://chengmatt.github.io/SPoRC/reference/get_osa.md),
[`get_retrospective_plot()`](https://chengmatt.github.io/SPoRC/reference/get_retrospective_plot.md),
[`get_retrospective_relative_difference()`](https://chengmatt.github.io/SPoRC/reference/get_retrospective_relative_difference.md),
[`plot_resids()`](https://chengmatt.github.io/SPoRC/reference/plot_resids.md)
