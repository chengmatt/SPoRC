# Run Jitter Analysis for Model Diagnostics

Tests how sensitive the optimization is to its starting values by
perturbing the parameter vector with additive normal noise, refitting,
and recording the resulting series and diagnostics. Each iteration
perturbs the fixed effects, and the random effects too under
`jitter_random = TRUE`, optimizes with
[`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html), optionally
takes extra Newton steps, and extracts the reported quantities. The
starting values of every iteration can be returned alongside the results
with `return_start_values = TRUE`, so any run can be inspected or
repeated. Iterations run sequentially or in parallel through `future`.

## Usage

``` r
do_jitter(
  data,
  parameters,
  mapping,
  random = NULL,
  sd,
  n_jitter,
  n_newton_loops = 0,
  do_par,
  n_cores,
  par_vec = NULL,
  jitter_random = FALSE,
  return_start_values = FALSE
)
```

## Arguments

- data:

  List of model data for the `RTMB` objective function.

- parameters:

  Named list of parameters for `RTMB::MakeADFun()`.

- mapping:

  Named list of parameter mappings.

- random:

  Character vector of random-effect parameters.

- sd:

  Standard deviation of the additive normal noise.

- n_jitter:

  Number of jittered runs.

- n_newton_loops:

  Extra Newton steps after
  [`nlminb()`](https://rdrr.io/r/stats/nlminb.html) converges. Default
  0.

- do_par:

  Logical, whether the iterations run in parallel.

- n_cores:

  Parallel workers used when `do_par = TRUE`.

- par_vec:

  Optional starting values to jitter, either the fixed-effect vector
  (`length(obj$par)`, such as `fit$optim$par`) or the joint fixed and
  random vector (`length(obj$env$par)`, such as
  `fit$env$last.par.best`). Any other length is an error, and `NULL`
  uses the model's own start.

- jitter_random:

  Logical, whether the random effects are perturbed alongside the fixed
  effects. Only the fixed effects are searched by
  [`nlminb()`](https://rdrr.io/r/stats/nlminb.html), so the random draws
  move the starting point of the inner Laplace solve and check whether
  it settles on the same modes. Either way the inner solve starts from
  the random values in `par_vec`, or from the model's own start when
  `par_vec` holds none. Default `FALSE`.

- return_start_values:

  Logical, whether the starting values of each run are returned with the
  results. Default `FALSE`. `TRUE` returns a list of `jitter_df` and
  `start_values`, described under Value.

## Value

A data frame of the jitter results: the spawning stock biomass and
recruitment series of each run, with its jitter index, whether the
Hessian was positive definite, the joint negative log-likelihood, and
the maximum absolute gradient of the fixed effects. Under
`return_start_values = TRUE`, a list of two elements:

- `jitter_df`:

  The data frame above.

- `start_values`:

  A matrix of the values each run started from, with one row per jitter
  (row `i` is `jitter == i`) and one column per element of the model's
  joint parameter vector (`obj$env$par`), named by parameter. Random
  effects sit in their own columns, the ones named in `random`, and hold
  the unperturbed start when `jitter_random = FALSE`. Passing row `i`
  back as `par_vec` with `sd = 0` and `n_jitter = 1` repeats run `i`.

## See also

Other Model Diagnostics:
[`do_likelihood_profile()`](https://chengmatt.github.io/SPoRC/dev/reference/do_likelihood_profile.md),
[`do_retrospective()`](https://chengmatt.github.io/SPoRC/dev/reference/do_retrospective.md),
[`do_runs_test()`](https://chengmatt.github.io/SPoRC/dev/reference/do_runs_test.md),
[`get_catch_fits_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_catch_fits_plot.md),
[`get_comp_prop()`](https://chengmatt.github.io/SPoRC/dev/reference/get_comp_prop.md),
[`get_idx_fits()`](https://chengmatt.github.io/SPoRC/dev/reference/get_idx_fits.md),
[`get_idx_fits_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_idx_fits_plot.md),
[`get_model_rep_from_mcmc()`](https://chengmatt.github.io/SPoRC/dev/reference/get_model_rep_from_mcmc.md),
[`get_nLL_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_nLL_plot.md),
[`get_osa()`](https://chengmatt.github.io/SPoRC/dev/reference/get_osa.md),
[`get_retrospective_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_retrospective_plot.md),
[`get_retrospective_relative_difference()`](https://chengmatt.github.io/SPoRC/dev/reference/get_retrospective_relative_difference.md),
[`plot_resids()`](https://chengmatt.github.io/SPoRC/dev/reference/plot_resids.md)
