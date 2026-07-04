# Run Jitter Analysis for Model Diagnostics

Performs a jitter analysis to evaluate sensitivity of model optimization
to starting parameter values. The function repeatedly perturbs the
parameter vector with additive normal noise, refits the model, and
records resulting time series and diagnostic metrics.

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
  par_vec = NULL
)
```

## Arguments

- data:

  A list of model data used to construct the `RTMB` objective function.

- parameters:

  A named list of model parameters used to initialize
  `RTMB::MakeADFun()`.

- mapping:

  A named list defining parameter mappings for `RTMB::MakeADFun()`.

- random:

  Character vector specifying random-effect parameters.

- sd:

  Numeric value specifying the standard deviation of the additive normal
  noise used to jitter parameters.

- n_jitter:

  Integer specifying the number of jittered optimization runs.

- n_newton_loops:

  Integer specifying the number of additional Newton optimization steps
  performed after [`nlminb()`](https://rdrr.io/r/stats/nlminb.html)
  convergence. Default = 0.

- do_par:

  Logical indicating whether jitter iterations should be executed in
  parallel.

- n_cores:

  Integer specifying the number of parallel workers to use when
  `do_par = TRUE`.

- par_vec:

  Optional numeric vector of parameter values used as the starting point
  for jittering. If `NULL`, the model's default starting parameter
  vector is jittered. If provided, jittering is applied to this vector
  (for example, the maximum likelihood estimates).

## Value

A `data.frame` containing jitter iteration results. The output includes
time series of spawning stock biomass (SSB) and recruitment, along with
diagnostic information for each jitter run, including:

- Jitter index

- Whether the Hessian is positive definite

- Joint negative log-likelihood

- Maximum absolute gradient of fixed effects

## Details

Each jitter iteration:

- Perturbs the starting parameter vector with random normal noise.

- Optimizes the objective function using
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html).

- Optionally performs additional Newton steps to refine the solution.

- Extracts reported quantities (e.g., spawning biomass and recruitment)
  and diagnostic statistics.

The analysis can be executed sequentially or in parallel using the
`future` framework.

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
