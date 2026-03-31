# Run retrospective analyses for RTMB models

Performs retrospective peels by truncating the input data, optionally
applying Francis reweighting and parallelization, and returns estimates
of spawning stock biomass (SSB) and recruitment for each peel.

## Usage

``` r
do_retrospective(
  n_retro,
  data,
  parameters,
  mapping,
  random = NULL,
  do_par,
  n_cores,
  newton_loops = 3,
  do_francis = FALSE,
  n_francis_iter = NULL,
  nlminb_control = list(iter.max = 1e+05, eval.max = 1e+05, rel.tol = 1e-15),
  do_sdrep = FALSE,
  fishidx_datalag = array(0, dim = c(data$n_regions, data$n_fish_fleets)),
  fishage_datalag = array(0, dim = c(data$n_regions, data$n_fish_fleets)),
  fishlen_datalag = array(0, dim = c(data$n_regions, data$n_fish_fleets)),
  srvidx_datalag = array(0, dim = c(data$n_regions, data$n_srv_fleets)),
  srvage_datalag = array(0, dim = c(data$n_regions, data$n_srv_fleets)),
  srvlen_datalag = array(0, dim = c(data$n_regions, data$n_srv_fleets)),
  tag_datalag = 0
)
```

## Arguments

- n_retro:

  Integer. Number of retrospective peels to perform.

- data:

  List. Data input for the RTMB model.

- parameters:

  List. Parameter values for the RTMB model.

- mapping:

  List. Mapping information for the RTMB model.

- random:

  Character vector. Names of random effects in the model. Default is
  `NULL`.

- do_par:

  Logical. Whether to run retrospective peels in parallel. Default is
  `FALSE`.

- n_cores:

  Integer. Number of cores to use for parallel execution if
  `do_par = TRUE`.

- newton_loops:

  Integer. Number of Newton loops to run during model fitting. Default
  is 3.

- do_francis:

  Logical. Whether to apply Francis reweighting within each
  retrospective peel. Default is `FALSE`.

- n_francis_iter:

  Integer. Number of Francis reweighting iterations. Required if
  `do_francis = TRUE`.

- nlminb_control:

  List. Control parameters passed to `nlminb` during model fitting.
  Default is `list(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-15)`.

- do_sdrep:

  Logical. Whether to return standard errors from `sdreport`. Default is
  `FALSE`.

- fishidx_datalag:

  Integer array. Lags for fishery index data \[regions x fleets\].
  Default is zeros.

- fishage_datalag:

  Integer array. Lags for fishery age composition data \[regions x
  fleets\]. Default is zeros.

- fishlen_datalag:

  Integer array. Lags for fishery length composition data \[regions x
  fleets\]. Default is zeros.

- srvidx_datalag:

  Integer array. Lags for survey index data \[regions x fleets\].
  Default is zeros.

- srvage_datalag:

  Integer array. Lags for survey age composition data \[regions x
  fleets\]. Default is zeros.

- srvlen_datalag:

  Integer array. Lags for survey length composition data \[regions x
  fleets\]. Default is zeros.

- tag_datalag:

  Integer. Lag for tagging data. Default is 0.

## Value

A `data.frame` containing retrospective estimates of SSB and
recruitment. Columns include:

- `Region`: Region index.

- `Year`: Year index.

- `Type`: "SSB" or "Recruitment".

- `peel`: Peel number (0 = full data, 1 = 1-year peel, etc.).

- `value`: Estimated value of SSB or recruitment.

- `pdHess` and `max_grad` (optional): Information from `sdreport` if
  `do_sdrep = TRUE`.

## See also

Other Model Diagnostics:
[`do_jitter()`](https://chengmatt.github.io/SPoRC/reference/do_jitter.md),
[`do_likelihood_profile()`](https://chengmatt.github.io/SPoRC/reference/do_likelihood_profile.md),
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
