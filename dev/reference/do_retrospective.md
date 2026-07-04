# Run Retrospective Diagnostics for RTMB Models

Conducts retrospective analyses by sequentially removing terminal years
("peels") from the dataset and refitting the model. For each peel, the
function truncates the model inputs, optionally applies data lags and
Francis composition reweighting, fits the model, and extracts estimates
of spawning stock biomass (SSB) and recruitment.

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
  fishage_discard_datalag = array(0, dim = c(data$n_regions, data$n_fish_fleets)),
  fishlen_discard_datalag = array(0, dim = c(data$n_regions, data$n_fish_fleets)),
  srvidx_datalag = array(0, dim = c(data$n_regions, data$n_srv_fleets)),
  srvage_datalag = array(0, dim = c(data$n_regions, data$n_srv_fleets)),
  srvlen_datalag = array(0, dim = c(data$n_regions, data$n_srv_fleets)),
  fishidx_pop_datalag = array(0, dim = c(data$n_pop, data$n_regions, data$n_fish_fleets)),
  fishage_pop_datalag = array(0, dim = c(data$n_pop, data$n_regions, data$n_fish_fleets)),
  fishlen_pop_datalag = array(0, dim = c(data$n_pop, data$n_regions, data$n_fish_fleets)),
  fishage_discard_pop_datalag = array(0, dim = c(data$n_pop, data$n_regions,
    data$n_fish_fleets)),
  fishlen_discard_pop_datalag = array(0, dim = c(data$n_pop, data$n_regions,
    data$n_fish_fleets)),
  srvidx_pop_datalag = array(0, dim = c(data$n_pop, data$n_regions, data$n_srv_fleets)),
  srvage_pop_datalag = array(0, dim = c(data$n_pop, data$n_regions, data$n_srv_fleets)),
  srvlen_pop_datalag = array(0, dim = c(data$n_pop, data$n_regions, data$n_srv_fleets)),
  conv_tag_datalag = 0,
  return_models = FALSE
)
```

## Arguments

- n_retro:

  Integer specifying the number of retrospective peels to perform. A
  value of `n_retro = 0` fits the model using the full dataset only.

- data:

  List containing the data supplied to the RTMB model.

- parameters:

  List containing the model parameters.

- mapping:

  List defining parameter mappings used during estimation.

- random:

  Character vector identifying random-effect parameters in the model.
  Default is `NULL`.

- do_par:

  Logical indicating whether retrospective peels should be run in
  parallel. Default is `FALSE`.

- n_cores:

  Integer specifying the number of cores to use when `do_par = TRUE`.

- newton_loops:

  Integer specifying the number of Newton optimization loops used during
  model fitting. Default is `3`.

- do_francis:

  Logical indicating whether Francis composition reweighting should be
  applied within each retrospective peel. Default is `FALSE`.

- n_francis_iter:

  Integer specifying the number of Francis reweighting iterations.
  Required if `do_francis = TRUE`.

- nlminb_control:

  List of control arguments passed to
  [`stats::nlminb`](https://rdrr.io/r/stats/nlminb.html) during model
  fitting. Default is
  `list(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-15)`.

- do_sdrep:

  Logical indicating whether standard errors should be calculated using
  `RTMB::sdreport`. Default is `FALSE`.

- fishidx_datalag:

  Integer array specifying lags applied to fishery index data \\\[region
  \times fleet\]\\. Default is zeros.

- fishage_datalag:

  Integer array specifying lags applied to fishery age-composition data
  \\\[region \times fleet\]\\. Default is zeros.

- fishlen_datalag:

  Integer array specifying lags applied to fishery length-composition
  data \\\[region \times fleet\]\\. Default is zeros.

- fishage_discard_datalag:

  Integer array specifying lags applied to fishery discard
  age-composition data \\\[region \times fleet\]\\. Default is zeros.

- fishlen_discard_datalag:

  Integer array specifying lags applied to fishery discard
  length-composition data \\\[region \times fleet\]\\. Default is zeros.

- srvidx_datalag:

  Integer array specifying lags applied to survey index data \\\[region
  \times fleet\]\\. Default is zeros.

- srvage_datalag:

  Integer array specifying lags applied to survey age-composition data
  \\\[region \times fleet\]\\. Default is zeros.

- srvlen_datalag:

  Integer array specifying lags applied to survey length-composition
  data \\\[region \times fleet\]\\. Default is zeros.

- fishidx_pop_datalag:

  Integer array specifying lags applied to population-specific fishery
  index data \\\[n\\pop \times region \times fleet\]\\. Default is
  zeros.

- fishage_pop_datalag:

  Integer array specifying lags applied to population-specific fishery
  age-composition data \\\[n\\pop \times region \times fleet\]\\.
  Default is zeros.

- fishlen_pop_datalag:

  Integer array specifying lags applied to population-specific fishery
  length-composition data \\\[n\\pop \times region \times fleet\]\\.
  Default is zeros.

- fishage_discard_pop_datalag:

  Integer array specifying lags applied to population-specific fishery
  discard age-composition data \\\[n\\pop \times region \times
  fleet\]\\. Default is zeros.

- fishlen_discard_pop_datalag:

  Integer array specifying lags applied to population-specific fishery
  discard length-composition data \\\[n\\pop \times region \times
  fleet\]\\. Default is zeros.

- srvidx_pop_datalag:

  Integer array specifying lags applied to population-specific survey
  index data \\\[n\\pop \times region \times fleet\]\\. Default is
  zeros.

- srvage_pop_datalag:

  Integer array specifying lags applied to population-specific survey
  age-composition data \\\[n\\pop \times region \times fleet\]\\.
  Default is zeros.

- srvlen_pop_datalag:

  Integer array specifying lags applied to population-specific survey
  length-composition data \\\[n\\pop \times region \times fleet\]\\.
  Default is zeros.

- conv_tag_datalag:

  Integer specifying the lag applied to conventional tagging data.
  Default is `0`.

- return_models:

  Logical indicating whether fitted model objects should be returned for
  each retrospective peel. Default is `FALSE`. When `TRUE`, the function
  returns a named list with two elements: `retro_df` (the long-format
  `data.frame` of SSB and recruitment estimates) and `retro_models` (a
  named list of fitted model objects, indexed as `peel_0`, `peel_1`,
  ..., `peel_n`). When `FALSE`, only the `data.frame` is returned.

## Value

A long-format `data.frame` containing retrospective estimates of
spawning stock biomass and recruitment. Columns include:

- `Pop` – Population index.

- `Region` – Region index.

- `Year` – Model year.

- `Type` – Quantity reported (`"SSB"` or `"Recruitment"`).

- `peel` – Retrospective peel number (0 = full data, 1 = one-year peel,
  etc.).

- `value` – Estimated value of the quantity.

- `pdHess` – Logical indicator of positive-definite Hessian (only
  present when `do_sdrep = TRUE`).

- `max_grad` – Maximum absolute gradient of fixed effects (only present
  when `do_sdrep = TRUE`).

## Details

Retrospective analyses are commonly used to evaluate the stability of
model estimates through time and to diagnose potential model
misspecification.

## See also

Other Model Diagnostics:
[`do_jitter()`](https://chengmatt.github.io/SPoRC/dev/reference/do_jitter.md),
[`do_likelihood_profile()`](https://chengmatt.github.io/SPoRC/dev/reference/do_likelihood_profile.md),
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
