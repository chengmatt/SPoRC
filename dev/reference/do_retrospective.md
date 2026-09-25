# Run Retrospective Diagnostics for RTMB Models

Refits the model with terminal years removed one peel at a time,
truncating the inputs, applying any data lags and Francis reweighting,
and extracting spawning stock biomass and recruitment from each fit.

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

  Number of peels. `0` fits the full dataset only.

- data:

  List of data supplied to the RTMB model.

- parameters:

  List of model parameters.

- mapping:

  List of parameter mappings used during estimation.

- random:

  Character vector of random-effect parameters. Default `NULL`.

- do_par:

  Logical, whether the peels run in parallel. Default `FALSE`.

- n_cores:

  Cores used when `do_par = TRUE`.

- newton_loops:

  Newton optimization loops per fit. Default `3`.

- do_francis:

  Logical, whether Francis reweighting runs within each peel. Default
  `FALSE`.

- n_francis_iter:

  Francis reweighting iterations, required when `do_francis = TRUE`.

- nlminb_control:

  Control list passed to
  [`stats::nlminb`](https://rdrr.io/r/stats/nlminb.html). Default
  `list(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-15)`.

- do_sdrep:

  Logical, whether standard errors are computed through
  `RTMB::sdreport`. Default `FALSE`.

- fishidx_datalag, fishage_datalag, fishlen_datalag,
  fishage_discard_datalag, fishlen_discard_datalag, srvidx_datalag,
  srvage_datalag, srvlen_datalag:

  Integer arrays \\\[region \times fleet\]\\ of the lag applied to each
  pooled data source. Default zeros.

- fishidx_pop_datalag, fishage_pop_datalag, fishlen_pop_datalag,
  fishage_discard_pop_datalag, fishlen_discard_pop_datalag,
  srvidx_pop_datalag, srvage_pop_datalag, srvlen_pop_datalag:

  The population-specific counterparts, \\\[n\\pop \times region \times
  fleet\]\\. Default zeros.

- conv_tag_datalag:

  Integer lag applied to the conventional tagging data. Default `0`.

- return_models:

  Logical, whether the fitted objects are returned per peel. Default
  `FALSE`. `TRUE` returns a list of `retro_df` and `retro_models`, the
  latter indexed `peel_0`, `peel_1` and so on.

## Value

A long-format data frame of retrospective spawning stock biomass and
recruitment, with columns `Pop`, `Region`, `Year`, `Type` (`"SSB"` or
`"Recruitment"`), `peel` (0 for the full data), and `value`. Under
`do_sdrep = TRUE` it also holds `pdHess` and `max_grad`, the maximum
absolute gradient of the fixed effects.

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
