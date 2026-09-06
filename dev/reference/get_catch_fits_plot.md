# Get Catch and Discard Fits Plot

Plots observed catch and discard time series alongside model-predicted
values for one or more SPoRC model runs, for both pooled (region-level)
and population-specific data sources.

## Usage

``` r
get_catch_fits_plot(data, rep, model_names)
```

## Arguments

- data:

  List of length `n_models`, where each element is a SPoRC data list.
  `ObsCatch` `[n_regions × n_yrs × n_seas × n_fish_fleets]` provides
  pooled observed catch, and `ObsDiscard` provides pooled observed
  discards. `Wt_Catch` and `Wt_Discard`, together with `ln_sigmaC`, are
  used to reconstruct observation-error standard deviations for
  confidence interval construction.

  When `UseCatch_pop` or `UseDiscard_pop` contain active elements,
  population-level data (`ObsCatch_pop`, `ObsDiscard_pop`) and
  corresponding weights are also used.

- rep:

  List of length `n_models`, where each element is a SPoRC model report
  (output of `obj$report()` after optimization). `PredCatch` and
  `PredDiscard` `[n_pop × n_regions × n_yrs × n_seas × n_fish_fleets]`
  are summed across populations for pooled trajectories and used
  directly for population-specific trajectories. `ln_sigmaC` and
  `ln_sigmaC_pop` are used for confidence interval construction.

- model_names:

  Character vector of length `n_models` giving display names for each
  model run. Used in the legend for predicted trajectories.

## Value

A list of `ggplot` objects:

- \[\[1\]\] catch_fit_rg_plot:

  Pooled catch fits (region-level). Produced when `UseCatch == 1`.
  Observations shown as `geom_pointrange` with 95% lognormal confidence
  intervals; predictions shown as lines colored by model. Faceted by
  (Season × Fleet) × Region with free y-scales.

- \[\[2\]\] catch_fit_pop_plot:

  Population-specific catch fits. Produced when `UseCatch_pop == 1`.
  Faceted by (Population × Season × Fleet) × Region. Returns `NULL` if
  not used.

- \[\[3\]\] discard_fit_rg_plot:

  Pooled discard fits (region-level). Produced when `UseDiscard == 1`.
  Same structure as catch plots.

- \[\[4\]\] discard_fit_pop_plot:

  Population-specific discard fits. Produced when `UseDiscard_pop == 1`.
  Returns `NULL` if not used.

## Details

Pooled predictions are summed across populations before comparison with
observed data. Years where observed values are zero are excluded from
both observed and predicted layers to avoid issues in lognormal
confidence interval construction.

This function produces diagnostic plots for both catch and discard data.
Observed and predicted time series are shown for each, with lognormal
confidence intervals derived from `ln_sigmaC` (or population-level
equivalents) and sampling weights.

Observations equal to zero are excluded prior to plotting to avoid
issues in log-space confidence interval construction.

Predicted catch and discard are aggregated across populations for pooled
diagnostics, while population-specific plots use unaggregated outputs.

## See also

Other Model Diagnostics:
[`do_jitter()`](https://chengmatt.github.io/SPoRC/dev/reference/do_jitter.md),
[`do_likelihood_profile()`](https://chengmatt.github.io/SPoRC/dev/reference/do_likelihood_profile.md),
[`do_retrospective()`](https://chengmatt.github.io/SPoRC/dev/reference/do_retrospective.md),
[`do_runs_test()`](https://chengmatt.github.io/SPoRC/dev/reference/do_runs_test.md),
[`get_comp_prop()`](https://chengmatt.github.io/SPoRC/dev/reference/get_comp_prop.md),
[`get_idx_fits()`](https://chengmatt.github.io/SPoRC/dev/reference/get_idx_fits.md),
[`get_idx_fits_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_idx_fits_plot.md),
[`get_model_rep_from_mcmc()`](https://chengmatt.github.io/SPoRC/dev/reference/get_model_rep_from_mcmc.md),
[`get_nLL_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_nLL_plot.md),
[`get_osa()`](https://chengmatt.github.io/SPoRC/dev/reference/get_osa.md),
[`get_retrospective_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_retrospective_plot.md),
[`get_retrospective_relative_difference()`](https://chengmatt.github.io/SPoRC/dev/reference/get_retrospective_relative_difference.md),
[`plot_resids()`](https://chengmatt.github.io/SPoRC/dev/reference/plot_resids.md)
