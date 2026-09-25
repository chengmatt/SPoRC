# Get Catch and Discard Fits Plot

Plots the observed catch and discard series against the predictions of
one or more SPoRC runs, pooled and population-specific. Pooled
predictions are summed across populations first, and years with a zero
observation are dropped from both layers so the lognormal intervals can
be built.

## Usage

``` r
get_catch_fits_plot(data, rep, model_names)
```

## Arguments

- data:

  List of length `n_models`, each a SPoRC data list. `ObsCatch`
  `[n_regions × n_yrs × n_seas × n_fish_fleets]` and `ObsDiscard` hold
  the pooled observations, and `Wt_Catch`, `Wt_Discard` and `ln_sigmaC`
  rebuild the observation error standard deviations for the intervals.
  `ObsCatch_pop` and `ObsDiscard_pop` with their weights are read when
  `UseCatch_pop` or `UseDiscard_pop` hold active elements.

- rep:

  List of length `n_models`, each a SPoRC report from `obj$report()`.
  `PredCatch` and `PredDiscard`
  `[n_pop × n_regions × n_yrs × n_seas × n_fish_fleets]` are summed
  across populations for the pooled series and read directly for the
  population-specific ones, with `ln_sigmaC` and `ln_sigmaC_pop` for the
  intervals.

- model_names:

  Character vector of length `n_models` of display names, used in the
  legend.

## Value

A list of four `ggplot` objects: the pooled catch fits, produced when
`UseCatch == 1` and faceted by season and fleet against region with free
y-scales; the population-specific catch fits, produced when
`UseCatch_pop == 1` and faceted by population, season and fleet against
region; and the two discard counterparts, produced when `UseDiscard` and
`UseDiscard_pop` are `1`. Observations are drawn as `geom_pointrange`
with 95% lognormal intervals and predictions as lines colored by model.
An element is `NULL` when its data source is unused.

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
