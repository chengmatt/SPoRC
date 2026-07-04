# Get Index Fits Plot

Plots observed survey and fishery indices alongside model-predicted
values for one or more SPoRC model runs. Observed values are shown as
points with approximate 95 lines coloured by model. Years where the
observed index is zero (i.e. `Use*Idx = 0`) are excluded from both the
points and lines.

## Usage

``` r
get_idx_fits_plot(data, rep, model_names)
```

## Arguments

- data:

  List of length `n_models`, where each element is a SPoRC data list.
  Passed to `get_idx_fits` along with `rep[[i]]`; year labels are taken
  from `data[[i]]$years`.

- rep:

  List of length `n_models`, where each element is a SPoRC report list
  (i.e. the output of `obj$report()` after optimisation). Predicted
  index values are extracted internally via `get_idx_fits`.

- model_names:

  Character vector of length `n_models` giving display names for each
  model run. Used as the colour legend label on predicted trajectories.

## Value

A single `ggplot` object. Observed indices are shown as
`geom_pointrange` (black) with lower and upper confidence interval
bounds from `get_idx_fits`. Predicted indices are shown as `geom_line`
coloured by model.

## See also

Other Model Diagnostics:
[`do_jitter()`](https://chengmatt.github.io/SPoRC/dev/reference/do_jitter.md),
[`do_likelihood_profile()`](https://chengmatt.github.io/SPoRC/dev/reference/do_likelihood_profile.md),
[`do_retrospective()`](https://chengmatt.github.io/SPoRC/dev/reference/do_retrospective.md),
[`do_runs_test()`](https://chengmatt.github.io/SPoRC/dev/reference/do_runs_test.md),
[`get_catch_fits_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_catch_fits_plot.md),
[`get_comp_prop()`](https://chengmatt.github.io/SPoRC/dev/reference/get_comp_prop.md),
[`get_idx_fits()`](https://chengmatt.github.io/SPoRC/dev/reference/get_idx_fits.md),
[`get_model_rep_from_mcmc()`](https://chengmatt.github.io/SPoRC/dev/reference/get_model_rep_from_mcmc.md),
[`get_nLL_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_nLL_plot.md),
[`get_osa()`](https://chengmatt.github.io/SPoRC/dev/reference/get_osa.md),
[`get_retrospective_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_retrospective_plot.md),
[`get_retrospective_relative_difference()`](https://chengmatt.github.io/SPoRC/dev/reference/get_retrospective_relative_difference.md),
[`plot_resids()`](https://chengmatt.github.io/SPoRC/dev/reference/plot_resids.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  get_idx_fits_plot(
    data        = list(data1, data2),
    rep         = list(rep1, rep2),
    model_names = c("Base", "Sensitivity")
  )
} # }
```
