# Get Plot of Negative Log Likelihood Values

Extracts, weights, and visualises all negative log-likelihood (nLL)
components from one or more SPoRC model runs. Likelihood weights stored
in the data list are applied to the relevant components (catch, indices,
recruitment, tagging, fishing mortality) before plotting, so reported
values reflect the actual contribution of each component to the joint
nLL. Components with a value of zero are silently excluded from both the
plot and table.

## Usage

``` r
get_nLL_plot(data, rep, model_names)
```

## Arguments

- data:

  List of length `n_models`, where each element is a SPoRC data list.
  Likelihood weights are extracted from `Wt_Catch`, `Wt_F`, `Wt_Rec`,
  `Wt_Tagging`, `Wt_SrvIdx`, and `Wt_FishIdx`.

- rep:

  List of length `n_models`, where each element is a SPoRC report list
  (i.e. the output of `obj$report()` after optimisation). The following
  nLL components are extracted: `jnLL`, `h_nLL`, `M_nLL`,
  `rec_region_prop_nLL`, `Rec_nLL`, `Init_Rec_nLL`, `sel_nLL`,
  `conv_fish_tag_nLL`, `Catch_nLL`, `Fmort_nLL`, `srv_q_nLL`,
  `fish_q_nLL`, `SrvIdx_nLL`, `FishIdx_nLL`, `TagRep_nLL`,
  `Movement_nLL`, `SrvAgeComps_nLL`, `FishAgeComps_nLL`,
  `SrvLenComps_nLL`, `FishLenComps_nLL`, and the population-specific
  counterparts `Catch_pop_nLL`, `FishIdx_pop_nLL`, `SrvIdx_pop_nLL`,
  `FishAgeComps_pop_nLL`, `FishLenComps_pop_nLL`, `SrvAgeComps_pop_nLL`,
  `SrvLenComps_pop_nLL`. Missing components are handled via
  `safe_extract`.

- model_names:

  Character vector of length `n_models` giving display names for each
  model run. Used as facet labels on the bar chart.

## Value

A list of two objects:

- \[\[1\]\] nLL_plot:

  A stacked bar chart (`ggplot`) with one facet per model. Bars show the
  weighted nLL contribution of each component, coloured by component
  type (Prior, Penalty, Catch, Index, Age, Length, Tagging, jnLL).
  Components with value 0 are excluded.

- \[\[2\]\] table_plot:

  A `ggdraw` table (via `gridExtra` and `cowplot`) showing weighted nLL
  values in wide format, with one column per non-zero component and one
  row per model.

## See also

Other Model Diagnostics:
[`do_jitter()`](https://chengmatt.github.io/SPoRC/dev/reference/do_jitter.md),
[`do_likelihood_profile()`](https://chengmatt.github.io/SPoRC/dev/reference/do_likelihood_profile.md),
[`do_retrospective()`](https://chengmatt.github.io/SPoRC/dev/reference/do_retrospective.md),
[`do_runs_test()`](https://chengmatt.github.io/SPoRC/dev/reference/do_runs_test.md),
[`get_catch_fits_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_catch_fits_plot.md),
[`get_comp_prop()`](https://chengmatt.github.io/SPoRC/dev/reference/get_comp_prop.md),
[`get_idx_fits()`](https://chengmatt.github.io/SPoRC/dev/reference/get_idx_fits.md),
[`get_idx_fits_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_idx_fits_plot.md),
[`get_model_rep_from_mcmc()`](https://chengmatt.github.io/SPoRC/dev/reference/get_model_rep_from_mcmc.md),
[`get_osa()`](https://chengmatt.github.io/SPoRC/dev/reference/get_osa.md),
[`get_retrospective_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_retrospective_plot.md),
[`get_retrospective_relative_difference()`](https://chengmatt.github.io/SPoRC/dev/reference/get_retrospective_relative_difference.md),
[`plot_resids()`](https://chengmatt.github.io/SPoRC/dev/reference/plot_resids.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  out <- get_nLL_plot(
    data        = list(data1, data2),
    rep         = list(rep1, rep2),
    model_names = c("Base", "Sensitivity")
  )
  out[[1]]  # bar chart
  out[[2]]  # table
} # }
```
