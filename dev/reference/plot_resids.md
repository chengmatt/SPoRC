# Plot OSA residuals from outputs of get_osa

Generates diagnostic plots for one-step-ahead (OSA) residuals. Includes
QQ-plots with SDNR annotations and bubble plots showing residual
magnitude and sign.

## Usage

``` r
plot_resids(osa_results)
```

## Arguments

- osa_results:

  List obtained from get_osa() containing residuals dataframe.

## Value

List of plots: `sdnr_plot` (QQ-plot) and a second element that is a
`bubble_plot` (composition/tag residual magnitude and sign) or a
`resid_plot` (index-type residual vs. year).

## Details

Panels are faceted by every structural dimension the residual data frame
actually spans. Composition plots always facet by `region`/`sex` and
additionally facet by `fleet`, `pop`, and `seas` whenever
`osa_results$res` contains more than one of each (`seas` matters because
`year` + bin alone don't uniquely place a bubble-plot point when
compositions are collected in more than one season). Tagging plots facet
by `region`, `recovery_season`, `fleet`, and `pop_pool` (movement/tag
population-pooling group, when more than one is present) whenever those
span more than one level. The release-conditioned "tail" (non-recapture)
row has no `pop_pool` of its own and gets its own `"Tail (non-recap)"`
panel rather than being dropped or blended in. `recovery_year` +
`years_at_liberty` alone don't uniquely place a bubble-plot point –
cohorts released in different regions/seasons of the same release year
can share both – so the bubble plot (but not the QQ-plot/SDNR grouping)
jitters point positions slightly to keep coincidentally co-located
residuals visible rather than fully overplotted. Index-type residuals
(from `get_osa(..., index_source = ...)`, carrying an `idx_type` column
`%in% c("Catch","Discard","FishIdx","SrvIdx")` instead of `comp_type`)
facet by `region`, `season`, `fleet`, and `pop` whenever those span more
than one level, and pair the QQ-plot with a residual-vs-year point plot
instead of a bubble plot (there is no bin/age/length dimension to plot
against). Note: these are one-step-ahead residuals; for the simpler raw
log-scale (Pearson-style) index residual and the observed-vs-predicted
index fit, see
[`get_idx_fits`](https://chengmatt.github.io/SPoRC/dev/reference/get_idx_fits.md)
/
[`get_idx_fits_plot`](https://chengmatt.github.io/SPoRC/dev/reference/get_idx_fits_plot.md)
instead.

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
[`get_nLL_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_nLL_plot.md),
[`get_osa()`](https://chengmatt.github.io/SPoRC/dev/reference/get_osa.md),
[`get_retrospective_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_retrospective_plot.md),
[`get_retrospective_relative_difference()`](https://chengmatt.github.io/SPoRC/dev/reference/get_retrospective_relative_difference.md)
