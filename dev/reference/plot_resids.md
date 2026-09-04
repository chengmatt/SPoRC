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

List of two plots, named and also safe to index by position: `sdnr_plot`
(QQ-plot) first, then `bubble_plot` (composition/tag residual magnitude
and sign) or `resid_plot` (index-type residual vs. year).

## Details

Panels are faceted by every structural dimension the residual data frame
actually spans. Composition plots always facet by `region`/`sex` and
additionally facet by `fleet`, `pop`, and `seas` whenever
`osa_results$res` contains more than one of each (`seas` matters because
`year` + bin alone don't uniquely place a bubble-plot point when
compositions are collected in more than one season). Tagging plots only
show QQ plots given the number of dimensions in tagging data. Index-type
residuals (from `get_osa(..., index_source = ...)`, carrying an
`idx_type` column `%in% c("Catch","Discard","FishIdx","SrvIdx")` instead
of `comp_type`) facet by `region`, `season`, `fleet`, and `pop` whenever
those span more than one level, and pair the QQ-plot with a
residual-vs-year point plot instead of a bubble plot (there is no
bin/age/length dimension to plot against). Note: these are
one-step-ahead residuals; for the simpler raw log-scale (Pearson-style)
index residual and the observed-vs-predicted index fit, see
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
