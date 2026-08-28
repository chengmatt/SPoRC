# Get Time Series Plots

Generates a suite of time series plots for key population dynamics
quantities (spawning stock biomass, dynamic unfished SSB, total biomass,
recruitment, and fishing mortality) across one or more SPoRC model runs.
Optionally overlays approximate 95 sdreport.

## Usage

``` r
get_ts_plot(rep, sd_rep, model_names, do_ci = TRUE)
```

## Arguments

- rep:

  List of length `n_models`, where each element is a SPoRC report list
  (i.e. the output of `obj$report()` after optimization).

- sd_rep:

  List of length `n_models`, where each element is a SPoRC sdreport list
  (i.e. the output of `sdreport(obj)`). Used to extract delta-method
  standard errors on log-scale quantities (`log_SSB`,
  `log_Dynamic_SSB0`, `log_Total_Biom`, `log_Rec`).

- model_names:

  Character vector of length `n_models` giving display names for each
  model run. Used as legend labels across all plots.

- do_ci:

  Logical. If `TRUE` (default), approximate 95 ribbons are added to SSB,
  total biomass, dynamic SSB0, and recruitment plots. Ribbons are
  computed on the natural scale as `exp(log(value) ± 1.96 * se)`.

## Value

A named list of seven `ggplot` objects:

- \[\[1\]\] comb_ts_plot:

  All quantities faceted by Type × Region.

- \[\[2\]\] f_ts_plot:

  Fishing mortality by fleet and season, faceted by Region.

- \[\[3\]\] rec_ts_plot:

  Recruitment by population, faceted by Region.

- \[\[4\]\] ssb_ts_plot:

  Spawning stock biomass by population, faceted by Region. Dynamic SSB0
  is excluded from this panel.

- \[\[5\]\] total_biom_plot:

  Total biomass by population, faceted by Region.

- \[\[6\]\] ssb0_plot:

  Dynamic unfished SSB (B0) by population, faceted by Region.

- \[\[7\]\] ssb_ssb0_plot:

  SSB and dynamic SSB0 overlaid on the same panel, distinguished by
  linetype, faceted by Region. Useful for visualizing depletion.

## See also

Other Plotting:
[`get_biological_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_biological_plot.md),
[`get_data_fitted_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_data_fitted_plot.md),
[`get_selex_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_selex_plot.md),
[`plot_all_basic()`](https://chengmatt.github.io/SPoRC/dev/reference/plot_all_basic.md),
[`theme_sablefish()`](https://chengmatt.github.io/SPoRC/dev/reference/theme_sablefish.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  plots <- get_ts_plot(
    rep        = list(rep1, rep2),
    sd_rep     = list(sdrep1, sdrep2),
    model_names = c("Base", "Sensitivity"),
    do_ci      = TRUE
  )
  plots[[4]] # SSB time series
  plots[[7]] # SSB vs dynamic B0
} # }
```
