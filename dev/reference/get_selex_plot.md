# Get Fishery and Survey Selectivity Plots

Generates age- or length-based selectivity plots for all fishery and
survey fleets across one or more SPoRC model runs. Selectivity curves
are faceted by sex, region, and fleet, and coloured by year when
multiple years are requested, allowing visualisation of time-varying
selectivity alongside model comparisons.

## Usage

``` r
get_selex_plot(rep, model_names, Selex_Type = "age", year_indx = NULL)
```

## Arguments

- rep:

  List of length `n_models`, where each element is a SPoRC report list
  (i.e. the output of `obj$report()` after optimisation). Must contain
  `fish_sel` and `srv_sel` (age-based) and/or `fish_sel_l` and
  `srv_sel_l` (length-based) depending on `Selex_Type`.

- model_names:

  Character vector of length `n_models` giving display names for each
  model run. Used as the linetype legend label.

- Selex_Type:

  Character string specifying whether to plot age- or length-based
  selectivity. One of `"age"` (default) or `"length"`. Must match the
  `Selex_Type` used when fitting the model (i.e. `fish_sel_l` /
  `srv_sel_l` are only populated when `Selex_Type = 1` in the model).

- year_indx:

  Integer or integer vector of year indices to include in the plot. When
  multiple years are supplied, curves are coloured by year on a
  continuous viridis scale, making time-variation in selectivity
  visible. If `NULL` (default), only the terminal year is plotted.

## Value

A list of two `ggplot` objects:

- \[\[1\]\] fish_sel_plot:

  Fishery selectivity curves faceted by Sex × (Region + Fleet). Lines
  are coloured by year and distinguished by linetype across models.

- \[\[2\]\] srv_sel_plot:

  Survey selectivity curves with identical faceting and aesthetic
  structure as the fishery plot.

## See also

Other Plotting:
[`get_biological_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_biological_plot.md),
[`get_data_fitted_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_data_fitted_plot.md),
[`get_ts_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_ts_plot.md),
[`plot_all_basic()`](https://chengmatt.github.io/SPoRC/dev/reference/plot_all_basic.md),
[`theme_sablefish()`](https://chengmatt.github.io/SPoRC/dev/reference/theme_sablefish.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  # Terminal year only
  plots <- get_selex_plot(list(rep1, rep2), c("Base", "Sensitivity"))

  # All years to visualise time-varying selectivity
  plots <- get_selex_plot(list(rep1), "Base", Selex_Type = "age", year_indx = 1:40)

  plots[[1]] # fishery selectivity
  plots[[2]] # survey selectivity
} # }
```
