# Get Data Fitted to Plot

Produces a dot-plot showing which data types were fitted to in each
year, region, and model. Each point represents a year in which a given
data source was active (i.e. the corresponding `Use*` indicator equals
1). This is useful for quickly auditing data availability and comparing
data structures across model configurations.

## Usage

``` r
get_data_fitted_plot(data, model_names)
```

## Arguments

- data:

  List of length `n_models`, where each element is a SPoRC data list.
  Pooled `Use*` indicator arrays are extracted with dimensions
  `[n_regions × n_years × n_seas × n_fleets]`: `UseSrvLenComps`,
  `UseSrvAgeComps`, `UseFishLenComps`, `UseFishAgeComps`, `UseCatch`,
  `UseFishIdx`, `UseSrvIdx`. Population-specific `Use*_pop` indicator
  arrays with dimensions
  `[n_pop × n_regions × n_years × n_seas × n_fleets]` are included when
  any element equals 1: `UseFishAgeComps_pop`, `UseFishLenComps_pop`,
  `UseFishIdx_pop`, `UseSrvAgeComps_pop`, `UseSrvLenComps_pop`,
  `UseSrvIdx_pop`. If `use_conv_fish_tagging` contains any 1s,
  `conv_tag_release_indicator` is also used to construct a tagging
  activity indicator array.

- model_names:

  Character vector of length `n_models` giving display names for each
  model run. Used as row facet labels.

## Value

A single `ggplot` object: a dot-plot with Year on the x-axis and data
type (labelled by source, population where applicable, season, and
fleet) on the y-axis, faceted by Model × Region. Points appear only in
years where the corresponding `Use*` indicator is 1. The legend is
suppressed; data types are distinguished by y-axis position and fill
colour.

## See also

Other Plotting:
[`get_biological_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_biological_plot.md),
[`get_selex_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_selex_plot.md),
[`get_ts_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_ts_plot.md),
[`plot_all_basic()`](https://chengmatt.github.io/SPoRC/dev/reference/plot_all_basic.md),
[`theme_sablefish()`](https://chengmatt.github.io/SPoRC/dev/reference/theme_sablefish.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  get_data_fitted_plot(
    data        = list(data1, data2),
    model_names = c("Base", "Sensitivity")
  )
} # }
```
