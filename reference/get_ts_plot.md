# Get Time Series Plots

Get Time Series Plots

## Usage

``` r
get_ts_plot(rep, sd_rep, model_names, do_ci = TRUE)
```

## Arguments

- rep:

  List of n_models of \`SPoRC\` report lists

- sd_rep:

  List of n_models of \`SPoRC\` sdreport lists

- model_names:

  Vector of model names

- do_ci:

  Boolean for whether confidence intervals are plotted

## Value

Plots of spawning biomass, dynamic b0, total biomass, recruitment, and
fishing mortality time-series across models

## See also

Other Plotting:
[`get_biological_plot()`](https://chengmatt.github.io/SPoRC/reference/get_biological_plot.md),
[`get_data_fitted_plot()`](https://chengmatt.github.io/SPoRC/reference/get_data_fitted_plot.md),
[`get_selex_plot()`](https://chengmatt.github.io/SPoRC/reference/get_selex_plot.md),
[`plot_all_basic()`](https://chengmatt.github.io/SPoRC/reference/plot_all_basic.md),
[`theme_sablefish()`](https://chengmatt.github.io/SPoRC/reference/theme_sablefish.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  get_ts_plot(list(rep1, rep2), list(sd_rep1, sd_rep2), c("Model1", "Model2"), do_ci = TRUE)
} # }
```
