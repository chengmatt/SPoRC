# Get Fishery and Survey Selectivity Plots

Get Fishery and Survey Selectivity Plots

## Usage

``` r
get_selex_plot(rep, model_names, Selex_Type = "age", year_indx = NULL)
```

## Arguments

- rep:

  List of n_models of \`SPoRC\` report lists

- model_names:

  Vector of model names

- Selex_Type:

  Character vector specifying whether to output age or length-based
  selectivity (age, length)

- year_indx:

  Year index for which selectivity year to plot (can be an integer or
  vector)

## Value

Plots of terminal year fishery and survey selectivity by fleet, region,
and sex across models

## See also

Other Plotting:
[`get_biological_plot()`](https://chengmatt.github.io/SPoRC/reference/get_biological_plot.md),
[`get_data_fitted_plot()`](https://chengmatt.github.io/SPoRC/reference/get_data_fitted_plot.md),
[`get_ts_plot()`](https://chengmatt.github.io/SPoRC/reference/get_ts_plot.md),
[`plot_all_basic()`](https://chengmatt.github.io/SPoRC/reference/plot_all_basic.md),
[`theme_sablefish()`](https://chengmatt.github.io/SPoRC/reference/theme_sablefish.md)

## Examples

``` r
if (FALSE) { # \dontrun{
get_selex_plot(list(rep1, rep2), c("Model1", "Model2"), year_indx = c(1:30))
} # }
```
