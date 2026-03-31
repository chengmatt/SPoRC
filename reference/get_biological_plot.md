# Get Plots of Biological Quantities

Get Plots of Biological Quantities

## Usage

``` r
get_biological_plot(data, rep, model_names)
```

## Arguments

- data:

  List of n_models of \`SPoRC\` data lists

- rep:

  List of n_models of \`SPoRC\` report lists

- model_names:

  Vector of model names

## Value

A list of plots for terminal year movement, natural mortality,
weight-at-age, and maturity at age across models

## See also

Other Plotting:
[`get_data_fitted_plot()`](https://chengmatt.github.io/SPoRC/reference/get_data_fitted_plot.md),
[`get_selex_plot()`](https://chengmatt.github.io/SPoRC/reference/get_selex_plot.md),
[`get_ts_plot()`](https://chengmatt.github.io/SPoRC/reference/get_ts_plot.md),
[`plot_all_basic()`](https://chengmatt.github.io/SPoRC/reference/plot_all_basic.md),
[`theme_sablefish()`](https://chengmatt.github.io/SPoRC/reference/theme_sablefish.md)

## Examples

``` r
if (FALSE) { # \dontrun{
get_biological_plot(list(data1, data2), list(rep1, rep2), c("Model1", "Model2"))
} # }
```
