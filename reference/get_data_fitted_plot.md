# Get Data Fitted to Plot

Get Data Fitted to Plot

## Usage

``` r
get_data_fitted_plot(data, model_names)
```

## Arguments

- data:

  List of n_models of \`SPoRC\` data lists

- model_names:

  Character vector of model names

## Value

A plot of data that were fitted to across models

## See also

Other Plotting:
[`get_biological_plot()`](https://chengmatt.github.io/SPoRC/reference/get_biological_plot.md),
[`get_selex_plot()`](https://chengmatt.github.io/SPoRC/reference/get_selex_plot.md),
[`get_ts_plot()`](https://chengmatt.github.io/SPoRC/reference/get_ts_plot.md),
[`plot_all_basic()`](https://chengmatt.github.io/SPoRC/reference/plot_all_basic.md),
[`theme_sablefish()`](https://chengmatt.github.io/SPoRC/reference/theme_sablefish.md)

## Examples

``` r
if (FALSE) { # \dontrun{
get_data_fitted_plot(list(data1, data2), c("Model1", "Model2"))
} # }
```
