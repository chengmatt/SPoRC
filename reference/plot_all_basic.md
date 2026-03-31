# Plotting function for all basic quantities

Plotting function for all basic quantities

## Usage

``` r
plot_all_basic(data, rep, sd_rep, model_names, out_path)
```

## Arguments

- data:

  List of n_models of \`SPoRC\` data lists

- rep:

  List of n_models of \`SPoRC\` report lists

- sd_rep:

  List of n_models of sd report lists from \`SPoRC\`

- model_names:

  Character vector of model names

- out_path:

  Path to the output directory. Users only need to specify the path.

## Value

A series of plots compared across models outputted as a pdf in the
specified directory

## See also

Other Plotting:
[`get_biological_plot()`](https://chengmatt.github.io/SPoRC/reference/get_biological_plot.md),
[`get_data_fitted_plot()`](https://chengmatt.github.io/SPoRC/reference/get_data_fitted_plot.md),
[`get_selex_plot()`](https://chengmatt.github.io/SPoRC/reference/get_selex_plot.md),
[`get_ts_plot()`](https://chengmatt.github.io/SPoRC/reference/get_ts_plot.md),
[`theme_sablefish()`](https://chengmatt.github.io/SPoRC/reference/theme_sablefish.md)

## Examples

``` r
if (FALSE) { # \dontrun{
plot_all_basic(
  data = list(data1, data2),
  rep = list(rep1, rep2),
  sd_rep = list(sd_rep1, sd_rep2),
  model_names = c("Model1", "Model2"),
  out_path = here::here()
)
} # }
```
