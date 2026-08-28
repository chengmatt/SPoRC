# Plotting Function for All Basic Quantities

Convenience wrapper that calls all core SPoRC plotting functions and
writes their output to a single PDF file. Equivalent to calling
`get_biological_plot`, `get_data_fitted_plot`, `get_ts_plot`,
`get_selex_plot`, and `get_nLL_plot` in sequence and printing each to
the same device.

## Usage

``` r
plot_all_basic(data, rep, sd_rep, model_names, out_path)
```

## Arguments

- data:

  List of length `n_models`, where each element is a SPoRC data list.
  Passed to `get_biological_plot`, `get_data_fitted_plot`, and
  `get_nLL_plot`.

- rep:

  List of length `n_models`, where each element is a SPoRC report list
  (i.e. the output of `obj$report()` after optimization). Passed to
  `get_biological_plot`, `get_ts_plot`, `get_selex_plot`, and
  `get_nLL_plot`.

- sd_rep:

  List of length `n_models`, where each element is a SPoRC sdreport list
  (i.e. the output of `sdreport(obj)`). Passed to `get_ts_plot` for
  delta-method confidence intervals.

- model_names:

  Character vector of length `n_models` giving display names for each
  model run. Passed to all sub-functions as legend and facet labels.

- out_path:

  Path to the output directory. The PDF is written to
  `here::here(out_path, "plot_results.pdf")` at 25 × 13 inches per page.

## Value

Called for its side effect. Writes `plot_results.pdf` to `out_path` and
returns `NULL` invisibly.

## See also

Other Plotting:
[`get_biological_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_biological_plot.md),
[`get_data_fitted_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_data_fitted_plot.md),
[`get_selex_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_selex_plot.md),
[`get_ts_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_ts_plot.md),
[`theme_sablefish()`](https://chengmatt.github.io/SPoRC/dev/reference/theme_sablefish.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  plot_all_basic(
    data        = list(data1, data2),
    rep         = list(rep1, rep2),
    sd_rep      = list(sd_rep1, sd_rep2),
    model_names = c("Base", "Sensitivity"),
    out_path    = here::here()
  )
} # }
```
