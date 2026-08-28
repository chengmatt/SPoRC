# Get Plots of Biological Quantities

Generates plots of movement probabilities, natural mortality,
weight-at-age, and maturity-at-age across one or more SPoRC model runs.
All plots use a user-specified year index; movement plots additionally
show all seasons and are returned as a separate plot per population.

## Usage

``` r
get_biological_plot(data, rep, model_names)
```

## Arguments

- data:

  List of length `n_models`, where each element is a SPoRC data list.
  Used to extract `WAA`, `MatAA`, `ages`, `n_pop`, and
  `do_recruits_move`.

- rep:

  List of length `n_models`, where each element is a SPoRC report list
  (i.e. the output of `obj$report()` after optimization). Must contain
  `Movement` and `natmort`.

- model_names:

  Character vector of length `n_models` giving display names for each
  model run. Used as color legend labels across all plots.

## Value

A list of four elements:

- \[\[1\]\] move_plot:

  A list of `n_pop` `ggplot` objects, one per population. Each plot
  shows age-specific movement probabilities for all seasons at
  `year_indx`, faceted by Region_To × (Region_From + Season), with lines
  colored by model and distinguished by sex linetype. If
  `do_recruits_move = 0`, age-1 fish are excluded. Y-axis is fixed to
  \[0, 1\].

- \[\[2\]\] natmort_plot:

  Natural mortality at age at `year_indx`, faceted by Region × Sex.
  Lines colored by model and distinguished by population linetype.

- \[\[3\]\] waa_plot:

  Spawning weight-at-age at `year_indx`, faceted by Region × (Sex +
  Season). Lines colored by model and distinguished by population
  linetype.

- \[\[4\]\] mataa_plot:

  Maturity-at-age at `year_indx`, faceted by Region × (Sex + Season).
  Lines colored by model and distinguished by population linetype.

## See also

Other Plotting:
[`get_data_fitted_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_data_fitted_plot.md),
[`get_selex_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_selex_plot.md),
[`get_ts_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_ts_plot.md),
[`plot_all_basic()`](https://chengmatt.github.io/SPoRC/dev/reference/plot_all_basic.md),
[`theme_sablefish()`](https://chengmatt.github.io/SPoRC/dev/reference/theme_sablefish.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  plots <- get_biological_plot(
    data        = list(data1, data2),
    rep         = list(rep1, rep2),
    model_names = c("Base", "Sensitivity"),
    year_indx   = 40
  )

  plots[[1]][[1]]  # movement plot for population 1
  plots[[1]][[2]]  # movement plot for population 2, etc
  plots[[2]]       # natural mortality
  plots[[3]]       # weight-at-age
  plots[[4]]       # maturity-at-age
} # }
```
