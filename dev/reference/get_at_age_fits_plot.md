# Plot observed and predicted age-disaggregated observations

Companion to
[`get_catch_fits_plot`](https://chengmatt.github.io/SPoRC/dev/reference/get_catch_fits_plot.md)
for the at-age streams: retained catch, discards, and the fishery and
survey indices at age. Each is a lognormal observation with its own
standard deviation, so intervals are built the same way, and the panel
is faceted by age.

## Usage

``` r
get_at_age_fits_plot(data, rep, model_names, stream = "CatchAA")
```

## Arguments

- data:

  List of length `n_models` of SPoRC data lists, read for the at-age
  observation and use arrays.

- rep:

  List of length `n_models` of SPoRC model reports, read for the
  predicted at-age values and their standard deviations.

- model_names:

  Character vector naming each model run.

- stream:

  Which stream to plot: `"CatchAA"` (default), `"DiscardAA"`,
  `"FishIdxAA"` or `"SrvIdxAA"`.

## Value

A `ggplot`, or `NULL` when no model fits that stream.
