# A fleet's deviations and their levels extended past the fit

The fitted years keep the deviations given and their levels. Each year
past the fit starts at zero and takes the sharing of the last fitted
year that estimates any deviation, under levels of its own, so
[`draw_sel_dev_surface`](https://chengmatt.github.io/SPoRC/dev/reference/draw_sel_dev_surface.md)
draws it fresh rather than copying a fitted year's value.

## Usage

``` r
extend_devs_past_fit(devs, map, n_fit, n_total)
```

## Arguments

- devs:

  Deviations `[region, year, bin, sex, fleet]` over at least the fitted
  years.

- map:

  Their levels, shaped as `devs`, `NA` where fixed.

- n_fit:

  Fitted years.

- n_total:

  Years the operating model runs.

## Value

List of `devs` and `map` over `n_total` years.
