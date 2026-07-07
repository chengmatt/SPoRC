# Resolve a modular selectivity penalty weight vector

Converts a user-supplied selectivity penalty weight specification into
the complete named weight vector consumed by
[`Get_Selex_Smoothness_Penalty`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Selex_Smoothness_Penalty.md).

## Usage

``` r
resolve_sel_pen_wts(pen_wts)
```

## Arguments

- pen_wts:

  `NULL`, or a named numeric vector/list giving independent weights for
  any subset of `"smooth_bin_curve"`, `"smooth_bin_diff"`,
  `"smooth_yr_diff"`, `"smooth_yr_curve"`, `"smooth_dome"`,
  `"smooth_mean_center"`. Any name not supplied defaults to `0`.

## Value

Named numeric vector of length 6:
`c(smooth_bin_curve, smooth_bin_diff, smooth_yr_diff, smooth_yr_curve, smooth_dome, smooth_mean_center)`.
