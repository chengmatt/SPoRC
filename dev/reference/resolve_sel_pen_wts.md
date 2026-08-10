# Resolve a modular selectivity penalty weight vector

Converts a user-supplied selectivity penalty weight specification into
the complete named weight vector consumed by
[`Get_Selex_Smoothness_Penalty`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Selex_Smoothness_Penalty.md).

## Usage

``` r
resolve_sel_pen_wts(pen_wts, n_fleets = 1)
```

## Arguments

- pen_wts:

  `NULL`, or a named numeric vector/list giving independent weights for
  any subset of `"smooth_bin_curve"`, `"smooth_bin_diff"`,
  `"smooth_yr_diff"`, `"smooth_yr_curve"`, `"smooth_dome"`,
  `"smooth_mean_center"`. Any name not supplied defaults to `0`. Each
  weight is either a scalar applied to every year, or a vector with one
  value per model year, which lets a penalty act only in some years or
  act with a different strength in each. The list may also carry
  `"bin_range"`, a length-two vector giving the first and last age or
  length bin the penalties apply over; the default of every bin is what
  a missing `bin_range` reproduces.

- n_fleets:

  Integer. Number of fleets the specification covers. A single named
  specification is applied to every fleet; an unnamed list of length
  `n_fleets` gives each fleet its own.

## Value

List of length `n_fleets`. Each element is a named list of length 7: the
six penalty terms, each a scalar or a per-year vector, plus `bin_range`
(`NULL` for all bins).
