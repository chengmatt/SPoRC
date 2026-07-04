# Extend an array along its year dimension

Appends additional year slices to an array along a specified dimension,
using one of several fill strategies. Used in SPoRC to extend
biological, selectivity, and sample-size arrays from the conditioning
period into projection years before running closed-loop MSE simulations.

## Usage

``` r
extend_years(arr, n_years, yr_dim, fill = "zeros")
```

## Arguments

- arr:

  Array of any dimensionality to extend.

- n_years:

  Integer. Total number of years in the extended array (i.e., the new
  size of dimension `yr_dim`). Must be greater than `dim(arr)[yr_dim]`.

- yr_dim:

  Integer. Index of the dimension in `arr` corresponding to years.

- fill:

  Character string or numeric. Fill strategy for the appended year
  slices:

  `"zeros"`

  :   Fill with zeros.

  `"last"`

  :   Repeat the last year slice that contains at least one non-`NA`,
      non-`NaN` value. If no valid slice exists, fills with `NA`.

  `"mean"`

  :   Fill with the per-element mean across years, excluding zeros,
      `NA`, and `NaN` values. Elements with no valid values are set to
      zero.

  `"F_pattern"`

  :   Fill with zeros; signals to downstream functions (e.g.,
      [`predict_sim_fish_iss_fmort`](https://chengmatt.github.io/SPoRC/dev/reference/predict_sim_fish_iss_fmort.md))
      that sample sizes should be dynamically updated based on projected
      fishing mortality during the closed-loop simulation.

  Numeric scalar or array

  :   Fill all appended slices with the supplied constant value,
      recycled via [`array()`](https://rdrr.io/r/base/array.html).

## Value

Array with the same dimensions as `arr` except that
`dim(result)[yr_dim] == dim(arr)[yr_dim] + n_years`, formed by binding
`arr` and the fill array along `yr_dim` via
[`abind::abind`](https://rdrr.io/pkg/abind/man/abind.html).
