# One fleet's fitted bins, or NULL when that fleet is fit over all of them

The same decision
[`bins_or_null`](https://chengmatt.github.io/SPoRC/dev/reference/bins_or_null.md)
makes, taken one fleet at a time and answered with the bin numbers
themselves. The fitting likelihoods subset a single fleet's observations
and need those numbers; the OSA packers walk every fleet, so they take
the whole array.

## Usage

``` r
fleet_bins_or_null(bins_arr, f)
```

## Arguments

- bins_arr:

  A `[n_obs_bins x n_fleets]` 0/1 array.

- f:

  Fleet index.

## Value

Integer vector of the bins fleet `f` is fit over, or `NULL` if it is fit
over all of them.
