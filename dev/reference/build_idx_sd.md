# Build the total index standard deviation from a per-fleet estimated component

Applies
[`combine_idx_sd`](https://chengmatt.github.io/SPoRC/dev/reference/combine_idx_sd.md)
fleet by fleet, where the fleet is the last dimension of `se`. Handles
both the aggregated arrays, indexed region by year by season by fleet,
and the population-specific arrays, which have a leading population
dimension.

## Usage

``` r
build_idx_sd(se, ln_sigma, form)
```

## Arguments

- se:

  Reported standard errors, with fleet as the last dimension.

- ln_sigma:

  Log-scale estimated component, one value per fleet.

- form:

  Integer form code, see
  [`combine_idx_sd`](https://chengmatt.github.io/SPoRC/dev/reference/combine_idx_sd.md).

## Value

An array shaped like `se`.
