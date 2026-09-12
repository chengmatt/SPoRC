# An unused at-age input moved onto the observed ages

Fits made before the at-age data sources sat on the observed ages hold
their unused ones on the model ages. With nothing observed there is
nothing to lose, so a placeholder on the observed ages stands in; a data
source the fit observes is returned as it is, for the shape checks to
judge.

## Usage

``` r
unused_at_age_on_obs_ages(x, used, age_dim, n_obs_ages, fill = 0)
```

## Arguments

- x:

  The array, or `NULL`.

- used:

  Logical. `TRUE` when the fitted model observes this data source.

- age_dim:

  Integer. Position of the age dim in `x`.

- n_obs_ages:

  Integer. Number of observed ages the operating model draws on.

- fill:

  Value the placeholder holds.

## Value

`x`, or an array of `fill` shaped like `x` but on `n_obs_ages` ages.
