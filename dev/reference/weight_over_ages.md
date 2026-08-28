# Weight an at-age likelihood component down to the weight's own dimensions

An at-age component is reported over age and sex, and its weight is not,
so both margins are summed within a cell before the weight is applied.
This is what the objective does when it forms `jnLL`, and the profile
has to match it or the components it reports will not add up to the
total.

## Usage

``` r
weight_over_ages(component, weight)
```

## Arguments

- component:

  An at-age likelihood array, or `NULL`.

- weight:

  The stream's weight.

## Value

The weighted component with age and sex summed out, or `NULL`.
