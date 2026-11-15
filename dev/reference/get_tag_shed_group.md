# Tag shedding groups

Cohorts mapped to the same estimated shedding parameter, or fixed at the
same shedding value, are in one group, numbered by the first cohort in
it. Under `move_timing = 2` the cohorts of a group share one seasonal
operator.

## Usage

``` r
get_tag_shed_group(conv_tag_shed_map, ln_conv_tag_shed)
```

## Arguments

- conv_tag_shed_map:

  Numeric vector `[cohort]`, the map of `ln_conv_tag_shed` (`NA` where
  it is fixed), or `NULL`, which makes every cohort its own group.

- ln_conv_tag_shed:

  Vector `[cohort]` of log shedding rates.

## Value

Integer vector `[cohort]` of group numbers.
