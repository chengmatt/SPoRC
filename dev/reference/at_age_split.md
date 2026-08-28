# Decode an at-age aggregation type into its split margins

The Type codes follow the composition vocabulary: `"agg"` sums over both
regions and sexes, `"spltRaggS"` keeps regions apart and sums over
sexes, `"aggRspltS"` does the reverse, and `"spltRspltS"` keeps both
apart.

## Usage

``` r
at_age_split(code)
```

## Arguments

- code:

  Integer, `0` to `3` in the order above.

## Value

A list with logical `region` and `sex`, `TRUE` where that margin is
split.
