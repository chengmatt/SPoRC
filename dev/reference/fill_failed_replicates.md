# Give failed replicates the shape of the ones that refit

A replicate whose refit fails is left empty or a single `NA`, and
`simplify2array` returns a list rather than an array as soon as one
replicate does. Each failed replicate becomes `NA` in the shape of a
replicate that refit, so the results keep their replicate dim and
arithmetic against the truth runs, `NA` where a refit failed.

## Usage

``` r
fill_failed_replicates(x)
```

## Arguments

- x:

  List of one result per replicate.

## Value

`x`, failed replicates filled.
