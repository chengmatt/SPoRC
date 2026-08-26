# Combine reported index standard errors with an estimated component

An index observation carries a standard error from its own survey
design, and an assessment may additionally estimate a component covering
everything that design does not. This returns the total standard
deviation the index likelihood should use.

## Usage

``` r
combine_idx_sd(se, extra, form)
```

## Arguments

- se:

  Reported standard errors, any shape.

- extra:

  Estimated component on the natural scale, conformable with `se`.

- form:

  Integer. `0` returns `se` unchanged, `1` additive, `2` in quadrature,
  `3` replacement.

## Value

The total standard deviation, shaped like `se`.

## Details

Additive represents a parameter added to an existing estiamte of the
standard error. Quadrature treats the two as independent variances.
Replacement discards the reported errors and estimates a rate.
