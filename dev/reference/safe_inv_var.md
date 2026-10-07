# Inverse variance for Francis reweighting

Inverse of the sample variance of `x` after dropping `NA`, used to
weight each year's standardized residual. A vector with fewer than two
values, or a variance that is zero or not finite, returns `1` rather
than `NA` or `Inf`.

## Usage

``` r
safe_inv_var(x)
```

## Arguments

- x:

  Numeric vector.

## Value

Numeric scalar, `1 / var(x)` or `1` where that is undefined.
