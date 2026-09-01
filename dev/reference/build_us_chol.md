# Cholesky factor of an unstructured correlation matrix

The factor behind
[`build_us_corr`](https://chengmatt.github.io/SPoRC/dev/reference/build_us_corr.md),
returned in its own right for callers that need to whiten rather than to
form the correlation. Whitening a margin, \\z = L^{-1} x\\, is what lets
a correlation over one margin compose with an arbitrary structure over
the others: the transformed margin is independent, so each slice can
then be scored by whatever density the remaining margins call for,
separable or not.

## Usage

``` r
build_us_chol(pars, n)
```

## Arguments

- pars:

  Numeric vector of length \\n(n-1)/2\\, unconstrained.

- n:

  Dimension of the margin.

## Value

An \\n imes n\\ lower triangular matrix with unit-length rows.

## Details

Rows are normalized to unit length, which makes the lower triangular
matrix a Cholesky factor of a correlation matrix for any parameter
values, so no positive-definiteness constraint is needed on the
parameters themselves.
