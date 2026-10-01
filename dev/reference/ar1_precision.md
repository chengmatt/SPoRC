# Precision of a stationary AR1 with unit marginal variance

Precision of a stationary AR1 with unit marginal variance

## Usage

``` r
ar1_precision(n, rho)
```

## Arguments

- n:

  Length of the series.

- rho:

  Correlation between neighbors.

## Value

Tridiagonal `n` by `n` matrix, the inverse of `rho^abs(i - j)`.
