# Evaluate a Dirichlet log-density

Computes the Dirichlet log-density \\\log \Gamma(\sum \alpha) - \sum
\log \Gamma(\alpha_k) + \sum (\alpha_k - 1) \log x_k\\ for a
compositional vector `x` and concentration parameter vector `alpha`.
Used in SPoRC as a prior for movement rates and recruitment regional
apportionment.

## Usage

``` r
ddirichlet(x, alpha, log = TRUE)
```

## Arguments

- x:

  Numeric vector of compositional values summing to 1; all elements must
  be strictly positive.

- alpha:

  Numeric vector of Dirichlet concentration parameters of the same
  length as `x`. Larger values of \\\sum \alpha\\ concentrate mass near
  \\\alpha / \sum \alpha\\.

- log:

  Logical. If `TRUE` (default), returns the log-density; otherwise
  returns the density.

## Value

Numeric. Log-density (or density if `log = FALSE`).
