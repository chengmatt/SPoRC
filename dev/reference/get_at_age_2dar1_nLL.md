# Evaluate an at-age observation block correlated over both age and year

A separable first-order autoregression treats the residual surface as an
AR(1) across ages and an AR(1) across years, with the covariance the
Kronecker product of the two. It is defined over a complete grid, so the
caller supplies a rectangular block and the whole block's density is
returned as one number.

## Usage

``` r
get_at_age_2dar1_nLL(resid, sigma, trans_rho_age, trans_rho_year)
```

## Arguments

- resid:

  Matrix of residuals, years by ages.

- sigma:

  Matrix of standard deviations, shaped like `resid`.

- trans_rho_age:

  Unconstrained correlation between adjacent ages.

- trans_rho_year:

  Unconstrained correlation between adjacent years.

## Value

The negative log likelihood of the whole block, a scalar.
