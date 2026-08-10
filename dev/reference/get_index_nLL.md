# Evaluate an index negative log-likelihood under a chosen error structure

Calls an abundance or biomass index to use one of three error
structures. `like_type = 0` is lognormal, where `sigma` is a log-scale
standard deviation and `const` is added inside both logarithms.
`like_type = 1` is normal on the arithmetic scale, where `sigma` is an
arithmetic standard deviation. `like_type = 2` is multivariate normal on
the arithmetic scale with a fixed covariance supplied through `Sigma`,
which places the whole series in a single density and is what a
survey-provided covariance across years calls for, since a diagonal
likelihood would treat correlated residuals as independent information.

## Usage

``` r
get_index_nLL(obs, pred, sigma, like_type, Sigma = NULL, const = 0)
```

## Arguments

- obs:

  Numeric vector of observed index values.

- pred:

  Numeric vector of predicted index values.

- sigma:

  Numeric vector of standard deviations. Ignored when `like_type = 2`.

- like_type:

  Integer. `0` lognormal, `1` normal, `2` multivariate normal.

- Sigma:

  Covariance matrix. Required when `like_type = 2`. Must be symmetric
  and positive definite;
  [`dmvnorm`](https://rdrr.io/pkg/RTMB/man/MVgauss.html) reads only the
  lower triangle and returns `NaN` without warning otherwise, so this is
  checked during setup rather than here.

- const:

  Numeric constant added inside the logarithms of the lognormal form.
  Default `0`.

## Value

Numeric vector of negative log-likelihood contributions the same length
as `obs`. The multivariate normal cannot be split across observations,
so its total is returned in the first element with zeros elsewhere.
