# Evaluate a logistic-normal log-likelihood

Applies the additive log-ratio (ALR) transformation to observed and
predicted compositions (removing the last reference bin) and evaluates a
multivariate normal log-density
([`RTMB::dmvnorm`](https://rdrr.io/pkg/RTMB/man/MVgauss.html), using a
covariance matrix). The ALR mean vector is \\\mu_k = \log(\hat{p}\_k /
\hat{p}\_K)\\, \\k = 1,\ldots,K-1\\.

## Usage

``` r
dlogistnormal(obs, pred, Sigma, give_log = TRUE)
```

## Arguments

- obs:

  Numeric vector of observed composition values of length \\K\\. Need
  not sum to 1; the last element is the ALR reference.

- pred:

  Numeric vector of predicted proportions of length \\K\\; the last
  element is the ALR reference.

- Sigma:

  Covariance matrix \\\Sigma\\

- give_log:

  Logical. If `TRUE` (default), returns the log-likelihood; otherwise
  returns the likelihood.

## Value

Numeric. Log-likelihood (or likelihood if `give_log = FALSE`).
