# Logistic Normal Likelihood

Logistic Normal Likelihood

## Usage

``` r
dlogistnormal(obs, pred, Sigma_or_Q, type = "dgmrf", give_log = TRUE)
```

## Arguments

- obs:

  Vector of observed values in numbers (can be integers or
  non-integers - vector length of n_bins)

- pred:

  Vector of predicted values in proportions (vector length of n_bins)

- Sigma_or_Q:

  Sigma or Precision Matrix

- type:

  Whether to fit using dgmrf (default; uses precision) or dmvnorm (uses
  sigma)

- give_log:

  whether or not to use log space likelihood

## Value

Returns likelihood values from a logistic normal
