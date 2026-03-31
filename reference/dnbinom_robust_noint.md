# Negative Binomial Likelihood

Approximates negative binomial when non-integer values are provided

## Usage

``` r
dnbinom_robust_noint(x, log_mu, log_var_minus_mu, give_log = TRUE)
```

## Arguments

- x:

  observations

- log_mu:

  log mu

- log_var_minus_mu:

  log var minus mu - reparameterize negbin

- give_log:

  whether to give log

## Value

Returns likelihood values from a robust negative binomial
