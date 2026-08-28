# Evaluate a robust negative binomial log-likelihood

Computes the negative binomial log-likelihood using a \\(\mu, \sigma^2 -
\mu)\\ reparameterization that remains valid for non-integer
observations via `lgamma`. The overdispersion parameter is recovered as
\\k = \mu^2 / (\sigma^2 - \mu)\\.

## Usage

``` r
dnbinom_robust_noint(x, log_mu, log_var_minus_mu, give_log = TRUE)
```

## Arguments

- x:

  Numeric. Observed count (may be non-integer).

- log_mu:

  Numeric. Log of the mean parameter \\\mu\\.

- log_var_minus_mu:

  Numeric. Log of the excess variance \\\sigma^2 - \mu\\. Must satisfy
  \\\sigma^2 \> \mu\\ (i.e., overdispersion); the implied size parameter
  is \\k = \mu^2 / (\sigma^2 - \mu)\\.

- give_log:

  Logical. If `TRUE` (default), returns the log-likelihood; otherwise
  returns the likelihood.

## Value

Numeric. Log-likelihood (or likelihood if `give_log = FALSE`).
