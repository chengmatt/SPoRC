# Evaluate a Dirichlet-multinomial log-likelihood

Computes the Dirichlet-multinomial log-likelihood following the
parameterisation of Thorson et al. (CCSRA). The concentration parameters
are \\\alpha_k = \exp(\ln\theta) \times N \times \hat{p}\_k\\, so
\\\exp(\ln\theta)\\ is the per-observation overdispersion scalar: values
near zero approach the multinomial and larger values increase variance.
Non-integer observed counts are supported via `lgamma`.

## Usage

``` r
ddirmult(obs, pred, Ntotal, ln_theta, give_log = TRUE)
```

## Arguments

- obs:

  Numeric vector of observed proportions of length \\K\\ (need not sum
  exactly to 1 after `addtotag` offsets).

- pred:

  Numeric vector of predicted proportions of length \\K\\; must sum to
  1.

- Ntotal:

  Numeric. Total count (input sample size \\N\\).

- ln_theta:

  Numeric. Log overdispersion parameter. The Dirichlet concentration is
  \\\exp(\ln\theta) \times N \times \hat{p}\_k\\.

- give_log:

  Logical. If `TRUE` (default), returns the log-likelihood; otherwise
  returns the likelihood.

## Value

Numeric. Log-likelihood (or likelihood if `give_log = FALSE`).
