# Evaluate a Poisson log-likelihood for non-integer counts

Computes \\-\lambda + x \log \lambda - \log \Gamma(x+1)\\, which reduces
to the standard Poisson log-likelihood for integer `x` and extends it
continuously to non-integer values via `lgamma`.

## Usage

``` r
dpois_noint(x, pred, give_log = TRUE)
```

## Arguments

- x:

  Numeric. Observed count (may be non-integer).

- pred:

  Numeric. Predicted mean \\\lambda \> 0\\.

- give_log:

  Logical. If `TRUE` (default), returns the log-likelihood; otherwise
  returns the likelihood.

## Value

Numeric. Log-likelihood (or likelihood if `give_log = FALSE`).
