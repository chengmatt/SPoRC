# Compute the corrected marginal AIC (AICc)

Calculates the corrected Akaike Information Criterion \\\mathrm{AICc} =
p k + 2 \ell + 2k(k+1)/(n-k-1)\\, where \\k\\ is the number of estimated
parameters, \\\ell\\ is the negative log-likelihood at the optimum,
\\p\\ is the penalty multiplier, and \\n\\ is the sample size. When \\n
= \infty\\ the small-sample correction term vanishes and the result
reduces to standard marginal AIC. Compatible with output from both
`optim` (`$value`) and `nlminb` (`$objective`).

## Usage

``` r
marg_AIC(opt, p = 2, n = Inf)
```

## Arguments

- opt:

  Named list of optimiser output. Must contain `$par` and either
  `$objective` (e.g., `nlminb`) or `$value` (e.g., `optim`).

- p:

  Numeric. Penalty multiplier per parameter. Default `2` (standard AIC).

- n:

  Numeric. Sample size used for the small-sample correction. Default
  `Inf` (no correction, equivalent to AIC).

## Value

Numeric. The AICc value.
