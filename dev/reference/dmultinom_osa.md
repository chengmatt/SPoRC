# Keep-aware multinomial log-density for OSA residuals (cdf-capable)

Conditional-binomial decomposition of an \\A\\-bin multinomial for
[`oneStepPredict`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html),
following Trijoulet et al. (2023). Each of the first \\A-1\\ bins is a
binomial conditional on the running remainder, gated by its `keep`
element; the final bin is fixed by the sum-to-\\N\\ constraint. In
addition to the density term, the analytic conditional binomial CDF is
accumulated through the `cdf_lower` / `cdf_upper` indicators, so this
density supports **both** `method = "cdf"` and
`method = "oneStepGeneric"`.

## Usage

``` r
dmultinom_osa(xobs, p, log = TRUE)
```

## Arguments

- xobs:

  An `"osa"` object from `oneStepPredict`, or a plain numeric count
  vector (length \\A\\) during fitting.

- p:

  Predicted proportions (length \\A\\); normalized internally.

- log:

  Logical; return the log-density (default) or the density.

## Value

Scalar (log-)density contribution.

## Details

The running remainder is frozen to the observed total so peeling a late
bin cannot drive the remaining count negative.
