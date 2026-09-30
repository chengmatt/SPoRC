# Keep-aware multinomial log-density for OSA residuals (cdf-capable)

Conditional-binomial decomposition of an \\A\\-bin multinomial for
[`oneStepPredict`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html),
following Trijoulet et al. (2023). Each of the first \\A-1\\ bins is a
binomial conditional on the running remainder, selected by its `keep`
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

The conditional trial count is the total left after the earlier bins,
frozen at the observed counts. Given what has already been peeled it is
a constant, so it must not move with the candidate value
`oneStepPredict` sweeps over: tying it to the candidate turns the
binomial into a size-weighted negative binomial and mis-calibrates the
residuals, most visibly at a small sample size across many bins. RTMB's
own `dmultinom` OSA method and WHAM's `age_comp_osa.hpp` hold it fixed
the same way. Nothing after the peeled bin moves either, so no remaining
count can be driven negative.
