# Keep-aware Dirichlet-multinomial log-density for OSA residuals (cdf-capable)

Conditional beta-binomial decomposition of an \\A\\-bin
Dirichlet-multinomial for
[`oneStepPredict`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html),
following Trijoulet et al. (2023). Each of the first \\A-1\\ bins is a
beta-binomial conditional on the running remainder, gated by its `keep`
element; the analytic conditional beta-binomial CDF is accumulated
through `cdf_lower` / `cdf_upper`, so this density supports **both**
`method = "cdf"` and `method = "oneStepGeneric"`.

## Usage

``` r
ddirmult_osa(xobs, alpha, log = TRUE)
```

## Arguments

- xobs:

  An `"osa"` object from `oneStepPredict`, or a plain numeric count
  vector (length \\A\\) during fitting.

- alpha:

  Concentration parameters (length \\A\\). Typically \\\alpha =
  \hat{p}\\\exp(\ln\theta)\\N\_{total}\\; must match the fitting
  parameterization.

- log:

  Logical; return the log-density (default) or the density.

## Value

Scalar (log-)density contribution.
