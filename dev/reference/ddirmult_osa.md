# Keep-aware Dirichlet-multinomial log-density for OSA residuals

Computes the Dirichlet-multinomial log-density of a single composition
using the conditional decomposition required by
[`oneStepPredict`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html). An
\\A\\-bin Dirichlet-multinomial is written as \\A - 1\\ conditional
two-category Dirichlet-multinomials (beta-binomials), each gated by its
`keep` element. The final bin is fixed by the sum-to-\\N\\ constraint
and contributes nothing (its residual is undefined and reported as
`NA`).

## Usage

``` r
ddirmult_osa(xobs, alpha, log = TRUE)
```

## Arguments

- xobs:

  Either an object of class `"osa"` supplied by `oneStepPredict`, or a
  plain numeric vector of observed counts (length \\A\\) during ordinary
  fitting.

- alpha:

  Dirichlet concentration parameters (length \\A\\). Typically \\\alpha
  = \hat{p} \times \exp(\ln\theta) \times N\_{total}\\; this must match
  the parameterization used by the fitting likelihood.

- log:

  Boolean on whether to return nLL

## Value

Scalar log-density contribution for the composition.

## Details

With all `keep` equal to one, the sum of the conditional log-densities
equals the joint Dirichlet-multinomial log-density. The running
remainder is frozen to the observed total (see
[`osa_extract_values`](https://chengmatt.github.io/SPoRC/dev/reference/osa_extract_values.md))
so that peeling a late bin cannot drive the remaining count negative.

Intended for use with `method = "oneStepGeneric"`.
