# Keep-aware multinomial log-density for OSA residuals

Computes the multinomial log-density of a single composition using the
conditional decomposition required by
[`oneStepPredict`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html).
Following Trijoulet et al. (2023) and WHAM (`src/age_comp_osa.hpp`), an
\\A\\-bin multinomial is written as \\A - 1\\ conditional two-category
multinomials, each gated by its `keep` element. The final bin is fixed
by the sum-to-\\N\\ constraint and contributes nothing (its residual is
undefined and reported as `NA`).

## Usage

``` r
dmultinom_osa(xobs, p, log = TRUE)
```

## Arguments

- xobs:

  Either an object of class `"osa"` supplied by `oneStepPredict`, or a
  plain numeric vector of observed counts (length \\A\\) during ordinary
  fitting.

- p:

  Predicted proportions (length \\A\\); normalized internally.

- log:

  Boolean on whether to return nLL

## Value

Scalar log-density contribution for the composition.

## Details

With all `keep` equal to one, the sum of the conditional log-densities
equals the joint multinomial log-density; enabling OSA therefore changes
how the likelihood is decomposed, not its value. The running remainder
is frozen to the observed total (see
[`osa_extract_values`](https://chengmatt.github.io/SPoRC/dev/reference/osa_extract_values.md))
so that peeling a late bin cannot drive the remaining count negative.

Intended for use with `method = "oneStepGeneric"`; the two-category
conditionals carry no analytic CDF hooks, so the `"cdf"` method is not
supported (consistent with WHAM, which omits it for compositions).
