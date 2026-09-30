# Covariate observation density by family and link

The observations of one covariate given its grid cells. The cell goes
through the link to the mean (identity, exp, inverse logit or inverse
cloglog), and the family's density is taken about that mean. Fixed (0)
has no density. Normal (1) has an estimated sd, gaussian_fixed_sd (5) a
known sd per observation, bernoulli (2) a coin flip at the mean, poisson
(3) that mean, gamma (4) shape \\1/CV^2\\ and that mean, lognormal (6)
that mean as its median, tweedie (7) that mean with a dispersion and a
power in (1, 2).

## Usage

``` r
get_dsem_obs_nLL(
  y,
  x,
  family,
  link,
  obs_sd,
  tweedie_p,
  fixed_sd = NULL,
  keep = 1
)
```

## Arguments

- y:

  Observed values, no NA. The objective reads these off the vector
  [`pack_dsem_cov_osa`](https://chengmatt.github.io/SPoRC/dev/reference/pack_dsem_cov_osa.md)
  builds, so that OSA residuals are available for them.

- x:

  Grid cells for those years, on the link scale.

- family:

  Family code.

- link:

  Link code.

- obs_sd:

  Measurement sd (normal), CV (gamma), sd of the log (lognormal) or
  dispersion (tweedie). Unused otherwise.

- tweedie_p:

  Tweedie power in (1, 2). Unused otherwise.

- fixed_sd:

  Known sd per observation for gaussian_fixed_sd. Unused otherwise.

- keep:

  Indicator per observation, which
  [`oneStepPredict`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html)
  switches off for every observation it has not reached yet. One during
  an ordinary fit, which leaves the likelihood as it was.

## Value

Scalar negative log likelihood.
