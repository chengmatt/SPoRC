# Covariate observation density by family and link

The observations of one covariate given its grid cells. The cell goes
through the link to the mean (identity, exp, inverse logit or inverse
cloglog), and the family's density is taken about that mean. Fixed (0)
has no density. Normal (1) has an estimated sd, gaussian_fixed_sd (5) a
known sd per observation, bernoulli (2) a coin flip at the mean, poisson
(3) that mean, Gamma (4) shape \\1/CV^2\\ and that mean, lognormal (6)
that mean as its median, tweedie (7) that mean with a dispersion and a
power in (1, 2).

## Usage

``` r
get_dsem_obs_nLL(y, x, family, link, obs_sd, tweedie_p, fixed_sd = NULL)
```

## Arguments

- y:

  Observed values, no NA.

- x:

  Grid cells for those years, on the link scale.

- family:

  Family code.

- link:

  Link code.

- obs_sd:

  Measurement sd (normal), CV (Gamma), sd of the log (lognormal) or
  dispersion (tweedie). Unused otherwise.

- tweedie_p:

  Tweedie power in (1, 2). Unused otherwise.

- fixed_sd:

  Known sd per observation for gaussian_fixed_sd. Unused otherwise.

## Value

Scalar negative log likelihood.
