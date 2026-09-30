# Pack the covariate observations that OSA residuals are available for

The observed cells of every covariate a dsem observes with error, in one
vector that the objective registers with
[`OBS`](https://rdrr.io/pkg/RTMB/man/TMB-interface.html) and
[`oneStepPredict`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html)
peels one observation at a time. Each observation keeps its own family's
density, weighted by the indicator the peel switches, so nothing is put
on a shared scale. One vector holds the families whose observations are
continuous and the other the counts, since the two need different
arguments out of `oneStepPredict`. The covariates run in the order they
were given and each one's observed years ascend, which is the order each
residual is conditioned in.

## Usage

``` r
pack_dsem_cov_osa(dsem_cov_obs, dsem_cov_family, family = "continuous")
```

## Arguments

- dsem_cov_obs:

  Matrix `[grid year, covariate]` of observations, NA in a year a
  covariate is not observed.

- dsem_cov_family:

  Family code of each covariate.

- family:

  `"continuous"` for normal, gamma, the fixed sd normal and lognormal,
  or `"bernoulli"`, `"poisson"` or `"tweedie"`.

## Value

List with `vec`, the observations, and `map`, the grid year and
covariate behind each one. `NULL` when no covariate is on a family of
the requested kind.

## Details

The counts and the tweedie each have a vector to themselves. A support
holds for every observation in the vector it is given, so a bernoulli's
outcomes of 0 and 1 and a poisson's counts cannot share one, and a
tweedie needs a range as well.
