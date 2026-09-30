# Pack the covariate observations that OSA residuals are available for

The observed cells of every covariate a dsem observes with error, in one
vector per kind that the objective registers with
[`OBS`](https://rdrr.io/pkg/RTMB/man/TMB-interface.html). Each keeps its
own family's density, weighted by the peel's indicator. The kinds are
packed apart because `oneStepPredict`'s support and range apply to every
observation in a call.

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
