# Run internal (model-based) OSA residuals for a dsem's covariate observations

Internal counterpart to
[`run_internal_index_osa`](https://chengmatt.github.io/SPoRC/dev/reference/run_internal_index_osa.md)
for the covariates a dsem observes with error, packed into a tracked
vector by
[`pack_dsem_cov_osa`](https://chengmatt.github.io/SPoRC/dev/reference/pack_dsem_cov_osa.md).
Every family has residuals but `"fixed"`, where the covariate is the
grid cell itself and there is no observation to peel.

## Usage

``` r
run_internal_dsem_osa(
  model,
  data,
  family = "continuous",
  osa_method = NULL,
  parallel = FALSE
)
```

## Arguments

- model:

  A fitted RTMB model object from
  [`fit_model`](https://chengmatt.github.io/SPoRC/dev/reference/fit_model.md).

- data:

  The model `data` list (e.g. `input_list$data`) used to build `model`.

- family:

  `"continuous"` for the normal, gamma, fixed sd normal and lognormal
  covariates, or `"bernoulli"`, `"poisson"` or `"tweedie"`.

- osa_method:

  Optional override for
  [`RTMB::oneStepPredict`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html)'s
  `method`. `"oneStepGaussian"` is the default when every covariate in
  the call is `"normal"` or `"gaussian_fixed_sd"` on the identity link,
  where the observation's density is normal in the value tracked, so it
  is exact and about five times quicker. Every other case defaults to
  `"oneStepGeneric"`, which integrates each observation out, and refuses
  a Gaussian method: on a skewed observation the curvature it reads as a
  variance is not one, and the residuals come back too large in the tail
  however right the model is.

- parallel:

  Whether or not to parallelize OSA computation. Defaults to `FALSE`.

## Value

A list with one element `res`: columns `covariate`, `year`, `family`,
`resid` and `idx_type` (set to `"DsemCov"`), or `NULL` when the model
has no dsem or no covariate on the requested family.

## Details

Each kind of observation takes its own arguments out of
[`oneStepPredict`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html), and
those arguments hold for every observation in the vector they are given,
which is why the families are packed apart. A continuous covariate takes
a range wide enough that no residual is lost to the integrator. A
bernoulli's outcomes are 0 and 1 and a poisson's counts run past the
largest one observed, so neither support can cover the other. A
tweedie's zeros are named as a point mass with a range for the positive
part, which puts it on the mixed discrete/continuous path.
