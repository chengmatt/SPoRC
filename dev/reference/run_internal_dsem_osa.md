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
  parallel = FALSE,
  seed = 123
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

- parallel:

  Whether or not to parallelize OSA computation. Defaults to `FALSE`.

- seed:

  Seed for the discrete residuals' uniform draw, as
  [`get_osa`](https://chengmatt.github.io/SPoRC/dev/reference/get_osa.md)
  takes it.

## Value

A list with one element `res`: columns `covariate`, `year`, `family`,
`resid` and `idx_type` (set to `"DsemCov"`), or `NULL` when the model
has no dsem or no covariate on the requested family.
