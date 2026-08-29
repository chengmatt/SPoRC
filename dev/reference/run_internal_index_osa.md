# Run internal (model-based) OSA residuals for an index-type data source

Internal counterpart to
[`run_internal_comp_osa`](https://chengmatt.github.io/SPoRC/dev/reference/run_internal_comp_osa.md)
for catch/discard/index data (`ObsCatch`, `ObsDiscard`, `ObsFishIdx`,
`ObsSrvIdx`, and their `_pop` variants). These are always continuous
observations.

## Usage

``` r
run_internal_index_osa(
  model,
  data,
  index_source,
  pop = FALSE,
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

- index_source:

  One of `"Catch"`, `"Discard"`, `"FishIdx"`, `"SrvIdx"`, or their
  at-age forms `"CatchAA"`, `"DiscardAA"`, `"SrvIdxAA"`. At-age sources
  return extra `age` and `sex` columns.

- pop:

  Logical; population-specific index source. Default `FALSE`.

- osa_method:

  Optional override for
  [`RTMB::oneStepPredict`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html)'s
  `method`. Must be one of `"oneStepGeneric"`,
  `"oneStepGaussianOffMode"`, or `"oneStepGaussian"` (the `"cdf"` method
  is not permitted). Defaults to `"oneStepGeneric"`.

- parallel:

  Whether or not to parallelize OSA computation. Defaults to `FALSE`.

## Value

A list with one element `res`: columns `fleet`, `region`, `year`,
`season`, `pop`, `age`, `sex`, `resid`, and `idx_type` (set to
`index_source`; named `idx_type` rather than `comp_type` because
index-type sources are not compositions), or `NULL` if no data of the
requested source is present. `age` and `sex` are `NA` for the aggregated
sources.
