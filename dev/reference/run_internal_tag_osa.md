# Run internal (model-based) OSA residuals for conventional tagging data

Internal counterpart to
[`run_internal_comp_osa`](https://chengmatt.github.io/SPoRC/dev/reference/run_internal_comp_osa.md)
for conventional tag recapture data packed via
[`pack_tag_osa`](https://chengmatt.github.io/SPoRC/dev/reference/pack_tag_osa.md)
(requires `do_internal_conv_tag_osa = TRUE` in
[`Setup_Mod_Dim`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Dim.md)).
Called by
[`get_osa`](https://chengmatt.github.io/SPoRC/dev/reference/get_osa.md)
when `tag = TRUE` and a fitted `model` is supplied.

## Usage

``` r
run_internal_tag_osa(model, data, osa_method = NULL, parallel = FALSE)
```

## Arguments

- model:

  A fitted RTMB model object from
  [`fit_model`](https://chengmatt.github.io/SPoRC/dev/reference/fit_model.md),
  built with `do_internal_conv_tag_osa = TRUE`.

- data:

  The model `data` list (e.g. `input_list$data`) used to build `model`.

- osa_method:

  Optional override for
  [`RTMB::oneStepPredict`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html)'s
  `method`. Must be one of `"oneStepGeneric"`,
  `"oneStepGaussianOffMode"`, or `"oneStepGaussian"` (the `"cdf"` method
  is not permitted). Defaults to `"oneStepGeneric"` (tag recapture data
  are always discrete/count-valued).

- parallel:

  Whether or not to parallelize OSA computation. Defaults to `FALSE`.

## Value

A list with one element `res`: columns `fleet`, `region`, `pop_pool`,
`age_pool`, `sex_pool`, `cohort`, `release_year`, `release_region`,
`release_season`, `recovery_year`, `recovery_season`,
`years_at_liberty`, `is_tail`, `resid`, `family`, `comp_type = "Tag"`,
or `NULL` if no tagging data is present.
