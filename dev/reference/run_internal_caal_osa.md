# Run internal (model-based) OSA residuals for conditional age-at-length data

Internal counterpart to
[`run_internal_comp_osa`](https://chengmatt.github.io/SPoRC/dev/reference/run_internal_comp_osa.md)
for CAAL data packed via
[`pack_caal_osa`](https://chengmatt.github.io/SPoRC/dev/reference/pack_caal_osa.md)
(requires `do_internal_comp_osa = TRUE` in
[`Setup_Mod_Dim`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Dim.md)).
Called by
[`get_osa`](https://chengmatt.github.io/SPoRC/dev/reference/get_osa.md)
when `comp_source` is `"Fish_caal"` or `"Srv_caal"`. CAAL carries only
the discrete families, so there is no family argument; the residuals
come back with an extra `len` column giving the length bin each age
composition was conditioned on.

## Usage

``` r
run_internal_caal_osa(
  model,
  data,
  comp_source,
  bins,
  bin_label,
  osa_method = NULL,
  parallel = FALSE
)
```

## Arguments

- model:

  Fitted model object.

- data:

  The data list the model was built from.

- comp_source:

  Either `"Fish_caal"` or `"Srv_caal"`.

- bins:

  Age bins, used to label the residuals.

- bin_label:

  Label for the bin axis, typically `"Age"`.

- osa_method:

  Optional override for the `oneStepPredict` method.

- parallel:

  Logical, passed to `oneStepPredict`.

## Value

A list with one element `res`, matching
[`run_internal_comp_osa`](https://chengmatt.github.io/SPoRC/dev/reference/run_internal_comp_osa.md)'s
schema plus a `len` column.
