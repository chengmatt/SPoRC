# Run internal (model-based) OSA residuals for a composition data source

Internal counterpart to
[`run_external_comp_osa`](https://chengmatt.github.io/SPoRC/dev/reference/run_external_comp_osa.md)'s
external (post-hoc, compResidual-based) path, called by
[`get_osa`](https://chengmatt.github.io/SPoRC/dev/reference/get_osa.md)
when a fitted `model` is supplied. Calls
[`RTMB::oneStepPredict()`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html)
directly on the model's internally tracked OSA vector (built via
`do_internal_comp_osa = TRUE` in
[`Setup_Mod_Dim`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Dim.md)),
and relabels the resulting residuals using
[`pack_comp_osa`](https://chengmatt.github.io/SPoRC/dev/reference/pack_comp_osa.md)'s
`return_labels = TRUE` output so the result matches the same `res`
schema produced by the external path.

## Usage

``` r
run_internal_comp_osa(
  model,
  data,
  comp_source,
  family,
  pop = FALSE,
  discard = FALSE,
  parallel = FALSE,
  bins,
  bin_label,
  osa_method = NULL
)
```

## Arguments

- model:

  A fitted RTMB model object from
  [`fit_model`](https://chengmatt.github.io/SPoRC/dev/reference/fit_model.md),
  built with `do_internal_comp_osa = TRUE`.

- data:

  The model `data` list (e.g. `input_list$data`) used to build `model`.

- comp_source:

  One of `"FishAge"`, `"FishLen"`, `"SrvAge"`, `"SrvLen"`.

- family:

  Character, `"discrete"` or `"continuous"`.

- pop:

  Logical; population-specific composition source. Default `FALSE`.

- discard:

  Logical; discard composition source. Default `FALSE`.

- parallel:

  Whether or not to parallelize OSA computation. Defaults to `FALSE`.

- bins:

  Vector of age or length bin labels for display. Must span every
  observed bin of the stream, not just the ones a `*_bins` restriction
  fits: residuals are labelled by true observed bin number, so a subset
  here shifts every label.

- bin_label:

  Character label describing whether bins represent ages or lengths.

- osa_method:

  Optional override for
  [`RTMB::oneStepPredict`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html)'s
  `method`. Must be one of `"oneStepGeneric"`,
  `"oneStepGaussianOffMode"`, or `"oneStepGaussian"`; the `"cdf"` method
  is not permitted (it is numerically fragile for the discrete
  likelihoods used here). Defaults to `"oneStepGeneric"` for discrete
  data and `"oneStepGaussianOffMode"` for continuous (logistic-normal)
  data. See
  [`TMB::oneStepPredict`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html)
  for further details. Note that if data are discrete, the only valid
  option is `"oneStepGeneric"`.

## Value

A list with one element `res`, matching
[`get_osa`](https://chengmatt.github.io/SPoRC/dev/reference/get_osa.md)'s
external-mode schema (columns `fleet`, `index_label`, `year`, `index`,
`resid`, `region`, `sex`, `seas`, `comp_type`) plus a `pop` column
(population index; always 1 for `pop = FALSE` sources), or `NULL` if no
data of the requested family/source is present.
