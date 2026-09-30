# Call `RTMB::oneStepPredict()` with the model's TMB DLL resolved

Used by the internal OSA routines, and works around two quirks of
[`oneStepPredict`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html):

- Its `parallel` branch calls
  [`TMB::openmp()`](https://rdrr.io/pkg/TMB/man/openmp.html) without a
  `DLL` argument, so TMB falls back to guessing the DLL and errors with
  "Multiple TMB models loaded" whenever a session has more than one TMB
  DLL loaded (e.g. RTMB alongside compResidual, which
  [`run_external_comp_osa`](https://chengmatt.github.io/SPoRC/dev/reference/run_external_comp_osa.md)
  loads).

- `discreteSupport` and `range` are detected with
  [`missing()`](https://rdrr.io/r/base/missing.html), so giving either a
  `NULL` is not the same as omitting it: a `NULL` support puts a
  continuous family on the mixed discrete/continuous path, which then
  errors under every Gaussian method and demands a `range` under
  `oneStepGeneric`. Both are forwarded here only when they are
  non-`NULL`.

## Usage

``` r
osa_one_step_predict(
  model,
  ...,
  discreteSupport = NULL,
  range = NULL,
  parallel = FALSE
)
```

## Arguments

- model:

  A fitted RTMB model object from
  [`fit_model`](https://chengmatt.github.io/SPoRC/dev/reference/fit_model.md).

- ...:

  Further arguments passed to
  [`oneStepPredict`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html).

- discreteSupport:

  Values a discrete observation can take, or `NULL` (the default) to
  omit the argument.

- range:

  Interval a part discrete observation's continuous part is integrated
  over (a tweedie's is `c(0, Inf)`), or `NULL` (the default) to omit it.

- parallel:

  Whether or not to parallelize OSA computation. Defaults to `FALSE`.

## Value

The [`oneStepPredict`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html)
result.
