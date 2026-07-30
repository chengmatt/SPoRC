# Call `RTMB::oneStepPredict()` with the model's TMB DLL resolved

Wrapper used by the internal OSA runners, which works around two quirks
of [`oneStepPredict`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html):

- Its `parallel` branch calls
  [`TMB::openmp()`](https://rdrr.io/pkg/TMB/man/openmp.html) without a
  `DLL` argument, so TMB falls back to guessing the DLL and errors with
  "Multiple TMB models loaded" whenever a session has more than one TMB
  DLL loaded (e.g. RTMB alongside compResidual, which
  [`run_external_comp_osa`](https://chengmatt.github.io/SPoRC/dev/reference/run_external_comp_osa.md)
  loads). `oneStepPredict()` itself accepts no `DLL` argument, so the
  model's own DLL (`model$env$DLL`) is instead bound as the default of
  [`TMB::openmp()`](https://rdrr.io/pkg/TMB/man/openmp.html) for the
  duration of the call and restored on exit. Serial calls never reach
  that code and are forwarded untouched.

- `discreteSupport` is detected with
  [`missing()`](https://rdrr.io/r/base/missing.html), so supplying it as
  `NULL` is not the same as omitting it: a `NULL` sends continuous
  families down the mixed discrete/continuous branch, which rejects the
  Gaussian methods. It is forwarded here only when it is non-`NULL`.

## Usage

``` r
osa_one_step_predict(model, ..., discreteSupport = NULL, parallel = FALSE)
```

## Arguments

- model:

  A fitted RTMB model object from
  [`fit_model`](https://chengmatt.github.io/SPoRC/dev/reference/fit_model.md).

- ...:

  Further arguments passed to
  [`oneStepPredict`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html).

- discreteSupport:

  Support of the discrete observations, or `NULL` (the default) to omit
  the argument entirely.

- parallel:

  Whether or not to parallelize OSA computation. Defaults to `FALSE`.

## Value

The [`oneStepPredict`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html)
result.
