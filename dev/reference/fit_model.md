# Fit a SPoRC RTMB model

Constructs an RTMB automatic differentiation function via
`RTMB::MakeADFun`, optimizes it with
[`stats::nlminb`](https://rdrr.io/r/stats/nlminb.html), and optionally
refines the solution with Newton steps using the analytic Hessian. The
best parameter vector (`obj$env$last.par.best`), optimizer output, and
model report are attached to the returned object.

## Usage

``` r
fit_model(
  data,
  parameters,
  mapping,
  random = NULL,
  newton_loops = 3,
  silent = FALSE,
  do_optim = TRUE,
  nlminb_control = list(iter.max = 1e+05, eval.max = 1e+05, rel.tol = 1e-15),
  lower = NULL,
  upper = NULL,
  model = SPoRC_rtmb,
  ...
)
```

## Arguments

- data:

  Named list of model data as constructed by the `Setup_Mod_*` family of
  functions.

- parameters:

  Named list of parameter starting values.

- mapping:

  Named list of factor maps controlling parameter sharing and fixing.

- random:

  Character vector of parameter names to integrate out as random
  effects. `NULL` (default) fits a fixed-effects-only model.

- newton_loops:

  Integer. Number of Newton refinement steps applied after `nlminb`
  convergence to reduce gradient magnitudes. Each step solves
  \\\Delta\theta = -H^{-1} g\\ and updates the objective. Default `3`.
  Errors and warnings are caught silently via `tryCatch`, so a step that
  fails leaves the `nlminb` solution in place without a message. \\H\\
  comes from the AD tape (`obj$he`) for fixed-effects models, which is
  exact and costs a single call. Random-effects models fall back to
  [`optimHess`](https://rdrr.io/r/stats/optim.html) differencing the
  gradient, since RTMB does not implement a tape Hessian when random
  effects are present. Refinement stops early if \\H\\ comes back
  non-finite, which happens on models that have not converged, where
  second derivatives can be undefined at parameter values the objective
  and gradient still evaluate at. The `nlminb` solution is kept in that
  case.

- silent:

  Logical. If `TRUE`, suppresses RTMB and optimizer console output.
  Default `FALSE`.

- do_optim:

  Logical. If `TRUE` (default), runs `nlminb` and Newton refinement. If
  `FALSE`, returns the un-optimized `MakeADFun` object only.

- nlminb_control:

  Named list of control parameters passed to
  [`stats::nlminb`](https://rdrr.io/r/stats/nlminb.html). Default
  `list(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-15)`.

- lower:

  Numeric vector of lower bounds for `obj$par` (the estimated parameter
  vector, i.e. after mapping and random-effects marginalization), passed
  to [`stats::nlminb`](https://rdrr.io/r/stats/nlminb.html) and used to
  clamp each Newton refinement step. `NULL` (default) is unbounded
  (`-Inf` for every element).

- upper:

  Numeric vector of upper bounds for `obj$par`, passed to
  [`stats::nlminb`](https://rdrr.io/r/stats/nlminb.html) and used to
  clamp each Newton refinement step. `NULL` (default) is unbounded
  (`Inf` for every element).

- model:

  Function with signature `function(pars, data)` passed to
  `RTMB::MakeADFun` via
  [`cmb`](https://chengmatt.github.io/SPoRC/dev/reference/cmb.md).
  Default
  [`SPoRC_rtmb`](https://chengmatt.github.io/SPoRC/dev/reference/SPoRC_rtmb.md).
  Allows non-SPoRC RTMB models to be fit with the same optimization and
  Newton-refinement machinery.

- ...:

  Additional arguments forwarded to `RTMB::MakeADFun`.

## Value

The RTMB `ADFun` object with additional fields: `$optim` (the `nlminb`
output list, with `$lower`/`$upper` recording the bounds used), `$rep`
(the model report evaluated at `obj$env$last.par.best`), and `$data`,
`$parameters`, `$mapping`, `$random`.

## See also

Other Utility:
[`get_logistN_Sigma()`](https://chengmatt.github.io/SPoRC/dev/reference/get_logistN_Sigma.md),
[`get_par_est_info()`](https://chengmatt.github.io/SPoRC/dev/reference/get_par_est_info.md),
[`post_optim_sanity_checks()`](https://chengmatt.github.io/SPoRC/dev/reference/post_optim_sanity_checks.md),
[`rho_trans()`](https://chengmatt.github.io/SPoRC/dev/reference/rho_trans.md),
[`set_data_indicator_unused()`](https://chengmatt.github.io/SPoRC/dev/reference/set_data_indicator_unused.md)

## Examples

``` r
if (FALSE) { # \dontrun{
obj <- fit_model(data, parameters, mapping, random = NULL, newton_loops = 3)
obj$rep$SSB
} # }
```
