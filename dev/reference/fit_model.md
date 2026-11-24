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
  newton_hessian = "exact",
  silent = FALSE,
  do_optim = TRUE,
  nlminb_control = list(iter.max = 1e+05, eval.max = 1e+05, rel.tol = 1e-12),
  lower = NULL,
  upper = NULL,
  model = SPoRC_rtmb,
  do_sdrep = FALSE,
  post_optim_sanity_checks_control = list(gradient_tol = 0.001, se_tol = 100, corr_tol =
    0.99),
  tape_config = NULL,
  tmb_config = NULL,
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
  \\\Delta\theta = -H^{-1} g\\ and updates the objective. Default `3`. A
  step that would raise the largest absolute gradient is not taken and
  refinement stops there, as happens when the Hessian is not positive
  definite. Errors and warnings are caught silently via `tryCatch`, so a
  step that fails leaves the `nlminb` solution in place without a
  message. \\H\\ comes from the AD tape (`obj$he`) for fixed-effects
  models, which is exact and costs a single call. Random-effects models
  fall back to [`optimHess`](https://rdrr.io/r/stats/optim.html)
  differencing the gradient, since RTMB does not implement a tape
  Hessian when random effects are present. Refinement stops early if
  \\H\\ comes back non-finite, which happens on models that have not
  converged, where second derivatives can be undefined at parameter
  values the objective and gradient still evaluate at. The `nlminb`
  solution is kept in that case.

- newton_hessian:

  Character. Where each refinement step's Hessian comes from. `"exact"`
  (default) is the behavior above. `"diff"` differences the exact
  gradient ([`optimHess`](https://rdrr.io/r/stats/optim.html), steps of
  1e-5), for models whose tape Hessian is too large to build, and in
  place of the Newton step reruns `nlminb` with each parameter in its
  own standard errors (the parameters scaled by the Cholesky factor of
  the inverse Hessian). A differenced Hessian is not accurate enough for
  a raw Newton step when the curvature spans many orders of magnitude,
  but it makes the surface round for `nlminb`, which then converges in a
  few iterations. The last Hessian is returned as `$hessian`, for
  `RTMB::sdreport(obj, hessian.fixed = obj$hessian)`.

- silent:

  Logical. If `TRUE`, suppresses RTMB and optimizer console output.
  Default `FALSE`.

- do_optim:

  Logical. If `TRUE` (default), runs `nlminb` and Newton refinement. If
  `FALSE`, returns the un-optimized `MakeADFun` object only.

- nlminb_control:

  Named list of control parameters passed to
  [`stats::nlminb`](https://rdrr.io/r/stats/nlminb.html). Default
  `list(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-12)`.

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
  Newton-refinement routines. A `data` list without the dims
  [`Setup_Mod_Dim`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Dim.md)
  writes (`years`, `n_regions`, `ages`, `lens`, `n_sexes`,
  `n_fish_fleets` and `n_srv_fleets`) is taken as such a model, and the
  SPoRC setup steps are skipped: filling in movement and tagging when
  they were not set up, and refreshing the copies of the map kept in the
  data list. A SPoRC data list with a model of your own, such as an
  edited copy of `SPoRC_rtmb`, still takes them.

- do_sdrep:

  Logical. If `TRUE`, runs `RTMB::sdreport(obj)` at the fitted
  parameters, returns it as `$sd_rep`, and runs
  [`post_optim_sanity_checks`](https://chengmatt.github.io/SPoRC/dev/reference/post_optim_sanity_checks.md)
  on it. Default `FALSE`, which skips both. `sdreport` differences the
  gradient for its Hessian, so on an ill-conditioned fit it can report a
  Hessian that is not positive definite where the tape Hessian is; rerun
  it as `RTMB::sdreport(obj, hessian.fixed = obj$he())` before reading
  that check as a real failure.

- post_optim_sanity_checks_control:

  Named list of the tolerances passed to
  [`post_optim_sanity_checks`](https://chengmatt.github.io/SPoRC/dev/reference/post_optim_sanity_checks.md)
  under `do_sdrep = TRUE`: `gradient_tol` (largest absolute gradient),
  `se_tol` (largest standard error) and `corr_tol` (largest absolute
  correlation between two parameters). Default
  `list(gradient_tol = 1e-3, se_tol = 100, corr_tol = 0.99)`. A list
  naming only some of them, such as `list(gradient_tol = 1e-4)`, takes
  the defaults for the rest. Any other element, named or not, stops the
  fit with an error before it starts. The checks print a message for
  each one that fails, or a single message when all pass, whatever
  `silent` is set to.

- tape_config:

  Named list of how RTMB records the model, passed to
  [`RTMB::TapeConfig()`](https://rdrr.io/pkg/RTMB/man/Tape.html) before
  it is built. `NULL` (default) leaves the session's settings as they
  are. Elements, each a single string:

  `vectorize`

  :   `"enable"` records arithmetic, math functions, sums and sparse
      matrix products as one operation per vector rather than one per
      element. `"disable"` is RTMB's default.

  `atomic`

  :   `"enable"` records each matrix product as one operation.
      `"disable"` records matrix products and multivariate normal
      densities element by element; RTMB's default records the densities
      as one operation and the products element by element, so
      `"disable"` also changes the densities.

  `comparison`

  :   What happens when the model compares two quantities it is
      differentiating: `"forbid"` (RTMB's default) stops with an error,
      `"tape"` records the comparison as an indicator, and `"allow"`,
      for testing only, decides it once at the starting values.

  `"NA"` leaves an element as set. Each recording gives the same
  objective and gradient, to rounding. A vectorized or atomic operation
  is taken as making every output depend on every input, so with
  `random` set, `vectorize = "enable"` makes the Hessian of the random
  effects dense, and `atomic = "enable"` does too when one product spans
  many random effects; each Laplace gradient is then several times
  slower. Time a fit both ways before keeping either.

- tmb_config:

  Named list of TMB settings for the models RTMB builds, passed to
  `TMB::config(..., DLL = "RTMB")` before the model is built. `NULL`
  (default) leaves them as they are. Any name
  `TMB::config(DLL = "RTMB")` returns is allowed, each a single 0 or 1
  (or a count for `nthreads`). The ones that change a fit:

  `tmbad.sparse_hessian_compress`

  :   `1` compresses the sparse Hessian of the random effects, which
      takes less memory and can be slower. Read only when `random` is
      set. Default `0`.

  `tmbad.atomic_sparse_log_determinant`

  :   `0` computes the Laplace log determinant from the plain
      factorization rather than TMB's single operation, which has no
      second derivative.
      [`Setup_Mod_DSEM`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_DSEM.md)
      sets `0` itself when a moderated arrow needs it, and setting `1`
      then gives NaN. Default `1`.

  `optimize.instantly`

  :   `0` skips removing repeated and unused operations from the
      recorded model, which is slower but keeps the operations in the
      order they were recorded, so each can be traced to the line that
      recorded it. Default `1`.

  Both lists are in force for the builds, Newton steps and `sdreport`
  inside `fit_model`, and the earlier settings are restored when it
  returns or stops. The diagnostics that refit the model
  ([`do_retrospective`](https://chengmatt.github.io/SPoRC/dev/reference/do_retrospective.md),
  [`do_jitter`](https://chengmatt.github.io/SPoRC/dev/reference/do_jitter.md),
  [`do_likelihood_profile`](https://chengmatt.github.io/SPoRC/dev/reference/do_likelihood_profile.md))
  do not record / save these settings

- ...:

  Additional arguments forwarded to `RTMB::MakeADFun`.

## Value

The RTMB `ADFun` object, a list of class `sporc_fit` with `print`,
`summary`, `logLik` and `plot` methods (see
[`sporc_fit`](https://chengmatt.github.io/SPoRC/dev/reference/sporc_fit.md)),
with additional fields: `$optim` (the `nlminb` output list, with
`$lower`/`$upper` recording the bounds used and `$gradient` the gradient
at the estimates), `$rep` (the model report evaluated at
`obj$env$last.par.best`; with `random` set, `jnLL` is at the random
effect mode and `$rep` also has `marginal_nLL`, the marginal the fit
minimized, and `laplace_nLL`, the marginal minus `jnLL`), `$data`,
`$parameters`, `$mapping`, `$random`, `$provenance` (the package version
and commit that built the fit, see
[`fit_provenance`](https://chengmatt.github.io/SPoRC/dev/reference/fit_provenance.md)),
`$hessian` (the last refinement step's Hessian; `NULL` if none was
taken), and under `do_sdrep = TRUE` `$sd_rep` (the `sdreport`) and
`$post_optim_sanity_checks` (`TRUE` if every check passed, `FALSE`
otherwise). Both are `NULL` under `do_sdrep = FALSE`.

## Details

With a projection set up by
[`Setup_Mod_Projection`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Projection.md),
the model is fitted with the projection off, since it adds nothing to
the likelihood, then rebuilt at the estimates with the projection on.
RTMB records a model once, at its starting values, and a rule or cap
that compares two numbers keeps whichever side held then, so the rebuild
is what puts the estimated stock on the right side. `do_sdrep = TRUE`
then covers the projected quantities.

## See also

Other Utility:
[`get_logistN_Sigma()`](https://chengmatt.github.io/SPoRC/dev/reference/get_logistN_Sigma.md),
[`get_par_est_info()`](https://chengmatt.github.io/SPoRC/dev/reference/get_par_est_info.md),
[`label_dims()`](https://chengmatt.github.io/SPoRC/dev/reference/label_dims.md),
[`post_optim_sanity_checks()`](https://chengmatt.github.io/SPoRC/dev/reference/post_optim_sanity_checks.md),
[`rho_trans()`](https://chengmatt.github.io/SPoRC/dev/reference/rho_trans.md),
[`set_data_indicator_unused()`](https://chengmatt.github.io/SPoRC/dev/reference/set_data_indicator_unused.md),
[`sporc_fit`](https://chengmatt.github.io/SPoRC/dev/reference/sporc_fit.md)

## Examples

``` r
if (FALSE) { # \dontrun{
obj <- fit_model(data, parameters, mapping, random = NULL, newton_loops = 3)
obj$rep$SSB
} # }
```
