# Extract and tabulate parameter estimates and metadata from a fitted SPoRC model

Joins the parameter list, factor map, and `sdreport` to produce two tidy
data frames: one for estimated parameters (with initial values,
posterior estimates, standard errors, and absolute gradients) and one
for fixed (non-estimated) parameters. Parameter indices are re-sequenced
to match the sequential numbering used internally by RTMB's `sdreport`,
and map entries of `NA` (i.e., parameters fixed via `mapping`) are
labeled `"NE"` (not estimated).

## Usage

``` r
get_par_est_info(parameters, mapping, sd_rep)
```

## Arguments

- parameters:

  Named list of parameter starting values passed to `RTMB::MakeADFun`.

- mapping:

  Named list of factor maps passed to `RTMB::MakeADFun`. Parameters
  absent from `mapping` are treated as freely estimated. Elements with
  `NA` factor levels are treated as fixed.

- sd_rep:

  `sdreport` object returned by `RTMB::sdreport`. Must contain
  `$par.fixed`, `$par.random`, `$cov.fixed`, `$diag.cov.random`, and
  `$gradient.fixed`.

## Value

A named list with two elements:

- `est_pars`:

  Data frame of estimated parameters, with columns `Par` (parameter
  name), `Est_Val` (estimated value on the native scale), `SE_Val`
  (standard error), `Abs_Grad_Val` (absolute gradient; `NA` for random
  effects), `Init_Val` (starting value), and dimension columns `Dim1`,
  `Dim2`, ... indicating the array indices of each element.

- `non_est_pars`:

  Data frame of non-estimated (fixed) parameters with the same dimension
  and `Init_Val` columns, plus `map = "NE"` indicating they were
  excluded from estimation.

## See also

Other Utility:
[`fit_model()`](https://chengmatt.github.io/SPoRC/dev/reference/fit_model.md),
[`get_logistN_Sigma()`](https://chengmatt.github.io/SPoRC/dev/reference/get_logistN_Sigma.md),
[`post_optim_sanity_checks()`](https://chengmatt.github.io/SPoRC/dev/reference/post_optim_sanity_checks.md),
[`rho_trans()`](https://chengmatt.github.io/SPoRC/dev/reference/rho_trans.md),
[`set_data_indicator_unused()`](https://chengmatt.github.io/SPoRC/dev/reference/set_data_indicator_unused.md)
