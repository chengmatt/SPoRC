# Print, summarize, compare and plot a fitted SPoRC model

[`fit_model`](https://chengmatt.github.io/SPoRC/dev/reference/fit_model.md)
returns a list of class `sporc_fit`. Typing its name prints a short
description of the fit in place of the whole list,
[`summary()`](https://rdrr.io/r/base/summary.html) tabulates the
estimates by parameter, and
[`logLik()`](https://rdrr.io/r/stats/logLik.html) returns the maximized
log likelihood, so `AIC(fit_1, fit_2)` compares fits.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws one of
the package's plots for this fit alone; to compare fits, call the
plotting function itself with a list of them. The fit is still a list:
`fit$rep`, `fit$sd_rep` and `RTMB::sdreport(fit)` work as before, and
`unclass(fit)` prints every element.

## Usage

``` r
# S3 method for class 'sporc_fit'
print(x, digits = 6, ...)

# S3 method for class 'sporc_fit'
summary(object, ...)

# S3 method for class 'sporc_fit'
logLik(object, ...)

# S3 method for class 'sporc_fit'
plot(x, type = "ts", model_name = "Model 1", ...)
```

## Arguments

- x:

  A fit returned by
  [`fit_model`](https://chengmatt.github.io/SPoRC/dev/reference/fit_model.md).

- digits:

  Number of significant digits printed for the objective and the largest
  absolute gradient. Default `6`.

- ...:

  For `plot`, further arguments of the plotting function, such as
  `do_ci`, `Selex_Type`, `year_indx` or `data_source`. Not used by the
  other methods.

- object:

  A fit returned by
  [`fit_model`](https://chengmatt.github.io/SPoRC/dev/reference/fit_model.md).

- type:

  Which plot to draw, each from the function named: `"ts"` (the default,
  [`get_ts_plot`](https://chengmatt.github.io/SPoRC/dev/reference/get_ts_plot.md),
  which needs `do_sdrep = TRUE`), `"selex"`
  ([`get_selex_plot`](https://chengmatt.github.io/SPoRC/dev/reference/get_selex_plot.md)),
  `"biological"`
  ([`get_biological_plot`](https://chengmatt.github.io/SPoRC/dev/reference/get_biological_plot.md)),
  `"data_fitted"`
  ([`get_data_fitted_plot`](https://chengmatt.github.io/SPoRC/dev/reference/get_data_fitted_plot.md)),
  `"nLL"`
  ([`get_nLL_plot`](https://chengmatt.github.io/SPoRC/dev/reference/get_nLL_plot.md)),
  `"idx_fits"`
  ([`get_idx_fits_plot`](https://chengmatt.github.io/SPoRC/dev/reference/get_idx_fits_plot.md)),
  `"catch_fits"`
  ([`get_catch_fits_plot`](https://chengmatt.github.io/SPoRC/dev/reference/get_catch_fits_plot.md))
  or `"at_age_fits"`
  ([`get_at_age_fits_plot`](https://chengmatt.github.io/SPoRC/dev/reference/get_at_age_fits_plot.md)).

- model_name:

  Name the fit is given in the plot legend. Default `"Model 1"`.

## Value

`plot` returns what the plotting function returns, a `ggplot` or a list
of them. `print` returns the fit invisibly. `summary` returns a data
frame with one row per parameter name: `Par`, `Type` (`"fixed"` or
`"random"`), `N` (number of estimated values), `Est_Min` and `Est_Max`
(range of the estimates, on the scale the parameter is estimated on),
`SE_Max` (largest standard error) and `Abs_Grad_Max` (largest absolute
gradient, `NA` for random effects).
[`get_par_est_info`](https://chengmatt.github.io/SPoRC/dev/reference/get_par_est_info.md)
gives one row per estimated value. `logLik` returns minus the objective
the fit minimized, with `df` the number of estimated fixed effects. That
objective includes every penalty and prior in the model, and with
`random` set it is the marginal. No number of observations is set, so
[`BIC()`](https://rdrr.io/r/stats/AIC.html) is not available. A fit run
with `do_optim = FALSE` has no estimates, and `summary` and `logLik`
stop on it.

## Details

`print` says whether `nlminb` converged and gives the largest absolute
gradient at the estimates, both stored with every optimized fit. The
Hessian check and the standard errors are read from `$sd_rep`, so they
are shown only for a fit run with `do_sdrep = TRUE`.

## See also

Other Utility:
[`fit_model()`](https://chengmatt.github.io/SPoRC/dev/reference/fit_model.md),
[`get_logistN_Sigma()`](https://chengmatt.github.io/SPoRC/dev/reference/get_logistN_Sigma.md),
[`get_par_est_info()`](https://chengmatt.github.io/SPoRC/dev/reference/get_par_est_info.md),
[`label_dims()`](https://chengmatt.github.io/SPoRC/dev/reference/label_dims.md),
[`post_optim_sanity_checks()`](https://chengmatt.github.io/SPoRC/dev/reference/post_optim_sanity_checks.md),
[`rho_trans()`](https://chengmatt.github.io/SPoRC/dev/reference/rho_trans.md),
[`set_data_indicator_unused()`](https://chengmatt.github.io/SPoRC/dev/reference/set_data_indicator_unused.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- fit_model(data, parameters, mapping, do_sdrep = TRUE)
fit
summary(fit)
AIC(fit)
plot(fit, type = "selex")
} # }
```
