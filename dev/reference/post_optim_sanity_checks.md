# Run post-optimisation convergence checks on a fitted SPoRC model

Evaluates four convergence criteria on the output of a fitted RTMB
model: (1) finiteness of the joint negative log-likelihood, (2) maximum
absolute gradient of fixed-effect parameters, (3) positive-definiteness
of the Hessian, and (4) finiteness and magnitude of parameter standard
errors and pairwise correlations. A diagnostic message is printed for
each failed check identifying the offending parameter or parameter pair.

## Usage

``` r
post_optim_sanity_checks(
  sd_rep,
  rep,
  gradient_tol = 0.001,
  se_tol = 100,
  corr_tol = 0.99
)
```

## Arguments

- sd_rep:

  `sdreport` object returned by `RTMB::sdreport`. Must contain
  `$gradient.fixed`, `$par.fixed`, `$pdHess`, and `$cov.fixed`.

- rep:

  Report list returned by `obj$rep` after fitting a SPoRC model. Must
  contain `$jnLL`.

- gradient_tol:

  Numeric. Maximum tolerated absolute gradient for any fixed-effect
  parameter. Default `1e-3`; values above this threshold suggest the
  optimiser did not reach a local minimum.

- se_tol:

  Numeric. Maximum tolerated parameter standard error. Default `100`;
  very large SEs indicate poorly identified parameters.

- corr_tol:

  Numeric. Maximum tolerated absolute pairwise parameter correlation.
  Default `0.99`; values approaching 1 indicate near-redundant
  parameters.

## Value

Logical. `TRUE` if all four checks pass; `FALSE` if any check fails. In
either case, informative messages are printed via `message`.

## See also

Other Utility:
[`fit_model()`](https://chengmatt.github.io/SPoRC/dev/reference/fit_model.md),
[`get_logistN_Sigma()`](https://chengmatt.github.io/SPoRC/dev/reference/get_logistN_Sigma.md),
[`get_par_est_info()`](https://chengmatt.github.io/SPoRC/dev/reference/get_par_est_info.md),
[`rho_trans()`](https://chengmatt.github.io/SPoRC/dev/reference/rho_trans.md),
[`set_data_indicator_unused()`](https://chengmatt.github.io/SPoRC/dev/reference/set_data_indicator_unused.md)
