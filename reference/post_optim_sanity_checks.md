# Post Optimization Model Convergence Checks

Post Optimization Model Convergence Checks

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

  sd report list from a \`SPoRC\` model

- rep:

  report list from a \`SPoRC\` model

- gradient_tol:

  Value for maximum gradient tolerance to use

- se_tol:

  Value for maximum standard error tolerance to use

- corr_tol:

  Value for maximum correlation tolerance to use

## See also

Other Utility:
[`fit_model()`](https://chengmatt.github.io/SPoRC/reference/fit_model.md),
[`get_logistN_Sigma()`](https://chengmatt.github.io/SPoRC/reference/get_logistN_Sigma.md),
[`get_par_est_info()`](https://chengmatt.github.io/SPoRC/reference/get_par_est_info.md),
[`rho_trans()`](https://chengmatt.github.io/SPoRC/reference/rho_trans.md),
[`set_data_indicator_unused()`](https://chengmatt.github.io/SPoRC/reference/set_data_indicator_unused.md)
