# Helper function for extracting parameter information and names from TMB

Helper function for extracting parameter information and names from TMB

## Usage

``` r
get_par_est_info(parameters, mapping, sd_rep)
```

## Arguments

- parameters:

  Parameter list from setting up TMB object

- mapping:

  Mapping list from setting up TMB object

- sd_rep:

  SD Report from TMB obj

## Value

A list of dataframes for estimated and non-estimated parameter values.

## See also

Other Utility:
[`fit_model()`](https://chengmatt.github.io/SPoRC/reference/fit_model.md),
[`get_logistN_Sigma()`](https://chengmatt.github.io/SPoRC/reference/get_logistN_Sigma.md),
[`post_optim_sanity_checks()`](https://chengmatt.github.io/SPoRC/reference/post_optim_sanity_checks.md),
[`rho_trans()`](https://chengmatt.github.io/SPoRC/reference/rho_trans.md),
[`set_data_indicator_unused()`](https://chengmatt.github.io/SPoRC/reference/set_data_indicator_unused.md)
