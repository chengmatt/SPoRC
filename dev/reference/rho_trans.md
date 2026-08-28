# Transform a real-valued parameter to the interval (-1, 1)

Applies the scaled logistic transformation \\2 / (1 + e^{-2x}) - 1\\ to
map an unconstrained real value to \\(-1, 1)\\, suitable for
parameterizing correlation coefficients. Used in SPoRC to back-transform
raw correlation parameters before constructing AR(1) and constant
covariance matrices for logistic-normal composition likelihoods.

## Usage

``` r
rho_trans(x)
```

## Arguments

- x:

  Numeric. Unconstrained real-valued parameter.

## Value

Numeric. Transformed value in \\(-1, 1)\\.

## See also

Other Utility:
[`fit_model()`](https://chengmatt.github.io/SPoRC/dev/reference/fit_model.md),
[`get_logistN_Sigma()`](https://chengmatt.github.io/SPoRC/dev/reference/get_logistN_Sigma.md),
[`get_par_est_info()`](https://chengmatt.github.io/SPoRC/dev/reference/get_par_est_info.md),
[`post_optim_sanity_checks()`](https://chengmatt.github.io/SPoRC/dev/reference/post_optim_sanity_checks.md),
[`set_data_indicator_unused()`](https://chengmatt.github.io/SPoRC/dev/reference/set_data_indicator_unused.md)
