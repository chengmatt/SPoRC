# Run RTMB model

Run RTMB model

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
  ...
)
```

## Arguments

- data:

  Data list

- parameters:

  Parameter list

- mapping:

  Mapping list

- random:

  Character of random effects to integrate out

- newton_loops:

  Number of newton loops to run to get gradients down

- silent:

  Boolean on whether or not model run is silent

- do_optim:

  Boolean on whether or not model is optimized

- nlminb_control:

  List argument controls by nlminb

- ...:

  Additional arguments taken by MakeADFun

## Value

Returns a list object that is optimized, with results outputted from the
RTMB model

## See also

Other Utility:
[`get_logistN_Sigma()`](https://chengmatt.github.io/SPoRC/reference/get_logistN_Sigma.md),
[`get_par_est_info()`](https://chengmatt.github.io/SPoRC/reference/get_par_est_info.md),
[`post_optim_sanity_checks()`](https://chengmatt.github.io/SPoRC/reference/post_optim_sanity_checks.md),
[`rho_trans()`](https://chengmatt.github.io/SPoRC/reference/rho_trans.md),
[`set_data_indicator_unused()`](https://chengmatt.github.io/SPoRC/reference/set_data_indicator_unused.md)

## Examples

``` r
if (FALSE) { # \dontrun{
model <- fit_model(data,
                 parameters,
                 mapping,
                 random = NULL,
                 newton_loops = 3)
} # }
```
