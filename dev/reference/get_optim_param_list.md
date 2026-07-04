# Populate a parameter list with optimised values from sdreport

Replaces starting values in `parameters` with the corresponding
optimised estimates from `sd_rep`, respecting the factor-map sharing
structure in `mapping`. For each parameter: if a map exists, factor
level integers are used to index into the estimated value vector so that
shared elements receive the same optimised value and `NA`-mapped (fixed)
elements are left unchanged. Parameters absent from `mapping` are
treated as fully estimated and filled in sequentially. Random effects
are sourced from `sd_rep$par.random`; all other estimated parameters
from `sd_rep$par.fixed`.

## Usage

``` r
get_optim_param_list(parameters, mapping, sd_rep, random)
```

## Arguments

- parameters:

  Named list of parameter starting values passed to `RTMB::MakeADFun`.

- mapping:

  Named list of factor maps passed to `RTMB::MakeADFun`.

- sd_rep:

  `sdreport` object returned by `RTMB::sdreport`. Must contain
  `$par.fixed`, `$par.random`, and `$cov.fixed`.

- random:

  Character vector of parameter names declared as random effects in
  `RTMB::MakeADFun`. Used to route extraction to `sd_rep$par.random`
  rather than `sd_rep$par.fixed`.

## Value

The `parameters` list with all estimated elements replaced by their
optimised values. Fixed (`NA`-mapped) elements retain their original
starting values.
