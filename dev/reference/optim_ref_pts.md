# Optimize Reference Point Models

Constructs and optimizes an RTMB automatic differentiation objective
function for estimating fisheries reference points (e.g., SPR-based or
Fmsy-based biological reference points). After optimization, retrieves
the model report and standard deviation report from the best parameter
estimates.

## Usage

``` r
optim_ref_pts(model_name, data_list, pars_list)
```

## Arguments

- model_name:

  Function. An RTMB-compatible model function (e.g., `global_SPR`,
  `global_Fmsy`, `local_Fmsy_sglpop`) that defines the objective
  function for the reference point calculation.

- data_list:

  List. A named list of data inputs passed to `model_name` via
  `MakeADFun`. Should contain all quantities treated as fixed data
  within the reference point model (e.g., biological parameters,
  selectivity, movement matrices).

- pars_list:

  List. A named list of initial parameter values passed to `MakeADFun`.
  These are the parameters over which the objective function is
  optimized (e.g., `ln_Fmsy`, `ln_F_spr`).

## Value

An RTMB AD function object (list) with the following additional elements
appended after optimization:

- `$optim`:

  Output from [`nlminb`](https://rdrr.io/r/stats/nlminb.html), including
  convergence code, final objective value, and optimized parameter
  estimates.

- `$rep`:

  Named list of reported quantities from the model (e.g., equilibrium
  SSB, yield, reference point values), evaluated at the best parameter
  estimates via `$env$last.par.best`.

- `$sd_rep`:

  Output from `sdreport`, containing standard errors and summary
  statistics for all estimated and derived quantities.
