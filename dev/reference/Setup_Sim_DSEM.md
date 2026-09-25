# Set up dsem series in a simulation list

Two ways in. From a fit: pass `data` and `pars` and the operating model
draws from the fitted arrows at their estimates, either every year or
with the fitted years kept and only later years drawn. From scratch:
pass `dsem_arrows` and `dsem_values` and it draws from arrows you wrote
at values you chose, which is how a self test states its truth.

## Usage

``` r
Setup_Sim_DSEM(
  sim_list,
  data = NULL,
  pars = NULL,
  rep = NULL,
  dsem_cov_use = NULL,
  condition_on_fit = FALSE,
  dsem_arrows = NULL,
  dsem_values = NULL,
  dsem_processes = "rec",
  dsem_cov_mu = NULL,
  dsem_cov_obs_sd = NULL,
  dsem_cov_family = NULL,
  dsem_cov_tweedie_p = NULL,
  dsem_cov_link = NULL,
  mod_var_logscale = FALSE,
  dsem_variance = "conditional"
)
```

## Arguments

- sim_list:

  Simulation list, for example from `condition_closed_loop_simulations`
  or the `Setup_Sim_*` calls.

- data:

  Data list of the fitted model, set up by `Setup_Mod_DSEM`. `NULL` from
  scratch.

- pars:

  Parameter list at the fitted values, as `get_optim_param_list` or
  `obj$env$parList()` returns it. `NULL` from scratch.

- rep:

  Report of the fit. Needed for a numbers at age series, whose
  `NAA_pred` turns the fitted log state into innovations, and for a
  growth series linked under selectivity at length, which the rebuilt
  keys are read through.

- dsem_cov_use:

  Optional matrix `[year, covariate]` of 0/1 flags for the years each
  covariate is observed in the simulation. Defaults to the fitted
  pattern for fitted years and every year after, or every year from
  scratch.

- condition_on_fit:

  `TRUE` keeps every fitted year at its estimate and draws later years
  given them. `FALSE` (default) draws every year.

- dsem_arrows:

  From scratch: the arrow lines, one per element or one per line of a
  string, naming covariates and deviation series.

- dsem_values:

  From scratch: named vector giving every arrow parameter its value, sds
  on the natural scale, for example
  `c(b_env = 0.5, rho_env = 0.6, sd_env = 1, sd_rec = 0.8)`.

- dsem_processes:

  From scratch: processes the arrows may name, `"rec"` (default) or
  `"NAA"`, or both.

- dsem_cov_mu:

  From scratch: named vector of covariate means, one per covariate the
  arrows name. Default zero for each.

- dsem_cov_obs_sd:

  From scratch: named vector, one per covariate, of the spread of its
  observations on its family's own scale: the sd for `"normal"`, the CV
  for `"gamma"`, the sd of the log for `"lognormal"`, the dispersion for
  `"tweedie"`; `NA` for a covariate observed without error. Default
  `NA`.

- dsem_cov_family:

  From scratch: named character vector, one per covariate, of the family
  the replicates' observations are drawn from, as in
  [`Setup_Mod_DSEM`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_DSEM.md).
  Default `"normal"` where `dsem_cov_obs_sd` is given and `"fixed"`
  otherwise. A link-scale family reads `dsem_cov_mu` and the arrows on
  that scale.

- dsem_cov_tweedie_p:

  From scratch: named vector of tweedie powers, strictly between 1
  and 2. Default 1.5.

- dsem_cov_link:

  From scratch: named character vector of links, as in
  [`Setup_Mod_DSEM`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_DSEM.md).
  Default each family's own. A `"gaussian_fixed_sd"` covariate from
  scratch takes its one sd from `dsem_cov_obs_sd` for every year.

- mod_var_logscale:

  From scratch: passed to `read_dsem_arrows`.

- dsem_variance:

  From scratch: passed to `read_dsem_arrows` as `variance`,
  `"conditional"` (default) or `"diagonal"`. From a fit the data list's
  own setting is used.

## Value

`sim_list` with the dsem arrows, values and observation pattern added,
plus the deviation arrays and the growth or movement arguments the
operating model needs to consume a drawn series.

## What a drawn series feeds

Each linked series is written into the array the operating model reads,
at the cells the fit linked, so those cells are never copied from the
report:

- recruitment into `ln_RecDevs`, which `generate_recruitment` reads.

- catchability into `ln_fish_q_devs` or `ln_srv_q_devs`, and
  [`draw_sim_q_devs`](https://chengmatt.github.io/SPoRC/dev/reference/draw_sim_q_devs.md)
  scales that fleet's catchability by them whatever the operating
  model's own `fish_q_model` or `srv_q_model` says. The catchability the
  deviations scale is the block mean, so the series drives the
  conditioning years and every year past them alike.

- the numbers at age into `naa_eta_all`, as the log state less its
  deterministic prediction, replacing those cells of the replicate's own
  innovations.

- growth and movement into arrays built here, since the operating model
  has none of its own: fitted years at their estimates and later years
  at zero. The fit's own `Get_Growth` or `Get_Movement` is then rerun
  per replicate to rebuild weight at age, the size-age keys and the
  movement matrix, reading every argument from `data` and `pars` by
  name. Under selectivity at length a fleet's selectivity at age is its
  curve read through the rebuilt key, as in the fit.

Growth propagated cohort by cohort (`growth_tv_type = 1`) runs the way
the fit runs it: the years before `growth_cohort_styr` are built up
front, and from there `run_annual_cycle` advances each replicate one
year at a time from its own start of year numbers. From scratch only
recruitment and the numbers at age can be linked, since those are the
arrays the operating model has itself.

## See also

[`Setup_Mod_DSEM`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_DSEM.md)
