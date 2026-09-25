# A data list and parameter list for a dsem written from scratch

Stands in for a fit when an operating model states its own dsem: the
arrows are read against the operating model's arrays (every cell
estimated, no map), and the values given become the parameters. Only
recruitment and the numbers at age can be named, since those are the
arrays the operating model has.

## Usage

``` r
scratch_dsem_fit(
  sim_list,
  dsem_arrows,
  dsem_values,
  dsem_processes,
  dsem_cov_mu,
  dsem_cov_obs_sd,
  mod_var_logscale,
  dsem_variance = "conditional",
  dsem_cov_family = NULL,
  dsem_cov_tweedie_p = NULL,
  dsem_cov_link = NULL
)
```

## Arguments

- sim_list:

  Simulation list holding the dims and the arrays.

- dsem_arrows, dsem_values, dsem_processes, dsem_cov_mu,
  dsem_cov_obs_sd, mod_var_logscale, dsem_variance, dsem_cov_family,
  dsem_cov_tweedie_p, dsem_cov_link:

  As in
  [`Setup_Sim_DSEM`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_DSEM.md).

## Value

List with `data` and `pars` holding the dsem fields `Setup_Sim_DSEM`
reads from a fit.
