# Fill one stratum's growth output for a set of years from one start-of-year state

Evaluates growth at every timing a fleet, a composition or spawning
needs, builds the size-age transitions and the derived weight and
maturity, and writes them into the years named, which all start from the
same state.

## Usage

``` r
growth_fill_year(
  out,
  p,
  r,
  s,
  fill_yrs,
  growth_pars,
  L_beg,
  L1_birth,
  cv_ref_fn,
  a_prop,
  cohort,
  ages,
  growth_A1,
  growth_A2,
  growth_L0,
  growth_len_lower,
  growth_cv_type,
  growth_sd_type,
  growth_dist,
  derive_waa,
  w_mid,
  seasdur,
  spawn_seas,
  t_spawn,
  n_seas,
  t_fish,
  t_srv,
  len_devs = NULL,
  growth_L2_asymptote = 0
)
```
