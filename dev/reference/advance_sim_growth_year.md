# Advance one replicate's cohort growth by a year

The fit's `Get_Growth_Year` at this replicate's deviations and its own
start of year numbers at age, which blend the plus group, called from
`run_annual_cycle` before anything else in the year is formed, as the
fit's population loop does.

## Usage

``` r
advance_sim_growth_year(y, sim, sim_env)
```

## Arguments

- y:

  Year.

- sim:

  Replicate.

- sim_env:

  Simulation environment holding `growth_state` from
  `derive_sim_growth`.

## Value

`invisible(NULL)`; `sim_env` is modified in place.
