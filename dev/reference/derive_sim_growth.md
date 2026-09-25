# Rebuild growth for every replicate from its deviation arrays

Runs the fit's own `Get_Growth` at each replicate's `ln_growth_devs` and
`ln_growth_semipar_devs` and writes what it returns over that
replicate's slices through `take_sim_growth_years`, so a drawn growth
series reaches the population the way it does in the fit. Under cohort
growth only the years before `growth_cohort_styr` are built here; the
rest come one year at a time from `advance_sim_growth_year` inside the
annual cycle, and the growth list each replicate advances from is kept
in `growth_state`.

## Usage

``` r
derive_sim_growth(sim_env)
```

## Arguments

- sim_env:

  Simulation environment holding `dsem_growth_args` from
  `Setup_Sim_DSEM` and the deviation arrays with the replicate dim last.

## Value

`invisible(NULL)`; `sim_env` is modified in place.
