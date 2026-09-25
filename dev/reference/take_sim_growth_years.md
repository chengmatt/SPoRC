# Write years of one replicate's growth into the operating model

Weight at age, the size-age keys, and then what the fit forms from them
in `compute_mortality_year`: a fleet asking for it gets the weight of
what it selects, and with selectivity at length each fleet's selectivity
at age is its curve read through the key. Weight at age is left alone
when the fit takes it as data (`derive_waa = 0`).

## Usage

``` r
take_sim_growth_years(sim_env, sim, growth, yrs)
```

## Arguments

- sim_env:

  Simulation environment holding `dsem_length_sel` from
  `Setup_Sim_DSEM`.

- sim:

  Replicate.

- growth:

  Growth list from `Get_Growth` or `Get_Growth_Year`.

- yrs:

  Years to write.

## Value

`invisible(NULL)`; `sim_env` is modified in place.
