# The growth arguments the operating model rebuilds growth with

`growth_args` from `Setup_Sim_Growth_RE` when the fit penalizes its own
growth deviations, else `dsem_growth_args` from `Setup_Sim_DSEM`.

## Usage

``` r
sim_growth_args(sim_env)
```

## Arguments

- sim_env:

  Simulation environment.

## Value

Named list of `Get_Growth` arguments, or `NULL`.
