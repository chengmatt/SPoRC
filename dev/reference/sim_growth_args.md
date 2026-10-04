# The growth arguments the operating model rebuilds growth with

`growth_args` from `Setup_Sim_Growth_RE` when the fit penalizes its own
growth deviations, else `dsem_growth_args` from `Setup_Sim_DSEM`. A
replicate on its own parameter draw takes its own.

## Usage

``` r
sim_growth_args(sim_env, sim = NULL)
```

## Arguments

- sim_env:

  Simulation environment.

- sim:

  Replicate, or `NULL` for the arguments every replicate shares.

## Value

Named list of `Get_Growth` arguments, or `NULL`.
