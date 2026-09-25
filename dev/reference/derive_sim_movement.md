# Rebuild movement for every replicate from its deviation array

Runs the fit's own `Get_Movement` at each replicate's `move_devs` and
writes the movement matrix (and the rate matrix under continuous
movement) over that replicate's slice.

## Usage

``` r
derive_sim_movement(sim_env)
```

## Arguments

- sim_env:

  Simulation environment holding `dsem_move_args` from `Setup_Sim_DSEM`
  and `move_devs` with the replicate dim last.

## Value

`invisible(NULL)`; `sim_env` is modified in place.
