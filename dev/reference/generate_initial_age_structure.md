# Initialise age structure for a simulation replicate

Simulates or reads in initial age deviations and calls
[`Get_Init_NAA`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Init_NAA.md)
to compute both the fished and unfished equilibrium numbers-at-age for
year 1 and season 1. Results are written directly into the simulation
environment arrays `NAA` and `NAA0`. This function is called once per
simulation replicate at `y = 1` by
[`run_annual_cycle`](https://chengmatt.github.io/SPoRC/dev/reference/run_annual_cycle.md).

## Usage

``` r
generate_initial_age_structure(y, sim, sim_env)
```

## Arguments

- y:

  Integer. Year index (must be `1`).

- sim:

  Integer. Simulation replicate index.

- sim_env:

  Simulation environment created by
  [`Setup_sim_env`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_sim_env.md).
  Modified in place: `$ln_InitDevs`, `$NAA[,,1,1,,,sim]`, and
  `$NAA0[,,1,1,,,sim]` are updated.

## Value

`invisible(NULL)`. All modifications are made by reference within
`sim_env`.

## Details

Initial deviation sharing follows the same logic as the estimation
model: deviations are drawn once per population when `n_pop > 1`, or
once per region when `n_pop = 1` and `init_dd = 0` (local
density-dependence). If `ln_InitDevs_input` exists in the simulation
environment, those values are used directly rather than simulating new
draws. Populations with `R0 = 0` receive zero deviations. The
equilibrium solver uses `init_iter = n_ages × 5` iterations.
