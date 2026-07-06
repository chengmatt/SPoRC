# Compute spawning-time biomass quantities for one simulation year/season

Computes Total_Biom, SSB, Dynamic_SSB0, and eff_SSB for year `y` at
season `seas` (always called with `seas == spawn_seas`) from the current
`NAA`/`NAA0` state in `sim_env`. Factored out of
[`apply_pop_dy`](https://chengmatt.github.io/SPoRC/dev/reference/apply_pop_dy.md)
so it can be evaluated either before or after that season's
mortality/ageing step depending on `rec_lag`, without duplicating the
underlying math. Pure/read-only: returns a list rather than modifying
`sim_env`.

## Usage

``` r
compute_spawn_biomass_sim(y, seas, sim, sim_env)
```

## Arguments

- y:

  Year integer

- seas:

  Season integer

- sim:

  Simulation integer

- sim_env:

  Simulation environment
