# Compute spawning-time biomass quantities for one simulation year/season

The operating model's state at year `y` and season `seas`, always the
spawning season, sliced at replicate `sim` and given to
[`biom_at_spawn`](https://chengmatt.github.io/SPoRC/dev/reference/biom_at_spawn.md),
which the estimation model and the forward projection also run.

## Usage

``` r
compute_biom_y_sim(y, seas, sim, sim_env, NAA_s = NULL, ZAA_s = NULL)
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

- NAA_s, ZAA_s:

  Arrays `[n_pop, n_regions, 1, 1, n_ages, n_sexes]` of numbers and
  total mortality to use instead of the stored ones, for a trial F.
  `NULL` reads `sim_env`.
