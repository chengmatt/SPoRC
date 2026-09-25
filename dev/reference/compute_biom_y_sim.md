# Compute spawning-time biomass quantities for one simulation year/season

The operating model's state at year `y` and season `seas`, always the
spawning season, sliced at replicate `sim` and given to
[`biom_at_spawn`](https://chengmatt.github.io/SPoRC/dev/reference/biom_at_spawn.md),
which the estimation model and the forward projection also run.

## Usage

``` r
compute_biom_y_sim(y, seas, sim, sim_env)
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
