# Compute Biomass for Population Projections

Compute Biomass for Population Projections

## Usage

``` r
derive_proj_biom(
  y,
  seas,
  proj_NAA,
  proj_NAA0,
  WAA,
  MatAA,
  proj_ZAA,
  natmort,
  t_spawn,
  seasdur,
  n_seas,
  n_pop,
  n_regions,
  n_ages,
  n_sexes,
  sgl_seas_spawning_movement,
  natal_region,
  stray_rate,
  Movement = NULL,
  Mrate = NULL,
  move_timing = 0,
  do_recruits_move = 1
)
```

## Arguments

- y:

  Projection year integer

- seas:

  Season integer (always spawn_seas)
