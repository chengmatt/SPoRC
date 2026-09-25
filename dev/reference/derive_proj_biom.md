# Compute Biomass for Population Projections

Spawning-time biomass for projection year `y` at the spawning season,
from the current projected state. Plain R, so it can run either side of
mortality depending on `rec_lag`.

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
  do_recruits_move = 1,
  expm_nsub = 0
)
```

## Arguments

- y:

  Projection year integer

- seas:

  Season integer (always spawn_seas)

- proj_NAA, proj_NAA0:

  Projected abundance arrays, fished and unfished.

- proj_ZAA:

  Projected total mortality.

- t_spawn:

  Fraction of the season elapsed at spawning.

- seasdur:

  Numeric vector of season durations.

- n_seas, n_pop, n_regions, n_ages, n_sexes:

  Model dimensions.

- sgl_seas_spawning_movement:

  Natal homing movement `[pop, region, region, year, age, sex]`.

- natal_region:

  Integer vector of each population's natal region.

- stray_rate:

  Matrix `[pop, year]` of stray rates.

- Movement, Mrate:

  Movement fractions and generator
  `[pop, region, region, year, season, age, sex]`, or `NULL`.

- move_timing, do_recruits_move, expm_nsub:

  Movement timing, whether age one moves, and substeps for the implicit
  matrix exponential.

## Value

List of `Total_Biom_y`, `SSB_y`, `Dynamic_SSB0_y` and `eff_SSB_y`.
