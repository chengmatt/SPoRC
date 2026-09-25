# Compute Biomass

Spawning-time biomass quantities for year `y` from the annual cycle's
current state at season `seas`, always the spawning season.

## Usage

``` r
compute_biom_y(
  y,
  seas,
  NAA,
  NAA0,
  WAA,
  MatAA,
  ZAA,
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

  Year integer

- seas:

  Season integer

- NAA, NAA0:

  Abundance arrays `[pop, region, year, season, age, sex]`, fished and
  unfished.

- WAA, MatAA, ZAA, natmort:

  Weight, maturity, total mortality and the natural mortality rate on
  the same dims.

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
