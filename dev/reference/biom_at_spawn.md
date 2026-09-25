# Spawning-time biomass for one year and season

The calculation doing the annual cycle, the forward projection and the
operating model.

## Usage

``` r
biom_at_spawn(
  NAA_s,
  NAA0_s,
  WAA_s,
  MatAA_s,
  ZAA_s,
  natmort_s,
  spawn_move_s,
  stray_rate_y,
  t_spawn,
  seasdur_seas,
  n_seas,
  n_pop,
  n_regions,
  n_ages,
  n_sexes,
  natal_region,
  Movement_s = NULL,
  Mrate_s = NULL,
  move_timing = 0,
  do_recruits_move = 1,
  expm_nsub = 0
)
```

## Arguments

- NAA_s:

  Array `[pop, region, 1, 1, age, sex]` of abundance at the start of the
  season.

- NAA0_s:

  Array on the same dims of unfished abundance.

- WAA_s:

  Array on the same dims of weight at age.

- MatAA_s:

  Array on the same dims of maturity at age.

- ZAA_s:

  Array on the same dims of total mortality over the season.

- natmort_s:

  Array on the same dims of the natural mortality rate.

- spawn_move_s:

  Array `[pop, region, region, 1, age, sex]` of natal homing movement,
  read only by a single-season multi-population model, `NULL` otherwise.

- stray_rate_y:

  Numeric `[pop]` of this year's stray rates, read only by a
  multi-population model, `NULL` otherwise.

- t_spawn:

  Fraction of the season elapsed at spawning.

- seasdur_seas:

  Duration of this season as a fraction of the year.

- n_seas, n_pop, n_regions, n_ages, n_sexes:

  Model dimensions.

- natal_region:

  Integer vector of each population's natal region.

- Movement_s:

  Array `[pop, region, region, 1, 1, age, sex]` of movement fractions,
  or `NULL` when movement does not run at spawning.

- Mrate_s:

  Array on the same dims of the movement generator, or `NULL`.

- move_timing:

  `0` movement before mortality, `1` and `2` movement within the season,
  which folds the mortality discount into the propagation.

- do_recruits_move:

  Whether age one moves.

- expm_nsub:

  Substeps for the implicit matrix exponential.

## Value

List of `Total_Biom_y`, `SSB_y`, `Dynamic_SSB0_y` and `eff_SSB_y`.
