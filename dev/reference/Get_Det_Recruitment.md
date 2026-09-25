# Deterministic Recruitment

Recruitment by population and region under mean, Beverton-Holt or Ricker
recruitment, spread over regions and seasons by the recruitment
proportions. Unfished spawning biomass per recruit is computed
internally by projecting one recruit through every age and season.
Equations are in the model equations vignette.

## Usage

``` r
Get_Det_Recruitment(
  recruitment_model,
  rec_dd,
  y,
  rec_lag,
  R0,
  rec_region_prop,
  rec_seas_prop,
  h,
  n_pop,
  n_regions,
  n_ages,
  n_fish_fleets,
  WAA,
  MatAA,
  natmort,
  SSB_vals,
  Movement,
  sgl_seas_spawning_movement,
  stray_rate,
  do_recruits_move,
  t_spawn,
  init_F,
  dmr,
  fish_sel,
  ret_sel,
  n_seas,
  spawn_seas,
  natal_region,
  seasdur,
  sexratio_f,
  Mrate = NULL,
  move_timing = 0,
  expm_nsub = 0
)
```

## Arguments

- recruitment_model:

  Integer. 0 = mean recruitment, 1 = Beverton-Holt, 2 = Ricker.

- rec_dd:

  Integer. 0 = density dependence within each population or region, 1 =
  shared across regions, valid only when `n_pop = 1`.

- y:

  Current model year index.

- rec_lag:

  Lag in seasons between spawning and recruitment. 1 uses `SSB_vals`
  from that many seasons prior, 0 uses the same year's SSB, which the
  caller must supply computed from survivors only.

- R0:

  Numeric vector (`n_pop`) of unfished recruitment by population.

- rec_region_prop:

  Matrix (`n_pop × n_regions`) of the proportion of recruitment
  allocated to each region.

- rec_seas_prop:

  Matrix (`n_pop × n_seas`) of seasonal recruitment proportions. Must be
  zero before `spawn_seas` when `rec_lag = 0`.

- h:

  Matrix (`n_pop × n_regions`) of steepness values.

- n_pop:

  Number of populations.

- n_regions:

  Number of spatial regions.

- n_ages:

  Number of age classes (including the plus group).

- n_fish_fleets:

  Integer. Number of fishery fleets.

- WAA:

  Array (`n_pop × n_regions × n_seas × n_ages`) of weight-at-age.

- MatAA:

  Array (`n_pop × n_regions × n_seas × n_ages`) of maturity-at-age.

- natmort:

  Array (`n_pop × n_regions × n_seas × n_ages`) of natural mortality, a
  rate per year in each season.

- SSB_vals:

  Array (`n_pop × n_regions × n_years`) of spawning biomass.

- Movement:

  Array (`n_pop × origin × destination × n_seas × n_ages`) of seasonal
  movement probabilities.

- sgl_seas_spawning_movement:

  Array (`n_pop × origin × destination × n_ages`) of spawning movement
  when a single season is used and `n_pop > 1`.

- stray_rate:

  Numeric vector of stray rates by population, scaling each other
  population's contribution to recruitment in a natal region.

- do_recruits_move:

  Indicator for whether recruits move in their first year.

- t_spawn:

  Fraction of the spawning season that occurs before spawning.

- init_F:

  Array (`n_regions × n_seas × n_fish_fleets`) of initial fishing
  mortality.

- dmr:

  Array (`n_regions × n_seas × n_fish_fleets`) of initial (first year)
  discard mortality.

- fish_sel:

  Array (`n_pop × n_regions × n_seas × n_ages × n_fish_fleets`) of total
  fishery selectivity.

- ret_sel:

  Array (`n_pop × n_regions × n_seas × n_ages × n_fish_fleets`) of
  retained fishery selectivity.

- n_seas:

  Number of seasons per year.

- spawn_seas:

  Season index in which spawning occurs.

- natal_region:

  Integer vector (`n_pop`) mapping each population to its natal region.

- seasdur:

  Numeric vector (`n_seas`) of seasonal durations as fractions of a
  year.

- sexratio_f:

  Matrix (`n_pop × n_regions`) of female recruitment proportions.
