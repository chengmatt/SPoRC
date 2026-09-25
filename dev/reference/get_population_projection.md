# Population projection (numbers-at-age dynamics)

Advances numbers at age through every modeled year and season: inserts
recruitment at the timing `rec_lag` sets, applies movement, computes the
biomass quantities through `compute_biom_y`, and applies mortality and
ageing. `ZAA` is an input here, not derived from `NAA`, so it must
already be computed. Called once from the population projection section
of `SPoRC_rtmb.R`.

## Usage

``` r
get_population_projection(
  n_pop,
  n_regions,
  n_seas,
  n_ages,
  n_sexes,
  n_yrs,
  n_fish_fleets,
  n_est_rec_devs,
  rec_lag,
  rec_model,
  rec_dd,
  R0,
  rec_region_prop,
  rec_seas_prop,
  h_trans,
  R0_yr = NULL,
  natal_region,
  t_spawn,
  spawn_seas,
  seasdur,
  init_F,
  ln_RecDevs,
  sexratio,
  WAA,
  MatAA,
  natmort,
  Movement,
  stray_rate,
  sgl_seas_spawning_movement,
  do_recruits_move,
  fish_sel,
  ret_sel,
  dmr,
  ZAA,
  NAA,
  NAA0,
  NAA_bef,
  NAA_aft,
  Rec,
  SSB,
  Total_Biom,
  Dynamic_SSB0,
  eff_SSB,
  Mrate = NULL,
  move_timing = 0,
  SR_ref_yr = 1,
  sr_penalty = 0,
  sr_R0 = NULL,
  growth_mortality_year_fn = NULL,
  growth_mortality_state = NULL,
  expm_nsub = 0,
  n_est_naa_re = 0,
  ln_NAA = NULL,
  naa_re_ages = NULL,
  naa_re_yrs = NULL,
  naa_re_seas = NULL
)
```

## Arguments

- n_pop, n_regions, n_seas, n_ages, n_sexes, n_yrs, n_fish_fleets:

  Dimension sizes.

- n_est_rec_devs:

  Number of estimated recruitment deviations.

- rec_lag:

  Integer. `0` inserts recruitment inside the spawning-season biomass
  computation; non-zero inserts it once a year ahead of the seasonal
  loop.

- rec_model, rec_dd, R0, rec_region_prop, rec_seas_prop, h_trans,
  natal_region, t_spawn, spawn_seas, seasdur, init_F:

  Recruitment and timing arguments passed through to
  `Get_Det_Recruitment`.

- R0_yr:

  Matrix `[n_pop x n_yrs]` of R0 by year when R0 has time blocks, or
  `NULL` for the single `R0`. Only the recruitment computed each year
  reads it; everything needing one value still uses `R0`.

- ln_RecDevs:

  Array `[pop, region, year]` of log recruitment deviations, applied
  multiplicatively to deterministic recruitment for
  `y <= n_est_rec_devs`.

- sexratio:

  Array `[pop, region, year, sex]` of the recruitment sex ratio.

- WAA, MatAA, natmort:

  Arrays `[pop, region, year, season, age, sex]` of weight at age,
  maturity at age and natural mortality.

- Movement:

  Array `[pop, region_from, region_to, year, season, age, sex]` of
  movement rates.

- stray_rate:

  Array `[pop, year]` of stray rate.

- sgl_seas_spawning_movement:

  Array `[pop, region_from, region_to, year, age, sex]` of single-season
  spawning movement.

- do_recruits_move:

  Integer (0/1) for whether age-1 recruits move.

- fish_sel, ret_sel:

  Arrays `[pop, region, year, season, age, sex, fish_fleet]` of total
  and retained fishery selectivity.

- dmr:

  Array `[region, year, season, fish_fleet]` of discard mortality rate.

- ZAA:

  Array `[pop, region, year, season, age, sex]` of total mortality at
  age, precomputed.

- NAA, NAA0, NAA_bef, NAA_aft:

  Arrays `[pop, region, year+1, season, age, sex]`, the output
  containers for the fished and unfished numbers at age and for the
  numbers immediately before and after movement.

- Rec:

  Array `[pop, region, year]`, the output container for total
  recruitment before seasonal apportionment.

- SSB, Total_Biom, Dynamic_SSB0:

  Arrays `[pop, region, year]`, output containers.

- eff_SSB:

  Array `[pop, year]`, the output container for effective
  (natal-homing-adjusted) SSB.

- SR_ref_yr:

  Integer year index supplying the biological inputs (weight at age,
  maturity, natural mortality and movement) to unfished spawning biomass
  per recruit, and so to `S0` and the curve's scale. Default `1`. Set it
  to `n_yrs` to condition the curve on terminal weight at age, as
  several ADMB assessments do. It is an index, not a calendar year, so a
  caller that truncates the year dim must clamp it.

- growth_mortality_year_fn:

  Optional function of `(y, NAA_y, growth_mortality_state)` called at
  the top of every year with the numbers at age at the start of that
  year, array `[pop, region, age, sex]`, and the state kept from the
  previous year. It returns a list with `state`, advanced to the next
  call and returned to the caller, and `ZAA_y`, `WAA_y` and `MatAA_y`,
  which replace that year's slices. Passing the state in and out keeps
  the per-year step a function of its arguments.

- growth_mortality_state:

  Initial state for `growth_mortality_year_fn`, passed through the year
  loop and returned. `NULL` (default) uses the arrays as given. This is
  how cohort growth, whose plus group blends by numbers, is evaluated
  inside the year loop.

- n_est_naa_re:

  Number of estimated state-space numbers at age. Zero leaves the
  numbers deterministic. Never inferred from `dim(ln_NAA)`, which is
  non-zero once the setup function has run at all.

- ln_NAA:

  Array `[pop, region, year, season, age, sex]` of log numbers at the
  start of a season, overwriting the deterministic prediction wherever
  the state is active. Season one is the year boundary, after ageing and
  the plus group; later seasons are states on the within-year survival
  step.

- naa_re_ages, naa_re_yrs, naa_re_seas:

  Integer index vectors the state is active over.

## Value

List with `NAA`, `NAA0`, `NAA_bef`, `NAA_aft`, `Rec`, `SSB`,
`Total_Biom`, `Dynamic_SSB0`, `eff_SSB`, `Aggregated_SSB` and
`Dynamic_Aggregated_SSB0` (arrays `[year]`, summed across population and
region), and `NAA_int` `[pop, region, year, season, age, sex]`, the
season-integrated abundance the spatial Baranov equation needs,
populated only under `move_timing = 2` and all zeros otherwise.

## Details

Every array argument matching an output name is passed in already
dimensioned, usually all zero aside from any initial-year values
inserted upstream, and comes back filled over `1:n_yrs`.
