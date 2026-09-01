# Population projection (numbers-at-age dynamics)

Advances numbers-at-age forward through all modeled years and seasons:
inserts recruitment (timing controlled by `rec_lag`), applies movement,
computes SSB/biomass quantities via `compute_biom_y`, and applies
mortality/ageing. Called once from the "Population Projection" section
of `SPoRC_rtmb.R`. `ZAA` (total mortality at age) must already be
computed before calling this, since it is treated as an input here
rather than derived from `NAA`.

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
  naa_re_yrs = NULL
)
```

## Arguments

- n_pop, n_regions, n_seas, n_ages, n_sexes, n_yrs, n_fish_fleets:

  Dimension sizes.

- n_est_rec_devs:

  Number of estimated recruitment deviations.

- rec_lag:

  Integer. Recruitment timing: `0` inserts recruitment within the
  spawning-season biomass computation; non-zero inserts recruitment once
  per year ahead of the seasonal loop.

- rec_model, rec_dd, R0, rec_region_prop, rec_seas_prop, h_trans,
  natal_region, t_spawn, spawn_seas, seasdur, init_F:

  Recruitment and timing arguments passed through to
  `Get_Det_Recruitment`.

- ln_RecDevs:

  Array `[pop, region, year]` of log recruitment deviations; applied
  multiplicatively to deterministic recruitment for
  `y <= n_est_rec_devs`.

- sexratio:

  Array `[pop, region, year, sex]` of recruitment sex ratio.

- WAA, MatAA:

  Arrays `[pop, region, year, season, age, sex]` of weight-at-age and
  maturity-at-age.

- natmort:

  Array `[pop, region, year, age, sex]` of natural mortality at age.

- Movement:

  Array `[pop, region_from, region_to, year, season, age, sex]` of
  movement rates.

- stray_rate:

  Array `[pop, year]` of stray rate.

- sgl_seas_spawning_movement:

  Array `[pop, region_from, region_to, year, age, sex]` of
  single-season-spawning movement rates.

- do_recruits_move:

  Integer (0/1) switch for whether age-1 recruits are subject to
  movement.

- fish_sel, ret_sel:

  Arrays `[pop, region, year, season, age, sex, fish_fleet]` of
  total/retained fishery selectivity.

- dmr:

  Array `[region, year, season, fish_fleet]` of discard mortality rate.

- ZAA:

  Array `[pop, region, year, season, age, sex]` of total mortality at
  age (precomputed).

- NAA, NAA0:

  Arrays `[pop, region, year+1, season, age, sex]`, output containers
  for fished/unfished numbers at age.

- NAA_bef, NAA_aft:

  Arrays `[pop, region, year+1, season, age, sex]`, output containers
  for numbers at age immediately before/after movement.

- Rec:

  Array `[pop, region, year]`, output container for total recruitment
  before seasonal apportionment.

- SSB, Total_Biom, Dynamic_SSB0:

  Arrays `[pop, region, year]`, output containers.

- eff_SSB:

  Array `[pop, year]`, output container for effective
  (natal-homing-adjusted) SSB.

- SR_ref_yr:

  Integer year index supplying the biological inputs, weight at age,
  maturity, natural mortality and movement, to unfished spawning biomass
  per recruit, and so to `S0` and the scale of the stock-recruit curve.
  Default `1`, the first model year, which is what the function used to
  hardcode. Set to `n_yrs` to condition the curve on terminal weight at
  age, which is what several ADMB assessments do; with time-varying
  weight at age the two differ and the whole curve shifts with them. It
  is a year INDEX, not a calendar year, so callers that truncate the
  year dimension (retrospectives) must clamp it.

- growth_mortality_year_fn:

  Optional function of `(y, NAA_y, growth_mortality_state)` called at
  the top of every year with the numbers at age at the start of that
  year, array `[pop, region, age, sex]`, and the state carried from the
  previous year. It returns a list with `state`, carried forward to the
  next call and returned to the caller, and `ZAA_y`, `WAA_y` and
  `MatAA_y`, the year's slices of total mortality, weight and maturity
  at age, which replace those handed in for that year. Passing the state
  in and out keeps the per-year step a function of its arguments.

- growth_mortality_state:

  Initial state for `growth_mortality_year_fn`, carried through the year
  loop and returned as `growth_mortality_state`. Ignored when
  `growth_mortality_year_fn` is `NULL`. This is how cohort growth, whose
  plus group blends by numbers, is evaluated inside the year loop.
  `NULL` (the default) uses the arrays as given.

- n_est_naa_re:

  Number of estimated state-space numbers at age. Zero leaves the
  numbers deterministic. Never inferred from `dim(ln_NAA)`, which is
  non-zero once the setup function has run at all.

- ln_NAA:

  Array `[pop, region, year, age, sex]` of log numbers at age,
  overwriting the deterministic prediction wherever the state is active.

- naa_re_ages, naa_re_yrs:

  Integer index vectors the state is active over.

## Value

List with elements `NAA`, `NAA0`, `NAA_bef`, `NAA_aft`, `Rec`, `SSB`,
`Total_Biom`, `Dynamic_SSB0`, `eff_SSB`, `Aggregated_SSB` (array
`[year]`, SSB summed across pop/region), `Dynamic_Aggregated_SSB0`
(array `[year]`, likewise for `Dynamic_SSB0`), and `NAA_int` (array
`[pop, region, year, season, age, sex]`). `NAA_int` holds the
season-integrated abundance needed by the spatial Baranov catch equation
and is populated only when `move_timing = 2`; it is all zeros otherwise.

## Details

All array arguments matching an output name (`NAA`, `NAA0`, `NAA_bef`,
`NAA_aft`, `Rec`, `SSB`, `Total_Biom`, `Dynamic_SSB0`, `eff_SSB`) are
passed in already dimensioned (typically all-zero, aside from any
initial-year values already inserted upstream) and returned fully
populated over `1:n_yrs`.
