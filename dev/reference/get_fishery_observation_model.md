# Fishery observation model

Converts realized fishing mortality and population state into predicted
catch-at-age/length, discards, and fishery indices. Called once from the
"Fishery Observation Model" section of `SPoRC_rtmb.R`. All array
arguments matching an output name (e.g. `fish_q`, `CAA`) are passed in
already dimensioned (typically all-zero) and returned fully populated
over `1:n_yrs`.

## Usage

``` r
get_fishery_observation_model(
  n_pop,
  n_regions,
  n_yrs,
  n_seas,
  n_fish_fleets,
  n_sexes,
  fish_q_blocks,
  ln_fish_q,
  fish_q,
  ret_FAA,
  disc_FAA,
  ZAA,
  NAA,
  CAA,
  DAA,
  CAL,
  DAL,
  PredCatch,
  PredDiscard,
  PredFishIdx,
  fit_lengths,
  SizeAgeTrans,
  catch_units,
  discard_units,
  WAA_fish,
  dmr,
  fish_idx_type,
  fish_sel,
  ret_sel,
  Mrate = NULL,
  move_timing = 0,
  seasdur = rep(1, n_seas),
  NAA_int = NULL,
  t_fish = NULL,
  fish_idx_ages = NULL,
  fish_q_type = NULL,
  do_fish_q_cov = 0,
  fish_q_cov = NULL,
  fish_q_coeff = NULL,
  ObsFishIdx = NULL,
  UseFishIdx = NULL
)
```

## Arguments

- n_pop, n_regions, n_yrs, n_seas, n_fish_fleets, n_sexes:

  Dimension sizes.

- fish_q_blocks:

  Array `[region, year, fish_fleet]` of time-block catchability indices.

- ln_fish_q:

  Array `[region, block, fish_fleet]` of log fishery catchability.

- fish_q:

  Array `[region, year, fish_fleet]`, output container.

- ret_FAA, disc_FAA:

  Arrays `[pop, region, year, season, age, sex, fish_fleet]` of
  retained/discarded fishing mortality at age.

- ZAA:

  Array `[pop, region, year, season, age, sex]` of total mortality at
  age.

- NAA:

  Array `[pop, region, year, season, age, sex]` of numbers at age.

- CAA, DAA:

  Arrays `[pop, region, year, season, age, sex, fish_fleet]`, output
  containers for retained/discarded catch at age.

- CAL, DAL:

  Arrays `[pop, region, year, season, len, sex, fish_fleet]`, output
  containers for retained/discarded catch at length.

- PredCatch, PredDiscard, PredFishIdx:

  Arrays `[pop, region, year, season, fish_fleet]`, output containers.

- fit_lengths:

  Integer (0/1) switch for computing length compositions.

- SizeAgeTrans:

  Array `[pop, region, year, season, len, age, sex]` age-to-length
  transition matrix.

- catch_units, discard_units:

  Integer vectors `[fish_fleet]` selecting abundance/biomass/fraction
  units.

- WAA_fish:

  Array `[pop, region, year, season, age, sex, fish_fleet]` of fishery
  weight at age.

- dmr:

  Array `[region, year, season, fish_fleet]` of discard mortality rate.

- fish_idx_type:

  Integer vector `[fish_fleet]` selecting abundance/biomass fishery
  index type.

- fish_sel, ret_sel:

  Arrays `[pop, region, year, season, age, sex, fish_fleet]` of
  total/retained fishery selectivity.

- t_fish:

  Array `[region, season, fish_fleet]` of fishery index timing, the
  fraction of the season elapsed when the index is observed. Numbers at
  age are decayed by `exp(-t_fish * ZAA)` before the index is formed,
  mirroring `t_srv` for surveys. `NULL` (the default) skips the decay
  entirely and reproduces a start-of-season index.

## Value

List with elements `fish_q`, `CAA`, `DAA`, `CAL`, `DAL`, `PredCatch`,
`PredDiscard`, `PredFishIdx`.
