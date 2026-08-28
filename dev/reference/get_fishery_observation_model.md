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
  FishIAA = NULL,
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
  UseFishIdx = NULL,
  do_caal = 0,
  Fish_caal = NULL,
  Fish_caal_discard = NULL,
  SizeAgeTrans_fish = NULL,
  fish_len_comp_sel = NULL,
  fish_selex_type = 0,
  ret_selex_type = 0,
  fish_sel_l = NULL,
  ret_sel_l = NULL,
  Fmort = NULL,
  expm_nsub = 0
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

- FishIAA:

  Container for the fishery index numbers at age, defaulting to a zero
  array shaped like `CAA`. Filled with the fleet's timing and movement
  treatment applied so the at-age index likelihood reads the same
  quantity the aggregated one does.

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

- do_caal:

  Integer (0/1) switch for building the joint catch at length and age
  arrays. Off by default.

- Fish_caal, Fish_caal_discard:

  Arrays `[pop, region, year, season, len, age, sex, fish_fleet]`,
  output containers for retained/discarded catch at length and age. Only
  written when `do_caal == 1` and `fit_lengths == 1`.

- SizeAgeTrans_fish:

  Array `[pop, region, year, season, len, age, sex, fish_fleet]`, each
  fleet's own key from the growth module, or `NULL` to read the shared
  data key.

- fish_len_comp_sel:

  Integer vector `[fish_fleet]`. `0` (default) selects the catch at age
  and spreads it over lengths afterwards, so every fish of an age is
  equally catchable. `1` selects at length instead: the fish at each age
  are spread over the lengths of the composition key and then selected
  length by length, \\C_l = s(l) \sum_a P(l \mid a)\\ N_a (1 - e^{-Z_a})
  F / Z_a\\, so the long fish of an age are taken more often. Requires
  `fish_selex_type == 1`.

- fish_selex_type, ret_selex_type:

  Integer (0 age, 1 length), the scale total and retention selectivity
  are defined on.

- fish_sel_l, ret_sel_l:

  Arrays `[region, year, len, sex, fish_fleet]` of selectivity at
  length, read under `fish_len_comp_sel == 1`.

- Fmort:

  Array `[region, year, season, fish_fleet]` of fishing mortality, read
  under `fish_len_comp_sel == 1`.

## Value

List with elements `fish_q`, `CAA`, `DAA`, `CAL`, `DAL`, `Fish_caal`,
`Fish_caal_discard`, `PredCatch`, `PredDiscard`, `PredFishIdx`.
