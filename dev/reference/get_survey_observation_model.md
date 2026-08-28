# Survey observation model

Converts population state into predicted survey index-at-age/length and
survey indices. Called once from the "Survey Observation Model" section
of `SPoRC_rtmb.R`. All array arguments matching an output name (e.g.
`srv_q`, `SrvIAA`) are passed in already dimensioned and returned fully
populated (or, for `srv_sel`, further updated) over `1:n_yrs`.

## Usage

``` r
get_survey_observation_model(
  n_pop,
  n_regions,
  n_yrs,
  n_seas,
  n_srv_fleets,
  n_sexes,
  srv_q_blocks,
  ln_srv_q,
  srv_q,
  do_srv_q_cov,
  srv_q_cov,
  srv_q_coeff,
  srv_selex_type,
  srv_sel,
  srv_sel_l,
  SizeAgeTrans,
  NAA,
  ZAA,
  t_srv,
  SrvIAA,
  fit_lengths,
  SrvIAL,
  srv_idx_type,
  WAA_srv,
  PredSrvIdx,
  Mrate = NULL,
  move_timing = 0,
  seasdur = rep(1, n_seas),
  srv_idx_ages = NULL,
  srv_q_type = NULL,
  ObsSrvIdx = NULL,
  UseSrvIdx = NULL,
  RecDev_anom = NULL,
  do_caal = 0,
  Srv_caal = NULL,
  SizeAgeTrans_srv = NULL,
  srv_len_comp_sel = NULL,
  expm_nsub = 0
)
```

## Arguments

- n_pop, n_regions, n_yrs, n_seas, n_srv_fleets, n_sexes:

  Dimension sizes.

- srv_q_blocks:

  Array `[region, year, srv_fleet]` of time-block catchability indices.

- ln_srv_q:

  Array `[region, block, srv_fleet]` of log survey catchability.

- srv_q:

  Array `[region, year, srv_fleet]`, output container.

- do_srv_q_cov:

  Integer (0/1) switch for a catchability covariate effect.

- srv_q_cov:

  Array `[region, year, srv_fleet, covariate]` of covariate values.

- srv_q_coeff:

  Array `[region, srv_fleet, covariate]` of covariate coefficients.

- srv_selex_type:

  Integer (0 = age-based, 1 = length-based) switch.

- srv_sel:

  Array `[pop, region, year, season, age, sex, srv_fleet]` of survey
  selectivity at age; when `srv_selex_type == 1` this is derived here
  from `srv_sel_l` via `SizeAgeTrans` for `1:n_yrs`.

- srv_sel_l:

  Array `[region, year, len, sex, srv_fleet]` of survey selectivity at
  length.

- SizeAgeTrans:

  Array `[pop, region, year, season, len, age, sex]` age-to-length
  transition matrix.

- NAA:

  Array `[pop, region, year, season, age, sex]` of numbers at age.

- ZAA:

  Array `[pop, region, year, season, age, sex]` of total mortality at
  age.

- t_srv:

  Array `[region, season, srv_fleet]` of survey timing (fraction of
  season elapsed).

- SrvIAA:

  Array `[pop, region, year, season, age, sex, srv_fleet]`, output
  container for survey index at age.

- fit_lengths:

  Integer (0/1) switch for computing length compositions.

- SrvIAL:

  Array `[pop, region, year, season, len, sex, srv_fleet]`, output
  container for survey index at length.

- srv_idx_type:

  Integer vector `[srv_fleet]` selecting abundance/biomass survey index
  type.

- WAA_srv:

  Array `[pop, region, year, season, age, sex, srv_fleet]` of survey
  weight at age.

- PredSrvIdx:

  Array `[pop, region, year, season, srv_fleet]`, output container.

- srv_idx_ages:

  Array `[age, srv_fleet]` of 0/1 weights selecting which ages
  contribute to each fleet's index total, or `NULL` for all ages.
  Restricting to a single age turns that fleet into an index of that age
  alone, and the compositions keep using the full age range because the
  restriction is applied to the index sum rather than to selectivity.

- srv_q_type:

  Integer vector `[srv_fleet]` selecting how catchability is obtained:
  `0` estimated as `exp(ln_srv_q)`, `1` solved analytically as the ratio
  of mean observed to mean predicted, `2` solved analytically on the log
  scale as `exp(mean(log(obs) - log(pred)))`. `NULL` means all
  estimated.

- ObsSrvIdx, UseSrvIdx:

  Arrays `[region, year, season, srv_fleet]` of observed index values
  and their use flags. Required only when a fleet solves catchability
  analytically.

- RecDev_anom:

  Array `[pop, region, deviation]` of recruitment deviations measured
  from the centre their penalty asserts, or `NULL` when no fleet
  observes them. Read only by fleets with `srv_idx_type == 2`.

- do_caal:

  Integer (0/1) switch for building the joint survey index at length and
  age array. Off by default.

- Srv_caal:

  Array `[pop, region, year, season, len, age, sex, srv_fleet]`, output
  container for the survey index at length and age. Only written when
  `do_caal == 1` and `fit_lengths == 1`.

- SizeAgeTrans_srv:

  Array `[pop, region, year, season, len, age, sex, srv_fleet]`, each
  survey's own key from the growth module, or `NULL` to read the shared
  data key.

- srv_len_comp_sel:

  Integer vector `[srv_fleet]`. `0` (default) expands the index at age
  through the key; `1` applies the length selectivity at length to the
  numbers spread over the composition own key, \\I_l = s(l) \sum_a P(l
  \mid a)\\ N_a e^{-t Z_a}\\. Needs `srv_selex_type == 1`.

## Value

List with elements `srv_q`, `srv_sel`, `SrvIAA`, `SrvIAL`, `Srv_caal`,
`PredSrvIdx`.
