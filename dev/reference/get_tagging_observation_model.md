# Tagging observation model

Forward-simulates conventional tag cohort availability (release,
movement, mortality/ageing) and predicted recaptures across all
recapture years, seasons, and cohorts. Called once from the "Tagging
Observation Model" section of `SPoRC_rtmb.R`, guarded by
`any(use_conv_fish_tagging == 1)`. Depends on population state
(`fish_sel`, `ret_sel`, `srv_sel`, `NAA_bef`, `Movement`, `natmort`,
`Fmort`, `dmr`) already computed upstream by the population projection
and selectivity sections.

## Usage

``` r
get_tagging_observation_model(
  n_fish_fleets,
  n_regions,
  n_conv_tag_cohorts,
  n_yrs,
  n_seas,
  n_pop,
  n_ages,
  n_sexes,
  conv_tag_fish_reporting_blocks,
  conv_tag_fish_reporting_pars,
  conv_tag_fish_reporting,
  conv_tag_release_indicator,
  conv_tag_max_liberty,
  use_conv_fish_tagging,
  Fmort,
  fish_sel,
  ret_sel,
  dmr,
  natmort,
  seasdur,
  ln_conv_tag_shed,
  conv_tag_t_tagging,
  conv_tagged_fish,
  conv_fish_tag_attr,
  conv_tag_release_platform,
  srv_sel,
  NAA_bef,
  ln_init_conv_tag_mort,
  do_recruits_move,
  Movement,
  conv_tag_fish_avail,
  pred_conv_tag_fish_recap,
  Mrate = NULL,
  move_timing = 0,
  expm_nsub = 0,
  NAA_scalar = NULL
)
```

## Arguments

- n_fish_fleets, n_regions, n_conv_tag_cohorts, n_yrs, n_seas, n_pop,
  n_ages, n_sexes:

  Dimension sizes.

- conv_tag_fish_reporting_blocks:

  Array `[region, block, fish_fleet]` mapping years to reporting-rate
  blocks.

- conv_tag_fish_reporting_pars:

  Array `[region, block, fish_fleet]` of tag reporting rate parameters
  on the logit scale.

- conv_tag_fish_reporting:

  Array `[region, year, fish_fleet]`, output container for reporting
  rate on the natural scale.

- conv_tag_release_indicator:

  Matrix `[cohort, 3]` of release region/year/season per tag cohort.

- conv_tag_max_liberty:

  Integer maximum years at liberty tracked.

- use_conv_fish_tagging:

  Integer vector `[fish_fleet]` (0/1) flagging fleets with tagging data.

- Fmort:

  Array `[region, year, season, fish_fleet]` of fishing mortality.

- fish_sel, ret_sel:

  Arrays `[pop, region, year, season, age, sex, fish_fleet]` of
  total/retained fishery selectivity.

- dmr:

  Array `[region, year, season, fish_fleet]` of discard mortality rate.

- natmort:

  Array `[pop, region, year, season, age, sex]` of natural mortality at
  age.

- seasdur:

  Numeric vector `[season]` of season duration (fraction of year).

- ln_conv_tag_shed:

  Numeric vector `[cohort]` of log tag shedding rate.

- conv_tag_t_tagging:

  Numeric vector `[cohort]` of within-season timing of tag release (1 =
  start of season).

- conv_tagged_fish:

  Array `[cohort, age, sex, ...]` of raw tag release numbers, passed
  through to `release_conv_tag_attr`.

- conv_fish_tag_attr, conv_tag_release_platform:

  Arguments passed through to `release_conv_tag_attr` controlling how
  releases are apportioned across dimensions.

- srv_sel:

  Array `[pop, region, year, season, age, sex, srv_fleet]` of survey
  selectivity, used by `release_conv_tag_attr` when releases are
  platform-apportioned via survey gear.

- NAA_bef:

  Array `[pop, region, year, season, age, sex]` of numbers at age before
  movement, used by `release_conv_tag_attr` to apportion releases.

- ln_init_conv_tag_mort:

  Numeric vector `[cohort]` of log initial tag-induced mortality rate.

- do_recruits_move:

  Integer (0/1) switch for whether age-1 fish are subject to movement.

- Movement:

  Array `[pop, region_from, region_to, year, season, age, sex]` of
  movement rates.

- conv_tag_fish_avail:

  Array `[liberty+1, season, cohort, pop, region, age, sex]`, output
  container for tags available for recapture.

- pred_conv_tag_fish_recap:

  Array `[liberty, season, cohort, pop, region, age, sex, fish_fleet]`,
  output container for predicted recaptures.

- NAA_scalar:

  Array `[pop, region, year, season, age, sex]` of the factor the
  state-space numbers at age applied to the deterministic prediction,
  one wherever it did not apply. Tagged fish are a subset of the
  population and the innovation reads as unmodelled mortality, so the
  cohorts take the same factor at every boundary the state acts on,
  within a year as well as across one. `NULL` (the default) leaves them
  on the deterministic trajectory, which is correct only when the state
  is off.

## Value

List with elements `conv_tag_fish_reporting`, `conv_tag_fish_avail`,
`pred_conv_tag_fish_recap`.

## Details

`conv_tag_fish_reporting`, `conv_tag_fish_avail`, and
`pred_conv_tag_fish_recap` are passed in already dimensioned (typically
all-zero) and returned fully populated.
