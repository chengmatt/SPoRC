# End of a season for one tag cohort: ageing and recaptures

Survival (under `move_timing = 0`, after movement), the shift to the
next season or the next year and age, the state-space numbers at age
factor, and the predicted recaptures. A season the cohort is at liberty
for none of (`conv_tag_t_tagging = 0` at release) has no recaptures.

## Usage

``` r
end_tag_cohort_season(
  avail_tc,
  recap_tc,
  seas_vals,
  tag_step,
  tag_int,
  y,
  rseas,
  n_seas,
  reporting_y,
  NAA_scalar,
  move_timing,
  n_pop,
  n_regions,
  n_ages,
  n_sexes,
  n_fish_fleets
)
```

## Arguments

- avail_tc, recap_tc:

  Arrays of the cohort's tags
  `[liberty + 1, season, pop, region, age, sex]` and predicted
  recaptures `[liberty, season, pop, region, age, sex, fish_fleet]`.

- seas_vals:

  The cohort's `seas_vals` from `start_tag_cohort_season`.

- tag_step, tag_int:

  The cohort's tags at the end of the season and integrated over it,
  from `move_tag_cohorts` (`NULL` under `move_timing = 0`).

- y, rseas, n_seas:

  Year, season and number of seasons.

- reporting_y:

  Array `[region, fish_fleet]` of tag reporting rates in year `y`.

- NAA_scalar:

  Array `[pop, region, year, season, age, sex]` of the state-space
  numbers at age factor, or `NULL`.

- move_timing:

  Where movement sits relative to mortality (0, 1 or 2).

- n_pop, n_regions, n_ages, n_sexes, n_fish_fleets:

  Dimension sizes.

## Value

List with `avail_tc` and `recap_tc`.
