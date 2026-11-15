# Movement for every tag cohort at liberty in a season

Every cohort at liberty goes through the season's movement matrix, or
its seasonal operator under `move_timing = 2`, through
`move_tag_vectors`. Under `move_timing = 0` the cohorts' tags are moved
in place, and mortality follows in `end_tag_cohort_season`. Under
`move_timing = 1` each cohort is survived and then moved. Under
`move_timing = 2` the cohorts of a shedding group with the same time at
liberty share one seasonal operator and its integral.

## Usage

``` r
move_tag_cohorts(
  tag_avail_cohort,
  cohort_seas_vals,
  cohorts_at_liberty,
  rseas,
  move_seas,
  rate_seas,
  tag_shed_group,
  move_timing,
  do_recruits_move,
  expm_nsub,
  n_pop,
  n_regions,
  n_ages,
  n_sexes
)
```

## Arguments

- tag_avail_cohort:

  List `[cohort]` of each cohort's tags.

- cohort_seas_vals:

  List `[cohort]` of `seas_vals` from `start_tag_cohort_season`.

- cohorts_at_liberty:

  Integer vector of the cohorts at liberty this season.

- rseas:

  Season.

- move_seas:

  Array `[pop, region_from, region_to, age, sex]`, the season's movement
  matrices.

- rate_seas:

  Array `[pop, region_from, region_to, age, sex]`, the season's movement
  rates under `move_timing = 2`, else `NULL`.

- tag_shed_group:

  Integer vector `[cohort]` from `get_tag_shed_group`.

- move_timing, do_recruits_move, expm_nsub:

  Movement settings.

- n_pop, n_regions, n_ages, n_sexes:

  Dimension sizes.

## Value

List with `tag_avail_cohort` (moved under `move_timing = 0`), and
`tag_step_cohort` and `tag_int_cohort`, each cohort's tags at the end of
the season and integrated over it under `move_timing` 1 and 2.
