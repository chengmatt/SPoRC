# Start of a season for one tag cohort: mortality and release

Total mortality with the cohort's shedding, the share of the season it
is at liberty for, the release in its release season, and survival.

## Usage

``` r
start_tag_cohort_season(
  avail_tc,
  tagged_fish,
  mort,
  ry,
  rseas,
  tr,
  tseas,
  ln_conv_tag_shed_tc,
  conv_tag_t_tagging_tc,
  ln_init_conv_tag_mort_tc,
  seasdur,
  move_timing,
  n_pop,
  n_ages,
  n_sexes
)
```

## Arguments

- avail_tc:

  Array `[liberty + 1, season, pop, region, age, sex]` of the cohort's
  tags.

- tagged_fish:

  Array `[pop, age, sex]` of fish released, read in the release season
  only.

- mort:

  Output of `get_tag_mort` for this year and season.

- ry, rseas:

  Years at liberty and season.

- tr, tseas:

  Release region and season.

- ln_conv_tag_shed_tc, conv_tag_t_tagging_tc, ln_init_conv_tag_mort_tc:

  The cohort's log shedding rate, share of the release season remaining
  at release, and log tag-induced mortality.

- seasdur:

  Numeric vector `[season]` of season durations.

- move_timing:

  Where movement sits relative to mortality (0, 1 or 2).

- n_pop, n_ages, n_sexes:

  Dimension sizes.

## Value

List with `avail_tc`, the cohort's tags with this season's release in,
and `seas_vals`, what the rest of the season reads: `ry`, `tmp_ZAA`,
`tmp_ret_FAA`, `tmp_SAA`, `tag_frac`, `tag_dur`, `tag_moves` and
`tag_moves_seas`.
