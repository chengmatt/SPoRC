# Fishing and natural mortality by tagged fish in one year and season

Builds the mortality-at-age components a tag cohort experiences over a
single season. These depend only on the year and season, not on which
cohort is at liberty, so both
[`get_tagging_observation_model`](https://chengmatt.github.io/SPoRC/dev/reference/get_tagging_observation_model.md)
and the simulation's `generate_fishery_conv_tags_recap` work them out
once per year and season and read them back for every cohort.

## Usage

``` r
get_tag_mort(
  y,
  rseas,
  n_pop,
  n_regions,
  n_ages,
  n_sexes,
  n_fish_fleets,
  use_conv_fish_tagging,
  Fmort,
  fish_sel,
  ret_sel,
  dmr,
  natmort,
  seasdur
)
```

## Arguments

- y:

  Integer year index.

- rseas:

  Integer season index.

- n_pop, n_regions, n_ages, n_sexes, n_fish_fleets:

  Dimension sizes.

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

## Value

List with elements `FAA`, `ret_FAA`, `disc_DAA`, each
`[pop, region, 1, age, sex, fish_fleet]`, and `Z_before_shed`,
`[pop, region, 1, age, sex]`.

## Details

Only fleets flagged in `use_conv_fish_tagging` contribute; the others
are left at zero, since tags are only ever returned by fleets that
report them.

`Z_before_shed` excludes tag shedding, which is a property of the
release cohort rather than of the year and season. Callers add
`exp(ln_conv_tag_shed[tc]) * seasdur[rseas]` to it.

Every array argument is expected without a trailing simulation
dimension, so the operating model slices by `sim` before calling.
