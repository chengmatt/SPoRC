# Enumerate valid conventional-tag recovery events

Enumerates all valid `(tc, ry, rseas)` recovery events applying the
exact skip logic used in the fitting loop (release season and mixing
period). Returns a `data.frame` with one row per event, in loop order.

## Usage

``` r
tag_grid(
  n_conv_tag_cohorts,
  conv_tag_release_indicator,
  conv_tag_max_liberty,
  n_yrs,
  n_seas,
  conv_tag_mixing_period
)
```

## Arguments

- n_conv_tag_cohorts:

  Number of conventional tag cohorts.

- conv_tag_release_indicator:

  Matrix giving release region, year, season for each cohort.

- conv_tag_max_liberty:

  Maximum years at liberty to evaluate.

- n_yrs:

  Total number of modeled years.

- n_seas:

  Number of seasons per year.

- conv_tag_mixing_period:

  Minimum seasons at liberty before tags are modeled.

## Value

A `data.frame` with columns `tc`, `ry`, `rseas`, `tr`, `ty`, and
`tseas`, in the same order the fitting loop would visit them. If no
valid events exist, returns an empty data.frame.
