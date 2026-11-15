# Tag cohorts at liberty in a year and season

The cohorts released at or before this season that are still within
their years at liberty.

## Usage

``` r
get_tag_cohorts_at_liberty(
  y,
  rseas,
  conv_tag_release_indicator,
  conv_tag_max_liberty
)
```

## Arguments

- y, rseas:

  Year and season.

- conv_tag_release_indicator:

  Matrix `[cohort, 3]` of release region, year and season.

- conv_tag_max_liberty:

  Integer, the most years a cohort is followed.

## Value

Integer vector of cohort indices.
