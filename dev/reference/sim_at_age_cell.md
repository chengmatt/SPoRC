# Draw one at-age data source for one region, year, season and fleet

The operating model states an at-age observation the way the estimation
model reads it: summed over whichever of regions and sexes the fleet
reports together, read through the fleet's ageing error onto the
observed ages, with the fleet's own density and its own standard
deviation. A data source summed over regions is one number, so it is
drawn once, when the region loop reaches region one.

## Usage

``` r
sim_at_age_cell(
  numbers,
  weight,
  use,
  se,
  ln_sigma,
  type_code,
  like_code,
  form_code,
  use_weight,
  r,
  ageing_error = NULL,
  bias_correct_oe = 0,
  std_resid = NULL
)
```

## Arguments

- numbers:

  Array `[n_pop, n_regions, n_ages, n_sexes]` of the quantity at model
  age for this year, season and fleet.

- weight:

  Array shaped like `numbers`, read when `use_weight`.

- use:

  Integer array `[n_regions, n_obs_ages, n_sexes]` of use flags.

- se:

  Reported standard errors shaped like `use`.

- ln_sigma:

  Log-scale observation error, `[n_obs_ages, n_sexes]`.

- type_code, like_code, form_code:

  The fleet's aggregation, density and error-source codes.

- use_weight:

  Logical, `TRUE` for an observation in weight.

- r:

  Region the loop is on.

- ageing_error:

  Matrix `[n_ages, n_obs_ages]` reading model ages as observed ages for
  this year and fleet, or `NULL` for the identity.

- bias_correct_oe:

  `1` draws a lognormal observation at its mean, as the estimation model
  reads it under the same setting; `0` (default) at its median.

- std_resid:

  Standardized residuals `[n_obs_ages, n_sexes]` drawn up front by
  [`draw_sim_at_age_corr`](https://chengmatt.github.io/SPoRC/dev/reference/draw_sim_at_age_corr.md)
  for a fleet whose ages are correlated, each scaled here by its age's
  sd. `NULL` (default) draws each age on its own.

## Value

A list with `true` and `obs`, both `[n_obs_ages, n_sexes]` and `NA`
wherever nothing was drawn.
