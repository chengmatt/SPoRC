# Draw one at-age observation stream for one region, year, season and fleet

The operating model states an at-age observation the way the estimation
model reads it: summed over whichever of regions and sexes the fleet
reports together, with the fleet's own density and its own standard
deviation. A stream summed over regions is one number, so it is drawn
once, when the region loop reaches region one.

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
  r
)
```

## Arguments

- numbers:

  Array `[n_pop, n_regions, n_ages, n_sexes]` of the quantity at age for
  this year, season and fleet.

- weight:

  Array shaped like `numbers`, read when `use_weight`.

- use:

  Integer array `[n_regions, n_ages, n_sexes]` of use flags.

- se:

  Reported standard errors shaped like `use`.

- ln_sigma:

  Log-scale observation error, `[n_ages, n_sexes]`.

- type_code, like_code, form_code:

  The fleet's aggregation, density and error-source codes.

- use_weight:

  Logical, `TRUE` for an observation in weight.

- r:

  Region the loop is on.

## Value

A list with `true` and `obs`, both `[n_ages, n_sexes]` and `NA` wherever
nothing was drawn.
