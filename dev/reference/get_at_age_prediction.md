# Predicted value for one age-disaggregated observation

Catch and discards are already at age and only need their units applied.
The indices apply an age-specific catchability to the numbers available
to that fleet.

## Usage

``` r
get_at_age_prediction(source, pop, p, r, y, seas, a, f, arrays)
```

## Arguments

- source:

  Character, one of `"catch"`, `"discard"`, `"fish_index"` or
  `"srv_index"`.

- pop:

  Logical. `TRUE` for the population-specific stream, in which case `p`
  indexes a single population rather than summing over all.

- p, r, y, seas, a, f:

  Population, region, year, season, age and fleet indices.

- arrays:

  Named list of the model arrays the prediction reads: `CAA`, `DAA`,
  `SrvIAA`, `FishIAA`, `WAA_fish` and `catch_units`. The two index
  arrays already carry their fleet's selectivity, timing and movement
  treatment, and the age shape of catchability lives in that
  selectivity: a fleet fit age by age uses the `"nonparfree"`
  selectivity form, whose values carry the height of the curve as well
  as its shape.

## Value

The predicted observation, a scalar.
