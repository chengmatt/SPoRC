# Predicted value for one age-disaggregated observation

Catch and discards are already at age and only need their units applied.
The discards are dead discards, so they are raised by the discard
mortality rate to the total the observation counts, exactly as the
aggregated data source does. The survey index applies an age-specific
catchability to the numbers available to that fleet.

## Usage

``` r
get_at_age_prediction(source, arrays, p_idx, r_idx, s_idx, y, seas, a, f)
```

## Arguments

- source:

  Character, one of `"catch"`, `"discard"` or `"srv_index"`.

- arrays:

  Named list of the model arrays the prediction reads: `CAA`, `DAA`,
  `SrvIAA`, `WAA_fish`, `dmr`, `catch_units` and `discard_units`. The
  two index arrays already have their fleet's selectivity, timing and
  movement treatment, and the age shape of catchability lives in that
  selectivity: a fleet fit age by age uses the `"nonparfree"`
  selectivity form, whose values hold the height of the curve as well as
  its shape.

- p_idx, r_idx, s_idx:

  Population, region and sex indices, each either one index or the whole
  extent of that dim.

- y, seas, a, f:

  Year, season, age and fleet indices.

## Value

The predicted observation, a scalar.

## Details

Population, region and sex arrive as index vectors rather than single
indices. A dim the fleet splits over is a single index, and a dim it
sums over is the whole extent, so one expression covers every
aggregation.
