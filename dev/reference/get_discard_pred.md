# Predicted discards over the populations and seasons an observation covers

Discards in numbers or weight add up, so the prediction is the sum of
`PredDiscard` over the populations and seasons the observation covers,
as
[`get_seas_pred`](https://chengmatt.github.io/SPoRC/dev/reference/get_seas_pred.md)
gives. A discard fraction does not add up: two populations each
discarding 30 percent of their catch discard 30 percent together, not
60. So a fraction is the discards over the total catch, each summed over
those populations and seasons first.

## Usage

``` r
get_discard_pred(
  PredDiscard,
  CAA,
  DAA,
  dmr,
  WAA_fish,
  units,
  pops,
  r,
  y,
  seasons,
  f
)
```

## Arguments

- PredDiscard:

  Prediction array `[pop, region, year, season, fleet]`.

- CAA, DAA:

  Retained and dead discarded catch at age
  `[pop, region, year, season, age, sex, fleet]`.

- dmr:

  Discard mortality rate `[region, year, season, fleet]`, which raises
  dead discards to all discards.

- WAA_fish:

  Fishery weight at age, shaped like `CAA`, read for a fraction of
  weight.

- units:

  The fleet's `discard_units`: 0 numbers, 1 weight, 2 fraction of
  numbers, 3 fraction of weight.

- pops:

  Populations summed: every one for a regional observation, one for a
  population-specific one.

- r, y, f:

  Region, year and fleet of the observation.

- seasons:

  Seasons summed: the observation's own, or every season for a year
  total.

## Value

Scalar prediction.
