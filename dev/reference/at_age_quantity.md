# Quantity at model age an at-age data source is drawn from

What the estimation model's prediction reads (`get_at_age_prediction`):
catch at age, dead discards raised by the discard mortality rate, or the
survey's available numbers, in weight where the fleet reports weight,
for the season given, or summed over every season for a year total.

## Usage

``` r
at_age_quantity(sim_env, source, y, seasons, f, sim)
```

## Arguments

- sim_env:

  Simulation environment, after the year's dynamics.

- source:

  `"catch"`, `"discard"` or `"srv_idx"`.

- y, f, sim:

  Year, fleet and replicate.

- seasons:

  Seasons summed, one for a seasonal data source.

## Value

Array `[n_pop, n_regions, n_ages, n_sexes]`.
