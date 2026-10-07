# Year's discards of a fleet in one region

Total discards, or the discarded fraction of the total catch, summed
over every season and the populations given before any ratio is taken,
so a fraction is the year's discards over the year's catch rather than a
sum of seasonal fractions. Mirrors `get_discard_pred` in the estimation
model.

## Usage

``` r
year_discard_total(sim_env, r, y, f, sim, pops)
```

## Arguments

- sim_env:

  Simulation environment, after the year's dynamics.

- r, y, f, sim:

  Region, year, fleet and replicate.

- pops:

  Populations summed: all of them for the regional data source, one for
  the population-specific one.

## Value

Scalar, in the fleet's `discard_units`.
