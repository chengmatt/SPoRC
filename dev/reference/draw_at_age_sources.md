# Draw one fleet's at-age data sources for a region, year and season

The aggregated data source, then the population-specific one, each from
[`at_age_quantity`](https://chengmatt.github.io/SPoRC/dev/reference/at_age_quantity.md):
the season's, or every season summed when the fleet's observation is a
year total (`_seas_Type`), which is drawn in the last season, once every
season is formed, into the season its use flags name. A
population-specific source draws each population from its own numbers.

## Usage

``` r
draw_at_age_sources(sim_env, source, y, seas, r, f, sim, bias_correct_oe = 0)
```

## Arguments

- sim_env:

  Simulation environment.

- source:

  `"catch"`, `"discard"` or `"srv_idx"`.

- y, seas, r, f, sim:

  Year, season, region, fleet and replicate.

- bias_correct_oe:

  Whether lognormal draws sit at their mean.

## Value

`invisible(NULL)`; `sim_env` is modified in place.
