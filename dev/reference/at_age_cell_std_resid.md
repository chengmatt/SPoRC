# One cell's standardized at-age residuals, or NULL to draw them in the cycle

One cell's standardized at-age residuals, or NULL to draw them in the
cycle

## Usage

``` r
at_age_cell_std_resid(sim_env, data_source, r, y, seas, f, sim, p = NULL)
```

## Arguments

- sim_env:

  Simulation environment.

- data_source:

  `"CatchAA"`, `"DiscardAA"` or `"SrvIdxAA"`.

- r, y, seas, f, sim:

  Region, year, season, fleet and replicate.

- p:

  Population, for a population-specific data source. `NULL` (default)
  for the aggregated one.

## Value

`[n_obs_ages, n_sexes]`, or `NULL` for a fleet drawn age by age.
