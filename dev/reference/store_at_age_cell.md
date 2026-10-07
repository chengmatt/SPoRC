# Write a simulated at-age cell into its true and observed containers

Write a simulated at-age cell into its true and observed containers

## Usage

``` r
store_at_age_cell(sim_env, data_source, drawn, r, y, seas, f, sim, p = NULL)
```

## Arguments

- sim_env:

  Environment holding the simulation containers.

- data_source:

  Data source tag, e.g. `"CatchAA"`.

- drawn:

  List returned by
  [`sim_at_age_cell`](https://chengmatt.github.io/SPoRC/dev/reference/sim_at_age_cell.md).

- r, y, seas, f, sim:

  Region, year, season, fleet and replicate.

- p:

  Population, for a population-specific data source, whose containers
  have a leading population dim. `NULL` (default) for the aggregated
  one.

## Value

`invisible(NULL)`, called for its side effect.
