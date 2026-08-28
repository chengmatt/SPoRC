# Write a simulated at-age cell into its true and observed containers

Write a simulated at-age cell into its true and observed containers

## Usage

``` r
store_at_age_cell(sim_env, stream, drawn, r, y, seas, f, sim)
```

## Arguments

- sim_env:

  Environment holding the simulation containers.

- stream:

  Stream tag, e.g. `"CatchAA"`.

- drawn:

  List returned by
  [`sim_at_age_cell`](https://chengmatt.github.io/SPoRC/dev/reference/sim_at_age_cell.md).

- r, y, seas, f, sim:

  Region, year, season, fleet and replicate.

## Value

`invisible(NULL)`, called for its side effect.
