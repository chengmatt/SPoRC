# Work out the time-series where the sd is fixed at zero

A series with an sd of zero has no innovation, so its value each year is
exactly what the arrows pointing into it give: its own mean plus, for
every path into it, that arrow's coefficient times how far the series it
reads sits from its mean in the year it reads. Series are visited in an
order that sets each one before anything reads it (`series_order`), one
year at a time.

## Usage

``` r
fill_dsem_derived(x_grid, mu_grid, arrow_value, dsem_model, only = NULL)
```

## Arguments

- x_grid:

  Matrix `[year, series]`, the derived columns overwritten.

- mu_grid:

  Matrix `[year, series]` of series means.

- arrow_value:

  Output of
  [`get_dsem_arrow_values`](https://chengmatt.github.io/SPoRC/dev/reference/get_dsem_arrow_values.md).

- dsem_model:

  Output of
  [`read_dsem_arrows`](https://chengmatt.github.io/SPoRC/dev/reference/read_dsem_arrows.md).

- only:

  Integer vector of series to fill, or `NULL` (default) for every
  derived series. Used to fill the catchability series before the
  observations are built, since the rest need the population first.

## Value

`x_grid` with every derived cell set.
