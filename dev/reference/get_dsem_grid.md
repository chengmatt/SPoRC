# Assemble the dsem grid

Builds the dsem matrices, shifts the means by a first year offset, and
sets the cells of a series with an sd of zero from the cells that keep
an innovation.

## Usage

``` r
get_dsem_grid(
  dsem_beta,
  ln_dsem_sd,
  x_grid,
  mu_grid,
  dsem_model,
  dsem_cells,
  delta0 = NULL
)
```

## Arguments

- dsem_beta, ln_dsem_sd, dsem_model:

  As in
  [`get_dsem_arrow_values`](https://chengmatt.github.io/SPoRC/dev/reference/get_dsem_arrow_values.md).

- x_grid:

  Matrix `[year, series]` of the grid, needed when an arrow is
  moderated.

- mu_grid:

  Matrix `[year, series]` of series means.

- dsem_cells:

  Output of
  [`get_dsem_cells`](https://chengmatt.github.io/SPoRC/dev/reference/get_dsem_cells.md).

- delta0:

  Optional numeric vector, one first-year offset per series.

## Value

List with `x_grid` (solved cells set), `mu_grid` (offsets propagated),
`parts` and `solve_mat`.
