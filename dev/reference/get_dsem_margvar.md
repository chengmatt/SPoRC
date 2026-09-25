# The variance of each grid cell given the cells the model is handed

A lognormal deviation with variance \\v\\ has \\E\[\exp(x)\] =
\exp(\mu + v/2)\\, so recruitment keeps its mean at \\R_0\\ only if the
cell's mean drops by \\v/2\\. The \\v\\ that does it is the cell's
variance given what the model is handed (the observed covariate values,
and the rows before a series starts), not the sd line's square: under a
self path \\\rho\\ a settled year has \\\sigma^2 / (1 - \rho^2)\\, and
every lagged path adds to it. For unknown cells \\U\\ and known cells
\\K\\, \\\mathrm{Var}(x_U \mid x_K) = (Q\_{UU})^{-1}\\, the inverse of
the unknown block of the precision, so each cell's variance is a
diagonal entry of that inverse.

## Usage

``` r
get_dsem_margvar(
  dsem_beta,
  ln_dsem_sd,
  x_grid,
  dsem_model,
  dsem_cells,
  known_cell
)
```

## Arguments

- dsem_beta, ln_dsem_sd, dsem_model:

  As in
  [`get_dsem_arrow_values`](https://chengmatt.github.io/SPoRC/dev/reference/get_dsem_arrow_values.md).

- x_grid:

  Matrix `[year, series]` of the grid, needed when an arrow is
  moderated.

- dsem_cells:

  Output of
  [`get_dsem_cells`](https://chengmatt.github.io/SPoRC/dev/reference/get_dsem_cells.md).

- known_cell:

  Logical over the `n_grid_yrs * n_series` cells (years within series),
  `TRUE` where the value is handed to the model.

## Value

Matrix `[year, series]` of variances, zero on the known cells.

## Details

A series with an sd of zero (derived) has no innovation, so its rows
leave the quadratic form; an unknown cell pointing into one is refused,
since the derived value would then be random while the model treats it
as fixed.
