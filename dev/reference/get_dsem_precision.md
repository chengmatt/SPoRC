# The precision of the whole grid

\\Q = (I - B)^{\top} V^{-1} (I - B)\\, the inverse covariance of the
stacked grid: if the innovations have covariance \\V\\ and the grid is
\\x - \mu = (I - B)^{-1}\varepsilon\\, then the grid's covariance is
\\(I - B)^{-1} V (I - B)^{-\top}\\ and its inverse is this product. It
is the same form
[`Get_3d_precision`](https://chengmatt.github.io/SPoRC/dev/reference/Get_3d_precision.md)
builds for the numbers at age field.
[`RTMB::dgmrf`](https://rdrr.io/pkg/RTMB/man/MVgauss.html) evaluates the
density from it without forming the covariance.

## Usage

``` r
get_dsem_precision(
  dsem_beta,
  ln_dsem_sd,
  dsem_model,
  dsem_cells,
  x_grid = NULL
)
```

## Arguments

- dsem_beta, ln_dsem_sd, dsem_model:

  As in
  [`get_dsem_arrow_values`](https://chengmatt.github.io/SPoRC/dev/reference/get_dsem_arrow_values.md).

- dsem_cells:

  Output of
  [`get_dsem_cells`](https://chengmatt.github.io/SPoRC/dev/reference/get_dsem_cells.md).

- x_grid:

  Matrix `[year, series]` of the grid, needed when an arrow is
  moderated.

## Value

Sparse precision matrix over the `n_grid_yrs * n_series` cells.
