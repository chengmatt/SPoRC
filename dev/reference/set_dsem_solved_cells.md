# Set the solved cells from the cells that keep an innovation

Set the solved cells from the cells that keep an innovation

## Usage

``` r
set_dsem_solved_cells(x_grid, mu_grid, solve_mat, dsem_cells)
```

## Arguments

- x_grid:

  Matrix `[year, series]`, the solved cells overwritten.

- mu_grid:

  Matrix `[year, series]` of series means.

- solve_mat:

  Output of
  [`get_dsem_solve_mat`](https://chengmatt.github.io/SPoRC/dev/reference/get_dsem_solve_mat.md).

- dsem_cells:

  Output of
  [`get_dsem_cells`](https://chengmatt.github.io/SPoRC/dev/reference/get_dsem_cells.md).

## Value

`x_grid` with every solved cell set.
