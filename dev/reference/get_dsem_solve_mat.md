# Solve the cells with no innovation out of the ones that keep it

A series with an sd of zero has no innovation, so its rows of \\(I -
B)(x - \mu) = \varepsilon\\ read zero. Writing \\u\\ for its cells
(`unobs_idx`) and \\o\\ for the rest (`obs_idx`), those rows give every
one of those cells from the others, \$\$x_u - \mu_u = -(I -
B)\_{uu}^{-1}(I - B)\_{uo}(x_o - \mu_o),\$\$ and this returns the matrix
in the middle. It is kept sparse: a solved cell is set only by the cells
an arrow reaches it from, so most entries are zero. When no arrow runs
between two solved cells, \\(I - B)\_{uu}\\ is the identity and the
matrix is \\(I - B)\_{uo}\\ itself.

## Usage

``` r
get_dsem_solve_mat(IminusB, dsem_cells)
```

## Arguments

- IminusB:

  Sparse \\I - B\\, from
  [`get_dsem_matrices`](https://chengmatt.github.io/SPoRC/dev/reference/get_dsem_matrices.md).

- dsem_cells:

  Output of
  [`get_dsem_cells`](https://chengmatt.github.io/SPoRC/dev/reference/get_dsem_cells.md).

## Value

Sparse matrix, one row per solved cell and one column per cell that
keeps an innovation, or `NULL` when nothing is solved out.
