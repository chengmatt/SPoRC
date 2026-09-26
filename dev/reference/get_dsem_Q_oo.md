# The precision of the cells that have an innovation

With nothing solved out this is \\(I - B)^{\top}V^{-1}(I - B)\\ over
every cell. Once a series has an sd of zero \\V\\ is singular and that
product does not exist. Substituting the solved cells into the remaining
rows leaves \\S(x_o - \mu_o) = \varepsilon_o\\ with \$\$S = (I -
B)\_{oo} - (I - B)\_{ou}(I - B)\_{uu}^{-1}(I - B)\_{uo},\$\$ so those
cells are normal with precision \\S^{\top}V\_{oo}^{-1}S\\. Building it
rather than reading the density off the innovations is what accounts for
the determinant of \\S\\, which is not one once the arrows loop within a
year.

## Usage

``` r
get_dsem_Q_oo(IminusB, parts, solve_mat, dsem_cells)
```

## Arguments

- IminusB:

  Sparse \\I - B\\, from
  [`get_dsem_matrices`](https://chengmatt.github.io/SPoRC/dev/reference/get_dsem_matrices.md).

- parts:

  The rest of that output, read for \\V^{-1}\\, \\V\\ or the sds.

- solve_mat:

  Output of
  [`get_dsem_solve_mat`](https://chengmatt.github.io/SPoRC/dev/reference/get_dsem_solve_mat.md).

- dsem_cells:

  Output of
  [`get_dsem_cells`](https://chengmatt.github.io/SPoRC/dev/reference/get_dsem_cells.md).

## Value

Sparse precision over the cells that keep an innovation.

## Details

The middle term is \\\tilde{V}\_{oo} = V\_{oo} + C V\_{uo} +
V\_{ou}C^{\top}\\ with \\C = M\_{ou}M\_{uu}^{-1}\\, the two cross terms
covering a covariance line between a solved cell and one that keeps its
innovation. That line is refused by
[`read_dsem_arrows`](https://chengmatt.github.io/SPoRC/dev/reference/read_dsem_arrows.md),
so both terms are always zero and only \\V\_{oo}\\ is formed.
