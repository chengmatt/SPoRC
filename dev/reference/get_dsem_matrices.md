# Get DSEM matrices

Stack the grid into one single vector, years within series. Every path
arrow puts its coefficient into \\B\\ at the cell it points to (row) and
the cell it reads (column, the same series `lag` rows earlier), so
\\(I - B)(x - \mu)\\ turns the grid into its innovations. Every sd line
puts its value on the diagonal of \\\Gamma\\ and every covariance line
off it, so \\V = \Gamma^{\top}\Gamma\\ is the covariance of the
innovations in one year. The positions are worked out once, outside of
the tape, by
[`get_dsem_cells`](https://chengmatt.github.io/SPoRC/dev/reference/get_dsem_cells.md);
this function only writes the numbers into the stored slots, which is
what lets the matrices be built on the tape.

## Usage

``` r
get_dsem_matrices(
  dsem_beta,
  ln_dsem_sd,
  dsem_model,
  dsem_cells,
  x_grid = NULL,
  need_Vinv = TRUE
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

- need_Vinv:

  Whether to form \\V^{-1}\\; the derived-series route only needs the
  sds.

## Value

List with `IminusB` (sparse, \\I - B\\), `Vinv` (sparse \\V^{-1}\\, or
`NULL` when not asked for) and `sd_cell` (each cell's innovation sd, or
`NULL` when a covariance line couples them).

## Details

A moderated arrow has no single value: its coefficient in year \\t\\ is
the moderating series' value in that year, read from `x_grid`, and a
moderated sd is that value or its exponential (`mod_var_logscale`).
Also, under `variance = "diagonal"` or `"marginal"` (set in
[`read_dsem_arrows`](https://chengmatt.github.io/SPoRC/dev/reference/read_dsem_arrows.md))
an sd line is the marginal sd of its series rather than the innovation
sd.
