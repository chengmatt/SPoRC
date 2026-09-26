# Grid cells each dsem arrow lands on

Cells of the year by series grid are numbered series by series, years
within a series, which is the order
[`as.vector()`](https://rdrr.io/r/base/vector.html) lays out a
`[year, series]` matrix. A path from series `j` to series `k` at lag `L`
fills cell `(t, k)` of row and cell `(t - L, j)` of column in \\B\\, for
every year `t > L`.

## Usage

``` r
get_dsem_cells(dsem_model, n_grid_yrs)
```

## Arguments

- dsem_model:

  Output of `read_dsem_arrows`.

- n_grid_yrs:

  Number of years in the grid.

## Value

List with `n_cells`, `n_grid_yrs` and, for \\I - B\\ and \\\Gamma\\, a
numbered sparse template (`m`), the entry filling each stored slot
(`slot_entry`), the arrow behind each entry (`entry_arrow`, 0 for the
diagonal of \\I - B\\), the row and column of each entry (`entry_row`,
`entry_col`) and its year (`entry_yr`), which a moderated arrow reads
its value at. `det_is_one` says whether the cells can be ordered so that
everything an arrow comes from is set before what it points to,
`needs_dense_logdet` whether dgmrf will hold a random effect in the
precision, and, for a reduced rank model, `project_k` with the cells of
a solved series, `obs_idx` and `unobs_idx` splitting the cells the way
dsem does, and `Q_oo` the pattern
[`get_dsem_Q_oo`](https://chengmatt.github.io/SPoRC/dev/reference/get_dsem_Q_oo.md)
fills.
