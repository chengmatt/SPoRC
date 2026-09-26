# Figure out which entries of the dsem precision can be non-zero

[`RTMB::dgmrf`](https://rdrr.io/pkg/RTMB/man/MVgauss.html) takes the
precision (the inverse of the covariance) as a sparse matrix, and on the
tape only the entries it already stores can be filled, so the positions
are fixed here at setup from the arrows alone.

## Usage

``` r
get_dsem_Q_oo_pattern(IminusB_m, Gamma_m, obs_idx, unobs_idx, has_cov)
```

## Arguments

- IminusB_m:

  Numbered sparse template of \\I - B\\, from
  [`get_dsem_cells`](https://chengmatt.github.io/SPoRC/dev/reference/get_dsem_cells.md),
  read for where the paths sit.

- Gamma_m:

  Numbered sparse template of \\\Gamma\\. Not read: the covariance
  entries only matter here through `has_cov`.

- obs_idx, unobs_idx:

  Cells that keep a deviation of their own, and cells of a series with
  an sd of zero, which are solved out of them.

- has_cov:

  Whether any arrow is a covariance.

## Value

`NULL` when nothing is solved out, since the precision then has \\I -
B\\'s own entries. Otherwise a list with `m`, the numbered sparse
template, and `slot_lin`, where each stored entry sits in the dense
matrix counted down the columns.
