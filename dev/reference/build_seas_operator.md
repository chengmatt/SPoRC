# Build a seasonal transition operator combining movement and survival

Returns the linear operator that advances a numbers-at-age vector across
one season, folding together movement and total mortality according to
`move_timing`. The operator is returned in the package's row-vector
convention, i.e. `N_new = t(N) %*% T` where `T[from, to]`, matching the
storage convention of `Movement`.

## Usage

``` r
build_seas_operator(Move, Z, Q = NULL, dur = 1, move_timing = 0, expm_nsub = 0)
```

## Arguments

- Move:

  Square `[n_regions x n_regions]` movement matrix in row convention
  (`Move[from, to]`), as stored in `Movement`. Required for
  `move_timing` 0 and 1; ignored for `move_timing = 2`.

- Z:

  Numeric vector of length `n_regions` giving total mortality for this
  season, already scaled by season duration (i.e. the same quantity
  stored in `ZAA`). Pass zeros for a movement-only step.

- Q:

  Square `[n_regions x n_regions]` instantaneous rate matrix (generator)
  in row convention, as stored in `Mrate`. Required when
  `move_timing = 2`; ignored otherwise.

- dur:

  Season duration used to scale `Q`. Should be `seasdur[seas]`. Only
  used when `move_timing = 2`.

- move_timing:

  Integer flag for the movement/mortality ordering: `0` = movement then
  mortality (default, historical SPoRC behaviour), `1` = mortality then
  movement, `2` = continuous (simultaneous) movement and mortality.

- expm_nsub:

  Integer controlling how the matrix exponential is evaluated under
  `move_timing = 2`: `0` (default) uses
  [`Matrix::expm`](https://rdrr.io/pkg/Matrix/man/expm-methods.html), a
  value \\n \ge 1\\ uses the implicit backward Euler scheme \\(I -
  A/n)^{-n}\\. See
  [`mat_exp`](https://chengmatt.github.io/SPoRC/dev/reference/mat_exp.md).

## Value

A square `[n_regions x n_regions]` matrix in row convention.

## Details

Writing \\s = \exp(-Z)\\ and letting \\M\\ be the row-convention
movement matrix, the three operators are \$\$T_0 = M \\
\mathrm{diag}(s), \qquad T_1 = \mathrm{diag}(s) \\ M, \qquad T_2 =
\left\[\exp\left(Q^\top \Delta -
\mathrm{diag}(Z)\right)\right\]^\top\$\$ where \\\Delta\\ is `dur`.
\\Q^\top\\ converts the stored row-convention generator back to the
column convention the exponential is taken in, which is
[`Matrix::expm`](https://rdrr.io/pkg/Matrix/man/expm-methods.html) or
the implicit solve of
[`mat_exp`](https://chengmatt.github.io/SPoRC/dev/reference/mat_exp.md)
according to `expm_nsub`.

The three agree exactly when `Z` is constant across regions, because a
scalar multiple of the identity commutes with the generator. They also
agree when `Z` is all zero (movement only) and when movement is absent
(`Move` the identity, `Q` all zero).
