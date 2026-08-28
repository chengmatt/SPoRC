# Advance a numbers-at-region vector across one season

Thin convenience wrapper around `build_seas_operator` that applies the
seasonal transition to a single numbers-at-region vector. For
`move_timing` 0 and 1 the operator is never formed explicitly, which
keeps the AD tape smaller than the equivalent matrix product.

## Usage

``` r
advance_seas(N, Move, Z, Q = NULL, dur = 1, move_timing = 0, expm_nsub = 0)
```

## Arguments

- N:

  Numeric vector of length `n_regions`.

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

Numeric vector of length `n_regions`.
