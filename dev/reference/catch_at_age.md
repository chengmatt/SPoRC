# Catch-at-age over one season, consistent with the movement timing

Applies the catch equation to a numbers-at-region vector in the way each
`move_timing` requires. Catch is taken where the fish actually are
during the season: at their post-movement locations under
`move_timing = 0` (movement happens at the start of the season), at
their pre-movement locations under `move_timing = 1` (movement happens
at the end). Under `move_timing = 2` fish redistribute among regions
while they are being caught, so the region-local \\F/Z\\(1 -
e^{-Z})\\N\\ form is invalid and the season-integrated (spatial Baranov)
abundance is used instead.

## Usage

``` r
catch_at_age(
  N,
  Move,
  Z,
  Q = NULL,
  dur = 1,
  F_landed,
  move_timing = 0,
  expm_nsub = 0
)
```

## Arguments

- N:

  Numeric vector of length `n_regions` at the start of the season.

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

- F_landed:

  Numeric vector of length `n_regions` giving the seasonal fishing
  mortality to apply (retained, discarded or landed, per the caller).

- move_timing:

  Integer flag for the movement/mortality ordering: `0` = movement then
  mortality (default, historical SPoRC behavior), `1` = mortality then
  movement, `2` = continuous (simultaneous) movement and mortality.

- expm_nsub:

  Integer controlling how the matrix exponential is evaluated under
  `move_timing = 2`: `0` (default) uses
  [`Matrix::expm`](https://rdrr.io/pkg/Matrix/man/expm-methods.html), a
  value \\n \ge 1\\ uses the implicit backward Euler scheme \\(I -
  A/n)^{-n}\\. See
  [`mat_exp`](https://chengmatt.github.io/SPoRC/dev/reference/mat_exp.md).

## Value

Numeric vector of length `n_regions` giving catch by region.
