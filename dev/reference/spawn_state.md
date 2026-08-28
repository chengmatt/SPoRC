# Numbers at the spawning point within a season

Propagates a numbers-at-region vector from the start of a season to the
spawning point `t_spawn` within it, consistently with `move_timing`.
Each timing has its own natural spawning state, so no additional
convention is imposed:

- `move_timing = 0`:

  Fish move at the start of the season and then experience `t_spawn`
  worth of mortality, so spawning happens at the post-movement location.
  This reproduces the historical SPoRC calculation.

- `move_timing = 1`:

  Movement occurs at the end of the season, so spawning happens at the
  pre-movement location after `t_spawn` worth of mortality.

- `move_timing = 2`:

  The population is propagated a fraction `t_spawn` of the way through
  the season under the combined movement-mortality generator, so
  spawners are partially redistributed.

## Usage

``` r
spawn_state(
  N,
  Move,
  Z,
  Q = NULL,
  dur = 1,
  t_spawn = 0,
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

- t_spawn:

  Fraction of the season elapsed before spawning, in `[0, 1]`.

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
