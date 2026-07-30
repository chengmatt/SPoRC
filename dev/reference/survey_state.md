# Numbers at the survey point within a season

Propagates a numbers-at-region vector from the start of a season to the
survey timing `t_srv` within it, consistently with `move_timing`. A
survey index is a snapshot at a point inside the season, not an
accumulation over it, so this uses the same partial propagation as
[`spawn_state`](https://chengmatt.github.io/SPoRC/dev/reference/spawn_state.md)
rather than the season integral used for catch.

## Usage

``` r
survey_state(N, Move, Z, Q = NULL, dur = 1, t_srv = 0, move_timing = 0)
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

- t_srv:

  Numeric vector of length `n_regions` giving the fraction of the season
  elapsed before the survey, per region. A scalar is recycled.

- move_timing:

  Integer flag for the movement/mortality ordering: `0` = movement then
  mortality (default, historical SPoRC behaviour), `1` = mortality then
  movement, `2` = continuous (simultaneous) movement and mortality.

## Value

Numeric vector of length `n_regions`.

## Details

Under `move_timing` 0 and 1 the fish occupy one region for the whole
season (post-movement and pre-movement respectively), so the historical
\\N \exp(-t\_{srv} Z)\\ form is exact and is used unchanged.

Under `move_timing = 2` the population is propagated by \\\exp(A \\
t\_{srv})\\ with \\A = Q^\top \Delta - \mathrm{diag}(Z)\\.

`t_srv` may differ by region, which a single propagation operator cannot
represent: fish observed in region \\r\\ arrived from regions whose
elapsed times differ. The convention adopted here is that the survey in
region \\r\\ observes the population propagated to *that region's*
survey time, i.e. element \\r\\ of \\\exp(A \\ t\_{srv,r}) N\\. When
`t_srv` is constant across regions (the usual case) this reduces to a
single propagation and one matrix exponential.
