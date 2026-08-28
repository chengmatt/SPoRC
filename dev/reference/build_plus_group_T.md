# Build annual transition matrices for the plus-group analytical solution

Constructs the four annual transition matrices needed to solve for the
equilibrium plus-group abundance analytically. Each matrix accumulates
survival and movement across all seasons for either the penultimate age
or the plus-group age, under either unfished or fished conditions.

## Usage

``` r
build_plus_group_T(
  M_penult,
  M_plus,
  F_penult,
  F_plus,
  Mov_penult,
  Mov_plus,
  n_regions,
  n_seas,
  seasdur,
  Mrate_penult = NULL,
  Mrate_plus = NULL,
  move_timing = 0,
  expm_nsub = 0
)
```

## Arguments

- M_penult:

  Numeric vector `[n_regions]`. Natural mortality for the penultimate
  age class, used as an annual rate (scaled by `seasdur` internally).

- M_plus:

  Numeric vector `[n_regions]`. Natural mortality for the plus-group age
  class.

- F_penult:

  Numeric matrix `[n_regions, n_seas]`. Total fishing mortality per
  season for the penultimate age, already summed across fleets.

- F_plus:

  Numeric matrix `[n_regions, n_seas]`. Total fishing mortality per
  season for the plus-group age.

- Mov_penult:

  Numeric array `[n_regions, n_regions, n_seas]`. Movement transition
  matrices for the penultimate age. Entry `[r1, r2, s]` is the
  probability of moving from region `r1` to region `r2` in season `s`.

- Mov_plus:

  Numeric array `[n_regions, n_regions, n_seas]`. Movement transition
  matrices for the plus-group age.

- n_regions:

  Integer. Number of spatial regions.

- n_seas:

  Integer. Number of seasons.

- seasdur:

  Numeric vector `[n_seas]`. Fractional duration of each season (must
  sum to one).

- Mrate_penult:

  Numeric array `[n_regions, n_regions, n_seas]`. Instantaneous movement
  rate matrices (generator \\Q\\) for the penultimate age, in the same
  `[from, to]` convention as `Mov_penult`. Required when
  `move_timing = 2`, ignored otherwise.

- Mrate_plus:

  Numeric array `[n_regions, n_regions, n_seas]`. Instantaneous movement
  rate matrices for the plus-group age. Required when `move_timing = 2`,
  ignored otherwise.

- move_timing:

  Integer flag for movement/mortality sequencing: `0` = movement then
  mortality (default), `1` = mortality then movement, `2` = continuous.
  See
  [`build_seas_operator`](https://chengmatt.github.io/SPoRC/dev/reference/build_seas_operator.md).

- expm_nsub:

  Integer controlling how the matrix exponential is evaluated under
  `move_timing = 2`: `0` uses
  [`Matrix::expm`](https://rdrr.io/pkg/Matrix/man/expm-methods.html),
  \\n \ge 1\\ uses \\n\\ implicit backward Euler substeps. See
  [`mat_exp`](https://chengmatt.github.io/SPoRC/dev/reference/mat_exp.md).

## Value

A named list with four transition matrices, each of dimension
`[n_regions, n_regions]`:

- `T_penult_unfished`:

  Annual transition for the penultimate age under unfished conditions.

- `T_plus_unfished`:

  Annual transition for the plus-group age under unfished conditions.

- `T_penult_fished`:

  Annual transition for the penultimate age under fished conditions.

- `T_plus_fished`:

  Annual transition for the plus-group age under fished conditions.

## Details

The equilibrium plus-group vector \\N\_+\\ satisfies \$\$N\_+ = T\_+
N\_+ + T\_{n-1} N\_{n-1}\$\$ which rearranges to \\(I - T\_+) N\_+ =
T\_{n-1} N\_{n-1}\\, solved in
[`solve_plus_group`](https://chengmatt.github.io/SPoRC/dev/reference/solve_plus_group.md).

All arguments are sliced by the caller to remove the population
dimension, so this helper works identically for the single-population
spatial case (`global_Fmsy`, `local_Fmsy_sglpop`) and the
multi-population case (`global_SPR`, `local_Fmsy_multipop`).
