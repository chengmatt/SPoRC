# Solve for equilibrium plus-group numbers given transition matrices

Given the annual transition matrices produced by
[`build_plus_group_T`](https://chengmatt.github.io/SPoRC/dev/reference/build_plus_group_T.md)
and the penultimate-age abundance vector, solves the linear system
\\(I - T\_+) N\_+ = T\_{n-1} N\_{n-1}\\ for the equilibrium plus-group
abundance under both unfished and fished conditions.

## Usage

``` r
solve_plus_group(Ts, N_penult_u, N_penult_f, n_regions)
```

## Arguments

- Ts:

  Named list returned by
  [`build_plus_group_T`](https://chengmatt.github.io/SPoRC/dev/reference/build_plus_group_T.md).

- N_penult_u:

  Numeric vector `[n_regions]`. Unfished penultimate-age abundance
  (per-recruit) at the start of the year.

- N_penult_f:

  Numeric vector `[n_regions]`. Fished penultimate-age abundance
  (per-recruit) at the start of the year.

- n_regions:

  Integer. Number of spatial regions.

## Value

A named list:

- `unfished`:

  Numeric vector `[n_regions]`. Equilibrium plus-group abundance per
  recruit under unfished conditions.

- `fished`:

  Numeric vector `[n_regions]`. Equilibrium plus-group abundance per
  recruit under fished conditions.
