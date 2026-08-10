# Resolve the factor scale and loading for one simulated index cell

Looks a cell up in a fleet's factor decomposition from
[`build_idx_factor`](https://chengmatt.github.io/SPoRC/dev/reference/build_idx_factor.md).
A cell inside the covariance gets its own marginal scale and loading; a
cell outside it (a projection year, or a cell the fleet never fit) gets
the mean scale and loading, so closed loop can keep drawing the series
past the end of the data.

## Usage

``` r
resolve_idx_factor(mvn, r, y, seas)
```

## Arguments

- mvn:

  One fleet's element from
  [`build_idx_factor`](https://chengmatt.github.io/SPoRC/dev/reference/build_idx_factor.md).

- r, y, seas:

  Integer indices of the simulated cell.

## Value

List with scalars `d` and `lambda`.
