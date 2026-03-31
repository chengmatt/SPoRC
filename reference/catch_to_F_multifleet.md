# Solve for fishing mortality rates that achieve target catches for multiple fleets

Solve for fishing mortality rates that achieve target catches for
multiple fleets

## Usage

``` r
catch_to_F_multifleet(
  target_catch,
  NAA,
  WAA,
  natmort,
  fish_sel,
  f_init = 0.05,
  control = list(btol = 1e-06)
)
```

## Arguments

- target_catch:

  Numeric vector of target catch values for each fleet

- NAA:

  Matrix of numbers-at-age (ages x sexes)

- WAA:

  Matrix of weight-at-age (ages x sexes)

- natmort:

  Matrix of natural mortality (ages x sexes)

- fish_sel:

  3D array of fishery selectivity (ages x sexes x fleets)

- f_init:

  Initial guess for F values (scalar or vector)

- control:

  List of control parameters for nleqslv

## Value

Numeric vector of F values for each fleet

## See also

Other Closed Loop Simulations:
[`catch_to_F_singlefleet()`](https://chengmatt.github.io/SPoRC/reference/catch_to_F_singlefleet.md),
[`condition_closed_loop_simulations()`](https://chengmatt.github.io/SPoRC/reference/condition_closed_loop_simulations.md),
[`get_closed_loop_reference_points()`](https://chengmatt.github.io/SPoRC/reference/get_closed_loop_reference_points.md)
