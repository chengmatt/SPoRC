# Go from TAC to Fishing Mortality using bisection for when a single fishery fleet exists

Go from TAC to Fishing Mortality using bisection for when a single
fishery fleet exists

## Usage

``` r
catch_to_F_singlefleet(
  f_guess,
  catch,
  NAA,
  WAA,
  natmort,
  fish_sel,
  n.iter = 20,
  lb = 0,
  ub = 2
)
```

## Arguments

- f_guess:

  Initial guess of F

- catch:

  Provided catch values

- NAA:

  Numbers, dimensioned by ages, and sexes

- WAA:

  Weight, dimensioned by ages and sexes

- natmort:

  Natural mortality dimensioned by ages and sex

- fish_sel:

  Fishery selectivity, dimesnioned by ages and sex

- n.iter:

  Number of iterations for bisection

- lb:

  Lower bound of F

- ub:

  Upper bound of F

## Value

Fishing mortality values for a single fleet

## See also

Other Closed Loop Simulations:
[`catch_to_F_multifleet()`](https://chengmatt.github.io/SPoRC/reference/catch_to_F_multifleet.md),
[`condition_closed_loop_simulations()`](https://chengmatt.github.io/SPoRC/reference/condition_closed_loop_simulations.md),
[`get_closed_loop_reference_points()`](https://chengmatt.github.io/SPoRC/reference/get_closed_loop_reference_points.md)
