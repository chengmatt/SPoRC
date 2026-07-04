# Convert a target catch to fishing mortality for a single fleet via bisection

Uses interval bisection to find the scalar fishing mortality \\F\\ that
produces a predicted catch (Baranov catch equation, summed over ages and
sexes in biomass) equal to `catch`. Retained and discarded catch
components are handled separately via retention selectivity and discard
mortality. Intended for closed-loop MSE harvest control rules where a
TAC must be translated into an \\F\\ for the operating model.

## Usage

``` r
catch_to_F_singlefleet(
  f_guess,
  catch,
  NAA,
  WAA,
  natmort,
  fish_sel,
  ret_sel = {
     tmp = fish_sel[]
tmp[] <- 1
     tmp
 },
  dmr = 0,
  n.iter = 20,
  lb = 0,
  ub = 2
)
```

## Arguments

- f_guess:

  Numeric. Initial \\F\\ guess (not used directly by the bisection
  algorithm but retained for API consistency).

- catch:

  Numeric. Target catch in biomass units.

- NAA:

  Numeric matrix `[n_ages x n_sexes]`. Numbers-at-age at the start of
  the time step.

- WAA:

  Numeric matrix `[n_ages x n_sexes]`. Weight-at-age.

- natmort:

  Numeric matrix `[n_ages x n_sexes]`. Instantaneous natural mortality
  rate.

- fish_sel:

  Numeric matrix `[n_ages x n_sexes]`. Fishery selectivity (scaled to a
  maximum of 1).

- ret_sel:

  Numeric matrix `[n_ages x n_sexes]`. Retention selectivity, i.e. the
  proportion of selected catch that is retained (vs. discarded). Default
  is an array of 1s matching `dim(fish_sel)` (full retention).

- dmr:

  Numeric scalar. Discard mortality rate, i.e. the fraction of discarded
  fish that die. Default `0` (all discards survive).

- n.iter:

  Integer. Number of bisection iterations. Default `20`; approximately
  \\\log_2((ub - lb) / \epsilon)\\ iterations are required for tolerance
  \\\epsilon\\.

- lb:

  Numeric. Lower bound of the \\F\\ search interval. Default `0`.

- ub:

  Numeric. Upper bound of the \\F\\ search interval. Default `2`.

## Value

Scalar numeric. The \\F\\ value at the final bisection midpoint that
most closely produces `catch`.

## See also

Other Closed Loop Simulations:
[`catch_to_F_multifleet()`](https://chengmatt.github.io/SPoRC/dev/reference/catch_to_F_multifleet.md),
[`condition_closed_loop_simulations()`](https://chengmatt.github.io/SPoRC/dev/reference/condition_closed_loop_simulations.md),
[`get_closed_loop_reference_points()`](https://chengmatt.github.io/SPoRC/dev/reference/get_closed_loop_reference_points.md)
