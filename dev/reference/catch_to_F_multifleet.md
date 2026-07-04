# Convert target catches to fishing mortality rates for multiple fleets

Solves for the vector of fleet-specific fishing mortality rates
\\\mathbf{F} = (F_1, \ldots, F_k)\\ that simultaneously satisfy the
Baranov catch equations for all fleets, given a vector of target
catches. Total mortality at age accounts for retained and discard
mortality contributions from all fleets: \\Z_a = M_a + \sum_f (F_f \\
s\_{a,f} \\ r\_{a,f} + F_f \\ s\_{a,f} \\ (1 - r\_{a,f}) \\ d_f)\\,
where \\r\_{a,f}\\ is retention selectivity and \\d_f\\ is
fleet-specific discard mortality. The system of equations is solved via
[`nleqslv`](https://bertcarnell.github.io/nleqslv/reference/nleqslv.html).
Intended for closed-loop MSE harvest control rules with multiple
interacting fishery fleets.

## Usage

``` r
catch_to_F_multifleet(
  target_catch,
  NAA,
  WAA,
  natmort,
  fish_sel,
  ret_sel = array(1, dim = dim(fish_sel)),
  dmr = rep(0, length(target_catch)),
  f_init = 0.05,
  control = list(btol = 1e-06)
)
```

## Arguments

- target_catch:

  Numeric vector `[n_fleets]`. Target catch in biomass units for each
  fleet.

- NAA:

  Numeric matrix `[n_ages x n_sexes]`. Numbers-at-age.

- WAA:

  Numeric matrix `[n_ages x n_sexes]`. Weight-at-age.

- natmort:

  Numeric matrix `[n_ages x n_sexes]`. Instantaneous natural mortality
  rate.

- fish_sel:

  Numeric array `[n_ages x n_sexes x n_fleets]`. Fishery selectivity for
  each fleet.

- ret_sel:

  Numeric array `[n_ages x n_sexes x n_fleets]`. Retention selectivity
  for each fleet, i.e. the proportion of selected catch that is
  retained. Default is an array of 1s matching `dim(fish_sel)` (full
  retention for all fleets).

- dmr:

  Numeric vector `[n_fleets]`. Fleet-specific discard mortality rates,
  i.e. the fraction of discarded fish that die. Default
  `rep(0, length(target_catch))` (all discards survive).

- f_init:

  Numeric scalar or vector `[n_fleets]`. Starting values for the \\F\\
  solver. If a scalar is supplied it is recycled across all fleets.
  Default `0.05`.

- control:

  Named list of control parameters passed to
  [`nleqslv`](https://bertcarnell.github.io/nleqslv/reference/nleqslv.html).
  Default `list(btol = 1e-6)`.

## Value

Numeric vector `[n_fleets]` of solved fishing mortality rates, one per
fleet.

## See also

Other Closed Loop Simulations:
[`catch_to_F_singlefleet()`](https://chengmatt.github.io/SPoRC/dev/reference/catch_to_F_singlefleet.md),
[`condition_closed_loop_simulations()`](https://chengmatt.github.io/SPoRC/dev/reference/condition_closed_loop_simulations.md),
[`get_closed_loop_reference_points()`](https://chengmatt.github.io/SPoRC/dev/reference/get_closed_loop_reference_points.md)
