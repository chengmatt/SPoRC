# Map fishery composition correlation parameters

Constructs factor maps for a composition correlation-parameter array
(e.g. `FishAge_corr_pars`, region- and sex-specific AR1/sex correlation)
and its aggregated counterpart (e.g. `FishAge_corr_pars_agg`) for 1D and
2D logistic-normal composition likelihoods. Parameters are activated
only when the corresponding `LikeType` is in `c(3, 4)` (1D / 2D
logistic-normal); all other likelihoods, or fleets with no observed
compositions, map correlation parameters to `NA`. For the 2D
logistic-normal (`LikeType == 4`), both trailing elements of the
`[...,2]` slice are activated: element 1 for the AR1 coefficient and
element 2 for the sex correlation (skipped when `n_sexes == 1`).

## Usage

``` r
do_comp_corr_pars_mapping(
  input_list,
  comp_prefix,
  discard = FALSE,
  has_pop = FALSE,
  fleet_field = "n_fish_fleets"
)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map` sublists.

- comp_prefix:

  Character, either `"FishAge"` or `"FishLen"`.

- discard:

  Logical. If `TRUE`, maps the discard variant. Default `FALSE`.

- has_pop:

  Logical. If `TRUE`, maps the population-specific variant. Default
  `FALSE`.

- fleet_field:

  Character. Name of the `$data` field giving the number of fleets to
  loop over. Default `"n_fish_fleets"`; pass `"n_srv_fleets"` for survey
  composition types (`"SrvAge"`, `"SrvLen"`), which have no `discard`
  variant.

## Value

The input `input_list` with the corresponding `$map$<...>_corr_pars` and
`$map$<...>_corr_pars_agg` set to factor vectors.

## Details

One function serves every composition block. `comp_prefix` selects age
versus length and fishery versus survey, `discard` selects the retained
or discarded data source, `has_pop` selects the aggregated or population
specific data source, and `fleet_field` names the fleet count to size
the map by. Same parameterization as
[`do_comp_theta_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_comp_theta_mapping.md).
