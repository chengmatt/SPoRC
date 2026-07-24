# Map fishery composition overdispersion (theta) parameters

Constructs factor maps for a composition overdispersion parameter (e.g.
`ln_FishAge_theta`) and its aggregated counterpart (e.g.
`ln_FishAge_theta_agg`) based on the composition type and likelihood
specified in the corresponding `$data` fields. Parameters are mapped to
`NA` for fleets using multinomial likelihoods (`LikeType == 0`), with no
observed compositions, or (per-region/per-pop-region) with no active
comps for that cell.

## Usage

``` r
do_comp_theta_mapping(
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

  Logical. If `TRUE`, maps the discard variant (`..._discard_...`
  fields). Default `FALSE`.

- has_pop:

  Logical. If `TRUE`, maps the population-specific variant
  (`..._pop_...` fields, with an added population dimension and loop).
  Default `FALSE`.

- fleet_field:

  Character. Name of the `$data` field giving the number of fleets to
  loop over. Default `"n_fish_fleets"`; pass `"n_srv_fleets"` for survey
  composition types (`"SrvAge"`, `"SrvLen"`), which have no `discard`
  variant.

## Value

The input `input_list` with the corresponding `$map$ln_<...>_theta` and
`$map$ln_<...>_theta_agg` set to factor vectors. Active parameters
receive sequential integer indices; inactive parameters are `NA`.
