# Set up a conditional age-at-length data source

Shared by
[`Setup_Mod_FishIdx_and_Comps`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_FishIdx_and_Comps.md)
and
[`Setup_Mod_SrvIdx_and_Comps`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_SrvIdx_and_Comps.md).
Validates the CAAL arrays, converts the likelihood and composition type
specifications, fills in a default input sample size when none is
supplied, and creates the Dirichlet-multinomial overdispersion
parameters together with their maps. The likelihood weights
(`Wt_Fish_caal`, `Wt_Srv_caal`) belong to
[`Setup_Mod_Weighting`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Weighting.md)
with every other weight.

## Usage

``` r
setup_caal_source(
  input_list,
  ObsCAAL,
  UseCAAL,
  ISS_CAAL,
  CAAL_LikeType,
  CAAL_Type,
  fleet_type
)
```

## Arguments

- input_list:

  Named list with `$data`, `$par` and `$map`.

- ObsCAAL:

  Observed CAAL array
  `[n_regions x n_years x n_seas x n_lens x n_ages x n_sexes x n_fleets]`,
  or `NULL` for none.

- UseCAAL:

  Use flags `[n_regions x n_years x n_seas x n_lens x n_fleets]`. A
  length bin with no aged fish has a zero and is skipped.

- ISS_CAAL:

  Input sample sizes
  `[n_regions x n_years x n_seas x n_lens x n_sexes x n_fleets]`, the
  number aged within each length bin. When `NULL` it is summed from
  `ObsCAAL`.

- CAAL_LikeType:

  Character vector of length `n_fleets`. One of `"none"`,
  `"Multinomial"` or `"Dirichlet-Multinomial"`.

- CAAL_Type:

  Character vector of composition type specifications, using the same
  `"CompType_Year_x-y_Fleet_z"` convention as the marginal compositions.

- fleet_type:

  Character, either `"Fish"` or `"Srv"`.

## Value

The updated `input_list`.

## Details

Only the multinomial and Dirichlet-multinomial families are available. A
CAAL row is the age composition of the otoliths taken from one length
bin, so it is usually a small sample that is mostly zeros, and the
logistic-normal forms need an additive log-ratio transform that such a
row cannot support.

Supplying any CAAL data switches on `do_caal`, since the likelihood
reads the joint arrays at length and age that flag builds.
