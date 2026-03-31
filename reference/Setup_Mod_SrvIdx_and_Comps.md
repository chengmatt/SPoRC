# Setup observed survey indices and composition data (age and length comps)

Setup observed survey indices and composition data (age and length
comps)

## Usage

``` r
Setup_Mod_SrvIdx_and_Comps(
  input_list,
  ObsSrvIdx,
  ObsSrvIdx_SE,
  UseSrvIdx,
  srv_idx_type,
  ObsSrvAgeComps,
  UseSrvAgeComps,
  ObsSrvLenComps,
  UseSrvLenComps,
  ISS_SrvAgeComps,
  ISS_SrvLenComps,
  SrvAgeComps_LikeType,
  SrvLenComps_LikeType,
  SrvAgeComps_Type,
  SrvLenComps_Type,
  ...
)
```

## Arguments

- input_list:

  List containing a data list, parameter list, and map list

- ObsSrvIdx:

  Observed survey index data as a numeric array with dimensions
  `[n_regions, n_years, n_srv_fleets]`.

- ObsSrvIdx_SE:

  Standard errors associated with `ObsSrvIdx`, also dimensioned
  `[n_regions, n_years, n_srv_fleets]`.

- UseSrvIdx:

  Logical or binary indicator array
  (`[n_regions, n_years, n_srv_fleets]`) specifying whether to include a
  survey index in the likelihood (`1`) or ignore it (`0`).

- srv_idx_type:

  Character vector of length `n_srv_fleets` specifying the type of index
  data. Options are `"abd"` for abundance, `"biom"` for biomass, and
  `"none"` if no index is available.

- ObsSrvAgeComps:

  Observed survey age composition data as a numeric array with
  dimensions `[n_regions, n_years, n_ages, n_sexes, n_srv_fleets]`.
  Values should reflect counts or proportions (not required to sum to 1,
  but should be on a comparable scale).

- UseSrvAgeComps:

  Indicator array (`[n_regions, n_years, n_srv_fleets]`) specifying
  whether to fit survey age composition data (`1`) or ignore it (`0`).

- ObsSrvLenComps:

  Observed survey length composition data as a numeric array with
  dimensions `[n_regions, n_years, n_lens, n_sexes, n_srv_fleets]`.
  Values should reflect counts or proportions.

- UseSrvLenComps:

  Indicator array (`[n_regions, n_years, n_srv_fleets]`) specifying
  whether to fit survey length composition data (`1`) or ignore it
  (`0`).

- ISS_SrvAgeComps:

  Input sample size for age compositions, array dimensioned
  `[n_regions, n_years, n_sexes, n_srv_fleets]`. Required if observed
  age comps are normalized (i.e., sum to 1), to correctly scale the
  contribution to the likelihood.

- ISS_SrvLenComps:

  Same as `ISS_SrvAgeComps`, but for length compositions.

- SrvAgeComps_LikeType:

  Character vector of length `n_srv_fleets` specifying the likelihood
  type used for survey age composition data. Options include
  `"Multinomial"`, `"Dirichlet-Multinomial"`, and
  `"iid-Logistic-Normal"`. Use `"none"` to omit the likelihood.

- SrvLenComps_LikeType:

  Same as `SrvAgeComps_LikeType`, but for survey length composition
  data.

- SrvAgeComps_Type:

  Character vector specifying how age compositions are structured by
  fleet and year range. Options include:

  - `"agg"`: Aggregated across regions and sexes.

  - `"spltRspltS"`: Split by region and by sex (compositions sum to 1
    within region-sex group).

  - `"spltRjntS"`: Split by region but summed jointly across sexes.

  - `"none"`: No composition data used.

  Format each element as
  `"<type>_Year_<start>-<end>_Fleet_<fleet number>"` (e.g.,
  `"agg_Year_1-10_Fleet_1"`).

- SrvLenComps_Type:

  Same as `SrvAgeComps_Type`, but for length compositions.

- ...:

  Additional arguments specifying starting values for overdispersion
  parameters (e.g., `ln_SrvAge_theta`, `ln_SrvLen_theta`,
  `ln_SrvAge_theta_agg`, `ln_SrvLen_theta_agg`).

## See also

Other Model Setup:
[`Setup_Mod_Biologicals()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Biologicals.md),
[`Setup_Mod_Catch_and_F()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Catch_and_F.md),
[`Setup_Mod_Dim()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Dim.md),
[`Setup_Mod_FishIdx_and_Comps()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_FishIdx_and_Comps.md),
[`Setup_Mod_Fishsel_and_Q()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Fishsel_and_Q.md),
[`Setup_Mod_Movement()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Movement.md),
[`Setup_Mod_Rec()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Rec.md),
[`Setup_Mod_Srvsel_and_Q()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Srvsel_and_Q.md),
[`Setup_Mod_Tagging()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Tagging.md),
[`Setup_Mod_Weighting()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Weighting.md)
