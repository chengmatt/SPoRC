# Setup observed fishery indices and composition data (age and length comps)

Setup observed fishery indices and composition data (age and length
comps)

## Usage

``` r
Setup_Mod_FishIdx_and_Comps(
  input_list,
  ObsFishIdx,
  ObsFishIdx_SE,
  fish_idx_type,
  UseFishIdx,
  ObsFishAgeComps,
  UseFishAgeComps,
  ISS_FishAgeComps,
  ObsFishLenComps,
  UseFishLenComps,
  ISS_FishLenComps,
  FishAgeComps_LikeType,
  FishLenComps_LikeType,
  FishAgeComps_Type,
  FishLenComps_Type,
  ...
)
```

## Arguments

- input_list:

  List containing a data list, parameter list, and map list

- ObsFishIdx:

  Observed fishery index data as a numeric array with dimensions
  `[n_regions, n_years, n_fish_fleets]`.

- ObsFishIdx_SE:

  Standard errors associated with `ObsFishIdx`, also dimensioned
  `[n_regions, n_years, n_fish_fleets]`.

- fish_idx_type:

  Character vector of length `n_fish_fleets` specifying the type of
  index data. Options are `"abd"` for abundance, `"biom"` for biomass,
  and `"none"` if no index is available.

- UseFishIdx:

  Logical or binary indicator array
  (`[n_regions, n_years, n_fish_fleets]`) specifying whether to include
  a fishery index in the likelihood (`1`) or ignore it (`0`).

- ObsFishAgeComps:

  Observed fishery age composition data as a numeric array with
  dimensions `[n_regions, n_years, n_ages, n_sexes, n_fish_fleets]`.
  Values should reflect counts or proportions (not required to sum to 1,
  but should be on a comparable scale).

- UseFishAgeComps:

  Indicator array (`[n_regions, n_years, n_fish_fleets]`) specifying
  whether to fit fishery age composition data (`1`) or ignore it (`0`).

- ISS_FishAgeComps:

  Input sample size for age compositions, array dimensioned
  `[n_regions, n_years, n_sexes, n_fish_fleets]`. Required if observed
  age comps are normalized (i.e., sum to 1), to correctly scale the
  contribution to the likelihood.

- ObsFishLenComps:

  Observed fishery length composition data as a numeric array with
  dimensions `[n_regions, n_years, n_lens, n_sexes, n_fish_fleets]`.
  Values should reflect counts or proportions.

- UseFishLenComps:

  Indicator array (`[n_regions, n_years, n_fish_fleets]`) specifying
  whether to fit fishery length composition data (`1`) or ignore it
  (`0`).

- ISS_FishLenComps:

  Same as `ISS_FishAgeComps`, but for length compositions.

- FishAgeComps_LikeType:

  Character vector of length `n_fish_fleets` specifying the likelihood
  type used for fishery age composition data. Options include
  `"Multinomial"`, `"Dirichlet-Multinomial"`, and
  `"iid-Logistic-Normal"`. Use `"none"` to omit the likelihood.

- FishLenComps_LikeType:

  Same as `FishAgeComps_LikeType`, but for fishery length composition
  data.

- FishAgeComps_Type:

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

- FishLenComps_Type:

  Same as `FishAgeComps_Type`, but for length compositions.

- ...:

  Additional arguments specifying starting values for overdispersion
  parameters (e.g., `ln_FishAge_theta`, `ln_FishLen_theta`,
  `ln_FishAge_theta_agg`, `ln_FishLen_theta_agg`).

## See also

Other Model Setup:
[`Setup_Mod_Biologicals()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Biologicals.md),
[`Setup_Mod_Catch_and_F()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Catch_and_F.md),
[`Setup_Mod_Dim()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Dim.md),
[`Setup_Mod_Fishsel_and_Q()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Fishsel_and_Q.md),
[`Setup_Mod_Movement()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Movement.md),
[`Setup_Mod_Rec()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Rec.md),
[`Setup_Mod_SrvIdx_and_Comps()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_SrvIdx_and_Comps.md),
[`Setup_Mod_Srvsel_and_Q()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Srvsel_and_Q.md),
[`Setup_Mod_Tagging()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Tagging.md),
[`Setup_Mod_Weighting()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Weighting.md)
