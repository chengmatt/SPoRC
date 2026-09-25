# Set up discard age and length composition inputs

Sets the observed discard compositions, pooled and population-specific,
with their use flags, input sample sizes, likelihood and composition
types, and the overdispersion and correlation parameters and maps.

## Usage

``` r
Setup_Mod_Discard_Comps(
  input_list,
  ObsFishAgeComps_discard,
  UseFishAgeComps_discard,
  ISS_FishAgeComps_discard,
  ObsFishLenComps_discard,
  UseFishLenComps_discard,
  ISS_FishLenComps_discard,
  FishAgeComps_discard_LikeType,
  FishLenComps_discard_LikeType,
  FishAgeComps_discard_Type,
  FishLenComps_discard_Type,
  ObsFishAgeComps_discard_pop,
  UseFishAgeComps_discard_pop,
  ISS_FishAgeComps_discard_pop,
  ObsFishLenComps_discard_pop,
  UseFishLenComps_discard_pop,
  ISS_FishLenComps_discard_pop,
  FishAgeComps_discard_pop_LikeType,
  FishLenComps_discard_pop_LikeType,
  FishAgeComps_discard_pop_Type,
  FishLenComps_discard_pop_Type,
  FishAgeComps_discard_bins = NULL,
  FishLenComps_discard_bins = NULL,
  FishAgeComps_discard_pop_bins = NULL,
  FishLenComps_discard_pop_bins = NULL,
  ...
)
```

## Arguments

- input_list:

  Main model input list with data, parameter and mapping lists.

- ObsFishAgeComps_discard:

  Observed discard age compositions
  `[n_regions, n_years, n_seas, n_sexes, n_fish_fleets]`.

- UseFishAgeComps_discard:

  Use flags `[n_regions, n_years, n_seas, n_fish_fleets]`.

- ISS_FishAgeComps_discard:

  Optional input sample sizes
  `[n_regions, n_years, n_seas, n_sexes, n_fish_fleets]`.

- ObsFishLenComps_discard:

  Observed discard length compositions
  `[n_regions, n_years, n_seas, n_lens, n_sexes, n_fish_fleets]`.

- UseFishLenComps_discard:

  Use flags `[n_regions, n_years, n_seas, n_fish_fleets]`.

- ISS_FishLenComps_discard:

  Optional input sample sizes
  `[n_regions, n_years, n_seas, n_sexes, n_fish_fleets]`.

- FishAgeComps_discard_LikeType, FishLenComps_discard_LikeType:

  Character vectors `[n_fish_fleets]`, one of `"none"`, `"Multinomial"`,
  `"Dirichlet-Multinomial"`, the three logistic-normal forms or their
  three `-miss0` counterparts.

- FishAgeComps_discard_Type, FishLenComps_discard_Type:

  Encoded composition structure by year and fleet.

- ObsFishAgeComps_discard_pop:

  Observed population-specific discard age compositions
  `[n_pop, n_regions, n_years, n_seas, n_sexes, n_fish_fleets]`.

- UseFishAgeComps_discard_pop:

  Use flags `[n_pop, n_regions, n_years, n_seas, n_fish_fleets]`.

- ISS_FishAgeComps_discard_pop:

  Optional input sample sizes
  `[n_pop, n_regions, n_years, n_seas, n_sexes, n_fish_fleets]`.

- ObsFishLenComps_discard_pop:

  Observed population-specific discard length compositions
  `[n_pop, n_regions, n_years, n_seas, n_lens, n_sexes, n_fish_fleets]`.

- UseFishLenComps_discard_pop:

  Use flags `[n_pop, n_regions, n_years, n_seas, n_fish_fleets]`.

- ISS_FishLenComps_discard_pop:

  Optional input sample sizes
  `[n_pop, n_regions, n_years, n_seas, n_sexes, n_fish_fleets]`.

- FishAgeComps_discard_pop_LikeType, FishLenComps_discard_pop_LikeType:

  Character vectors `[n_fish_fleets]` for the population-specific
  discard compositions.

- FishAgeComps_discard_pop_Type, FishLenComps_discard_pop_Type:

  Encoded composition structure for the population-specific discard
  compositions.

- FishAgeComps_discard_bins:

  Which age bins each fleet's discard age composition is fitted over,
  either a list with one element per fleet (bin indices, or `NULL` for
  all) or an `[n_obs_ages x n_fish_fleets]` array of 0/1 weights.
  Indices are observed bins, after any ageing error. Excluded bins leave
  the likelihood rather than being forced to be explained. Default
  `NULL`, all bins.

- FishLenComps_discard_bins, FishAgeComps_discard_pop_bins,
  FishLenComps_discard_pop_bins:

  Which bins the discard length, population-specific discard age and
  population-specific discard length compositions are fitted over, in
  the same format as `FishAgeComps_discard_bins`.

- ...:

  Optional starting values for the overdispersion and correlation
  parameters.

## Value

`input_list` with the discard composition arrays, the computed or
supplied input sample sizes, the integer-coded likelihood types, the
composition type matrices by year and fleet, the population-specific
discard structures, and the overdispersion and correlation parameters
with their maps.

## See also

Other Model Setup:
[`Setup_Mod_Biologicals()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Biologicals.md),
[`Setup_Mod_Catch_and_F()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Catch_and_F.md),
[`Setup_Mod_Dim()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Dim.md),
[`Setup_Mod_FishIdx_and_Comps()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_FishIdx_and_Comps.md),
[`Setup_Mod_Fishsel_and_Q()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Fishsel_and_Q.md),
[`Setup_Mod_Movement()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Movement.md),
[`Setup_Mod_Rec()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Rec.md),
[`Setup_Mod_Retsel()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Retsel.md),
[`Setup_Mod_SrvIdx_and_Comps()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_SrvIdx_and_Comps.md),
[`Setup_Mod_Srvsel_and_Q()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Srvsel_and_Q.md),
[`Setup_Mod_Tagging()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Tagging.md),
[`Setup_Mod_Weighting()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Weighting.md)
