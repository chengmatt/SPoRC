# Set up discard age and length composition inputs

Set up discard age and length composition inputs

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
  ...
)
```

## Arguments

- input_list:

  Main model input list containing data, parameters, and mapping
  structures.

- ObsFishAgeComps_discard:

  5D array of observed discard age compositions:
  `[n_regions, n_years, n_seas, n_sexes, n_fish_fleets]`

- UseFishAgeComps_discard:

  Matrix indicating whether discard age comps are used:
  `[n_regions, n_years, n_seas, n_fish_fleets]`

- ISS_FishAgeComps_discard:

  Optional ISS (effective sample size) array for discard age comps:
  `[n_regions, n_years, n_seas, n_sexes, n_fish_fleets]`

- ObsFishLenComps_discard:

  6D array of observed discard length compositions:
  `[n_regions, n_years, n_seas, n_lens, n_sexes, n_fish_fleets]`

- UseFishLenComps_discard:

  Matrix indicating whether discard length comps are used:
  `[n_regions, n_years, n_seas, n_fish_fleets]`

- ISS_FishLenComps_discard:

  Optional ISS array for discard length comps:
  `[n_regions, n_years, n_seas, n_sexes, n_fish_fleets]`

- FishAgeComps_discard_LikeType:

  Character vector (length n_fish_fleets) specifying likelihood type:
  one of
  `c("none","Multinomial","Dirichlet-Multinomial","iid-Logistic-Normal","1d-Logistic-Normal","2d-Logistic-Normal")`

- FishLenComps_discard_LikeType:

  Character vector (length n_fish_fleets) specifying likelihood type

- FishAgeComps_discard_Type:

  List/encoded character strings defining composition structure by year
  and fleet.

- FishLenComps_discard_Type:

  List/encoded character strings defining composition structure by year
  and fleet.

- ObsFishAgeComps_discard_pop:

  6D array of population-specific discard age comps:
  `[n_pop, n_regions, n_years, n_seas, n_sexes, n_fish_fleets]`

- UseFishAgeComps_discard_pop:

  5D array indicating use of population age comps:
  `[n_pop, n_regions, n_years, n_seas, n_fish_fleets]`

- ISS_FishAgeComps_discard_pop:

  Optional ISS array for population discard age comps:
  `[n_pop, n_regions, n_years, n_seas, n_sexes, n_fish_fleets]`

- ObsFishLenComps_discard_pop:

  7D array of population-specific discard length comps:
  `[n_pop, n_regions, n_years, n_seas, n_lens, n_sexes, n_fish_fleets]`

- UseFishLenComps_discard_pop:

  5D array indicating use of population length comps:
  `[n_pop, n_regions, n_years, n_seas, n_fish_fleets]`

- ISS_FishLenComps_discard_pop:

  Optional ISS array for population length comps:
  `[n_pop, n_regions, n_years, n_seas, n_sexes, n_fish_fleets]`

- FishAgeComps_discard_pop_LikeType:

  Character vector (length n_fish_fleets) for population age likelihood
  types.

- FishLenComps_discard_pop_LikeType:

  Character vector (length n_fish_fleets) for population length
  likelihood types.

- FishAgeComps_discard_pop_Type:

  Encoded structure definitions for population age comps by year/fleet.

- FishLenComps_discard_pop_Type:

  Encoded structure definitions for population length comps by
  year/fleet.

- ...:

  Optional starting values for parameters (e.g., dispersion,
  correlation)

## Value

The input `input_list` updated with:

- discard age and length composition data structures

- ISS (effective sample size) arrays (computed or supplied)

- likelihood type mappings (integer-coded)

- composition type matrices by year and fleet

- population-specific discard composition structures

- parameter arrays for dispersion and correlation

- mapping configurations for estimation

## Details

All composition arrays follow consistent indexing conventions:

- Age compositions: `[region, year, season, sex, fleet]`

- Length compositions: `[region, year, season, length, sex, fleet]`

- Population age comps: `[pop, region, year, season, sex, fleet]`

- Population length comps:
  `[pop, region, year, season, length, sex, fleet]`

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
