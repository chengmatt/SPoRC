# Extract simulation outputs into SPoRC estimation model format

Subsets and reshapes biological, tagging, fishery, and survey arrays
from a simulation environment or output list to cover years `1:y` and
simulation replicate `sim`, producing a named list ready for direct use
in
[`Setup_Mod_Biologicals`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Biologicals.md),
[`Setup_Mod_Catch_and_F`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Catch_and_F.md),
[`Setup_Mod_SrvIdx_and_Comps`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_SrvIdx_and_Comps.md),
and
[`Setup_Mod_Tagging`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Tagging.md).
Binary `Use*` indicator arrays are derived automatically from the
extracted observation arrays (non-NA, positive values set to 1).

## Usage

``` r
simulation_data_to_SPoRC(sim_env, y, sim)
```

## Arguments

- sim_env:

  Simulation environment or list (e.g., output from
  [`Simulate_Pop_Static`](https://chengmatt.github.io/SPoRC/dev/reference/Simulate_Pop_Static.md)
  or a
  [`Setup_sim_env`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_sim_env.md)
  environment) containing all operating model arrays.

- y:

  Integer. Last year to include; years `1:y` are retained.

- sim:

  Integer. Simulation replicate index to extract.

## Value

Named list with the following elements (all arrays have `y` in the year
dimension unless noted): `WAA`
`[n_pop x n_regions x y x n_seas x n_ages x n_sexes]`, `WAA_fish`
`[... x n_fish_fleets]`, `WAA_srv` `[... x n_srv_fleets]`, `MatAA`,
`SizeAgeTrans` (or `NULL`), `AgeingError` `[y x n_obs_ages x n_ages]`,
`use_conv_fish_tagging`, `conv_tag_release_indicator`,
`obs_conv_tag_fish_recap`, `conv_tagged_fish`, `conv_tagged_fish_attr`,
`n_tag_cohorts` (all `NULL` when tagging inactive), `ObsCatch`,
`ln_sigmaC`, `UseCatch`, `ObsCatch_pop`, `ln_sigmaC_pop`,
`UseCatch_pop`, `ObsDiscard`, `ln_sigmaD`, `UseDiscard`,
`ObsDiscard_pop`, `ln_sigmaD_pop`, `UseDiscard_pop`, `ObsFishIdx`,
`ObsFishIdx_SE`, `UseFishIdx`, `ObsFishIdx_pop`, `ObsFishIdx_pop_SE`,
`UseFishIdx_pop`, `ObsFishAgeComps`, `ISS_FishAgeComps`,
`UseFishAgeComps`, `ObsFishAgeComps_pop`, `ISS_FishAgeComps_pop`,
`UseFishAgeComps_pop`, `ObsFishLenComps` (or `NULL`), `ISS_FishLenComps`
(or `NULL`), `UseFishLenComps`, `ObsFishLenComps_pop` (or `NULL`),
`ISS_FishLenComps_pop` (or `NULL`), `UseFishLenComps_pop`,
`ObsFishAgeComps_discard`, `ISS_FishAgeComps_discard`,
`UseFishAgeComps_discard`, `ObsFishAgeComps_discard_pop`,
`ISS_FishAgeComps_discard_pop`, `UseFishAgeComps_discard_pop`,
`ObsFishLenComps_discard` (or `NULL`), `ISS_FishLenComps_discard` (or
`NULL`), `UseFishLenComps_discard`, `ObsFishLenComps_discard_pop` (or
`NULL`), `ISS_FishLenComps_discard_pop` (or `NULL`),
`UseFishLenComps_discard_pop`, `ObsSrvIdx`, `ObsSrvIdx_SE`, `UseSrvIdx`,
`ObsSrvIdx_pop`, `ObsSrvIdx_pop_SE`, `UseSrvIdx_pop`, `ObsSrvAgeComps`,
`ISS_SrvAgeComps`, `UseSrvAgeComps`, `ObsSrvAgeComps_pop`,
`ISS_SrvAgeComps_pop`, `UseSrvAgeComps_pop`, `ObsSrvLenComps` (or
`NULL`), `ISS_SrvLenComps` (or `NULL`), `UseSrvLenComps`,
`ObsSrvLenComps_pop` (or `NULL`), `ISS_SrvLenComps_pop` (or `NULL`),
`UseSrvLenComps_pop`.

## Details

Population-specific arrays (`ObsCatch_pop`, `ObsFishIdx_pop`,
`ObsFishAgeComps_pop`, `ObsFishLenComps_pop`, `ObsSrvIdx_pop`,
`ObsSrvAgeComps_pop`, `ObsSrvLenComps_pop`) are extracted when present
in `sim_env`; corresponding `Use*_PopSpec` flags are derived
automatically.

Input sample sizes for age and length compositions are extracted for
both aggregate and population-specific data sources, including retained
(`ISS_FishAgeComps`, `ISS_FishLenComps`, `ISS_FishAgeComps_pop`,
`ISS_FishLenComps_pop`, `ISS_SrvAgeComps`, `ISS_SrvLenComps`,
`ISS_SrvAgeComps_pop`, `ISS_SrvLenComps_pop`) and discard
(`ISS_FishAgeComps_discard`, `ISS_FishLenComps_discard`,
`ISS_FishAgeComps_discard_pop`, `ISS_FishLenComps_discard_pop`) data
sources.

Length composition outputs (`ObsFishLenComps`, `ObsSrvLenComps`, and
their population-specific and discard variants) and `SizeAgeTrans` are
`NULL` when no size-age transition matrix is present in `sim_env`.
Tagging outputs are `NULL` when `use_conv_fish_tagging = 0`; otherwise,
only cohorts with release years in `1:y` are retained.

## See also

[`Setup_Mod_Biologicals`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Biologicals.md),
[`Setup_Mod_Catch_and_F`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Catch_and_F.md),
[`Setup_Mod_SrvIdx_and_Comps`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_SrvIdx_and_Comps.md),
[`Setup_Mod_Tagging`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Tagging.md),
[`Simulate_Pop_Static`](https://chengmatt.github.io/SPoRC/dev/reference/Simulate_Pop_Static.md),
[`Setup_sim_env`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_sim_env.md)
