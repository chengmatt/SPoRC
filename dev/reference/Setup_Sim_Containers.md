# Initialise output containers for the operating model simulation

Allocates and appends zero-initialised arrays to `sim_list` for all
biological, fishery, and survey quantities tracked during simulation.
This function should be called after
[`Setup_Sim_Dim`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Dim.md)
and before any operating model dynamics are run. All arrays are
pre-allocated and populated in subsequent simulation steps.

## Usage

``` r
Setup_Sim_Containers(sim_list)
```

## Arguments

- sim_list:

  A simulation list returned by
  [`Setup_Sim_Dim`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Dim.md).
  Dimension elements (e.g., `n_pop`, `n_regions`, `n_yrs`, `n_seas`,
  `n_ages`, `n_sexes`, `n_sims`, `n_fish_fleets`, `n_srv_fleets`,
  `n_obs_ages`, `n_lens`) are used to size all output containers.

## Value

The input `sim_list` with the following zero-initialised arrays added:

\*\*Biological containers\*\*

- `$NAA`:

  Numbers-at-age
  `[n_pop × n_regions × (n_yrs+1) × n_seas × n_ages × n_sexes × n_sims]`.
  The `+1` year dimension stores initial conditions and propagates
  population state through the final year.

- `$NAA_bef`, `$NAA_aft`:

  Numbers-at-age immediately before and after fishing mortality is
  applied; same dimensions as `$NAA`.

- `$NAA0`:

  Unfished numbers-at-age; same dimensions as `$NAA`. Used to compute
  dynamic \\B_0\\ reference quantities.

- `$ZAA`:

  Total instantaneous mortality-at-age
  `[n_pop × n_regions × n_yrs × n_seas × n_ages × n_sexes × n_sims]`.

- `$Rec`:

  Recruitment `[n_pop × n_regions × n_yrs × n_sims]`.

- `$SSB`:

  Spawning stock biomass `[n_pop × n_regions × n_yrs × n_sims]`.

- `$eff_SSB`:

  Effective (population-aggregated) spawning stock biomass
  `[n_pop × n_yrs × n_sims]`.

- `$Dynamic_SSB0`:

  Dynamic unfished spawning stock biomass
  `[n_pop × n_regions × n_yrs × n_sims]`.

- `$Total_Biom`:

  Total biomass `[n_pop × n_regions × n_yrs × n_sims]`.

- `$ln_RecDevs`:

  Log-scale recruitment deviations
  `[n_pop × n_regions × n_yrs × n_sims]`.

- `$ln_InitDevs`:

  Log-scale initial age-structure deviations
  `[n_pop × n_regions × (n_ages - 1) × n_sexes × n_sims]`.

\*\*Fishery containers\*\*

\*Retained catch and indices (aggregated)\*

- `$ObsCatch`, `$TrueCatch`:

  Observed and true catch
  `[n_regions × n_yrs × n_seas × n_fish_fleets × n_sims]`.

- `$ObsFishIdx`, `$TrueFishIdx`:

  Observed and true fishery CPUE index; same dimensions as `$ObsCatch`.

- `$ObsFishAgeComps`:

  Observed fishery age compositions
  `[n_regions × n_yrs × n_seas × n_obs_ages × n_sexes × n_fish_fleets × n_sims]`.

- `$ObsFishLenComps`:

  Observed fishery length compositions
  `[n_regions × n_yrs × n_seas × n_lens × n_sexes × n_fish_fleets × n_sims]`.

\*Discards (aggregated)\*

- `$ObsDiscard`, `$TrueDiscard`:

  Observed and true discard
  `[n_regions × n_yrs × n_seas × n_fish_fleets × n_sims]`.

- `$ObsFishAgeComps_discard`:

  Observed discard age compositions; same dimensions as
  `$ObsFishAgeComps`.

- `$ObsFishLenComps_discard`:

  Observed discard length compositions; same dimensions as
  `$ObsFishLenComps`.

\*Population-specific quantities\*

- `$ObsCatch_pop`, `$TrueCatch_pop`:

  Observed and true catch
  `[n_pop × n_regions × n_yrs × n_seas × n_fish_fleets × n_sims]`.

- `$ObsFishIdx_pop`, `$TrueFishIdx_pop`:

  Observed and true fishery index; same dimensions as `$ObsCatch_pop`.

- `$ObsFishAgeComps_pop`:

  Observed fishery age compositions
  `[n_pop × n_regions × n_yrs × n_seas × n_obs_ages × n_sexes × n_fish_fleets × n_sims]`.

- `$ObsFishLenComps_pop`:

  Observed fishery length compositions
  `[n_pop × n_regions × n_yrs × n_seas × n_lens × n_sexes × n_fish_fleets × n_sims]`.

\*Discards (population-specific)\*

- `$ObsDiscard_pop`, `$TrueDiscard_pop`:

  Observed and true discard; same dimensions as `$ObsCatch_pop`.

- `$ObsFishAgeComps_discard_pop`:

  Observed discard age compositions; same dimensions as
  `$ObsFishAgeComps_pop`.

- `$ObsFishLenComps_discard_pop`:

  Observed discard length compositions; same dimensions as
  `$ObsFishLenComps_pop`.

\*True catch/discard at age and length\*

- `$CAA`, `$DAA`:

  Catch- and discard-at-age (true)
  `[n_pop × n_regions × n_yrs × n_seas × n_ages × n_sexes × n_fish_fleets × n_sims]`.

- `$CAL`, `$DAL`:

  Catch- and discard-at-length (true)
  `[n_pop × n_regions × n_yrs × n_seas × n_lens × n_sexes × n_fish_fleets × n_sims]`.

\*\*Survey containers\*\*

\*Aggregated\*

- `$ObsSrvIdx`, `$TrueSrvIdx`:

  Observed and true survey index
  `[n_regions × n_yrs × n_seas × n_srv_fleets × n_sims]`.

- `$ObsSrvAgeComps`:

  Observed survey age compositions
  `[n_regions × n_yrs × n_seas × n_obs_ages × n_sexes × n_srv_fleets × n_sims]`.

- `$ObsSrvLenComps`:

  Observed survey length compositions
  `[n_regions × n_yrs × n_seas × n_lens × n_sexes × n_srv_fleets × n_sims]`.

\*Population-specific\*

- `$ObsSrvIdx_pop`, `$TrueSrvIdx_pop`:

  Observed and true survey index
  `[n_pop × n_regions × n_yrs × n_seas × n_srv_fleets × n_sims]`.

- `$ObsSrvAgeComps_pop`:

  Observed survey age compositions
  `[n_pop × n_regions × n_yrs × n_seas × n_obs_ages × n_sexes × n_srv_fleets × n_sims]`.

- `$ObsSrvLenComps_pop`:

  Observed survey length compositions
  `[n_pop × n_regions × n_yrs × n_seas × n_lens × n_sexes × n_srv_fleets × n_sims]`.

\*True index at age and length\*

- `$SrvIAA`:

  Survey index-at-age (true)
  `[n_pop × n_regions × n_yrs × n_seas × n_ages × n_sexes × n_srv_fleets × n_sims]`.

- `$SrvIAL`:

  Survey index-at-length (true), same structure as `$SrvIAA` with
  `n_lens` replacing `n_ages`.

## See also

Other Simulation Setup:
[`Setup_Sim_Biologicals()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Biologicals.md),
[`Setup_Sim_Dim()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Dim.md),
[`Setup_Sim_Fishing()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Fishing.md),
[`Setup_Sim_Rec()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Rec.md),
[`Setup_Sim_Survey()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Survey.md),
[`Setup_Sim_Tagging()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Tagging.md),
[`Setup_sim_env()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_sim_env.md),
[`Simulate_Pop_Static()`](https://chengmatt.github.io/SPoRC/dev/reference/Simulate_Pop_Static.md),
[`run_annual_cycle()`](https://chengmatt.github.io/SPoRC/dev/reference/run_annual_cycle.md),
[`simulation_self_test()`](https://chengmatt.github.io/SPoRC/dev/reference/simulation_self_test.md)
