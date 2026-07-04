# Generate fishery catches, compositions, and indices in simulation

Applies Baranov's catch equation to compute retained catch-at-age
(`CAA`) and dead discard catch-at-age (`DAA`) for all populations,
regions, seasons, and fleets, derives catch-at-length (`CAL` and `DAL`)
when a size-age transition matrix is available, and generates observed
catch and discard indices (with lognormal error), fishery abundance or
biomass indices, and age and length composition samples for both
retained and discarded catch. Composition sampling calls
[`simulate_comps`](https://chengmatt.github.io/SPoRC/dev/reference/simulate_comps.md)
and respects the likelihood type (`comp_fishage_like`,
`comp_fishlen_like`) and aggregation structure (`FishAgeComps_Type`,
`FishLenComps_Type`) specified in `sim_env`.

## Usage

``` r
generate_fishery_catch_comp_idx(y, sim, sim_env)
```

## Arguments

- y:

  Integer. Year index.

- sim:

  Integer. Simulation replicate index.

- sim_env:

  Simulation environment created by
  [`Setup_sim_env`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_sim_env.md).
  Modified in place. The following elements are updated:

  `CAA`, `DAA`

  :   Retained and dead discard catch-at-age for all populations,
      regions, seasons, and fleets.

  `CAL` and `DAL`

  :   Retained and dead discard catch-at-length if `SizeAgeTrans` is
      present.

  `TrueCatch`, `ObsCatch`

  :   Regional retained catch indices (abundance or biomass).

  `TrueCatch_pop`, `ObsCatch_pop`

  :   Population-specific retained catch indices.

  `TrueDiscard`, `ObsDiscard`

  :   Regional discard indices (abundance, biomass, or fraction).

  `TrueDiscard_pop`, `ObsDiscard_pop`

  :   Population-specific discard indices.

  `TrueFishIdx`, `ObsFishIdx`

  :   Regional fishery indices (abundance or biomass).

  `TrueFishIdx_pop`, `ObsFishIdx_pop`

  :   Population-specific fishery indices.

  `ObsFishAgeComps`, `ObsFishAgeComps_pop`

  :   Observed retained fishery age compositions.

  `ObsFishLenComps`, `ObsFishLenComps_pop`

  :   Observed retained fishery length compositions if `SizeAgeTrans` is
      available.

  `ObsFishAgeComps_discard`, `ObsFishAgeComps_discard_pop`

  :   Observed discard fishery age compositions.

  `ObsFishLenComps_discard`, `ObsFishLenComps_discard_pop`

  :   Observed discard fishery length compositions if `SizeAgeTrans` is
      available.

  `ISS_FishAgeComps`, `ISS_FishAgeComps_pop`, `ISS_FishLenComps`, `ISS_FishLenComps_pop`, `ISS_FishAgeComps_discard`, `ISS_FishAgeComps_discard_pop`, `ISS_FishLenComps_discard`, `ISS_FishLenComps_discard_pop`

  :   Effective sample sizes for retained and discard age and length
      compositions.

## Value

`invisible(NULL)`. All modifications are made by reference within
`sim_env`.

## Details

Composition draws are skipped for fleet-season cells with zero fishing
mortality (`Fmort = 0`). Discard composition draws are additionally
skipped when retention selectivity is fully 1 for the fleet-region-year-
season cell (i.e., no discarding occurs). Discard indices support four
unit types: abundance (`discard_units = 0`), biomass (`1`), abundance
fraction (`2`), and biomass fraction (`3`).

When `ISS_FishAgeComps_fill = "F_pattern"` and feedback is active,
sample sizes for retained and discard compositions in the current and
prior years are updated via
[`predict_sim_fish_iss_fmort`](https://chengmatt.github.io/SPoRC/dev/reference/predict_sim_fish_iss_fmort.md)
(scaled by fishing mortality) before sampling.

For each combination of season, region, and fleet, the function:

1.  Applies Baranov's catch equation to compute retained and dead
    discard catch-at-age.

2.  Converts catch-at-age to catch-at-length if `SizeAgeTrans` is
    available.

3.  Calculates true regional and population-specific catch, discard, and
    fishery indices.

4.  Applies lognormal observation error to generate observed indices.

5.  Simulates retained age and length compositions using
    [`simulate_comps`](https://chengmatt.github.io/SPoRC/dev/reference/simulate_comps.md),
    skipping fleet-season cells with zero fishing mortality.

6.  Simulates discard age and length compositions when retention
    selectivity is not fully 1.
