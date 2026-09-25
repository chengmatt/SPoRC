# Generate fishery catches, compositions, and indices in simulation

Takes retained and dead discard catch at age from Baranov's equation for
every population, region, season and fleet, converts them to catch at
length when a size-age key is available, and draws the observed catch,
discard and fishery indices under lognormal error together with the age
and length compositions of both retained and discarded catch. The
composition draws go through
[`simulate_comps`](https://chengmatt.github.io/SPoRC/dev/reference/simulate_comps.md)
and follow the likelihood and aggregation type set in `sim_env`.

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

  Simulation environment from
  [`Setup_sim_env`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_sim_env.md),
  modified in place. It gains the retained and dead discard catch at age
  `CAA` and `DAA`, their at-length counterparts `CAL` and `DAL` when a
  size-age key is present, the true and observed regional catch, discard
  and fishery indices with their population-specific counterparts, the
  observed retained and discard age and length compositions with theirs,
  and the input sample sizes of all eight composition data sources.

## Value

`invisible(NULL)`; everything is modified by reference within `sim_env`.

## Details

Composition draws are skipped where `Fmort = 0`, and the discard ones
also where retention selectivity is fully 1, so nothing is discarded.
Discard indices come in four units: abundance, biomass, abundance
fraction and biomass fraction. Under
`ISS_FishAgeComps_fill = "F_pattern"` with feedback active, the sample
sizes for the current and prior years are rescaled by
[`predict_sim_fish_iss_fmort`](https://chengmatt.github.io/SPoRC/dev/reference/predict_sim_fish_iss_fmort.md)
before sampling.
