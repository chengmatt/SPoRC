# Generate survey indices and compositions in a simulation

Computes survey index-at-age (`SrvIAA`) for all populations using the
mid-survey abundance formula \\N \cdot s \cdot e^{-t\_{\text{srv}} Z}\\,
derives index-at-length (`SrvIAL`) when a size-age transition matrix is
available, generates observed survey indices (with lognormal error) as
abundance, biomass, or the recruitment deviations depending on
`srv_idx_type`, and draws age and length composition samples via
[`simulate_comps`](https://chengmatt.github.io/SPoRC/dev/reference/simulate_comps.md).

## Usage

``` r
generate_survey_comp_idx(y, sim, sim_env)
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

  `SrvIAA`

  :   Survey index-at-age for all populations.

  `SrvIAL`

  :   Survey index-at-length if `SizeAgeTrans` is present.

  `TrueSrvIdx`, `ObsSrvIdx`

  :   Aggregated survey index values.

  `TrueSrvIdx_pop`, `ObsSrvIdx_pop`

  :   Population-specific survey index values.

  `ObsSrvAgeComps`, `ObsSrvAgeComps_pop`

  :   Observed survey age compositions.

  `ObsSrvLenComps`, `ObsSrvLenComps_pop`

  :   Observed survey length compositions if `SizeAgeTrans` is
      available.

## Value

`invisible(NULL)`. All modifications are performed by reference within
`sim_env`.

## Details

This function loops over seasons, regions, and survey fleets for all
populations and replicates. It computes mid-period abundance, applies
survey selectivity, calculates true survey indices (abundance or
biomass), applies lognormal observation error, and simulates age and
length composition samples. Population-specific compositions are also
generated if requested.

## See also

[`simulate_comps`](https://chengmatt.github.io/SPoRC/dev/reference/simulate_comps.md),
[`Setup_sim_env`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_sim_env.md)
