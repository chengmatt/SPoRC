# Deterministic recruitment for a simulation year

Recruitment from the stock-recruit curve, before deviations, at the
spawning biomass in `SSB_vals`, so the recruits left by a trial F's
spawning biomass can be worked out without running the year.

## Usage

``` r
sim_det_recruitment(y, sim, sim_env, SSB_vals)
```

## Arguments

- y:

  Integer. Year index.

- sim:

  Integer. Simulation replicate index.

- sim_env:

  Simulation environment from
  [`Setup_sim_env`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_sim_env.md),
  modified in place: `$sbpr_table_cache`.

- SSB_vals:

  Array `[n_pop, n_regions, n_yrs]` of spawning biomass.

## Value

Array `[n_pop, n_regions]` of deterministic recruitment.
