# Generate conventional tag recaptures from fisheries in simulation

For each tag cohort and recovery season in year `y`, advances the
available tagged fish through movement and mortality, applies Baranov's
equation for the predicted recaptures, and draws the observed ones
through
[`simulate_conv_tag_fish_recaptures`](https://chengmatt.github.io/SPoRC/dev/reference/simulate_conv_tag_fish_recaptures.md).
Cohorts not yet released, already at `conv_tag_max_liberty`, or released
in a future year are skipped.

## Usage

``` r
generate_fishery_conv_tags_recap(y, sim, sim_env)
```

## Arguments

- y:

  Integer. Year index.

- sim:

  Integer. Simulation replicate index.

- sim_env:

  Simulation environment from
  [`Setup_sim_env`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_sim_env.md),
  modified in place. It gains `conv_tag_fish_avail`, the tagged fish
  available by age, region, season and fleet for each cohort;
  `pred_conv_tag_fish_recap` and `obs_conv_tag_fish_recap`, the
  predicted recaptures and the observed ones after sampling error;
  `conv_tag_fish_surv`, the survivors after mortality and movement; and
  `conv_tag_fish_reported`, the reporting-adjusted counts by fleet and
  region.

## Value

`invisible(NULL)`; everything is modified by reference within `sim_env`.

## Details

Total fishing mortality entering \\Z\\ splits into retained, \\F \cdot
s\_{\text{fish}} \cdot s\_{\text{ret}}\\, and dead discards, \\F \cdot
s\_{\text{fish}} \cdot (1 - s\_{\text{ret}}) \cdot \text{dmr}\\, as in
[`apply_pop_dy`](https://chengmatt.github.io/SPoRC/dev/reference/apply_pop_dy.md),
and the Baranov numerator takes the retained component alone, since tags
come back from retained catch.

At release, tags enter `conv_tag_fish_avail[1, rseas, tc, ...]`
discounted by the initial tag-induced mortality
`ln_init_conv_tag_mort[tc]`, and when `conv_tag_t_tagging[tc] < 1` total
mortality is scaled by the fraction of the season remaining, for that
cell alone. Chronic shedding `ln_conv_tag_shed[tc]` joins natural and
fishing mortality in the total rate. All three are vectors of length
`n_tag_rel_events` indexed by release event, so timing, initial
mortality and shedding may differ across cohorts. At the end of a season
the survivors advance to the next season, or to the next year's first
season with plus group accumulation, and the reporting rates in
`conv_tag_fish_reporting` are applied by fleet and region.
