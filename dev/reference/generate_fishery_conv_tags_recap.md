# Generate conventional tag recaptures from fisheries in simulation

For each tag cohort (`tc`) and recovery season (`rseas`) in year `y`,
advances available tagged fish through movement and mortality, applies
Baranov's equation to compute predicted recaptures
(`pred_conv_tag_fish_recap`), and draws observed recaptures
(`obs_conv_tag_fish_recap`) via
[`simulate_conv_tag_fish_recaptures`](https://chengmatt.github.io/SPoRC/dev/reference/simulate_conv_tag_fish_recaptures.md).
Cohorts not yet released, already at `conv_tag_max_liberty`, or with
release year in the future are silently skipped.

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

  Simulation environment created by
  [`Setup_sim_env`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_sim_env.md).
  Modified in place. The following elements are updated:

  `conv_tag_fish_avail`

  :   Available tagged fish at age/region/season/fleet for each cohort.

  `pred_conv_tag_fish_recap`

  :   Predicted conventional tag recaptures by fleet, region, season,
      and cohort.

  `obs_conv_tag_fish_recap`

  :   Observed conventional tag recaptures after sampling error.

  `conv_tag_fish_surv`

  :   Surviving tagged fish after mortality and movement.

  `conv_tag_fish_reported`

  :   Reporting-adjusted recapture counts by fleet and region.

## Value

`invisible(NULL)`. All modifications are made by reference within
`sim_env`.

## Details

Total fishing mortality entering Z is decomposed into retained (\\F
\cdot s\_{\text{fish}} \cdot s\_{\text{ret}}\\) and dead discard (\\F
\cdot s\_{\text{fish}} \cdot (1 - s\_{\text{ret}}) \cdot \text{dmr}\\)
components, consistent with
[`apply_pop_dy`](https://chengmatt.github.io/SPoRC/dev/reference/apply_pop_dy.md).
Predicted recaptures use only the retained component in the Baranov
numerator, reflecting that tags are recovered from retained catch only.

At initial release (`ry = 1`, `rseas = tseas`), tags are placed into
`conv_tag_fish_avail[1, rseas, tc, ...]` after discounting for initial
tag-induced mortality (`ln_init_conv_tag_mort[tc]`). When
`conv_tag_t_tagging[tc] < 1`, total mortality is scaled by the fraction
of the season remaining at release for that cell only. Chronic shedding
(`ln_conv_tag_shed[tc]`) enters the total mortality rate alongside
natural and fishing mortality. `conv_tag_t_tagging`,
`ln_init_conv_tag_mort`, and `ln_conv_tag_shed` are each vectors of
length `n_tag_rel_events`, indexed by release event (`tc`), so timing,
initial mortality, and shedding can differ across release cohorts. At
the end of each season, survivors advance to the next season or the next
year's first season with plus-group accumulation. Tag reporting rates
from `conv_tag_fish_reporting` are applied fleet- and
region-specifically.

Tagged fish dynamics follow the same seasonal progression logic as the
population projection, including natural mortality, fishing mortality,
movement, and tag shedding. Recaptures are computed only from the
retained catch component, consistent with tag return processes. Cohorts
exceeding `conv_tag_max_liberty` are removed from the active tracking
pool.
