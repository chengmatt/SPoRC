# Set up conventional tagging dynamics for the operating model simulation

Sets the release cohorts, release platform, tagging timing, tag-induced
mortality, chronic shedding, the recapture likelihood and reporting
rates, and allocates the recapture containers. Call after
[`Setup_Sim_Dim`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Dim.md).

## Usage

``` r
Setup_Sim_Tagging(
  sim_list,
  n_tags = NULL,
  n_tags_rel_input = NULL,
  use_conv_fish_tagging = 0,
  conv_tag_max_liberty = sim_list$n_ages/2,
  conv_tag_release_indicator = expand.grid(regions = 1:sim_list$n_regions, tag_years =
    1:sim_list$n_yrs, tag_seas = 1:sim_list$n_seas),
  conv_tag_release_platform = matrix(c("survey", "1"), nrow =
    nrow(conv_tag_release_indicator), ncol = 2, byrow = TRUE, dimnames = list(NULL,
    c("platform", "fleet"))),
  conv_tag_t_tagging = 1,
  ln_init_conv_tag_mort = -1000,
  ln_conv_tag_shed = -1000,
  conv_fish_tag_attr = "p_a_s",
  conv_tag_fish_reporting_input = array(0.5, dim = c(sim_list$n_regions, sim_list$n_yrs,
    sim_list$n_fish_fleets, sim_list$n_sims)),
  conv_fish_tag_like = 0,
  ln_conv_fish_tag_theta = log(1)
)
```

## Arguments

- sim_list:

  Simulation list returned by
  [`Setup_Sim_Dim`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Dim.md).

- n_tags:

  Constant number of tags released per release event. Cannot be given
  alongside `n_tags_rel_input`. Default `NULL`.

- n_tags_rel_input:

  Array of cohort-specific release numbers, overriding `n_tags` and not
  given alongside it. Default `NULL`.

- use_conv_fish_tagging:

  Integer (0/1) for whether fishery tag recaptures are simulated.
  Default `0`.

- conv_tag_max_liberty:

  Integer maximum years at liberty tracked per cohort, which sizes the
  recapture containers. Default `n_ages / 2`.

- conv_tag_release_indicator:

  Data frame of release cohorts with columns `regions`, `tag_years` and
  `tag_seas`. Default every combination of the region, year and season
  in `sim_list`.

- conv_tag_release_platform:

  Character matrix `[n_conv_tag_cohorts × 2]` of the release platform
  and fleet index per cohort, aligned row by row with
  `conv_tag_release_indicator`. `"population"` releases tags into the
  population independently of any sampling process, with `fleet` set to
  `NA`; `"fishery"` and `"survey"` release them through the fleet the
  `fleet` column names. Default `"survey"` with fleet `1`.

- conv_tag_t_tagging:

  Numeric scalar or vector `[n_tag_rel_events]` in \\\[0, 1\]\\, the
  fraction of the season remaining at release: `1` the start of the
  season, `0.5` mid-season, `0` the end. A scalar is recycled. Default
  `1`.

- ln_init_conv_tag_mort:

  Log-scale tag-induced mortality applied at release, a scalar or a
  vector `[n_tag_rel_events]` with a scalar recycled. Default `-1000`,
  effectively none.

- ln_conv_tag_shed:

  Log-scale annual chronic shedding rate, a scalar or a vector
  `[n_tag_rel_events]` with a scalar recycled. Default `-1000`,
  effectively none.

- conv_fish_tag_attr:

  Character scalar or vector `[n_tag_rel_events]` naming which dims are
  resolved at release and so kept at recapture, built from `"p"`
  (population), `"a"` (age) and `"s"` (sex) joined by underscores:
  `"p_a_s"` (default), `"a_s"`, `"p_a"`, `"p_s"`, `"a"`, `"s"`, `"p"` or
  `"none"`. A scalar is recycled; a vector lets events resolve different
  dims. Region and fleet are always kept.

- conv_tag_fish_reporting_input:

  Fishery reporting rate array
  `[n_regions × n_yrs × n_fish_fleets × n_sims]` in \\\[0, 1\]\\, the
  probability a recaptured tag is reported. Default `0.5`.

- conv_fish_tag_like:

  Tag recapture likelihood: `0`/`"Poisson"` (default), `1`/`"NegBin"`,
  `2`/`"Multinomial_Release"`, `3`/`"Multinomial_Recapture"`,
  `4`/`"Dirichlet-Multinomial_Release"` or
  `5`/`"Dirichlet-Multinomial_Recapture"`.

- ln_conv_fish_tag_theta:

  Log-scale overdispersion for the negative binomial and the two
  Dirichlet-multinomial likelihoods, ignored by the others. Default
  `log(1)`.

## Value

`sim_list` with the tagging fields appended: `$n_tags` or
`$n_tags_rel_input`, `$conv_tag_max_liberty`, `$conv_tag_t_tagging`,
`$ln_init_conv_tag_mort` and `$ln_conv_tag_shed` (each of length
`n_tag_rel_events`), `$conv_tag_release_indicator`,
`$conv_tag_release_platform`, `$n_tag_rel_events`,
`$use_conv_fish_tagging`, `$conv_fish_tag_like`, `$conv_fish_tag_attr`,
`$ln_conv_fish_tag_theta` and `$conv_tag_fish_reporting`, plus the
zero-filled `$conv_tagged_fish`, `$conv_tagged_fish_attr`,
`$conv_tag_fish_avail`, `$obs_conv_tag_fish_recap` and
`$pred_conv_tag_fish_recap`.

## See also

Other Simulation Setup:
[`Setup_Sim_Biologicals()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Biologicals.md),
[`Setup_Sim_Containers()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Containers.md),
[`Setup_Sim_Dim()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Dim.md),
[`Setup_Sim_Fishing()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Fishing.md),
[`Setup_Sim_NAA_state()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_NAA_state.md),
[`Setup_Sim_Rec()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Rec.md),
[`Setup_Sim_Survey()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Survey.md),
[`Setup_sim_env()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_sim_env.md),
[`Simulate_Pop_Static()`](https://chengmatt.github.io/SPoRC/dev/reference/Simulate_Pop_Static.md),
[`run_annual_cycle()`](https://chengmatt.github.io/SPoRC/dev/reference/run_annual_cycle.md),
[`simulation_self_test()`](https://chengmatt.github.io/SPoRC/dev/reference/simulation_self_test.md)
