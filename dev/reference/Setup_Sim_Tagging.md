# Set up conventional tagging dynamics for the operating model simulation

Populates `sim_list` with all conventional tagging inputs needed by the
operating model: tag release cohort definitions, release platform,
tagging timing, tag-induced mortality, chronic shedding, recapture
likelihood settings, fishery reporting rates, and zero-filled containers
for predicted and observed recaptures. Must be called after
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

  Numeric scalar. Constant number of tags released per release event.
  Cannot be specified simultaneously with `n_tags_rel_input`. Default
  `NULL`.

- n_tags_rel_input:

  Array specifying cohort-specific tag release numbers (e.g., varying by
  length, age, or other structure). Overrides `n_tags` when supplied.
  Cannot be specified simultaneously with `n_tags`. Default `NULL`.

- use_conv_fish_tagging:

  Integer (0/1). Whether conventional tag recaptures from fisheries are
  simulated. Default `0`.

- conv_tag_max_liberty:

  Integer. Maximum number of years-at-liberty tracked per tag cohort.
  Controls the size of recapture containers. Default: `n_ages / 2`.

- conv_tag_release_indicator:

  Data frame defining tag release cohorts. Required columns: `regions`
  (release region), `tag_years` (release year), `tag_seas` (release
  season). Default: all combinations of regions, years, and seasons from
  `sim_list`.

- conv_tag_release_platform:

  Character matrix `[n_conv_tag_cohorts × 2]` specifying the release
  platform and associated fleet index for each cohort. Rows must align
  with `conv_tag_release_indicator`. Column `platform` must be one of:

  `"population"`

  :   Tags released directly into the population, independent of any
      fleet or survey sampling process. `fleet` column should be `NA`.

  `"fishery"`

  :   Tags released through a fishery fleet. `fleet` column gives the
      fleet index.

  `"survey"`

  :   Tags released through a survey fleet. `fleet` column gives the
      fleet index.

  Default: `"survey"` platform with fleet `1` for every cohort.

- conv_tag_t_tagging:

  Numeric scalar or vector of length `n_tag_rel_events` (one value per
  row of `conv_tag_release_indicator`), each in \\\[0, 1\]\\. Fraction
  of the season remaining at the time of tag release for that release
  event. `1` = start of season; `0.5` = mid-season; `0` = end of season.
  A scalar is recycled to all release events. Default `1`.

- ln_init_conv_tag_mort:

  Log-scale initial tag-induced mortality applied at the moment of
  release. Numeric scalar or vector of length `n_tag_rel_events`, one
  value per release event. A scalar is recycled to all release events.
  Default `-1000` (approximately zero mortality).

- ln_conv_tag_shed:

  Log-scale annual chronic tag shedding rate. Numeric scalar or vector
  of length `n_tag_rel_events`, one value per release event. A scalar is
  recycled to all release events. Default `-1000` (approximately no
  shedding).

- conv_fish_tag_attr:

  Character scalar or character vector of length `n_tag_rel_events`
  specifying which biological dimensions are resolved at release (and
  hence retained at recapture) for each tag release event. A scalar is
  recycled to every event; a vector lets different events resolve
  different dimensions (e.g. event 1 known at `"p_a_s"`, event 2 only at
  `"p_a"`). Each element is built from any combination of `"p"`
  (population), `"a"` (age), and `"s"` (sex), joined by underscores.
  Region and fleet are always retained. Valid values: `"p_a_s"`,
  `"a_s"`, `"p_a"`, `"p_s"`, `"a"`, `"s"`, `"p"`, `"none"`. Default
  `"p_a_s"`.

- conv_tag_fish_reporting_input:

  Fishery tag reporting rate array
  `[n_regions × n_yrs × n_fish_fleets × n_sims]`. Values represent the
  probability that a recaptured tag is reported and must be in \\\[0,
  1\]\\. Default: `0.5` for all cells.

- conv_fish_tag_like:

  Integer or character scalar specifying the tag recapture likelihood.
  Default `0` (Poisson). Options:

  `0`/`"Poisson"`

  :   Poisson.

  `1`/`"NegBin"`

  :   Negative binomial.

  `2`/`"Multinomial_Release"`

  :   Multinomial conditioned on releases.

  `3`/`"Multinomial_Recapture"`

  :   Multinomial conditioned on recaptures.

  `4`/`"Dirichlet-Multinomial_Release"`

  :   Dirichlet-multinomial conditioned on releases.

  `5`/`"Dirichlet-Multinomial_Recapture"`

  :   Dirichlet-multinomial conditioned on recaptures.

- ln_conv_fish_tag_theta:

  Log-scale overdispersion parameter for negative-binomial (`1`) or
  Dirichlet-multinomial (`4`, `5`) likelihoods. Ignored for Poisson and
  multinomial likelihoods. Default `log(1)`.

## Value

The input `sim_list` with tagging-related fields appended: `$n_tags` or
`$n_tags_rel_input` (depending on which is provided),
`$conv_tag_max_liberty`, `$conv_tag_t_tagging` (length
`n_tag_rel_events`), `$ln_init_conv_tag_mort` (length
`n_tag_rel_events`), `$ln_conv_tag_shed` (length `n_tag_rel_events`),
`$conv_tag_release_indicator`, `$conv_tag_release_platform`,
`$n_tag_rel_events`, `$use_conv_fish_tagging`, `$conv_fish_tag_like`,
`$conv_fish_tag_attr`, `$ln_conv_fish_tag_theta`,
`$conv_tag_fish_reporting`, and zero-filled containers
`$conv_tagged_fish`, `$conv_tagged_fish_attr`, `$conv_tag_fish_avail`,
`$obs_conv_tag_fish_recap`, `$pred_conv_tag_fish_recap`.

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
