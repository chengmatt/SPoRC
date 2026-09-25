# Set up the conventional tagging module for model fitting

Sets the release cohorts and recapture data, the tag likelihood, the
mixing period, the release platform, which dims are attended and how
they pool, and the reporting rate blocks, sharing and priors. Call after
[`Setup_Mod_Dim`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Dim.md).

## Usage

``` r
Setup_Mod_Tagging(
  input_list,
  use_conv_fish_tagging = rep(0, input_list$data$n_fish_fleets),
  conv_tag_release_indicator = NULL,
  conv_tag_max_liberty = 0,
  conv_tagged_fish = NA,
  obs_conv_tag_fish_recap = NA,
  conv_fish_tag_like = NA,
  conv_tag_mixing_period = 1,
  conv_tag_t_tagging = 1,
  use_conv_tag_fishrep_prior = 0,
  conv_tag_fishrep_prior = NULL,
  conv_tag_pop_pool = as.list(1:input_list$data$n_pop),
  conv_tag_age_pool = as.list(seq_along(input_list$data$ages)),
  conv_tag_sex_pool = as.list(1:input_list$data$n_sexes),
  init_conv_tag_mort_spec = NULL,
  conv_tag_shed_spec = NULL,
  conv_tagrep_spec = "fix",
  conv_tag_fish_reporting_blocks = NULL,
  conv_fish_tag_attr = "p_a_s",
  conv_tag_release_platform = NULL,
  ...
)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, `$map` and `$verbose`.

- use_conv_fish_tagging:

  Integer vector `[n_fish_fleets]` (0/1) of whether each fleet's tagging
  data enter the likelihood. Default `0`.

- conv_tag_release_indicator:

  Integer matrix `[n_conv_tag_cohorts × 3]` of the release region, year
  and season of each cohort. Required when any
  `use_conv_fish_tagging = 1`. Default `NULL`.

- conv_tag_max_liberty:

  Integer maximum years at liberty in the likelihood; later recaptures
  are ignored. Must exceed `0` when tagging is active. Default `0`.

- conv_tagged_fish:

  Array `[n_conv_tag_cohorts × n_pop × n_ages × n_sexes]` of fish
  released per cohort. Dims not attended in `conv_fish_tag_attr` take
  all their fish in index 1 and zero elsewhere. Default `NA`.

- obs_conv_tag_fish_recap:

  Observed recaptures
  `[conv_tag_max_liberty × n_seas × n_conv_tag_cohorts × n_pop × n_regions × n_ages × n_sexes × n_fish_fleets]`.
  Default `NA`.

- conv_fish_tag_like:

  Tag recapture likelihood: `"Poisson"`, `"NegBin"`,
  `"Multinomial_Release"`, `"Multinomial_Recapture"`,
  `"Dirichlet-Multinomial_Release"` or
  `"Dirichlet-Multinomial_Recapture"`, stored as `0`-`5`. Default `NA`.

- conv_tag_mixing_period:

  Integer years, or seasons in a seasonal model, after release before
  recaptures contribute to the likelihood, allowing the tags to mix
  before they inform movement. Default `1`.

- conv_tag_t_tagging:

  Numeric scalar or vector `[n_conv_tag_cohorts]` in \\\[0, 1\]\\, the
  fraction of the season remaining at release: `1` the start of the
  season, `0.5` mid-season, `0` the end. A scalar is recycled. Default
  `1`.

- use_conv_tag_fishrep_prior:

  Integer (0/1) for priors on the reporting rates. Default `0`.

- conv_tag_fishrep_prior:

  Data frame with columns `region`, `block`, `fleet`, `mu`, `sd` and
  `type`. Default `NULL`.

- conv_tag_pop_pool, conv_tag_age_pool, conv_tag_sex_pool:

  Lists of integer vectors defining the population, age and sex pooling
  groups for the tagging likelihood. Use `list(1:n)` for a dim that is
  not attended; custom groupings such as `list(1:5, 6:10)` work for an
  attended one. A structure inconsistent with `conv_fish_tag_attr` warns
  and is overridden to a single group. Default one group per level.

- init_conv_tag_mort_spec, conv_tag_shed_spec:

  `"fix"`, `"est_shared"` or `"est_all"`: whether the initial
  tag-induced mortality and the chronic shedding rate are held at their
  starting values, estimated as one value across release events, or
  estimated per event. See
  [`do_conv_init_tag_mort_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_conv_init_tag_mort_mapping.md)
  and
  [`do_conv_tag_shed_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_conv_tag_shed_mapping.md).
  Default `NULL`.

- conv_tagrep_spec:

  Sharing structure for `conv_tag_fish_reporting_pars`
  `[n_regions × max_tagrep_blocks × n_fish_fleets]`. Default `"fix"`,
  which warns. See
  [`do_conv_tag_fish_reporting_pars_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_conv_tag_fish_reporting_pars_mapping.md).

- conv_tag_fish_reporting_blocks:

  Character vector of time blocks for the fishery reporting rates, each
  `"none_Region_r_Fleet_f"` or `"Block_b_Year_y1-y2_Region_r_Fleet_f"`
  with `"terminal"` allowed as the end year. Parsed into an
  `[n_regions × n_years × n_fish_fleets]` array. `NULL` (default) gives
  one constant block per region and fleet.

- conv_fish_tag_attr:

  Character scalar or vector `[n_conv_tag_cohorts]` naming which dims
  are resolved at release, built from `"p"` (population), `"a"` (age)
  and `"s"` (sex) joined by underscores: `"p_a_s"` (default), `"a_s"`,
  `"p_a"`, `"p_s"`, `"a"`, `"s"`, `"p"` or `"none"`. A scalar is
  recycled; a vector lets events resolve different dims. Region and
  fleet are always kept. An unattended dim takes all that event's fish
  in index 1 and is apportioned to full resolution by the release
  platform. A dim may only be split into several pooling groups if it is
  attended in every release event; otherwise its pooling argument is
  overridden to a single group with a warning.

- conv_tag_release_platform:

  Character matrix `[n_conv_tag_cohorts × 2]` of the release platform
  and fleet index per cohort, in the format
  [`Setup_Sim_Tagging`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Tagging.md)
  uses. Default `NULL`.

- ...:

  Optional starting values: `ln_init_conv_tag_mort` and
  `ln_conv_tag_shed`, each a scalar or a vector `[n_conv_tag_cohorts]`
  with a scalar recycled, default `-1000`; `ln_conv_fish_tag_theta`,
  default `0`; and `conv_tag_fish_reporting_pars`
  `[n_regions × max_tagrep_blocks × n_fish_fleets]`, default `0` on the
  logit scale, about a 0.5 reporting probability, with inactive fleet
  slots overwritten to `-1000`.

## Value

`input_list` with the tagging configuration in `$data`, the starting
values in `$par` for `ln_init_conv_tag_mort` and `ln_conv_tag_shed`
(each of length `n_conv_tag_cohorts`, or 1 when tagging is inactive),
`ln_conv_fish_tag_theta` and `conv_tag_fish_reporting_pars`, and the
factor maps for all four in `$map`.

## See also

Other Model Setup:
[`Setup_Mod_Biologicals()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Biologicals.md),
[`Setup_Mod_Catch_and_F()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Catch_and_F.md),
[`Setup_Mod_Dim()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Dim.md),
[`Setup_Mod_Discard_Comps()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Discard_Comps.md),
[`Setup_Mod_FishIdx_and_Comps()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_FishIdx_and_Comps.md),
[`Setup_Mod_Fishsel_and_Q()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Fishsel_and_Q.md),
[`Setup_Mod_Movement()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Movement.md),
[`Setup_Mod_Rec()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Rec.md),
[`Setup_Mod_Retsel()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Retsel.md),
[`Setup_Mod_SrvIdx_and_Comps()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_SrvIdx_and_Comps.md),
[`Setup_Mod_Srvsel_and_Q()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Srvsel_and_Q.md),
[`Setup_Mod_Weighting()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Weighting.md)
