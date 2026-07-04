# Set up the conventional tagging module for model fitting

Configures all conventional tagging components of the estimation model:
tag release cohort definitions, recapture data, tag recapture
likelihood, mixing period, release platform, dimension-attendance and
pooling structure, reporting rate time blocks and sharing, and optional
reporting rate priors. Delegates parameter mapping to four internal
helpers
([`do_conv_init_tag_mort_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_conv_init_tag_mort_mapping.md),
[`do_conv_tag_shed_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_conv_tag_shed_mapping.md),
[`do_conv_tag_theta_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_conv_tag_theta_mapping.md),
[`do_conv_tag_fish_reporting_pars_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_conv_tag_fish_reporting_pars_mapping.md)).
Must be called after
[`Setup_Mod_Dim`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Dim.md)
and before model compilation.

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
  conv_tag_age_pool = as.list(1:length(input_list$data$ages)),
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

  Named list with `$data`, `$par`, `$map`, and `$verbose` sublists, as
  returned by upstream setup functions.

- use_conv_fish_tagging:

  Integer vector `[n_fish_fleets]` (0/1). Whether conventional tagging
  data are included in the likelihood for each fishery fleet. Default:
  `0` for all fleets.

- conv_tag_release_indicator:

  Integer matrix `[n_conv_tag_cohorts × 3]` giving the release region,
  year, and season for each tag cohort. Required when any
  `use_conv_fish_tagging = 1`. Default `NULL`.

- conv_tag_max_liberty:

  Integer. Maximum years-at-liberty included in the likelihood;
  recaptures beyond this horizon are ignored. Must be \\\> 0\\ when
  tagging is active. Default `0`.

- conv_tagged_fish:

  Array `[n_conv_tag_cohorts × n_pop × n_ages × n_sexes]` of tagged fish
  released per cohort. Required when any `use_conv_fish_tagging = 1`.
  Dimensions not attended in `conv_fish_tag_attr` should have all fish
  placed into index 1 with remaining indices set to zero. Default `NA`.

- obs_conv_tag_fish_recap:

  Array of observed recaptures
  `[conv_tag_max_liberty × n_seas × n_conv_tag_cohorts × n_pop × n_regions × n_ages × n_sexes × n_fish_fleets]`.
  Required when any `use_conv_fish_tagging = 1`. Default `NA`.

- conv_fish_tag_like:

  Character string specifying the tag recapture likelihood. One of
  `"Poisson"`, `"NegBin"`, `"Multinomial_Release"`,
  `"Multinomial_Recapture"`, `"Dirichlet-Multinomial_Release"`,
  `"Dirichlet-Multinomial_Recapture"`. Converted to integer codes
  (`0`–`5`) before storage. Default `NA`.

- conv_tag_mixing_period:

  Integer. Minimum number of years (or seasons in seasonal models)
  post-release before recaptures contribute to the likelihood. Allows
  time for tags to mix within the population before informing movement
  estimation. Default `1`.

- conv_tag_t_tagging:

  Numeric scalar in \\\[0, 1\]\\. Fraction of the season remaining at
  tag release. `1` = start of season; `0.5` = mid-season; `0` = end of
  season. Default `1`.

- use_conv_tag_fishrep_prior:

  Integer (0/1). Whether priors are applied to reporting rate
  parameters. Default `0`.

- conv_tag_fishrep_prior:

  Data frame of prior specifications for reporting rates. Required
  columns: `region`, `block`, `fleet`, `mu`, `sd`, `type`. Ignored when
  `use_conv_tag_fishrep_prior = 0`. Default `NULL`.

- conv_tag_pop_pool:

  List of integer vectors defining population pooling groups for the
  tagging likelihood. When `"p"` is not attended in
  `conv_fish_tag_attr`, use `list(1:n_pop)`. If the pooling structure is
  inconsistent with `conv_fish_tag_attr`, a warning is issued and the
  structure is automatically overridden to a single group. Default:
  `as.list(1:n_pop)` (population-specific).

- conv_tag_age_pool:

  List of integer vectors defining age pooling groups. When `"a"` is not
  attended, use `list(1:n_ages)`. Custom groupings (e.g.,
  `list(1:5, 6:10)`) are supported when `"a"` is attended. Default:
  `as.list(1:n_ages)`.

- conv_tag_sex_pool:

  List of integer vectors defining sex pooling groups. When `"s"` is not
  attended, use `list(1:n_sexes)`. Default: `as.list(1:n_sexes)`.

- init_conv_tag_mort_spec:

  Character string (`"fix"` or `"est"`). Whether initial tag-induced
  mortality is fixed at its starting value or estimated. See
  [`do_conv_init_tag_mort_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_conv_init_tag_mort_mapping.md).
  Default `NULL`.

- conv_tag_shed_spec:

  Character string (`"fix"` or `"est"`). Whether chronic tag shedding is
  fixed or estimated. See
  [`do_conv_tag_shed_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_conv_tag_shed_mapping.md).
  Default `NULL`.

- conv_tagrep_spec:

  Character string. Sharing structure for reporting rate parameters
  `conv_tag_fish_reporting_pars`
  `[n_regions × max_tagrep_blocks × n_fish_fleets]`. See
  [`do_conv_tag_fish_reporting_pars_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_conv_tag_fish_reporting_pars_mapping.md)
  for full option descriptions. Default `"fix"` (a warning is issued if
  this was unintentional).

- conv_tag_fish_reporting_blocks:

  Character vector defining time blocks for fishery tag reporting rates.
  Each element follows one of:

  `"none_Region_r_Fleet_f"`

  :   Constant reporting rate for region `r` and fleet `f`.

  `"Block_b_Year_y1-y2_Region_r_Fleet_f"`

  :   Block `b` applies to years `y1`–`y2`. Use `"terminal"` for the end
      year to extend to the final model year.

  Parsed into an array `[n_regions × n_years × n_fish_fleets]`. If
  `NULL`, a single constant block is used for all region–fleet
  combinations. Default `NULL`.

- conv_fish_tag_attr:

  Character string specifying which biological dimensions are attended
  (resolved) in the recapture likelihood. Built from any combination of
  `"p"` (population), `"a"` (age), and `"s"` (sex), joined by
  underscores. Region and fleet are always retained. When a dimension is
  not attended, all released and recaptured fish must be placed into
  index 1 of that dimension, and the corresponding pooling argument must
  pool all indices into a single group. The likelihood then marginalises
  over the unattended dimension. Valid values: `"p_a_s"`, `"a_s"`,
  `"p_a"`, `"p_s"`, `"a"`, `"s"`, `"p"`, `"none"`. Default `"p_a_s"`.

- conv_tag_release_platform:

  Character matrix `[n_conv_tag_cohorts × 2]` specifying the release
  platform and fleet index per cohort. Same format as in
  [`Setup_Sim_Tagging`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Tagging.md).
  Default `NULL`.

- ...:

  Optional named starting values for parameters. Supported names and
  defaults: `ln_init_conv_tag_mort` (scalar, default `-1000`),
  `ln_conv_tag_shed` (scalar, default `-1000`), `ln_conv_fish_tag_theta`
  (scalar, default `0`), `conv_tag_fish_reporting_pars`
  `[n_regions × max_tagrep_blocks × n_fish_fleets]`, default `0` (logit
  scale ≈ 0.5 reporting probability; inactive fleet slots overwritten to
  `-1000`).

## Value

The input `input_list` with tagging configuration stored in `$data`
(`use_conv_fish_tagging`, `conv_tag_release_indicator`,
`n_conv_tag_cohorts`, `conv_tag_max_liberty`, `conv_tagged_fish`,
`obs_conv_tag_fish_recap`, `conv_fish_tag_like`,
`conv_tag_mixing_period`, `conv_tag_t_tagging`,
`use_conv_tag_fishrep_prior`, `conv_tag_fishrep_prior`,
`conv_tag_pop_pool`, `conv_tag_age_pool`, `conv_tag_sex_pool`,
`conv_tag_fish_reporting_blocks`, `conv_fish_tag_attr`,
`conv_tag_release_platform`); starting values in `$par` for
`ln_init_conv_tag_mort`, `ln_conv_tag_shed`, `ln_conv_fish_tag_theta`,
and `conv_tag_fish_reporting_pars`; and factor maps in `$map` for all
four parameter arrays.

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
