# Set up the recruitment module and associated processes

Configures all recruitment-related components of the estimation model:
stock-recruit relationship type and density-dependence structure,
Beverton-Holt steepness and priors, recruitment variability
(\\\sigma_R\\), annual and initial age-structure deviations, regional
and seasonal recruitment apportionment, spawning movement and stray
rates, sex ratio dynamics, equilibrium initialisation method, and the
recruitment bias ramp. Delegates parameter mapping to a family of
internal helpers
([`do_sigmaR_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_sigmaR_mapping.md),
[`do_RecDevs_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_RecDevs_mapping.md),
[`do_InitDevs_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_InitDevs_mapping.md),
[`do_h_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_h_mapping.md),
[`do_sexratio_pars_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_sexratio_pars_mapping.md),
[`do_rec_region_prop_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_rec_region_prop_mapping.md),
[`do_rec_seas_prop_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_rec_seas_prop_mapping.md)).
Must be called after
[`Setup_Mod_Dim`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Dim.md)
and
[`Setup_Mod_Biologicals`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Biologicals.md),
and before model compilation.

## Usage

``` r
Setup_Mod_Rec(
  input_list,
  rec_model,
  rec_dd = "global",
  rec_lag = 1,
  Use_h_prior = 0,
  h_prior = NULL,
  rec_region_prop_spec = NULL,
  use_rec_region_prop_prior = 0,
  rec_region_prop_prior = NULL,
  rec_seas_prop_spec = "fix",
  use_rec_seas_prop_prior = 0,
  rec_seas_prop_prior = NULL,
  use_fixed_rec_seas_prop = 1,
  fixed_rec_seas_prop = {
     rec_seas_prop = array(0, dim = c(input_list$data$n_pop,
    input_list$data$n_seas))
rec_seas_prop[, 1] <- 1
     rec_seas_prop
 },
  do_rec_bias_ramp = 0,
  bias_year = NA,
  max_bias_ramp_fct = 1,
  sigmaR_switch = 1,
  dont_est_recdev_last = 0,
  init_age_strc = 2,
  equil_init_age_strc = 1,
  init_F_prop = array(0, dim = c(input_list$data$n_regions, input_list$data$n_seas,
    input_list$data$n_fish_fleets)),
  sigmaR_spec = "est_all",
  InitDevs_spec = NULL,
  RecDevs_spec = NULL,
  h_spec = NULL,
  sgl_seas_spawning_movement = NA,
  t_spawn = 0,
  stray_rate_spec = "fix",
  stray_rate_blocks = paste0("none_Pop_", seq_len(input_list$data$n_pop)),
  use_fixed_stray_rate = if (stray_rate_spec != "fix") 0 else 1,
  fixed_stray_rate = array(0, dim = c(input_list$data$n_pop,
    length(input_list$data$years))),
  use_stray_rate_prior = 0,
  stray_rate_prior = NULL,
  spawn_seas = 1,
  sexratio_spec = "fix",
  sexratio_blocks = {
     grid <- expand.grid(region = 1:input_list$data$n_regions, pop
    = 1:input_list$data$n_pop)
     blks <- paste0("none_Pop_", grid$pop, "_Region_",
    grid$region)
     blks
 },
  use_rinit = 0,
  init_age_devs_shared = NULL,
  use_r0_prior = 0,
  r0_prior = NULL,
  ...
)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, `$map`, and `$verbose` sublists, as
  returned by upstream setup functions. Population, region, age, year,
  and season dimensions must already be defined in `$data`.

- rec_model:

  Character string (required). Stock-recruit relationship:

  `"mean_rec"`

  :   Fixed mean recruitment; no stock-recruit relationship. Steepness
      is automatically fixed and not estimated.

  `"bh_rec"`

  :   Beverton-Holt stock-recruit relationship.

- rec_dd:

  Density-dependence structure. Default `"global"`.

  `"local"`

  :   Independent stock-recruit relationship per population. Required
      when `n_pop > 1`.

  `"global"`

  :   Single pooled spawner-recruit relationship across all regions.
      Constrains `h_spec`, `RecDevs_spec`, and `InitDevs_spec` to shared
      or fixed options when `n_regions > 1`.

- rec_lag:

  Integer. Lag between spawning biomass and recruitment (in seasons).
  `1` (default) is the classic lagged case: recruitment uses SSB from
  `rec_lag` seasons prior and may enter in any season. `0` is age-0
  recruitment: recruitment uses the SAME year's SSB, and because that
  SSB isn't known until `spawn_seas` is reached, recruits may only enter
  in `spawn_seas` itself or a later season – when
  `use_fixed_rec_seas_prop = 1`, `fixed_rec_seas_prop` must be zero
  before `spawn_seas`; when estimated, this is enforced structurally via
  a restricted softmax (see
  [`do_rec_seas_prop_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_rec_seas_prop_mapping.md)).

- Use_h_prior:

  Integer (0/1). Whether normal priors on steepness are applied. Only
  relevant for `rec_model = "bh_rec"`. Default `0`.

- h_prior:

  Data frame of steepness prior parameters. Required columns: `pop`,
  `region`, `mu`, `sd`. Ignored when `Use_h_prior = 0`. Default `NULL`.

- rec_region_prop_spec:

  Character or `NULL`. Regional recruitment dispersal structure. Default
  `NULL` (estimate all proportions freely). See
  [`do_rec_region_prop_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_rec_region_prop_mapping.md)
  for full option descriptions including `"no_dispersal"`. Stored as
  `$data$rec_region_prop_spec`: `0` = full dispersal, `1` = no
  dispersal.

- use_rec_region_prop_prior:

  Integer (0/1). Whether Dirichlet priors are applied to regional
  recruitment proportions. Not valid when `n_regions = 1`. Default `0`.

- rec_region_prop_prior:

  Data frame of Dirichlet prior concentration parameters. Required
  columns: `pop` and `alpha`, where `alpha` is a list-column of
  length-`n_regions` vectors. Ignored when
  `use_rec_region_prop_prior = 0`. Default `NULL`.

- rec_seas_prop_spec:

  Character or `NULL`. Seasonal recruitment apportionment structure.
  Default `"fix"`. See
  [`do_rec_seas_prop_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_rec_seas_prop_mapping.md)
  for full option descriptions including `"est_shared_pop"`.

- use_rec_seas_prop_prior:

  Integer (0/1). Whether Dirichlet priors are applied to seasonal
  recruitment proportions. Not valid when `n_seas = 1`. When
  `rec_lag = 0` and `spawn_seas > 1`, the prior is evaluated only over
  seasons `spawn_seas:n_seas` (the seasons before `spawn_seas` are
  structurally zero, not estimated). Default `0`.

- rec_seas_prop_prior:

  Data frame of Dirichlet prior concentration parameters for seasonal
  proportions. Required columns: `pop` and `alpha`. Ignored when
  `use_rec_seas_prop_prior = 0`. Default `NULL`.

- use_fixed_rec_seas_prop:

  Integer (0/1). Whether fixed (non-estimated) seasonal proportions from
  `fixed_rec_seas_prop` are used. Automatically reset to `0` with a
  warning if `rec_seas_prop_spec` requests estimation. Default `1`.

- fixed_rec_seas_prop:

  Array `[n_pop x n_seas]`. Fixed seasonal recruitment proportions used
  when `use_fixed_rec_seas_prop = 1`. Default: all recruitment assigned
  to season 1. When `rec_lag = 0` and `spawn_seas > 1`, must be zero for
  every season before `spawn_seas` – an error is raised otherwise.

- do_rec_bias_ramp:

  Integer (0/1). Whether a recruitment bias ramp is applied to
  `ln_RecDevs` to account for reduced information in early and terminal
  years. Default `0`.

- bias_year:

  Numeric. Calendar year at which the bias ramp reaches its maximum
  correction. Only used when `do_rec_bias_ramp = 1`. Default `NA`.

- max_bias_ramp_fct:

  Numeric in \\\[0, 1\]\\. Maximum bias correction factor applied at
  `bias_year`. Default `1`.

- sigmaR_switch:

  Integer. Year index at which \\\sigma_R\\ switches from the
  early-period value (index 1) to the late-period value (index 2). If
  \\\leq 1\\, a single \\\sigma_R\\ is applied throughout. Default `1`.

- dont_est_recdev_last:

  Non-negative integer. Number of terminal years for which recruitment
  deviations are not estimated. Automatically overridden to `0` if
  `n_proj_yrs_devs > 0`, since projected deviation years are penalised
  toward the mean and are effectively estimated regardless. Default `0`.

- init_age_strc:

  Equilibrium initialisation method. Default `2`.

  `0`/`"iterative"`

  :   Iterates the population to approximate equilibrium. Slowest but
      most general.

  `1`/`"scalar_no_move"`

  :   Scalar geometric series assuming no movement.

  `2`/`"matrix"`

  :   Matrix geometric series incorporating movement. Recommended
      default for spatial models.

  `3`/`"scalar_plus_only"`

  :   Scalar geometric series with movement only in the plus group.

- equil_init_age_strc:

  Plus-group treatment during stochastic initialisation. Default `1`.

  `0`/`"equil"`

  :   Deterministic equilibrium; no `ln_InitDevs` are estimated.

  `1`/`"stoch_no_plus"`

  :   Stochastic deviations for all ages except the plus group.

  `2`/`"stoch_all"`

  :   Stochastic deviations for all ages including the plus group.

  `3`/`"stoch_shared_ages"`

  :   Stochastic deviations with user-defined age sharing via
      `init_age_devs_shared`. Deviations are estimated independently
      across all populations and regions, but ages sharing the same
      value in `init_age_devs_shared` are constrained to a single
      parameter. The plus group is not automatically fixed; include an
      `NA` in `init_age_devs_shared` at the plus-group position to fix
      it, or share it with the preceding age by repeating that index
      (e.g. `c(1:42, rep(42, 9))`). Requires `init_age_devs_shared` to
      be non-`NULL`.

- init_F_prop:

  Numeric array `[n_regions x n_seas x n_fish_fleets]`. Seasonal
  distribution of fishing mortality applied during equilibrium
  initialisation by fleet. Default: zero for all seasons and fleets.

- sigmaR_spec:

  Character. Estimation structure for \\\sigma_R\\, stored in
  `ln_sigmaR` `[2 x n_pop x n_regions]`, where index 1 = initial
  deviation period and index 2 = annual deviation period. Default
  `"est_all"`. See
  [`do_sigmaR_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_sigmaR_mapping.md)
  for full option descriptions.

- InitDevs_spec:

  Character or `NULL`. Sharing structure for initial age-structure
  deviations `ln_InitDevs` `[n_pop x n_regions x (n_ages - 1)]`. Default
  `NULL` (estimate all independently). See
  [`do_InitDevs_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_InitDevs_mapping.md)
  for full option descriptions.

- RecDevs_spec:

  Character or `NULL`. Sharing structure for annual recruitment
  deviations `ln_RecDevs` `[n_pop x n_regions x n_years]`. Default
  `NULL` (estimate all independently). See
  [`do_RecDevs_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_RecDevs_mapping.md)
  for full option descriptions.

- h_spec:

  Character or `NULL`. Sharing structure for Beverton-Holt steepness
  `steepness_h` `[n_pop x n_regions]`, parameterised in bounded logit
  space \\(0.2, 1)\\. Default `NULL` (estimate by population when
  `n_pop > 1`, by region when `n_pop = 1`). Ignored when
  `rec_model = "mean_rec"`. See
  [`do_h_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_h_mapping.md)
  for full option descriptions.

- sgl_seas_spawning_movement:

  Spawning movement array
  `[n_pop x n_regions x n_regions x n_years x n_ages x n_sexes]`. Each
  `[p, , r, y, a, s]` slice is a row-stochastic movement matrix giving
  the probability of fish from each origin region spawning in region
  `r`. If `NA` (default), 100% natal homing is assumed and the array is
  constructed internally.

- t_spawn:

  Numeric. Spawn timing as a fraction of the season elapsed before
  spawning. `0` (default) = spawning before any mortality; `1` =
  spawning after all mortality.

- stray_rate_spec:

  Character string. Estimation structure for `stray_rate_pars`
  `[n_pop x max_stray_blocks]`, parameterised on the logit scale.
  Ignored when `use_fixed_stray_rate = 1` or `n_pop = 1`. Default
  `"fix"`. Options:

  `"fix"`

  :   All parameters fixed at starting values (mapped to `NA`). Use this
      alongside `use_fixed_stray_rate = 0` to hold stray rates at a
      specified value without estimating.

  `"est_all"`

  :   Estimate independently per population x block. Produces one
      parameter per population per unique block.

  `"est_shared_pop"`

  :   Single parameter per block, shared across all populations.
      Requires identical block structures across all populations – an
      error is raised if block indices differ.

- stray_rate_blocks:

  Character vector of length `n_pop` defining the temporal block
  structure for stray rate parameters. Valid formats:

  `"none_Pop_x"`

  :   Constant stray rate for population `x` across all years (single
      block).

  `"Block_k_Year_a-b_Pop_x"`

  :   Block `k` applies to years `a` through `b` for population `x`. Use
      `"terminal"` in place of the end year to extend through the final
      model year.

  Default: a single constant block for every population. **Note:** stray
  rate is generally unidentifiable from fisheries data alone.
  Time-blocking is provided for completeness but regularisation via
  `use_stray_rate_prior` in the penalty setup is strongly recommended
  whenever `stray_rate_spec != "fix"`.

- use_fixed_stray_rate:

  Integer (0/1). Whether stray rates are supplied as a fixed external
  array (`fixed_stray_rate`) rather than estimated as model parameters.
  Default `1` (fixed), preserving existing behaviour. Set to `0` to
  estimate stray rates via `stray_rate_pars`.

- fixed_stray_rate:

  Array `[n_pop x n_years]`. Fixed stray rate values used when
  `use_fixed_stray_rate = 1`. Values should be in \\\[0, 1\]\\. Default:
  `0` (no straying) for all populations and years. Ignored when
  `use_fixed_stray_rate = 0`.

- use_stray_rate_prior:

  Integer (0/1). Whether Beta priors are applied to estimated stray rate
  parameters. Only relevant when `use_fixed_stray_rate = 0` and
  `n_pop > 1`. An error is raised if `use_stray_rate_prior = 1`
  alongside `use_fixed_stray_rate = 1` since `stray_rate_pars` would not
  be estimated. Default `0`.

- stray_rate_prior:

  Data frame of Beta prior parameters for stray rates. Required columns:
  `pop` (population index), `block` (block index matching
  `stray_rate_blocks`), `mu` (prior mean, in \\(0,1)\\), `sd` (prior
  standard deviation). One row per population x block combination to
  penalise. Ignored when `use_stray_rate_prior = 0`. Default `NULL`.

- spawn_seas:

  Integer. Season index in which spawning occurs. Default `1`.

- sexratio_spec:

  Character. Estimation structure for sex ratio parameters
  `sexratio_pars` `[n_pop x n_regions x n_blocks]`. Default `"fix"`. See
  [`do_sexratio_pars_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_sexratio_pars_mapping.md)
  for full option descriptions. Must be `"fix"` when `n_sexes = 1`.

- sexratio_blocks:

  Character vector defining temporal block structure for sex ratio
  parameters. One entry per population-region combination. Valid
  formats:

  `"none_Pop_x_Region_x"`

  :   Constant sex ratio for population `x` and region `x` (single block
      across all years).

  `"Block_k_Year_a-b_Pop_x_Region_x"`

  :   Block `k` applies to years `a` through `b`. Use `"terminal"` in
      place of the end year to extend through the final model year.

  Default: a single constant block for every population-region
  combination.

- use_rinit:

  Integer (0/1). Whether a separate initial recruitment scalar
  `ln_rinit` is used to initialise the population independently of the
  recruitment `ln_global_R0`. When `0` (default), `ln_rinit` is fixed
  and `ln_global_R0` governs both initialisation and recruitment. When
  `1`, both `ln_rinit` and `ln_global_R0` are estimated, with `ln_rinit`
  used exclusively for equilibrium initialisation and `ln_global_R0`
  used for the stock-recruit relationship.

- init_age_devs_shared:

  Integer vector of length `n_ages - 1` specifying an explicit
  parameter-sharing structure for `ln_InitDevs` along the age dimension.
  Each element gives the factor level assigned to that age position;
  positions sharing the same integer value are constrained to a single
  estimated parameter. Used in conjunction with
  `equil_init_age_strc = 3` (`"stoch_shared_ages"`), which activates
  user-defined age sharing while still estimating deviations
  independently across populations and regions. The sharing structure is
  also respected by `InitDevs_spec` options: `"est_shared_r"` applies
  the vector per population (with a population-level offset so pops
  remain independent), and `"est_shared_pop_r"` applies it globally (no
  offset, all pops and regions share the same parameters). A typical use
  case is replicating ADMB models where ages beyond the data plus group
  share the last estimated deviation, e.g. `c(1:42, rep(42, 9))` for a
  52-age model with 43 data ages, giving 42 free parameters. When `NULL`
  (default), age sharing follows the standard behaviour determined by
  `equil_init_age_strc` alone.

- use_r0_prior:

  Integer (0/1). Whether to apply a lognormal prior on R0 for any
  populations. Default 0.

- r0_prior:

  Data frame with columns `pop` (population index), `mu` (prior mean on
  natural scale), and `sd` (prior SD on log scale). Required when
  `use_r0_prior = 1`.

- ...:

  Optional named starting values for parameters. Any of: `ln_global_R0`
  `[n_pop]`, `ln_rinit` `[n_pop]`, `rec_region_prop_pars`
  `[n_pop x (n_regions - 1)]`, `rec_seas_prop_pars`
  `[n_pop x (n_seas - 1)]`, `steepness_h` `[n_pop x n_regions]` (bounded
  logit scale), `ln_InitDevs` `[n_pop x n_regions x (n_ages - 1)]`,
  `ln_RecDevs` `[n_pop x n_regions x n_years]`, `ln_sigmaR`
  `[2 x n_pop x n_regions]`, `sexratio_pars`
  `[n_pop x n_regions x n_blocks]`. Unspecified parameters use internal
  defaults.

## Value

The input `input_list` with all recruitment-related fields populated in
`$data` and `$par`, and factor maps constructed in `$map` for:
`rec_region_prop_pars`, `rec_seas_prop_pars`, `ln_sigmaR`,
`ln_InitDevs`, `ln_RecDevs`, `steepness_h`, `sexratio_pars`, and
`stray_rate_pars`. Character-coded inputs for `init_age_strc` and
`equil_init_age_strc` are converted to integer codes before storage.

## See also

Other Model Setup:
[`Setup_Mod_Biologicals()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Biologicals.md),
[`Setup_Mod_Catch_and_F()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Catch_and_F.md),
[`Setup_Mod_Dim()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Dim.md),
[`Setup_Mod_Discard_Comps()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Discard_Comps.md),
[`Setup_Mod_FishIdx_and_Comps()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_FishIdx_and_Comps.md),
[`Setup_Mod_Fishsel_and_Q()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Fishsel_and_Q.md),
[`Setup_Mod_Movement()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Movement.md),
[`Setup_Mod_Retsel()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Retsel.md),
[`Setup_Mod_SrvIdx_and_Comps()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_SrvIdx_and_Comps.md),
[`Setup_Mod_Srvsel_and_Q()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Srvsel_and_Q.md),
[`Setup_Mod_Tagging()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Tagging.md),
[`Setup_Mod_Weighting()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Weighting.md)
