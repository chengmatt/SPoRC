# Set up the recruitment module and associated processes

Sets the stock-recruit form and density dependence, steepness and its
prior, \\\sigma_R\\, annual and initial deviations, regional and
seasonal apportionment, spawning movement, stray rates, sex ratio, the
equilibrium initialization and the bias ramp. Call after
[`Setup_Mod_Dim`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Dim.md)
and
[`Setup_Mod_Biologicals`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Biologicals.md).

## Usage

``` r
Setup_Mod_Rec(
  input_list,
  rec_model,
  rec_dd = "global",
  rec_lag = 1,
  SR_ref_yr = 1,
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
  dont_pen_recdev_first = 0,
  init_age_strc = 2,
  equil_init_age_strc = 1,
  init_F_prop = array(0, dim = c(input_list$data$n_regions, input_list$data$n_seas,
    input_list$data$n_fish_fleets)),
  init_F_form = "prop",
  init_F_spec = "fix",
  sigmaR_spec = "est_all",
  InitDevs_spec = NULL,
  InitDevs_sex_spec = "est_shared_s",
  RecDevs_spec = NULL,
  RecDevs_model = "iid",
  RecDevs_rho_spec = "fix",
  RecDevs_rw_init_sigma = 5,
  RecDevs_pen_center = "fixed",
  Use_rec_level_pen = 0,
  rec_level_pen_sigma = 1,
  rec_level_pen_center = "own_mean",
  rec_level_pen_yrs = NULL,
  Use_init_sex_pen = 0,
  init_sex_pen_sigma = 1,
  sr_penalty = "none",
  sr_pen_sigma = 1,
  sr_pen_yrs = NULL,
  sr_R0_spec = "shared",
  InitDevs_pen_center = "fixed",
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
  R0_blocks = NULL,
  R0_ref_block = 1,
  use_r0_prior = 0,
  r0_prior = NULL,
  Use_rinit_pen = 0,
  rinit_pen_sd = 1,
  ...,
  ln_global_R0_spec = "est"
)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, `$map` and `$verbose`. Dimensions
  must already be set in `$data`.

- rec_model:

  Character, required. `"mean_rec"` (steepness fixed and not estimated),
  `"bh_rec"`, or `"ricker_rec"`. Steepness is not interchangeable
  between the last two, see
  [`Get_Det_Recruitment`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Det_Recruitment.md).

- rec_dd:

  Density dependence. `"local"` is one stock-recruit relationship per
  population, required when `n_pop > 1`. `"global"` (default) pools
  across regions and restricts `h_spec`, `RecDevs_spec` and
  `InitDevs_spec` to shared or fixed options when `n_regions > 1`.

- rec_lag:

  Integer lag in seasons between spawning biomass and recruitment. `1`
  (default) uses SSB from that many seasons prior. `0` is age-0
  recruitment on the same year's SSB, so recruits may only enter in
  `spawn_seas` or later.

- SR_ref_yr:

  Integer year index supplying every input to unfished spawning biomass
  per recruit, and so to `S0` and the curve's scale: weight-at-age,
  maturity, natural mortality, movement, stray rate, sex ratio, and what
  enters through `init_F`. `R0` is the exception and is always the
  year's own value. Default `1`. Ignored under `rec_model = "mean_rec"`.

- Use_h_prior:

  Integer (0/1). Normal priors on steepness. Default `0`.

- h_prior:

  Data frame with columns `pop`, `region`, `mu` and `sd`. Default
  `NULL`.

- rec_region_prop_spec:

  Character or `NULL`. Regional recruitment dispersal structure, default
  `NULL` (all estimated freely). Stored as `$data$rec_region_prop_spec`,
  `0` = full dispersal, `1` = none. See
  [`do_rec_region_prop_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_rec_region_prop_mapping.md).

- use_rec_region_prop_prior:

  Integer (0/1). Dirichlet priors on regional recruitment proportions.
  Not valid when `n_regions = 1`. Default `0`.

- rec_region_prop_prior:

  Data frame of Dirichlet concentrations with columns `pop` and `alpha`,
  a list-column of length-`n_regions` vectors. Default `NULL`.

- rec_seas_prop_spec:

  Character or `NULL`. Seasonal apportionment structure, default
  `"fix"`. See
  [`do_rec_seas_prop_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_rec_seas_prop_mapping.md).

- use_rec_seas_prop_prior:

  Integer (0/1). Dirichlet priors on seasonal proportions. Not valid
  when `n_seas = 1`. Evaluated only over `spawn_seas:n_seas` when
  `rec_lag = 0` and `spawn_seas > 1`. Default `0`.

- rec_seas_prop_prior:

  Data frame of Dirichlet concentrations with columns `pop` and `alpha`.
  Default `NULL`.

- use_fixed_rec_seas_prop:

  Integer (0/1). Whether `fixed_rec_seas_prop` is used. Reset to `0`
  with a warning if `rec_seas_prop_spec` estimates. Default `1`.

- fixed_rec_seas_prop:

  Array `[n_pop x n_seas]` of fixed seasonal proportions, default all
  recruitment in season 1. Must be zero before `spawn_seas` when
  `rec_lag = 0` and `spawn_seas > 1`.

- do_rec_bias_ramp:

  Integer (0/1). Whether a bias ramp is applied to `ln_RecDevs`. Under
  `0` every penalty is centered on the full \\-\sigma_R^2/2\\, under `1`
  the center follows the ramp. Default `0`.

- bias_year:

  Numeric calendar year at which the ramp reaches its maximum
  correction. Default `NA`.

- max_bias_ramp_fct:

  Numeric in \\\[0, 1\]\\, the maximum correction applied at
  `bias_year`. Default `1`.

- sigmaR_switch:

  Integer year index at which \\\sigma_R\\ switches from the early to
  the late value. \\\leq 1\\ (default) uses one value throughout.

- dont_est_recdev_last:

  Non-negative integer. Terminal years for which recruitment deviations
  are not estimated. Forced to `0` when `n_proj_yrs_devs > 0`, refused
  under `RecDevs_model = "dsem"`. Default `0`.

- dont_pen_recdev_first:

  Integer. How many leading years of recruitment deviations are
  estimated but left out of the penalty. `0` (default) penalizes every
  year.

- init_age_strc:

  Initialization method, default `2`. `0`/`"iterative"` iterates to
  approximate equilibrium, `1`/`"scalar_no_move"` is a scalar geometric
  series without movement, `2`/`"matrix"` is the matrix series with
  movement, and `3`/`"scalar_plus_only"` moves only the plus group; all
  four treat `ln_InitDevs` as multiplicative deviations from the
  equilibrium. `4`/`"free"` projects no equilibrium: ages 2 and older
  are `exp(ln_InitDevs)` apportioned by sex ratio, so the deviations are
  on the scale of numbers and `equil_init_age_strc` becomes a prior on
  log abundance.

- equil_init_age_strc:

  Plus-group treatment during stochastic initialization, default `1`.
  `0`/`"equil"` estimates no `ln_InitDevs`, `1`/`"stoch_no_plus"`
  estimates every age but the plus group, `2`/`"stoch_all"` estimates
  every age, `3`/`"stoch_shared_ages"` shares ages through
  `init_age_devs_shared` (non-`NULL` required; the plus group is not
  fixed automatically), and `4`/`"stoch_all_no_pen"` estimates every age
  and penalizes none, for pairing with `init_age_strc = "free"`.

- init_F_prop:

  Numeric array `[n_regions x n_seas x n_fish_fleets]`. Legacy
  interface. A non-zero value without `init_F_par` is converted to
  `ln_init_F = log(init_F_prop)` with `init_F_form = "prop"`. Prefer
  `init_F_par`. Default zero.

- init_F_form:

  What `init_F_par` means. `"prop"` (default) gives
  `init_F = exp(ln_init_F) * exp(ln_F_mean)`, a proportion of the
  estimated mean F. `"abs"` gives `init_F = exp(ln_init_F)`, independent
  of `ln_F_mean`. Use `"abs"` when bridging an assessment with a
  separate historical F; under `"prop"` the two collapse into one
  parameter and catch constrains only their product.

- init_F_spec:

  `"fix"` (default) or `"est"`, whether `init_F_par` is estimated. Sets
  only the mapping, so it combines freely with `init_F_form`. Refused
  under `init_age_strc = "free"`. The value comes from `init_F_par`
  `[n_regions x n_seas x n_fish_fleets]`, passed through `...`, on the
  logit scale under `"prop"` and the log scale under `"abs"`.

- sigmaR_spec:

  Character. Estimation structure for \\\sigma_R\\, stored in
  `ln_sigmaR` `[2 x n_pop x n_regions]` with index 1 the initial period
  and 2 the annual period. Default `"est_all"`, see
  [`do_sigmaR_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_sigmaR_mapping.md).

- InitDevs_spec:

  Character or `NULL`. Sharing structure for `ln_InitDevs`
  `[n_pop x n_regions x (n_ages - 1) x n_sexes]`. Default `NULL` (all
  independent), see
  [`do_InitDevs_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_InitDevs_mapping.md).

- InitDevs_sex_spec:

  `"est_shared_s"` (default) estimates one initial age deviation curve
  read by every sex; `"est_all"` gives each sex its own, pooling the
  level across sexes under an `"own_mean"` `InitDevs_pen_center`.
  Requires `n_sexes > 1`.

- RecDevs_spec:

  Character or `NULL`. Sharing structure for `ln_RecDevs`
  `[n_pop x n_regions x n_years]`. Default `NULL` (all independent), see
  [`do_RecDevs_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_RecDevs_mapping.md).

- RecDevs_model:

  Process error on `ln_RecDevs`. `"iid"` (default) is independent about
  the prior mean, `"rw"` centers each deviation on the previous one with
  a diffuse first year, `"ar1"` reverts toward zero at rate
  `RecDevs_rho` with a stationary first year, and `"dsem"` takes the
  density from the arrows in
  [`Setup_Mod_DSEM`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_DSEM.md).
  Steps span estimated years, not calendar years. `"rw"` and `"ar1"` are
  refused alongside `do_rec_bias_ramp = 1` or
  `RecDevs_pen_center = "own_mean"`; `"dsem"` reads `sigmaR` off the
  arrows, so `sigmaR_spec` other than `"fix"`,
  `dont_est_recdev_last > 0` and a nonzero ramp are refused.

- RecDevs_rho_spec:

  Sharing structure for `RecDevs_rho` `[n_pop x n_regions]`:
  `"est_all"`, `"est_shared_pop"`, `"est_shared_r"`,
  `"est_shared_pop_r"` or `"fix"` (default). Only read under
  `RecDevs_model = "ar1"`, see
  [`do_RecDevs_rho_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_RecDevs_rho_mapping.md).

- RecDevs_rw_init_sigma:

  Standard deviation given to year one of a random walk, which sets the
  level of the series. Default `5`; `NA` instead starts the walk at zero
  under its own sigma. Only read under `RecDevs_model = "rw"`.

- RecDevs_pen_center, InitDevs_pen_center:

  Where the recruitment and initial age deviation penalties are
  centered. `"fixed"` (default) centers on zero or the bias-corrected
  \\-\sigma_R^2/2\\; `"own_mean"` centers on the estimated deviations'
  own mean, leaving their level free to be set elsewhere. Cannot be
  combined with `do_rec_bias_ramp = 1`.

- Use_rec_level_pen:

  Integer (0/1). Whether the log recruitment series itself is penalized,
  separately from the deviation penalty. Default `0`.

- rec_level_pen_sigma:

  Numeric standard deviation of that penalty. A sum of squares with
  weight \\w\\ corresponds to \\1/\sqrt{2w}\\. Default `1`.

- rec_level_pen_center:

  `"own_mean"` (default) centers on the mean of the log recruitment
  series, `"fixed"` centers on zero.

- rec_level_pen_yrs:

  Years the penalty applies over, or `NULL` (default) for every year.

- Use_init_sex_pen:

  Integer (0/1). Whether each later sex's initial age deviations are
  tied to the first sex's by a Gaussian on their difference. Requires
  `n_sexes > 1` and `InitDevs_sex_spec = "est_all"`. Enters the
  objective unweighted. Default `0`.

- init_sex_pen_sigma:

  Numeric standard deviation of that tie. Default `1`.

- sr_penalty:

  Character. `"none"` (default), `"bh"` or `"ricker"`. Only valid under
  `rec_model = "mean_rec"`. Fits the curve as a likelihood on \\\log
  R_y - \log\widehat{R}\_y\\ without letting it generate recruitment.

- sr_pen_sigma:

  Numeric standard deviation of that residual.

- sr_pen_yrs:

  Years the stock-recruit penalty applies over, or `NULL` (default) for
  every year with a lagged spawning biomass. Naming a year without one
  is an error.

- sr_R0_spec:

  Character. `"shared"` (default) takes the curve's scale from
  `ln_global_R0`, `"est"` gives it its own `ln_sr_R0`, and `"rinit"`
  takes it from `ln_rinit` and requires `use_rinit = 1`.

- h_spec:

  Character or `NULL`. Sharing structure for `steepness_h`
  `[n_pop x n_regions]`, on a logit scale bounded to \\(0.2, 1)\\.
  Default `NULL`, ignored under `rec_model = "mean_rec"`. See
  [`do_h_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_h_mapping.md).

- sgl_seas_spawning_movement:

  Array `[n_pop x n_regions x n_regions x n_years x n_ages x n_sexes]`,
  each `[p, , r, y, a, s]` slice row-stochastic. `NA` (default) assumes
  complete natal homing and builds the array internally.

- t_spawn:

  Numeric fraction of the spawning season elapsed before spawning. `0`
  (default) spawns before any mortality, `1` after all of it.

- stray_rate_spec:

  Estimation structure for `stray_rate_pars`
  `[n_pop x max_stray_blocks]` on the logit scale. `"fix"` (default)
  holds every value, `"est_all"` estimates per population and block, and
  `"est_shared_pop"` gives one parameter per block shared across
  populations, which requires identical block structures. Ignored when
  `use_fixed_stray_rate = 1` or `n_pop = 1`.

- stray_rate_blocks:

  Character vector of length `n_pop`, either `"none_Pop_x"` for one
  block or `"Block_k_Year_a-b_Pop_x"`, with `"terminal"` allowed as the
  end year. Default one block each. Stray rate is generally
  unidentifiable from fisheries data alone, so use
  `use_stray_rate_prior` whenever `stray_rate_spec != "fix"`.

- use_fixed_stray_rate:

  Integer (0/1). Whether stray rates come from `fixed_stray_rate` rather
  than being estimated. Default `1`.

- fixed_stray_rate:

  Array `[n_pop x n_years]` of stray rates in \\\[0, 1\]\\. Default `0`.

- use_stray_rate_prior:

  Integer (0/1). Beta priors on estimated stray rates. Only relevant
  when `use_fixed_stray_rate = 0` and `n_pop > 1`; an error is raised
  alongside `use_fixed_stray_rate = 1`. Default `0`.

- stray_rate_prior:

  Data frame with columns `pop`, `block`, `mu` in \\(0,1)\\ and `sd`,
  one row per population and block. Default `NULL`.

- spawn_seas:

  Integer season index in which spawning occurs. Default `1`.

- sexratio_spec:

  Estimation structure for `sexratio_pars`
  `[n_pop x n_regions x n_blocks]`. Default `"fix"`, which is required
  when `n_sexes = 1`. See
  [`do_sexratio_pars_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_sexratio_pars_mapping.md).

- sexratio_blocks:

  Character vector, one entry per population and region, either
  `"none_Pop_x_Region_x"` or `"Block_k_Year_a-b_Pop_x_Region_x"`, with
  `"terminal"` allowed as the end year. Default one block each.

- use_rinit:

  Integer (0/1). Whether `ln_rinit` initializes the population
  separately from `ln_global_R0`. Under `0` (default) `ln_rinit` is
  fixed and `ln_global_R0` does both jobs.

- init_age_devs_shared:

  Integer vector of length `n_ages - 1` giving the factor level of each
  age position for `ln_InitDevs`; positions sharing a value share one
  parameter. Read under `equil_init_age_strc = 3`, and also respected by
  `InitDevs_spec = "est_shared_r"` (applied per population, with a
  population offset) and `"est_shared_pop_r"` (applied globally, no
  offset). `c(1:42, rep(42, 9))` gives 42 free parameters for a 52-age
  model with 43 data ages. Default `NULL`.

- R0_blocks:

  Character vector of time blocks for `R0`, one per population, in the
  selectivity block vocabulary: `"none_Pop_<p>"` or
  `"Block_<b>_Year_<a>-<e>_Pop_<p>"` on 1-based year indices, with
  `"terminal"` allowed. Under `"mean_rec"` a block is a productivity
  regime; under a stock-recruit form it makes the curve time-varying, so
  `S0`, depletion and any reference point step at the boundary. Default
  `NULL`.

- R0_ref_block:

  Integer, the block whose `R0` is used wherever a single value is
  needed: the initial age structure, the regional apportionment, the
  `R0` prior, the `ln_rinit` penalty, and the stock-recruit scale under
  `sr_R0_spec = "shared"`. Default 1.

- use_r0_prior:

  Integer (0/1). Lognormal prior on `R0`. Default 0.

- r0_prior:

  Data frame with columns `pop`, `mu` on the natural scale and `sd` on
  the log scale. Required when `use_r0_prior = 1`.

- Use_rinit_pen:

  Integer (0/1). Whether \\\log(R\_{init} / R_0)\\ is penalized under
  `use_rinit = 1`. An equilibrium recruitment stands for an average of
  several years, so \\\sigma_R / (1 / M - 0.5)\\ is a reasonable sd.
  Default 0.

- rinit_pen_sd:

  Standard deviation of that penalty, log scale. Default 1.

- ...:

  Optional named starting values: `ln_global_R0` `[n_pop]`, `ln_rinit`
  `[n_pop]`, `rec_region_prop_pars` `[n_pop x (n_regions - 1)]`,
  `rec_seas_prop_pars` `[n_pop x (n_seas - 1)]`, `steepness_h`
  `[n_pop x n_regions]` on the bounded logit scale, `ln_InitDevs`
  `[n_pop x n_regions x (n_ages - 1) x n_sexes]` (a 3-D array is
  expanded across sexes), `ln_RecDevs` `[n_pop x n_regions x n_years]`,
  `ln_sigmaR` `[2 x n_pop x n_regions]`, `sexratio_pars`
  `[n_pop x n_regions x n_blocks]`.

- ln_global_R0_spec:

  `"est"` (default) or `"fix"`. `"fix"` maps `ln_global_R0` off at its
  starting value, so the deviations hold log recruitment outright. Under
  `rec_model = "mean_rec"` a random walk with
  `dont_pen_recdev_first >= 1` leaves the level unidentified and is
  refused; a walk with the first year still penalized is accepted with a
  warning.

## Value

`input_list` with recruitment fields set in `$data` and `$par`, and maps
built in `$map` for `rec_region_prop_pars`, `rec_seas_prop_pars`,
`ln_sigmaR`, `ln_InitDevs`, `ln_RecDevs`, `RecDevs_rho`, `steepness_h`,
`sexratio_pars` and `stray_rate_pars`. Character codes for
`init_age_strc` and `equil_init_age_strc` are converted to integers
before storage.

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
