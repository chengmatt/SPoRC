# Changelog

## version 1.2.0.9000 (dev-popn-seasons)

- Incorporates population-specific, seasonal, and discarding dynamics.

### Major changes

- Incorporated ability to simulate and estimate both population-specific
  (natal homing) and seasonal dynamics.
- Most model and simulation dimensions now include population- and
  season-specific indices, following the general dimension order:
  population, region, year, season, age, sex, fleet.
- Recoded tagging module to allow fleet-specific tag reporting rates, as
  well as missing attributes in tagged fish.
- `conv_tag_t_tagging`, `ln_init_conv_tag_mort`, and `ln_conv_tag_shed`
  (in both `Setup_Sim_Tagging`/`Simulate_Population` and
  `Setup_Mod_Tagging`/the RTMB model) are now per-tag-release-event
  vectors (length `n_conv_tag_cohorts`/`n_tag_rel_events`) rather than
  global scalars, so time-of-tagging, tag-induced mortality, and chronic
  shedding can differ across release cohorts. Scalars are still accepted
  and recycled to all release events, so existing calls are unaffected.
  `init_conv_tag_mort_spec`/`conv_tag_shed_spec` in `Setup_Mod_Tagging`
  now accept `"fix"`, `"est_shared"` (one value estimated across all
  events, replacing the old `"est"`), or `"est_all"` (an independent
  value estimated per release event).
- Reference points for local density-dependence in both meta-population
  and natal homing contexts.
- Added ability to simulate and fit to population-specific catches,
  indices, compositions, and tagging data.
- Added ability to simulate and fit to population-specific discards and
  discarded compositions.
- Added in ability to internally construct OSA residuals for catch,
  indices, compositions, and tagging data.
- Added `Fdev_model` option (`"iid"`, `"rw"`, `"ar1"`) to
  `Setup_Mod_Catch_and_F` for fishing mortality deviations
  (`ln_F_devs`), with a new AR1 correlation parameter `Fdev_rho`
  (shared/fixed via `Fdev_rho_spec`). Previously only IID deviations
  were supported.
- Fishing mortality deviations now distinguish a genuinely missing
  aggregate catch observation from a true recorded zero. If `ObsCatch`
  is `NA` at a cell where `UseCatch == 0` (and no population-specific
  catch is used), fishing is assumed to have continued and
  `Fmort`/`ln_F_devs` are estimated as an ordinary active year; if
  `ObsCatch` holds a real value (typically `0`), the cell is treated as
  a genuine closure (`Fmort` forced to zero, no deviation estimated),
  matching prior behavior. See `@param ObsCatch` in
  `Setup_Mod_Catch_and_F` and `Get_Fdev_PE_loglik`.
- Added age-0 (`rec_lag = 0`) Beverton-Holt recruitment, set via
  `rec_lag` in `Setup_Mod_Rec`/`Setup_Sim_Rec`. Previously `rec_lag` had
  to be `>= 1` (recruitment driven by SSB from `rec_lag` seasons prior,
  entering in any season). With `rec_lag = 0`, recruitment for a year is
  driven by that *same* year’s own SSB.
- Incorporated additional movement ordering options including: movement
  after mortality and continuous movement dynamics.
- Added the ability to fit a stock-recruit curve as a penalty rather
  than as the recruitment process, via `sr_penalty` (`"none"`, `"bh"`,
  `"ricker"`) in `Setup_Mod_Rec`, valid only with
  `rec_model = "mean_rec"`. Recruitment stays a mean with deviations and
  the curve is evaluated alongside the dynamics without advancing them,
  with a normal density placed on the log residual over `sr_pen_yrs` at
  `sr_pen_sigma`. This is the arrangement several ADMB templates use,
  where a mean recruitment parameter generates recruitment and a
  separate unfished recruitment parameter carries the curve. The curve’s
  scale is set by `sr_R0_spec`: `"shared"` reuses `ln_global_R0`,
  `"est"` estimates its own `ln_sr_R0`, and `"rinit"` takes it from
  `ln_rinit` so one parameter sets both the unfished age structure and
  the curve (requires `use_rinit = 1`). `SR_pred`, `SR_pen_nLL` and
  `sr_R0` are reported. MSY reference points remain unavailable under
  mean recruitment, since the curve does not govern the stock;
  `Get_Reference_Points` names `sr_R0` and `h_trans` for callers who
  want to supply that curve through `srr_opt` deliberately.
- Added `fish_sel_norm_bins` and `srv_sel_norm_bins` to
  `Setup_Mod_Fishsel_and_Q` and `Setup_Mod_Srvsel_and_Q`, naming the
  bins the mean-one standardization averages over for non-parametric
  log-scale selectivity (`"nonparlog"`). The default is every bin. A
  gear whose catchability is defined against part of the bin range
  standardizes over that part; the difference between windows is a
  constant that catchability absorbs when it is free, but not when it
  carries an informative prior.
- Added a BSAI Atka mackerel case study bridging the 2024 assessment
  (AMAK Model 16.0b), with `sgl_rg_bsai_atka_data`.

### Minor changes

- Indices can now be restricted to a subset of ages, via `srv_idx_ages`
  in `Setup_Mod_SrvIdx_and_Comps` and `fish_idx_ages` in
  `Setup_Mod_FishIdx_and_Comps`. Either a list with one element per
  fleet (a vector of ages, or `NULL` for all ages) or an array
  `[n_ages x n_fleets]` of 0/1 weights. The restriction applies to the
  index sum only, so the fleet’s compositions still use the full age
  range and its selectivity is untouched; a fleet restricted to a single
  age becomes an index of that age alone, which is how an age-1 index
  could be specified. Existing input lists sum over every age, which
  `NULL` (the default) reproduces.
- Index catchability can now be concentrated out of the likelihood
  rather than estimated, via `srv_q_type` in `Setup_Mod_Srvsel_and_Q`.
  `"est"` (default) estimates `ln_srv_q` as before; `"arith"` solves it
  as the ratio of mean observed to mean predicted index; `"geo"` solves
  it on the log scale as `exp(mean(log(obs) - log(pred)))`. Both
  analytic forms solve one catchability per region and fleet using only
  the years with observations, ignore block structure, and fix that
  fleet’s `ln_srv_q` regardless of `srv_q_spec`.
- Index likelihoods are now selectable per fleet, via
  `SrvIdx_LikeType`/`FishIdx_LikeType`. Options are `"lognormal"`
  (default, the previous and only behaviour), `"normal"` on the
  arithmetic scale, and `"mvn"`, a multivariate normal on the arithmetic
  scale using a fixed covariance matrix supplied through
  `SrvIdx_Cov`/`FishIdx_Cov`. This is appropraite when a covariance is
  availiable across years and a diagonal likelihood would overweight the
  series. The covariance is inverted and its log-determinant taken once
  at setup, so only the quadratic form is on the AD tape.
- The simulator now draws index observations under the fleet’s chosen
  error structure, via `FishIdx_LikeType`/`FishIdx_Cov`/`UseFishIdx` in
  `Setup_Sim_Fishing` and `SrvIdx_LikeType`/`SrvIdx_Cov`/`UseSrvIdx` in
  `Setup_Sim_Survey`, which `simulation_self_test` fills from the input
  list automatically. Previously every index stream was drawn lognormal
  regardless of what the estimation model would fit, so a self-test of a
  normal or mvn index generated data under the wrong error structure and
  reported the mismatch as estimation bias. An mvn fleet takes its scale
  from the covariance rather than the SE array (the two are not
  interchangeable; the pollock trawl covariance’s diagonal is about
  twice the reported SEs) through a one-factor decomposition
  (`cov_to_factor`): one shared draw per fleet and replicate scales the
  whole series, with per-observation loadings resolved by each
  observation’s position among the fleet’s fitted cells, and cells
  outside the covariance (projection years in closed loop) reuse the
  mean loading and scale. The approximation suits a survey-wide scaling
  error, not a banded or autoregressive covariance.
- Population-specific index likelihoods now follow the fleet’s
  `FishIdx_LikeType`/`SrvIdx_LikeType` for `"lognormal"` and `"normal"`;
  previously they were silently lognormal whatever the fleet’s choice.
  Under `"mvn"` the population stream stays lognormal, since the
  covariance describes the regional series only. The simulator draws its
  population-specific streams under the same rules.
- Deviation penalties can now be centred on the deviations’ own mean
  rather than on an asserted prior mean, via `RecDevs_pen_center` and
  `InitDevs_pen_center` in `Setup_Mod_Rec` and `Fdev_pen_center` in
  `Setup_Mod_Catch_and_F`. `"fixed"` (default, previous behaviour)
  centres on zero or the bias-corrected `-sigmaR^2/2` and constrains
  both the level and the spread of the deviations; `"own_mean"` centres
  on the mean of the estimated deviations so only the spread is
  penalized. The second is what a sum of squares about a series’ own
  mean amounts to. Because the level is then unpenalized it must be
  pinned elsewhere, by a prior on `R0` or by fixing a deviation, or the
  likelihood is flat along it; for fishing mortality under a
  mean-plus-deviations parameterization the level is already carried by
  `ln_F_mean`, so `"own_mean"` is what avoids penalizing it twice.
  `RecDevs_pen_center = "own_mean"` cannot be combined with
  `do_rec_bias_ramp = 1`, whose offset is meaningless once the mean is
  estimated rather than asserted.
- Fishing mortality can now be parameterized as free annual log-F, via
  `ln_F_mean_spec = "fix"` in `Setup_Mod_Catch_and_F`: `ln_F_mean` is
  mapped off at its starting value (defaulting to 0 under this spec) so
  the deviations carry all of log F. It must be paired with
  `Fdev_pen_center = "own_mean"`, `Fdev_model = "rw"`, or
  `Use_F_pen = 0`; an `"iid"` or `"ar1"` penalty centred on the fixed
  zero mean would shrink the deviations toward F = 1, so that
  combination errors at setup. `"est"` (default) keeps the
  mean-plus-deviations form unchanged. Setup now also warns when
  `Fdev_pen_center = "own_mean"` is combined with an estimated
  `ln_F_mean` in configurations where nothing else reads the mean
  (absolute-rate initialization F, or a free initial age structure),
  since the mean and the deviations’ level then trade off along an
  exactly flat ridge.
- The steepness prior’s beta support can be configured through optional
  `lb` and `ub` columns on `h_prior`, defaulting to the previous
  hard-coded `0.2` and `1`.
- Added `init_age_strc = 4` (`"free"`), which projects no equilibrium
  age structure at all: numbers at age 2 and older are
  `exp(ln_InitDevs)` outright, apportioned by sex ratio, with age 1
  still taken from recruitment. The four existing options all project an
  equilibrium forward from `R0` and treat the deviations as
  multiplicative departures from it, which pulls the initial age
  structure toward `R0` and leaves the two partly confounded. Use `4`
  when the initial age structure carries no information about `R0`.
  Per-age sharing and fixing of the deviations is unchanged and still
  goes through `equil_init_age_strc` and `init_age_devs_shared`. Note
  that under `4` the deviations are log abundances rather than log
  ratios, so pair it with `equil_init_age_strc = 0` if no prior on that
  scale is wanted.
- Added selectivity form `"nonparlog"`, non-parametric on the log scale
  and standardized so each year’s selectivity averages to one across
  bins. It differs from the existing `"nonpar"` in both respects: that
  one bounds every raw value below one through `plogis` and standardizes
  over years and bins jointly, whereas this one leaves the scale free
  and centers within the year. Only the differences among the parameters
  within a year are identified, so their level is absorbed by
  catchability or fishing mortality. It reuses the existing
  `*_sel_nonpar_est_bins` bin-sharing, so a flat tail beyond the last
  estimated age needs no extra machinery.
- Added bin-override selectivity deviations, via `*_sel_bin_dev_bins`
  and `cont_tv_*sel_bin_devs` in `Setup_Mod_Fishsel_and_Q`,
  `Setup_Mod_Retsel`, and `Setup_Mod_Srvsel_and_Q`. Bins named for a
  fleet take a freely estimated annual value, `exp(ln_*sel_bin_devs)`,
  instead of whatever the fleet’s functional form produces, while the
  rest of the curve keeps its parametric shape. This covers a gear whose
  curve is logistic over most of its range but whose youngest bin is
  governed by availability rather than by the gear. The override is
  applied last, so it also overrides a form that standardizes
  internally. The deviations carry their own process error (`"none"`,
  `"iid"`, or `"rw"`) through the same `Get_sel_PE_loglik` the other
  deviations use, and live in their own parameter array rather than
  sharing `ln_*sel_devs`, whose third dimension already carries the
  functional form’s parameter deviations under `"iid"`/`"rw"` and would
  collide. Only the named bins are estimated; a fleet naming none adds
  no free parameters, which is what older input lists get.
- `Wt_Rec` is now applied inside the sum rather than to the whole
  recruitment penalty, so it accepts a numeric array
  `[n_pop x n_regions x n_est_rec_devs]` for deviation-specific
  weighting as every other `Wt_*` already did. A weight of zero leaves a
  deviation estimated but removes it from the penalty, which is how a
  stock-recruit relationship is fit over a window of years while
  recruitment stays free everywhere; this is distinct from
  `dont_est_recdev_last`, which deletes the deviations so those years
  revert to the deterministic prediction. The initial age penalty is
  dimensioned differently and so has moved to its own `Wt_Init_Rec`,
  which defaults to `Wt_Rec` when that is a scalar and must be supplied
  explicitly when it is an array. A scalar `Wt_Rec` is unchanged in
  every respect, since `sum(w * x)` and `w * sum(x)` agree; input lists
  built before this change get `Wt_Init_Rec` filled from `Wt_Rec` by
  [`maintain_backwards_compatibility()`](https://chengmatt.github.io/SPoRC/dev/reference/maintain_backwards_compatibility.md).
  `Setup_Mod_Rec` gained no new argument, and self-testing still
  requires a scalar to derive the operating model’s `ln_sigmaR`.
- Added a centering penalty on sets of selectivity fixed-effect
  parameters, via `Use_fish_selex_penalty`/`fish_selex_penalty` in
  `Setup_Mod_Fishsel_and_Q` and the retention and survey equivalents.
  Each row of the table penalizes `wt * (log(mean(exp(pars))))^2` over
  the set of parameters named in its `par` column, which may be a single
  index or a vector naming a whole set. This pins the scalar of a
  non-parametric selectivity curve that catchability or fishing
  mortality would otherwise absorb, and is a softer constraint than
  fixing a bin outright. Intended for parameter sets held on the log
  scale.
- Selectivity prior tables
  (`fish_selex_prior`/`ret_selex_prior`/`srv_selex_prior`) gained an
  optional `type` column, so a row may now constrain either a fixed
  selectivity parameter or a realized selectivity value. `"par"` (the
  default when the column is absent, and the previous and only
  behaviour) is the existing lognormal prior on one parameter, with `mu`
  on the natural scale and `sd` on the log scale. `"value"` is a normal
  prior on the realized selectivity value at one bin, with both
  hyperparameters on the natural scale: `par` instead names the bin on
  the grid the stream’s selectivity is parameterized on (ages or lengths
  per its selectivity type), and the value is read at the first model
  year of `block` (blocked and time-invariant selectivity are constant
  within a block). A `"value"` row constrains the derived selectivity
  value rather than the parameters, a penalty of the form
  `square(sel(age) - 1) / (2 * CV^2)` on a realized selectivity value
  being the motivating case: that statement’s Gauss-Newton Hessian in
  (log a50, log slope) is rank one, so no set of independent parameter
  priors can represent it, and bridging such a model through
  parameter-space priors alone changes the estimation problem it solves.
- Selectivity smoothness penalties are now specified per fleet and per
  term, rather than one specification shared by every fleet with one
  weight per term.
  `fish_sel_pen_wts`/`ret_sel_pen_wts`/`srv_sel_pen_wts` still accept a
  single named specification applied to all fleets, and now also accept
  an unnamed list with one specification per fleet. Within a
  specification, each weight may be a vector with one value per model
  year instead of a scalar, so a penalty can act only in the years
  selectivity is allowed to change, or act with a different strength in
  each year — a random walk with a year-specific standard deviation is
  `1 / (2 * sigma^2)` in that year with `normalize = FALSE`. A
  specification may also carry `bin_range`, the first and last bin the
  penalties act over, and `normalize`; both accept either one setting
  shared by all terms or a named list giving each term its own, so a
  shape penalty confined to the older ages can sit alongside a random
  walk spanning every age. Input lists built before this change hold a
  single shared named vector and are replicated across fleets by
  [`maintain_backwards_compatibility()`](https://chengmatt.github.io/SPoRC/dev/reference/maintain_backwards_compatibility.md),
  leaving their objective unchanged.
- Fishery indices can now be observed part-way through a season, via
  `t_fish` in `Setup_Mod_FishIdx_and_Comps` (estimation) and
  `Setup_Sim_Fishing` (simulation). `t_fish` is an array
  `[n_regions x n_seas x n_fish_fleets]` giving the fraction of the
  season elapsed at observation; numbers at age are decayed by
  `exp(-t_fish * ZAA)` before the index is formed, the same convention
  `t_srv` already used for surveys. Previously the fishery index was
  always formed from start-of-season numbers, which is what `t_fish = 0`
  (the default) reproduces, so existing input and simulation lists are
  unaffected — both `SPoRC_rtmb` and `Setup_sim_env` fill in zeros when
  `t_fish` is absent.
- Recruitment deviations that are mapped off are no longer penalized.
  `do_RecDevs_mapping` now mirrors the `ln_RecDevs` factor map into
  `data$map_ln_RecDevs`, `fit_model` refreshes it from the map handed to
  `RTMB::MakeADFun`, and `get_recruitment_penalty` evaluates the penalty
  only where it is not `NA` — the same convention `ln_F_devs` and
  `logit_dmr_devs` already followed. Previously the penalty ran over
  index ranges (`sigmaR_switch:n_est_rec_devs`) with no knowledge of the
  map, so a deviation fixed at zero still contributed `-dnorm(0, ...)`.
- Initialization fishing mortality is now its own parameter,
  `init_F_par` (`[n_regions x n_seas x n_fish_fleets]`), with two
  switches in `Setup_Mod_Rec`: `init_F_form` (`"prop"` = a proportion of
  the mean F via inverse-logit, bounded (0,1), so the initial age
  structure moves with mean F; `"abs"` = an absolute rate on the log
  scale, independent of mean F) and `init_F_spec` (`"fix"`/`"est"`,
  which sets only the mapping, so all four combinations are available).
  Previously `init_F_prop` held a fixed proportion as *data*, which
  meant the initialization F and the mean F were a single parameter;
  because catch constrains only the product of numbers-at-age and
  fishing mortality, that let the optimizer deplete the initial age
  structure and scale the F series together, fitting catch equally well
  with a smaller, harder-fished stock. `init_F_prop` is still accepted
  and is converted to the `"prop"` form, so existing calls are
  unaffected.
- Added `backwards_compatibility_guard()`, which fills in data and
  parameters that current `SPoRC_rtmb` calls expect but older input
  lists predate, so previously built objects keep evaluating unchanged.
- Changed parameter names of ln_srv_fixed_sel_pars and
  ln_fish_fixed_sel_pars to srv_fixed_sel_pars and fish_fixed_sel_pars
  for clarity.
- Included new options to estimate non-parametric selectivity, bicubic
  selectivity, logistic selectivity with an asymptote parameter, as well
  as provide fixed selectivity (fishery, retention, and survey) inputs.
- Changed dimensions of init_F_prop to be region, season, and
  fleet-specific (as opposed to just being based on the first fishery
  fleet).
- Coded in an argument in `do_retrospective` to return retrospective
  RTMB model objects (`return_models`).
- Added 95% confidence intervals for SDNR based on a Chi squared test
  for OSA residuals.
- Force non-parametric selectivity to be mean-standardized.
- Added OSA residuals and `oneStepPredict` functionality to time-series
  observations (indices and catch) as well as composition and tagging
  data. `get_osa` and `plot_resids` modified to accomodate plotting of
  OSA residuals for the aforementioned data sources. The internal-OSA
  `osa_method` is restricted to `"oneStepGeneric"`,
  `"oneStepGaussianOffMode"`, and `"oneStepGaussian"`; the `"cdf"`
  method is disallowed as it is numerically fragile for the discrete
  likelihoods used here and can silently return mis-calibrated
  residuals. Index-type OSA output (`get_osa(index_source = ...)`) now
  carries an `idx_type` discriminator column (values
  `"Catch"`/`"Discard"`/`"FishIdx"`/`"SrvIdx"`) instead of `comp_type`,
  since these sources are not compositions. Clarified in
  `get_idx_fits`’s documentation that its `resid` column is a raw
  log-scale (Pearson-style) residual, distinct from the one-step-ahead
  residuals produced by `get_osa`/`plot_resids`.
- Consolidated selectivity smoothness/regularization penalties into a
  single set of six weights (`smooth_bin_curve`, `smooth_bin_diff`,
  `smooth_yr_diff`, `smooth_yr_curve`, `smooth_dome`,
  `smooth_mean_center`) set via
  `fish_sel_pen_wts`/`ret_sel_pen_wts`/`srv_sel_pen_wts` in
  `Setup_Mod_Weighting`, evaluated directly on the realized selectivity
  surface so they apply uniformly to any selectivity functional form and
  fleet. Removed the old `cont_tv_*_sel_penalty` on/off flags and the
  legacy `bin_curve`/`yr_curve` terms; all six weights now default to
  `0` (off) unless explicitly set.

### Improvements

- Allowed both numeric and character codes for specifying dynamics in
  simulator, so as to be more consistent with how estimation models are
  specified, while maintaining backwards compatibility.
- Added a generic internal process-error sharing-spec helper
  (`build_pe_map`/`build_shared_spec_map`) and migrated `sigmaF_spec`,
  `sigmaC_spec`, `sigmaR_spec`, `Fdev_rho_spec`, and the fishing
  mortality deviation map (`do_Fmort_mapping`) onto it, replacing
  hand-enumerated per-combination branches with a single
  dimension-collapsing implementation (no change in behavior for
  existing models) (for developers).
- Refactored movement’s continuous process-error map
  (`do_cont_vary_move_mapping`) and log-likelihood
  (`Get_move_PE_loglik`) to use the same generic sharing-spec machinery
  and a single dimension-aware likelihood loop, replacing 10
  hand-written `iid_*` branches in each (no change in behavior for
  existing models; new dedicated tests added since this module
  previously had no coverage) (for developers).
- Newton refinement in `fit_model` now takes its Hessian from the AD
  tape (`obj$he`) instead of finite differencing the gradient with
  `optimHess`, which needed one gradient evaluation per parameter.
  Random-effects models continue to use `optimHess`, since RTMB does not
  implement a tape Hessian when random effects are present. Newton
  refinement now also stops early if the Hessian comes back non-finite,
  which can happen on models that have not converged, where second
  derivatives are undefined at parameter values the objective and
  gradient still evaluate at; previously a non-finite Hessian passed NaN
  through [`solve()`](https://rdrr.io/r/base/solve.html) without
  erroring and left `optim$par` and `optim$objective` as NaN.

### Bug Fixes

- Fixed bug on smoothing penalty of time-varying selectivity where it
  was not being applied in the first year.

## version 1.1.0

Release Date: 2026-3-31

### Major changes

- Added capability to estimate movement using continuous time Markov
  chains (CTMC) using preference functions.

### Bug Fixes

- Fixed bug to allow for size-based selectivity in the estimation model.

## version 1.0.0

Release Date: 2025-11-24

### Major changes

- First stable public release of SPoRC.
- Introduced a modular framework for estimating age-, sex-, season-, and
  region-structured population dynamics.
- Added closed-loop simulation and management strategy evaluation (MSE)
  capabilities.
