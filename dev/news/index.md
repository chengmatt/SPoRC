# Changelog

## version 1.2.0.9000 (dev-popn-seasons)

- Incorporates population-specific, seasonal, and discarding dynamics.

### Major changes

- Added population-specific (natal homing) and seasonal dynamics, in
  both simulation and estimation. Most dimensions now follow population,
  region, year, season, age, sex, fleet.
- Recoded tagging to allow fleet-specific reporting rates and missing
  attributes in tagged fish. `conv_tag_t_tagging`,
  `ln_init_conv_tag_mort`, and `ln_conv_tag_shed` are now
  per-release-event vectors rather than global scalars (scalars still
  work, recycled). `init_conv_tag_mort_spec`/`conv_tag_shed_spec` take
  `"fix"`, `"est_shared"`, or `"est_all"`.
- Reference points now support local density-dependence for
  meta-population and natal homing models.
- Can simulate and fit population-specific catches, indices,
  compositions, tagging data, discards, and discarded compositions.
- Can build internal OSA residuals for catch, indices, compositions, and
  tagging data.
- Added `Fdev_model` (`"iid"`, `"rw"`, `"ar1"`) to
  `Setup_Mod_Catch_and_F` for fishing mortality deviations, with
  `Fdev_rho`/`Fdev_rho_spec` for the AR1 case.
- Fishing mortality deviations now distinguish a missing catch
  observation from a true recorded zero: `NA` in `ObsCatch` with
  `UseCatch == 0` estimates F as an active year; an observed zero forces
  `Fmort` to zero.
- Added age-0 (`rec_lag = 0`) Beverton-Holt recruitment via `rec_lag` in
  `Setup_Mod_Rec`/`Setup_Sim_Rec`, driven by the same year’s own SSB.
- Added movement-after-mortality and continuous movement dynamics.
- Added `sr_penalty` (`"none"`, `"bh"`, `"ricker"`) in `Setup_Mod_Rec`
  to fit a stock-recruit curve as a penalty rather than the recruitment
  process itself, valid only with `rec_model = "mean_rec"`. `sr_R0_spec`
  sets the curve’s scale (`"shared"`, `"est"`, `"rinit"`). MSY reference
  points remain unavailable under mean recruitment.
- Added `fish_sel_norm_bins`/`srv_sel_norm_bins` to name which bins
  non-parametric log-scale selectivity (`"nonparlog"`) standardizes
  over.
- Added a BSAI Atka mackerel case study bridging the 2024 AMAK
  assessment, with `sgl_rg_bsai_atka_data`.
- Added a West Coast sablefish case study bridging the 2025 assessment,
  with `sgl_rg_wc_sablefish_data`.
- Selectivity sexes can be linked through offsets:
  `fish_sel_sex_offset`/`ret_sel_sex_offset`/`srv_sel_sex_offset`, with
  `"par"`, `"scale"`, `"par_scale"`, `"apical"`, and `"par_apical"`
  (double normal only) forms.
- Added an optional `_NSelBins_<n>` suffix on any selectivity model
  string, holding bins beyond `n` at bin `n`’s value (the plateau
  convention the bicubic form already had).
- Initial age deviations gained a sex dimension (`ln_InitDevs` is now
  `[n_pop, n_regions, n_ages - 1, n_sexes]`), via `InitDevs_sex_spec` in
  `Setup_Mod_Rec`/`Setup_Sim_Rec` (`"est_shared_s"` default,
  bit-identical to before; `"est_all"` gives each sex its own curve).
- Sex-specific initial age deviations can be penalized against each
  other via `Use_init_sex_pen`/`init_sex_pen_sigma` in `Setup_Mod_Rec`
  (requires `InitDevs_sex_spec = "est_all"`).
- Corrected the double normal (`"dbnrml"`) selectivity form: the
  plateau-start parameter is now on the bin scale, and both limbs are
  rescaled to hit `p5`/`p6` at the first/last bin exactly. Existing
  `"dbnrml"` fits should be re-run.
- Growth starting values now go through `starting_values` in
  `Setup_Mod_Biologicals`, like every other parameter: `growth_pars` is
  now `ln_growth_pars`, and
  `growth_tv_sigma`/`growth_semipar_sigma`/`growth_semipar_rho` are now
  one array, `growth_pe_pars`.
- `SizeAgeTrans_fish`/`SizeAgeTrans_srv` in `Setup_Mod_Biologicals` let
  a fleet read its own fixed size-age key under `growth_model = "none"`,
  mirroring `WAA_fish`/`WAA_srv`.
- `get_osa`’s external path now reads `N`/`DM_theta`/`LN_Sigma`
  consistently across every `comp_type`, and fills in `bins`/`bin_label`
  from the data when omitted.
- `plot_resids` facets conditional age-at-length residuals by length
  bin, so they no longer stack onto a single point.
- Francis reweighting now covers conditional age-at-length, via
  `new_fish_caal_wts`/`new_srv_caal_wts` from `do_francis_reweighting`.
- Added `get_caal_fits`, returning conditional age-at-length fits as
  observed/predicted mean age within each length bin.
- `Get_Reference_Points` now reads the growth-derived weight at age from
  the report instead of the data placeholder.
- Retrospective peels now truncate the growth deviations, their maps,
  and the conditional age-at-length arrays, which were previously left
  at full length and read incorrectly.
- Added `srv_waa_selected` in `Setup_Mod_SrvIdx_and_Comps`, the survey
  twin of `fish_waa_selected`.
- `growth_A2` in `Setup_Mod_Biologicals` accepts `"Linf"` in place of a
  second reference age.
- Growth parameters can vary over time via `growth_tv_model`
  (`"iid"`/`"rw"`, with `growth_tv_link`, `growth_tv_years`), reported
  as `growth_pars_y`.
- Time-varying growth can be carried cohort by cohort via
  `growth_tv_type = "cohort"`, evaluated inside the population dynamics
  year loop.
- Growth gained the Richards form (`growth_model = "richards"`), a sixth
  parameter `rho`; `rho = 1` is von Bertalanffy.
- Added semi-parametric growth (`growth_semipar`: `"iid"`, `"rw"`,
  `"2dar1"`, `"3dmarg"`, `"3dcond"`), a year-by-age deviation surface on
  mean length at age.
- Length compositions can be recorded on coarser bins than the
  population via `LenBinMap` in `Setup_Mod_Biologicals`.
- Length selectivity can be applied at length rather than at age when
  forming length compositions, via `FishLenComps_sel`/`SrvLenComps_sel`.
- The catch in biomass can use the selection-weighted weight at age via
  `fish_waa_selected`.
- The double normal’s ascending limb can be anchored at a bin other than
  the first, via `fish_sel_dbnrml_startbin`/`srv_sel_dbnrml_startbin`.
- The initial equilibrium recruitment’s offset from the recruitment
  level can be penalized via `Use_rinit_pen`/`rinit_pen_sd` in
  `Setup_Mod_Rec`.
- Added an EBS Pacific cod case study bridging the 2024 assessment,
  [`vignette("ae_ebs_pacific_cod_case_study")`](https://chengmatt.github.io/SPoRC/dev/articles/ae_ebs_pacific_cod_case_study.md),
  with `sgl_rg_ebs_pcod_data`.
- Fixed the 3D GMRF process error’s precision matrix, previously built
  over the full bin dimension rather than the bins actually evaluated.
- Survey fleets can observe recruitment deviations directly via
  `srv_idx_type = "recdev"`, reported as `RecDev_anom`.
- A conditioned operating model now records `ln_RecDevs` correctly
  instead of leaving it at zero.
- The double normal’s default starting values now place its peak
  mid-range instead of at zero, which previously divided by zero.

### Minor changes

- `likelihood_profile`/`get_nLL_plot` now carry `Init_Sex_nLL`,
  `Rec_level_nLL`, and `SR_pen_nLL`.
- Indices can be restricted to a subset of ages via
  `srv_idx_ages`/`fish_idx_ages`.
- Index catchability can be concentrated out of the likelihood via
  `srv_q_type` (`"est"`/`"arith"`/`"geo"`).
- Index likelihoods are now selectable per fleet via
  `SrvIdx_LikeType`/`FishIdx_LikeType` (`"lognormal"`, `"normal"`,
  `"mvn"`).
- The simulator now draws index observations under the fleet’s chosen
  error structure instead of always lognormal.
- Population-specific index likelihoods now follow the fleet’s
  `LikeType` for `"lognormal"`/`"normal"` instead of always lognormal.
- Deviation penalties can be centred on their own mean via
  `RecDevs_pen_center`/`InitDevs_pen_center`/`Fdev_pen_center`
  (`"own_mean"`).
- Fishing mortality can be parameterized as free annual log-F via
  `ln_F_mean_spec = "fix"`.
- The steepness prior’s beta support is configurable via `lb`/`ub`
  columns on `h_prior`.
- Added `init_age_strc = 4` (“free”), skipping the equilibrium age
  structure projection entirely.
- Added selectivity form `"nonparlog"`, non-parametric on the log scale,
  standardized to average one within a year.
- Added bin-override selectivity deviations via
  `*_sel_bin_dev_bins`/`cont_tv_*sel_bin_devs`.
- `Wt_Rec` now applies inside the sum, accepting a per-deviation weight
  array; the initial age penalty moved to its own `Wt_Init_Rec`.
- Added a centering penalty on selectivity fixed-effect parameter sets
  via `Use_fish_selex_penalty`/`fish_selex_penalty` and its
  retention/survey equivalents.
- Selectivity prior tables gained an optional `type` column (`"par"` or
  `"value"`, the latter constraining a realized selectivity value rather
  than a parameter).
- Selectivity smoothness penalties are now specified per fleet and per
  term rather than shared across all fleets.
- Fishery indices can be observed part-way through a season via
  `t_fish`, the same convention `t_srv` already used.
- Recruitment deviations that are mapped off are no longer penalized.
- Initialization fishing mortality is now its own parameter,
  `init_F_par`, with `init_F_form` (`"prop"`/`"abs"`) and `init_F_spec`.
- Added `backwards_compatibility_guard()` to fill in fields older input
  lists predate.
- Renamed `ln_srv_fixed_sel_pars`/`ln_fish_fixed_sel_pars` to
  `srv_fixed_sel_pars`/`fish_fixed_sel_pars`.
- Added non-parametric selectivity, bicubic selectivity, logistic
  selectivity with an asymptote, and fixed selectivity inputs.
- `init_F_prop` is now region/season/fleet-specific rather than tied to
  the first fishery fleet.
- `do_retrospective` can return the retrospective RTMB model objects via
  `return_models`.
- Added 95% confidence intervals for SDNR based on a chi-squared test
  for OSA residuals.
- Non-parametric selectivity is now forced to be mean-standardized.
- Added OSA residuals for time-series observations (indices, catch) as
  well as compositions and tagging data.
- Consolidated selectivity smoothness penalties into six weights
  (`smooth_bin_curve`, `smooth_bin_diff`, `smooth_yr_diff`,
  `smooth_yr_curve`, `smooth_dome`, `smooth_mean_center`) via
  `fish_sel_pen_wts`/`ret_sel_pen_wts`/`srv_sel_pen_wts`, evaluated on
  the realized selectivity surface for any form. Removed the old
  `cont_tv_*_sel_penalty` flags and the legacy `bin_curve`/`yr_curve`
  terms; all six weights default to off.

### Improvements

- Numeric and character codes are both accepted for simulator dynamics,
  matching the estimation model.
- Added a shared process-error sharing-spec helper
  (`build_pe_map`/`build_shared_spec_map`); migrated `sigmaF_spec`,
  `sigmaC_spec`, `sigmaR_spec`, `Fdev_rho_spec`, and `do_Fmort_mapping`
  onto it (for developers).
- Refactored movement’s continuous process-error map and log-likelihood
  onto the same shared machinery (for developers).
- Newton refinement in `fit_model` now takes its Hessian from the AD
  tape (`obj$he`) instead of finite-differencing, and stops early on a
  non-finite Hessian.

### Bug Fixes

- Fixed the time-varying selectivity smoothing penalty not being applied
  in the first year.

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
