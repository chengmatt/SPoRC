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
- Random walk and AR1 fishing mortality deviations no longer require
  catch-active years to be contiguous. A fishery closure spanning
  multiple years is now handled via the exact closed-form marginal
  transition over the elapsed gap `d` (random walk: variance inflates to
  `d*sigma^2`; AR1: mean decays by `rho^d`, variance becomes
  `sigma^2 * sum_{i=0}^{d-1} rho^(2i)`) – mathematically identical to
  estimating deviations for the closed years and integrating them out,
  without actually estimating them. Both reduce exactly to the standard
  single-step transition when `d = 1`.
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

### Minor changes

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
  observations (indices and catch).
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

### Bug Fixes

- Fixed bug on smoothing penalty of time-varying selectivity where it
  was not being applied in the first year.
- Fixed inconsistent logic between the fishing mortality and discard
  mortality rate zeroing conditions for population-specific catch:
  `Fmort` was forced to zero if *any* population’s catch was unused,
  while the parallel discard mortality rate condition correctly required
  *all* populations to be unused. Both now consistently use the latter.
- Fixed a data issue in the bundled `sgl_rg_sable_data` dataset where 3
  cells of `ObsCatch` (fleet 2, years 1-3, representing a fleet that had
  not yet started fishing) were `NA`; corrected to `0` to represent a
  true closure under the new missing-vs-true-zero convention above
  (prior fitted results for this dataset are unaffected).

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
