# version 1.2.0 (dev-popn-seasons)
- Incorporates population-specific, seasonal, and discarding dynamics.

## Major changes
- Incorporated ability to simulate and estimate both population-specific (natal homing) and seasonal dynamics.
- Most model and simulation dimensions now include population- and season-specific indices, following the general dimension order: population, region, year, season, age, sex, fleet.
- Recoded tagging module to allow fleet-specific tag reporting rates, as well as missing attributes in tagged fish. 
- Reference points for local density-dependence in both meta-population and natal homing contexts.
- Added ability to simulate and fit to population-specific catches, indices, compositions, and tagging data. 
- Added ability to simulate and fit to population-specific discards and discarded compositions. 

## Minor changes
- Changed parameter names of ln_srv_fixed_sel_pars and ln_fish_fixed_sel_pars to srv_fixed_sel_pars and fish_fixed_sel_pars for clarity.
- Included new options to estimate non-parametric selectivity, logistic selectivity with an asymptote parameter, as well as provide fixed selectivity (fishery, retention, and survey) inputs. 
- Changed dimensions of init_F_prop to be region, season, and fleet-specific (as opposed to just being based on the first fishery fleet).
- Coded in an argument in `do_retrospective` to return retrospective RTMB model objects (`return_models`).
- Added 95% confidence intervals for SDNR based on a Chi squared test for OSA residuals.
- Force non-parametric selectivity to be mean-standardized.

## Improvements
- Allowed both numeric and character codes for specifying dynamics in simulator, so as to be more consistent with how estimation models are specified, while maintaining backwards compatibility. 

## Bug Fixes
- Fixed bug on smoothing penalty of time-varying selectivity where it was not being applied in the first year. 

# version 1.1.0
Release Date: 2026-3-31

## Major changes
- Added capability to estimate movement using continuous time Markov chains (CTMC) using preference functions. 

## Bug Fixes
- Fixed bug to allow for size-based selectivity in the estimation model. 

# version 1.0.0
Release Date: 2025-11-24

## Major changes
- First stable public release of SPoRC.
- Introduced a modular framework for estimating age-, sex-, season-, and region-structured population dynamics.
- Added closed-loop simulation and management strategy evaluation (MSE) capabilities.
