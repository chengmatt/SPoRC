# version 1.2.0.9000 (dev-popn-seasons)
- Incorporates population-specific, seasonal, and discarding dynamics.

## Major changes

### Population structure and dynamics

* Added population-specific (natal homing) and seasonal dynamics to both simulation and estimation. Most dimensions now follow population, region, year, season, age, sex, and fleet.
* Added movement-after-mortality and continuous movement dynamics.
* Added continuous-time Markov chain (CTMC) movement estimation using preference functions.
* Added `move_expm_nsub` to `Setup_Mod_Movement`, allowing CTMC matrix exponentials to be evaluated using implicit backward Euler substeps, providing a faster alternative for differentiation while preserving non-negative fractions and exact survivor catch balance.
* Added support for local density dependence in reference points for meta-population and natal homing models.
* Added uncertainty estimation for reference points.

### Data and observation processes

* Added population-specific catches, indices, compositions, tagging data, discards, and discarded compositions to simulation and estimation.
* Recoded tagging to allow fleet-specific reporting rates and missing attributes in tagged fish. `conv_tag_t_tagging`, `ln_init_conv_tag_mort`, and `ln_conv_tag_shed` are now per-release-event vectors rather than global scalars; scalars remain supported and are recycled. `init_conv_tag_mort_spec`/`conv_tag_shed_spec` accept `"fix"`, `"est_shared"`, or `"est_all"`.
* Added internal OSA residuals for catch, indices, compositions, and tagging data.
* Added support for at-age data streams (e.g., ICES-style stock assessments).
* Survey fleets can observe recruitment deviations directly via `srv_idx_type = "recdev"`, reported as `RecDev_anom`.

### Recruitment and fishing mortality

* Added age-0 (`rec_lag = 0`) Beverton-Holt recruitment, driven by the same year's SSB.
* Added `sr_penalty` (`"none"`, `"bh"`, `"ricker"`) to fit a stock-recruit relationship as a penalty rather than as the recruitment process itself. `sr_R0_spec` controls the scale of the stock-recruit curve. MSY reference points remain unavailable under mean recruitment.
* Added `Fdev_model` (`"iid"`, `"rw"`, `"ar1"`) to `Setup_Mod_Catch_and_F`, with `Fdev_rho`/`Fdev_rho_spec` for the AR1 case.
* Fishing mortality deviations now distinguish missing catch observations from true recorded zeros: `NA` in `ObsCatch` with `UseCatch == 0` estimates F as an active year, whereas an observed zero forces `Fmort` to zero.
* Added the ability to parameterize fishing mortality as free annual log-F via `ln_F_mean_spec = "fix"`.

### Biological processes and growth

* Growth parameters can vary over time via `growth_tv_model` (`"iid"`/`"rw"`), with support for year- and cohort-specific variation.
* Added semi-parametric growth (`growth_semipar`: `"iid"`, `"rw"`, `"2dar1"`, `"3dmarg"`, `"3dcond"`) for year-by-age deviations in mean length at age.
* Added the Richards growth form (`growth_model = "richards"`), with von Bertalanffy growth recovered when `rho = 1`.
* Added flexible size-age and weight-at-age inputs, including fleet-specific size-age keys and selection-weighted weight at age.
* Added conditional age-at-length fitting and Francis reweighting.

### Selectivity

* Added non-parametric, bicubic, logistic-with-asymptote, and fixed selectivity forms.
* Added sex-specific selectivity offsets, allowing selectivity parameters to be linked across sexes using `"par"`, `"scale"`, `"par_scale"`, `"apical"`, and `"par_apical"` forms.
* Added time-, bin-, and surface-based selectivity deviations and expanded selectivity smoothness penalties.
* Corrected the double-normal (`"dbnrml"`) selectivity parameterization and added flexible anchoring of its ascending limb.
* Added `_NSelBins_<n>` suffixes to selectivity model strings to hold bins beyond `n` at the value of bin `n`.
* Added `fish_sel_norm_bins`/`srv_sel_norm_bins` to specify the bins over which non-parametric log-scale selectivity is standardized.

### Reference points

* Reference points now support local density dependence for meta-population and natal homing models.
* Added uncertainty estimation for reference points.
* `Get_Reference_Points` now uses growth-derived weight at age from the model report rather than the data placeholder.

### Case studies and applications

* Added BSAI Atka mackerel, West Coast sablefish, EBS Pacific cod, North Sea sand eel, GOA Rex Sole, BSAI Northern Rock Sole case studies bridging recent assessments.

### Model fitting and diagnostics

* Improved Newton refinement in `fit_model` by obtaining the Hessian directly from the AD tape rather than finite differencing.
* Added 95% confidence intervals for SDNR based on a chi-squared test for OSA residuals.
* Added OSA residuals for time-series observations, compositions, and tagging data.

## Improvements
- Numeric and character codes are both accepted for simulator dynamics, matching the estimation model.
- Added a shared process-error sharing-spec helper (`build_pe_map`/`build_shared_spec_map`); migrated `sigmaF_spec`, `sigmaC_spec`, `sigmaR_spec`, `Fdev_rho_spec`, and `do_Fmort_mapping` onto it (for developers).
- Refactored movement's continuous process-error map and log-likelihood onto the same shared machinery (for developers).
- Newton refinement in `fit_model` now takes its Hessian from the AD tape (`obj$he`) instead of finite-differencing, and stops early on a non-finite Hessian.

## Bug Fixes
- Fixed the time-varying selectivity smoothing penalty not being applied in the first year.

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
