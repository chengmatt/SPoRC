# version 1.2.0.9000 (dev-popn-seasons)
- Incorporates population-specific, seasonal, and discarding dynamics.

## Major changes

### Population structure and dynamics

* Added population-specific (natal homing) and seasonal dynamics to both simulation and estimation. Most dimensions now follow population, region, year, season, age, sex, and fleet.
* Added movement-after-mortality and continuous movement dynamics.
* Added continuous-time Markov chain (CTMC) movement estimation using preference functions.
* Added `move_expm_nsub` to `Setup_Mod_Movement`, allowing CTMC matrix exponentials to be evaluated using implicit backward Euler substeps, providing a faster alternative for differentiation while preserving non-negative fractions and exact survivor catch balance.
* Added a state-space treatment of the numbers at age via `NAA_re` in `Setup_Mod_Biologicals` and `Setup_Sim_NAA_state`. The log numbers themselves become random effects for ages two and older including the plus group, with the deterministic mortality and ageing step as the prediction they are scored against. Forms are `"iid"`, `"1dar1_a"`, `"1dar1_y"`, `"2dar1"`, `"3dcond"` and `"3dmarg"` over the age-year grid, composable with unstructured correlations across populations, regions, seasons and sexes.
* The state runs at season boundaries as well as year boundaries. `NAA_re_seasons` takes `"annual"` (the default, a state at season one only, so the numbers within a year stay deterministic and the innovation is purely annual), `"all"`, or an arbitrary set of season indices, which need not be contiguous so that states can be placed only in the seasons carrying observations. `NAA_re_season` adds an unstructured correlation across the active seasons, and `NAA_sigma_seasblk_spec` blocks the process error standard deviation over them. `ln_NAA`, `ln_sigmaNAA`, `naa_sigma_blocks`, `NAA_pred` and `NAA_scalar` all carry the season margin between year and age; saved input lists without it are promoted and hold the state at season one.
* The state-space numbers at age are not projected. `Do_Population_Projection()` advances projected numbers deterministically from the terminal year with recruitment the only stochastic element, so a forecast from a state-space fit omits that process error. The closed loop operating model does carry the state forward.
* Fixed the tag cohort rescale under a state-space fit: it indexed the state factor by years at liberty rather than by calendar year, so it was correct only for cohorts released in the first model year.
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
* `Setup_Mod_Rec` gained `R0_blocks`, which takes the same string forms as the selectivity blocks (`"none_Pop_p"`, `"Block_b_Year_a-e_Pop_p"`) and gives `ln_global_R0` a block dimension, together with `R0_ref_block` for the single value the initial age structure, the R0 prior, the `ln_rinit` penalty and the stock-recruit scale all read. Under mean recruitment an R0 block is a productivity regime; under a stock-recruit form it is the curve's scale. `R0_yr` is reported. A model without R0 blocks keeps one column and is byte-identical to before.
* An analytically solved catchability is now solved **within each catchability time block** rather than pooled over the whole series. `srv_q_type` / `fish_q_type` of `"arith"` or `"geo"` previously computed one value per region and fleet and discarded `srv_q_blocks` / `fish_q_blocks` silently, so a blocked analytic q looked configured and did nothing. Each block is now solved from the observations it owns, a block with none takes the pooled value, and a fleet with a single block is unchanged to machine precision.

### Reference points

* Reference points now support local density dependence for meta-population and natal homing models.
* Added uncertainty estimation for reference points.
* `Get_Reference_Points` now uses growth-derived weight at age from the model report rather than the data placeholder.
* `Do_Population_Projection()` gained the hooks needed to drive it as a management strategy loop. `rec_devs` supplies multiplicative recruitment deviations from outside, so a deterministic recruitment option becomes a stochastic one under deviations the caller draws and owns, and a set of replicates can share them across management procedures. A `HCR_function` that declares a `state` argument is handed the year's numbers at age, spawning biomass, total biomass and catch, so a rule can be written on mean weight or age structure rather than spawning biomass alone; rules without that argument are called exactly as before. Catch at age is now returned as `proj_CAA`.
* `Do_Population_Projection()` can now be differentiated through. The projection call chain carries RTMB's replacement operators, so a projection can be taped with `MakeTape()` or `MakeADFun()` and handed to an optimizer with an exact gradient, which allows an F schedule to be solved against an objective rather than scanned over a grid. `recruitment_opt = "inv_gauss"` and `fmort_opt = "Catch"` are refused when the inputs are AD types, since one draws random recruitment and the other inverts the catch target with a numerical solve, and both would otherwise return a gradient that is wrong rather than an error.

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
- Fixed the projection dropping the AD class when joining assessment and projected spawning biomass for the stock-recruit curve, which took recruitment off the tape.

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
