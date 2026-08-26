# Package index

## Model Setup

Functions for setting up a `SPoRC` model

- [`Setup_Mod_Biologicals()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Biologicals.md)
  : Set up biological inputs for the estimation model
- [`Setup_Mod_Catch_and_F()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Catch_and_F.md)
  : Set up fishing mortality, discard mortality, and catch observation
  inputs
- [`Setup_Mod_FishIdx_and_Comps()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_FishIdx_and_Comps.md)
  : Set up discards, fishery index, age composition, and length
  composition inputs
- [`Setup_Mod_Fishsel_and_Q()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Fishsel_and_Q.md)
  : Set up total and retained fishery selectivity and catchability
  specifications
- [`Setup_Mod_Movement()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Movement.md)
  : Set up movement model inputs and parameter structures
- [`Setup_Mod_Weighting()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Weighting.md)
  : Set likelihood and penalty weights for the estimation model
- [`Setup_Mod_Rec()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Rec.md)
  : Set up the recruitment module and associated processes
- [`Setup_Mod_SrvIdx_and_Comps()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_SrvIdx_and_Comps.md)
  : Set up observed survey indices and composition data
- [`Setup_Mod_Srvsel_and_Q()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Srvsel_and_Q.md)
  : Set up survey selectivity and catchability specifications
- [`Setup_Mod_Tagging()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Tagging.md)
  : Set up the conventional tagging module for model fitting
- [`Setup_Mod_Dim()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Dim.md)
  : Initialise model dimension settings

## Reference Points and Projections

Functions for deriving reference points and catch projections

- [`Do_Population_Projection()`](https://chengmatt.github.io/SPoRC/dev/reference/Do_Population_Projection.md)
  : Do Population Projections
- [`Get_Reference_Points()`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Reference_Points.md)
  : Compute fishing and biological reference points from an assessment
  or simulation
- [`get_key_quants()`](https://chengmatt.github.io/SPoRC/dev/reference/get_key_quants.md)
  : Generate Key Projection Quantities and Table Plot

## Simulation Setup

Functions for setting up simulations

- [`Setup_Sim_Biologicals()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Biologicals.md)
  : Set up biological parameter inputs for closed-loop simulation
- [`Setup_Sim_Containers()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Containers.md)
  : Initialise output containers for the operating model simulation
- [`Setup_Sim_Fishing()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Fishing.md)
  : Setup Simulation Fishing Inputs
- [`Setup_Sim_Rec()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Rec.md)
  : Set up recruitment dynamics for the operating model simulation
- [`Setup_Sim_Survey()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Survey.md)
  : Set up survey parameterisation for the operating model simulation
- [`Setup_Sim_Tagging()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Tagging.md)
  : Set up conventional tagging dynamics for the operating model
  simulation
- [`Setup_Sim_Dim()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Dim.md)
  : Initialise simulation dimension settings
- [`Setup_sim_env()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_sim_env.md)
  : Construct and populate a simulation execution environment
- [`run_annual_cycle()`](https://chengmatt.github.io/SPoRC/dev/reference/run_annual_cycle.md)
  : Run the annual cycle for a single simulation year
- [`Simulate_Pop_Static()`](https://chengmatt.github.io/SPoRC/dev/reference/Simulate_Pop_Static.md)
  : Simulate a static (open-loop) spatial age- and sex-structured
  population
- [`simulation_self_test()`](https://chengmatt.github.io/SPoRC/dev/reference/simulation_self_test.md)
  : Run a simulation self-test of a fitted RTMB estimation model
- [`simulation_data_to_SPoRC()`](https://chengmatt.github.io/SPoRC/dev/reference/simulation_data_to_SPoRC.md)
  : Extract simulation outputs into SPoRC estimation model format

## Closed-Loop Simulation

Functions for setting up and running closed-loop population simulations

- [`condition_closed_loop_simulations()`](https://chengmatt.github.io/SPoRC/dev/reference/condition_closed_loop_simulations.md)
  : Construct and Condition Closed-Loop Simulation Inputs
- [`get_closed_loop_reference_points()`](https://chengmatt.github.io/SPoRC/dev/reference/get_closed_loop_reference_points.md)
  : Get Closed Loop Reference Points
- [`catch_to_F_singlefleet()`](https://chengmatt.github.io/SPoRC/dev/reference/catch_to_F_singlefleet.md)
  : Convert a target catch to fishing mortality for a single fleet via
  bisection
- [`catch_to_F_multifleet()`](https://chengmatt.github.io/SPoRC/dev/reference/catch_to_F_multifleet.md)
  : Convert target catches to fishing mortality rates for multiple
  fleets

## Francis Reweighting

Functions for setting up Francis Reweighting

- [`run_francis()`](https://chengmatt.github.io/SPoRC/dev/reference/run_francis.md)
  : Run Iterative Francis Reweighting Procedure
- [`do_francis_reweighting()`](https://chengmatt.github.io/SPoRC/dev/reference/do_francis_reweighting.md)
  : Get Francis Weights

## Model Diagnostics

Functions for model diagnostics

- [`do_jitter()`](https://chengmatt.github.io/SPoRC/dev/reference/do_jitter.md)
  : Run Jitter Analysis for Model Diagnostics
- [`do_retrospective()`](https://chengmatt.github.io/SPoRC/dev/reference/do_retrospective.md)
  : Run Retrospective Diagnostics for RTMB Models
- [`get_retrospective_relative_difference()`](https://chengmatt.github.io/SPoRC/dev/reference/get_retrospective_relative_difference.md)
  : Derive relative difference from terminal year from a retrospective
  analysis.
- [`do_runs_test()`](https://chengmatt.github.io/SPoRC/dev/reference/do_runs_test.md)
  : Runs Test for Residual Randomness
- [`get_idx_fits()`](https://chengmatt.github.io/SPoRC/dev/reference/get_idx_fits.md)
  : Extract Index Fit Results
- [`get_comp_prop()`](https://chengmatt.github.io/SPoRC/dev/reference/get_comp_prop.md)
  : Get Composition Proportions from RTMB Output
- [`get_caal_fits()`](https://chengmatt.github.io/SPoRC/dev/reference/get_caal_fits.md)
  : Conditional age-at-length fits as mean age within each length bin
- [`get_osa()`](https://chengmatt.github.io/SPoRC/dev/reference/get_osa.md)
  : Compute OSA residuals for composition data
- [`plot_resids()`](https://chengmatt.github.io/SPoRC/dev/reference/plot_resids.md)
  : Plot OSA residuals from outputs of get_osa
- [`get_nLL_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_nLL_plot.md)
  : Get Plot of Negative Log Likelihood Values
- [`get_idx_fits_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_idx_fits_plot.md)
  : Get Index Fits Plot
- [`get_catch_fits_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_catch_fits_plot.md)
  : Get Catch and Discard Fits Plot
- [`get_at_age_fits_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_at_age_fits_plot.md)
  : Plot observed and predicted age-disaggregated observations
- [`get_retrospective_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_retrospective_plot.md)
  : Get Retrospective Plot
- [`do_likelihood_profile()`](https://chengmatt.github.io/SPoRC/dev/reference/do_likelihood_profile.md)
  : Run Likelihood Profile
- [`get_model_rep_from_mcmc()`](https://chengmatt.github.io/SPoRC/dev/reference/get_model_rep_from_mcmc.md)
  : Extract model report quantities from MCMC posterior samples
- [`marg_AIC()`](https://chengmatt.github.io/SPoRC/dev/reference/marg_AIC.md)
  : Compute the corrected marginal AIC (AICc)

## Utility

Functions for convenience

- [`fit_model()`](https://chengmatt.github.io/SPoRC/dev/reference/fit_model.md)
  : Fit a SPoRC RTMB model
- [`set_data_indicator_unused()`](https://chengmatt.github.io/SPoRC/dev/reference/set_data_indicator_unused.md)
  : Set data indicators to unused for specified years
- [`post_optim_sanity_checks()`](https://chengmatt.github.io/SPoRC/dev/reference/post_optim_sanity_checks.md)
  : Run post-optimisation convergence checks on a fitted SPoRC model
- [`get_par_est_info()`](https://chengmatt.github.io/SPoRC/dev/reference/get_par_est_info.md)
  : Extract and tabulate parameter estimates and metadata from a fitted
  SPoRC model
- [`rho_trans()`](https://chengmatt.github.io/SPoRC/dev/reference/rho_trans.md)
  : Transform a real-valued parameter to the interval (-1, 1)
- [`get_logistN_Sigma()`](https://chengmatt.github.io/SPoRC/dev/reference/get_logistN_Sigma.md)
  : Construct a logistic-normal covariance matrix

## Plotting

Functions for plotting

- [`get_ts_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_ts_plot.md)
  : Get Time Series Plots
- [`get_selex_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_selex_plot.md)
  : Get Fishery and Survey Selectivity Plots
- [`get_biological_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_biological_plot.md)
  : Get Plots of Biological Quantities
- [`get_data_fitted_plot()`](https://chengmatt.github.io/SPoRC/dev/reference/get_data_fitted_plot.md)
  : Get Data Fitted to Plot
- [`plot_all_basic()`](https://chengmatt.github.io/SPoRC/dev/reference/plot_all_basic.md)
  : Plotting Function for All Basic Quantities
- [`theme_sablefish()`](https://chengmatt.github.io/SPoRC/dev/reference/theme_sablefish.md)
  : ggplot2 theme for SPoRC plots

## Data

Assessment data inputs provided by SPoRC

- [`sgl_rg_sable_data`](https://chengmatt.github.io/SPoRC/dev/reference/sgl_rg_sable_data.md)
  : Sablefish data for single region case study
- [`mlt_rg_sable_data`](https://chengmatt.github.io/SPoRC/dev/reference/mlt_rg_sable_data.md)
  : Sablefish data for multi region (5 area) case study
- [`three_rg_sable_data`](https://chengmatt.github.io/SPoRC/dev/reference/three_rg_sable_data.md)
  : Sablefish data for multi region (3 area) case study
- [`sgl_rg_ebswp_data`](https://chengmatt.github.io/SPoRC/dev/reference/sgl_rg_ebswp_data.md)
  : EBS Walleye Pollock data for single region case study
- [`sgl_rg_dusky_data`](https://chengmatt.github.io/SPoRC/dev/reference/sgl_rg_dusky_data.md)
  : Dusky data for single region assessment case study
- [`sgl_rg_goa_nork_data`](https://chengmatt.github.io/SPoRC/dev/reference/sgl_rg_goa_nork_data.md)
  : GOA Northern Rockfish data for single region case study
- [`sgl_rg_bsai_nork_data`](https://chengmatt.github.io/SPoRC/dev/reference/sgl_rg_bsai_nork_data.md)
  : BSAI Northern Rockfish data for single region case study
- [`sgl_rg_bsai_pop_data`](https://chengmatt.github.io/SPoRC/dev/reference/sgl_rg_bsai_pop_data.md)
  : BSAI Pacific Ocean Perch data for single region case study
- [`sgl_rg_bsai_atka_data`](https://chengmatt.github.io/SPoRC/dev/reference/sgl_rg_bsai_atka_data.md)
  : BSAI Atka Mackerel data for single region case study
- [`sgl_rg_wc_sablefish_data`](https://chengmatt.github.io/SPoRC/dev/reference/sgl_rg_wc_sablefish_data.md)
  : West Coast Sablefish data for single region case study
- [`sgl_rg_bsai_nrs_data`](https://chengmatt.github.io/SPoRC/dev/reference/sgl_rg_bsai_nrs_data.md)
  : BSAI northern rock sole data for single region assessment case study
- [`sgl_rg_rebs_data`](https://chengmatt.github.io/SPoRC/dev/reference/sgl_rg_rebs_data.md)
  : BSAI Blackspotted and Rougheye Rockfish data for single region case
  study
- [`mlt_rg_goa_rex_data`](https://chengmatt.github.io/SPoRC/dev/reference/mlt_rg_goa_rex_data.md)
  : GOA rex sole bridge data (Model 25.1)
- [`sgl_rg_ebs_pcod_data`](https://chengmatt.github.io/SPoRC/dev/reference/sgl_rg_ebs_pcod_data.md)
  : EBS Pacific cod bridge data (Model 24.1)

## Fitted Model Objects

Fitted model output provided by SPoRC

- [`sgl_rg_sable_rep`](https://chengmatt.github.io/SPoRC/dev/reference/sgl_rg_sable_rep.md)
  : Sablefish report for single region case study
- [`mlt_rg_sable_rep`](https://chengmatt.github.io/SPoRC/dev/reference/mlt_rg_sable_rep.md)
  : Sablefish report for 5 region case study
- [`three_rg_sable_rep`](https://chengmatt.github.io/SPoRC/dev/reference/three_rg_sable_rep.md)
  : Sablefish report for 3 region case study
- [`dusky_rtmb_model`](https://chengmatt.github.io/SPoRC/dev/reference/dusky_rtmb_model.md)
  : Dusky model outputs from single regino model
