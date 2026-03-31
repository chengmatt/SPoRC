# Package index

## Model Setup

Functions for setting up a `SPoRC` model

- [`Setup_Mod_Biologicals()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Biologicals.md)
  : Setup biological inputs for estimation model
- [`Setup_Mod_Catch_and_F()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Catch_and_F.md)
  : Setup fishing mortality and catch observations
- [`Setup_Mod_FishIdx_and_Comps()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_FishIdx_and_Comps.md)
  : Setup observed fishery indices and composition data (age and length
  comps)
- [`Setup_Mod_Fishsel_and_Q()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Fishsel_and_Q.md)
  : Setup fishery selectivity and catchability specifications
- [`Setup_Mod_Movement()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Movement.md)
  : Setup Movement Processes for SPoRC
- [`Setup_Mod_Weighting()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Weighting.md)
  : Set up SPoRC model weighting
- [`Setup_Mod_Rec()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Rec.md)
  : Setup model objects for specifying recruitment module and associated
  processes
- [`Setup_Mod_SrvIdx_and_Comps()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_SrvIdx_and_Comps.md)
  : Setup observed survey indices and composition data (age and length
  comps)
- [`Setup_Mod_Srvsel_and_Q()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Srvsel_and_Q.md)
  : Setup survey selectivity and catchability specifications
- [`Setup_Mod_Tagging()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Tagging.md)
  : Setup tagging processes and parameters
- [`Setup_Mod_Dim()`](https://chengmatt.github.io/SPoRC/reference/Setup_Mod_Dim.md)
  : Set up model dimensions

## Reference Points and Projections

Functions for deriving reference points and catch projections

- [`Do_Population_Projection()`](https://chengmatt.github.io/SPoRC/reference/Do_Population_Projection.md)
  : Do Population Projections
- [`Get_Reference_Points()`](https://chengmatt.github.io/SPoRC/reference/Get_Reference_Points.md)
  : Wrapper function to get reference points
- [`get_key_quants()`](https://chengmatt.github.io/SPoRC/reference/get_key_quants.md)
  : Generate Key Projection Quantities and Table Plot

## Simulation Setup

Functions for setting up simulations

- [`Setup_Sim_Biologicals()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Biologicals.md)
  : Set up simulation containers and inputs for biological parameters
- [`Setup_Sim_Containers()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Containers.md)
  : Setup containers for simulation and output
- [`Setup_Sim_Fishing()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Fishing.md)
  : Setup values and dimensions of fishing processes
- [`Setup_Sim_Rec()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Rec.md)
  : Set up recruitment dynamics for simulation
- [`Setup_Sim_Survey()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Survey.md)
  : Setup values for survey parameterization
- [`Setup_Sim_Tagging()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Tagging.md)
  : Set up simulated tagging dynamics
- [`Setup_Sim_Dim()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Dim.md)
  : Initialize Simulation Dimension Settings
- [`Setup_sim_env()`](https://chengmatt.github.io/SPoRC/reference/Setup_sim_env.md)
  : Constructs simulation objects in a new simulation environment for
  use in simulation functions
- [`run_annual_cycle()`](https://chengmatt.github.io/SPoRC/reference/run_annual_cycle.md)
  : Run Annual Cycle in Simulation Environment
- [`Simulate_Pop_Static()`](https://chengmatt.github.io/SPoRC/reference/Simulate_Pop_Static.md)
  : Simulates a static spatial, sex, and age-structured population (no
  feedback loop)
- [`simulation_self_test()`](https://chengmatt.github.io/SPoRC/reference/simulation_self_test.md)
  : Conduct a Simulation Self Test
- [`simulation_data_to_SPoRC()`](https://chengmatt.github.io/SPoRC/reference/simulation_data_to_SPoRC.md)
  : Extract simulation data into SPoRC format

## Closed-Loop Simulation

Functions for setting up and running closed-loop population simulations

- [`condition_closed_loop_simulations()`](https://chengmatt.github.io/SPoRC/reference/condition_closed_loop_simulations.md)
  : Set up simulation list for closed-loop projections
- [`get_closed_loop_reference_points()`](https://chengmatt.github.io/SPoRC/reference/get_closed_loop_reference_points.md)
  : Get Closed Loop Reference Points
- [`catch_to_F_singlefleet()`](https://chengmatt.github.io/SPoRC/reference/catch_to_F_singlefleet.md)
  : Go from TAC to Fishing Mortality using bisection for when a single
  fishery fleet exists
- [`catch_to_F_multifleet()`](https://chengmatt.github.io/SPoRC/reference/catch_to_F_multifleet.md)
  : Solve for fishing mortality rates that achieve target catches for
  multiple fleets

## Francis Reweighting

Functions for setting up Francis Reweighting

- [`run_francis()`](https://chengmatt.github.io/SPoRC/reference/run_francis.md)
  : Run Iterative Francis Reweighting Procedure
- [`do_francis_reweighting()`](https://chengmatt.github.io/SPoRC/reference/do_francis_reweighting.md)
  : Get Francis Weights

## Model Diagnostics

Functions for model diagnostics

- [`do_jitter()`](https://chengmatt.github.io/SPoRC/reference/do_jitter.md)
  : Run Jitter Analysis
- [`do_retrospective()`](https://chengmatt.github.io/SPoRC/reference/do_retrospective.md)
  : Run retrospective analyses for RTMB models
- [`get_retrospective_relative_difference()`](https://chengmatt.github.io/SPoRC/reference/get_retrospective_relative_difference.md)
  : Derive relative difference from terminal year from a retrospective
  analysis.
- [`do_runs_test()`](https://chengmatt.github.io/SPoRC/reference/do_runs_test.md)
  : Runs test function taken from SS3 diags.
- [`get_idx_fits()`](https://chengmatt.github.io/SPoRC/reference/get_idx_fits.md)
  : Gets index fits results
- [`get_comp_prop()`](https://chengmatt.github.io/SPoRC/reference/get_comp_prop.md)
  : Gets composition data proportions normalized according to the
  assessment specifications from RTMB
- [`get_osa()`](https://chengmatt.github.io/SPoRC/reference/get_osa.md)
  : Compute OSA residuals for composition data
- [`plot_resids()`](https://chengmatt.github.io/SPoRC/reference/plot_resids.md)
  : Plots OSA residuals from outputs from get_osa. Much of this code is
  taken from the afscOM package, but with modificaitons to plot
  features.
- [`get_nLL_plot()`](https://chengmatt.github.io/SPoRC/reference/get_nLL_plot.md)
  : Get plot of negative log likelihood values
- [`get_idx_fits_plot()`](https://chengmatt.github.io/SPoRC/reference/get_idx_fits_plot.md)
  : Get Index Fits Plot
- [`get_catch_fits_plot()`](https://chengmatt.github.io/SPoRC/reference/get_catch_fits_plot.md)
  : Title Get Catch Fits Plot
- [`get_retrospective_plot()`](https://chengmatt.github.io/SPoRC/reference/get_retrospective_plot.md)
  : Get Retrospective Plot
- [`do_likelihood_profile()`](https://chengmatt.github.io/SPoRC/reference/do_likelihood_profile.md)
  : Run Likelihood Profile
- [`get_model_rep_from_mcmc()`](https://chengmatt.github.io/SPoRC/reference/get_model_rep_from_mcmc.md)
  : Extract model report from MCMC posterior samples
- [`marg_AIC()`](https://chengmatt.github.io/SPoRC/reference/marg_AIC.md)
  : Calculate the Corrected marginal AIC (AICc) from Optimization
  Results

## Utility

Functions for convenience

- [`fit_model()`](https://chengmatt.github.io/SPoRC/reference/fit_model.md)
  : Run RTMB model
- [`set_data_indicator_unused()`](https://chengmatt.github.io/SPoRC/reference/set_data_indicator_unused.md)
  : Set Data Indicators to Unused for Specified Years
- [`post_optim_sanity_checks()`](https://chengmatt.github.io/SPoRC/reference/post_optim_sanity_checks.md)
  : Post Optimization Model Convergence Checks
- [`get_par_est_info()`](https://chengmatt.github.io/SPoRC/reference/get_par_est_info.md)
  : Helper function for extracting parameter information and names from
  TMB
- [`rho_trans()`](https://chengmatt.github.io/SPoRC/reference/rho_trans.md)
  : Title Constrains value between -1 and 1
- [`get_logistN_Sigma()`](https://chengmatt.github.io/SPoRC/reference/get_logistN_Sigma.md)
  : Construct logistic-normal covariance matrix

## Plotting

Functions for plotting

- [`get_ts_plot()`](https://chengmatt.github.io/SPoRC/reference/get_ts_plot.md)
  : Get Time Series Plots
- [`get_selex_plot()`](https://chengmatt.github.io/SPoRC/reference/get_selex_plot.md)
  : Get Fishery and Survey Selectivity Plots
- [`get_biological_plot()`](https://chengmatt.github.io/SPoRC/reference/get_biological_plot.md)
  : Get Plots of Biological Quantities
- [`get_data_fitted_plot()`](https://chengmatt.github.io/SPoRC/reference/get_data_fitted_plot.md)
  : Get Data Fitted to Plot
- [`plot_all_basic()`](https://chengmatt.github.io/SPoRC/reference/plot_all_basic.md)
  : Plotting function for all basic quantities
- [`theme_sablefish()`](https://chengmatt.github.io/SPoRC/reference/theme_sablefish.md)
  : ggplot theme for sablefish

## Data

Data objects provided by SPoRC

- [`sgl_rg_sable_data`](https://chengmatt.github.io/SPoRC/reference/sgl_rg_sable_data.md)
  : Sablefish data for single region case study
- [`mlt_rg_sable_data`](https://chengmatt.github.io/SPoRC/reference/mlt_rg_sable_data.md)
  : Sablefish data for multi region (5 area) case study
- [`sgl_rg_ebswp_data`](https://chengmatt.github.io/SPoRC/reference/sgl_rg_ebswp_data.md)
  : EBS Walleye Pollock data for single region case study
- [`three_rg_sable_data`](https://chengmatt.github.io/SPoRC/reference/three_rg_sable_data.md)
  : Sablefish data for multi region (3 area) case study
- [`sgl_rg_sable_rep`](https://chengmatt.github.io/SPoRC/reference/sgl_rg_sable_rep.md)
  : Sablefish report for single region case study
- [`mlt_rg_sable_rep`](https://chengmatt.github.io/SPoRC/reference/mlt_rg_sable_rep.md)
  : Sablefish report for 5 region case study
- [`three_rg_sable_rep`](https://chengmatt.github.io/SPoRC/reference/three_rg_sable_rep.md)
  : Sablefish report for 3 region case study
- [`sgl_rg_dusky_data`](https://chengmatt.github.io/SPoRC/reference/sgl_rg_dusky_data.md)
  : Dusky data for single region assessment case study
- [`dusky_rtmb_model`](https://chengmatt.github.io/SPoRC/reference/dusky_rtmb_model.md)
  : Dusky model outputs from single regino model
