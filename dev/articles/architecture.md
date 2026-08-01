# Architecture

This document is for people editing the `SPoRC` codebase: contributors,
successors, future you. The goal is to keep the package navigable as it
grows, especially as spatial model development adds more moving parts.
Update it whenever the pipeline below changes shape (a new `Setup_Mod_*`
stage, a new section of the objective function, a new post fit
diagnostic). A stale architecture document should be treated as a bug.

## Mental model

`SPoRC` runs in three phases, in this order.

1.  **Setup.** A chain of R functions builds up three plain lists:
    `data`, `parameters` and `mapping`. No likelihood is evaluated yet.
    This phase is data wrangling and input validation only.
2.  **Fit.** `data` and `parameters` go to `RTMB::MakeADFun` wrapped
    around a single objective function,
    [`SPoRC_rtmb()`](https://chengmatt.github.io/SPoRC/dev/reference/SPoRC_rtmb.md),
    which is optimized with `nlminb` plus Newton refinement.
3.  **Post fit.** The fitted object (`obj$rep`, `obj$sdrep`,
    `obj$optim`) is passed to diagnostics, reference points, plotting
    and resampling routines.

&nbsp;

    Setup_Mod_Dim()
       -> Setup_Mod_Rec()
       -> Setup_Mod_Biologicals()
       -> Setup_Mod_Movement()
       -> Setup_Mod_Tagging()
       -> Setup_Mod_Catch_and_F()
       -> Setup_Mod_FishIdx_and_Comps()   [+ Setup_Mod_Discard_Comps()]
       -> Setup_Mod_Fishsel_and_Q()       [+ Setup_Mod_Retsel()]
       -> Setup_Mod_SrvIdx_and_Comps()
       -> Setup_Mod_Srvsel_and_Q()
       -> Setup_Mod_Weighting()
            |
            v  input_list$data, input_list$par, input_list$map
       fit_model(data, parameters, mapping)
            |  RTMB::MakeADFun(cmb(SPoRC_rtmb, data), parameters, map = mapping)
            |  nlminb() + Newton steps
            v
       obj  (obj$rep, obj$sdrep, obj$optim)
            |
            +--> Get_Reference_Points()          SPR and MSY
            +--> get_idx_fits() / get_osa()      fit extraction and OSA residuals
            +--> plot_* / get_*_plot             figures and tables
            +--> do_retrospective(), do_jitter(), do_likelihood_profile(), run_francis()
            +--> Do_Population_Projection(), condition_closed_loop_simulations()

### The `input_list` accumulator

Every `Setup_Mod_*` function takes an `input_list` (initialized by
[`Setup_Mod_Dim()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Dim.md)
as `list(data = list(), par = list(), map = list())`) and returns it
with more keys filled in. This is the central structure of the package:
setup functions are pipeline stages that thread one growing list through
each other, not independent builders.

Order matters. Later stages read dimensions (`input_list$data$n_regions`
and friends) set by earlier ones, and default arguments frequently
reference `input_list$data$...` in their own definitions. See
[`Setup_Mod_Rec()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Rec.md)’s
defaults for an example.

At the end of the chain `input_list$data`, `input_list$par` and
`input_list$map` are unpacked and passed to
[`fit_model()`](https://chengmatt.github.io/SPoRC/dev/reference/fit_model.md).
The `Setup_Sim_*` functions live in the same files as their
`Setup_Mod_*` counterparts and follow the identical accumulator pattern,
but build inputs for the operating model instead of for fitting.

### The objective function

`R/model_objective.R` holds the one RTMB objective function for every
model configuration the package supports. There is no per model
branching at the file level. The objective always runs end to end and is
regulated by the switches set during setup, such as `rec_model`,
`move_type`, and the likelihood and OSA flags.

It moves through, in order: parameter transforms (movement, natural
mortality, selectivity), then mortality, recruitment, initial age
structure, population projection, the observation models (fishery,
survey, tagging), the likelihood components (retained and discarded
catch, indices, compositions, tags), and finally priors and penalties.

Each section hands off to a module in one of the other `model_*.R`
files, so the objective reads as an order of operations rather than as
the arithmetic itself. If you are adding a new data source or process,
this is where its likelihood contribution gets wired in, and the
matching `Setup_Mod_*` function is where its data, parameters and
mapping get prepared beforehand.

Within `model_objective.R` the `## Section` and `### Subsection` comment
banner style is used consistently. Keep using it if you extend the file,
since it is the only navigation aid inside a function of this size.

Note that `R/model_transition.R` is worth knowing about before anything
else. It computes what happens to a vector of abundance over one season,
given movement, mortality, season duration and a movement timing. Eleven
other files go through it, including the objective function’s dynamics,
the reference point solvers and the operating model.

## Naming conventions

### Files

The prefix on a file name says which phase it belongs to.

| Prefix | Phase |
|----|----|
| `setup_` | Builds `data`, `parameters` and `mapping`, one process or data source per file |
| `model_` | The objective function and the modules it delegates to |
| `refpts_`, `projection` | Reference points and forward projection off a fitted model |
| `sim_` | The operating model: self testing and closed loop simulation |
| `diag_` | Post fit diagnostics and residuals |
| `plot_` | Figures and summary tables |
| `utils_` | Helpers shared across phases |
| `data_` | Roxygen documentation for the bundled datasets |

### Functions

| Prefix | Role |
|----|----|
| `Setup_Mod_*` | Builds inputs for a fitted object |
| `Setup_Sim_*` | Same, for the operating model. Lives beside its `Setup_Mod_*` counterpart |
| `do_*_mapping` | Internal helper that builds the RTMB `map` for one parameter block. This is how parameters get fixed, shared across blocks, or turned into random effect deviations |
| `Get_*` | Pulls a derived quantity out of set up data, parameters or a fitted object |
| `do_*` (top level) | A post fit procedure that refits or perturbs a fitted model |
| `get_*_plot`, `get_*_fits` | Post fit plotting and tabulation |

A `do_*_mapping` helper lives in `setup_mapping.R` only if more than one
`setup_*` file calls it. Single caller helpers stay in the file that
uses them.

## File map

Every script, what it defines, and what it talks to. This table is
generated from the actual call graph, and shows “what interacts with
what”. “Calls into” lists the files whose functions this one calls;
“called from” is the reverse.

![Directed graph of SPoRC pipeline stages, with arrows pointing from
each stage to the stages it calls
into.](figures/arch_stage_diagram-1.png)

Edge labels count how many file to file dependencies sit behind each
arrow. `Helpers` is a pure sink and `Plotting` a pure source, which is
the shape you want: the objective function never reaches up into the
post fit machinery.

Below the same information at file resolution. A filled cell means the
file on that row calls a function defined in the file in that column,
and cells are grouped by stage. Empty rows call nothing else in the
package, empty columns are called by nothing.

![Dependency matrix of the 52 files in R/, grouped by pipeline stage,
where a filled cell means the row file calls a function defined in the
column file.](figures/arch_dependency_matrix-1.png)

The `model_transition` column is the densest in the matrix, which is the
point made above: everything that moves fish through a season goes
through that one file.

### Per stage detail

#### Setup: build data, parameters and mapping

| Script | Defines | Calls into | Called from |
|----|----|----|----|
| `setup_biologicals.R` | `do_natmort_mapping`, `Setup_Mod_Biologicals`, `Setup_Sim_Biologicals` | `setup_checks.R`, `utils_setup.R` | `sim_closed_loop.R`, `sim_self_test.R` |
| `setup_checks.R` | `check_data_dimensions`, `check_sim_dimensions` | nothing | `setup_biologicals.R`, `setup_fishery_catch.R`, `setup_fishery_comps.R`, `setup_fishery_selectivity.R`, `setup_movement.R`, `setup_recruitment.R`, `setup_survey.R`, `setup_tagging.R` |
| `setup_dimensions.R` | `Setup_Mod_Dim`, `Setup_Sim_Dim` | `utils_setup.R` | `sim_closed_loop.R`, `sim_self_test.R` |
| `setup_fishery_catch.R` | `do_dmr_dev_mapping`, `do_dmr_mean_mapping`, `do_Fdev_rho_mapping`, `do_Fmort_mapping`, `do_sigma_dmr_mapping`, `do_sigmaC_mapping`, `do_sigmaC_pop_mapping`, `do_sigmaD_mapping`, `do_sigmaD_pop_mapping`, `do_sigmaF_mapping`, `Setup_Mod_Catch_and_F`, `Setup_Sim_Fishing` | `setup_checks.R`, `setup_mapping.R`, `utils_setup.R` | `sim_closed_loop.R`, `sim_self_test.R` |
| `setup_fishery_comps.R` | `Setup_Mod_Discard_Comps`, `Setup_Mod_FishIdx_and_Comps` | `setup_checks.R`, `setup_mapping.R`, `utils_setup.R` |  |
| `setup_fishery_selectivity.R` | `Setup_Mod_Fishsel_and_Q`, `Setup_Mod_Retsel` | `setup_checks.R`, `setup_mapping.R`, `utils_math.R`, `utils_setup.R` |  |
| `setup_mapping.R` | `build_pe_map`, `build_shared_spec_map`, `do_comp_corr_pars_mapping`, `do_comp_theta_mapping`, `do_fixed_sel_pars_mapping`, `do_q_mapping`, `do_sel_devs_mapping`, `do_sel_pe_pars_mapping`, `sync_dev_map_data` | `utils_setup.R` | `model_fit.R`, `setup_fishery_catch.R`, `setup_fishery_comps.R`, `setup_fishery_selectivity.R`, `setup_movement.R`, `setup_recruitment.R`, `setup_survey.R` |
| `setup_movement.R` | `do_cont_vary_move_mapping`, `do_move_pars_mapping`, `Setup_Mod_Movement` | `model_movement.R`, `setup_checks.R`, `setup_mapping.R`, `utils_setup.R` |  |
| `setup_recruitment.R` | `do_h_mapping`, `do_InitDevs_mapping`, `do_rec_region_prop_mapping`, `do_rec_seas_prop_mapping`, `do_RecDevs_mapping`, `do_sexratio_pars_mapping`, `do_sigmaR_mapping`, `do_stray_rate_mapping`, `Setup_Mod_Rec`, `Setup_Sim_Rec` | `setup_checks.R`, `setup_mapping.R`, `utils_setup.R` | `sim_closed_loop.R`, `sim_self_test.R` |
| `setup_sim_containers.R` | `Setup_Sim_Containers` | nothing | `sim_closed_loop.R`, `sim_self_test.R` |
| `setup_survey.R` | `Setup_Mod_SrvIdx_and_Comps`, `Setup_Mod_Srvsel_and_Q`, `Setup_Sim_Survey` | `setup_checks.R`, `setup_mapping.R`, `utils_math.R`, `utils_setup.R` | `sim_closed_loop.R`, `sim_self_test.R` |
| `setup_tagging.R` | `do_conv_init_tag_mort_mapping`, `do_conv_tag_fish_reporting_pars_mapping`, `do_conv_tag_shed_mapping`, `do_conv_tag_theta_mapping`, `recycle_tag_event_par`, `Setup_Mod_Tagging`, `Setup_Sim_Tagging` | `setup_checks.R`, `utils_setup.R` | `sim_closed_loop.R`, `sim_self_test.R` |
| `setup_weighting.R` | `Setup_Mod_Weighting` | `utils_setup.R` |  |

#### Objective function and its modules

| Script | Defines | Calls into | Called from |
|----|----|----|----|
| `model_biomass.R` | `compute_biom_y`, `derive_proj_biom` | `model_transition.R` | `model_population_dynamics.R`, `projection.R` |
| `model_distributions.R` | `dbeta_symmetric`, `ddirichlet`, `ddirmult`, `dlogistnormal`, `dnbinom_robust_noint`, `dpois_noint`, `get_beta_scaled_pars` | nothing | `model_lik_comps.R`, `model_lik_tags.R`, `model_priors_penalties.R` |
| `model_fit.R` | `cmb`, `fit_model` | `setup_mapping.R` | `diag_francis.R`, `diag_jitter.R`, `diag_likelihood_profile.R`, `diag_retrospective.R`, `refpts_main.R`, `sim_self_test.R` |
| `model_init_naa.R` | `Get_Init_NAA` | `model_transition.R` | `model_objective.R`, `sim_population.R` |
| `model_lik_comps.R` | `eval_comp_osa`, `Get_Comp_Likelihoods`, `Get_Comp_Likelihoods_OSA`, `pack_comp_osa` | `model_distributions.R`, `model_osa.R`, `utils_math.R` | `diag_osa_residuals.R`, `model_objective.R` |
| `model_lik_tags.R` | `eval_tag_osa`, `get_conv_tag_likelihoods`, `pack_tag_osa`, `tag_fam_of`, `tag_grid` | `model_distributions.R`, `model_osa.R` | `diag_osa_residuals.R`, `model_objective.R` |
| `model_movement.R` | `Get_Movement`, `get_movement_dp_design_matrix` | nothing | `model_objective.R`, `setup_movement.R` |
| `model_objective.R` | `maintain_backwards_compatibility`, `SPoRC_rtmb` | `model_init_naa.R`, `model_lik_comps.R`, `model_lik_tags.R`, `model_movement.R`, `model_obs_fishery_survey.R`, `model_obs_tagging.R`, `model_population_dynamics.R`, `model_priors_penalties.R`, `model_selectivity.R`, `utils_setup.R` |  |
| `model_obs_fishery_survey.R` | `get_fishery_observation_model`, `get_survey_observation_model` | `model_transition.R` | `model_objective.R` |
| `model_obs_tagging.R` | `get_tag_mort`, `get_tagging_observation_model`, `release_conv_tag_attr` | `model_transition.R` | `model_objective.R`, `sim_observations.R` |
| `model_osa.R` | `ddirmult_osa`, `ddirmult2`, `dmultinom_osa`, `osa_extract_cdf`, `osa_extract_keep`, `osa_extract_values`, `osa_extract_x`, `osa_pbetabinom`, `osa_pbinom`, `osa_squeeze` | nothing | `model_lik_comps.R`, `model_lik_tags.R` |
| `model_population_dynamics.R` | `get_population_projection` | `model_biomass.R`, `model_recruitment.R`, `model_transition.R` | `model_objective.R` |
| `model_precision.R` | `Get_3d_precision` | nothing | `model_priors_penalties.R` |
| `model_priors_penalties.R` | `get_dmr_penalty`, `Get_Fdev_PE_loglik`, `Get_move_PE_loglik`, `get_movement_dirichlet_prior`, `get_natmort_prior`, `get_q_prior`, `get_r0_prior`, `get_recruitment_penalty`, `get_recruitment_proportion_priors`, `Get_sel_PE_loglik`, `get_selex_fixed_prior`, `Get_Selex_Smoothness_Penalty`, `get_steepness_prior`, `get_tagrep_prior` | `model_distributions.R`, `model_precision.R`, `utils_math.R` | `model_objective.R` |
| `model_recruitment.R` | `Get_Det_Recruitment` | `model_transition.R` | `model_population_dynamics.R`, `projection.R`, `sim_population.R` |
| `model_selectivity.R` | `Get_Selex`, `Get_Selex_Array` | nothing | `model_objective.R` |
| `model_transition.R` | `advance_seas`, `build_seas_operator`, `catch_at_age`, `integrate_seas_abundance`, `seas_operator_and_integral`, `spawn_state`, `survey_state` | nothing | `model_biomass.R`, `model_init_naa.R`, `model_obs_fishery_survey.R`, `model_obs_tagging.R`, `model_population_dynamics.R`, `model_recruitment.R`, `projection.R`, `refpts_main.R`, `refpts_msy.R`, `refpts_spr.R`, `sim_observations.R`, `sim_population.R` |

#### Reference points and projection

| Script | Defines | Calls into | Called from |
|----|----|----|----|
| `projection.R` | `build_proj_F`, `Do_Population_Projection`, `proj_catch_at_F`, `proj_log_catch_resid`, `proj_target_catch`, `run_proj_year`, `solve_proj_F_catch`, `solve_proj_year_F` | `model_biomass.R`, `model_recruitment.R`, `model_transition.R`, `sim_random_variates.R` | `plot_figures_tables.R` |
| `refpts_main.R` | `build_plus_group_T`, `Get_Reference_Points`, `optim_ref_pts`, `solve_plus_group` | `model_fit.R`, `model_transition.R` | `plot_figures_tables.R`, `refpts_msy.R`, `refpts_spr.R`, `sim_closed_loop.R` |
| `refpts_msy.R` | `global_BH_Fmsy`, `local_BH_Fmsy_multipop`, `local_BH_Fmsy_sglpop`, `single_region_BH_Fmsy` | `model_transition.R`, `refpts_main.R` |  |
| `refpts_spr.R` | `global_SPR`, `single_region_SPR` | `model_transition.R`, `refpts_main.R` |  |

#### Operating model

| Script | Defines | Calls into | Called from |
|----|----|----|----|
| `sim_closed_loop.R` | `catch_to_F_multifleet`, `catch_to_F_singlefleet`, `condition_closed_loop_simulations`, `get_closed_loop_reference_points` | `refpts_main.R`, `setup_biologicals.R`, `setup_dimensions.R`, `setup_fishery_catch.R`, `setup_recruitment.R`, `setup_sim_containers.R`, `setup_survey.R`, `setup_tagging.R`, `utils_postfit.R`, `utils_setup.R` |  |
| `sim_observations.R` | `generate_fishery_catch_comp_idx`, `generate_fishery_conv_tags_recap`, `generate_survey_comp_idx`, `marginalize_conv_fish_tags`, `predict_sim_fish_iss_fmort`, `release_conv_tags`, `simulate_comps`, `simulate_conv_tag_fish_recaptures` | `model_obs_tagging.R`, `model_transition.R`, `sim_random_variates.R`, `utils_math.R` | `sim_population.R` |
| `sim_population.R` | `apply_pop_dy`, `compute_biom_y_sim`, `generate_initial_age_structure`, `generate_recruitment`, `run_annual_cycle`, `Simulate_Pop_Static` | `model_init_naa.R`, `model_recruitment.R`, `model_transition.R`, `sim_observations.R`, `sim_setup.R` | `sim_self_test.R` |
| `sim_random_variates.R` | `rdirM`, `rinvgauss_rec`, `rlogistnormal` | `utils_math.R` | `projection.R`, `sim_observations.R` |
| `sim_self_test.R` | `simulation_data_to_SPoRC`, `simulation_self_test` | `model_fit.R`, `setup_biologicals.R`, `setup_dimensions.R`, `setup_fishery_catch.R`, `setup_recruitment.R`, `setup_sim_containers.R`, `setup_survey.R`, `setup_tagging.R`, `sim_population.R`, `utils_postfit.R` |  |
| `sim_setup.R` | `Setup_sim_env` | nothing | `sim_population.R` |

#### Diagnostics

| Script | Defines | Calls into | Called from |
|----|----|----|----|
| `diag_fits.R` | `get_comp_prop`, `get_idx_fits`, `Restrc_Comps` | nothing | `diag_francis.R`, `plot_figures_tables.R` |
| `diag_francis.R` | `do_francis_reweighting`, `get_francis_weights`, `run_francis`, `safe_inv_var` | `diag_fits.R`, `model_fit.R` | `diag_retrospective.R` |
| `diag_jitter.R` | `do_jitter` | `model_fit.R` |  |
| `diag_likelihood_profile.R` | `do_likelihood_profile` | `model_fit.R` |  |
| `diag_osa_residuals.R` | `comp_osa_field_map`, `get_osa`, `index_osa_field_map`, `osa_keep_subset`, `osa_one_step_predict`, `plot_resids`, `run_external_comp_osa`, `run_internal_comp_osa`, `run_internal_index_osa`, `run_internal_tag_osa`, `validate_osa_method` | `model_lik_comps.R`, `model_lik_tags.R` |  |
| `diag_retrospective.R` | `do_retrospective`, `get_retrospective_relative_difference`, `truncate_yr` | `diag_francis.R`, `model_fit.R` | `plot_figures_tables.R` |
| `diag_runs_test.R` | `do_runs_test` | nothing |  |

#### Plotting

| Script | Defines | Calls into | Called from |
|----|----|----|----|
| `plot_figures_tables.R` | `get_biological_plot`, `get_catch_fits_plot`, `get_data_fitted_plot`, `get_idx_fits_plot`, `get_key_quants`, `get_nLL_plot`, `get_retrospective_plot`, `get_selex_plot`, `get_ts_plot`, `plot_all_basic`, `theme_sablefish` | `diag_fits.R`, `diag_retrospective.R`, `projection.R`, `refpts_main.R`, `utils_setup.R` |  |

#### Shared helpers

| Script | Defines | Calls into | Called from |
|----|----|----|----|
| `utils_math.R` | `get_AR1_CorrMat`, `get_Constant_CorrMat`, `get_logistN_Sigma`, `Get_Natural_Cubic_Spline_Weights`, `rho_trans` | nothing | `model_lik_comps.R`, `model_priors_penalties.R`, `setup_fishery_selectivity.R`, `setup_survey.R`, `sim_observations.R`, `sim_random_variates.R` |
| `utils_postfit.R` | `get_model_rep_from_mcmc`, `get_optim_param_list`, `get_par_est_info`, `marg_AIC`, `post_optim_sanity_checks` | nothing | `sim_closed_loop.R`, `sim_self_test.R` |
| `utils_setup.R` | `collect_message`, `convert_to_numeric`, `extend_years`, `resolve_sel_pen_wts`, `safe_extract`, `set_data_indicator_unused` | nothing | `model_objective.R`, `plot_figures_tables.R`, `setup_biologicals.R`, `setup_dimensions.R`, `setup_fishery_catch.R`, `setup_fishery_comps.R`, `setup_fishery_selectivity.R`, `setup_mapping.R`, `setup_movement.R`, `setup_recruitment.R`, `setup_survey.R`, `setup_tagging.R`, `setup_weighting.R`, `sim_closed_loop.R` |

#### Package data

| Script                 | Defines              | Calls into | Called from |
|------------------------|----------------------|------------|-------------|
| `data_documentation.R` | (documentation only) | nothing    |             |

## How the code is layered

The code is structured as follows and interacts in a hopefully intuitive
fashion. Reading it as layers, lowest first:

1.  `utils_math.R`, `utils_setup.R`, `utils_postfit.R`,
    `setup_checks.R`, `model_transition.R`, `model_distributions.R`,
    `model_osa.R`, `model_precision.R`. These call nothing else in the
    package.
2.  The model modules (`model_recruitment.R`, `model_biomass.R`,
    `model_init_naa.R`, `model_movement.R`, `model_selectivity.R`,
    `model_obs_*.R`, `model_lik_*.R`, `model_priors_penalties.R`,
    `model_population_dynamics.R`).
3.  `model_objective.R`, which calls only downward into layer 2.
4.  The post fit and simulation module: reference points, projection,
    diagnostics, plotting, operating model.

## Tests

Test files in `tests/testthat/` carry the same prefixes as `R/`, so the
tests for a given source file sort next to each other and the name says
what is under test rather than which stock the fixture happens to use.

| Prefix | Files | What it covers |
|----|----|----|
| `test-setup_*` | 9 | Input building: the `map` factor builders and the `Setup_Mod_*` validation |
| `test-model_*` | 21 | One objective function module each: selectivity, movement, transition, observation models, likelihoods, distributions |
| `test-utils_*` | 2 | Shared numerical helpers |
| `test-sim_*` | 4 | Operating model, including simulate then refit self tests |
| `test-refpts_*` | 8 | SPR and MSY solvers, one file per spatial structure |
| `test-projection_*` | 3 | Forward projection off a fitted model |
| `test-diag_*` | 8 | Post fit diagnostics: retrospectives and OSA residuals |
| `test-integration_*` | 4 | Cross cutting agreement between the objective, the reference points and the operating model |
| `test-regression_*` | 5 | End to end fits pinning `obj$rep` and `nll` for known configurations |

That is 64 test files in total.

Two groups are worth calling out. The `test-regression_*` files pin
`obj$rep` and `nll` values for known configurations against bundled
example fits, so they are the ones most likely to catch an accidental
change in the objective function’s numerics. The
`test-integration_move_timing_*` files are the tightest constraint on
`model_transition.R` and everything built on it: they check that the
reference points, the forward projection and the operating model all
agree under each movement timing, to a tolerance that leaves no room for
a sequencing mistake.

## Extending the model

Adding a new process, fleet type or data source interacts with a
predictable set of places.

1.  **Setup.** Add or extend a `Setup_Mod_*` function to accept the new
    spec, validate it, and append the results to `input_list$data`,
    `$par` and `$map`. If parameters need fixing or sharing, add a
    `do_*_mapping` helper in the same file, or reach for the shared
    builders in `setup_mapping.R` if the pattern already exists.
2.  **Objective.** Add a `## Section` in `model_objective.R` that
    consumes the new entries and contributes to `jnll` or `nll`. Put it
    in the section it belongs to (transform, observation model, or
    likelihood and prior) rather than appending at the end. If the logic
    runs to more than a few lines, extract it into a `model_*.R` module
    instead of growing the objective inline.
3.  **Simulation.** If the process should be simulateable, add the
    matching `Setup_Sim_*` inputs and the generation logic in
    `sim_population.R` or `sim_observations.R`.
4.  **Post fit.** Expose derived output through a `Get_*` function if
    other code needs to reuse it, and through `plot_figures_tables.R` or
    `diag_fits.R` if it needs to be seen.
5.  **Tests.** Add a targeted test near the feature, and check whether
    the regression tests need updated reference values.
6.  **Vignette.** If the feature is user facing, document it in the
    relevant lettered vignette rather than only in code comments.

For the mechanics of getting a change in (branch naming, running tests
locally, opening a PR) see the
[Contributing](https://chengmatt.github.io/SPoRC/dev/articles/contributing.md)
guide.
