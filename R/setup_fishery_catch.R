# Stage 1 of 3: model setup
#
# Catch and fishing mortality inputs, plus discard mortality rate.
# Setup_Mod_Catch_and_F is the entry point and calls every do_*_mapping helper in
# this file to build the parameter maps for F, its deviations, the catch and
# discard observation error terms, and the discard mortality rate.
# Setup_Sim_Fishing is the operating model counterpart.

#' Setup Simulation Fishing Inputs
#'
#' Initializes and validates fishing-related inputs for a simulation list (`sim_list`).
#' This includes fishing mortality, selectivity, catchability, observation error,
#' and age- and length-composition parameters for both aggregate and population-specific data.
#'
#' @param sim_list A list containing simulation settings, including the number of populations
#'   (`n_pop`), regions (`n_regions`), years (`n_yrs`), seasons (`n_seas`), ages (`n_ages`),
#'   sexes (`n_sexes`), fishing fleets (`n_fish_fleets`), and simulations (`n_sims`).
#'
#' @param ln_sigmaC Numeric array. Log-scale observation SD for total catch,
#'   dimensions `n_regions x n_yrs x n_seas x n_fish_fleets`. Default: log(0.02).
#' @param ln_sigmaC_pop Numeric array. Log-scale observation SD for population-specific catch,
#'   dimensions `n_pop x n_regions x n_yrs x n_seas x n_fish_fleets`. Default: log(0.02).
#' @param catch_units Numeric vector. Catch units (0 = abundance, 1 = biomass),
#'   length `n_fish_fleets`. Default: 1.
#'
#' @param init_F_val Numeric array. Initial fishing mortality,
#'   dimensions `n_regions x n_seas x n_fish_fleets`. Default: 0.
#' @param Fmort_input Numeric array. Fishing mortality,
#'   dimensions `n_regions x n_yrs x n_seas x n_fish_fleets x n_sims`. Default: 0.1.
#' @param fish_sel_input Numeric array. Fishery selectivity,
#'   dimensions `n_pop x n_regions x n_yrs x n_seas x n_ages x n_sexes x n_fish_fleets x n_sims`.
#' @param fish_q_input Numeric array. Catchability,
#'   dimensions `n_regions x n_yrs x n_fish_fleets x n_sims`. Default: 1.
#'
#' @param ObsFishIdx_SE Numeric array. Observation SD for fishery indices,
#'   dimensions `n_regions x n_yrs x n_seas x n_fish_fleets`. Default: 0.2.
#' @param ObsFishIdx_pop_SE Numeric array. Observation SD for population-specific fishery indices,
#'   dimensions `n_pop x n_regions x n_yrs x n_seas x n_fish_fleets`. Default: 0.2.
#' @param fish_idx_type Numeric array. Index type (0 = abundance, 1 = biomass),
#'   dimensions `n_regions x n_fish_fleets`. Default: 1.
#'
#' @param comp_fishage_like Numeric vector. Likelihood for age composition
#'   (0 = Multinomial, 1 = Dirichlet-Multinomial, 2–4 = Logistic-Normal variants),
#'   length `n_fish_fleets`. Default: 0.
#' @param ISS_FishAgeComps Numeric array. Effective sample sizes for age compositions,
#'   dimensions `n_regions x n_yrs x n_seas x n_sexes x n_fish_fleets x n_sims`. Default: 100.
#' @param ln_FishAge_theta Numeric array. Log-scale overdispersion for fishery age compositions,
#'   dimensions `n_regions x n_sexes x n_fish_fleets`. Default: log(1).
#' @param ln_FishAge_theta_agg Numeric vector. Aggregated log-scale overdispersion for fishery age compositions,
#'   length `n_fish_fleets`. Default: log(1).
#' @param FishAge_corr_pars Numeric array. Correlation parameters for fishery age compositions,
#'   dimensions `n_regions x n_sexes x n_fish_fleets x 2`. Default: 0.01.
#' @param FishAge_corr_pars_agg Numeric vector. Aggregated correlation parameters for fishery age compositions,
#'   length `n_fish_fleets`. Default: 0.01.
#' @param FishAgeComps_Type Numeric array. Composition structure for fishery age compositions
#'   (0 = aggregated, 1 = split region/sex, 2 = split region joint sex, 999 = none),
#'   dimensions `n_yrs x n_fish_fleets`. Default: 2.
#'
#' @param comp_fishlen_like Numeric vector. Likelihood for length composition
#'   (0 = Multinomial, 1 = Dirichlet-Multinomial, 2–4 = Logistic-Normal variants),
#'   length `n_fish_fleets`. Default: 0.
#' @param ISS_FishLenComps Numeric array. Effective sample sizes for length compositions,
#'   dimensions `n_regions x n_yrs x n_seas x n_sexes x n_fish_fleets x n_sims`. Default: 100.
#' @param ln_FishLen_theta Numeric array. Log-scale overdispersion for fishery length compositions,
#'   dimensions `n_regions x n_sexes x n_fish_fleets`. Default: log(1).
#' @param ln_FishLen_theta_agg Numeric vector. Aggregated log-scale overdispersion for fishery length compositions,
#'   length `n_fish_fleets`. Default: log(1).
#' @param FishLen_corr_pars Numeric array. Correlation parameters for fishery length compositions,
#'   dimensions `n_regions x n_sexes x n_fish_fleets x 2`. Default: 0.01.
#' @param FishLen_corr_pars_agg Numeric vector. Aggregated correlation parameters for fishery length compositions,
#'   length `n_fish_fleets`. Default: 0.01.
#' @param FishLenComps_Type Numeric array. Composition structure for fishery length compositions
#'   (0 = aggregated, 1 = split region/sex, 2 = split region joint sex, 999 = none),
#'   dimensions `n_yrs x n_fish_fleets`. Default: 2.
#'
#' @param comp_fishage_pop_like Numeric vector. Likelihood for population-specific fishery age composition
#'   (0 = Multinomial, 1 = Dirichlet-Multinomial, 2–4 = Logistic-Normal variants),
#'   length `n_fish_fleets`. Default: 0.
#' @param ISS_FishAgeComps_pop Numeric array. Effective sample sizes for population-specific fishery age compositions,
#'   dimensions `n_pop x n_regions x n_yrs x n_seas x n_sexes x n_fish_fleets x n_sims`. Default: 100.
#' @param ln_FishAge_pop_theta Numeric array. Log-scale overdispersion for population-specific fishery age compositions,
#'   dimensions `n_pop x n_regions x n_sexes x n_fish_fleets`. Default: log(1).
#' @param ln_FishAge_pop_theta_agg Numeric array. Aggregated log-scale overdispersion for population-specific fishery age compositions,
#'   dimensions `n_pop x n_fish_fleets`. Default: log(1).
#' @param FishAge_pop_corr_pars Numeric array. Correlation parameters for population-specific fishery age compositions,
#'   dimensions `n_pop x n_regions x n_sexes x n_fish_fleets x 2`. Default: 0.01.
#' @param FishAge_pop_corr_pars_agg Numeric array. Aggregated correlation parameters for population-specific fishery age compositions,
#'   dimensions `n_pop x n_fish_fleets`. Default: 0.01.
#' @param FishAgeComps_pop_Type Numeric array. Composition structure for population-specific fishery age compositions
#'   (0 = aggregated, 1 = split region/sex, 2 = split region joint sex, 999 = none),
#'   dimensions `n_yrs x n_fish_fleets`. Default: 2.
#'
#' @param comp_fishlen_pop_like Numeric vector. Likelihood for population-specific fishery length composition
#'   (0 = Multinomial, 1 = Dirichlet-Multinomial, 2–4 = Logistic-Normal variants),
#'   length `n_fish_fleets`. Default: 0.
#' @param ISS_FishLenComps_pop Numeric array. Effective sample sizes for population-specific fishery length compositions,
#'   dimensions `n_pop x n_regions x n_yrs x n_seas x n_sexes x n_fish_fleets x n_sims`. Default: 100.
#' @param ln_FishLen_pop_theta Numeric array. Log-scale overdispersion for population-specific fishery length compositions,
#'   dimensions `n_pop x n_regions x n_sexes x n_fish_fleets`. Default: log(1).
#' @param ln_FishLen_pop_theta_agg Numeric array. Aggregated log-scale overdispersion for population-specific fishery length compositions,
#'   dimensions `n_pop x n_fish_fleets`. Default: log(1).
#' @param FishLen_pop_corr_pars Numeric array. Correlation parameters for population-specific fishery length compositions,
#'   dimensions `n_pop x n_regions x n_sexes x n_fish_fleets x 2`. Default: 0.01.
#' @param FishLen_pop_corr_pars_agg Numeric array. Aggregated correlation parameters for population-specific fishery length compositions,
#'   dimensions `n_pop x n_fish_fleets`. Default: 0.01.
#' @param FishLenComps_pop_Type Numeric array. Composition structure for population-specific fishery length compositions
#'   (0 = aggregated, 1 = split region/sex, 2 = split region joint sex, 999 = none),
#'   dimensions `n_yrs x n_fish_fleets`. Default: 2.
#'
#' @param ret_sel_input Numeric array. Retained selectivity at age,
#'   dimensions `n_pop x n_regions x n_yrs x n_seas x n_ages x n_sexes x n_fish_fleets x n_sims`.
#'   Default: 1.
#' @param dmr_input Numeric array. Discard mortality rate,
#'   dimensions `n_regions x n_yrs x n_seas x n_fish_fleets x n_sims`. Default: 0.
#' @param discard_units Numeric vector. Discard units
#'   (0 = abundance, 1 = biomass, 2 = abundance fraction, 3 = biomass fraction),
#'   length `n_fish_fleets`. Default: 3.
#' @param ln_sigmaD Numeric array. Log-scale observation SD for discards,
#'   dimensions `n_regions x n_yrs x n_seas x n_fish_fleets`. Default: log(0.02).
#' @param ln_sigmaD_pop Numeric array. Log-scale observation SD for population-specific discards,
#'   dimensions `n_pop x n_regions x n_yrs x n_seas x n_fish_fleets`. Default: log(0.02).
#'
#' @param comp_fishage_discard_like Numeric vector. Likelihood for discard age composition
#'   (0 = Multinomial, 1 = Dirichlet-Multinomial, 2–4 = Logistic-Normal variants, 999 = none),
#'   length `n_fish_fleets`. Default: 0.
#' @param ISS_FishAgeComps_discard Numeric array. Effective sample sizes for discard age compositions,
#'   dimensions `n_regions x n_yrs x n_seas x n_sexes x n_fish_fleets x n_sims`. Default: 100.
#' @param ln_FishAge_discard_theta Numeric array. Log-scale overdispersion for discard age compositions,
#'   dimensions `n_regions x n_sexes x n_fish_fleets`. Default: log(1).
#' @param ln_FishAge_discard_theta_agg Numeric vector. Aggregated log-scale overdispersion for discard age compositions,
#'   length `n_fish_fleets`. Default: log(1).
#' @param FishAge_discard_corr_pars Numeric array. Correlation parameters for discard age compositions,
#'   dimensions `n_regions x n_sexes x n_fish_fleets x 2`. Default: 0.01.
#' @param FishAge_discard_corr_pars_agg Numeric vector. Aggregated correlation parameters for discard age compositions,
#'   length `n_fish_fleets`. Default: 0.01.
#' @param FishAgeComps_discard_Type Numeric array. Composition structure for discard age compositions
#'   (0 = aggregated, 1 = split region/sex, 2 = split region joint sex, 999 = none),
#'   dimensions `n_yrs x n_fish_fleets`. Default: 2.
#'
#' @param comp_fishlen_discard_like Numeric vector. Likelihood for discard length composition
#'   (0 = Multinomial, 1 = Dirichlet-Multinomial, 2–4 = Logistic-Normal variants, 999 = none),
#'   length `n_fish_fleets`. Default: 0.
#' @param ISS_FishLenComps_discard Numeric array. Effective sample sizes for discard length compositions,
#'   dimensions `n_regions x n_yrs x n_seas x n_sexes x n_fish_fleets x n_sims`. Default: 100.
#' @param ln_FishLen_discard_theta Numeric array. Log-scale overdispersion for discard length compositions,
#'   dimensions `n_regions x n_sexes x n_fish_fleets`. Default: log(1).
#' @param ln_FishLen_discard_theta_agg Numeric vector. Aggregated log-scale overdispersion for discard length compositions,
#'   length `n_fish_fleets`. Default: log(1).
#' @param FishLen_discard_corr_pars Numeric array. Correlation parameters for discard length compositions,
#'   dimensions `n_regions x n_sexes x n_fish_fleets x 2`. Default: 0.01.
#' @param FishLen_discard_corr_pars_agg Numeric vector. Aggregated correlation parameters for discard length compositions,
#'   length `n_fish_fleets`. Default: 0.01.
#' @param FishLenComps_discard_Type Numeric array. Composition structure for discard length compositions
#'   (0 = aggregated, 1 = split region/sex, 2 = split region joint sex, 999 = none),
#'   dimensions `n_yrs x n_fish_fleets`. Default: 2.
#'
#' @param comp_fishage_discard_pop_like Numeric vector. Likelihood for population-specific discard age composition
#'   (0 = Multinomial, 1 = Dirichlet-Multinomial, 2–4 = Logistic-Normal variants, 999 = none),
#'   length `n_fish_fleets`. Default: 0.
#' @param ISS_FishAgeComps_discard_pop Numeric array. Effective sample sizes for population-specific discard age compositions,
#'   dimensions `n_pop x n_regions x n_yrs x n_seas x n_sexes x n_fish_fleets x n_sims`. Default: 100.
#' @param ln_FishAge_discard_pop_theta Numeric array. Log-scale overdispersion for population-specific discard age compositions,
#'   dimensions `n_pop x n_regions x n_sexes x n_fish_fleets`. Default: log(1).
#' @param ln_FishAge_discard_pop_theta_agg Numeric array. Aggregated log-scale overdispersion for population-specific discard age compositions,
#'   dimensions `n_pop x n_fish_fleets`. Default: log(1).
#' @param FishAge_discard_pop_corr_pars Numeric array. Correlation parameters for population-specific discard age compositions,
#'   dimensions `n_pop x n_regions x n_sexes x n_fish_fleets x 2`. Default: 0.01.
#' @param FishAge_discard_pop_corr_pars_agg Numeric array. Aggregated correlation parameters for population-specific discard age compositions,
#'   dimensions `n_pop x n_fish_fleets`. Default: 0.01.
#' @param FishAgeComps_discard_pop_Type Numeric array. Composition structure for population-specific discard age compositions
#'   (0 = aggregated, 1 = split region/sex, 2 = split region joint sex, 999 = none),
#'   dimensions `n_yrs x n_fish_fleets`. Default: 2.
#'
#' @param comp_fishlen_discard_pop_like Numeric vector. Likelihood for population-specific discard length composition
#'   (0 = Multinomial, 1 = Dirichlet-Multinomial, 2–4 = Logistic-Normal variants, 999 = none),
#'   length `n_fish_fleets`. Default: 0.
#' @param ISS_FishLenComps_discard_pop Numeric array. Effective sample sizes for population-specific discard length compositions,
#'   dimensions `n_pop x n_regions x n_yrs x n_seas x n_sexes x n_fish_fleets x n_sims`. Default: 100.
#' @param ln_FishLen_discard_pop_theta Numeric array. Log-scale overdispersion for population-specific discard length compositions,
#'   dimensions `n_pop x n_regions x n_sexes x n_fish_fleets`. Default: log(1).
#' @param ln_FishLen_discard_pop_theta_agg Numeric array. Aggregated log-scale overdispersion for population-specific discard length compositions,
#'   dimensions `n_pop x n_fish_fleets`. Default: log(1).
#' @param FishLen_discard_pop_corr_pars Numeric array. Correlation parameters for population-specific discard length compositions,
#'   dimensions `n_pop x n_regions x n_sexes x n_fish_fleets x 2`. Default: 0.01.
#' @param FishLen_discard_pop_corr_pars_agg Numeric array. Aggregated correlation parameters for population-specific discard length compositions,
#'   dimensions `n_pop x n_fish_fleets`. Default: 0.01.
#' @param FishLenComps_discard_pop_Type Numeric array. Composition structure for population-specific discard length compositions
#'   (0 = aggregated, 1 = split region/sex, 2 = split region joint sex, 999 = none),
#'   dimensions `n_yrs x n_fish_fleets`. Default: 2.
#'
#' @return A modified `sim_list` with validated fishing-related inputs.
#'
#' @export Setup_Sim_Fishing
#' @family Simulation Setup
Setup_Sim_Fishing <- function(sim_list,

                              # Retained / total fishery dynamics
                              ln_sigmaC = array(log(0.02), dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets)),
                              ln_sigmaC_pop = array(log(0.02), dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets)),
                              catch_units = array(1, dim = c(sim_list$n_fish_fleets)),
                              init_F_val = array(0, dim = c(sim_list$n_regions, sim_list$n_seas, sim_list$n_fish_fleets)),
                              Fmort_input = array(0.1, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets, sim_list$n_sims)),
                              fish_sel_input,
                              fish_q_input = array(1, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_fish_fleets, sim_list$n_sims)),
                              ObsFishIdx_SE = array(0.2, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets)),
                              ObsFishIdx_pop_SE = array(0.2, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets)),
                              fish_idx_type = array(1, dim = c(sim_list$n_regions, sim_list$n_fish_fleets)),

                              # Retained age compositions
                              comp_fishage_like = rep(0, sim_list$n_fish_fleets),
                              ISS_FishAgeComps = array(100, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims)),
                              ln_FishAge_theta = array(log(1), dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets)),
                              ln_FishAge_theta_agg = rep(log(1), sim_list$n_fish_fleets),
                              FishAge_corr_pars_agg = rep(0.01, sim_list$n_fish_fleets),
                              FishAge_corr_pars = array(0.01, dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets, 2)),
                              FishAgeComps_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_fish_fleets)),

                              # Retained length compositions
                              comp_fishlen_like = rep(0, sim_list$n_fish_fleets),
                              ISS_FishLenComps = array(100, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims)),
                              ln_FishLen_theta = array(log(1), dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets)),
                              ln_FishLen_theta_agg = rep(log(1), sim_list$n_fish_fleets),
                              FishLen_corr_pars_agg = rep(0.01, sim_list$n_fish_fleets),
                              FishLen_corr_pars = array(0.01, dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets, 2)),
                              FishLenComps_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_fish_fleets)),

                              # Retained age compositions (population-specific)
                              comp_fishage_pop_like = rep(0, sim_list$n_fish_fleets),
                              ISS_FishAgeComps_pop = array(100, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims)),
                              ln_FishAge_pop_theta = array(log(1), dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets)),
                              ln_FishAge_pop_theta_agg = array(log(1), dim = c(sim_list$n_pop, sim_list$n_fish_fleets)),
                              FishAge_pop_corr_pars = array(0.01, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets, 2)),
                              FishAge_pop_corr_pars_agg = array(0.01, dim = c(sim_list$n_pop, sim_list$n_fish_fleets)),
                              FishAgeComps_pop_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_fish_fleets)),

                              # Retained length compositions (population-specific)
                              comp_fishlen_pop_like = rep(0, sim_list$n_fish_fleets),
                              ISS_FishLenComps_pop = array(100, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims)),
                              ln_FishLen_pop_theta = array(log(1), dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets)),
                              ln_FishLen_pop_theta_agg = array(log(1), dim = c(sim_list$n_pop, sim_list$n_fish_fleets)),
                              FishLen_pop_corr_pars = array(0.01, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets, 2)),
                              FishLen_pop_corr_pars_agg = array(0.01, dim = c(sim_list$n_pop, sim_list$n_fish_fleets)),
                              FishLenComps_pop_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_fish_fleets)),

                              # Retention and discards
                              ret_sel_input = array(1, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims)),
                              dmr_input = array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets, sim_list$n_sims)),
                              discard_units = array(3, dim = c(sim_list$n_fish_fleets)),
                              ln_sigmaD = array(log(0.02), dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets)),
                              ln_sigmaD_pop = array(log(0.02), dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets)),

                              # Discard age compositions (non-population specific)
                              comp_fishage_discard_like = rep(0, sim_list$n_fish_fleets),
                              ISS_FishAgeComps_discard = array(100, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims)),
                              ln_FishAge_discard_theta = array(log(1), dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets)),
                              ln_FishAge_discard_theta_agg = rep(log(1), sim_list$n_fish_fleets),
                              FishAge_discard_corr_pars = array(0.01, dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets, 2)),
                              FishAge_discard_corr_pars_agg = rep(0.01, sim_list$n_fish_fleets),
                              FishAgeComps_discard_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_fish_fleets)),

                              # Discard length compositions (non-population specific)
                              comp_fishlen_discard_like = rep(0, sim_list$n_fish_fleets),
                              ISS_FishLenComps_discard = array(100, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims)),
                              ln_FishLen_discard_theta = array(log(1), dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets)),
                              ln_FishLen_discard_theta_agg = rep(log(1), sim_list$n_fish_fleets),
                              FishLen_discard_corr_pars = array(0.01, dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets, 2)),
                              FishLen_discard_corr_pars_agg = rep(0.01, sim_list$n_fish_fleets),
                              FishLenComps_discard_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_fish_fleets)),

                              # Discard age compositions (population specific)
                              comp_fishage_discard_pop_like = rep(0, sim_list$n_fish_fleets),
                              ISS_FishAgeComps_discard_pop = array(100, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims)),
                              ln_FishAge_discard_pop_theta = array(log(1), dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets)),
                              ln_FishAge_discard_pop_theta_agg = array(log(1), dim = c(sim_list$n_pop, sim_list$n_fish_fleets)),
                              FishAge_discard_pop_corr_pars = array(0.01, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets, 2)),
                              FishAge_discard_pop_corr_pars_agg = array(0.01, dim = c(sim_list$n_pop, sim_list$n_fish_fleets)),
                              FishAgeComps_discard_pop_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_fish_fleets)),

                              # Discard length compositions (population specific)
                              comp_fishlen_discard_pop_like = rep(0, sim_list$n_fish_fleets),
                              ISS_FishLenComps_discard_pop = array(100, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims)),
                              ln_FishLen_discard_pop_theta = array(log(1), dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets)),
                              ln_FishLen_discard_pop_theta_agg = array(log(1), dim = c(sim_list$n_pop, sim_list$n_fish_fleets)),
                              FishLen_discard_pop_corr_pars = array(0.01, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets, 2)),
                              FishLen_discard_pop_corr_pars_agg = array(0.01, dim = c(sim_list$n_pop, sim_list$n_fish_fleets)),
                              FishLenComps_discard_pop_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_fish_fleets))

                              ) {

  # Convert character inputs to numeric codes
  catch_units <- convert_to_numeric(catch_units,  list(abd = 0, biom = 1))
  fish_idx_type <- convert_to_numeric(fish_idx_type, list(abd = 0, biom = 1))
  comp_fishage_like <- convert_to_numeric(comp_fishage_like, list(Multinomial = 0,  `Dirichlet-Multinomial` = 1, `iid-Logistic-Normal` = 2, `1d-Logistic-Normal` = 3, `2d-Logistic-Normal` = 4))
  comp_fishlen_like <- convert_to_numeric(comp_fishlen_like, list(Multinomial = 0, `Dirichlet-Multinomial` = 1, `iid-Logistic-Normal` = 2, `1d-Logistic-Normal` = 3, `2d-Logistic-Normal` = 4))
  FishAgeComps_Type <- convert_to_numeric(FishAgeComps_Type,  list(agg = 0, spltRspltS = 1, spltRjntS = 2, none = 999))
  FishLenComps_Type <- convert_to_numeric(FishLenComps_Type,  list(agg = 0, spltRspltS = 1, spltRjntS = 2, none = 999))
  FishAgeComps_pop_Type <- convert_to_numeric(FishAgeComps_pop_Type,  list(agg = 0, spltRspltS = 1, spltRjntS = 2, none = 999))
  FishLenComps_pop_Type <- convert_to_numeric(FishLenComps_pop_Type,  list(agg = 0, spltRspltS = 1, spltRjntS = 2, none = 999))
  comp_fishage_discard_like <- convert_to_numeric(comp_fishage_discard_like, list(Multinomial = 0, `Dirichlet-Multinomial` = 1, `iid-Logistic-Normal` = 2, `1d-Logistic-Normal` = 3, `2d-Logistic-Normal` = 4, none = 999))
  comp_fishlen_discard_like <- convert_to_numeric(comp_fishlen_discard_like, list(Multinomial = 0, `Dirichlet-Multinomial` = 1, `iid-Logistic-Normal` = 2, `1d-Logistic-Normal` = 3, `2d-Logistic-Normal` = 4, none = 999))
  FishAgeComps_discard_Type <- convert_to_numeric(FishAgeComps_discard_Type, list(agg = 0, spltRspltS = 1, spltRjntS = 2, none = 999))
  FishLenComps_discard_Type <- convert_to_numeric(FishLenComps_discard_Type, list(agg = 0, spltRspltS = 1, spltRjntS = 2, none = 999))
  comp_fishage_discard_pop_like <- convert_to_numeric(comp_fishage_discard_pop_like, list(Multinomial = 0, `Dirichlet-Multinomial` = 1, `iid-Logistic-Normal` = 2, `1d-Logistic-Normal` = 3, `2d-Logistic-Normal` = 4, none = 999))
  comp_fishlen_discard_pop_like <- convert_to_numeric(comp_fishlen_discard_pop_like, list(Multinomial = 0, `Dirichlet-Multinomial` = 1, `iid-Logistic-Normal` = 2, `1d-Logistic-Normal` = 3, `2d-Logistic-Normal` = 4, none = 999))
  FishAgeComps_discard_pop_Type <- convert_to_numeric(FishAgeComps_discard_pop_Type, list(agg = 0, spltRspltS = 1, spltRjntS = 2, none = 999))
  FishLenComps_discard_pop_Type <- convert_to_numeric(FishLenComps_discard_pop_Type, list(agg = 0, spltRspltS = 1, spltRjntS = 2, none = 999))
  discard_units <- convert_to_numeric(discard_units, list(abd = 0, biom = 1, abd_frac = 2, biom_frac = 3))

  # Validate dimensions of all input parameters
  check_sim_dimensions(ln_sigmaC, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs, n_seas = sim_list$n_seas,
                       n_fish_fleets = sim_list$n_fish_fleets, what = "ln_sigmaC")
  check_sim_dimensions(ln_sigmaC_pop, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs, n_seas = sim_list$n_seas,
                       n_fish_fleets = sim_list$n_fish_fleets, n_pop = sim_list$n_pop, what = "ln_sigmaC_pop")
  check_sim_dimensions(Fmort_input, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs, n_seas = sim_list$n_seas,
                       n_fish_fleets = sim_list$n_fish_fleets, n_sims = sim_list$n_sims, what = "Fmort_input")
  check_sim_dimensions(fish_sel_input, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs,
                       n_ages = sim_list$n_ages, n_sexes = sim_list$n_sexes, n_pop = sim_list$n_pop, n_seas = sim_list$n_seas,
                       n_fish_fleets = sim_list$n_fish_fleets, n_sims = sim_list$n_sims, what = "fish_sel_input")
  check_sim_dimensions(fish_q_input, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs,
                       n_fish_fleets = sim_list$n_fish_fleets, n_sims = sim_list$n_sims, what = "fish_q_input")
  check_sim_dimensions(ObsFishIdx_SE, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs, n_seas = sim_list$n_seas,
                       n_fish_fleets = sim_list$n_fish_fleets, what = "ObsFishIdx_SE")
  check_sim_dimensions(ObsFishIdx_pop_SE, n_pop = sim_list$n_pop, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs, n_seas = sim_list$n_seas,
                       n_fish_fleets = sim_list$n_fish_fleets, what = "ObsFishIdx_pop_SE")

  # Validate fishery age composition parameters
  check_sim_dimensions(comp_fishage_like, n_fish_fleets = sim_list$n_fish_fleets, what = "comp_fishage_like")
  check_sim_dimensions(ISS_FishAgeComps, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs, n_seas = sim_list$n_seas,
                       n_sexes = sim_list$n_sexes, n_fish_fleets = sim_list$n_fish_fleets,
                       n_sims = sim_list$n_sims, what = "ISS_FishAgeComps")
  check_sim_dimensions(ln_FishAge_theta, n_regions = sim_list$n_regions, n_sexes = sim_list$n_sexes,
                       n_fish_fleets = sim_list$n_fish_fleets, what = "ln_FishAge_theta")
  check_sim_dimensions(ln_FishAge_theta_agg, n_fish_fleets = sim_list$n_fish_fleets, what = "ln_FishAge_theta_agg")
  check_sim_dimensions(FishAge_corr_pars_agg, n_fish_fleets = sim_list$n_fish_fleets, what = "FishAge_corr_pars_agg")
  check_sim_dimensions(FishAge_corr_pars, n_regions = sim_list$n_regions, n_sexes = sim_list$n_sexes,
                       n_fish_fleets = sim_list$n_fish_fleets, what = "FishAge_corr_pars")
  check_sim_dimensions(FishAgeComps_Type, n_years = sim_list$n_yrs, n_fish_fleets = sim_list$n_fish_fleets,
                       what = "FishAgeComps_Type")
  check_sim_dimensions(comp_fishage_pop_like, n_fish_fleets = sim_list$n_fish_fleets, what = "comp_fishage_pop_like")
  check_sim_dimensions(ISS_FishAgeComps_pop, n_pop = sim_list$n_pop, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs, n_seas = sim_list$n_seas, n_sexes = sim_list$n_sexes, n_fish_fleets = sim_list$n_fish_fleets, n_sims = sim_list$n_sims, what = "ISS_FishAgeComps_pop")
  check_sim_dimensions(ln_FishAge_pop_theta, n_pop = sim_list$n_pop, n_regions = sim_list$n_regions, n_sexes = sim_list$n_sexes, n_fish_fleets = sim_list$n_fish_fleets, what = "ln_FishAge_pop_theta")
  check_sim_dimensions(ln_FishAge_pop_theta_agg, n_pop = sim_list$n_pop, n_fish_fleets = sim_list$n_fish_fleets, what = "ln_FishAge_pop_theta_agg")
  check_sim_dimensions(FishAge_pop_corr_pars_agg, n_pop = sim_list$n_pop, n_fish_fleets = sim_list$n_fish_fleets, what = "FishAge_pop_corr_pars_agg")
  check_sim_dimensions(FishAge_pop_corr_pars, n_pop = sim_list$n_pop, n_regions = sim_list$n_regions, n_sexes = sim_list$n_sexes, n_fish_fleets = sim_list$n_fish_fleets, what = "FishAge_pop_corr_pars")
  check_sim_dimensions(FishAgeComps_pop_Type, n_years = sim_list$n_yrs, n_fish_fleets = sim_list$n_fish_fleets, what = "FishAgeComps_pop_Type")


  # Validate fishery length composition parameters
  check_sim_dimensions(comp_fishlen_like, n_fish_fleets = sim_list$n_fish_fleets, what = "comp_fishlen_like")
  check_sim_dimensions(ISS_FishLenComps, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs, n_seas = sim_list$n_seas,
                       n_sexes = sim_list$n_sexes, n_fish_fleets = sim_list$n_fish_fleets,
                       n_sims = sim_list$n_sims, what = "ISS_FishLenComps")
  check_sim_dimensions(ln_FishLen_theta, n_regions = sim_list$n_regions, n_sexes = sim_list$n_sexes,
                       n_fish_fleets = sim_list$n_fish_fleets, what = "ln_FishLen_theta")
  check_sim_dimensions(ln_FishLen_theta_agg, n_fish_fleets = sim_list$n_fish_fleets, what = "ln_FishLen_theta_agg")
  check_sim_dimensions(FishLen_corr_pars_agg, n_fish_fleets = sim_list$n_fish_fleets, what = "FishLen_corr_pars_agg")
  check_sim_dimensions(FishLen_corr_pars, n_regions = sim_list$n_regions, n_sexes = sim_list$n_sexes,
                       n_fish_fleets = sim_list$n_fish_fleets, what = "FishLen_corr_pars")
  check_sim_dimensions(FishLenComps_Type, n_years = sim_list$n_yrs, n_fish_fleets = sim_list$n_fish_fleets,
                       what = "FishLenComps_Type")
  check_sim_dimensions(comp_fishlen_pop_like, n_fish_fleets = sim_list$n_fish_fleets, what = "comp_fishlen_pop_like")
  check_sim_dimensions(ISS_FishLenComps_pop, n_pop = sim_list$n_pop, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs, n_seas = sim_list$n_seas, n_sexes = sim_list$n_sexes, n_fish_fleets = sim_list$n_fish_fleets, n_sims = sim_list$n_sims, what = "ISS_FishLenComps_pop")
  check_sim_dimensions(ln_FishLen_pop_theta, n_pop = sim_list$n_pop, n_regions = sim_list$n_regions, n_sexes = sim_list$n_sexes, n_fish_fleets = sim_list$n_fish_fleets, what = "ln_FishLen_pop_theta")
  check_sim_dimensions(ln_FishLen_pop_theta_agg, n_pop = sim_list$n_pop, n_fish_fleets = sim_list$n_fish_fleets, what = "ln_FishLen_pop_theta_agg")
  check_sim_dimensions(FishLen_pop_corr_pars_agg, n_pop = sim_list$n_pop, n_fish_fleets = sim_list$n_fish_fleets, what = "FishLen_pop_corr_pars_agg")
  check_sim_dimensions(FishLen_pop_corr_pars, n_pop = sim_list$n_pop, n_regions = sim_list$n_regions, n_sexes = sim_list$n_sexes, n_fish_fleets = sim_list$n_fish_fleets, what = "FishLen_pop_corr_pars")
  check_sim_dimensions(FishLenComps_pop_Type, n_years = sim_list$n_yrs, n_fish_fleets = sim_list$n_fish_fleets, what = "FishLenComps_pop_Type")


  # Validate retention and discard inputs
  check_sim_dimensions(ret_sel_input, n_pop = sim_list$n_pop, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs,
                       n_seas = sim_list$n_seas, n_ages = sim_list$n_ages, n_sexes = sim_list$n_sexes,
                       n_fish_fleets = sim_list$n_fish_fleets, n_sims = sim_list$n_sims, what = "ret_sel_input")
  check_sim_dimensions(dmr_input, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs, n_seas = sim_list$n_seas,
                       n_fish_fleets = sim_list$n_fish_fleets, n_sims = sim_list$n_sims, what = "dmr_input")
  check_sim_dimensions(ln_sigmaD, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs, n_seas = sim_list$n_seas,
                       n_fish_fleets = sim_list$n_fish_fleets, what = "ln_sigmaD")
  check_sim_dimensions(ln_sigmaD_pop, n_pop = sim_list$n_pop, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs,
                       n_seas = sim_list$n_seas, n_fish_fleets = sim_list$n_fish_fleets, what = "ln_sigmaD_pop")

  # Validate discard age composition parameters
  check_sim_dimensions(comp_fishage_discard_like, n_fish_fleets = sim_list$n_fish_fleets, what = "comp_fishage_discard_like")
  check_sim_dimensions(ISS_FishAgeComps_discard, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs, n_seas = sim_list$n_seas,
                       n_sexes = sim_list$n_sexes, n_fish_fleets = sim_list$n_fish_fleets, n_sims = sim_list$n_sims, what = "ISS_FishAgeComps_discard")
  check_sim_dimensions(ln_FishAge_discard_theta, n_regions = sim_list$n_regions, n_sexes = sim_list$n_sexes,
                       n_fish_fleets = sim_list$n_fish_fleets, what = "ln_FishAge_discard_theta")
  check_sim_dimensions(ln_FishAge_discard_theta_agg, n_fish_fleets = sim_list$n_fish_fleets, what = "ln_FishAge_discard_theta_agg")
  check_sim_dimensions(FishAge_discard_corr_pars_agg, n_fish_fleets = sim_list$n_fish_fleets, what = "FishAge_discard_corr_pars_agg")
  check_sim_dimensions(FishAge_discard_corr_pars, n_regions = sim_list$n_regions, n_sexes = sim_list$n_sexes,
                       n_fish_fleets = sim_list$n_fish_fleets, what = "FishAge_discard_corr_pars")
  check_sim_dimensions(FishAgeComps_discard_Type, n_years = sim_list$n_yrs, n_fish_fleets = sim_list$n_fish_fleets, what = "FishAgeComps_discard_Type")

  # Validate discard length composition parameters
  check_sim_dimensions(comp_fishlen_discard_like, n_fish_fleets = sim_list$n_fish_fleets, what = "comp_fishlen_discard_like")
  check_sim_dimensions(ISS_FishLenComps_discard, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs, n_seas = sim_list$n_seas,
                       n_sexes = sim_list$n_sexes, n_fish_fleets = sim_list$n_fish_fleets, n_sims = sim_list$n_sims, what = "ISS_FishLenComps_discard")
  check_sim_dimensions(ln_FishLen_discard_theta, n_regions = sim_list$n_regions, n_sexes = sim_list$n_sexes,
                       n_fish_fleets = sim_list$n_fish_fleets, what = "ln_FishLen_discard_theta")
  check_sim_dimensions(ln_FishLen_discard_theta_agg, n_fish_fleets = sim_list$n_fish_fleets, what = "ln_FishLen_discard_theta_agg")
  check_sim_dimensions(FishLen_discard_corr_pars_agg, n_fish_fleets = sim_list$n_fish_fleets, what = "FishLen_discard_corr_pars_agg")
  check_sim_dimensions(FishLen_discard_corr_pars, n_regions = sim_list$n_regions, n_sexes = sim_list$n_sexes,
                       n_fish_fleets = sim_list$n_fish_fleets, what = "FishLen_discard_corr_pars")
  check_sim_dimensions(FishLenComps_discard_Type, n_years = sim_list$n_yrs, n_fish_fleets = sim_list$n_fish_fleets, what = "FishLenComps_discard_Type")

  # Validate population-specific discard age composition parameters
  check_sim_dimensions(comp_fishage_discard_pop_like, n_fish_fleets = sim_list$n_fish_fleets, what = "comp_fishage_discard_pop_like")
  check_sim_dimensions(ISS_FishAgeComps_discard_pop, n_pop = sim_list$n_pop, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs,
                       n_seas = sim_list$n_seas, n_sexes = sim_list$n_sexes, n_fish_fleets = sim_list$n_fish_fleets,
                       n_sims = sim_list$n_sims, what = "ISS_FishAgeComps_discard_pop")
  check_sim_dimensions(ln_FishAge_discard_pop_theta, n_pop = sim_list$n_pop, n_regions = sim_list$n_regions,
                       n_sexes = sim_list$n_sexes, n_fish_fleets = sim_list$n_fish_fleets, what = "ln_FishAge_discard_pop_theta")
  check_sim_dimensions(ln_FishAge_discard_pop_theta_agg, n_pop = sim_list$n_pop, n_fish_fleets = sim_list$n_fish_fleets, what = "ln_FishAge_discard_pop_theta_agg")
  check_sim_dimensions(FishAge_discard_pop_corr_pars_agg, n_pop = sim_list$n_pop, n_fish_fleets = sim_list$n_fish_fleets, what = "FishAge_discard_pop_corr_pars_agg")
  check_sim_dimensions(FishAge_discard_pop_corr_pars, n_pop = sim_list$n_pop, n_regions = sim_list$n_regions,
                       n_sexes = sim_list$n_sexes, n_fish_fleets = sim_list$n_fish_fleets, what = "FishAge_discard_pop_corr_pars")
  check_sim_dimensions(FishAgeComps_discard_pop_Type, n_years = sim_list$n_yrs, n_fish_fleets = sim_list$n_fish_fleets, what = "FishAgeComps_discard_pop_Type")

  # Validate population-specific discard length composition parameters
  check_sim_dimensions(comp_fishlen_discard_pop_like, n_fish_fleets = sim_list$n_fish_fleets, what = "comp_fishlen_discard_pop_like")
  check_sim_dimensions(ISS_FishLenComps_discard_pop, n_pop = sim_list$n_pop, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs,
                       n_seas = sim_list$n_seas, n_sexes = sim_list$n_sexes, n_fish_fleets = sim_list$n_fish_fleets,
                       n_sims = sim_list$n_sims, what = "ISS_FishLenComps_discard_pop")
  check_sim_dimensions(ln_FishLen_discard_pop_theta, n_pop = sim_list$n_pop, n_regions = sim_list$n_regions,
                       n_sexes = sim_list$n_sexes, n_fish_fleets = sim_list$n_fish_fleets, what = "ln_FishLen_discard_pop_theta")
  check_sim_dimensions(ln_FishLen_discard_pop_theta_agg, n_pop = sim_list$n_pop, n_fish_fleets = sim_list$n_fish_fleets, what = "ln_FishLen_discard_pop_theta_agg")
  check_sim_dimensions(FishLen_discard_pop_corr_pars_agg, n_pop = sim_list$n_pop, n_fish_fleets = sim_list$n_fish_fleets, what = "FishLen_discard_pop_corr_pars_agg")
  check_sim_dimensions(FishLen_discard_pop_corr_pars, n_pop = sim_list$n_pop, n_regions = sim_list$n_regions,
                       n_sexes = sim_list$n_sexes, n_fish_fleets = sim_list$n_fish_fleets, what = "FishLen_discard_pop_corr_pars")
  check_sim_dimensions(FishLenComps_discard_pop_Type, n_years = sim_list$n_yrs, n_fish_fleets = sim_list$n_fish_fleets, what = "FishLenComps_discard_pop_Type")

  # output variables into list
  sim_list$Fmort <- Fmort_input # input fishing mortality pattern
  sim_list$catch_units <- catch_units # catch units
  sim_list$ln_sigmaC <- ln_sigmaC # Observation sd for catch
  sim_list$ln_sigmaC_pop <- ln_sigmaC_pop # observation sd for pop-specific catch
  sim_list$init_F <- init_F_val # initial F value
  sim_list$fish_sel <- fish_sel_input # fishery selectivity
  sim_list$fish_q <- fish_q_input # fishery catchability
  sim_list$ObsFishIdx_SE <- ObsFishIdx_SE # fishery index SE
  sim_list$ObsFishIdx_pop_SE <- ObsFishIdx_pop_SE # fishery index SE pop-specific
  sim_list$fish_idx_type <- fish_idx_type # fishery index type

  # Fishery age compositions
  sim_list$comp_fishage_like <- comp_fishage_like
  sim_list$ISS_FishAgeComps <- ISS_FishAgeComps
  sim_list$ln_FishAge_theta <- ln_FishAge_theta
  sim_list$ln_FishAge_theta_agg <- ln_FishAge_theta_agg
  sim_list$FishAge_corr_pars_agg <- FishAge_corr_pars_agg
  sim_list$FishAge_corr_pars <- FishAge_corr_pars
  sim_list$FishAgeComps_Type <- FishAgeComps_Type

  # Fishery length compositions
  sim_list$comp_fishlen_like <- comp_fishlen_like
  sim_list$ISS_FishLenComps <- ISS_FishLenComps
  sim_list$ln_FishLen_theta <- ln_FishLen_theta
  sim_list$ln_FishLen_theta_agg <- ln_FishLen_theta_agg
  sim_list$FishLen_corr_pars_agg <- FishLen_corr_pars_agg
  sim_list$FishLen_corr_pars <- FishLen_corr_pars
  sim_list$FishLenComps_Type <- FishLenComps_Type

  # Population-specific stuff
  sim_list$comp_fishage_pop_like <- comp_fishage_pop_like
  sim_list$ISS_FishAgeComps_pop <- ISS_FishAgeComps_pop
  sim_list$ln_FishAge_pop_theta <- ln_FishAge_pop_theta
  sim_list$ln_FishAge_pop_theta_agg <- ln_FishAge_pop_theta_agg
  sim_list$FishAge_pop_corr_pars <- FishAge_pop_corr_pars
  sim_list$FishAge_pop_corr_pars_agg <- FishAge_pop_corr_pars_agg
  sim_list$FishAgeComps_pop_Type <- FishAgeComps_pop_Type

  sim_list$comp_fishlen_pop_like <- comp_fishlen_pop_like
  sim_list$ISS_FishLenComps_pop <- ISS_FishLenComps_pop
  sim_list$ln_FishLen_pop_theta <- ln_FishLen_pop_theta
  sim_list$ln_FishLen_pop_theta_agg <- ln_FishLen_pop_theta_agg
  sim_list$FishLen_pop_corr_pars <- FishLen_pop_corr_pars
  sim_list$FishLen_pop_corr_pars_agg <- FishLen_pop_corr_pars_agg
  sim_list$FishLenComps_pop_Type <- FishLenComps_pop_Type

  # Retention and discards
  sim_list$ret_sel <- ret_sel_input
  sim_list$dmr <- dmr_input
  sim_list$discard_units <- discard_units
  sim_list$ln_sigmaD <- ln_sigmaD
  sim_list$ln_sigmaD_pop <- ln_sigmaD_pop

  # Discard age compositions
  sim_list$comp_fishage_discard_like <- comp_fishage_discard_like
  sim_list$ISS_FishAgeComps_discard <- ISS_FishAgeComps_discard
  sim_list$ln_FishAge_discard_theta <- ln_FishAge_discard_theta
  sim_list$ln_FishAge_discard_theta_agg <- ln_FishAge_discard_theta_agg
  sim_list$FishAge_discard_corr_pars <- FishAge_discard_corr_pars
  sim_list$FishAge_discard_corr_pars_agg <- FishAge_discard_corr_pars_agg
  sim_list$FishAgeComps_discard_Type <- FishAgeComps_discard_Type

  # Discard length compositions
  sim_list$comp_fishlen_discard_like <- comp_fishlen_discard_like
  sim_list$ISS_FishLenComps_discard <- ISS_FishLenComps_discard
  sim_list$ln_FishLen_discard_theta <- ln_FishLen_discard_theta
  sim_list$ln_FishLen_discard_theta_agg <- ln_FishLen_discard_theta_agg
  sim_list$FishLen_discard_corr_pars <- FishLen_discard_corr_pars
  sim_list$FishLen_discard_corr_pars_agg <- FishLen_discard_corr_pars_agg
  sim_list$FishLenComps_discard_Type <- FishLenComps_discard_Type

  # Discard age compositions (population specific)
  sim_list$comp_fishage_discard_pop_like <- comp_fishage_discard_pop_like
  sim_list$ISS_FishAgeComps_discard_pop <- ISS_FishAgeComps_discard_pop
  sim_list$ln_FishAge_discard_pop_theta <- ln_FishAge_discard_pop_theta
  sim_list$ln_FishAge_discard_pop_theta_agg <- ln_FishAge_discard_pop_theta_agg
  sim_list$FishAge_discard_pop_corr_pars <- FishAge_discard_pop_corr_pars
  sim_list$FishAge_discard_pop_corr_pars_agg <- FishAge_discard_pop_corr_pars_agg
  sim_list$FishAgeComps_discard_pop_Type <- FishAgeComps_discard_pop_Type

  # Discard length compositions (population specific)
  sim_list$comp_fishlen_discard_pop_like <- comp_fishlen_discard_pop_like
  sim_list$ISS_FishLenComps_discard_pop <- ISS_FishLenComps_discard_pop
  sim_list$ln_FishLen_discard_pop_theta <- ln_FishLen_discard_pop_theta
  sim_list$ln_FishLen_discard_pop_theta_agg <- ln_FishLen_discard_pop_theta_agg
  sim_list$FishLen_discard_pop_corr_pars <- FishLen_discard_pop_corr_pars
  sim_list$FishLen_discard_pop_corr_pars_agg <- FishLen_discard_pop_corr_pars_agg
  sim_list$FishLenComps_discard_pop_Type <- FishLenComps_discard_pop_Type

  return(sim_list)
}


#' Map sigmaF (fishing mortality process error SD) parameters
#'
#' Constructs the \code{ln_sigmaF} factor map used by the TMB/RTMB objective
#' function to share or fix the log-scale standard deviation of fishing mortality
#' process error across regions, seasons, and fleets. All cells within a shared
#' group are assigned the same estimation index.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists, as constructed by upstream setup functions.
#' @param sigmaF_spec Character string controlling the sharing and estimation
#'   structure for \code{ln_sigmaF}. One of:
#'   \describe{
#'     \item{\code{"est_all"}}{Unique parameter per region × season × fleet combination.}
#'     \item{\code{"est_shared_r"}}{Shared across regions; unique per season × fleet.}
#'     \item{\code{"est_shared_seas"}}{Shared across seasons; unique per region × fleet.}
#'     \item{\code{"est_shared_f"}}{Shared across fleets; unique per region × season.}
#'     \item{\code{"est_shared_r_seas"}}{Shared across regions and seasons; unique per fleet.}
#'     \item{\code{"est_shared_r_f"}}{Shared across regions and fleets; unique per season.}
#'     \item{\code{"est_shared_seas_f"}}{Shared across seasons and fleets; unique per region.}
#'     \item{\code{"est_shared_r_seas_f"}}{Single parameter shared across all dimensions.}
#'     \item{\code{"fix"}}{All \code{ln_sigmaF} parameters fixed at starting values (mapped to \code{NA}).}
#'   }
#'
#' @return The input \code{input_list} with \code{$map$ln_sigmaF} set to a factor
#'   vector of length \code{prod(dim(par$ln_sigmaF))}. Each element is an integer
#'   estimation index for shared or estimated configurations, or \code{NA} when
#'   \code{sigmaF_spec = "fix"}.
#'
#' @keywords internal
do_sigmaF_mapping <- function(input_list, sigmaF_spec) {

  # Sigma F -----------------------------------------------------------------
  dims <- c(region = input_list$data$n_regions,
            season = input_list$data$n_seas,
            fleet  = input_list$data$n_fish_fleets)

  input_list$map$ln_sigmaF <- build_shared_spec_map(
    dims = dims, spec = sigmaF_spec,
    dim_abbrev = c(r = "region", seas = "season", f = "fleet")
  )

  # Print Message
  collect_message("sigmaF is specified as: ", sigmaF_spec)

  return(input_list)
}

#' Map AR1 correlation parameter for fishing mortality deviations
#'
#' Constructs the \code{Fdev_rho} factor map. \code{Fdev_rho} is only
#' meaningful when \code{Fdev_model = "ar1"} (see
#' \code{\link{Setup_Mod_Catch_and_F}}); for any other \code{Fdev_model}, all
#' \code{Fdev_rho} parameters are mapped to \code{NA} regardless of
#' \code{Fdev_rho_spec}, since they are unused by \code{\link{Get_Fdev_PE_loglik}}.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists, as constructed by upstream setup functions.
#' @param Fdev_rho_spec Character string controlling the sharing and
#'   estimation structure for \code{Fdev_rho}, following the same convention
#'   as \code{\link{do_sigmaF_mapping}}'s \code{sigmaF_spec}: one of
#'   \code{"est_all"}, \code{"est_shared_r"}, \code{"est_shared_seas"},
#'   \code{"est_shared_f"}, \code{"est_shared_r_seas"}, \code{"est_shared_r_f"},
#'   \code{"est_shared_seas_f"}, \code{"est_shared_r_seas_f"}, or \code{"fix"}.
#'
#' @return The input \code{input_list} with \code{$map$Fdev_rho} set to a
#'   factor vector of length \code{prod(dim(par$Fdev_rho))}.
#'
#' @keywords internal
do_Fdev_rho_mapping <- function(input_list, Fdev_rho_spec) {

  dims <- c(region = input_list$data$n_regions,
            season = input_list$data$n_seas,
            fleet  = input_list$data$n_fish_fleets)

  if(input_list$data$Fdev_model != 3) { # only AR1 uses Fdev_rho
    input_list$map$Fdev_rho <- factor(rep(NA, prod(dims)))
  } else {
    input_list$map$Fdev_rho <- build_shared_spec_map(
      dims = dims, spec = Fdev_rho_spec,
      dim_abbrev = c(r = "region", seas = "season", f = "fleet")
    )
  }

  collect_message("Fdev_rho is specified as: ", Fdev_rho_spec)

  return(input_list)
}

#' Map sigma_C (catch observation error SD) parameters
#'
#' Constructs the \code{ln_sigmaC} factor map used by the TMB/RTMB objective
#' function to share or fix the log-scale standard deviation of catch observation
#' error across regions, years, seasons, and fleets. All cells within a shared
#' group are assigned the same estimation index.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists, as constructed by upstream setup functions.
#' @param sigmaC_spec Character string controlling the sharing and estimation
#'   structure for \code{ln_sigmaC}. One of:
#'   \describe{
#'     \item{\code{"est_all"}}{Unique parameter per region × year × season × fleet.}
#'     \item{\code{"est_shared_r"}}{Shared across regions.}
#'     \item{\code{"est_shared_y"}}{Shared across years.}
#'     \item{\code{"est_shared_seas"}}{Shared across seasons.}
#'     \item{\code{"est_shared_f"}}{Shared across fleets.}
#'     \item{\code{"est_shared_r_y"}}{Shared across regions and years.}
#'     \item{\code{"est_shared_r_seas"}}{Shared across regions and seasons.}
#'     \item{\code{"est_shared_r_f"}}{Shared across regions and fleets.}
#'     \item{\code{"est_shared_y_seas"}}{Shared across years and seasons.}
#'     \item{\code{"est_shared_y_f"}}{Shared across years and fleets.}
#'     \item{\code{"est_shared_seas_f"}}{Shared across seasons and fleets.}
#'     \item{\code{"est_shared_r_y_seas"}}{Shared across regions, years, and seasons.}
#'     \item{\code{"est_shared_r_y_f"}}{Shared across regions, years, and fleets.}
#'     \item{\code{"est_shared_r_seas_f"}}{Shared across regions, seasons, and fleets.}
#'     \item{\code{"est_shared_y_seas_f"}}{Shared across years, seasons, and fleets.}
#'     \item{\code{"est_shared_r_y_seas_f"}}{Single parameter shared across all dimensions.}
#'     \item{\code{"fix"}}{All \code{ln_sigmaC} parameters fixed at starting values (mapped to \code{NA}).}
#'   }
#'
#' @return The input \code{input_list} with \code{$map$ln_sigmaC} set to a factor
#'   vector of length \code{prod(dim(par$ln_sigmaC))}. Each element is an integer
#'   estimation index, or \code{NA} when \code{sigmaC_spec = "fix"}.
#'
#' @keywords internal
do_sigmaC_mapping <- function(input_list, sigmaC_spec) {

  # Sigma C -----------------------------------------------------------------
  dims <- c(region = input_list$data$n_regions,
            year   = length(input_list$data$years),
            season = input_list$data$n_seas,
            fleet  = input_list$data$n_fish_fleets)

  input_list$map$ln_sigmaC <- build_shared_spec_map(
    dims = dims, spec = sigmaC_spec,
    dim_abbrev = c(r = "region", y = "year", seas = "season", f = "fleet")
  )

  # Print Message
  collect_message("sigmaC is specified as: ", sigmaC_spec)

  return(input_list)
}

#' Map fishing mortality parameters
#'
#' Constructs the \code{ln_F_devs} and \code{ln_F_mean} factor maps, assigning
#' unique estimation indices to cells where catch data are used
#' (\code{UseCatch == 1}) and mapping cells without catch data to \code{NA}.
#' This ensures that fishing mortality parameters are only estimated for
#' dimensions with observed catch. \code{ln_F_devs} is resolved per
#' region–year–season–fleet cell, while \code{ln_F_mean} is resolved per
#' region–season–fleet cell and is estimated whenever that cell is fished in
#' at least one year.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists. Requires \code{$data$UseCatch} and \code{$data$UseCatch_pop} to be populated by
#'   \code{\link{Setup_Mod_Catch_and_F}}.
#'
#' @return The input \code{input_list} with \code{$map$ln_F_devs} and
#'   \code{$map$ln_F_mean} set to factor vectors. Cells with catch are assigned
#'   sequential integer indices; cells without catch are \code{NA}.
#'
#' @keywords internal
do_Fmort_mapping <- function(input_list) {

  dims <- c(region = input_list$data$n_regions,
            year   = length(input_list$data$years),
            season = input_list$data$n_seas,
            fleet  = input_list$data$n_fish_fleets)

  # Estimate F devs if aggregated catch OR any pop-specific catch is used, or
  # if the aggregate catch observation is missing (NA) rather than a true
  # recorded zero -- fishing is assumed to have continued through a missing
  # observation, whereas a recorded zero (or no catch data used at all) with
  # no missing observation indicates a true closure
  has_catch <- input_list$data$UseCatch == 1 |
    apply(input_list$data$UseCatch_pop == 1, c(2,3,4,5), any) |
    is.na(input_list$data$ObsCatch)

  F_dev_map <- build_pe_map(dims, share_over = character(0))
  F_dev_map[!has_catch] <- NA

  input_list$map$ln_F_devs <- factor(as.vector(F_dev_map))

  # ln_F_mean only enters the objective through cells that are fished in at
  # least one year -- the same condition that keeps Fmort free rather than
  # pinned to zero. A region-season-fleet cell that is closed in every year
  # (e.g. a fleet that only operates in its own region/season) contributes
  # nothing to the likelihood, so estimating it leaves a flat gradient
  F_mean_dims <- dims[c("region", "season", "fleet")]
  F_mean_active <- apply(has_catch, c(1,3,4), any)

  F_mean_map <- build_pe_map(F_mean_dims, share_over = character(0))
  F_mean_map[!F_mean_active] <- NA

  input_list$map$ln_F_mean <- factor(as.vector(F_mean_map))
  return(input_list)
}

#' Map population-specific catch observation error SD parameters
#'
#' Constructs the \code{ln_sigmaC_pop} factor map used by the TMB/RTMB
#' objective function to share or fix the log-scale standard deviation of
#' population-specific catch observation error across populations, regions,
#' years, seasons, and fleets. All cells within a shared group are assigned
#' the same estimation index.
#'
#' The sharing specification encodes which dimensions are collapsed into a
#' single parameter via underscore-separated tokens. For example,
#' \code{"est_shared_pop_r"} shares across populations and regions (one
#' parameter per year × season × fleet combination), while
#' \code{"est_shared_r_f"} shares across regions and fleets (one parameter
#' per population × year × season combination).
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists, as constructed by upstream setup functions. Requires
#'   \code{$data$n_pop}, \code{$data$n_regions}, \code{$data$years},
#'   \code{$data$n_seas}, \code{$data$n_fish_fleets}, and
#'   \code{$par$ln_sigmaC_pop} to be populated before calling.
#' @param sigmaC_pop_spec Character string controlling the sharing and
#'   estimation structure for \code{ln_sigmaC_pop}. One of:
#'   \describe{
#'     \item{\code{"fix"}}{All parameters fixed at starting values
#'       (mapped to \code{NA}).}
#'     \item{\code{"est_all"}}{Unique parameter per population × region ×
#'       year × season × fleet cell.}
#'     \item{\code{"est_shared_pop"}}{Shared across populations; unique per
#'       region × year × season × fleet.}
#'     \item{\code{"est_shared_r"}}{Shared across regions; unique per
#'       population × year × season × fleet.}
#'     \item{\code{"est_shared_y"}}{Shared across years; unique per
#'       population × region × season × fleet.}
#'     \item{\code{"est_shared_seas"}}{Shared across seasons; unique per
#'       population × region × year × fleet.}
#'     \item{\code{"est_shared_f"}}{Shared across fleets; unique per
#'       population × region × year × season.}
#'     \item{\code{"est_shared_pop_r"}}{Shared across populations and regions.}
#'     \item{\code{"est_shared_pop_y"}}{Shared across populations and years.}
#'     \item{\code{"est_shared_pop_seas"}}{Shared across populations and seasons.}
#'     \item{\code{"est_shared_pop_f"}}{Shared across populations and fleets.}
#'     \item{\code{"est_shared_r_y"}}{Shared across regions and years.}
#'     \item{\code{"est_shared_r_seas"}}{Shared across regions and seasons.}
#'     \item{\code{"est_shared_r_f"}}{Shared across regions and fleets.}
#'     \item{\code{"est_shared_y_seas"}}{Shared across years and seasons.}
#'     \item{\code{"est_shared_y_f"}}{Shared across years and fleets.}
#'     \item{\code{"est_shared_seas_f"}}{Shared across seasons and fleets.}
#'     \item{\code{"est_shared_pop_r_y"}}{Shared across populations, regions,
#'       and years.}
#'     \item{\code{"est_shared_pop_r_seas"}}{Shared across populations, regions,
#'       and seasons.}
#'     \item{\code{"est_shared_pop_r_f"}}{Shared across populations, regions,
#'       and fleets.}
#'     \item{\code{"est_shared_pop_y_seas"}}{Shared across populations, years,
#'       and seasons.}
#'     \item{\code{"est_shared_pop_y_f"}}{Shared across populations, years,
#'       and fleets.}
#'     \item{\code{"est_shared_pop_seas_f"}}{Shared across populations, seasons,
#'       and fleets.}
#'     \item{\code{"est_shared_r_y_seas"}}{Shared across regions, years,
#'       and seasons.}
#'     \item{\code{"est_shared_r_y_f"}}{Shared across regions, years,
#'       and fleets.}
#'     \item{\code{"est_shared_r_seas_f"}}{Shared across regions, seasons,
#'       and fleets.}
#'     \item{\code{"est_shared_y_seas_f"}}{Shared across years, seasons,
#'       and fleets.}
#'     \item{\code{"est_shared_pop_r_y_seas"}}{Shared across populations,
#'       regions, years, and seasons.}
#'     \item{\code{"est_shared_pop_r_y_f"}}{Shared across populations, regions,
#'       years, and fleets.}
#'     \item{\code{"est_shared_pop_r_seas_f"}}{Shared across populations,
#'       regions, seasons, and fleets.}
#'     \item{\code{"est_shared_pop_y_seas_f"}}{Shared across populations, years,
#'       seasons, and fleets.}
#'     \item{\code{"est_shared_r_y_seas_f"}}{Shared across regions, years,
#'       seasons, and fleets.}
#'     \item{\code{"est_shared_pop_r_y_seas_f"}}{Single parameter shared across
#'       all dimensions.}
#'   }
#'
#'
#' @keywords internal
do_sigmaC_pop_mapping <- function(input_list, sigmaC_pop_spec) {

  dims <- c(pop    = input_list$data$n_pop,
            region = input_list$data$n_regions,
            year   = length(input_list$data$years),
            season = input_list$data$n_seas,
            fleet  = input_list$data$n_fish_fleets)

  input_list$map$ln_sigmaC_pop <- build_shared_spec_map(
    dims = dims, spec = sigmaC_pop_spec,
    dim_abbrev = c(pop = "pop", r = "region", y = "year", seas = "season", f = "fleet")
  )

  collect_message("sigmaC_pop is specified as: ", sigmaC_pop_spec)
  return(input_list)
}

#' Map sigma_dmr (discard mortality process error SD) parameters
#'
#' Constructs the \code{ln_sigma_dmr} factor map used by the TMB/RTMB objective
#' function to share or fix the logit-scale standard deviation of discard mortality
#' process error across regions, seasons, and fleets. All cells within a shared
#' group are assigned the same estimation index.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists, as constructed by upstream setup functions.
#' @param sigma_dmr_spec Character string controlling the sharing and estimation
#'   structure for \code{ln_sigma_dmr}. One of:
#'   \describe{
#'     \item{\code{"est_all"}}{Unique parameter per region × season × fleet combination.}
#'     \item{\code{"est_shared_r"}}{Shared across regions; unique per season × fleet.}
#'     \item{\code{"est_shared_seas"}}{Shared across seasons; unique per region × fleet.}
#'     \item{\code{"est_shared_f"}}{Shared across fleets; unique per region × season.}
#'     \item{\code{"est_shared_r_seas"}}{Shared across regions and seasons; unique per fleet.}
#'     \item{\code{"est_shared_r_f"}}{Shared across regions and fleets; unique per season.}
#'     \item{\code{"est_shared_seas_f"}}{Shared across seasons and fleets; unique per region.}
#'     \item{\code{"est_shared_r_seas_f"}}{Single parameter shared across all dimensions.}
#'     \item{\code{"fix"}}{All \code{ln_sigma_dmr} parameters fixed at starting values (mapped to \code{NA}).}
#'   }
#'
#' @return The input \code{input_list} with \code{$map$ln_sigma_dmr} set to a factor
#'   vector of length \code{prod(dim(par$ln_sigma_dmr))}. Each element is an integer
#'   estimation index for shared or estimated configurations, or \code{NA} when
#'   \code{sigma_dmr_spec = "fix"}.
#'
#' @keywords internal
do_sigma_dmr_mapping <- function(input_list, sigma_dmr_spec) {

  dims <- c(region = input_list$data$n_regions,
            season = input_list$data$n_seas,
            fleet  = input_list$data$n_fish_fleets)

  input_list$map$ln_sigma_dmr <- build_shared_spec_map(
    dims = dims, spec = sigma_dmr_spec,
    dim_abbrev = c(r = "region", seas = "season", f = "fleet")
  )

  # Print Message
  collect_message("sigma_dmr is specified as: ", sigma_dmr_spec)

  return(input_list)
}

#' Map sigmaD (discard mortality rate observation error SD) parameters
#'
#' Constructs the \code{ln_sigmaD} factor map used by the TMB/RTMB objective
#' function to share or fix the log-scale standard deviation of discard mrotality observation
#' error across regions, years, seasons, and fleets. All cells within a shared
#' group are assigned the same estimation index.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists, as constructed by upstream setup functions.
#' @param sigmaD_spec Character string controlling the sharing and estimation
#'   structure for \code{ln_sigmaD}. One of:
#'   \describe{
#'     \item{\code{"est_all"}}{Unique parameter per region × year × season × fleet.}
#'     \item{\code{"est_shared_r"}}{Shared across regions.}
#'     \item{\code{"est_shared_y"}}{Shared across years.}
#'     \item{\code{"est_shared_seas"}}{Shared across seasons.}
#'     \item{\code{"est_shared_f"}}{Shared across fleets.}
#'     \item{\code{"est_shared_r_y"}}{Shared across regions and years.}
#'     \item{\code{"est_shared_r_seas"}}{Shared across regions and seasons.}
#'     \item{\code{"est_shared_r_f"}}{Shared across regions and fleets.}
#'     \item{\code{"est_shared_y_seas"}}{Shared across years and seasons.}
#'     \item{\code{"est_shared_y_f"}}{Shared across years and fleets.}
#'     \item{\code{"est_shared_seas_f"}}{Shared across seasons and fleets.}
#'     \item{\code{"est_shared_r_y_seas"}}{Shared across regions, years, and seasons.}
#'     \item{\code{"est_shared_r_y_f"}}{Shared across regions, years, and fleets.}
#'     \item{\code{"est_shared_r_seas_f"}}{Shared across regions, seasons, and fleets.}
#'     \item{\code{"est_shared_y_seas_f"}}{Shared across years, seasons, and fleets.}
#'     \item{\code{"est_shared_r_y_seas_f"}}{Single parameter shared across all dimensions.}
#'     \item{\code{"fix"}}{All \code{ln_sigmaD} parameters fixed at starting values (mapped to \code{NA}).}
#'   }
#'
#' @return The input \code{input_list} with \code{$map$ln_sigmaD} set to a factor
#'   vector of length \code{prod(dim(par$ln_sigmaD))}. Each element is an integer
#'   estimation index, or \code{NA} when \code{sigmaD_spec = "fix"}.
#'
#' @keywords internal
do_sigmaD_mapping <- function(input_list, sigmaD_spec) {

  dims <- c(region = input_list$data$n_regions,
            year   = length(input_list$data$years),
            season = input_list$data$n_seas,
            fleet  = input_list$data$n_fish_fleets)

  input_list$map$ln_sigmaD <- build_shared_spec_map(
    dims = dims, spec = sigmaD_spec,
    dim_abbrev = c(r = "region", y = "year", seas = "season", f = "fleet")
  )

  # Print Message
  collect_message("sigmaD is specified as: ", sigmaD_spec)

  return(input_list)
}

#' Map discard mortality deviation parameters
#'
#' Constructs the \code{logit_dmr_devs} factor map, assigning unique estimation
#' indices to region–year–season–fleet cells where discard data are used
#' (\code{UseDiscard == 1}) and mapping cells without discard data to \code{NA}.
#' This ensures that \code{logit_dmr_devs} parameters are only estimated for
#' dimensions with observed discard.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists. Requires \code{$data$UseDiscard} and \code{$data$UseDiscard_pop} to be populated by
#'   \code{\link{Setup_Mod_Catch_and_F}}.
#' @param dmr_dev_spec Character string specifying whether to estimate or fix
#'   deviations. Currently supports \code{"est_all"} and \code{"fix"}.
#'
#' @return The input \code{input_list} with \code{$map$logit_dmr_devs} set to a
#'   factor vector. Cells with discard are assigned sequential integer indices;
#'   cells without discard are \code{NA}.
#'
#' @keywords internal
do_dmr_dev_mapping <- function(input_list, dmr_dev_spec) {

  valid_specs <- c("fix", "est_all")
  if(!dmr_dev_spec %in% valid_specs)
    stop("dmr_dev_spec '", dmr_dev_spec, "' not recognized. Valid options: ",
         paste(valid_specs, collapse = ", "))

  if(dmr_dev_spec == "fix") {
    input_list$map$logit_dmr_devs <- factor(rep(NA, length(input_list$par$logit_dmr_devs)))
    collect_message("dmr_devs is specified as: fix")
    return(input_list)
  }

  # est: estimate devs where discard data exist
  dmr_dev_map <- input_list$par$logit_dmr_devs
  dmr_dev_map[] <- NA
  dmr_dev_counter <- 1

  for(r in 1:input_list$data$n_regions) {
    for(y in 1:length(input_list$data$years)) {
      for(seas in 1:input_list$data$n_seas) {
        for(f in 1:input_list$data$n_fish_fleets) {

          has_agg_discard <- input_list$data$UseCatch[r,y,seas,f] == 1
          has_pop_discard <- any(input_list$data$UseCatch_pop[,r,y,seas,f] == 1)

          if(has_agg_discard || has_pop_discard) {
            dmr_dev_map[r,y,seas,f] <- dmr_dev_counter
            dmr_dev_counter <- dmr_dev_counter + 1
          }

        } # end f
      } # end seas
    } # end y
  } # end r

  input_list$map$logit_dmr_devs <- factor(dmr_dev_map)
  collect_message("dmr_devs is specified as: est")
  return(input_list)
}

#' Map population-specific discard observation error SD parameters
#'
#' Constructs the \code{ln_sigmaD_pop} factor map used by the TMB/RTMB
#' objective function to share or fix the log-scale standard deviation of
#' population-specific discard observation error across populations, regions,
#' years, seasons, and fleets. All cells within a shared group are assigned
#' the same estimation index.
#'
#' The sharing specification encodes which dimensions are collapsed into a
#' single parameter via underscore-separated tokens. For example,
#' \code{"est_shared_pop_r"} shares across populations and regions (one
#' parameter per year × season × fleet combination), while
#' \code{"est_shared_r_f"} shares across regions and fleets (one parameter
#' per population × year × season combination).
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists, as constructed by upstream setup functions. Requires
#'   \code{$data$n_pop}, \code{$data$n_regions}, \code{$data$years},
#'   \code{$data$n_seas}, \code{$data$n_fish_fleets}, and
#'   \code{$par$ln_sigmaD_pop} to be populated before calling.
#' @param sigmaD_pop_spec Character string controlling the sharing and
#'   estimation structure for \code{ln_sigmaD_pop}. One of:
#'   \describe{
#'     \item{\code{"fix"}}{All parameters fixed at starting values
#'       (mapped to \code{NA}).}
#'     \item{\code{"est_all"}}{Unique parameter per population × region ×
#'       year × season × fleet cell.}
#'     \item{\code{"est_shared_pop"}}{Shared across populations; unique per
#'       region × year × season × fleet.}
#'     \item{\code{"est_shared_r"}}{Shared across regions; unique per
#'       population × year × season × fleet.}
#'     \item{\code{"est_shared_y"}}{Shared across years; unique per
#'       population × region × season × fleet.}
#'     \item{\code{"est_shared_seas"}}{Shared across seasons; unique per
#'       population × region × year × fleet.}
#'     \item{\code{"est_shared_f"}}{Shared across fleets; unique per
#'       population × region × year × season.}
#'     \item{\code{"est_shared_pop_r"}}{Shared across populations and regions.}
#'     \item{\code{"est_shared_pop_y"}}{Shared across populations and years.}
#'     \item{\code{"est_shared_pop_seas"}}{Shared across populations and seasons.}
#'     \item{\code{"est_shared_pop_f"}}{Shared across populations and fleets.}
#'     \item{\code{"est_shared_r_y"}}{Shared across regions and years.}
#'     \item{\code{"est_shared_r_seas"}}{Shared across regions and seasons.}
#'     \item{\code{"est_shared_r_f"}}{Shared across regions and fleets.}
#'     \item{\code{"est_shared_y_seas"}}{Shared across years and seasons.}
#'     \item{\code{"est_shared_y_f"}}{Shared across years and fleets.}
#'     \item{\code{"est_shared_seas_f"}}{Shared across seasons and fleets.}
#'     \item{\code{"est_shared_pop_r_y"}}{Shared across populations, regions,
#'       and years.}
#'     \item{\code{"est_shared_pop_r_seas"}}{Shared across populations, regions,
#'       and seasons.}
#'     \item{\code{"est_shared_pop_r_f"}}{Shared across populations, regions,
#'       and fleets.}
#'     \item{\code{"est_shared_pop_y_seas"}}{Shared across populations, years,
#'       and seasons.}
#'     \item{\code{"est_shared_pop_y_f"}}{Shared across populations, years,
#'       and fleets.}
#'     \item{\code{"est_shared_pop_seas_f"}}{Shared across populations, seasons,
#'       and fleets.}
#'     \item{\code{"est_shared_r_y_seas"}}{Shared across regions, years,
#'       and seasons.}
#'     \item{\code{"est_shared_r_y_f"}}{Shared across regions, years,
#'       and fleets.}
#'     \item{\code{"est_shared_r_seas_f"}}{Shared across regions, seasons,
#'       and fleets.}
#'     \item{\code{"est_shared_y_seas_f"}}{Shared across years, seasons,
#'       and fleets.}
#'     \item{\code{"est_shared_pop_r_y_seas"}}{Shared across populations,
#'       regions, years, and seasons.}
#'     \item{\code{"est_shared_pop_r_y_f"}}{Shared across populations, regions,
#'       years, and fleets.}
#'     \item{\code{"est_shared_pop_r_seas_f"}}{Shared across populations,
#'       regions, seasons, and fleets.}
#'     \item{\code{"est_shared_pop_y_seas_f"}}{Shared across populations, years,
#'       seasons, and fleets.}
#'     \item{\code{"est_shared_r_y_seas_f"}}{Shared across regions, years,
#'       seasons, and fleets.}
#'     \item{\code{"est_shared_pop_r_y_seas_f"}}{Single parameter shared across
#'       all dimensions.}
#'   }
#'
#'
#' @keywords internal
do_sigmaD_pop_mapping <- function(input_list, sigmaD_pop_spec) {

  dims <- c(pop    = input_list$data$n_pop,
            region = input_list$data$n_regions,
            year   = length(input_list$data$years),
            season = input_list$data$n_seas,
            fleet  = input_list$data$n_fish_fleets)

  input_list$map$ln_sigmaD_pop <- build_shared_spec_map(
    dims = dims, spec = sigmaD_pop_spec,
    dim_abbrev = c(pop = "pop", r = "region", y = "year", seas = "season", f = "fleet")
  )

  collect_message("sigmaD_pop is specified as: ", sigmaD_pop_spec)
  return(input_list)
}

#' Map discard mortality process error SD parameters
#'
#' Constructs the \code{ln_sigma_dmr} factor map used in the TMB/RTMB
#' objective function to share or fix log-scale standard deviations of
#' discard mortality process error across regions, seasons, and fleets.
#'
#' The mapping assigns integer indices to parameter groups defined by
#' \code{dmr_mean_spec}. Cells within the same group share a single
#' estimated parameter.
#'
#' @param input_list Named list containing \code{$data}, \code{$par}, and
#'   \code{$map}.
#'
#' @param dmr_mean_spec Character string specifying sharing structure.
#'   Options include \code{"est_all"}, \code{"est_shared_r"},
#'   \code{"est_shared_seas"}, \code{"est_shared_f"}, and combinations thereof,
#'   or \code{"fix"}.
#'
#' @return Updated \code{input_list} with \code{$map$ln_sigma_dmr}.
#'
#' @keywords internal
do_dmr_mean_mapping <- function(input_list, dmr_mean_spec) {

  dims <- c(region = input_list$data$n_regions,
            season = input_list$data$n_seas,
            fleet  = input_list$data$n_fish_fleets)

  input_list$map$logit_dmr_mean <- build_shared_spec_map(
    dims = dims, spec = dmr_mean_spec,
    dim_abbrev = c(r = "region", seas = "season", f = "fleet")
  )

  collect_message("dmr_mean is specified as: ", dmr_mean_spec)
  return(input_list)
}

#' Set up fishing mortality, discard mortality, and catch observation inputs
#'
#' Populates \code{input_list} with observed catch, catch usage indicators,
#' fishing mortality parameters (\code{ln_F_mean}, \code{ln_F_devs}), and
#' observation/process error structures (\code{ln_sigmaC}, \code{ln_sigmaC_pop},
#' \code{ln_sigmaF}). Also populates discard observations, discard mortality
#' rate parameters (\code{logit_dmr_mean}, \code{logit_dmr_devs}), and
#' discard observation/process error structures (\code{ln_sigmaD},
#' \code{ln_sigmaD_pop}, \code{ln_sigma_dmr}).
#' Must be called after \code{\link{Setup_Mod_Biologicals}}.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, \code{$map},
#'   and \code{$verbose} sublists, as returned by upstream setup functions.
#' @param ObsCatch Observed aggregated catch array
#'   \code{[n_regions x n_years x n_seas x n_fish_fleets]}.
#'   Values should be in the units specified by \code{catch_units}. For a
#'   cell with \code{UseCatch == 0} (and no population-specific catch
#'   used), an \code{NA} entry here is treated as a genuinely missing
#'   observation -- fishing is assumed to have continued and \code{Fmort}/
#'   \code{ln_F_devs} are estimated normally for that year -- whereas a
#'   true recorded value (typically \code{0}) is treated as a real
#'   closure: \code{Fmort} is forced to zero and no deviation is estimated.
#'   See \code{\link{Get_Fdev_PE_loglik}}.
#' @param ObsCatch_pop Observed population-specific catch array
#'   \code{[n_pop x n_regions x n_years x n_seas x n_fish_fleets]}.
#'   Values should be in the units specified by \code{catch_units}.
#' @param UseCatch Binary indicator array
#'   \code{[n_regions x n_years x n_seas x n_fish_fleets]} controlling which
#'   aggregated catch observations enter the likelihood and whether
#'   \code{ln_F_devs} are estimated for each cell. \code{1} = use;
#'   \code{0} = exclude, unless \code{ObsCatch} is \code{NA} at that cell
#'   (see \code{ObsCatch} above), in which case \code{ln_F_devs} is still
#'   estimated as an ordinary active year despite not being fit against an
#'   observation.
#' @param UseCatch_pop Binary indicator array
#'   \code{[n_pop x n_regions x n_years x n_seas x n_fish_fleets]} controlling
#'   which population-specific catch observations enter the likelihood.
#'   \code{1} = use; \code{0} = exclude.
#' @param catch_units Character array \code{[n_fish_fleets]} specifying catch
#'   units per fleet. \code{"biom"} = biomass (default); \code{"abd"} =
#'   abundance. Converted internally to \code{0}/\code{1} integer codes.
#' @param Use_F_pen Integer flag for applying a fishing mortality penalty to
#'   penalise large deviations in \code{ln_F_devs}. \code{1} = apply
#'   (default); \code{0} = do not apply.
#' @param sigmaC_spec Character string specifying the sharing structure for
#'   \code{ln_sigmaC} (aggregated catch observation error SD). Default
#'   \code{"fix"} holds \code{ln_sigmaC} at its starting value
#'   (\code{log(0.01)} unless overridden via \code{...}). Sharing options
#'   follow the convention \code{"est_shared_<dims>"} where \code{<dims>} is
#'   an underscore-separated list of dimensions to collapse: \code{"r"}
#'   (regions), \code{"y"} (years), \code{"seas"} (seasons), \code{"f"}
#'   (fleets), or any combination (e.g., \code{"est_shared_r_y"},
#'   \code{"est_shared_r_y_seas_f"}). Use \code{"est_all"} for a fully
#'   independent parameter per cell. A warning is issued if \code{"fix"} is
#'   selected without providing a starting value in \code{...}.
#' @param sigmaC_pop_spec Character string specifying the sharing structure for
#'   \code{ln_sigmaC_pop} (population-specific catch observation error SD).
#'   Default \code{"fix"} holds \code{ln_sigmaC_pop} at its starting value
#'   (\code{log(0.01)} unless overridden via \code{...}). Sharing options
#'   follow the same convention as \code{sigmaC_spec} but with an additional
#'   population dimension: e.g., \code{"est_shared_pop"} shares across
#'   populations, \code{"est_shared_pop_r"} shares across populations and
#'   regions, and \code{"est_shared_pop_r_y_seas_f"} collapses all dimensions
#'   into a single parameter. A warning is issued if \code{"fix"} is selected
#'   without providing a starting value in \code{...}.
#' @param sigmaF_spec Character string specifying the sharing structure for
#'   \code{ln_sigmaF} (fishing mortality process error SD). Default
#'   \code{"fix"} holds \code{ln_sigmaF} at its starting value (\code{log(1)},
#'   i.e., \eqn{\sigma_F = 1}, unless overridden via \code{...}). A warning is
#'   issued if \code{"fix"} is selected without providing a starting value
#'   in \code{...}.
#' @param Fdev_model Character string specifying the process error structure
#'   for \code{ln_F_devs}. One of \code{"iid"} (default; independent
#'   deviations), \code{"rw"} (random walk; the first catch-active year per
#'   region/season/fleet is initialized with a diffuse \eqn{N(0,5)} prior),
#'   or \code{"ar1"} (first-order autoregressive; the first catch-active year
#'   is drawn from its stationary marginal distribution, and \code{Fdev_rho_spec}
#'   controls the AR1 correlation parameter). Catch-active years do not need
#'   to be contiguous for \code{"rw"} or \code{"ar1"}: the transition between
#'   two active years spanning a gap of \eqn{d} closed years is taken over
#'   the elapsed gap directly (the same marginal transition as estimating
#'   deviations for the closed years and integrating them out, without
#'   actually estimating them) -- see \code{\link{Get_Fdev_PE_loglik}}.
#'   A warning is issued if \code{"rw"} or \code{"ar1"} is selected but
#'   \code{Use_F_pen = 0} (the penalty is never evaluated, so the process
#'   structure has no effect), \code{sigmaF_spec = "fix"} (the process error
#'   SD is not estimated), or (for \code{"ar1"}) \code{Fdev_rho_spec =
#'   "fix"} (the correlation is not estimated) -- any of these may be
#'   intentional, but are common oversights when switching away from
#'   \code{"iid"}.
#' @param Fdev_rho_spec Character string specifying the sharing structure for
#'   the AR1 correlation parameter \code{Fdev_rho}, following the same
#'   convention as \code{sigmaF_spec}. Only used when \code{Fdev_model =
#'   "ar1"}; ignored (and mapped entirely to \code{NA}) otherwise.
#' @param ObsDiscard Observed aggregated discard array
#'   \code{[n_regions x n_years x n_seas x n_fish_fleets]}.
#'   Values should be in the units specified by \code{discard_units}.
#'   Default: \code{NULL} (no discard observations).
#' @param UseDiscard Binary indicator array
#'   \code{[n_regions x n_years x n_seas x n_fish_fleets]} controlling which
#'   aggregated discard observations enter the likelihood. \code{1} = use;
#'   \code{0} = exclude. Default: all zeros.
#' @param discard_units Character array \code{[n_fish_fleets]} specifying
#'   discard units per fleet. \code{"abd"} = abundance (\code{0}),
#'   \code{"biom"} = biomass (\code{1}), \code{"abd_frac"} = abundance
#'   fraction (\code{2}), \code{"biom_frac"} = biomass fraction (\code{3},
#'   default). Converted internally to integer codes.
#' @param UseDiscard_pop Binary indicator array
#'   \code{[n_pop x n_regions x n_years x n_seas x n_fish_fleets]} controlling
#'   which population-specific discard observations enter the likelihood.
#'   \code{1} = use; \code{0} = exclude. Default: all zeros.
#' @param ObsDiscard_pop Observed population-specific discard array
#'   \code{[n_pop x n_regions x n_years x n_seas x n_fish_fleets]}.
#'   Values should be in the units specified by \code{discard_units}.
#'   Default: \code{NULL} (no population-specific discard observations).
#' @param Use_dmr_pen Integer flag for applying a discard mortality rate
#'   penalty to penalise large deviations in \code{logit_dmr_devs}.
#'   \code{1} = apply; \code{0} = do not apply (default). Must be \code{1}
#'   when \code{dmr_dev_spec = "est_all"} and \code{0} when
#'   \code{dmr_dev_spec = "fix"}.
#' @param sigmaD_spec Character string specifying the sharing structure for
#'   \code{ln_sigmaD} (aggregated discard observation error SD). Default
#'   \code{"fix"} holds \code{ln_sigmaD} at its starting value
#'   (\code{log(0.01)} unless overridden via \code{...}). Sharing options
#'   follow the same convention as \code{sigmaC_spec}. A warning is issued
#'   if \code{"fix"} is selected without providing a starting value in
#'   \code{...}.
#' @param sigmaD_pop_spec Character string specifying the sharing structure for
#'   \code{ln_sigmaD_pop} (population-specific discard observation error SD).
#'   Default \code{"fix"} holds \code{ln_sigmaD_pop} at its starting value
#'   (\code{log(0.01)} unless overridden via \code{...}). Sharing options
#'   follow the same convention as \code{sigmaC_pop_spec}. A warning is
#'   issued if \code{"fix"} is selected without providing a starting value
#'   in \code{...}.
#' @param sigma_dmr_spec Character string specifying the sharing structure for
#'   \code{ln_sigma_dmr} (discard mortality rate process error SD). Default
#'   \code{"fix"} holds \code{ln_sigma_dmr} at its starting value
#'   (\code{log(1)} unless overridden via \code{...}). Sharing options
#'   follow the same convention as \code{sigmaF_spec}. A warning is issued
#'   if \code{"fix"} is selected without providing a starting value in
#'   \code{...}.
#' @param dmr_mean_spec Character string specifying the sharing/estimation
#'   structure for \code{logit_dmr_mean} (logit-scale mean discard mortality
#'   rate). Default \code{"fix"} holds at its starting value (\code{0},
#'   i.e., DMR = 0.5 on the natural scale, unless overridden via \code{...}).
#'   See \code{\link{do_dmr_mean_mapping}} for sharing options.
#' @param dmr_dev_spec Character string specifying the sharing/estimation
#'   structure for \code{logit_dmr_devs} (logit-scale annual discard mortality
#'   rate deviations). Default \code{"fix"} holds deviations at zero
#'   (unless overridden via \code{...}). Use \code{"est"} to estimate
#'   deviations; requires \code{Use_dmr_pen = 1}. See
#'   \code{\link{do_dmr_dev_mapping}} for sharing options.
#' @param ... Optional starting value overrides for catch and discard related parameters.
#'
#' @return The input \code{input_list} with \code{$data}, \code{$par}, and
#'   \code{$map} updated. Key additions:
#'   \describe{
#'     \item{\code{$data}}{
#'       \code{ObsCatch}, \code{ObsCatch_pop}, \code{UseCatch},
#'       \code{UseCatch_pop}, \code{Use_F_pen}, \code{catch_units},
#'       \code{Fdev_model},
#'       \code{ObsDiscard}, \code{ObsDiscard_pop}, \code{UseDiscard},
#'       \code{UseDiscard_pop}, \code{Use_dmr_pen}, \code{discard_units}.}
#'     \item{\code{$par}}{
#'       \code{ln_sigmaC}, \code{ln_sigmaC_pop}, \code{ln_sigmaF},
#'       \code{Fdev_rho},
#'       \code{ln_F_mean}, \code{ln_F_devs},
#'       \code{ln_sigmaD}, \code{ln_sigmaD_pop}, \code{ln_sigma_dmr},
#'       \code{logit_dmr_mean}, \code{logit_dmr_devs}.}
#'     \item{\code{$map}}{
#'       \code{ln_sigmaC}, \code{ln_sigmaC_pop}, \code{ln_sigmaF},
#'       \code{Fdev_rho},
#'       \code{ln_F_mean}, \code{ln_F_devs},
#'       \code{ln_sigmaD}, \code{ln_sigmaD_pop}, \code{ln_sigma_dmr},
#'       \code{logit_dmr_mean}, \code{logit_dmr_devs}.}
#'   }
#'
#' @export Setup_Mod_Catch_and_F
#' @family Model Setup
Setup_Mod_Catch_and_F <- function(input_list,

                                  # Retained Catch Stuff
                                  ObsCatch,
                                  UseCatch,
                                  catch_units = array("biom", dim = c(input_list$data$n_fish_fleets)),
                                  UseCatch_pop = array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions,
                                                                  length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets )),
                                  ObsCatch_pop = NULL,
                                  Use_F_pen = 1,
                                  sigmaC_spec = "fix",
                                  sigmaC_pop_spec = 'fix',
                                  sigmaF_spec = "fix",
                                  Fdev_model = "iid",
                                  Fdev_rho_spec = "fix",

                                  # Discarded Catch Stuff
                                  ObsDiscard = NULL,
                                  UseDiscard = array(0, dim = c(input_list$data$n_regions,
                                                                length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets )),
                                  discard_units = array("biom_frac", dim = c(input_list$data$n_fish_fleets)),
                                  UseDiscard_pop = array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions,
                                                                    length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets )),
                                  ObsDiscard_pop = NULL,
                                  Use_dmr_pen = 0,
                                  sigmaD_spec = "fix",
                                  sigmaD_pop_spec = 'fix',
                                  sigma_dmr_spec = "fix",
                                  dmr_mean_spec = 'fix',
                                  dmr_dev_spec = 'fix',
                                  ...
) {

  messages_list <<- character(0) # string to attach to for printing messages
  starting_values <- list(...)
  if(input_list$store_config) input_list$config$Setup_Mod_Catch_and_F <- mget(names(formals()))[-1]

  # Input Validation --------------------------------------------------------

  # Catch objects
  check_data_dimensions(ObsCatch, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsCatch')
  check_data_dimensions(UseCatch, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'UseCatch')

  if(any(UseCatch_pop == 1)) {
    check_data_dimensions(ObsCatch_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsCatch_pop')
    check_data_dimensions(UseCatch_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'UseCatch_pop')
  }

  if(any(UseDiscard == 1)) {
    check_data_dimensions(ObsDiscard, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsDiscard')
    check_data_dimensions(UseDiscard, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'UseDiscard')
  }

  if(any(UseDiscard_pop == 1)) {
    check_data_dimensions(ObsDiscard_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsDiscard_pop')
    check_data_dimensions(UseDiscard_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'UseDiscard_pop')
  }

  # Fishing Mortality checking
  if(!Use_F_pen %in% c(0,1)) stop("Use_F_pen incorrectly specified. Either set at 0 (don't use F penalty) or 1 (use F penalty)")
  else collect_message("Fishing mortality penalty is: ", ifelse(Use_F_pen == 0, 'Not Used', "Used"))
  if(sigmaC_spec == "fix" && !("ln_sigmaC" %in% names(starting_values))) warning("sigmaC is specified as fix, but no starting values / fixed values are provided. Either do this post-hoc, or use the ... argument if you do not want to use default values")
  if(sigmaC_pop_spec == "fix" && !("ln_sigmaC_pop" %in% names(starting_values))) warning("sigmaC_pop is specified as fix, but no starting values / fixed values are provided. Either do this post-hoc, or use the ... argument if you do not want to use default values")
  if(sigmaF_spec == "fix" && !("ln_sigmaF" %in% names(starting_values))) warning("sigmaF_spec is specified as fix, but no starting values / fixed values are provided. Either do this post-hoc, or use the ... argument if you do not want to use default values")

  # Fdev_model checking
  if(!Fdev_model %in% c("iid", "rw", "ar1")) stop("Fdev_model incorrectly specified. Must be one of 'iid', 'rw', or 'ar1'")
  else collect_message("Fdev_model is specified as: ", Fdev_model)

  if(Fdev_model %in% c("rw", "ar1") && Use_F_pen == 0)
    warning("Fdev_model = '", Fdev_model, "' but Use_F_pen = 0 -- the fishing mortality deviation penalty is never evaluated, so the ", Fdev_model, " process error structure has no effect on the model. Set Use_F_pen = 1 to actually apply it.")

  if(Fdev_model %in% c("rw", "ar1") && sigmaF_spec == "fix")
    warning("Fdev_model = '", Fdev_model, "' but sigmaF_spec = 'fix' -- the process error standard deviation (ln_sigmaF) driving the ", Fdev_model, " process is not being estimated. This may be intentional (e.g. fixing sigma at a known value), but if not, consider estimating ln_sigmaF via sigmaF_spec.")

  if(Fdev_model == "ar1" && Fdev_rho_spec == "fix")
    warning("Fdev_model = 'ar1' but Fdev_rho_spec = 'fix' -- the AR1 correlation parameter (Fdev_rho) is not being estimated. This may be intentional (e.g. fixing rho at a known value), but if not, consider estimating Fdev_rho via Fdev_rho_spec.")

  # Discard Mortality checking
  if(!Use_dmr_pen %in% c(0,1)) stop("Use_dmr_pen incorrectly specified. Either set at 0 (don't use D penalty) or 1 (use D penalty)")
  else collect_message("Discard mortality penalty is: ", ifelse(Use_dmr_pen == 0, 'Not Used', "Used"))
  if(sigmaD_spec == "fix" && !("ln_sigmaD" %in% names(starting_values))) warning("sigmaD is specified as fix, but no starting values / fixed values are provided. Either do this post-hoc, or use the ... argument if you do not want to use default values")
  if(sigmaD_pop_spec == "fix" && !("ln_sigmaD_pop" %in% names(starting_values))) warning("sigmaD_pop is specified as fix, but no starting values / fixed values are provided. Either do this post-hoc, or use the ... argument if you do not want to use default values")
  if(sigma_dmr_spec == "fix" && !("ln_sigma_dmr" %in% names(starting_values))) warning("sigma_dmr_spec is specified as fix, but no starting values / fixed values are provided. Either do this post-hoc, or use the ... argument if you do not want to use default values")

  # Validation checks for dmr_dev_spec and Use_dmr_pen consistency
  if(dmr_dev_spec == "est" && Use_dmr_pen == 0)
    warning("dmr_dev_spec is 'est' but Use_dmr_pen is 0. Estimating dmr deviations without a penalty will likely cause convergence issues.")

  if(dmr_dev_spec == "fix" && Use_dmr_pen == 1)
    stop("dmr_dev_spec is 'fix' but Use_dmr_pen is 1. Cannot apply a penalty on deviations that are not estimated.")

  # Catch units
  catch_units[catch_units == 'abd'] <- 0
  catch_units[catch_units == 'biom'] <- 1
  catch_units <- array(as.numeric(catch_units), dim = c(input_list$data$n_fish_fleets)) # convert to numeric array

  # Discard units
  discard_units[discard_units == 'abd'] <- 0
  discard_units[discard_units == 'biom'] <- 1
  discard_units[discard_units == 'abd_frac'] <- 2
  discard_units[discard_units == 'biom_frac'] <- 3
  discard_units <- array(as.numeric(discard_units), dim = c(input_list$data$n_fish_fleets)) # convert to numeric array

  # Populate Data List ------------------------------------------------------

  # Retained Catch Stuff
  input_list$data$ObsCatch <- ObsCatch
  input_list$data$UseCatch <- UseCatch
  input_list$data$ObsCatch_pop <- ObsCatch_pop
  input_list$data$UseCatch_pop <- UseCatch_pop
  input_list$data$Use_F_pen <- Use_F_pen
  input_list$data$catch_units <- catch_units
  input_list$data$Fdev_model <- match(Fdev_model, c("iid", "rw", "ar1")) # 1 = iid, 2 = rw, 3 = ar1

  # Discarded Catch Stuff
  input_list$data$ObsDiscard <- ObsDiscard
  input_list$data$UseDiscard <- UseDiscard
  input_list$data$ObsDiscard_pop <- ObsDiscard_pop
  input_list$data$UseDiscard_pop <- UseDiscard_pop
  input_list$data$Use_dmr_pen <- Use_dmr_pen
  input_list$data$discard_units <- discard_units

  # Populate Parameter List -------------------------------------------------

  # Catch observation error
  if("ln_sigmaC" %in% names(starting_values)) input_list$par$ln_sigmaC <- starting_values$ln_sigmaC
  else input_list$par$ln_sigmaC <- array(log(0.01), dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets))

  if("ln_sigmaC_pop" %in% names(starting_values)) input_list$par$ln_sigmaC_pop <- starting_values$ln_sigmaC_pop
  else input_list$par$ln_sigmaC_pop <- array(log(0.01), dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets))

  # Process error fishing deviations
  if("ln_sigmaF" %in% names(starting_values)) input_list$par$ln_sigmaF <- starting_values$ln_sigmaF
  else input_list$par$ln_sigmaF <- array(log(1), dim = c(input_list$data$n_regions, input_list$data$n_seas, input_list$data$n_fish_fleets))

  # AR1 correlation for fishing mortality deviations (only used when Fdev_model = "ar1")
  if("Fdev_rho" %in% names(starting_values)) input_list$par$Fdev_rho <- starting_values$Fdev_rho
  else input_list$par$Fdev_rho <- array(0, dim = c(input_list$data$n_regions, input_list$data$n_seas, input_list$data$n_fish_fleets))

  # Log mean fishing mortality
  if("ln_F_mean" %in% names(starting_values)) input_list$par$ln_F_mean <- starting_values$ln_F_mean
  else input_list$par$ln_F_mean <- array(log(0.1), dim = c(input_list$data$n_regions, input_list$data$n_seas, input_list$data$n_fish_fleets))

  # Log fishing deviations
  if("ln_F_devs" %in% names(starting_values)) input_list$par$ln_F_devs <- starting_values$ln_F_devs
  else input_list$par$ln_F_devs <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets))

  # Discard observation error
  if("ln_sigmaD" %in% names(starting_values)) input_list$par$ln_sigmaD <- starting_values$ln_sigmaD
  else input_list$par$ln_sigmaD <- array(log(0.01), dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets))

  if("ln_sigmaD_pop" %in% names(starting_values)) input_list$par$ln_sigmaD_pop <- starting_values$ln_sigmaD_pop
  else input_list$par$ln_sigmaD_pop <- array(log(0.01), dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets))

  # Process error discard deviations
  if("ln_sigma_dmr" %in% names(starting_values)) input_list$par$ln_sigma_dmr <- starting_values$ln_sigma_dmr
  else input_list$par$ln_sigma_dmr <- array(log(1), dim = c(input_list$data$n_regions, input_list$data$n_seas, input_list$data$n_fish_fleets))

  # Logit mean discard mortality
  if("logit_dmr_mean" %in% names(starting_values)) input_list$par$logit_dmr_mean <- starting_values$logit_dmr_mean
  else input_list$par$logit_dmr_mean <- array(0, dim = c(input_list$data$n_regions, input_list$data$n_seas, input_list$data$n_fish_fleets))

  # Logit discard mortality deviations
  if("logit_dmr_devs" %in% names(starting_values)) input_list$par$logit_dmr_devs <- starting_values$logit_dmr_devs
  else input_list$par$logit_dmr_devs <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets))

  # Mapping Options ---------------------------------------------------------

  # Retained Catch Stuff
  input_list <- do_sigmaC_mapping(input_list, sigmaC_spec)
  input_list <- do_sigmaC_pop_mapping(input_list, sigmaC_pop_spec)
  input_list <- do_sigmaF_mapping(input_list, sigmaF_spec)
  input_list <- do_Fdev_rho_mapping(input_list, Fdev_rho_spec)
  input_list <- do_Fmort_mapping(input_list)

  # Discard Catch Stuff
  input_list <- do_sigmaD_mapping(input_list, sigmaD_spec)
  input_list <- do_sigmaD_pop_mapping(input_list, sigmaD_pop_spec)
  input_list <- do_sigma_dmr_mapping(input_list, sigma_dmr_spec)
  input_list <- do_dmr_mean_mapping(input_list, dmr_mean_spec)
  input_list <- do_dmr_dev_mapping(input_list, dmr_dev_spec)

  # Print Messages ----------------------------------------------------------
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)
}

