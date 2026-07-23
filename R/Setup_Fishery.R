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

#' Map fishing mortality deviation parameters
#'
#' Constructs the \code{ln_F_devs} factor map, assigning unique estimation
#' indices to region–year–season–fleet cells where catch data are used
#' (\code{UseCatch == 1}) and mapping cells without catch data to \code{NA}.
#' This ensures that \code{ln_F_devs} parameters are only estimated for
#' dimensions with observed catch.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists. Requires \code{$data$UseCatch} and \code{$data$UseCatch_pop} to be populated by
#'   \code{\link{Setup_Mod_Catch_and_F}}.
#'
#' @return The input \code{input_list} with \code{$map$ln_F_devs} set to a
#'   factor vector. Cells with catch are assigned sequential integer indices;
#'   cells without catch are \code{NA}.
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

  # All valid sharing specs
  valid_specs <- c("fix", "est_all",
                   "est_shared_pop", "est_shared_r", "est_shared_y",
                   "est_shared_seas", "est_shared_f",
                   "est_shared_pop_r", "est_shared_pop_y", "est_shared_pop_seas",
                   "est_shared_pop_f", "est_shared_r_y", "est_shared_r_seas",
                   "est_shared_r_f", "est_shared_y_seas", "est_shared_y_f",
                   "est_shared_seas_f",
                   "est_shared_pop_r_y", "est_shared_pop_r_seas", "est_shared_pop_r_f",
                   "est_shared_pop_y_seas", "est_shared_pop_y_f", "est_shared_pop_seas_f",
                   "est_shared_r_y_seas", "est_shared_r_y_f", "est_shared_r_seas_f",
                   "est_shared_y_seas_f",
                   "est_shared_pop_r_y_seas", "est_shared_pop_r_y_f",
                   "est_shared_pop_r_seas_f", "est_shared_pop_y_seas_f",
                   "est_shared_r_y_seas_f",
                   "est_shared_pop_r_y_seas_f")

  if (!sigmaC_pop_spec %in% valid_specs)
    stop("sigmaC_pop_spec '", sigmaC_pop_spec,
         "' not recognized. Valid options: ", paste(valid_specs, collapse = ", "))

  # Fixing everything
  if (sigmaC_pop_spec == "fix") {
    input_list$map$ln_sigmaC_pop <- factor(rep(NA, length(input_list$par$ln_sigmaC_pop)))
    collect_message("sigmaC_pop is specified as: fix")
    return(input_list)
  }

  # Independently estimate everything
  if (sigmaC_pop_spec == "est_all") {
    input_list$map$ln_sigmaC_pop <- factor(seq_along(input_list$par$ln_sigmaC_pop))
    collect_message("sigmaC_pop is specified as: est_all")
    return(input_list)
  }

  # Pull dimensions out of the data list for readability
  n_pop         <- input_list$data$n_pop
  n_regions     <- input_list$data$n_regions
  n_years       <- length(input_list$data$years)
  n_seas        <- input_list$data$n_seas
  n_fish_fleets <- input_list$data$n_fish_fleets

  map_sigmaC_pop <- input_list$par$ln_sigmaC_pop
  counter <- 1  # running integer index for unique parameters

  # Walk every cell of the 5D array [p, r, y, seas, f]
  for (p in 1:n_pop) {
    for (r in 1:n_regions) {
      for (y in 1:n_years) {
        for (seas in 1:n_seas) {
          for (f in 1:n_fish_fleets) {

            # Check specification and see if dimensions > 1
            share_p    <- grepl("p",    sigmaC_pop_spec) && p    > 1
            share_r    <- grepl("_r",   sigmaC_pop_spec) && r    > 1
            share_y    <- grepl("_y",   sigmaC_pop_spec) && y    > 1
            share_seas <- grepl("seas", sigmaC_pop_spec) && seas > 1
            share_f    <- grepl("_f",   sigmaC_pop_spec) && f    > 1

            # Check to see if we are sharing a dimension, and if that dimesnion = 1
            is_first <- (!share_p    || p    == 1) &&
                        (!share_r    || r    == 1) &&
                        (!share_y    || y    == 1) &&
                        (!share_seas || seas == 1) &&
                        (!share_f    || f    == 1)

            if (is_first) {

              # figure out indexing, if sharing dimension, then use full vector, otherwise,
              # population unique index
              p_idx    <- if (share_p)    1:n_pop         else p
              r_idx    <- if (share_r)    1:n_regions     else r
              y_idx    <- if (share_y)    1:n_years       else y
              seas_idx <- if (share_seas) 1:n_seas        else seas
              f_idx    <- if (share_f)    1:n_fish_fleets else f
              map_sigmaC_pop[p_idx, r_idx, y_idx, seas_idx, f_idx] <- counter
              counter <- counter + 1  # update counter
            }  # end if
          } # end f
        } # end seas
      } # end y
    } # end r
  } # end p

  input_list$map$ln_sigmaC_pop <- factor(map_sigmaC_pop)
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

  # Sigma F -----------------------------------------------------------------
  map_sigma_dmr <- input_list$par$ln_sigma_dmr
  n_regions <- input_list$data$n_regions
  n_seas <- input_list$data$n_seas
  n_fish_fleets <- input_list$data$n_fish_fleets

  shared_specs <- c("est_shared_r", "est_shared_seas", "est_shared_f",
                    "est_shared_r_seas", "est_shared_r_f", "est_shared_seas_f",
                    "est_shared_r_seas_f")

  # In do_sigma_dmr_mapping, after defining shared_specs:
  valid_specs <- c(shared_specs, "fix", "est_all")
  if(!sigma_dmr_spec %in% valid_specs) stop("sigma_dmr_spec '", sigma_dmr_spec, "' not recognized. Valid options: ", paste(valid_specs, collapse = ", "))

  if(sigma_dmr_spec %in% shared_specs) {

    counter <- 1
    for(r in 1:n_regions) {
      for(seas in 1:n_seas) {
        for(f in 1:n_fish_fleets) {

          # --- Single dimension sharing ---
          if(sigma_dmr_spec == "est_shared_r" && r == 1) {
            map_sigma_dmr[, seas, f] <- counter; counter <- counter + 1
          }
          if(sigma_dmr_spec == "est_shared_seas" && seas == 1) {
            map_sigma_dmr[r, , f] <- counter; counter <- counter + 1
          }
          if(sigma_dmr_spec == "est_shared_f" && f == 1) {
            map_sigma_dmr[r, seas, ] <- counter; counter <- counter + 1
          }

          # --- Two dimension sharing ---
          if(sigma_dmr_spec == "est_shared_r_seas" && r == 1 && seas == 1) {
            map_sigma_dmr[, , f] <- counter; counter <- counter + 1
          }
          if(sigma_dmr_spec == "est_shared_r_f" && r == 1 && f == 1) {
            map_sigma_dmr[, seas, ] <- counter; counter <- counter + 1
          }
          if(sigma_dmr_spec == "est_shared_seas_f" && seas == 1 && f == 1) {
            map_sigma_dmr[r, , ] <- counter; counter <- counter + 1
          }

          # --- Three dimension sharing (single parameter) ---
          if(sigma_dmr_spec == "est_shared_r_seas_f" && r == 1 && seas == 1 && f == 1) {
            map_sigma_dmr[, , ] <- counter; counter <- counter + 1
          }

        } # end f
      } # end seas
    } # end r

    input_list$map$ln_sigma_dmr <- factor(map_sigma_dmr)
  }

  # Fixing sigma_dmr
  if(sigma_dmr_spec == "fix") input_list$map$ln_sigma_dmr <- factor(rep(NA, length(input_list$par$ln_sigma_dmr)))
  # Estimating all sigma_dmr
  if(sigma_dmr_spec == "est_all") input_list$map$ln_sigma_dmr <- factor(1:length(input_list$par$ln_sigma_dmr))

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

  # Sigma C -----------------------------------------------------------------
  map_sigmaD <- input_list$par$ln_sigmaD
  n_regions <- input_list$data$n_regions
  n_years <- length(input_list$data$years)
  n_seas <- input_list$data$n_seas
  n_fish_fleets <- input_list$data$n_fish_fleets

  shared_specs <- c("est_shared_r", "est_shared_y", "est_shared_seas", "est_shared_f",
                    "est_shared_r_y", "est_shared_r_seas", "est_shared_r_f",
                    "est_shared_y_seas", "est_shared_y_f", "est_shared_seas_f",
                    "est_shared_r_y_seas", "est_shared_r_y_f", "est_shared_r_seas_f", "est_shared_y_seas_f",
                    "est_shared_r_y_seas_f")

  # In do_sigmaF_mapping, after defining shared_specs:
  valid_specs <- c(shared_specs, "fix", "est_all")
  if(!sigmaD_spec %in% valid_specs) stop("sigmaD_spec '", sigmaD_spec, "' not recognized. Valid options: ", paste(valid_specs, collapse = ", "))

  if(sigmaD_spec %in% shared_specs) {

    counter <- 1
    for(r in 1:n_regions) {
      for(y in 1:n_years) {
        for(seas in 1:n_seas) {
          for(f in 1:n_fish_fleets) {

            # --- Single dimension sharing ---
            if(sigmaD_spec == "est_shared_r" && r == 1) {
              map_sigmaD[, y, seas, f] <- counter; counter <- counter + 1
            }
            if(sigmaD_spec == "est_shared_y" && y == 1) {
              map_sigmaD[r, , seas, f] <- counter; counter <- counter + 1
            }
            if(sigmaD_spec == "est_shared_seas" && seas == 1) {
              map_sigmaD[r, y, , f] <- counter; counter <- counter + 1
            }
            if(sigmaD_spec == "est_shared_f" && f == 1) {
              map_sigmaD[r, y, seas, ] <- counter; counter <- counter + 1
            }

            # --- Two dimension sharing ---
            if(sigmaD_spec == "est_shared_r_y" && r == 1 && y == 1) {
              map_sigmaD[, , seas, f] <- counter; counter <- counter + 1
            }
            if(sigmaD_spec == "est_shared_r_seas" && r == 1 && seas == 1) {
              map_sigmaD[, y, , f] <- counter; counter <- counter + 1
            }
            if(sigmaD_spec == "est_shared_r_f" && r == 1 && f == 1) {
              map_sigmaD[, y, seas, ] <- counter; counter <- counter + 1
            }
            if(sigmaD_spec == "est_shared_y_seas" && y == 1 && seas == 1) {
              map_sigmaD[r, , , f] <- counter; counter <- counter + 1
            }
            if(sigmaD_spec == "est_shared_y_f" && y == 1 && f == 1) {
              map_sigmaD[r, , seas, ] <- counter; counter <- counter + 1
            }
            if(sigmaD_spec == "est_shared_seas_f" && seas == 1 && f == 1) {
              map_sigmaD[r, y, , ] <- counter; counter <- counter + 1
            }

            # --- Three dimension sharing ---
            if(sigmaD_spec == "est_shared_r_y_seas" && r == 1 && y == 1 && seas == 1) {
              map_sigmaD[, , , f] <- counter; counter <- counter + 1
            }
            if(sigmaD_spec == "est_shared_r_y_f" && r == 1 && y == 1 && f == 1) {
              map_sigmaD[, , seas, ] <- counter; counter <- counter + 1
            }
            if(sigmaD_spec == "est_shared_r_seas_f" && r == 1 && seas == 1 && f == 1) {
              map_sigmaD[, y, , ] <- counter; counter <- counter + 1
            }
            if(sigmaD_spec == "est_shared_y_seas_f" && y == 1 && seas == 1 && f == 1) {
              map_sigmaD[r, , , ] <- counter; counter <- counter + 1
            }

            # --- Four dimension sharing (single parameter) ---
            if(sigmaD_spec == "est_shared_r_y_seas_f" && r == 1 && y == 1 && seas == 1 && f == 1) {
              map_sigmaD[, , , ] <- counter; counter <- counter + 1
            }

          } # end f
        } # end seas
      } # end y
    } # end r

    input_list$map$ln_sigmaD <- factor(map_sigmaD)
  }

  # Fixing sigmaD
  if(sigmaD_spec == "fix") input_list$map$ln_sigmaD <- factor(rep(NA, length(input_list$par$ln_sigmaD)))
  # Estimating all sigmaD
  if(sigmaD_spec == "est_all") input_list$map$ln_sigmaD <- factor(1:length(input_list$par$ln_sigmaD))

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

  # All valid sharing specs
  valid_specs <- c("fix", "est_all",
                   "est_shared_pop", "est_shared_r", "est_shared_y",
                   "est_shared_seas", "est_shared_f",
                   "est_shared_pop_r", "est_shared_pop_y", "est_shared_pop_seas",
                   "est_shared_pop_f", "est_shared_r_y", "est_shared_r_seas",
                   "est_shared_r_f", "est_shared_y_seas", "est_shared_y_f",
                   "est_shared_seas_f",
                   "est_shared_pop_r_y", "est_shared_pop_r_seas", "est_shared_pop_r_f",
                   "est_shared_pop_y_seas", "est_shared_pop_y_f", "est_shared_pop_seas_f",
                   "est_shared_r_y_seas", "est_shared_r_y_f", "est_shared_r_seas_f",
                   "est_shared_y_seas_f",
                   "est_shared_pop_r_y_seas", "est_shared_pop_r_y_f",
                   "est_shared_pop_r_seas_f", "est_shared_pop_y_seas_f",
                   "est_shared_r_y_seas_f",
                   "est_shared_pop_r_y_seas_f")

  if (!sigmaD_pop_spec %in% valid_specs)
    stop("sigmaD_pop_spec '", sigmaD_pop_spec,
         "' not recognized. Valid options: ", paste(valid_specs, collapse = ", "))

  # Fixing everything
  if (sigmaD_pop_spec == "fix") {
    input_list$map$ln_sigmaD_pop <- factor(rep(NA, length(input_list$par$ln_sigmaD_pop)))
    collect_message("sigmaD_pop is specified as: fix")
    return(input_list)
  }

  # Independently estimate everything
  if (sigmaD_pop_spec == "est_all") {
    input_list$map$ln_sigmaD_pop <- factor(seq_along(input_list$par$ln_sigmaD_pop))
    collect_message("sigmaD_pop is specified as: est_all")
    return(input_list)
  }

  # Pull dimensions out of the data list for readability
  n_pop         <- input_list$data$n_pop
  n_regions     <- input_list$data$n_regions
  n_years       <- length(input_list$data$years)
  n_seas        <- input_list$data$n_seas
  n_fish_fleets <- input_list$data$n_fish_fleets

  map_sigmaD_pop <- input_list$par$ln_sigmaD_pop
  counter <- 1  # running integer index for unique parameters

  # Walk every cell of the 5D array [p, r, y, seas, f]
  for (p in 1:n_pop) {
    for (r in 1:n_regions) {
      for (y in 1:n_years) {
        for (seas in 1:n_seas) {
          for (f in 1:n_fish_fleets) {

            # Check specification and see if dimensions > 1
            share_p    <- grepl("p",    sigmaD_pop_spec) && p    > 1
            share_r    <- grepl("_r",   sigmaD_pop_spec) && r    > 1
            share_y    <- grepl("_y",   sigmaD_pop_spec) && y    > 1
            share_seas <- grepl("seas", sigmaD_pop_spec) && seas > 1
            share_f    <- grepl("_f",   sigmaD_pop_spec) && f    > 1

            # Check to see if we are sharing a dimension, and if that dimesnion = 1
            is_first <- (!share_p    || p    == 1) &&
              (!share_r    || r    == 1) &&
              (!share_y    || y    == 1) &&
              (!share_seas || seas == 1) &&
              (!share_f    || f    == 1)

            if (is_first) {

              # figure out indexing, if sharing dimension, then use full vector, otherwise,
              # population unique index
              p_idx    <- if (share_p)    1:n_pop         else p
              r_idx    <- if (share_r)    1:n_regions     else r
              y_idx    <- if (share_y)    1:n_years       else y
              seas_idx <- if (share_seas) 1:n_seas        else seas
              f_idx    <- if (share_f)    1:n_fish_fleets else f
              map_sigmaD_pop[p_idx, r_idx, y_idx, seas_idx, f_idx] <- counter
              counter <- counter + 1  # update counter
            }  # end if
          } # end f
        } # end seas
      } # end y
    } # end r
  } # end p

  input_list$map$ln_sigmaD_pop <- factor(map_sigmaD_pop)
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

  map_dmr_mean <- input_list$par$logit_dmr_mean
  map_dmr_mean[] <- NA
  n_regions <- input_list$data$n_regions
  n_seas <- input_list$data$n_seas
  n_fish_fleets <- input_list$data$n_fish_fleets

  shared_specs <- c("est_shared_r", "est_shared_seas", "est_shared_f",
                    "est_shared_r_seas", "est_shared_r_f", "est_shared_seas_f",
                    "est_shared_r_seas_f")

  valid_specs <- c(shared_specs, "fix", "est_all")
  if(!dmr_mean_spec %in% valid_specs)
    stop("dmr_mean_spec '", dmr_mean_spec, "' not recognized. Valid options: ",
         paste(valid_specs, collapse = ", "))

  if(dmr_mean_spec == "fix") {
    input_list$map$logit_dmr_mean <- factor(rep(NA, length(input_list$par$logit_dmr_mean)))
    collect_message("dmr_mean is specified as: fix")
    return(input_list)
  }

  if(dmr_mean_spec == "est_all") {
    input_list$map$logit_dmr_mean <- factor(1:length(input_list$par$logit_dmr_mean))
    collect_message("dmr_mean is specified as: est_all")
    return(input_list)
  }

  counter <- 1
  for(r in 1:n_regions) {
    for(seas in 1:n_seas) {
      for(f in 1:n_fish_fleets) {

        if(dmr_mean_spec == "est_shared_r" && r == 1) {
          map_dmr_mean[, seas, f] <- counter; counter <- counter + 1
        }
        if(dmr_mean_spec == "est_shared_seas" && seas == 1) {
          map_dmr_mean[r, , f] <- counter; counter <- counter + 1
        }
        if(dmr_mean_spec == "est_shared_f" && f == 1) {
          map_dmr_mean[r, seas, ] <- counter; counter <- counter + 1
        }
        if(dmr_mean_spec == "est_shared_r_seas" && r == 1 && seas == 1) {
          map_dmr_mean[, , f] <- counter; counter <- counter + 1
        }
        if(dmr_mean_spec == "est_shared_r_f" && r == 1 && f == 1) {
          map_dmr_mean[, seas, ] <- counter; counter <- counter + 1
        }
        if(dmr_mean_spec == "est_shared_seas_f" && seas == 1 && f == 1) {
          map_dmr_mean[r, , ] <- counter; counter <- counter + 1
        }
        if(dmr_mean_spec == "est_shared_r_seas_f" && r == 1 && seas == 1 && f == 1) {
          map_dmr_mean[, , ] <- counter; counter <- counter + 1
        }

      }
    }
  }

  input_list$map$logit_dmr_mean <- factor(map_dmr_mean)
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
#'       \code{ln_F_devs},
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


#' Map fishery age composition overdispersion parameters
#'
#' Constructs factor maps for \code{ln_FishAge_theta} (fleet- region- and
#' sex-specific overdispersion) and \code{ln_FishAge_theta_agg} (aggregated
#' overdispersion) based on the composition type and likelihood specified in
#' \code{$data$FishAgeComps_Type} and \code{$data$FishAgeComps_LikeType}.
#' Parameters are mapped to \code{NA} for fleets using multinomial likelihoods
#' (\code{LikeType == 0}) or with no observed age compositions.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists. Requires \code{FishAgeComps_Type}, \code{FishAgeComps_LikeType},
#'   and \code{UseFishAgeComps} to be set by
#'   \code{\link{Setup_Mod_FishIdx_and_Comps}}.
#'
#' @return The input \code{input_list} with \code{$map$ln_FishAge_theta} and
#'   \code{$map$ln_FishAge_theta_agg} set to factor vectors. Active parameters
#'   receive sequential integer indices; inactive parameters are \code{NA}.
#'
#' @keywords internal
do_FishAge_theta_mapping <- function(input_list) {

  # setup counters
  counter_fishage_agg <- 1
  counter_fishage <- 1

  # initialize array to set up mapping
  map_FishAge_theta <- input_list$par$ln_FishAge_theta
  map_FishAge_theta_agg <- input_list$par$ln_FishAge_theta_agg
  map_FishAge_theta[] <- NA
  map_FishAge_theta_agg[] <- NA

  for(f in 1:input_list$data$n_fish_fleets) {

    # get unique fishery comp types
    fishage_comp_type <- unique(input_list$data$FishAgeComps_Type[,f])

    # If aggregated (ages)
    if(any(fishage_comp_type == 0) && input_list$data$FishAgeComps_LikeType[f] != 0) {
      map_FishAge_theta_agg[f] <- counter_fishage_agg
      counter_fishage_agg <- counter_fishage_agg + 1 # aggregated
    }

    # Loop through to make sure mapping stuff off correctly
    for(r in 1:input_list$data$n_regions) {

      # a region with no active comps for this fleet never contributes to the likelihood, so
      # its theta cell must stay NA (fixed) -- otherwise it's a free, unidentifiable parameter
      region_has_data <- sum(input_list$data$UseFishAgeComps[r,,,f]) > 0

      for(s in 1:input_list$data$n_sexes) {

        # if split by sex and region
        if(any(fishage_comp_type == 1) && input_list$data$FishAgeComps_LikeType[f] != 0 && region_has_data) {
          map_FishAge_theta[r,s,f] <- counter_fishage
          counter_fishage <- counter_fishage + 1 # split by sex and region
        }

        # joint by sex, split by region
        if(any(fishage_comp_type == 2) && input_list$data$FishAgeComps_LikeType[f] != 0 && s == 1 && region_has_data) {
          map_FishAge_theta[r,1,f] <- counter_fishage
          counter_fishage <- counter_fishage + 1 # joint by sex, split by region
        }

      } # end s loop
    } # end r loop

    # If we are using a multinomial or there aren't any age comps for a given fleet
    if(input_list$data$FishAgeComps_LikeType[f] == 0 || sum(input_list$data$UseFishAgeComps[,,,f]) == 0) {
      map_FishAge_theta[,,f] <- NA
      map_FishAge_theta_agg[f] <- NA
    }

  } # end f loop

  # Input into mapping list
  input_list$map$ln_FishAge_theta <- factor(map_FishAge_theta)
  input_list$map$ln_FishAge_theta_agg <- factor(map_FishAge_theta_agg)

  return(input_list)
}

#' Map population fishery age composition overdispersion parameters
#'
#' Constructs factor maps for \code{ln_FishAge_pop_theta} (fleet- region- and
#' sex-specific overdispersion) and \code{ln_FishAge_pop_theta_agg} (aggregated
#' overdispersion) based on the composition type and likelihood specified in
#' \code{$data$FishAgeComps_pop_Type} and \code{$data$FishAgeComps_pop_LikeType}.
#' Parameters are mapped to \code{NA} for fleets using multinomial likelihoods
#' (\code{LikeType == 0}) or with no observed age compositions.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists. Requires \code{FishAgeComps_pop_Type}, \code{FishAgeComps_pop_LikeType},
#'   and \code{UseFishAgeComps_pop} to be set by
#'   \code{\link{Setup_Mod_FishIdx_and_Comps}}.
#'
#' @return The input \code{input_list} with \code{$map$ln_FishAge_pop_theta} and
#'   \code{$map$ln_FishAge_pop_theta_agg} set to factor vectors. Active parameters
#'   receive sequential integer indices; inactive parameters are \code{NA}.
#'
#' @keywords internal
do_FishAge_pop_theta_mapping <- function(input_list) {

  # setup counters
  counter_fishage_agg <- 1
  counter_fishage <- 1

  # initialize array to set up mapping
  map_FishAge_pop_theta <- input_list$par$ln_FishAge_pop_theta
  map_FishAge_pop_theta_agg <- input_list$par$ln_FishAge_pop_theta_agg
  map_FishAge_pop_theta[] <- NA
  map_FishAge_pop_theta_agg[] <- NA

  for(p in 1:input_list$data$n_pop) {
    for(f in 1:input_list$data$n_fish_fleets) {

      # get unique fishery comp types
      fishage_comp_type <- unique(input_list$data$FishAgeComps_pop_Type[,f])

      # If aggregated (ages)
      if(any(fishage_comp_type == 0) && input_list$data$FishAgeComps_pop_LikeType[f] != 0) {
        map_FishAge_pop_theta_agg[p,f] <- counter_fishage_agg
        counter_fishage_agg <- counter_fishage_agg + 1 # aggregated
      }

      # Loop through to make sure mapping stuff off correctly
      for(r in 1:input_list$data$n_regions) {

        # a (pop, region) cell with no active comps for this fleet never contributes to the
        # likelihood, so its theta cell must stay NA (fixed) -- otherwise it's a free,
        # unidentifiable parameter
        region_has_data <- sum(input_list$data$UseFishAgeComps_pop[p,r,,,f]) > 0

        for(s in 1:input_list$data$n_sexes) {

          # if split by sex and region
          if(any(fishage_comp_type == 1) && input_list$data$FishAgeComps_pop_LikeType[f] != 0 && region_has_data) {
            map_FishAge_pop_theta[p,r,s,f] <- counter_fishage
            counter_fishage <- counter_fishage + 1 # split by sex and region
          }

          # joint by sex, split by region
          if(any(fishage_comp_type == 2) && input_list$data$FishAgeComps_pop_LikeType[f] != 0 && s == 1 && region_has_data) {
            map_FishAge_pop_theta[p,r,1,f] <- counter_fishage
            counter_fishage <- counter_fishage + 1 # joint by sex, split by region
          }

        } # end s loop
      } # end r loop

      # If we are using a multinomial or there aren't any age comps for a given fleet
      if(input_list$data$FishAgeComps_pop_LikeType[f] == 0 || sum(input_list$data$UseFishAgeComps_pop[p,,,,f]) == 0) {
        map_FishAge_pop_theta[p,,,f] <- NA
        map_FishAge_pop_theta_agg[p,f] <- NA
      }

    } # end f loop
  } # end p loop

  # Input into mapping list
  input_list$map$ln_FishAge_pop_theta <- factor(map_FishAge_pop_theta)
  input_list$map$ln_FishAge_pop_theta_agg <- factor(map_FishAge_pop_theta_agg)

  return(input_list)
}

#' Map fishery length composition overdispersion parameters
#'
#' Analogous to \code{\link{do_FishAge_theta_mapping}} but for length
#' compositions. Constructs factor maps for \code{ln_FishLen_theta} and
#' \code{ln_FishLen_theta_agg} based on \code{FishLenComps_Type},
#' \code{FishLenComps_LikeType}, and \code{UseFishLenComps}.
#' Parameters are mapped to \code{NA} for fleets using multinomial likelihoods
#' or with no observed length compositions.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists. Requires \code{FishLenComps_Type}, \code{FishLenComps_LikeType},
#'   and \code{UseFishLenComps} to be set by
#'   \code{\link{Setup_Mod_FishIdx_and_Comps}}.
#'
#' @return The input \code{input_list} with \code{$map$ln_FishLen_theta} and
#'   \code{$map$ln_FishLen_theta_agg} set to factor vectors.
#'
#' @keywords internal
do_FishLen_theta_mapping <- function(input_list) {

  # setup counters
  counter_fishlen_agg <- 1
  counter_fishlen <- 1

  # initialize array to set up mapping
  map_FishLen_theta <- input_list$par$ln_FishLen_theta
  map_FishLen_theta_agg <- input_list$par$ln_FishLen_theta_agg
  map_FishLen_theta[] <- NA
  map_FishLen_theta_agg[] <- NA

  for(f in 1:input_list$data$n_fish_fleets) {

    # get unique fishery comp types
    fishlen_comp_type <- unique(input_list$data$FishLenComps_Type[,f])

    # If aggregated (ages)
    if(any(fishlen_comp_type == 0) && input_list$data$FishLenComps_LikeType[f] != 0) {
      map_FishLen_theta_agg[f] <- counter_fishlen_agg
      counter_fishlen_agg <- counter_fishlen_agg + 1 # aggregated
    }

    # Loop through to make sure mapping stuff off correctly
    for(r in 1:input_list$data$n_regions) {
      for(s in 1:input_list$data$n_sexes) {

        # if split by sex and region
        if(any(fishlen_comp_type == 1) && input_list$data$FishLenComps_LikeType[f] != 0) {
          map_FishLen_theta[r,s,f] <- counter_fishlen
          counter_fishlen <- counter_fishlen + 1 # split by sex and region
        }

        # joint by sex, split by region
        if(any(fishlen_comp_type == 2) && input_list$data$FishLenComps_LikeType[f] != 0 && s == 1) {
          map_FishLen_theta[r,1,f] <- counter_fishlen
          counter_fishlen <- counter_fishlen + 1 # joint by sex, split by region
        }

      } # end s loop
    } # end r loop

    # If we are using a multinomial or there aren't any lenght comps for a given fleet
    if(input_list$data$FishLenComps_LikeType[f] == 0 || sum(input_list$data$UseFishLenComps[,,,f]) == 0) {
      map_FishLen_theta[,,f] <- NA
      map_FishLen_theta_agg[f] <- NA
    }

  } # end f loop

  # Input into mapping list
  input_list$map$ln_FishLen_theta <- factor(map_FishLen_theta)
  input_list$map$ln_FishLen_theta_agg <- factor(map_FishLen_theta_agg)

  return(input_list)
}

#' Map population fishery length composition overdispersion parameters
#'
#' Constructs factor maps for \code{ln_FishLen_pop_theta} (fleet- region- and
#' sex-specific overdispersion) and \code{ln_FishLen_pop_theta_agg} (aggregated
#' overdispersion) based on the composition type and likelihood specified in
#' \code{$data$FishLenComps_pop_Type} and \code{$data$FishLenComps_pop_LikeType}.
#' Parameters are mapped to \code{NA} for fleets using multinomial likelihoods
#' (\code{LikeType == 0}) or with no observed len compositions.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists. Requires \code{FishLenComps_pop_Type}, \code{FishLenComps_pop_LikeType},
#'   and \code{UseFishLenComps_pop} to be set by
#'   \code{\link{Setup_Mod_FishIdx_and_Comps}}.
#'
#' @return The input \code{input_list} with \code{$map$ln_FishLen_pop_theta} and
#'   \code{$map$ln_FishLen_pop_theta_agg} set to factor vectors. Active parameters
#'   receive sequential integer indices; inactive parameters are \code{NA}.
#'
#' @keywords internal
do_FishLen_pop_theta_mapping <- function(input_list) {

  # setup counters
  counter_fishlen_agg <- 1
  counter_fishlen <- 1

  # initialize array to set up mapping
  map_FishLen_pop_theta <- input_list$par$ln_FishLen_pop_theta
  map_FishLen_pop_theta_agg <- input_list$par$ln_FishLen_pop_theta_agg
  map_FishLen_pop_theta[] <- NA
  map_FishLen_pop_theta_agg[] <- NA

  for(p in 1:input_list$data$n_pop) {
    for(f in 1:input_list$data$n_fish_fleets) {

      # get unique fishery comp types
      fishlen_comp_type <- unique(input_list$data$FishLenComps_pop_Type[,f])

      # If aggregated (lens)
      if(any(fishlen_comp_type == 0) && input_list$data$FishLenComps_pop_LikeType[f] != 0) {
        map_FishLen_pop_theta_agg[p,f] <- counter_fishlen_agg
        counter_fishlen_agg <- counter_fishlen_agg + 1 # aggregated
      }

      # Loop through to make sure mapping stuff off correctly
      for(r in 1:input_list$data$n_regions) {
        for(s in 1:input_list$data$n_sexes) {

          # if split by sex and region
          if(any(fishlen_comp_type == 1) && input_list$data$FishLenComps_pop_LikeType[f] != 0) {
            map_FishLen_pop_theta[p,r,s,f] <- counter_fishlen
            counter_fishlen <- counter_fishlen + 1 # split by sex and region
          }

          # joint by sex, split by region
          if(any(fishlen_comp_type == 2) && input_list$data$FishLenComps_pop_LikeType[f] != 0 && s == 1) {
            map_FishLen_pop_theta[p,r,1,f] <- counter_fishlen
            counter_fishlen <- counter_fishlen + 1 # joint by sex, split by region
          }

        } # end s loop
      } # end r loop

      # If we are using a multinomial or there aren't any len comps for a given fleet
      if(input_list$data$FishLenComps_pop_LikeType[f] == 0 || sum(input_list$data$UseFishLenComps_pop[p,,,,f]) == 0) {
        map_FishLen_pop_theta[p,,,f] <- NA
        map_FishLen_pop_theta_agg[p,f] <- NA
      }

    } # end f loop
  } # end p loop

  # Input into mapping list
  input_list$map$ln_FishLen_pop_theta <- factor(map_FishLen_pop_theta)
  input_list$map$ln_FishLen_pop_theta_agg <- factor(map_FishLen_pop_theta_agg)

  return(input_list)
}

#' Map fishery age composition correlation parameters
#'
#' Constructs factor maps for \code{FishAge_corr_pars} (region- and sex-specific
#' AR1 and sex correlation parameters) and \code{FishAge_corr_pars_agg}
#' (aggregated correlation parameters) for 1D and 2D logistic-normal age
#' composition likelihoods.
#'
#' Parameters are activated only when
#' \code{FishAgeComps_LikeType} is in \code{c(3, 4)}.
#' These correspond to the 1D and 2D logistic-normal likelihoods.
#' All other likelihoods map correlation parameters to \code{NA}.
#'
#' For the 2D logistic-normal (\code{LikeType == 4}), both trailing elements of
#' the \code{[,,,,2]} slice are activated: element 1 for the age AR1 coefficient
#' and element 2 for the sex correlation (skipped when \code{n_sexes == 1}).
#'
#' @param input_list Named list containing \code{data}, \code{par}, and
#'   \code{map} components.
#'
#' @return The input \code{input_list} with elements
#'   \code{map\$FishAge_corr_pars} and
#'   \code{map\$FishAge_corr_pars_agg} set to factor vectors.
#'
#' @keywords internal
do_FishAge_corr_pars_mapping <- function(input_list) {

  # setup counters
  counter_fishage_corr <- 1
  counter_fishage_corr_agg <- 1

  # initialize array to set up mapping
  map_FishAge_corr_pars <- input_list$par$FishAge_corr_pars
  map_FishAge_corr_pars_agg <- input_list$par$FishAge_corr_pars_agg
  map_FishAge_corr_pars[] <- NA
  map_FishAge_corr_pars_agg[] <- NA

  for(f in 1:input_list$data$n_fish_fleets) {

    # No overdispersion parameters estimated
    if(input_list$data$FishAgeComps_LikeType[f] == 0 || sum(input_list$data$UseFishAgeComps[,,,f]) == 0) {
      map_FishAge_corr_pars[,,f,] <- NA
      map_FishAge_corr_pars_agg[f] <- NA
      next # skip if none
    }

    # get unique fishery comp types
    fishage_comp_type <- unique(input_list$data$FishAgeComps_Type[,f])

    # Aggregated Correlation Parameters
    if(any(fishage_comp_type == 0) && input_list$data$FishAgeComps_LikeType[f] != 0) {
      if(input_list$data$FishAgeComps_LikeType[f] == 3) {
        map_FishAge_corr_pars_agg[f] <- counter_fishage_corr_agg
        counter_fishage_corr_agg <- counter_fishage_corr_agg + 1 # aggregated
      }
    }

    # Loop through to make sure mapping stuff off correctly
    for(r in 1:input_list$data$n_regions) {
      for(s in 1:input_list$data$n_sexes) {

        # Split by region and sex
        if(any(fishage_comp_type == 1) && input_list$data$FishAgeComps_LikeType[f] != 0) {
          if(input_list$data$FishAgeComps_LikeType[f] == 3) {
            map_FishAge_corr_pars[r,s,f,1] <- counter_fishage_corr
            counter_fishage_corr <- counter_fishage_corr + 1
          }
        }

        # Joint by sex, split by region
        if(any(fishage_comp_type == 2) && input_list$data$FishAgeComps_LikeType[f] != 0 && s == 1) {

          # 1dar1 correlation
          if(input_list$data$FishAgeComps_LikeType[f] == 3) {
            map_FishAge_corr_pars[r,1,f,1] <- counter_fishage_corr
            counter_fishage_corr <- counter_fishage_corr + 1
          }

          # 2dar1 correlation
          if(input_list$data$FishAgeComps_LikeType[f] == 4) {
            for(i in 1:2) {
              if(i == 2 && input_list$data$n_sexes == 1) next # skip if we only have 1 sex
              map_FishAge_corr_pars[r,1,f,i] <- counter_fishage_corr
              counter_fishage_corr <- counter_fishage_corr + 1
            } # end i
          } # end if

        }
      } # end s loop
    } # end r loop

  } # end f loop

  # Input into mapping list
  input_list$map$FishAge_corr_pars_agg <- factor(map_FishAge_corr_pars_agg)
  input_list$map$FishAge_corr_pars <- factor(map_FishAge_corr_pars)

  return(input_list)
}

#' Map population fishery age composition correlation parameters
#'
#' Constructs factor maps for \code{FishAge_pop_corr_pars} (region- and sex-specific
#' AR1 and sex correlation parameters) and \code{FishAge_pop_corr_pars_agg}
#' (aggregated correlation parameters) for 1D and 2D logistic-normal age
#' composition likelihoods.
#'
#' Parameters are activated only when
#' \code{FishAgeComps_pop_LikeType} is in \code{c(3, 4)}.
#' These correspond to the 1D and 2D logistic-normal likelihoods.
#' All other likelihoods map correlation parameters to \code{NA}.
#'
#' For the 2D logistic-normal (\code{LikeType == 4}), both trailing elements of
#' the \code{[,,,,2]} slice are activated: element 1 for the age AR1 coefficient
#' and element 2 for the sex correlation (skipped when \code{n_sexes == 1}).
#'
#' @param input_list Named list containing \code{data}, \code{par}, and
#'   \code{map} components.
#'
#' @return The input \code{input_list} with elements
#'   \code{map\$FishAge_pop_corr_pars} and
#'   \code{map\$FishAge_pop_corr_pars_agg} set to factor vectors.
#'
#' @keywords internal
do_FishAge_pop_corr_pars_mapping <- function(input_list) {

  # setup counters
  counter_fishage_corr <- 1
  counter_fishage_corr_agg <- 1

  # initialize array to set up mapping
  map_FishAge_pop_corr_pars <- input_list$par$FishAge_pop_corr_pars
  map_FishAge_pop_corr_pars_agg <- input_list$par$FishAge_pop_corr_pars_agg
  map_FishAge_pop_corr_pars[] <- NA
  map_FishAge_pop_corr_pars_agg[] <- NA

  for(p in 1:input_list$data$n_pop) {
    for(f in 1:input_list$data$n_fish_fleets) {

      # No overdispersion parameters estimated
      if(input_list$data$FishAgeComps_pop_LikeType[f] == 0 || sum(input_list$data$UseFishAgeComps_pop[p,,,,f]) == 0) {
        map_FishAge_pop_corr_pars[p,,,f,] <- NA
        map_FishAge_pop_corr_pars_agg[p,f] <- NA
        next # skip if none
      }

      # get unique fishery comp types
      fishage_comp_type <- unique(input_list$data$FishAgeComps_pop_Type[,f])

      # Aggregated Correlation Parameters
      if(any(fishage_comp_type == 0) && input_list$data$FishAgeComps_pop_LikeType[f] != 0) {
        if(input_list$data$FishAgeComps_pop_LikeType[f] == 3) {
          map_FishAge_pop_corr_pars_agg[p,f] <- counter_fishage_corr_agg
          counter_fishage_corr_agg <- counter_fishage_corr_agg + 1 # aggregated
        }
      }

      # Loop through to make sure mapping stuff off correctly
      for(r in 1:input_list$data$n_regions) {
        for(s in 1:input_list$data$n_sexes) {

          # Split by region and sex
          if(any(fishage_comp_type == 1) && input_list$data$FishAgeComps_pop_LikeType[f] != 0) {
            if(input_list$data$FishAgeComps_pop_LikeType[f] == 3) {
              map_FishAge_pop_corr_pars[p,r,s,f,1] <- counter_fishage_corr
              counter_fishage_corr <- counter_fishage_corr + 1
            }
          }

          # Joint by sex, split by region
          if(any(fishage_comp_type == 2) && input_list$data$FishAgeComps_pop_LikeType[f] != 0 && s == 1) {

            # 1dar1 correlation
            if(input_list$data$FishAgeComps_pop_LikeType[f] == 3) {
              map_FishAge_pop_corr_pars[p,r,1,f,1] <- counter_fishage_corr
              counter_fishage_corr <- counter_fishage_corr + 1
            }

            # 2dar1 correlation
            if(input_list$data$FishAgeComps_pop_LikeType[f] == 4) {
              for(i in 1:2) {
                if(i == 2 && input_list$data$n_sexes == 1) next # skip if we only have 1 sex
                map_FishAge_pop_corr_pars[p,r,1,f,i] <- counter_fishage_corr
                counter_fishage_corr <- counter_fishage_corr + 1
              } # end i
            } # end if

          }
        } # end s loop
      } # end r loop

    } # end f loop
  } # end p loop

  # Input into mapping list
  input_list$map$FishAge_pop_corr_pars_agg <- factor(map_FishAge_pop_corr_pars_agg)
  input_list$map$FishAge_pop_corr_pars <- factor(map_FishAge_pop_corr_pars)

  return(input_list)
}

#' Map fishery length composition correlation parameters
#'
#' Analogous to \code{\link{do_FishAge_corr_pars_mapping}} but for length
#' compositions. Constructs factor maps for \code{FishLen_corr_pars} and
#' \code{FishLen_corr_pars_agg} for 1D and 2D logistic-normal length composition
#' likelihoods (\code{FishLenComps_LikeType} is in \code{c(3, 4)}).
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists.
#'
#' @return The input \code{input_list} with \code{$map$FishLen_corr_pars} and
#'   \code{$map$FishLen_corr_pars_agg} set to factor vectors.
#'
#' @keywords internal
do_FishLen_corr_pars_mapping <- function(input_list) {

  # setup counters
  counter_fishlen_corr <- 1
  counter_fishlen_corr_agg <- 1

  # initialize array to set up mapping
  map_FishLen_corr_pars <- input_list$par$FishLen_corr_pars
  map_FishLen_corr_pars_agg <- input_list$par$FishLen_corr_pars_agg
  map_FishLen_corr_pars[] <- NA
  map_FishLen_corr_pars_agg[] <- NA

  for(f in 1:input_list$data$n_fish_fleets) {

    # No overdispersion parameters estimated
    if(input_list$data$FishLenComps_LikeType[f] == 0 || sum(input_list$data$UseFishLenComps[,,,f]) == 0) {
      map_FishLen_corr_pars[,,f,] <- NA
      map_FishLen_corr_pars_agg[f] <- NA
      next # skip if none
    }

    # get unique fishery comp types
    fishlen_comp_type <- unique(input_list$data$FishLenComps_Type[,f])

    # Aggregated Correlation Parameters
    if(any(fishlen_comp_type == 0) && input_list$data$FishLenComps_LikeType[f] != 0) {
      if(input_list$data$FishLenComps_LikeType[f] == 3) {
        map_FishLen_corr_pars_agg[f] <- counter_fishlen_corr_agg
        counter_fishlen_corr_agg <- counter_fishlen_corr_agg + 1 # aggregated
      }
    }

    # Loop through to make sure mapping stuff off correctly
    for(r in 1:input_list$data$n_regions) {
      for(s in 1:input_list$data$n_sexes) {

        # Split by region and sex
        if(any(fishlen_comp_type == 1) && input_list$data$FishLenComps_LikeType[f] != 0) {
          if(input_list$data$FishLenComps_LikeType[f] == 3) {
            map_FishLen_corr_pars[r,s,f,1] <- counter_fishlen_corr
            counter_fishlen_corr <- counter_fishlen_corr + 1
          }
        }

        # Joint by sex, split by region
        if(any(fishlen_comp_type == 2) && input_list$data$FishLenComps_LikeType[f] != 0 && s == 1) {

          # 1dar1 correlation
          if(input_list$data$FishLenComps_LikeType[f] == 3) {
            map_FishLen_corr_pars[r,1,f,1] <- counter_fishlen_corr
            counter_fishlen_corr <- counter_fishlen_corr + 1
          }

          # 2dar1 correlation
          if(input_list$data$FishLenComps_LikeType[f] == 4) {
            for(i in 1:2) {
              if(i == 2 && input_list$data$n_sexes == 1) next # skip if we only have 1 sex
              map_FishLen_corr_pars[r,1,f,i] <- counter_fishlen_corr
              counter_fishlen_corr <- counter_fishlen_corr + 1
            } # end i
          } # end if

        }
      } # end s loop
    } # end r loop

  } # end f loop

  # Input into mapping list
  input_list$map$FishLen_corr_pars_agg <- factor(map_FishLen_corr_pars_agg)
  input_list$map$FishLen_corr_pars <- factor(map_FishLen_corr_pars)

  return(input_list)
}

#' Map population fishery length composition correlation parameters
#'
#' Constructs factor maps for \code{FishLen_pop_corr_pars} (region- and sex-specific
#' AR1 and sex correlation parameters) and \code{FishLen_pop_corr_pars_agg}
#' (aggregated correlation parameters) for 1D and 2D logistic-normal len
#' composition likelihoods.
#'
#' Parameters are activated only when
#' \code{FishLenComps_pop_LikeType} is in \code{c(3, 4)}.
#' These correspond to the 1D and 2D logistic-normal likelihoods.
#' All other likelihoods map correlation parameters to \code{NA}.
#'
#' For the 2D logistic-normal (\code{LikeType == 4}), both trailing elements of
#' the \code{[,,,,2]} slice are activated: element 1 for the len AR1 coefficient
#' and element 2 for the sex correlation (skipped when \code{n_sexes == 1}).
#'
#' @param input_list Named list containing \code{data}, \code{par}, and
#'   \code{map} components.
#'
#' @return The input \code{input_list} with elements
#'   \code{map\$FishLen_pop_corr_pars} and
#'   \code{map\$FishLen_pop_corr_pars_agg} set to factor vectors.
#'
#' @keywords internal
do_FishLen_pop_corr_pars_mapping <- function(input_list) {

  # setup counters
  counter_fishlen_corr <- 1
  counter_fishlen_corr_agg <- 1

  # initialize array to set up mapping
  map_FishLen_pop_corr_pars <- input_list$par$FishLen_pop_corr_pars
  map_FishLen_pop_corr_pars_agg <- input_list$par$FishLen_pop_corr_pars_agg
  map_FishLen_pop_corr_pars[] <- NA
  map_FishLen_pop_corr_pars_agg[] <- NA

  for(p in 1:input_list$data$n_pop) {
    for(f in 1:input_list$data$n_fish_fleets) {

      # No overdispersion parameters estimated
      if(input_list$data$FishLenComps_pop_LikeType[f] == 0 || sum(input_list$data$UseFishLenComps_pop[p,,,,f]) == 0) {
        map_FishLen_pop_corr_pars[p,,,f,] <- NA
        map_FishLen_pop_corr_pars_agg[p,f] <- NA
        next # skip if none
      }

      # get unique fishery comp types
      fishlen_comp_type <- unique(input_list$data$FishLenComps_pop_Type[,f])

      # Aggregated Correlation Parameters
      if(any(fishlen_comp_type == 0) && input_list$data$FishLenComps_pop_LikeType[f] != 0) {
        if(input_list$data$FishLenComps_pop_LikeType[f] == 3) {
          map_FishLen_pop_corr_pars_agg[p,f] <- counter_fishlen_corr_agg
          counter_fishlen_corr_agg <- counter_fishlen_corr_agg + 1 # aggregated
        }
      }

      # Loop through to make sure mapping stuff off correctly
      for(r in 1:input_list$data$n_regions) {
        for(s in 1:input_list$data$n_sexes) {

          # Split by region and sex
          if(any(fishlen_comp_type == 1) && input_list$data$FishLenComps_pop_LikeType[f] != 0) {
            if(input_list$data$FishLenComps_pop_LikeType[f] == 3) {
              map_FishLen_pop_corr_pars[p,r,s,f,1] <- counter_fishlen_corr
              counter_fishlen_corr <- counter_fishlen_corr + 1
            }
          }

          # Joint by sex, split by region
          if(any(fishlen_comp_type == 2) && input_list$data$FishLenComps_pop_LikeType[f] != 0 && s == 1) {

            # 1dar1 correlation
            if(input_list$data$FishLenComps_pop_LikeType[f] == 3) {
              map_FishLen_pop_corr_pars[p,r,1,f,1] <- counter_fishlen_corr
              counter_fishlen_corr <- counter_fishlen_corr + 1
            }

            # 2dar1 correlation
            if(input_list$data$FishLenComps_pop_LikeType[f] == 4) {
              for(i in 1:2) {
                if(i == 2 && input_list$data$n_sexes == 1) next # skip if we only have 1 sex
                map_FishLen_pop_corr_pars[p,r,1,f,i] <- counter_fishlen_corr
                counter_fishlen_corr <- counter_fishlen_corr + 1
              } # end i
            } # end if

          }
        } # end s loop
      } # end r loop

    } # end f loop
  } # end p loop

  # Input into mapping list
  input_list$map$FishLen_pop_corr_pars_agg <- factor(map_FishLen_pop_corr_pars_agg)
  input_list$map$FishLen_pop_corr_pars <- factor(map_FishLen_pop_corr_pars)

  return(input_list)
}

#' Map discarded fishery age composition overdispersion parameters
#'
#' Constructs factor maps for \code{ln_FishAge_discard_theta} (fleet- region- and
#' sex-specific overdispersion) and \code{ln_FishAge_discard_theta_agg} (aggregated
#' overdispersion) based on the composition type and likelihood specified in
#' \code{$data$FishAgeComps_discard_Type} and \code{$data$FishAgeComps_discard_LikeType}.
#' Parameters are mapped to \code{NA} for fleets using multinomial likelihoods
#' (\code{LikeType == 0}) or with no observed age compositions.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists. Requires \code{FishAgeComps_discard_Type}, \code{FishAgeComps_discard_LikeType},
#'   and \code{UseFishAgeComps_discard} to be set by
#'   \code{\link{Setup_Mod_FishIdx_and_Comps}}.
#'
#' @return The input \code{input_list} with \code{$map$ln_FishAge_discard_theta} and
#'   \code{$map$ln_FishAge_discard_theta_agg} set to factor vectors. Active parameters
#'   receive sequential integer indices; inactive parameters are \code{NA}.
#'
#' @keywords internal
do_FishAge_discard_theta_mapping <- function(input_list) {

  # setup counters
  counter_fishage_discard_agg <- 1
  counter_fishage_discard <- 1

  # initialize array to set up mapping
  map_FishAge_discard_theta <- input_list$par$ln_FishAge_discard_theta
  map_FishAge_discard_theta_agg <- input_list$par$ln_FishAge_discard_theta_agg
  map_FishAge_discard_theta[] <- NA
  map_FishAge_discard_theta_agg[] <- NA

  for(f in 1:input_list$data$n_fish_fleets) {

    # get unique fishery comp types
    fishage_discard_comp_type <- unique(input_list$data$FishAgeComps_discard_Type[,f])

    # If aggregated (ages)
    if(any(fishage_discard_comp_type == 0) && input_list$data$FishAgeComps_discard_LikeType[f] != 0) {
      map_FishAge_discard_theta_agg[f] <- counter_fishage_discard_agg
      counter_fishage_discard_agg <- counter_fishage_discard_agg + 1 # aggregated
    }

    # Loop through to make sure mapping stuff off correctly
    for(r in 1:input_list$data$n_regions) {
      for(s in 1:input_list$data$n_sexes) {

        # if split by sex and region
        if(any(fishage_discard_comp_type == 1) && input_list$data$FishAgeComps_discard_LikeType[f] != 0) {
          map_FishAge_discard_theta[r,s,f] <- counter_fishage_discard
          counter_fishage_discard <- counter_fishage_discard + 1 # split by sex and region
        }

        # joint by sex, split by region
        if(any(fishage_discard_comp_type == 2) && input_list$data$FishAgeComps_discard_LikeType[f] != 0 && s == 1) {
          map_FishAge_discard_theta[r,1,f] <- counter_fishage_discard
          counter_fishage_discard <- counter_fishage_discard + 1 # joint by sex, split by region
        }

      } # end s loop
    } # end r loop

    # If we are using a multinomial or there aren't any age comps for a given fleet
    if(input_list$data$FishAgeComps_discard_LikeType[f] == 0 || sum(input_list$data$UseFishAgeComps_discard[,,,f]) == 0) {
      map_FishAge_discard_theta[,,f] <- NA
      map_FishAge_discard_theta_agg[f] <- NA
    }

  } # end f loop

  # Input into mapping list
  input_list$map$ln_FishAge_discard_theta <- factor(map_FishAge_discard_theta)
  input_list$map$ln_FishAge_discard_theta_agg <- factor(map_FishAge_discard_theta_agg)

  return(input_list)
}

#' Map discarded population fishery age composition overdispersion parameters
#'
#' Constructs factor maps for \code{ln_FishAge_discard_pop_theta} (fleet- region- and
#' sex-specific overdispersion) and \code{ln_FishAge_discard_pop_theta_agg} (aggregated
#' overdispersion) based on the composition type and likelihood specified in
#' \code{$data$FishAgeComps_discard_pop_Type} and \code{$data$FishAgeComps_discard_pop_LikeType}.
#' Parameters are mapped to \code{NA} for fleets using multinomial likelihoods
#' (\code{LikeType == 0}) or with no observed age compositions.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists. Requires \code{FishAgeComps_discard_pop_Type}, \code{FishAgeComps_discard_pop_LikeType},
#'   and \code{UseFishAgeComps_discard_pop} to be set by
#'   \code{\link{Setup_Mod_FishIdx_and_Comps}}.
#'
#' @return The input \code{input_list} with \code{$map$ln_FishAge_discard_pop_theta} and
#'   \code{$map$ln_FishAge_discard_pop_theta_agg} set to factor vectors. Active parameters
#'   receive sequential integer indices; inactive parameters are \code{NA}.
#'
#' @keywords internal
do_FishAge_discard_pop_theta_mapping <- function(input_list) {

  # setup counters
  counter_fishage_discard_agg <- 1
  counter_fishage_discard <- 1

  # initialize array to set up mapping
  map_FishAge_discard_pop_theta <- input_list$par$ln_FishAge_discard_pop_theta
  map_FishAge_discard_pop_theta_agg <- input_list$par$ln_FishAge_discard_pop_theta_agg
  map_FishAge_discard_pop_theta[] <- NA
  map_FishAge_discard_pop_theta_agg[] <- NA

  for(p in 1:input_list$data$n_pop) {
    for(f in 1:input_list$data$n_fish_fleets) {

      # get unique fishery comp types
      fishage_discard_comp_type <- unique(input_list$data$FishAgeComps_discard_pop_Type[,f])

      # If aggregated (ages)
      if(any(fishage_discard_comp_type == 0) && input_list$data$FishAgeComps_discard_pop_LikeType[f] != 0) {
        map_FishAge_discard_pop_theta_agg[p,f] <- counter_fishage_discard_agg
        counter_fishage_discard_agg <- counter_fishage_discard_agg + 1 # aggregated
      }

      # Loop through to make sure mapping stuff off correctly
      for(r in 1:input_list$data$n_regions) {
        for(s in 1:input_list$data$n_sexes) {

          # if split by sex and region
          if(any(fishage_discard_comp_type == 1) && input_list$data$FishAgeComps_discard_pop_LikeType[f] != 0) {
            map_FishAge_discard_pop_theta[p,r,s,f] <- counter_fishage_discard
            counter_fishage_discard <- counter_fishage_discard + 1 # split by sex and region
          }

          # joint by sex, split by region
          if(any(fishage_discard_comp_type == 2) && input_list$data$FishAgeComps_discard_pop_LikeType[f] != 0 && s == 1) {
            map_FishAge_discard_pop_theta[p,r,1,f] <- counter_fishage_discard
            counter_fishage_discard <- counter_fishage_discard + 1 # joint by sex, split by region
          }

        } # end s loop
      } # end r loop

      # If we are using a multinomial or there aren't any age comps for a given fleet
      if(input_list$data$FishAgeComps_discard_pop_LikeType[f] == 0 || sum(input_list$data$UseFishAgeComps_discard_pop[p,,,,f]) == 0) {
        map_FishAge_discard_pop_theta[p,,,f] <- NA
        map_FishAge_discard_pop_theta_agg[p,f] <- NA
      }

    } # end f loop
  } # end p loop

  # Input into mapping list
  input_list$map$ln_FishAge_discard_pop_theta <- factor(map_FishAge_discard_pop_theta)
  input_list$map$ln_FishAge_discard_pop_theta_agg <- factor(map_FishAge_discard_pop_theta_agg)

  return(input_list)
}

#' Map discarded fishery length composition overdispersion parameters
#'
#' Constructs factor maps for \code{ln_FishLen_discard_theta} (fleet- region- and
#' sex-specific overdispersion) and \code{ln_FishLen_discard_theta_agg} (aggregated
#' overdispersion) based on the composition type and likelihood specified in
#' \code{$data$FishLenComps_discard_Type} and \code{$data$FishLenComps_discard_LikeType}.
#' Parameters are mapped to \code{NA} for fleets using multinomial likelihoods
#' (\code{LikeType == 0}) or with no observed len compositions.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists. Requires \code{FishLenComps_discard_Type}, \code{FishLenComps_discard_LikeType},
#'   and \code{UseFishLenComps_discard} to be set by
#'   \code{\link{Setup_Mod_FishIdx_and_Comps}}.
#'
#' @return The input \code{input_list} with \code{$map$ln_FishLen_discard_theta} and
#'   \code{$map$ln_FishLen_discard_theta_agg} set to factor vectors. Active parameters
#'   receive sequential integer indices; inactive parameters are \code{NA}.
#'
#' @keywords internal
do_FishLen_discard_theta_mapping <- function(input_list) {

  # setup counters
  counter_fishlen_discard_agg <- 1
  counter_fishlen_discard <- 1

  # initialize array to set up mapping
  map_FishLen_discard_theta <- input_list$par$ln_FishLen_discard_theta
  map_FishLen_discard_theta_agg <- input_list$par$ln_FishLen_discard_theta_agg
  map_FishLen_discard_theta[] <- NA
  map_FishLen_discard_theta_agg[] <- NA

  for(f in 1:input_list$data$n_fish_fleets) {

    # get unique fishery comp types
    fishlen_discard_comp_type <- unique(input_list$data$FishLenComps_discard_Type[,f])

    # If aggregated (lens)
    if(any(fishlen_discard_comp_type == 0) && input_list$data$FishLenComps_discard_LikeType[f] != 0) {
      map_FishLen_discard_theta_agg[f] <- counter_fishlen_discard_agg
      counter_fishlen_discard_agg <- counter_fishlen_discard_agg + 1 # aggregated
    }

    # Loop through to make sure mapping stuff off correctly
    for(r in 1:input_list$data$n_regions) {
      for(s in 1:input_list$data$n_sexes) {

        # if split by sex and region
        if(any(fishlen_discard_comp_type == 1) && input_list$data$FishLenComps_discard_LikeType[f] != 0) {
          map_FishLen_discard_theta[r,s,f] <- counter_fishlen_discard
          counter_fishlen_discard <- counter_fishlen_discard + 1 # split by sex and region
        }

        # joint by sex, split by region
        if(any(fishlen_discard_comp_type == 2) && input_list$data$FishLenComps_discard_LikeType[f] != 0 && s == 1) {
          map_FishLen_discard_theta[r,1,f] <- counter_fishlen_discard
          counter_fishlen_discard <- counter_fishlen_discard + 1 # joint by sex, split by region
        }

      } # end s loop
    } # end r loop

    # If we are using a multinomial or there aren't any len comps for a given fleet
    if(input_list$data$FishLenComps_discard_LikeType[f] == 0 || sum(input_list$data$UseFishLenComps_discard[,,,f]) == 0) {
      map_FishLen_discard_theta[,,f] <- NA
      map_FishLen_discard_theta_agg[f] <- NA
    }

  } # end f loop

  # Input into mapping list
  input_list$map$ln_FishLen_discard_theta <- factor(map_FishLen_discard_theta)
  input_list$map$ln_FishLen_discard_theta_agg <- factor(map_FishLen_discard_theta_agg)

  return(input_list)
}

#' Map discarded population fishery length composition overdispersion parameters
#'
#' Constructs factor maps for \code{ln_FishLen_discard_pop_theta} (fleet- region- and
#' sex-specific overdispersion) and \code{ln_FishLen_discard_pop_theta_agg} (aggregated
#' overdispersion) based on the composition type and likelihood specified in
#' \code{$data$FishLenComps_discard_pop_Type} and \code{$data$FishLenComps_discard_pop_LikeType}.
#' Parameters are mapped to \code{NA} for fleets using multinomial likelihoods
#' (\code{LikeType == 0}) or with no observed len compositions.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists. Requires \code{FishLenComps_discard_pop_Type}, \code{FishLenComps_discard_pop_LikeType},
#'   and \code{UseFishLenComps_discard_pop} to be set by
#'   \code{\link{Setup_Mod_FishIdx_and_Comps}}.
#'
#' @return The input \code{input_list} with \code{$map$ln_FishLen_discard_pop_theta} and
#'   \code{$map$ln_FishLen_discard_pop_theta_agg} set to factor vectors. Active parameters
#'   receive sequential integer indices; inactive parameters are \code{NA}.
#'
#' @keywords internal
do_FishLen_discard_pop_theta_mapping <- function(input_list) {

  # setup counters
  counter_fishlen_discard_agg <- 1
  counter_fishlen_discard <- 1

  # initialize array to set up mapping
  map_FishLen_discard_pop_theta <- input_list$par$ln_FishLen_discard_pop_theta
  map_FishLen_discard_pop_theta_agg <- input_list$par$ln_FishLen_discard_pop_theta_agg
  map_FishLen_discard_pop_theta[] <- NA
  map_FishLen_discard_pop_theta_agg[] <- NA

  for(p in 1:input_list$data$n_pop) {
    for(f in 1:input_list$data$n_fish_fleets) {

      # get unique fishery comp types
      fishlen_discard_comp_type <- unique(input_list$data$FishLenComps_discard_pop_Type[,f])

      # If aggregated (lens)
      if(any(fishlen_discard_comp_type == 0) && input_list$data$FishLenComps_discard_pop_LikeType[f] != 0) {
        map_FishLen_discard_pop_theta_agg[p,f] <- counter_fishlen_discard_agg
        counter_fishlen_discard_agg <- counter_fishlen_discard_agg + 1 # aggregated
      }

      # Loop through to make sure mapping stuff off correctly
      for(r in 1:input_list$data$n_regions) {
        for(s in 1:input_list$data$n_sexes) {

          # if split by sex and region
          if(any(fishlen_discard_comp_type == 1) && input_list$data$FishLenComps_discard_pop_LikeType[f] != 0) {
            map_FishLen_discard_pop_theta[p,r,s,f] <- counter_fishlen_discard
            counter_fishlen_discard <- counter_fishlen_discard + 1 # split by sex and region
          }

          # joint by sex, split by region
          if(any(fishlen_discard_comp_type == 2) && input_list$data$FishLenComps_discard_pop_LikeType[f] != 0 && s == 1) {
            map_FishLen_discard_pop_theta[p,r,1,f] <- counter_fishlen_discard
            counter_fishlen_discard <- counter_fishlen_discard + 1 # joint by sex, split by region
          }

        } # end s loop
      } # end r loop

      # If we are using a multinomial or there aren't any len comps for a given fleet
      if(input_list$data$FishLenComps_discard_pop_LikeType[f] == 0 || sum(input_list$data$UseFishLenComps_discard_pop[p,,,,f]) == 0) {
        map_FishLen_discard_pop_theta[p,,,f] <- NA
        map_FishLen_discard_pop_theta_agg[p,f] <- NA
      }

    } # end f loop
  } # end p loop

  # Input into mapping list
  input_list$map$ln_FishLen_discard_pop_theta <- factor(map_FishLen_discard_pop_theta)
  input_list$map$ln_FishLen_discard_pop_theta_agg <- factor(map_FishLen_discard_pop_theta_agg)

  return(input_list)
}


#' Map discarded fishery age composition correlation parameters
#'
#' Constructs factor maps for \code{FishAge_discard_corr_pars} (region- and sex-specific
#' AR1 and sex correlation parameters) and \code{FishAge_discard_corr_pars_agg}
#' (aggregated correlation parameters) for 1D and 2D logistic-normal age
#' composition likelihoods.
#'
#' Parameters are activated only when
#' \code{FishAgeComps_discard_LikeType} is in \code{c(3, 4)}.
#' These correspond to the 1D and 2D logistic-normal likelihoods.
#' All other likelihoods map correlation parameters to \code{NA}.
#'
#' For the 2D logistic-normal (\code{LikeType == 4}), both trailing elements of
#' the \code{[,,,,2]} slice are activated: element 1 for the age AR1 coefficient
#' and element 2 for the sex correlation (skipped when \code{n_sexes == 1}).
#'
#' @param input_list Named list containing \code{data}, \code{par}, and
#'   \code{map} components.
#'
#' @return The input \code{input_list} with elements
#'   \code{map\$FishAge_discard_corr_pars} and
#'   \code{map\$FishAge_discard_corr_pars_agg} set to factor vectors.
#'
#' @keywords internal
do_FishAge_discard_corr_pars_mapping <- function(input_list) {

  # setup counters
  counter_fishage_discard_corr <- 1
  counter_fishage_discard_corr_agg <- 1

  # initialize array to set up mapping
  map_FishAge_discard_corr_pars <- input_list$par$FishAge_discard_corr_pars
  map_FishAge_discard_corr_pars_agg <- input_list$par$FishAge_discard_corr_pars_agg
  map_FishAge_discard_corr_pars[] <- NA
  map_FishAge_discard_corr_pars_agg[] <- NA

  for(f in 1:input_list$data$n_fish_fleets) {

    # No overdispersion parameters estimated
    if(input_list$data$FishAgeComps_discard_LikeType[f] == 0 || sum(input_list$data$UseFishAgeComps_discard[,,,f]) == 0) {
      map_FishAge_discard_corr_pars[,,f,] <- NA
      map_FishAge_discard_corr_pars_agg[f] <- NA
      next # skip if none
    }

    # get unique fishery comp types
    fishage_discard_comp_type <- unique(input_list$data$FishAgeComps_discard_Type[,f])

    # Aggregated Correlation Parameters
    if(any(fishage_discard_comp_type == 0) && input_list$data$FishAgeComps_discard_LikeType[f] != 0) {
      if(input_list$data$FishAgeComps_discard_LikeType[f] == 3) {
        map_FishAge_discard_corr_pars_agg[f] <- counter_fishage_discard_corr_agg
        counter_fishage_discard_corr_agg <- counter_fishage_discard_corr_agg + 1 # aggregated
      }
    }

    # Loop through to make sure mapping stuff off correctly
    for(r in 1:input_list$data$n_regions) {
      for(s in 1:input_list$data$n_sexes) {

        # Split by region and sex
        if(any(fishage_discard_comp_type == 1) && input_list$data$FishAgeComps_discard_LikeType[f] != 0) {
          if(input_list$data$FishAgeComps_discard_LikeType[f] == 3) {
            map_FishAge_discard_corr_pars[r,s,f,1] <- counter_fishage_discard_corr
            counter_fishage_discard_corr <- counter_fishage_discard_corr + 1
          }
        }

        # Joint by sex, split by region
        if(any(fishage_discard_comp_type == 2) && input_list$data$FishAgeComps_discard_LikeType[f] != 0 && s == 1) {

          # 1dar1 correlation
          if(input_list$data$FishAgeComps_discard_LikeType[f] == 3) {
            map_FishAge_discard_corr_pars[r,1,f,1] <- counter_fishage_discard_corr
            counter_fishage_discard_corr <- counter_fishage_discard_corr + 1
          }

          # 2dar1 correlation
          if(input_list$data$FishAgeComps_discard_LikeType[f] == 4) {
            for(i in 1:2) {
              if(i == 2 && input_list$data$n_sexes == 1) next # skip if we only have 1 sex
              map_FishAge_discard_corr_pars[r,1,f,i] <- counter_fishage_discard_corr
              counter_fishage_discard_corr <- counter_fishage_discard_corr + 1
            } # end i
          } # end if

        }
      } # end s loop
    } # end r loop

  } # end f loop

  # Input into mapping list
  input_list$map$FishAge_discard_corr_pars_agg <- factor(map_FishAge_discard_corr_pars_agg)
  input_list$map$FishAge_discard_corr_pars <- factor(map_FishAge_discard_corr_pars)

  return(input_list)
}

#' Map discarded population fishery age composition correlation parameters
#'
#' Constructs factor maps for \code{FishAge_discard_pop_corr_pars} (region- and sex-specific
#' AR1 and sex correlation parameters) and \code{FishAge_discard_pop_corr_pars_agg}
#' (aggregated correlation parameters) for 1D and 2D logistic-normal age
#' composition likelihoods.
#'
#' Parameters are activated only when
#' \code{FishAgeComps_discard_pop_LikeType} is in \code{c(3, 4)}.
#' These correspond to the 1D and 2D logistic-normal likelihoods.
#' All other likelihoods map correlation parameters to \code{NA}.
#'
#' For the 2D logistic-normal (\code{LikeType == 4}), both trailing elements of
#' the \code{[,,,,2]} slice are activated: element 1 for the age AR1 coefficient
#' and element 2 for the sex correlation (skipped when \code{n_sexes == 1}).
#'
#' @param input_list Named list containing \code{data}, \code{par}, and
#'   \code{map} components.
#'
#' @return The input \code{input_list} with elements
#'   \code{map\$FishAge_discard_pop_corr_pars} and
#'   \code{map\$FishAge_discard_pop_corr_pars_agg} set to factor vectors.
#'
#' @keywords internal
do_FishAge_discard_pop_corr_pars_mapping <- function(input_list) {

  # setup counters
  counter_fishage_discard_corr <- 1
  counter_fishage_discard_corr_agg <- 1

  # initialize array to set up mapping
  map_FishAge_discard_pop_corr_pars <- input_list$par$FishAge_discard_pop_corr_pars
  map_FishAge_discard_pop_corr_pars_agg <- input_list$par$FishAge_discard_pop_corr_pars_agg
  map_FishAge_discard_pop_corr_pars[] <- NA
  map_FishAge_discard_pop_corr_pars_agg[] <- NA

  for(p in 1:input_list$data$n_pop) {
    for(f in 1:input_list$data$n_fish_fleets) {

      # No overdispersion parameters estimated
      if(input_list$data$FishAgeComps_discard_pop_LikeType[f] == 0 || sum(input_list$data$UseFishAgeComps_discard_pop[p,,,,f]) == 0) {
        map_FishAge_discard_pop_corr_pars[p,,,f,] <- NA
        map_FishAge_discard_pop_corr_pars_agg[p,f] <- NA
        next # skip if none
      }

      # get unique fishery comp types
      fishage_discard_comp_type <- unique(input_list$data$FishAgeComps_discard_pop_Type[,f])

      # Aggregated Correlation Parameters
      if(any(fishage_discard_comp_type == 0) && input_list$data$FishAgeComps_discard_pop_LikeType[f] != 0) {
        if(input_list$data$FishAgeComps_discard_pop_LikeType[f] == 3) {
          map_FishAge_discard_pop_corr_pars_agg[p,f] <- counter_fishage_discard_corr_agg
          counter_fishage_discard_corr_agg <- counter_fishage_discard_corr_agg + 1 # aggregated
        }
      }

      # Loop through to make sure mapping stuff off correctly
      for(r in 1:input_list$data$n_regions) {
        for(s in 1:input_list$data$n_sexes) {

          # Split by region and sex
          if(any(fishage_discard_comp_type == 1) && input_list$data$FishAgeComps_discard_pop_LikeType[f] != 0) {
            if(input_list$data$FishAgeComps_discard_pop_LikeType[f] == 3) {
              map_FishAge_discard_pop_corr_pars[p,r,s,f,1] <- counter_fishage_discard_corr
              counter_fishage_discard_corr <- counter_fishage_discard_corr + 1
            }
          }

          # Joint by sex, split by region
          if(any(fishage_discard_comp_type == 2) && input_list$data$FishAgeComps_discard_pop_LikeType[f] != 0 && s == 1) {

            # 1dar1 correlation
            if(input_list$data$FishAgeComps_discard_pop_LikeType[f] == 3) {
              map_FishAge_discard_pop_corr_pars[p,r,1,f,1] <- counter_fishage_discard_corr
              counter_fishage_discard_corr <- counter_fishage_discard_corr + 1
            }

            # 2dar1 correlation
            if(input_list$data$FishAgeComps_discard_pop_LikeType[f] == 4) {
              for(i in 1:2) {
                if(i == 2 && input_list$data$n_sexes == 1) next # skip if we only have 1 sex
                map_FishAge_discard_pop_corr_pars[p,r,1,f,i] <- counter_fishage_discard_corr
                counter_fishage_discard_corr <- counter_fishage_discard_corr + 1
              } # end i
            } # end if

          }
        } # end s loop
      } # end r loop

    } # end f loop
  } # end p loop

  # Input into mapping list
  input_list$map$FishAge_discard_pop_corr_pars_agg <- factor(map_FishAge_discard_pop_corr_pars_agg)
  input_list$map$FishAge_discard_pop_corr_pars <- factor(map_FishAge_discard_pop_corr_pars)

  return(input_list)
}


#' Map discarded fishery length composition correlation parameters
#'
#' Constructs factor maps for \code{FishLen_discard_corr_pars} (region- and sex-specific
#' AR1 and sex correlation parameters) and \code{FishLen_discard_corr_pars_agg}
#' (aggregated correlation parameters) for 1D and 2D logistic-normal len
#' composition likelihoods.
#'
#' Parameters are activated only when
#' \code{FishLenComps_discard_LikeType} is in \code{c(3, 4)}.
#' These correspond to the 1D and 2D logistic-normal likelihoods.
#' All other likelihoods map correlation parameters to \code{NA}.
#'
#' For the 2D logistic-normal (\code{LikeType == 4}), both trailing elements of
#' the \code{[,,,,2]} slice are activated: element 1 for the len AR1 coefficient
#' and element 2 for the sex correlation (skipped when \code{n_sexes == 1}).
#'
#' @param input_list Named list containing \code{data}, \code{par}, and
#'   \code{map} components.
#'
#' @return The input \code{input_list} with elements
#'   \code{map\$FishLen_discard_corr_pars} and
#'   \code{map\$FishLen_discard_corr_pars_agg} set to factor vectors.
#'
#' @keywords internal
do_FishLen_discard_corr_pars_mapping <- function(input_list) {

  # setup counters
  counter_fishlen_discard_corr <- 1
  counter_fishlen_discard_corr_agg <- 1

  # initialize array to set up mapping
  map_FishLen_discard_corr_pars <- input_list$par$FishLen_discard_corr_pars
  map_FishLen_discard_corr_pars_agg <- input_list$par$FishLen_discard_corr_pars_agg
  map_FishLen_discard_corr_pars[] <- NA
  map_FishLen_discard_corr_pars_agg[] <- NA

  for(f in 1:input_list$data$n_fish_fleets) {

    # No overdispersion parameters estimated
    if(input_list$data$FishLenComps_discard_LikeType[f] == 0 || sum(input_list$data$UseFishLenComps_discard[,,,f]) == 0) {
      map_FishLen_discard_corr_pars[,,f,] <- NA
      map_FishLen_discard_corr_pars_agg[f] <- NA
      next # skip if none
    }

    # get unique fishery comp types
    fishlen_discard_comp_type <- unique(input_list$data$FishLenComps_discard_Type[,f])

    # Aggregated Correlation Parameters
    if(any(fishlen_discard_comp_type == 0) && input_list$data$FishLenComps_discard_LikeType[f] != 0) {
      if(input_list$data$FishLenComps_discard_LikeType[f] == 3) {
        map_FishLen_discard_corr_pars_agg[f] <- counter_fishlen_discard_corr_agg
        counter_fishlen_discard_corr_agg <- counter_fishlen_discard_corr_agg + 1 # aggregated
      }
    }

    # Loop through to make sure mapping stuff off correctly
    for(r in 1:input_list$data$n_regions) {
      for(s in 1:input_list$data$n_sexes) {

        # Split by region and sex
        if(any(fishlen_discard_comp_type == 1) && input_list$data$FishLenComps_discard_LikeType[f] != 0) {
          if(input_list$data$FishLenComps_discard_LikeType[f] == 3) {
            map_FishLen_discard_corr_pars[r,s,f,1] <- counter_fishlen_discard_corr
            counter_fishlen_discard_corr <- counter_fishlen_discard_corr + 1
          }
        }

        # Joint by sex, split by region
        if(any(fishlen_discard_comp_type == 2) && input_list$data$FishLenComps_discard_LikeType[f] != 0 && s == 1) {

          # 1dar1 correlation
          if(input_list$data$FishLenComps_discard_LikeType[f] == 3) {
            map_FishLen_discard_corr_pars[r,1,f,1] <- counter_fishlen_discard_corr
            counter_fishlen_discard_corr <- counter_fishlen_discard_corr + 1
          }

          # 2dar1 correlation
          if(input_list$data$FishLenComps_discard_LikeType[f] == 4) {
            for(i in 1:2) {
              if(i == 2 && input_list$data$n_sexes == 1) next # skip if we only have 1 sex
              map_FishLen_discard_corr_pars[r,1,f,i] <- counter_fishlen_discard_corr
              counter_fishlen_discard_corr <- counter_fishlen_discard_corr + 1
            } # end i
          } # end if

        }
      } # end s loop
    } # end r loop

  } # end f loop

  # Input into mapping list
  input_list$map$FishLen_discard_corr_pars_agg <- factor(map_FishLen_discard_corr_pars_agg)
  input_list$map$FishLen_discard_corr_pars <- factor(map_FishLen_discard_corr_pars)

  return(input_list)
}

#' Map discarded population fishery length composition correlation parameters
#'
#' Constructs factor maps for \code{FishLen_discard_pop_corr_pars} (region- and sex-specific
#' AR1 and sex correlation parameters) and \code{FishLen_discard_pop_corr_pars_agg}
#' (aggregated correlation parameters) for 1D and 2D logistic-normal len
#' composition likelihoods.
#'
#' Parameters are activated only when
#' \code{FishLenComps_discard_pop_LikeType} is in \code{c(3, 4)}.
#' These correspond to the 1D and 2D logistic-normal likelihoods.
#' All other likelihoods map correlation parameters to \code{NA}.
#'
#' For the 2D logistic-normal (\code{LikeType == 4}), both trailing elements of
#' the \code{[,,,,2]} slice are activated: element 1 for the len AR1 coefficient
#' and element 2 for the sex correlation (skipped when \code{n_sexes == 1}).
#'
#' @param input_list Named list containing \code{data}, \code{par}, and
#'   \code{map} components.
#'
#' @return The input \code{input_list} with elements
#'   \code{map\$FishLen_discard_pop_corr_pars} and
#'   \code{map\$FishLen_discard_pop_corr_pars_agg} set to factor vectors.
#'
#' @keywords internal
do_FishLen_discard_pop_corr_pars_mapping <- function(input_list) {

  # setup counters
  counter_fishlen_discard_corr <- 1
  counter_fishlen_discard_corr_agg <- 1

  # initialize array to set up mapping
  map_FishLen_discard_pop_corr_pars <- input_list$par$FishLen_discard_pop_corr_pars
  map_FishLen_discard_pop_corr_pars_agg <- input_list$par$FishLen_discard_pop_corr_pars_agg
  map_FishLen_discard_pop_corr_pars[] <- NA
  map_FishLen_discard_pop_corr_pars_agg[] <- NA

  for(p in 1:input_list$data$n_pop) {
    for(f in 1:input_list$data$n_fish_fleets) {

      # No overdispersion parameters estimated
      if(input_list$data$FishLenComps_discard_pop_LikeType[f] == 0 || sum(input_list$data$UseFishLenComps_discard_pop[p,,,,f]) == 0) {
        map_FishLen_discard_pop_corr_pars[p,,,f,] <- NA
        map_FishLen_discard_pop_corr_pars_agg[p,f] <- NA
        next # skip if none
      }

      # get unique fishery comp types
      fishlen_discard_comp_type <- unique(input_list$data$FishLenComps_discard_pop_Type[,f])

      # Aggregated Correlation Parameters
      if(any(fishlen_discard_comp_type == 0) && input_list$data$FishLenComps_discard_pop_LikeType[f] != 0) {
        if(input_list$data$FishLenComps_discard_pop_LikeType[f] == 3) {
          map_FishLen_discard_pop_corr_pars_agg[p,f] <- counter_fishlen_discard_corr_agg
          counter_fishlen_discard_corr_agg <- counter_fishlen_discard_corr_agg + 1 # aggregated
        }
      }

      # Loop through to make sure mapping stuff off correctly
      for(r in 1:input_list$data$n_regions) {
        for(s in 1:input_list$data$n_sexes) {

          # Split by region and sex
          if(any(fishlen_discard_comp_type == 1) && input_list$data$FishLenComps_discard_pop_LikeType[f] != 0) {
            if(input_list$data$FishLenComps_discard_pop_LikeType[f] == 3) {
              map_FishLen_discard_pop_corr_pars[p,r,s,f,1] <- counter_fishlen_discard_corr
              counter_fishlen_discard_corr <- counter_fishlen_discard_corr + 1
            }
          }

          # Joint by sex, split by region
          if(any(fishlen_discard_comp_type == 2) && input_list$data$FishLenComps_discard_pop_LikeType[f] != 0 && s == 1) {

            # 1dar1 correlation
            if(input_list$data$FishLenComps_discard_pop_LikeType[f] == 3) {
              map_FishLen_discard_pop_corr_pars[p,r,1,f,1] <- counter_fishlen_discard_corr
              counter_fishlen_discard_corr <- counter_fishlen_discard_corr + 1
            }

            # 2dar1 correlation
            if(input_list$data$FishLenComps_discard_pop_LikeType[f] == 4) {
              for(i in 1:2) {
                if(i == 2 && input_list$data$n_sexes == 1) next # skip if we only have 1 sex
                map_FishLen_discard_pop_corr_pars[p,r,1,f,i] <- counter_fishlen_discard_corr
                counter_fishlen_discard_corr <- counter_fishlen_discard_corr + 1
              } # end i
            } # end if

          }
        } # end s loop
      } # end r loop

    } # end f loop
  } # end p loop

  # Input into mapping list
  input_list$map$FishLen_discard_pop_corr_pars_agg <- factor(map_FishLen_discard_pop_corr_pars_agg)
  input_list$map$FishLen_discard_pop_corr_pars <- factor(map_FishLen_discard_pop_corr_pars)

  return(input_list)
}

#' Set up discard age and length composition inputs
#'
#' @param input_list Main model input list containing data, parameters, and mapping structures.
#'
#' @param ObsFishAgeComps_discard 5D array of observed discard age compositions:
#'   \code{[n_regions, n_years, n_seas, n_sexes, n_fish_fleets]}
#'
#' @param UseFishAgeComps_discard Matrix indicating whether discard age comps are used:
#'   \code{[n_regions, n_years, n_seas, n_fish_fleets]}
#'
#' @param ISS_FishAgeComps_discard Optional ISS (effective sample size) array for discard age comps:
#'   \code{[n_regions, n_years, n_seas, n_sexes, n_fish_fleets]}
#'
#' @param ObsFishLenComps_discard 6D array of observed discard length compositions:
#'   \code{[n_regions, n_years, n_seas, n_lens, n_sexes, n_fish_fleets]}
#'
#' @param UseFishLenComps_discard Matrix indicating whether discard length comps are used:
#'   \code{[n_regions, n_years, n_seas, n_fish_fleets]}
#'
#' @param ISS_FishLenComps_discard Optional ISS array for discard length comps:
#'   \code{[n_regions, n_years, n_seas, n_sexes, n_fish_fleets]}
#'
#' @param FishAgeComps_discard_LikeType Character vector (length n_fish_fleets) specifying likelihood type:
#'   one of \code{c("none","Multinomial","Dirichlet-Multinomial","iid-Logistic-Normal","1d-Logistic-Normal","2d-Logistic-Normal")}
#'
#' @param FishLenComps_discard_LikeType Character vector (length n_fish_fleets) specifying likelihood type
#'
#' @param FishAgeComps_discard_Type List/encoded character strings defining composition structure by year and fleet.
#'
#' @param FishLenComps_discard_Type List/encoded character strings defining composition structure by year and fleet.
#'
#' @param ObsFishAgeComps_discard_pop 6D array of population-specific discard age comps:
#'   \code{[n_pop, n_regions, n_years, n_seas, n_sexes, n_fish_fleets]}
#'
#' @param UseFishAgeComps_discard_pop 5D array indicating use of population age comps:
#'   \code{[n_pop, n_regions, n_years, n_seas, n_fish_fleets]}
#'
#' @param ISS_FishAgeComps_discard_pop Optional ISS array for population discard age comps:
#'   \code{[n_pop, n_regions, n_years, n_seas, n_sexes, n_fish_fleets]}
#'
#' @param ObsFishLenComps_discard_pop 7D array of population-specific discard length comps:
#'   \code{[n_pop, n_regions, n_years, n_seas, n_lens, n_sexes, n_fish_fleets]}
#'
#' @param UseFishLenComps_discard_pop 5D array indicating use of population length comps:
#'   \code{[n_pop, n_regions, n_years, n_seas, n_fish_fleets]}
#'
#' @param ISS_FishLenComps_discard_pop Optional ISS array for population length comps:
#'   \code{[n_pop, n_regions, n_years, n_seas, n_sexes, n_fish_fleets]}
#'
#' @param FishAgeComps_discard_pop_LikeType Character vector (length n_fish_fleets) for population age likelihood types.
#'
#' @param FishLenComps_discard_pop_LikeType Character vector (length n_fish_fleets) for population length likelihood types.
#'
#' @param FishAgeComps_discard_pop_Type Encoded structure definitions for population age comps by year/fleet.
#'
#' @param FishLenComps_discard_pop_Type Encoded structure definitions for population length comps by year/fleet.
#'
#' @param ... Optional starting values for parameters (e.g., dispersion, correlation)
#'
#' @return The input \code{input_list} updated with:
#' \itemize{
#'   \item discard age and length composition data structures
#'   \item ISS (effective sample size) arrays (computed or supplied)
#'   \item likelihood type mappings (integer-coded)
#'   \item composition type matrices by year and fleet
#'   \item population-specific discard composition structures
#'   \item parameter arrays for dispersion and correlation
#'   \item mapping configurations for estimation
#' }
#'
#' @details
#' All composition arrays follow consistent indexing conventions:
#' \itemize{
#'   \item Age compositions: \code{[region, year, season, sex, fleet]}
#'   \item Length compositions: \code{[region, year, season, length, sex, fleet]}
#'   \item Population age comps: \code{[pop, region, year, season, sex, fleet]}
#'   \item Population length comps: \code{[pop, region, year, season, length, sex, fleet]}
#' }
#'
#'
#' @keywords internal
#' @importFrom stringr str_detect
#' @family Model Setup
Setup_Mod_Discard_Comps     <- function(input_list,
                                        ObsFishAgeComps_discard,
                                        UseFishAgeComps_discard,
                                        ISS_FishAgeComps_discard,
                                        ObsFishLenComps_discard,
                                        UseFishLenComps_discard,
                                        ISS_FishLenComps_discard,
                                        FishAgeComps_discard_LikeType,
                                        FishLenComps_discard_LikeType,
                                        FishAgeComps_discard_Type,
                                        FishLenComps_discard_Type,
                                        ObsFishAgeComps_discard_pop,
                                        UseFishAgeComps_discard_pop,
                                        ISS_FishAgeComps_discard_pop,
                                        ObsFishLenComps_discard_pop,
                                        UseFishLenComps_discard_pop,
                                        ISS_FishLenComps_discard_pop,
                                        FishAgeComps_discard_pop_LikeType ,
                                        FishLenComps_discard_pop_LikeType,
                                        FishAgeComps_discard_pop_Type,
                                        FishLenComps_discard_pop_Type,
                                        ...
) {

  starting_values <- list(...)
  if(input_list$store_config) input_list$config$Setup_Mod_FishIdx_and_Comps <- c(input_list$config$Setup_Mod_FishIdx_and_Comps, mget(names(formals()))[-1])

  # Input Validation ---------------------------------------------------------
  # Discard Fishery compositions
  check_data_dimensions(ObsFishAgeComps_discard, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsFishAgeComps_discard')
  check_data_dimensions(UseFishAgeComps_discard, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'UseFishAgeComps_discard')
  check_data_dimensions(UseFishLenComps_discard, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'UseFishLenComps_discard')
  if(input_list$data$fit_lengths == 1) check_data_dimensions(ObsFishLenComps_discard, n_regions = input_list$data$n_regions, n_seas = input_list$data$n_seas, n_years = length(input_list$data$years), n_lens = length(input_list$data$lens), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsFishLenComps_discard')
  if(!is.null(ISS_FishAgeComps_discard)) check_data_dimensions(ISS_FishAgeComps_discard, n_regions = input_list$data$n_regions, n_seas = input_list$data$n_seas, n_years = length(input_list$data$years), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ISS_FishAgeComps_discard')
  if(!is.null(ISS_FishLenComps_discard)) check_data_dimensions(ISS_FishLenComps_discard, n_regions = input_list$data$n_regions, n_seas = input_list$data$n_seas, n_years = length(input_list$data$years), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ISS_FishLenComps_discard')
  check_data_dimensions(FishAgeComps_discard_LikeType, n_fish_fleets = input_list$data$n_fish_fleets, what = 'FishAgeComps_discard_LikeType')
  check_data_dimensions(FishLenComps_discard_LikeType, n_fish_fleets = input_list$data$n_fish_fleets, what = 'FishLenComps_discard_LikeType')
  if(!all(FishAgeComps_discard_LikeType %in% c("none", "Multinomial", "Dirichlet-Multinomial", "iid-Logistic-Normal", "1d-Logistic-Normal", "2d-Logistic-Normal")))
    stop("Invalid specification for FishAgeComps_discard_LikeType Should be either none, Multinomial, Dirichlet-Multinomial, iid-Logistic-Normal, 1d-Logistic-Normal, 2d-Logistic-Normal")
  if(!all(FishLenComps_discard_LikeType %in% c("none", "Multinomial", "Dirichlet-Multinomial", "iid-Logistic-Normal", "1d-Logistic-Normal", "2d-Logistic-Normal")))
    stop("Invalid specification for FishLenComps_discard_LikeType Should be either none, Multinomial, Dirichlet-Multinomial, iid-Logistic-Normal, 1d-Logistic-Normal, 2d-Logistic-Normal")

  # Discard Fishery compositions (population-specific)
  if(any(UseFishAgeComps_discard_pop == 1)) check_data_dimensions(ObsFishAgeComps_discard_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsFishAgeComps_discard_pop')
  check_data_dimensions(UseFishAgeComps_discard_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'UseFishAgeComps_discard_pop')
  check_data_dimensions(UseFishLenComps_discard_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'UseFishLenComps_discard_pop')
  if(input_list$data$fit_lengths == 1 && any(UseFishLenComps_discard_pop == 1)) check_data_dimensions(ObsFishLenComps_discard_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_seas = input_list$data$n_seas, n_years = length(input_list$data$years), n_lens = length(input_list$data$lens), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsFishLenComps_discard_pop')
  if(!is.null(ISS_FishAgeComps_discard_pop)) check_data_dimensions(ISS_FishAgeComps_discard_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_seas = input_list$data$n_seas, n_years = length(input_list$data$years), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ISS_FishAgeComps_discard_pop')
  if(!is.null(ISS_FishLenComps_discard_pop)) check_data_dimensions(ISS_FishLenComps_discard_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_seas = input_list$data$n_seas, n_years = length(input_list$data$years), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ISS_FishLenComps_discard_pop')
  check_data_dimensions(FishAgeComps_discard_pop_LikeType, n_fish_fleets = input_list$data$n_fish_fleets, what = 'FishAgeComps_discard_pop_LikeType')
  check_data_dimensions(FishLenComps_discard_pop_LikeType, n_fish_fleets = input_list$data$n_fish_fleets, what = 'FishLenComps_discard_pop_LikeType')
  if(!all(FishAgeComps_discard_pop_LikeType %in% c("none", "Multinomial", "Dirichlet-Multinomial", "iid-Logistic-Normal", "1d-Logistic-Normal", "2d-Logistic-Normal")))
    stop("Invalid specification for FishAgeComps_discard_pop_LikeType Should be either none, Multinomial, Dirichlet-Multinomial, iid-Logistic-Normal, 1d-Logistic-Normal, 2d-Logistic-Normal")
  if(!all(FishLenComps_discard_pop_LikeType %in% c("none", "Multinomial", "Dirichlet-Multinomial", "iid-Logistic-Normal", "1d-Logistic-Normal", "2d-Logistic-Normal")))
    stop("Invalid specification for FishLenComps_discard_pop_LikeType Should be either none, Multinomial, Dirichlet-Multinomial, iid-Logistic-Normal, 1d-Logistic-Normal, 2d-Logistic-Normal")

  # checking to make sure defaults are not applied
  if(any(UseFishAgeComps_discard_pop == 1)) {
    if(is.null(ObsFishAgeComps_discard_pop)) stop("ObsFishAgeComps_discard_pop is NULL, but UseFishAgeComps_discard_pop contains 1s!")
    if(any(str_detect(FishAgeComps_discard_pop_LikeType, "none"))) warning("FishAgeComps_discard_pop_LikeType has nones, but UseFishAgeComps_discard_pop contains 1s! Please verify!")
    if(any(str_detect(FishAgeComps_discard_pop_Type, "none"))) warning("FishAgeComps_discard_pop_Type has nones, but UseFishAgeComps_discard_pop contains 1s! Please verify!")
  }

  if(any(UseFishLenComps_discard_pop == 1)) {
    if(is.null(ObsFishLenComps_discard_pop)) stop("ObsFishLenComps_discard_pop is NULL, but UseFishAgeComps_discard_pop contains 1s!")
    if(any(str_detect(FishLenComps_discard_pop_LikeType, "none"))) warning("FishLenComps_discard_pop_LikeType has nones, but UseFishAgeComps_discard_pop contains 1s! Please verify!")
    if(any(str_detect(FishLenComps_discard_pop_Type, "none"))) warning("FishLenComps_discard_pop_Type has nones, but UseFishAgeComps_discard_pop contains 1s! Please verify!")
  }

  # Discard Fishery Age Composition Options -----------------------------------------

  # Specifying composition likelihood
  comp_fishage_discard_like_vals <- vector()
  for(f in 1:input_list$data$n_fish_fleets) {
    if(FishAgeComps_discard_LikeType[f] == 'none') comp_fishage_discard_like_vals <- c(comp_fishage_discard_like_vals, 999)
    if(FishAgeComps_discard_LikeType[f] == "Multinomial") comp_fishage_discard_like_vals <- c(comp_fishage_discard_like_vals, 0)
    if(FishAgeComps_discard_LikeType[f] == "Dirichlet-Multinomial") comp_fishage_discard_like_vals <- c(comp_fishage_discard_like_vals, 1)
    if(FishAgeComps_discard_LikeType[f] == "iid-Logistic-Normal") comp_fishage_discard_like_vals <- c(comp_fishage_discard_like_vals, 2)
    if(FishAgeComps_discard_LikeType[f] == "1d-Logistic-Normal") comp_fishage_discard_like_vals <- c(comp_fishage_discard_like_vals, 3)
    if(FishAgeComps_discard_LikeType[f] == "2d-Logistic-Normal") comp_fishage_discard_like_vals <- c(comp_fishage_discard_like_vals, 4)
    collect_message(paste("Discard Fishery Age Composition Likelihoods", "for fishery fleet", f, "specified as:" , FishAgeComps_discard_LikeType[f]))
  } # end f loop

  # Specifying composition type
  FishAgeComps_discard_Type_Mat <- array(NA, dim = c(length(input_list$data$years), input_list$data$n_fish_fleets))
  for(i in 1:length(FishAgeComps_discard_Type)) {

    # Extract out components from list
    tmp <- FishAgeComps_discard_Type[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    comps_type_tmp <- tmp_vec[1] # get composition type
    fleet <- as.numeric(tmp_vec[5]) # extract fleet index

    # Checking character string
    if(!comps_type_tmp %in% c("agg", "spltRspltS", "spltRjntS", 'none')) stop("FishAgeComps_discard_Type not specified correctly. Must be one of: agg, spltRspltS, spltRjntS, none")
    if(!fleet %in% c(1:input_list$data$n_fish_fleets)) stop("Invalid fleet specified for FishAgeComps_discard_Type. This needs to be specified as CompType_Year_x-y_Fleet_x")

    # get year ranges
    if(!str_detect(tmp, "terminal")) { # if not terminal year
      year_range <- as.numeric(unlist(strsplit(tmp_vec[3], "-")))
      years <- year_range[1]:year_range[2] # get sequence of years
    } else { # if terminal year
      year_range <- unlist(strsplit(tmp_vec[3], '-'))[1] # get year range
      years <- as.numeric(year_range):length(input_list$data$years) # get sequence of years
    }

    # Composition type
    # define composition types
    if(comps_type_tmp == "agg") {
      if(comp_fishage_discard_like_vals[fleet] == 4) stop("Discard Age composition likelihood specified as 2d-Logistic-Normal, but composition type is aggregated. This is not valid.")
      comps_type_val <- 0
    }
    if(comps_type_tmp == "spltRspltS") comps_type_val <- 1
    if(comps_type_tmp == "spltRjntS") comps_type_val <- 2
    if(comps_type_tmp == "none") comps_type_val <- 999

    # input into matrix
    FishAgeComps_discard_Type_Mat[years,fleet] <- comps_type_val
  } # end i

  if(any(is.na(FishAgeComps_discard_Type_Mat))) stop("FishAgeComps_discard_Type is returning an NA. Did you update the year range of FishAgeComps_discard_Type?")

  # Specifying composition likelihood for population-specific data
  comp_fishage_discard_pop_like_vals <- vector()
  for(f in 1:input_list$data$n_fish_fleets) {
    if(FishAgeComps_discard_pop_LikeType[f] == 'none') comp_fishage_discard_pop_like_vals <- c(comp_fishage_discard_pop_like_vals, 999)
    if(FishAgeComps_discard_pop_LikeType[f] == "Multinomial") comp_fishage_discard_pop_like_vals <- c(comp_fishage_discard_pop_like_vals, 0)
    if(FishAgeComps_discard_pop_LikeType[f] == "Dirichlet-Multinomial") comp_fishage_discard_pop_like_vals <- c(comp_fishage_discard_pop_like_vals, 1)
    if(FishAgeComps_discard_pop_LikeType[f] == "iid-Logistic-Normal") comp_fishage_discard_pop_like_vals <- c(comp_fishage_discard_pop_like_vals, 2)
    if(FishAgeComps_discard_pop_LikeType[f] == "1d-Logistic-Normal") comp_fishage_discard_pop_like_vals <- c(comp_fishage_discard_pop_like_vals, 3)
    if(FishAgeComps_discard_pop_LikeType[f] == "2d-Logistic-Normal") comp_fishage_discard_pop_like_vals <- c(comp_fishage_discard_pop_like_vals, 4)
    collect_message(paste("Discard Population Fishery Age Composition Likelihoods", "for fishery fleet", f, "specified as:" , FishAgeComps_discard_pop_LikeType[f]))
  } # end f loop

  # Specifying composition type
  FishAgeComps_discard_pop_Type_Mat <- array(NA, dim = c(length(input_list$data$years), input_list$data$n_fish_fleets))
  for(i in 1:length(FishAgeComps_discard_pop_Type)) {

    # Extract out components from list
    tmp <- FishAgeComps_discard_pop_Type[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    comps_type_tmp <- tmp_vec[1] # get composition type
    fleet <- as.numeric(tmp_vec[5]) # extract fleet index

    # Checking character string
    if(!comps_type_tmp %in% c("agg", "spltRspltS", "spltRjntS", 'none')) stop("FishAgeComps_discard_pop_Type not specified correctly. Must be one of: agg, spltRspltS, spltRjntS, none")
    if(!fleet %in% c(1:input_list$data$n_fish_fleets)) stop("Invalid fleet specified for FishAgeComps_discard_pop_Type. This needs to be specified as CompType_Year_x-y_Fleet_x")

    # get year ranges
    if(!str_detect(tmp, "terminal")) { # if not terminal year
      year_range <- as.numeric(unlist(strsplit(tmp_vec[3], "-")))
      years <- year_range[1]:year_range[2] # get sequence of years
    } else { # if terminal year
      year_range <- unlist(strsplit(tmp_vec[3], '-'))[1] # get year range
      years <- as.numeric(year_range):length(input_list$data$years) # get sequence of years
    }

    # Composition type
    # define composition types
    if(comps_type_tmp == "agg") {
      if(comp_fishage_discard_pop_like_vals[fleet] == 4) stop("Discard Population Age composition likelihood specified as 2d-Logistic-Normal, but composition type is aggregated. This is not valid.")
      comps_type_val <- 0
    }
    if(comps_type_tmp == "spltRspltS") comps_type_val <- 1
    if(comps_type_tmp == "spltRjntS") comps_type_val <- 2
    if(comps_type_tmp == "none") comps_type_val <- 999

    # input into matrix
    FishAgeComps_discard_pop_Type_Mat[years,fleet] <- comps_type_val
  } # end i

  if(any(is.na(FishAgeComps_discard_pop_Type_Mat))) stop("FishAgeComps_discard_pop_Type is returning an NA. Did you update the year range of FishAgeComps_discard_pop_Type?")

  # Fishery Length Composition Options -----------------------------------------

  # Specifying composition likelihood
  comp_fishlen_discard_like_vals <- vector()
  for(f in 1:input_list$data$n_fish_fleets) {
    if(FishLenComps_discard_LikeType[f] == 'none') comp_fishlen_discard_like_vals <- c(comp_fishlen_discard_like_vals, 999)
    if(FishLenComps_discard_LikeType[f] == "Multinomial") comp_fishlen_discard_like_vals <- c(comp_fishlen_discard_like_vals, 0)
    if(FishLenComps_discard_LikeType[f] == "Dirichlet-Multinomial") comp_fishlen_discard_like_vals <- c(comp_fishlen_discard_like_vals, 1)
    if(FishLenComps_discard_LikeType[f] == "iid-Logistic-Normal") comp_fishlen_discard_like_vals <- c(comp_fishlen_discard_like_vals, 2)
    if(FishLenComps_discard_LikeType[f] == "1d-Logistic-Normal") comp_fishlen_discard_like_vals <- c(comp_fishlen_discard_like_vals, 3)
    if(FishLenComps_discard_LikeType[f] == "2d-Logistic-Normal") comp_fishlen_discard_like_vals <- c(comp_fishlen_discard_like_vals, 4)
    collect_message(paste("Discard Fishery Length Composition Likelihoods", "for fishery fleet", f, "specified as:" , FishLenComps_discard_LikeType[f]))
  } # end f loop

  # Specifying composition type
  FishLenComps_discard_Type_Mat <- array(NA, dim = c(length(input_list$data$years), input_list$data$n_fish_fleets))
  for(i in 1:length(FishLenComps_discard_Type)) {

    # Extract out components from list
    tmp <- FishLenComps_discard_Type[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    comps_type_tmp <- tmp_vec[1] # get composition type
    fleet <- as.numeric(tmp_vec[5]) # extract fleet index

    # define composition types
    if(!comps_type_tmp %in% c("agg", "spltRspltS", "spltRjntS", 'none')) stop("FishLenComps_discard_Type not specified correctly. Must be one of: agg, spltRspltS, spltRjntS, none")
    if(!fleet %in% c(1:input_list$data$n_fish_fleets)) stop("Invalid fleet specified for FishLenComps_discard_Type This needs to be specified as CompType_Year_x-y_Fleet_x")

    # get year ranges
    if(!str_detect(tmp, "terminal")) { # if not terminal year
      year_range <- as.numeric(unlist(strsplit(tmp_vec[3], "-")))
      years <- year_range[1]:year_range[2] # get sequence of years
    } else { # if terminal year
      year_range <- unlist(strsplit(tmp_vec[3], '-'))[1] # get year range
      years <- as.numeric(year_range):length(input_list$data$years) # get sequence of years
    }

    # define composition types
    if(comps_type_tmp == "agg") {
      if(comp_fishlen_discard_like_vals[fleet] == 4) stop("Discard Length composition likelihood specified as 2d-Logistic-Normal, but composition type is aggregated. This is not valid.")
      comps_type_val <- 0
    }
    if(comps_type_tmp == "spltRspltS") comps_type_val <- 1
    if(comps_type_tmp == "spltRjntS") comps_type_val <- 2
    if(comps_type_tmp == "none") comps_type_val <- 999

    # input into matrix
    FishLenComps_discard_Type_Mat[years,fleet] <- comps_type_val
  } # end i

  if(any(is.na(FishLenComps_discard_Type_Mat))) stop("FishLenComps_discard_Type_Mat is returning an NA. Did you update the year range of FishLenComps_discard_Type_Mat?")


  # Specifying composition likelihood for population-specific data
  comp_fishlen_discard_pop_like_vals <- vector()
  for(f in 1:input_list$data$n_fish_fleets) {
    if(FishLenComps_discard_pop_LikeType[f] == 'none') comp_fishlen_discard_pop_like_vals <- c(comp_fishlen_discard_pop_like_vals, 999)
    if(FishLenComps_discard_pop_LikeType[f] == "Multinomial") comp_fishlen_discard_pop_like_vals <- c(comp_fishlen_discard_pop_like_vals, 0)
    if(FishLenComps_discard_pop_LikeType[f] == "Dirichlet-Multinomial") comp_fishlen_discard_pop_like_vals <- c(comp_fishlen_discard_pop_like_vals, 1)
    if(FishLenComps_discard_pop_LikeType[f] == "iid-Logistic-Normal") comp_fishlen_discard_pop_like_vals <- c(comp_fishlen_discard_pop_like_vals, 2)
    if(FishLenComps_discard_pop_LikeType[f] == "1d-Logistic-Normal") comp_fishlen_discard_pop_like_vals <- c(comp_fishlen_discard_pop_like_vals, 3)
    if(FishLenComps_discard_pop_LikeType[f] == "2d-Logistic-Normal") comp_fishlen_discard_pop_like_vals <- c(comp_fishlen_discard_pop_like_vals, 4)
    collect_message(paste("Discard Population Fishery Length Composition Likelihoods", "for fishery fleet", f, "specified as:" , FishLenComps_discard_pop_LikeType[f]))
  } # end f loop

  # Specifying composition type
  FishLenComps_discard_pop_Type_Mat <- array(NA, dim = c(length(input_list$data$years), input_list$data$n_fish_fleets))
  for(i in 1:length(FishLenComps_discard_pop_Type)) {

    # Extract out components from list
    tmp <- FishLenComps_discard_pop_Type[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    comps_type_tmp <- tmp_vec[1] # get composition type
    fleet <- as.numeric(tmp_vec[5]) # extract fleet index

    # Checking character string
    if(!comps_type_tmp %in% c("agg", "spltRspltS", "spltRjntS", 'none')) stop("FishLenComps_discard_pop_Type not specified correctly. Must be one of: agg, spltRspltS, spltRjntS, none")
    if(!fleet %in% c(1:input_list$data$n_fish_fleets)) stop("Invalid fleet specified for FishLenComps_discard_pop_Type. This needs to be specified as CompType_Year_x-y_Fleet_x")

    # get year ranges
    if(!str_detect(tmp, "terminal")) { # if not terminal year
      year_range <- as.numeric(unlist(strsplit(tmp_vec[3], "-")))
      years <- year_range[1]:year_range[2] # get sequence of years
    } else { # if terminal year
      year_range <- unlist(strsplit(tmp_vec[3], '-'))[1] # get year range
      years <- as.numeric(year_range):length(input_list$data$years) # get sequence of years
    }

    # Composition type
    # define composition types
    if(comps_type_tmp == "agg") {
      if(comp_fishlen_discard_pop_like_vals[fleet] == 4) stop("Discard Population Len composition likelihood specified as 2d-Logistic-Normal, but composition type is aggregated. This is not valid.")
      comps_type_val <- 0
    }
    if(comps_type_tmp == "spltRspltS") comps_type_val <- 1
    if(comps_type_tmp == "spltRjntS") comps_type_val <- 2
    if(comps_type_tmp == "none") comps_type_val <- 999

    # input into matrix
    FishLenComps_discard_pop_Type_Mat[years,fleet] <- comps_type_val
  } # end i

  if(any(is.na(FishLenComps_discard_pop_Type_Mat))) stop("FishLenComps_discard_pop_Type is returning an NA. Did you update the year range of FishLenComps_discard_pop_Type?")

  # ISS Munging -------------------------------------------------------------

  # Discard Fishery Ages
  if(is.null(ISS_FishAgeComps_discard)) {
    collect_message("No ISS is specified for FishAgeComps_discard. ISS weighting is calculated by summing up values from ObsFishAgeComps_discard each year")
    ISS_FishAgeComps_discard <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_fish_fleets))
    for(y in 1:length(input_list$data$years)) {
      for(f in 1:input_list$data$n_fish_fleets) {
        for(seas in 1:input_list$data$n_seas) {
          # if aggregated across sexes and regions (0) or joint across sexes
          if(FishAgeComps_discard_Type_Mat[y,f] == 0) ISS_FishAgeComps_discard[1,y,seas,1,f] <- sum(ObsFishAgeComps_discard[,y,seas,,,f])
          # if split by region and sex
          if(FishAgeComps_discard_Type_Mat[y,f] == 1) ISS_FishAgeComps_discard[,y,seas,,f] <- apply(ObsFishAgeComps_discard[,y,seas,,,f, drop = FALSE], c(1,4), sum)
          # if split by region, joint by sex
          if(FishAgeComps_discard_Type_Mat[y,f] == 2) ISS_FishAgeComps_discard[,y,seas,1,f] <- apply(ObsFishAgeComps_discard[,y,seas,,,f, drop = FALSE], 1, sum)
        } # end seas loop
      } # end f loop
    } # end y loop
  }

  # Fishery Lengths
  if(is.null(ISS_FishLenComps_discard)) {
    collect_message("No ISS is specified for FishLenComps_discard. ISS weighting is calculated by summing up values from ObsFishLenComps_discard each year")
    ISS_FishLenComps_discard <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_fish_fleets))
    for(y in 1:length(input_list$data$years)) {
      for(f in 1:input_list$data$n_fish_fleets) {
        for(seas in 1:input_list$data$n_seas) {
          # if aggregated across sexes and regions (0)
          if(FishLenComps_discard_Type_Mat[y,f] == 0) ISS_FishLenComps_discard[1,y,seas,1,f] <- sum(ObsFishLenComps_discard[,y,seas,,,f])
          # if split by region and sex
          if(FishLenComps_discard_Type_Mat[y,f] == 1) ISS_FishLenComps_discard[,y,seas,,f] <- apply(ObsFishLenComps_discard[,y,seas,,,f, drop = FALSE], c(1,4), sum)
          # if split by region, joint by sex
          if(FishLenComps_discard_Type_Mat[y,f] == 2) ISS_FishLenComps_discard[,y,seas,1,f] <- apply(ObsFishLenComps_discard[,y,seas,,,f, drop = FALSE], 1, sum)
        } # end seas loop
      } # end f loop
    } # end y loop
  }

  # Fishery Ages
  if(is.null(ISS_FishAgeComps_discard_pop)) {
    collect_message("No ISS is specified for pop_FishAgeComps_discard. ISS weighting is calculated by summing up values from ObsFishAgeComps_discard_pop each year")
    ISS_FishAgeComps_discard_pop <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_fish_fleets))
    for(p in 1:input_list$data$n_pop) {
      for(y in 1:length(input_list$data$years)) {
        for(f in 1:input_list$data$n_fish_fleets) {
          for(seas in 1:input_list$data$n_seas) {
            # if aggregated across sexes and regions (0) or joint across sexes
            if(FishAgeComps_discard_pop_Type_Mat[y,f] == 0) ISS_FishAgeComps_discard_pop[p,1,y,seas,1,f] <- sum(ObsFishAgeComps_discard_pop[p,,y,seas,,,f])
            # if split by region and sex
            if(FishAgeComps_discard_pop_Type_Mat[y,f] == 1) ISS_FishAgeComps_discard_pop[p,,y,seas,,f] <- apply(ObsFishAgeComps_discard_pop[p,,y,seas,,,f, drop = FALSE], c(2,5), sum)
            # if split by region, joint by sex
            if(FishAgeComps_discard_pop_Type_Mat[y,f] == 2) ISS_FishAgeComps_discard_pop[p,,y,seas,1,f] <- apply(ObsFishAgeComps_discard_pop[p,,y,seas,,,f, drop = FALSE], 2, sum)
          } # end seas loop
        } # end f loop
      } # end y loop
    } # end p loop
  }

  # Fishery Lengths
  if(is.null(ISS_FishLenComps_discard_pop)) {
    collect_message("No ISS is specified for pop_FishLenComps_discard. ISS weighting is calculated by summing up values from ObsFishLenComps_discard_pop each year")
    ISS_FishLenComps_discard_pop <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_fish_fleets))
    for(p in 1:input_list$data$n_pop) {
      for(y in 1:length(input_list$data$years)) {
        for(f in 1:input_list$data$n_fish_fleets) {
          for(seas in 1:input_list$data$n_seas) {
            # if aggregated across sexes and regions (0)
            if(FishLenComps_discard_pop_Type_Mat[y,f] == 0) ISS_FishLenComps_discard_pop[p,1,y,seas,1,f] <- sum(ObsFishLenComps_discard_pop[p,,y,seas,,,f])
            # if split by region and sex
            if(FishLenComps_discard_pop_Type_Mat[y,f] == 1) ISS_FishLenComps_discard_pop[p,,y,seas,,f] <- apply(ObsFishLenComps_discard_pop[p,,y,seas,,,f, drop = FALSE], c(2,5), sum)
            # if split by region, joint by sex
            if(FishLenComps_discard_pop_Type_Mat[y,f] == 2) ISS_FishLenComps_discard_pop[p,,y,seas,1,f] <- apply(ObsFishLenComps_discard_pop[p,,y,seas,,,f, drop = FALSE], 2, sum)
          } # end seas loop
        } # end f loop
      } # end y loop
    } # end p loop
  }

  # Populate Data List ------------------------------------------------------

  input_list$data$ISS_FishAgeComps_discard <- ISS_FishAgeComps_discard
  input_list$data$ISS_FishLenComps_discard <- ISS_FishLenComps_discard
  input_list$data$ISS_FishAgeComps_discard_pop <- ISS_FishAgeComps_discard_pop
  input_list$data$ISS_FishLenComps_discard_pop <- ISS_FishLenComps_discard_pop
  input_list$data$ObsFishAgeComps_discard <- ObsFishAgeComps_discard
  input_list$data$UseFishAgeComps_discard <- UseFishAgeComps_discard
  input_list$data$ObsFishLenComps_discard <- ObsFishLenComps_discard
  input_list$data$UseFishLenComps_discard <- UseFishLenComps_discard
  input_list$data$ObsFishAgeComps_discard_pop <- ObsFishAgeComps_discard_pop
  input_list$data$UseFishAgeComps_discard_pop <- UseFishAgeComps_discard_pop
  input_list$data$ObsFishLenComps_discard_pop <- ObsFishLenComps_discard_pop
  input_list$data$UseFishLenComps_discard_pop <- UseFishLenComps_discard_pop
  input_list$data$FishAgeComps_discard_LikeType <- comp_fishage_discard_like_vals
  input_list$data$FishLenComps_discard_LikeType <- comp_fishlen_discard_like_vals
  input_list$data$FishAgeComps_discard_pop_LikeType <- comp_fishage_discard_pop_like_vals
  input_list$data$FishLenComps_discard_pop_LikeType <- comp_fishlen_discard_pop_like_vals
  input_list$data$FishAgeComps_discard_Type <- FishAgeComps_discard_Type_Mat
  input_list$data$FishLenComps_discard_Type <- FishLenComps_discard_Type_Mat
  input_list$data$FishAgeComps_discard_pop_Type <- FishAgeComps_discard_pop_Type_Mat
  input_list$data$FishLenComps_discard_pop_Type <- FishLenComps_discard_pop_Type_Mat

  # Populate Parameter List -------------------------------------------------

  # Dispersion parameters for the fishery age comps
  if("ln_FishAge_discard_theta" %in% names(starting_values)) input_list$par$ln_FishAge_discard_theta <- starting_values$ln_FishAge_discard_theta
  else input_list$par$ln_FishAge_discard_theta <- array(0, dim = c(input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_fish_fleets))

  # logistic normal correlation parameters for fishery age comps
  if("FishAge_discard_corr_pars" %in% names(starting_values)) input_list$par$FishAge_discard_corr_pars <- starting_values$FishAge_discard_corr_pars
  else input_list$par$FishAge_discard_corr_pars <- array(0.01, dim = c(input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_fish_fleets, 2))

  # aggregated
  if("ln_FishAge_discard_theta_agg" %in% names(starting_values)) input_list$par$ln_FishAge_discard_theta_agg <- starting_values$ln_FishAge_discard_theta_agg
  else input_list$par$ln_FishAge_discard_theta_agg <- array(0, dim = c(input_list$data$n_fish_fleets))

  # aggregated correlation parameters
  if("FishAge_discard_corr_pars_agg" %in% names(starting_values)) input_list$par$FishAge_discard_corr_pars_agg <- starting_values$FishAge_discard_corr_pars_agg
  else input_list$par$FishAge_discard_corr_pars_agg <- array(0.01, dim = c(input_list$data$n_fish_fleets))

  # Dispersion parameters for fishery length comps
  if("ln_FishLen_discard_theta" %in% names(starting_values)) input_list$par$ln_FishLen_discard_theta <- starting_values$ln_FishLen_discard_theta
  else input_list$par$ln_FishLen_discard_theta <- array(0, dim = c(input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_fish_fleets))

  # logistic normal correlation parameters for fishery length comps
  if("FishLen_discard_corr_pars" %in% names(starting_values)) input_list$par$FishLen_discard_corr_pars <- starting_values$FishLen_discard_corr_pars
  else input_list$par$FishLen_discard_corr_pars <- array(0.01, dim = c(input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_fish_fleets, 2))

  # aggregated
  if("ln_FishLen_discard_theta_agg" %in% names(starting_values)) input_list$par$ln_FishLen_discard_theta_agg <- starting_values$ln_FishLen_discard_theta_agg
  else input_list$par$ln_FishLen_discard_theta_agg <- array(0, dim = c(input_list$data$n_fish_fleets))

  if("FishLen_discard_corr_pars_agg" %in% names(starting_values)) input_list$par$FishLen_discard_corr_pars_agg <- starting_values$FishLen_discard_corr_pars_agg
  else input_list$par$FishLen_discard_corr_pars_agg <- array(0.01, dim = c(input_list$data$n_fish_fleets))

  # Dispersion parameters for the population fishery age comps
  if("ln_FishAge_discard_pop_theta" %in% names(starting_values)) input_list$par$ln_FishAge_discard_pop_theta <- starting_values$ln_FishAge_discard_pop_theta
  else input_list$par$ln_FishAge_discard_pop_theta <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_fish_fleets))

  # logistic normal correlation parameters for population fishery age comps
  if("FishAge_discard_pop_corr_pars" %in% names(starting_values)) input_list$par$FishAge_discard_pop_corr_pars <- starting_values$FishAge_discard_pop_corr_pars
  else input_list$par$FishAge_discard_pop_corr_pars <- array(0.01, dim = c(input_list$data$n_pop,input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_fish_fleets, 2))

  # aggregated population pars
  if("ln_FishAge_discard_pop_theta_agg" %in% names(starting_values)) input_list$par$ln_FishAge_discard_pop_theta_agg <- starting_values$ln_FishAge_discard_pop_theta_agg
  else input_list$par$ln_FishAge_discard_pop_theta_agg <- array(0, dim = c(input_list$data$n_pop,input_list$data$n_fish_fleets))

  # aggregated population correlation parameters
  if("FishAge_discard_pop_corr_pars_agg" %in% names(starting_values)) input_list$par$FishAge_discard_pop_corr_pars_agg <- starting_values$FishAge_discard_pop_corr_pars_agg
  else input_list$par$FishAge_discard_pop_corr_pars_agg <- array(0.01, dim = c(input_list$data$n_pop,input_list$data$n_fish_fleets))

  # Dispersion parameters for population fishery length comps
  if("ln_FishLen_discard_pop_theta" %in% names(starting_values)) input_list$par$ln_FishLen_discard_pop_theta <- starting_values$ln_FishLen_discard_pop_theta
  else input_list$par$ln_FishLen_discard_pop_theta <- array(0, dim = c(input_list$data$n_pop,input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_fish_fleets))

  # logistic normal correlation parameters for population fishery length comps
  if("FishLen_discard_pop_corr_pars" %in% names(starting_values)) input_list$par$FishLen_discard_pop_corr_pars <- starting_values$FishLen_discard_pop_corr_pars
  else input_list$par$FishLen_discard_pop_corr_pars <- array(0.01, dim = c(input_list$data$n_pop,input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_fish_fleets, 2))

  # aggregated population pars
  if("ln_FishLen_discard_pop_theta_agg" %in% names(starting_values)) input_list$par$ln_FishLen_discard_pop_theta_agg <- starting_values$ln_FishLen_discard_pop_theta_agg
  else input_list$par$ln_FishLen_discard_pop_theta_agg <- array(0, dim = c(input_list$data$n_pop,input_list$data$n_fish_fleets))

  if("FishLen_discard_pop_corr_pars_agg" %in% names(starting_values)) input_list$par$FishLen_discard_pop_corr_pars_agg <- starting_values$FishLen_discard_pop_corr_pars_agg
  else input_list$par$FishLen_discard_pop_corr_pars_agg <- array(0.01, dim = c(input_list$data$n_pop, input_list$data$n_fish_fleets))

  # Mapping Options ---------------------------------------------------------

  input_list <- do_FishAge_discard_theta_mapping(input_list)
  input_list <- do_FishLen_discard_theta_mapping(input_list)
  input_list <- do_FishAge_discard_corr_pars_mapping(input_list)
  input_list <- do_FishLen_discard_corr_pars_mapping(input_list)

  input_list <- do_FishAge_discard_pop_theta_mapping(input_list)
  input_list <- do_FishLen_discard_pop_theta_mapping(input_list)
  input_list <- do_FishAge_discard_pop_corr_pars_mapping(input_list)
  input_list <- do_FishLen_discard_pop_corr_pars_mapping(input_list)

  # Print Messages ----------------------------------------------------------
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)
}

#' Set up discards, fishery index, age composition, and length composition inputs
#'
#' Populates \code{input_list} with observed fishery indices, age compositions,
#' and length compositions (both pooled and population-specific) along with
#' their usage indicators, likelihood types, composition structure types, input
#' sample sizes, and overdispersion and correlation parameter starting values
#' and mappings. Must be called after \code{\link{Setup_Mod_Catch_and_F}}.
#'
#' When \code{ISS_FishAgeComps}, \code{ISS_FishLenComps},
#' \code{ISS_FishAgeComps_pop}, or \code{ISS_FishLenComps_pop} are \code{NULL},
#' input sample sizes are derived automatically by summing the observed
#' composition arrays within each year–fleet–season–region cell, consistent
#' with the specified composition type.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, \code{$map},
#'   and \code{$verbose} sublists, as returned by upstream setup functions.
#' @param ObsFishIdx Observed fishery CPUE or biomass index array
#'   \code{[n_regions × n_years × n_seas × n_fish_fleets]}.
#' @param ObsFishIdx_SE Standard errors of \code{ObsFishIdx} on the log scale,
#'   same dimensions as \code{ObsFishIdx}.
#' @param UseFishIdx Binary indicator array \code{[n_regions × n_years × n_seas × n_fish_fleets]}.
#'   \code{1} = include index in the likelihood; \code{0} = exclude.
#' @param fish_idx_type Character vector of length \code{n_fish_fleets} specifying
#'   the index type for each fleet. \code{"biom"} = biomass; \code{"abd"} =
#'   abundance; \code{"none"} = no index for this fleet.
#' @param ObsFishIdx_pop Observed population-specific fishery index array
#'   \code{[n_pop × n_regions × n_years × n_seas × n_fish_fleets]}.
#' @param ObsFishIdx_pop_SE Lognormal standard errors for \code{ObsFishIdx_pop},
#'   same dimensions \code{[n_pop × n_regions × n_years × n_seas × n_fish_fleets]}.
#' @param UseFishIdx_pop Binary indicator array
#'   \code{[n_pop × n_regions × n_years × n_seas × n_fish_fleets]}. \code{1} =
#'   include population-specific index in likelihood; \code{0} = exclude.
#'   Default: all zeros.
#' @param ObsFishAgeComps Observed fishery age composition array
#'   \code{[n_regions × n_years × n_seas × n_ages × n_sexes × n_fish_fleets]}.
#'   Values may be raw counts or proportions; if proportions, supply
#'   \code{ISS_FishAgeComps} explicitly.
#' @param UseFishAgeComps Binary indicator array
#'   \code{[n_regions × n_years × n_seas × n_fish_fleets]}.
#'   \code{1} = fit age compositions; \code{0} = exclude.
#' @param ISS_FishAgeComps Input sample size array
#'   \code{[n_regions × n_years × n_seas × n_sexes × n_fish_fleets]}.
#'   If \code{NULL} (default), computed automatically by summing
#'   \code{ObsFishAgeComps} within each year–fleet–season–region cell
#'   according to \code{FishAgeComps_Type}.
#' @param ObsFishLenComps Observed fishery length composition array
#'   \code{[n_regions × n_years × n_seas × n_lens × n_sexes × n_fish_fleets]}.
#'   Only required when \code{input_list$data$fit_lengths == 1}.
#' @param UseFishLenComps Binary indicator array
#'   \code{[n_regions × n_years × n_seas × n_fish_fleets]}.
#'   \code{1} = fit length compositions; \code{0} = exclude.
#' @param ISS_FishLenComps Input sample size array for length compositions
#'   \code{[n_regions × n_years × n_seas × n_sexes × n_fish_fleets]}.
#'   If \code{NULL} (default), derived automatically from \code{ObsFishLenComps}.
#' @param FishAgeComps_LikeType Character vector of length \code{n_fish_fleets}
#'   specifying the likelihood for fishery age compositions. Options:
#'   \code{"Multinomial"}, \code{"Dirichlet-Multinomial"},
#'   \code{"iid-Logistic-Normal"}, \code{"1d-Logistic-Normal"},
#'   \code{"2d-Logistic-Normal"}, \code{"none"}.
#' @param FishLenComps_LikeType Same as \code{FishAgeComps_LikeType} but for
#'   length compositions.
#' @param FishAgeComps_Type Character vector defining the age composition
#'   structure (aggregation level) for each fleet and time period. Each element
#'   must follow the format \code{"<type>_Year_<start>-<end>_Fleet_<f>"} or
#'   \code{"<type>_Year_<start>-terminal_Fleet_<f>"}. Valid types:
#'   \describe{
#'     \item{\code{"agg"}}{Aggregated across regions and sexes
#'       (incompatible with \code{"2d-Logistic-Normal"}).}
#'     \item{\code{"spltRspltS"}}{Split by region and sex.}
#'     \item{\code{"spltRjntS"}}{Split by region, summed jointly across sexes.}
#'     \item{\code{"none"}}{No composition data for this fleet and period.}
#'   }
#'   Example: \code{c("spltRjntS_Year_1-10_Fleet_1", "agg_Year_11-terminal_Fleet_1")}.
#' @param FishLenComps_Type Same format and options as \code{FishAgeComps_Type}
#'   but applied to length compositions.
#' @param ObsFishAgeComps_pop Observed population-specific fishery age
#'   composition array
#'   \code{[n_pop × n_regions × n_years × n_seas × n_ages × n_sexes × n_fish_fleets]}.
#'   Required when any element of \code{UseFishAgeComps_pop} is \code{1}.
#' @param UseFishAgeComps_pop Binary indicator array
#'   \code{[n_pop × n_regions × n_years × n_seas × n_fish_fleets]}.
#'   \code{1} = fit population-specific age compositions; \code{0} = exclude.
#'   Default: all zeros.
#' @param ISS_FishAgeComps_pop Input sample size array for population-specific
#'   age compositions
#'   \code{[n_pop × n_regions × n_years × n_seas × n_sexes × n_fish_fleets]}.
#'   If \code{NULL} (default), computed automatically by summing
#'   \code{ObsFishAgeComps_pop} within each population–year–fleet–season–region
#'   cell according to \code{FishAgeComps_pop_Type}.
#' @param ObsFishLenComps_pop Observed population-specific fishery length
#'   composition array
#'   \code{[n_pop × n_regions × n_years × n_seas × n_lens × n_sexes × n_fish_fleets]}.
#'   Required when \code{input_list$data$fit_lengths == 1} and any element of
#'   \code{UseFishLenComps_pop} is \code{1}.
#' @param UseFishLenComps_pop Binary indicator array
#'   \code{[n_pop × n_regions × n_years × n_seas × n_fish_fleets]}.
#'   \code{1} = fit population-specific length compositions; \code{0} = exclude.
#'   Default: all zeros.
#' @param ISS_FishLenComps_pop Input sample size array for population-specific
#'   length compositions
#'   \code{[n_pop × n_regions × n_years × n_seas × n_sexes × n_fish_fleets]}.
#'   If \code{NULL} (default), derived automatically from
#'   \code{ObsFishLenComps_pop}.
#' @param FishAgeComps_pop_LikeType Character vector of length
#'   \code{n_fish_fleets} specifying the likelihood for population-specific
#'   fishery age compositions. Same options as \code{FishAgeComps_LikeType}.
#'   Default: \code{"none"} for all fleets.
#' @param FishLenComps_pop_LikeType Character vector of length
#'   \code{n_fish_fleets} specifying the likelihood for population-specific
#'   fishery length compositions. Same options as \code{FishLenComps_LikeType}.
#'   Default: \code{"none"} for all fleets.
#' @param FishAgeComps_pop_Type Character vector defining the composition
#'   structure for population-specific age compositions. Same format and options
#'   as \code{FishAgeComps_Type}. Default: \code{"none"} for all fleets across
#'   all years.
#' @param FishLenComps_pop_Type Character vector defining the composition
#'   structure for population-specific length compositions. Same format and
#'   options as \code{FishLenComps_Type}. Default: \code{"none"} for all fleets
#'   across all years.
#' @param ... Optional starting value overrides for overdispersion and
#'   correlation parameters.
#' @param ObsFishAgeComps_discard Observed fishery age composition from discards
#'   \code{[n_regions × n_years × n_seas × n_ages × n_sexes × n_fish_fleets]}.
#'   Structure must match \code{ObsFishAgeComps}.
#'
#' @param UseFishAgeComps_discard Binary indicator array for discard age compositions
#'   \code{[n_regions × n_years × n_seas × n_fish_fleets]}.
#'   \code{1} = include discard age compositions in likelihood; \code{0} = exclude.
#'
#' @param ISS_FishAgeComps_discard Input sample size array for discard age compositions
#'   \code{[n_regions × n_years × n_seas × n_sexes × n_fish_fleets]}.
#'   If \code{NULL}, derived automatically from \code{ObsFishAgeComps_discard}
#'   using \code{FishAgeComps_discard_Type}.
#'
#' @param ObsFishLenComps_discard Observed fishery length composition from discards
#'   \code{[n_regions × n_years × n_seas × n_lens × n_sexes × n_fish_fleets]}.
#'   Required if \code{input_list$data$fit_lengths == 1}.
#'
#' @param UseFishLenComps_discard Binary indicator array for discard length compositions
#'   \code{[n_regions × n_years × n_seas × n_fish_fleets]}.
#'   \code{1} = include discard length compositions in likelihood; \code{0} = exclude.
#'
#' @param ISS_FishLenComps_discard Input sample size array for discard length compositions
#'   \code{[n_regions × n_years × n_seas × n_sexes × n_fish_fleets]}.
#'   If \code{NULL}, derived automatically from \code{ObsFishLenComps_discard}.
#'
#' @param FishAgeComps_discard_LikeType Character vector of length
#'   \code{n_fish_fleets} specifying likelihood type for discard age compositions.
#'   Options:
#'   \describe{
#'     \item{\code{"Multinomial"}}{Standard multinomial likelihood}
#'     \item{\code{"Dirichlet-Multinomial"}}{Overdispersed multinomial}
#'     \item{\code{"iid-Logistic-Normal"}}{Independent logistic-normal}
#'     \item{\code{"1d-Logistic-Normal"}}{1D correlated logistic-normal}
#'     \item{\code{"2d-Logistic-Normal"}}{2D correlated logistic-normal}
#'     \item{\code{"none"}}{No discard age composition likelihood}
#'   }
#'
#' @param FishLenComps_discard_LikeType Same specification as
#'   \code{FishAgeComps_discard_LikeType}, but for discard length compositions.
#'
#' @param FishAgeComps_discard_Type Character vector defining discard age composition
#'   structure by fleet and year block.
#'   Format:
#'   \code{"<type>_Year_<start>-<end>_Fleet_<f>"} or
#'   \code{"<type>_Year_<start>-terminal_Fleet_<f>"}.
#'   Valid types:
#'   \describe{
#'     \item{\code{"agg"}}{Aggregated across regions and sexes}
#'     \item{\code{"spltRspltS"}}{Split by region and sex}
#'     \item{\code{"spltRjntS"}}{Split by region, joint across sexes}
#'     \item{\code{"none"}}{No discard age composition}
#'   }
#'
#' @param FishLenComps_discard_Type Same format and options as
#'   \code{FishAgeComps_discard_Type}, applied to discard length compositions.
#'
#' @param ObsFishAgeComps_discard_pop Observed population-specific discard age
#'   composition array
#'   \code{[n_pop × n_regions × n_years × n_seas × n_ages × n_sexes × n_fish_fleets]}.
#'
#' @param UseFishAgeComps_discard_pop Binary indicator array for population-specific
#'   discard age compositions
#'   \code{[n_pop × n_regions × n_years × n_seas × n_fish_fleets]}.
#'
#' @param ISS_FishAgeComps_discard_pop Input sample size array for population-specific
#'   discard age compositions
#'   \code{[n_pop × n_regions × n_years × n_seas × n_sexes × n_fish_fleets]}.
#'   If \code{NULL}, computed from \code{ObsFishAgeComps_discard_pop}.
#'
#' @param ObsFishLenComps_discard_pop Observed population-specific discard length
#'   composition array
#'   \code{[n_pop × n_regions × n_years × n_seas × n_lens × n_sexes × n_fish_fleets]}.
#'
#' @param UseFishLenComps_discard_pop Binary indicator array for population-specific
#'   discard length compositions
#'   \code{[n_pop × n_regions × n_years × n_seas × n_fish_fleets]}.
#'
#' @param ISS_FishLenComps_discard_pop Input sample size array for population-specific
#'   discard length compositions
#'   \code{[n_pop × n_regions × n_years × n_seas × n_sexes × n_fish_fleets]}.
#'   If \code{NULL}, derived from \code{ObsFishLenComps_discard_pop}.
#'
#' @param FishAgeComps_discard_pop_LikeType Character vector of length
#'   \code{n_fish_fleets} specifying likelihood type for population-specific
#'   discard age compositions. Same options as \code{FishAgeComps_discard_LikeType}.
#'
#' @param FishLenComps_discard_pop_LikeType Same as above but for discard length compositions.
#'
#' @param FishAgeComps_discard_pop_Type Character vector defining structure for
#'   population-specific discard age compositions. Same format as
#'   \code{FishAgeComps_discard_Type}.
#'
#' @param FishLenComps_discard_pop_Type Character vector defining structure for
#'   population-specific discard length compositions. Same format as
#'   \code{FishLenComps_discard_Type}.
#'
#' @return The input \code{input_list} with \code{$data}, \code{$par}, and
#'   \code{$map} updated with all fishery index and composition fields, including
#'   pooled and population-specific observed arrays, computed or supplied ISS
#'   arrays, integer-coded likelihood and composition type matrices,
#'   overdispersion parameters, and their factor maps.
#'
#' @export Setup_Mod_FishIdx_and_Comps
#' @importFrom stringr str_detect
#' @family Model Setup
Setup_Mod_FishIdx_and_Comps <- function(input_list,
                                        ObsFishIdx,
                                        ObsFishIdx_SE,
                                        ObsFishIdx_pop = NULL,
                                        ObsFishIdx_pop_SE = NULL,
                                        UseFishIdx_pop = array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets)),
                                        fish_idx_type,
                                        UseFishIdx,

                                        # Retained Compositions
                                        ObsFishAgeComps,
                                        UseFishAgeComps,
                                        ISS_FishAgeComps = NULL,
                                        ObsFishLenComps,
                                        UseFishLenComps,
                                        ISS_FishLenComps = NULL,
                                        FishAgeComps_LikeType,
                                        FishLenComps_LikeType,
                                        FishAgeComps_Type,
                                        FishLenComps_Type,
                                        ObsFishAgeComps_pop = NULL,
                                        UseFishAgeComps_pop = array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets)),
                                        ISS_FishAgeComps_pop = NULL,
                                        ObsFishLenComps_pop = NULL,
                                        UseFishLenComps_pop = array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets)),
                                        ISS_FishLenComps_pop = NULL,
                                        FishAgeComps_pop_LikeType = rep("none", input_list$data$n_fish_fleets),
                                        FishLenComps_pop_LikeType = rep("none", input_list$data$n_fish_fleets),
                                        FishAgeComps_pop_Type = paste("none_Year_1-terminal_Fleet_", 1:input_list$data$n_fish_fleets, sep = ''),
                                        FishLenComps_pop_Type = paste("none_Year_1-terminal_Fleet_", 1:input_list$data$n_fish_fleets, sep = ''),

                                        # Discard Compositions (forwarded to Setup_Mod_Discard_Comps)
                                        ObsFishAgeComps_discard = array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, length(input_list$data$ages), input_list$data$n_sexes, input_list$data$n_fish_fleets)),
                                        UseFishAgeComps_discard = array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets)),
                                        ISS_FishAgeComps_discard = NULL,
                                        ObsFishLenComps_discard = array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, length(input_list$data$lens), input_list$data$n_sexes, input_list$data$n_fish_fleets)),
                                        UseFishLenComps_discard = array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets)),
                                        ISS_FishLenComps_discard = NULL,
                                        FishAgeComps_discard_LikeType = rep("none", input_list$data$n_fish_fleets),
                                        FishLenComps_discard_LikeType = rep("none", input_list$data$n_fish_fleets),
                                        FishAgeComps_discard_Type = paste("none_Year_1-terminal_Fleet_", 1:input_list$data$n_fish_fleets, sep = ''),
                                        FishLenComps_discard_Type = paste("none_Year_1-terminal_Fleet_", 1:input_list$data$n_fish_fleets, sep = ''),
                                        ObsFishAgeComps_discard_pop = NULL,
                                        UseFishAgeComps_discard_pop = array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets)),
                                        ISS_FishAgeComps_discard_pop = NULL,
                                        ObsFishLenComps_discard_pop = NULL,
                                        UseFishLenComps_discard_pop = array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets)),
                                        ISS_FishLenComps_discard_pop = NULL,
                                        FishAgeComps_discard_pop_LikeType = rep("none", input_list$data$n_fish_fleets),
                                        FishLenComps_discard_pop_LikeType = rep("none", input_list$data$n_fish_fleets),
                                        FishAgeComps_discard_pop_Type = paste("none_Year_1-terminal_Fleet_", 1:input_list$data$n_fish_fleets, sep = ''),
                                        FishLenComps_discard_pop_Type = paste("none_Year_1-terminal_Fleet_", 1:input_list$data$n_fish_fleets, sep = ''),
                                        ...
                                        ) {

  messages_list <<- character(0) # string to attach to for printing messages
  starting_values <- list(...)
  if(input_list$store_config) input_list$config$Setup_Mod_FishIdx_and_Comps <- mget(names(formals()))[-1]

  # Input Validation ---------------------------------------------------------

  # Fishery indices
  check_data_dimensions(ObsFishIdx, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsFishIdx')
  check_data_dimensions(ObsFishIdx_SE, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsFishIdx_SE')
  check_data_dimensions(UseFishIdx, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'UseFishIdx')

  if(any(UseFishIdx_pop == 1)) {
    check_data_dimensions(ObsFishIdx_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsFishIdx_pop')
    check_data_dimensions(ObsFishIdx_pop_SE, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsFishIdx_pop_SE')
    check_data_dimensions(UseFishIdx_pop, n_pop = input_list$data$n_pop,  n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'UseFishIdx_pop')
  }

  if(!all(fish_idx_type %in% c("biom", "abd", "none"))) stop("Invalid specification for fish_idx_type. Should be either abd, biom, or none")

  # Fishery compositions
  check_data_dimensions(ObsFishAgeComps, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsFishAgeComps')
  check_data_dimensions(UseFishAgeComps, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'UseFishAgeComps')
  check_data_dimensions(UseFishLenComps, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'UseFishLenComps')
  if(input_list$data$fit_lengths == 1) check_data_dimensions(ObsFishLenComps, n_regions = input_list$data$n_regions, n_seas = input_list$data$n_seas, n_years = length(input_list$data$years), n_lens = length(input_list$data$lens), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsFishLenComps')
  if(!is.null(ISS_FishAgeComps)) check_data_dimensions(ISS_FishAgeComps, n_regions = input_list$data$n_regions, n_seas = input_list$data$n_seas, n_years = length(input_list$data$years), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ISS_FishAgeComps')
  if(!is.null(ISS_FishLenComps)) check_data_dimensions(ISS_FishLenComps, n_regions = input_list$data$n_regions, n_seas = input_list$data$n_seas, n_years = length(input_list$data$years), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ISS_FishLenComps')
  check_data_dimensions(FishAgeComps_LikeType, n_fish_fleets = input_list$data$n_fish_fleets, what = 'FishAgeComps_LikeType')
  check_data_dimensions(FishLenComps_LikeType, n_fish_fleets = input_list$data$n_fish_fleets, what = 'FishLenComps_LikeType')
  if(!all(FishAgeComps_LikeType %in% c("none", "Multinomial", "Dirichlet-Multinomial", "iid-Logistic-Normal", "1d-Logistic-Normal", "2d-Logistic-Normal")))
    stop("Invalid specification for FishAgeComps_LikeType Should be either none, Multinomial, Dirichlet-Multinomial, iid-Logistic-Normal, 1d-Logistic-Normal, 2d-Logistic-Normal")
  if(!all(FishLenComps_LikeType %in% c("none", "Multinomial", "Dirichlet-Multinomial", "iid-Logistic-Normal", "1d-Logistic-Normal", "2d-Logistic-Normal")))
    stop("Invalid specification for FishLenComps_LikeType Should be either none, Multinomial, Dirichlet-Multinomial, iid-Logistic-Normal, 1d-Logistic-Normal, 2d-Logistic-Normal")

   # Fishery compositions (population-specific)
  if(any(UseFishAgeComps_pop == 1)) check_data_dimensions(ObsFishAgeComps_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsFishAgeComps_pop')
  check_data_dimensions(UseFishAgeComps_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'UseFishAgeComps_pop')
  check_data_dimensions(UseFishLenComps_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'UseFishLenComps_pop')
  if(input_list$data$fit_lengths == 1 && any(UseFishLenComps_pop == 1)) check_data_dimensions(ObsFishLenComps_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_seas = input_list$data$n_seas, n_years = length(input_list$data$years), n_lens = length(input_list$data$lens), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsFishLenComps_pop')
  if(!is.null(ISS_FishAgeComps_pop)) check_data_dimensions(ISS_FishAgeComps_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_seas = input_list$data$n_seas, n_years = length(input_list$data$years), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ISS_FishAgeComps_pop')
  if(!is.null(ISS_FishLenComps_pop)) check_data_dimensions(ISS_FishLenComps_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_seas = input_list$data$n_seas, n_years = length(input_list$data$years), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ISS_FishLenComps_pop')
  check_data_dimensions(FishAgeComps_pop_LikeType, n_fish_fleets = input_list$data$n_fish_fleets, what = 'FishAgeComps_pop_LikeType')
  check_data_dimensions(FishLenComps_pop_LikeType, n_fish_fleets = input_list$data$n_fish_fleets, what = 'FishLenComps_pop_LikeType')
  if(!all(FishAgeComps_pop_LikeType %in% c("none", "Multinomial", "Dirichlet-Multinomial", "iid-Logistic-Normal", "1d-Logistic-Normal", "2d-Logistic-Normal")))
    stop("Invalid specification for FishAgeComps_pop_LikeType Should be either none, Multinomial, Dirichlet-Multinomial, iid-Logistic-Normal, 1d-Logistic-Normal, 2d-Logistic-Normal")
  if(!all(FishLenComps_pop_LikeType %in% c("none", "Multinomial", "Dirichlet-Multinomial", "iid-Logistic-Normal", "1d-Logistic-Normal", "2d-Logistic-Normal")))
    stop("Invalid specification for FishLenComps_pop_LikeType Should be either none, Multinomial, Dirichlet-Multinomial, iid-Logistic-Normal, 1d-Logistic-Normal, 2d-Logistic-Normal")

  # checking to make sure defaults are not applied
  if(any(UseFishAgeComps_pop == 1)) {
    if(is.null(ObsFishAgeComps_pop)) stop("ObsFishAgeComps_pop is NULL, but UseFishAgeComps_pop contains 1s!")
    if(any(str_detect(FishAgeComps_pop_LikeType, "none"))) warning("FishAgeComps_pop_LikeType has nones, but UseFishAgeComps_pop contains 1s! Please verify!")
    if(any(str_detect(FishAgeComps_pop_Type, "none"))) warning("FishAgeComps_pop_Type has nones, but UseFishAgeComps_pop contains 1s! Please verify!")
  }

  if(any(UseFishLenComps_pop == 1)) {
    if(is.null(ObsFishLenComps_pop)) stop("ObsFishLenComps_pop is NULL, but UseFishLenComps_pop contains 1s!")
    if(any(str_detect(FishLenComps_pop_LikeType, "none"))) warning("FishLenComps_pop_LikeType has nones, but UseFishAgeComps_pop contains 1s! Please verify!")
    if(any(str_detect(FishLenComps_pop_Type, "none"))) warning("FishLenComps_pop_Type has nones, but UseFishAgeComps_pop contains 1s! Please verify!")
  }


  # Fishery Index Options ---------------------------------------------------

  fish_idx_type_vals <- array(NA, dim = c(input_list$data$n_fish_fleets))
  for(f in 1:input_list$data$n_fish_fleets) {
    if(fish_idx_type[f] == 'biom') fish_idx_type_vals[f] <- 1 # biomass
    if(fish_idx_type[f] == 'abd') fish_idx_type_vals[f] <- 0 # abundance
    if(fish_idx_type[f] == 'none') fish_idx_type_vals[f] <- 999 # none
    collect_message(paste("Fishery Index", "for fishery fleet", f, "specified as:" , fish_idx_type[f]))
  } # end f loop


  # Fishery Age Composition Options -----------------------------------------

  # Specifying composition likelihood
  comp_fishage_like_vals <- vector()
  for(f in 1:input_list$data$n_fish_fleets) {
    if(FishAgeComps_LikeType[f] == 'none') comp_fishage_like_vals <- c(comp_fishage_like_vals, 999)
    if(FishAgeComps_LikeType[f] == "Multinomial") comp_fishage_like_vals <- c(comp_fishage_like_vals, 0)
    if(FishAgeComps_LikeType[f] == "Dirichlet-Multinomial") comp_fishage_like_vals <- c(comp_fishage_like_vals, 1)
    if(FishAgeComps_LikeType[f] == "iid-Logistic-Normal") comp_fishage_like_vals <- c(comp_fishage_like_vals, 2)
    if(FishAgeComps_LikeType[f] == "1d-Logistic-Normal") comp_fishage_like_vals <- c(comp_fishage_like_vals, 3)
    if(FishAgeComps_LikeType[f] == "2d-Logistic-Normal") comp_fishage_like_vals <- c(comp_fishage_like_vals, 4)
    collect_message(paste("Fishery Age Composition Likelihoods", "for fishery fleet", f, "specified as:" , FishAgeComps_LikeType[f]))
  } # end f loop

  # Specifying composition type
  FishAgeComps_Type_Mat <- array(NA, dim = c(length(input_list$data$years), input_list$data$n_fish_fleets))
  for(i in 1:length(FishAgeComps_Type)) {

    # Extract out components from list
    tmp <- FishAgeComps_Type[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    comps_type_tmp <- tmp_vec[1] # get composition type
    fleet <- as.numeric(tmp_vec[5]) # extract fleet index

    # Checking character string
    if(!comps_type_tmp %in% c("agg", "spltRspltS", "spltRjntS", 'none')) stop("FishAgeComps_Type not specified correctly. Must be one of: agg, spltRspltS, spltRjntS, none")
    if(!fleet %in% c(1:input_list$data$n_fish_fleets)) stop("Invalid fleet specified for FishAgeComps_Type. This needs to be specified as CompType_Year_x-y_Fleet_x")

    # get year ranges
    if(!str_detect(tmp, "terminal")) { # if not terminal year
      year_range <- as.numeric(unlist(strsplit(tmp_vec[3], "-")))
      years <- year_range[1]:year_range[2] # get sequence of years
    } else { # if terminal year
      year_range <- unlist(strsplit(tmp_vec[3], '-'))[1] # get year range
      years <- as.numeric(year_range):length(input_list$data$years) # get sequence of years
    }

    # Composition type
    # define composition types
    if(comps_type_tmp == "agg") {
      if(comp_fishage_like_vals[fleet] == 4) stop("Age composition likelihood specified as 2d-Logistic-Normal, but composition type is aggregated. This is not valid.")
      comps_type_val <- 0
    }
    if(comps_type_tmp == "spltRspltS") comps_type_val <- 1
    if(comps_type_tmp == "spltRjntS") comps_type_val <- 2
    if(comps_type_tmp == "none") comps_type_val <- 999

    # input into matrix
    FishAgeComps_Type_Mat[years,fleet] <- comps_type_val
  } # end i

  if(any(is.na(FishAgeComps_Type_Mat))) stop("FishAgeComps_Type is returning an NA. Did you update the year range of FishAgeComps_Type?")

  # Specifying composition likelihood for population-specific data
  comp_fishage_pop_like_vals <- vector()
  for(f in 1:input_list$data$n_fish_fleets) {
    if(FishAgeComps_pop_LikeType[f] == 'none') comp_fishage_pop_like_vals <- c(comp_fishage_pop_like_vals, 999)
    if(FishAgeComps_pop_LikeType[f] == "Multinomial") comp_fishage_pop_like_vals <- c(comp_fishage_pop_like_vals, 0)
    if(FishAgeComps_pop_LikeType[f] == "Dirichlet-Multinomial") comp_fishage_pop_like_vals <- c(comp_fishage_pop_like_vals, 1)
    if(FishAgeComps_pop_LikeType[f] == "iid-Logistic-Normal") comp_fishage_pop_like_vals <- c(comp_fishage_pop_like_vals, 2)
    if(FishAgeComps_pop_LikeType[f] == "1d-Logistic-Normal") comp_fishage_pop_like_vals <- c(comp_fishage_pop_like_vals, 3)
    if(FishAgeComps_pop_LikeType[f] == "2d-Logistic-Normal") comp_fishage_pop_like_vals <- c(comp_fishage_pop_like_vals, 4)
    collect_message(paste("Population Fishery Age Composition Likelihoods", "for fishery fleet", f, "specified as:" , FishAgeComps_pop_LikeType[f]))
  } # end f loop

  # Specifying composition type
  FishAgeComps_pop_Type_Mat <- array(NA, dim = c(length(input_list$data$years), input_list$data$n_fish_fleets))
  for(i in 1:length(FishAgeComps_pop_Type)) {

    # Extract out components from list
    tmp <- FishAgeComps_pop_Type[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    comps_type_tmp <- tmp_vec[1] # get composition type
    fleet <- as.numeric(tmp_vec[5]) # extract fleet index

    # Checking character string
    if(!comps_type_tmp %in% c("agg", "spltRspltS", "spltRjntS", 'none')) stop("FishAgeComps_pop_Type not specified correctly. Must be one of: agg, spltRspltS, spltRjntS, none")
    if(!fleet %in% c(1:input_list$data$n_fish_fleets)) stop("Invalid fleet specified for FishAgeComps_pop_Type. This needs to be specified as CompType_Year_x-y_Fleet_x")

    # get year ranges
    if(!str_detect(tmp, "terminal")) { # if not terminal year
      year_range <- as.numeric(unlist(strsplit(tmp_vec[3], "-")))
      years <- year_range[1]:year_range[2] # get sequence of years
    } else { # if terminal year
      year_range <- unlist(strsplit(tmp_vec[3], '-'))[1] # get year range
      years <- as.numeric(year_range):length(input_list$data$years) # get sequence of years
    }

    # Composition type
    # define composition types
    if(comps_type_tmp == "agg") {
      if(comp_fishage_pop_like_vals[fleet] == 4) stop("Population Age composition likelihood specified as 2d-Logistic-Normal, but composition type is aggregated. This is not valid.")
      comps_type_val <- 0
    }
    if(comps_type_tmp == "spltRspltS") comps_type_val <- 1
    if(comps_type_tmp == "spltRjntS") comps_type_val <- 2
    if(comps_type_tmp == "none") comps_type_val <- 999

    # input into matrix
    FishAgeComps_pop_Type_Mat[years,fleet] <- comps_type_val
  } # end i

  if(any(is.na(FishAgeComps_pop_Type_Mat))) stop("FishAgeComps_pop_Type is returning an NA. Did you update the year range of FishAgeComps_pop_Type?")

  # Fishery Length Composition Options -----------------------------------------

  # Specifying composition likelihood
  comp_fishlen_like_vals <- vector()
  for(f in 1:input_list$data$n_fish_fleets) {
    if(FishLenComps_LikeType[f] == 'none') comp_fishlen_like_vals <- c(comp_fishlen_like_vals, 999)
    if(FishLenComps_LikeType[f] == "Multinomial") comp_fishlen_like_vals <- c(comp_fishlen_like_vals, 0)
    if(FishLenComps_LikeType[f] == "Dirichlet-Multinomial") comp_fishlen_like_vals <- c(comp_fishlen_like_vals, 1)
    if(FishLenComps_LikeType[f] == "iid-Logistic-Normal") comp_fishlen_like_vals <- c(comp_fishlen_like_vals, 2)
    if(FishLenComps_LikeType[f] == "1d-Logistic-Normal") comp_fishlen_like_vals <- c(comp_fishlen_like_vals, 3)
    if(FishLenComps_LikeType[f] == "2d-Logistic-Normal") comp_fishlen_like_vals <- c(comp_fishlen_like_vals, 4)
    collect_message(paste("Fishery Length Composition Likelihoods", "for fishery fleet", f, "specified as:" , FishLenComps_LikeType[f]))
  } # end f loop

  # Specifying composition type
  FishLenComps_Type_Mat <- array(NA, dim = c(length(input_list$data$years), input_list$data$n_fish_fleets))
  for(i in 1:length(FishLenComps_Type)) {

    # Extract out components from list
    tmp <- FishLenComps_Type[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    comps_type_tmp <- tmp_vec[1] # get composition type
    fleet <- as.numeric(tmp_vec[5]) # extract fleet index

    # define composition types
    if(!comps_type_tmp %in% c("agg", "spltRspltS", "spltRjntS", 'none')) stop("FishLenComps_Type not specified correctly. Must be one of: agg, spltRspltS, spltRjntS, none")
    if(!fleet %in% c(1:input_list$data$n_fish_fleets)) stop("Invalid fleet specified for FishLenComps_Type This needs to be specified as CompType_Year_x-y_Fleet_x")

    # get year ranges
    if(!str_detect(tmp, "terminal")) { # if not terminal year
      year_range <- as.numeric(unlist(strsplit(tmp_vec[3], "-")))
      years <- year_range[1]:year_range[2] # get sequence of years
    } else { # if terminal year
      year_range <- unlist(strsplit(tmp_vec[3], '-'))[1] # get year range
      years <- as.numeric(year_range):length(input_list$data$years) # get sequence of years
    }

    # define composition types
    if(comps_type_tmp == "agg") {
      if(comp_fishlen_like_vals[fleet] == 4) stop("Length composition likelihood specified as 2d-Logistic-Normal, but composition type is aggregated. This is not valid.")
      comps_type_val <- 0
    }
    if(comps_type_tmp == "spltRspltS") comps_type_val <- 1
    if(comps_type_tmp == "spltRjntS") comps_type_val <- 2
    if(comps_type_tmp == "none") comps_type_val <- 999

    # input into matrix
    FishLenComps_Type_Mat[years,fleet] <- comps_type_val
  } # end i

  if(any(is.na(FishLenComps_Type_Mat))) stop("FishLenComps_Type_Mat is returning an NA. Did you update the year range of FishLenComps_Type_Mat?")


  # Specifying composition likelihood for population-specific data
  comp_fishlen_pop_like_vals <- vector()
  for(f in 1:input_list$data$n_fish_fleets) {
    if(FishLenComps_pop_LikeType[f] == 'none') comp_fishlen_pop_like_vals <- c(comp_fishlen_pop_like_vals, 999)
    if(FishLenComps_pop_LikeType[f] == "Multinomial") comp_fishlen_pop_like_vals <- c(comp_fishlen_pop_like_vals, 0)
    if(FishLenComps_pop_LikeType[f] == "Dirichlet-Multinomial") comp_fishlen_pop_like_vals <- c(comp_fishlen_pop_like_vals, 1)
    if(FishLenComps_pop_LikeType[f] == "iid-Logistic-Normal") comp_fishlen_pop_like_vals <- c(comp_fishlen_pop_like_vals, 2)
    if(FishLenComps_pop_LikeType[f] == "1d-Logistic-Normal") comp_fishlen_pop_like_vals <- c(comp_fishlen_pop_like_vals, 3)
    if(FishLenComps_pop_LikeType[f] == "2d-Logistic-Normal") comp_fishlen_pop_like_vals <- c(comp_fishlen_pop_like_vals, 4)
    collect_message(paste("Population Fishery Length Composition Likelihoods", "for fishery fleet", f, "specified as:" , FishLenComps_pop_LikeType[f]))
  } # end f loop

  # Specifying composition type
  FishLenComps_pop_Type_Mat <- array(NA, dim = c(length(input_list$data$years), input_list$data$n_fish_fleets))
  for(i in 1:length(FishLenComps_pop_Type)) {

    # Extract out components from list
    tmp <- FishLenComps_pop_Type[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    comps_type_tmp <- tmp_vec[1] # get composition type
    fleet <- as.numeric(tmp_vec[5]) # extract fleet index

    # Checking character string
    if(!comps_type_tmp %in% c("agg", "spltRspltS", "spltRjntS", 'none')) stop("FishLenComps_pop_Type not specified correctly. Must be one of: agg, spltRspltS, spltRjntS, none")
    if(!fleet %in% c(1:input_list$data$n_fish_fleets)) stop("Invalid fleet specified for FishLenComps_pop_Type. This needs to be specified as CompType_Year_x-y_Fleet_x")

    # get year ranges
    if(!str_detect(tmp, "terminal")) { # if not terminal year
      year_range <- as.numeric(unlist(strsplit(tmp_vec[3], "-")))
      years <- year_range[1]:year_range[2] # get sequence of years
    } else { # if terminal year
      year_range <- unlist(strsplit(tmp_vec[3], '-'))[1] # get year range
      years <- as.numeric(year_range):length(input_list$data$years) # get sequence of years
    }

    # Composition type
    # define composition types
    if(comps_type_tmp == "agg") {
      if(comp_fishlen_pop_like_vals[fleet] == 4) stop("Population Len composition likelihood specified as 2d-Logistic-Normal, but composition type is aggregated. This is not valid.")
      comps_type_val <- 0
    }
    if(comps_type_tmp == "spltRspltS") comps_type_val <- 1
    if(comps_type_tmp == "spltRjntS") comps_type_val <- 2
    if(comps_type_tmp == "none") comps_type_val <- 999

    # input into matrix
    FishLenComps_pop_Type_Mat[years,fleet] <- comps_type_val
  } # end i

  if(any(is.na(FishLenComps_pop_Type_Mat))) stop("FishLenComps_pop_Type is returning an NA. Did you update the year range of FishLenComps_pop_Type?")

  # ISS Munging -------------------------------------------------------------

  # Fishery Ages
  if(is.null(ISS_FishAgeComps)) {
    collect_message("No ISS is specified for FishAgeComps. ISS weighting is calculated by summing up values from ObsFishAgeComps each year")
    ISS_FishAgeComps <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_fish_fleets))
    for(y in 1:length(input_list$data$years)) {
      for(f in 1:input_list$data$n_fish_fleets) {
        for(seas in 1:input_list$data$n_seas) {
          # if aggregated across sexes and regions (0) or joint across sexes
          if(FishAgeComps_Type_Mat[y,f] == 0) ISS_FishAgeComps[1,y,seas,1,f] <- sum(ObsFishAgeComps[,y,seas,,,f])
          # if split by region and sex
          if(FishAgeComps_Type_Mat[y,f] == 1) ISS_FishAgeComps[,y,seas,,f] <- apply(ObsFishAgeComps[,y,seas,,,f, drop = FALSE], c(1,4), sum)
          # if split by region, joint by sex
          if(FishAgeComps_Type_Mat[y,f] == 2) ISS_FishAgeComps[,y,seas,1,f] <- apply(ObsFishAgeComps[,y,seas,,,f, drop = FALSE], 1, sum)
        } # end seas loop
      } # end f loop
    } # end y loop
  }

  # Fishery Lengths
  if(is.null(ISS_FishLenComps)) {
    collect_message("No ISS is specified for FishLenComps. ISS weighting is calculated by summing up values from ObsFishLenComps each year")
    ISS_FishLenComps <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_fish_fleets))
    for(y in 1:length(input_list$data$years)) {
      for(f in 1:input_list$data$n_fish_fleets) {
        for(seas in 1:input_list$data$n_seas) {
          # if aggregated across sexes and regions (0)
          if(FishLenComps_Type_Mat[y,f] == 0) ISS_FishLenComps[1,y,seas,1,f] <- sum(ObsFishLenComps[,y,seas,,,f])
          # if split by region and sex
          if(FishLenComps_Type_Mat[y,f] == 1) ISS_FishLenComps[,y,seas,,f] <- apply(ObsFishLenComps[,y,seas,,,f, drop = FALSE], c(1,4), sum)
          # if split by region, joint by sex
          if(FishLenComps_Type_Mat[y,f] == 2) ISS_FishLenComps[,y,seas,1,f] <- apply(ObsFishLenComps[,y,seas,,,f, drop = FALSE], 1, sum)
        } # end seas loop
      } # end f loop
    } # end y loop
  }

  # Fishery Ages
  if(is.null(ISS_FishAgeComps_pop)) {
    collect_message("No ISS is specified for pop_FishAgeComps. ISS weighting is calculated by summing up values from ObsFishAgeComps_pop each year")
    ISS_FishAgeComps_pop <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_fish_fleets))
    for(p in 1:input_list$data$n_pop) {
      for(y in 1:length(input_list$data$years)) {
        for(f in 1:input_list$data$n_fish_fleets) {
          for(seas in 1:input_list$data$n_seas) {
            # if aggregated across sexes and regions (0) or joint across sexes
            if(FishAgeComps_pop_Type_Mat[y,f] == 0) ISS_FishAgeComps_pop[p,1,y,seas,1,f] <- sum(ObsFishAgeComps_pop[p,,y,seas,,,f])
            # if split by region and sex
            if(FishAgeComps_pop_Type_Mat[y,f] == 1) ISS_FishAgeComps_pop[p,,y,seas,,f] <- apply(ObsFishAgeComps_pop[p,,y,seas,,,f, drop = FALSE], c(2,5), sum)
            # if split by region, joint by sex
            if(FishAgeComps_pop_Type_Mat[y,f] == 2) ISS_FishAgeComps_pop[p,,y,seas,1,f] <- apply(ObsFishAgeComps_pop[p,,y,seas,,,f, drop = FALSE], 2, sum)
          } # end seas loop
        } # end f loop
      } # end y loop
    } # end p loop
  }

  # Fishery Lengths
  if(is.null(ISS_FishLenComps_pop)) {
    collect_message("No ISS is specified for pop_FishLenComps. ISS weighting is calculated by summing up values from ObsFishLenComps_pop each year")
    ISS_FishLenComps_pop <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_fish_fleets))
    for(p in 1:input_list$data$n_pop) {
      for(y in 1:length(input_list$data$years)) {
        for(f in 1:input_list$data$n_fish_fleets) {
          for(seas in 1:input_list$data$n_seas) {
            # if aggregated across sexes and regions (0)
            if(FishLenComps_pop_Type_Mat[y,f] == 0) ISS_FishLenComps_pop[p,1,y,seas,1,f] <- sum(ObsFishLenComps_pop[p,,y,seas,,,f])
            # if split by region and sex
            if(FishLenComps_pop_Type_Mat[y,f] == 1) ISS_FishLenComps_pop[p,,y,seas,,f] <- apply(ObsFishLenComps_pop[p,,y,seas,,,f, drop = FALSE], c(2,5), sum)
            # if split by region, joint by sex
            if(FishLenComps_pop_Type_Mat[y,f] == 2) ISS_FishLenComps_pop[p,,y,seas,1,f] <- apply(ObsFishLenComps_pop[p,,y,seas,,,f, drop = FALSE], 2, sum)
          } # end seas loop
        } # end f loop
      } # end y loop
    } # end p loop
  }

  # Populate Data List ------------------------------------------------------

  input_list$data$ISS_FishAgeComps <- ISS_FishAgeComps
  input_list$data$ISS_FishLenComps <- ISS_FishLenComps
  input_list$data$ISS_FishAgeComps_pop <- ISS_FishAgeComps_pop
  input_list$data$ISS_FishLenComps_pop <- ISS_FishLenComps_pop
  input_list$data$ObsFishIdx <- ObsFishIdx
  input_list$data$ObsFishIdx_SE <- ObsFishIdx_SE
  input_list$data$UseFishIdx <- UseFishIdx
  input_list$data$ObsFishIdx_pop <- ObsFishIdx_pop
  input_list$data$ObsFishIdx_pop_SE <- ObsFishIdx_pop_SE
  input_list$data$UseFishIdx_pop <- UseFishIdx_pop
  input_list$data$fish_idx_type <- fish_idx_type_vals
  input_list$data$ObsFishAgeComps <- ObsFishAgeComps
  input_list$data$UseFishAgeComps <- UseFishAgeComps
  input_list$data$ObsFishLenComps <- ObsFishLenComps
  input_list$data$UseFishLenComps <- UseFishLenComps
  input_list$data$ObsFishAgeComps_pop <- ObsFishAgeComps_pop
  input_list$data$UseFishAgeComps_pop <- UseFishAgeComps_pop
  input_list$data$ObsFishLenComps_pop <- ObsFishLenComps_pop
  input_list$data$UseFishLenComps_pop <- UseFishLenComps_pop
  input_list$data$FishAgeComps_LikeType <- comp_fishage_like_vals
  input_list$data$FishLenComps_LikeType <- comp_fishlen_like_vals
  input_list$data$FishAgeComps_pop_LikeType <- comp_fishage_pop_like_vals
  input_list$data$FishLenComps_pop_LikeType <- comp_fishlen_pop_like_vals
  input_list$data$FishAgeComps_Type <- FishAgeComps_Type_Mat
  input_list$data$FishLenComps_Type <- FishLenComps_Type_Mat
  input_list$data$FishAgeComps_pop_Type <- FishAgeComps_pop_Type_Mat
  input_list$data$FishLenComps_pop_Type <- FishLenComps_pop_Type_Mat

  # Populate Parameter List -------------------------------------------------

  # Dispersion parameters for the fishery age comps
  if("ln_FishAge_theta" %in% names(starting_values)) input_list$par$ln_FishAge_theta <- starting_values$ln_FishAge_theta
  else input_list$par$ln_FishAge_theta <- array(0, dim = c(input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_fish_fleets))

  # logistic normal correlation parameters for fishery age comps
  if("FishAge_corr_pars" %in% names(starting_values)) input_list$par$FishAge_corr_pars <- starting_values$FishAge_corr_pars
  else input_list$par$FishAge_corr_pars <- array(0.01, dim = c(input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_fish_fleets, 2))

  # aggregated
  if("ln_FishAge_theta_agg" %in% names(starting_values)) input_list$par$ln_FishAge_theta_agg <- starting_values$ln_FishAge_theta_agg
  else input_list$par$ln_FishAge_theta_agg <- array(0, dim = c(input_list$data$n_fish_fleets))

  # aggregated correlation parameters
  if("FishAge_corr_pars_agg" %in% names(starting_values)) input_list$par$FishAge_corr_pars_agg <- starting_values$FishAge_corr_pars_agg
  else input_list$par$FishAge_corr_pars_agg <- array(0.01, dim = c(input_list$data$n_fish_fleets))

  # Dispersion parameters for fishery length comps
  if("ln_FishLen_theta" %in% names(starting_values)) input_list$par$ln_FishLen_theta <- starting_values$ln_FishLen_theta
  else input_list$par$ln_FishLen_theta <- array(0, dim = c(input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_fish_fleets))

  # logistic normal correlation parameters for fishery length comps
  if("FishLen_corr_pars" %in% names(starting_values)) input_list$par$FishLen_corr_pars <- starting_values$FishLen_corr_pars
  else input_list$par$FishLen_corr_pars <- array(0.01, dim = c(input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_fish_fleets, 2))

  # aggregated
  if("ln_FishLen_theta_agg" %in% names(starting_values)) input_list$par$ln_FishLen_theta_agg <- starting_values$ln_FishLen_theta_agg
  else input_list$par$ln_FishLen_theta_agg <- array(0, dim = c(input_list$data$n_fish_fleets))

  if("FishLen_corr_pars_agg" %in% names(starting_values)) input_list$par$FishLen_corr_pars_agg <- starting_values$FishLen_corr_pars_agg
  else input_list$par$FishLen_corr_pars_agg <- array(0.01, dim = c(input_list$data$n_fish_fleets))

  # Dispersion parameters for the population fishery age comps
  if("ln_FishAge_pop_theta" %in% names(starting_values)) input_list$par$ln_FishAge_pop_theta <- starting_values$ln_FishAge_pop_theta
  else input_list$par$ln_FishAge_pop_theta <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_fish_fleets))

  # logistic normal correlation parameters for population fishery age comps
  if("FishAge_pop_corr_pars" %in% names(starting_values)) input_list$par$FishAge_pop_corr_pars <- starting_values$FishAge_pop_corr_pars
  else input_list$par$FishAge_pop_corr_pars <- array(0.01, dim = c(input_list$data$n_pop,input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_fish_fleets, 2))

  # aggregated population pars
  if("ln_FishAge_pop_theta_agg" %in% names(starting_values)) input_list$par$ln_FishAge_pop_theta_agg <- starting_values$ln_FishAge_pop_theta_agg
  else input_list$par$ln_FishAge_pop_theta_agg <- array(0, dim = c(input_list$data$n_pop,input_list$data$n_fish_fleets))

  # aggregated population correlation parameters
  if("FishAge_pop_corr_pars_agg" %in% names(starting_values)) input_list$par$FishAge_pop_corr_pars_agg <- starting_values$FishAge_pop_corr_pars_agg
  else input_list$par$FishAge_pop_corr_pars_agg <- array(0.01, dim = c(input_list$data$n_pop,input_list$data$n_fish_fleets))

  # Dispersion parameters for population fishery length comps
  if("ln_FishLen_pop_theta" %in% names(starting_values)) input_list$par$ln_FishLen_pop_theta <- starting_values$ln_FishLen_pop_theta
  else input_list$par$ln_FishLen_pop_theta <- array(0, dim = c(input_list$data$n_pop,input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_fish_fleets))

  # logistic normal correlation parameters for population fishery length comps
  if("FishLen_pop_corr_pars" %in% names(starting_values)) input_list$par$FishLen_pop_corr_pars <- starting_values$FishLen_pop_corr_pars
  else input_list$par$FishLen_pop_corr_pars <- array(0.01, dim = c(input_list$data$n_pop,input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_fish_fleets, 2))

  # aggregated population pars
  if("ln_FishLen_pop_theta_agg" %in% names(starting_values)) input_list$par$ln_FishLen_pop_theta_agg <- starting_values$ln_FishLen_pop_theta_agg
  else input_list$par$ln_FishLen_pop_theta_agg <- array(0, dim = c(input_list$data$n_pop,input_list$data$n_fish_fleets))

  if("FishLen_pop_corr_pars_agg" %in% names(starting_values)) input_list$par$FishLen_pop_corr_pars_agg <- starting_values$FishLen_pop_corr_pars_agg
  else input_list$par$FishLen_pop_corr_pars_agg <- array(0.01, dim = c(input_list$data$n_pop, input_list$data$n_fish_fleets))

  # Mapping Options ---------------------------------------------------------

  input_list <- do_FishAge_theta_mapping(input_list)
  input_list <- do_FishLen_theta_mapping(input_list)
  input_list <- do_FishAge_corr_pars_mapping(input_list)
  input_list <- do_FishLen_corr_pars_mapping(input_list)

  input_list <- do_FishAge_pop_theta_mapping(input_list)
  input_list <- do_FishLen_pop_theta_mapping(input_list)
  input_list <- do_FishAge_pop_corr_pars_mapping(input_list)
  input_list <- do_FishLen_pop_corr_pars_mapping(input_list)

  # Discard Compositions (forwarded to Setup_Mod_Discard_Comps) ---------------
  input_list <- Setup_Mod_Discard_Comps(
    input_list,
    ObsFishAgeComps_discard        = ObsFishAgeComps_discard,
    UseFishAgeComps_discard        = UseFishAgeComps_discard,
    ISS_FishAgeComps_discard       = ISS_FishAgeComps_discard,
    ObsFishLenComps_discard        = ObsFishLenComps_discard,
    UseFishLenComps_discard        = UseFishLenComps_discard,
    ISS_FishLenComps_discard       = ISS_FishLenComps_discard,
    FishAgeComps_discard_LikeType  = FishAgeComps_discard_LikeType,
    FishLenComps_discard_LikeType  = FishLenComps_discard_LikeType,
    FishAgeComps_discard_Type      = FishAgeComps_discard_Type,
    FishLenComps_discard_Type      = FishLenComps_discard_Type,
    ObsFishAgeComps_discard_pop    = ObsFishAgeComps_discard_pop,
    UseFishAgeComps_discard_pop    = UseFishAgeComps_discard_pop,
    ISS_FishAgeComps_discard_pop   = ISS_FishAgeComps_discard_pop,
    ObsFishLenComps_discard_pop    = ObsFishLenComps_discard_pop,
    UseFishLenComps_discard_pop    = UseFishLenComps_discard_pop,
    ISS_FishLenComps_discard_pop   = ISS_FishLenComps_discard_pop,
    FishAgeComps_discard_pop_LikeType = FishAgeComps_discard_pop_LikeType,
    FishLenComps_discard_pop_LikeType = FishLenComps_discard_pop_LikeType,
    FishAgeComps_discard_pop_Type     = FishAgeComps_discard_pop_Type,
    FishLenComps_discard_pop_Type     = FishLenComps_discard_pop_Type,
    ...
  )

  # Print Messages ----------------------------------------------------------
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)
}

#' Map fishery selectivity fixed-effect parameters
#'
#' Constructs the factor map for \code{fish_fixed_sel_pars} (e.g., \eqn{a_{50}}, \eqn{k}, \eqn{a_{max}}),
#' controlling whether selectivity shape parameters are estimated independently or shared across
#' regions, sexes, or fleets. Cells with no catch data (\code{UseCatch == 0}) are automatically
#' mapped to \code{NA}.
#'
#' Non-parametric selectivity (model \code{sel_model == 5}) is supported via
#' \code{fish_sel_nonpar_est_bins}, which defines bin groupings that are treated as
#' individual estimated parameters.
#'
#' Fleet sharing (\code{"est_shared_f_x"}) is handled in a second pass after all
#' base fleet mappings are established, copying the reference fleet's index
#' assignments into the sharing fleet.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists.
#'
#' @param fish_fixed_sel_pars_spec Character vector of length \code{n_fish_fleets}.
#' Each element specifies the estimation structure for one fleet. Options:
#' \describe{
#'   \item{\code{"est_all"}}{Separate parameters per region × sex × block × bin.}
#'   \item{\code{"est_shared_r"}}{Parameters shared across regions; unique per sex × block × bin.}
#'   \item{\code{"est_shared_s"}}{Parameters shared across sexes; unique per region × block × bin.}
#'   \item{\code{"est_shared_r_s"}}{Parameters shared across regions and sexes; unique per block × bin.}
#'   \item{\code{"est_shared_f_x"}}{Copy parameter mapping from fleet \code{x}. Fleet \code{x}
#'     must not itself use \code{"est_shared_f_y"}.}
#'   \item{\code{"fix"}}{All parameters fixed at starting values (mapped to \code{NA}).}
#' }
#'
#' @param bins Number of selectivity bins.
#'
#' @param fish_sel_nonpar_est_bins Optional list defining bin groupings for
#' non-parametric selectivity. Structure is \code{[[fleet]][[block]]}, where each
#' element is a list of bin index vectors representing grouped parameters.
#'
#' @return Updated \code{input_list} with \code{$map$fish_fixed_sel_pars}
#'   as a factor array.
#'
#' @keywords internal
do_fish_fixed_sel_pars_mapping <- function(input_list, fish_fixed_sel_pars_spec, bins, fish_sel_nonpar_est_bins) {

  # Initialize counter and mapping array for fixed effects fishery selectivity
  fish_fixed_sel_pars_counter <- 1
  map_fish_fixed_sel_pars <- input_list$par$fish_fixed_sel_pars
  map_fish_fixed_sel_pars[] <- NA

  for(f in 1:input_list$data$n_fish_fleets) {

    # Validate Options
    if(!fish_fixed_sel_pars_spec[f] %in% c("est_all", "est_shared_r", "est_shared_r_s", "fix", "est_shared_s", "fix_fish_sel_input") &&
       !stringr::str_detect(fish_fixed_sel_pars_spec[f], "est_shared_f_\\d+"))
      stop("fish_fixed_sel_pars_spec not correctly specfied. Should be one of these: est_all, est_shared_r, est_shared_r_s, est_shared_s, fix, or est_shared_f_# (where # is fleet number)")
    # checking fixed selex options
    if(input_list$data$use_fixed_fish_sel[f] == 1 && stringr::str_detect(fish_fixed_sel_pars_spec[f], 'est'))
      stop("use_fixed_fish_sel has 1s for a given fleet, but fish_fixed_sel_pars_spec is specified at an est variant.")
    if(input_list$data$use_fixed_fish_sel[f] == 0 && fish_fixed_sel_pars_spec[f] == 'fix_fish_sel_input')
      stop("use_fixed_fish_sel has 0s for a given fleet, but fish_fixed_sel_pars_spec is specified at fix_fish_sel_input")

    # Skip fleet sharing specs in first pass
    if(stringr::str_detect(fish_fixed_sel_pars_spec[f], "est_shared_f")) next

    for(r in 1:input_list$data$n_regions) {

      # Only add a counter if caatches are avaliable in some years for a given region and fleet combination
      if(sum(input_list$data$UseCatch[r,,,f]) > 0 || sum(input_list$data$UseCatch_pop[,r,,,f]) > 0) {

        # Extract number of fishery selectivity blocks
        fishsel_blocks_tmp <- unique(as.vector(input_list$data$fish_sel_blocks[r,,f]))

        for(s in 1:input_list$data$n_sexes) {
          for(b in 1:length(fishsel_blocks_tmp)) {

            block_years <- which(input_list$data$fish_sel_blocks[r,,f] == fishsel_blocks_tmp[b]) # figure out block years
            sel_model_this_block <- unique(input_list$data$fish_sel_model[r, block_years, f]) # get selectivity form for a given block
            if(length(sel_model_this_block) > 1) stop("Block ", fishsel_blocks_tmp[b], " for fleet ", f, " region ", r, " has multiple selectivity models assigned to it")

            # determine maximum selectivity parameters
            if(sel_model_this_block == 2) max_sel_pars <- 1 # exponential
            if(sel_model_this_block %in% c(0,1,3)) max_sel_pars <- 2 # logistic or gamma
            if(sel_model_this_block == 4) max_sel_pars <- 6 # double normal
            if(sel_model_this_block %in% c(6,7)) max_sel_pars <- 3 # logistic w/ asymptotic selectivity
            if(sel_model_this_block == 8) { # bicubic spline: flattened bin-node x year-node grid (group_bins below reduces to a plain 1:max_sel_pars mapping, same as other parametric forms)
              n_bin_nodes_this <- unique(input_list$data$fish_sel_bicubic_binnodes[r, block_years, f])
              n_yr_nodes_this <- unique(input_list$data$fish_sel_bicubic_yrnodes[r, block_years, f])
              max_sel_pars <- n_bin_nodes_this * n_yr_nodes_this
            }

            # non-parametric selectivity
            if(sel_model_this_block == 5) {

              if(is.null(fish_sel_nonpar_est_bins)) stop("Non-parametric fishery selectivtiy specified, but fish_sel_nonpar_est_bins is NULL. Please specify bins!")
              bin_groups <- fish_sel_nonpar_est_bins[[f]][[b]]
              max_sel_pars <- length(bin_groups)  # number of groups = number of estimated pars

              # validate
              all_bins <- unlist(bin_groups)
              if(any(all_bins < 1) || any(all_bins > bins))
                stop("fish_sel_nonpar_est_bins[[", f, "]][[", b, "]] contains indices outside 1:", bins)
              if(length(all_bins) != length(unique(all_bins)))
                stop("fish_sel_nonpar_est_bins[[", f, "]][[", b, "]] has duplicate bin indices")
            }

            for(i in 1:max_sel_pars) {

              # get non-parametric selectivity bins
              group_bins <- if(sel_model_this_block == 5) bin_groups[[i]] else i

              # Estimate all selectivity fixed effects parameters within the constraints of the defined blocks
              if(fish_fixed_sel_pars_spec[f] == "est_all") {
                for(bi in group_bins) map_fish_fixed_sel_pars[r,bi,b,s,f] <- fish_fixed_sel_pars_counter
                fish_fixed_sel_pars_counter <- fish_fixed_sel_pars_counter + 1
              } # end if

              # Estimating parameters shared across regions (but unique for each sex, fleet, parameter)
              if(fish_fixed_sel_pars_spec[f] == 'est_shared_r' && r == 1) {
                for(rr in 1:input_list$data$n_regions) {
                  if(fishsel_blocks_tmp[b] %in% input_list$data$fish_sel_blocks[rr,,f]) {
                    for(bi in group_bins) map_fish_fixed_sel_pars[rr, bi, b, s, f] <- fish_fixed_sel_pars_counter
                  } # end if
                } # end rr loop
                fish_fixed_sel_pars_counter <- fish_fixed_sel_pars_counter + 1
              } # end if

              # Estimating process error parameters shared across sexes (but unique for each region, fleet, parameter)
              if(fish_fixed_sel_pars_spec[f] == 'est_shared_s' && s == 1) {
                for(ss in 1:input_list$data$n_sexes) {
                  for(bi in group_bins) map_fish_fixed_sel_pars[r, bi, b, ss, f] <- fish_fixed_sel_pars_counter
                } # end ss loop
                fish_fixed_sel_pars_counter <- fish_fixed_sel_pars_counter + 1
              } # end if

              # Estimating process error parameters shared across regions and sexes (but unique for each fleet, parameter)
              if(fish_fixed_sel_pars_spec[f] == 'est_shared_r_s' && r == 1 && s == 1) {
                for(rr in 1:input_list$data$n_regions) {
                  for(ss in 1:input_list$data$n_sexes) {
                    if(fishsel_blocks_tmp[b] %in% input_list$data$fish_sel_blocks[rr,,f]) {
                      for(bi in group_bins) map_fish_fixed_sel_pars[rr, bi, b, ss, f] <- fish_fixed_sel_pars_counter
                    } # end if
                  } # end ss loop
                } #end rr loop
                fish_fixed_sel_pars_counter <- fish_fixed_sel_pars_counter + 1
              } # end if

            } # end i loop
          } # end b loop
        } # end s loop
      } # end if statement
    } # end r loop

    # fix all parameters
    if(fish_fixed_sel_pars_spec[f] %in% c("fix", "fixed_fish_sel_input")) map_fish_fixed_sel_pars[,,,,f] <- NA
    collect_message("fish_fixed_sel_pars_spec is specified as: ", fish_fixed_sel_pars_spec[f], " for fishery fleet ", f)

  } # end f loop

  # Handle fleet sharing after all base mappings are established
  for(f in 1:input_list$data$n_fish_fleets) {
    if(stringr::str_detect(fish_fixed_sel_pars_spec[f], "est_shared_f")) {
      # extract fleet sharing index
      flt_shared <- as.numeric(unlist(stringr::str_extract_all(fish_fixed_sel_pars_spec[f], "\\d+")))

      # Validate options here
      if(flt_shared > input_list$data$n_fish_fleets || flt_shared < 1) stop("Fleet sharing specification 'est_shared_f", flt_shared, "' for fleet ", f, " references invalid fleet number.")
      if(stringr::str_detect(fish_fixed_sel_pars_spec[flt_shared], "est_shared_f")) stop("Fleet ", f, " cannot share with fleet ", flt_shared, " because fleet ", flt_shared, " is self-sharing parameters, which does not make sense.")

      # Copy mapping from reference fleet
      map_fish_fixed_sel_pars[,,,,f] <- map_fish_fixed_sel_pars[,,,,flt_shared]
      collect_message("fish_fixed_sel_pars_spec is specified as: ", fish_fixed_sel_pars_spec[f], " for fishery fleet ", f, " (sharing with fleet ", flt_shared, ")")
    } # end if statement
  } # end f loop

  # input into mapping list
  input_list$map$fish_fixed_sel_pars <- factor(map_fish_fixed_sel_pars)
  return(input_list)
}


#' Map fishery catchability parameters
#'
#' Constructs the factor map for \code{ln_fish_q}, controlling whether
#' catchability parameters are estimated independently per region and time block
#' or shared across regions. Cells with no fishery index observations
#' (\code{UseFishIdx == 0}) are automatically mapped to \code{NA}.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists.
#' @param fish_q_spec Character vector of length \code{n_fish_fleets}. Options:
#'   \describe{
#'     \item{\code{"est_all"}}{Separate catchability per region × block × fleet.}
#'     \item{\code{"est_shared_r"}}{Single catchability shared across regions,
#'       unique per block × fleet.}
#'     \item{\code{"fix"}}{All catchability parameters fixed (mapped to \code{NA}).}
#'   }
#'
#' @return The input \code{input_list} with \code{$map$ln_fish_q} set to a
#'   factor vector.
#'
#' @keywords internal
do_fish_q_mapping <- function(input_list, fish_q_spec) {

  # Initialize counter and mapping array for fishery catchability
  fish_q_counter <- 1
  map_fish_q <- input_list$par$ln_fish_q
  map_fish_q[] <- NA

  for(f in 1:input_list$data$n_fish_fleets) {

    # Validate options
    if(!is.null(fish_q_spec)) {
      if(!fish_q_spec[f] %in% c("est_all", "est_shared_r", "fix"))
        stop("fish_q_spec not correctly specfied. Should be one of these: est_all, est_shared_r, fix")
    }

    for(r in 1:input_list$data$n_regions) {

      if(sum(input_list$data$UseFishIdx[r,,,f]) == 0 && sum(input_list$data$UseFishIdx_pop[,r,,,f]) == 0) {
        map_fish_q[r,,f] <- NA # fix parameters if we are not using fishery indices for these fleets and regions
      } else {

        # Extract number of fishery catchability blocks
        fishq_blocks_tmp <- unique(as.vector(input_list$data$fish_q_blocks[r,,f]))

        for(b in 1:length(fishq_blocks_tmp)) {

          # Estimate for all regions
          if(fish_q_spec[f] == 'est_all') {
            map_fish_q[r,b,f] <- fish_q_counter
            fish_q_counter <- fish_q_counter + 1
          } # end if

          # Estimate but share q across regions
          if(fish_q_spec[f] == 'est_shared_r' && r == 1) {
            for(rr in 1:input_list$data$n_regions) {
              if(fishq_blocks_tmp[b] %in% input_list$data$fish_q_blocks[rr,,f]) {
                map_fish_q[rr, b, f] <- fish_q_counter
              } # end if
            } # end rr loop
            fish_q_counter <- fish_q_counter + 1
          } # end if

        } # end b loop
      } # end else loop
    } # end r loop

    # fix all parameters
    if(fish_q_spec[f] == 'fix') map_fish_q[,,f] <- NA
    collect_message("fish_q_spec is specified as: ", fish_q_spec[f], " for fishery fleet ", f)
  } # end f loop

  # input into mapping list
  input_list$map$ln_fish_q <- factor(map_fish_q)

  return(input_list)
}

#' Map fishery selectivity process error hyperparameters
#'
#' Constructs the factor map for \code{fishsel_pe_pars}, which contains the
#' variance and correlation hyperparameters governing continuous time-varying
#' selectivity. The set of active parameters depends on the time-variation type
#' (\code{cont_tv_fish_sel}): iid/random-walk forms use up to 2 parameters
#' (log-sigma); 3D GMRF forms use up to 4 (partial correlations for age, year,
#' cohort dimensions plus log-sigma); the 2D AR1 form uses 3 (bin AR1, year AR1,
#' log-sigma). Correlation components can be selectively suppressed via
#' \code{corr_opt_semipar}.
#'
#' Fleet sharing (\code{"est_shared_f_x"}) is handled in a second pass.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists.
#' @param fishsel_pe_pars_spec Character vector of length \code{n_fish_fleets}.
#'   Options:
#'   \describe{
#'     \item{\code{"est_all"}}{Separate hyperparameters per region × sex.}
#'     \item{\code{"est_shared_r"}}{Shared across regions; unique per sex.}
#'     \item{\code{"est_shared_s"}}{Shared across sexes; unique per region.}
#'     \item{\code{"est_shared_r_s"}}{Shared across regions and sexes.}
#'     \item{\code{"est_shared_f_x"}}{Copy hyperparameters from fleet \code{x}.}
#'     \item{\code{"fix"} or \code{"none"}}{All parameters fixed (mapped to \code{NA}).}
#'   }
#' @param corr_opt_semipar Character vector of length \code{n_fish_fleets}
#'   specifying which correlation components to suppress for semi-parametric
#'   models. Valid values per fleet: \code{NA} (no suppression),
#'   \code{"corr_zero_y"}, \code{"corr_zero_b"}, \code{"corr_zero_y_b"},
#'   \code{"corr_zero_c"}, \code{"corr_zero_y_c"}, \code{"corr_zero_b_c"},
#'   \code{"corr_zero_y_b_c"}. Cohort options (\code{"corr_zero_c"}, etc.) are
#'   only valid for 3D GMRF forms and will error if applied to the 2D AR1
#'   (\code{cont_tv_fish_sel == 5}).
#' @param bins Number of selectivity bins
#'
#' @return The input \code{input_list} with \code{$map$fishsel_pe_pars} set to a
#'   factor vector. Index numbering is reset after any correlation suppression to
#'   maintain contiguous integer indices.
#'
#' @keywords internal
do_fishsel_pe_pars_mapping <- function(input_list, fishsel_pe_pars_spec, corr_opt_semipar, bins) {

  # Initialize counter and mapping array for fishery process errors
  fishsel_pe_pars_counter <- 1 # initalize counter
  map_fishsel_pe_pars <- input_list$par$fishsel_pe_pars # initalize array
  map_fishsel_pe_pars[] <- NA

  # Fishery process error parameters
  for(f in 1:input_list$data$n_fish_fleets) {

    # Validate options
    if(!is.null(fishsel_pe_pars_spec)) {
      if(!fishsel_pe_pars_spec[f] %in% c("fix", "none", "est_all", "est_shared_r", "est_shared_s", "est_shared_r_s") &&
         !stringr::str_detect(fishsel_pe_pars_spec[f], "est_shared_f_\\d+"))
        stop("fishsel_pe_pars_spec not correctly specfied. Should be one of these: est_all, est_shared_r, est_shared_r_s, est_shared_s, fix, or est_shared_f_# (where # is fleet number)")
    }

    # Skip fleet sharing specs in first pass
    if(!is.null(fishsel_pe_pars_spec)) if(stringr::str_detect(fishsel_pe_pars_spec[f], "est_shared_f")) next

    for(r in 1:input_list$data$n_regions) {

      # if no time-variation, then fix all parameters for this fleet
      if(input_list$data$cont_tv_fish_sel[r,f] == 0 || (sum(input_list$data$UseCatch[r,,,f]) == 0 && sum(input_list$data$UseCatch_pop[,r,,,f]) == 0)) {
        map_fishsel_pe_pars[r,,,f] <- NA
      } else { # if we have time-variation

        # Figure out max number of selectivity parameters for a given region and fleet
        if(unique(input_list$data$fish_sel_model[r,,f]) %in% 2) max_sel_pars <- 1 # exponential
        if(unique(input_list$data$fish_sel_model[r,,f]) %in% c(0,1,3)) max_sel_pars <- 2 # logistic or gamma
        if(unique(input_list$data$fish_sel_model[r,,f]) == 4) max_sel_pars <- 6 # double normal
        if(unique(input_list$data$fish_sel_model[r,,f]) == 5) max_sel_pars <- bins # non-parametric selectivity
        if(unique(input_list$data$fish_sel_model[r,,f]) %in% c(6,7)) max_sel_pars <- 3 # logistic selectivity w/ asmyptote

        for(s in 1:input_list$data$n_sexes) {

          # If iid time-variation or random walk for this fleet
          if(input_list$data$cont_tv_fish_sel[r,f] %in% c(1,2)) {

            for(i in 1:max_sel_pars) {

              # either fixing parameters or not used for a given fleet
              if(fishsel_pe_pars_spec[f] %in% c("none", "fix")) map_fishsel_pe_pars[r,i,s,f] <- NA

              # Estimating all parameters separately (unique for each region, sex, fleet, parameter)
              if(fishsel_pe_pars_spec[f] == "est_all") {
                map_fishsel_pe_pars[r,i,s,f] <- fishsel_pe_pars_counter
                fishsel_pe_pars_counter <- fishsel_pe_pars_counter + 1
              } # end est_all

              # Estimating process error parameters shared across regions (but unique for each sex, fleet, parameter)
              if(fishsel_pe_pars_spec[f] == 'est_shared_r' && r == 1) {
                map_fishsel_pe_pars[,i,s,f] <- fishsel_pe_pars_counter
                fishsel_pe_pars_counter <- fishsel_pe_pars_counter + 1
              }

              # Estimating process error parameters shared across sexes (but unique for each region, fleet, parameter)
              if(fishsel_pe_pars_spec[f] == 'est_shared_s' && s == 1) {
                map_fishsel_pe_pars[r,i,,f] <- fishsel_pe_pars_counter
                fishsel_pe_pars_counter <- fishsel_pe_pars_counter + 1
              }

              # Estimating process error parameters shared across regions and sexes (but unique for each fleet, parameter)
              if(fishsel_pe_pars_spec[f] == 'est_shared_r_s' && r == 1 && s == 1) {
                map_fishsel_pe_pars[,i,,f] <- fishsel_pe_pars_counter
                fishsel_pe_pars_counter <- fishsel_pe_pars_counter + 1
              }

            } # end i loop
          } # end iid or random walk variation

          # If 3d gmrf or 2dar1
          if(input_list$data$cont_tv_fish_sel[r,f] %in% c(3,4,5)) {

            # Set up indexing to loop through
            if(input_list$data$cont_tv_fish_sel[r,f] %in% c(3,4)) idx = 1:4 # 3dgmrf (1 = pcorr_age, 2 = pcorr_year, 3= pcorr_cohort, 4 = log_sigma)
            if(input_list$data$cont_tv_fish_sel[r,f] %in% c(5)) idx = c(1,2,4) # 2dar1 (1 = pcorr_bin, 2 = pcorr_year, 4 = log_sigma)
            if(input_list$data$cont_tv_fish_sel[r,f] %in% c(3,4) && input_list$data$fish_selex_type == 1) stop("Cohort-based selectivity deviations are specified, but selectivity is specified as length-based. Please choose another deviation form!")

            for(i in idx) {

              # either fixing parameters or not used for a given fleet
              if(fishsel_pe_pars_spec[f] %in% c("none", "fix")) map_fishsel_pe_pars[r,i,s,f] <- NA

              # Estimating all process error parameters
              if(fishsel_pe_pars_spec[f] == "est_all") {
                map_fishsel_pe_pars[r,i,s,f] <- fishsel_pe_pars_counter
                fishsel_pe_pars_counter <- fishsel_pe_pars_counter + 1
              } # end est_all

              # Estimating process error parameters shared across regions (but unique for each sex, fleet, parameter)
              if(fishsel_pe_pars_spec[f] == 'est_shared_r' && r == 1) {
                map_fishsel_pe_pars[,i,s,f] <- fishsel_pe_pars_counter
                fishsel_pe_pars_counter <- fishsel_pe_pars_counter + 1
              }

              # Estimating process error parameters shared across sexes (but unique for each region, fleet, parameter)
              if(fishsel_pe_pars_spec[f] == 'est_shared_s' && s == 1) {
                map_fishsel_pe_pars[r,i,,f] <- fishsel_pe_pars_counter
                fishsel_pe_pars_counter <- fishsel_pe_pars_counter + 1
              }

              # Estimating process error parameters shared across regions and sexes (but unique for each fleet, parameter)
              if(fishsel_pe_pars_spec[f] == 'est_shared_r_s' && r == 1 && s == 1) {
                map_fishsel_pe_pars[,i,,f] <- fishsel_pe_pars_counter
                fishsel_pe_pars_counter <- fishsel_pe_pars_counter + 1
              }

            } # end i loop

            # Options to set correaltions to 0 for 3dgmrf
            if(!is.null(corr_opt_semipar)) {

              opt <- input_list$data$cont_tv_fish_sel[r,f] # get random effects options

              # Validate options
              if(!corr_opt_semipar[f] %in% c(NA, "corr_zero_y", "corr_zero_b", "corr_zero_y_b", "corr_zero_c", "corr_zero_y_c", "corr_zero_b_c", "corr_zero_y_b_c"))
                stop("corr_opt_semipar not correctly specfied. Should be one of these: corr_zero_y, corr_zero_b, corr_zero_y_b, corr_zero_c, corr_zero_y_c, corr_zero_b_c, corr_zero_y_b_c, NA")
              if(opt == 5 && corr_opt_semipar[f] %in% c("corr_zero_c","corr_zero_y_c","corr_zero_b_c","corr_zero_y_b_c"))
                stop("Invalid corr_opt_semipar for 2dar1 (opt=5): cohort correlations are not allowed.")

              if (opt %in% c(3,4,5)) {
                # 2d and 3d options
                if (corr_opt_semipar[f] == "corr_zero_y")    map_fishsel_pe_pars[,2,,f]     <- NA
                if (corr_opt_semipar[f] == "corr_zero_b")    map_fishsel_pe_pars[,1,,f]     <- NA
                if (corr_opt_semipar[f] == "corr_zero_y_b")  map_fishsel_pe_pars[,1:2,,f]   <- NA
              }

              if(opt %in% c(3,4)) {
                # 3d gmrf options only (adds the cohort dimension)
                if (corr_opt_semipar[f] == "corr_zero_c")      map_fishsel_pe_pars[,3,,f]   <- NA
                if (corr_opt_semipar[f] == "corr_zero_y_c")    map_fishsel_pe_pars[,2:3,,f] <- NA
                if (corr_opt_semipar[f] == "corr_zero_b_c")    map_fishsel_pe_pars[,c(1,3),,f] <- NA
                if (corr_opt_semipar[f] == "corr_zero_y_b_c")  map_fishsel_pe_pars[,1:3,,f] <- NA
              }

              # Reset numbering for mapping off correlation parameters for clarity
              non_na_positions <- which(!is.na(map_fishsel_pe_pars))
              map_fishsel_pe_pars[non_na_positions] <- seq_along(non_na_positions)
              collect_message("corr_opt_semipar is specified as: ", corr_opt_semipar[f], "for fishery fleet", f)

            }
          } # end if 3d gmrf marginal or conditional variance

          # fix all parameters
          if(fishsel_pe_pars_spec[f] == "fix") map_fishsel_pe_pars[r,,s,f] <- NA

        } # end s loop
      } # end else
    } # end r loop

    if(!is.null(fishsel_pe_pars_spec)) collect_message("fishsel_pe_pars_spec is specified as: ", fishsel_pe_pars_spec[f], "for fishery fleet", f)

  } # end f loop

  # Handle fleet sharing after all base mappings are established
  for(f in 1:input_list$data$n_fish_fleets) {
    if(stringr::str_detect(fishsel_pe_pars_spec[f], "est_shared_f") && !is.null(fishsel_pe_pars_spec)) {
      # extract fleet sharing index
      flt_shared <- as.numeric(unlist(stringr::str_extract_all(fishsel_pe_pars_spec[f], "\\d+")))

      # Validate options here
      if(flt_shared > input_list$data$n_fish_fleets || flt_shared < 1) stop("Fleet sharing specification 'est_shared_f", flt_shared, "' for fleet ", f, " references invalid fleet number.")
      if(stringr::str_detect(fishsel_pe_pars_spec[flt_shared], "est_shared_f")) stop("Fleet ", f, " cannot share with fleet ", flt_shared, " because fleet ", flt_shared, " is self-sharing parameters, which does not make sense.")

      # Copy mapping from reference fleet
      map_fishsel_pe_pars[,,,f] <- map_fishsel_pe_pars[,,,flt_shared]
      collect_message("fishsel_pe_pars_spec is specified as: ", fishsel_pe_pars_spec[f], " for fishery fleet ", f, " (sharing with fleet ", flt_shared, ")")
    } # end if statement
  } # end f loop

  # input into mapping list
  input_list$map$fishsel_pe_pars <- factor(map_fishsel_pe_pars)

  return(input_list)
}

#' Map fishery selectivity deviation parameters
#'
#' Constructs the factor map for \code{ln_fishsel_devs}, the annual deviations in
#' continuous time-varying fishery selectivity. For iid and random-walk forms,
#' the deviation dimension corresponds to selectivity parameters (up to 6 for
#' double-normal); for semi-parametric forms (3D GMRF, 2D AR1), it corresponds
#' to age or length bins. Cells with no time-variation
#' (\code{cont_tv_fish_sel == 0}) or no catch data are mapped to \code{NA}.
#'
#' Bin-sharing (\code{"est_shared_b"} and related options) groups bins into
#' blocks defined by \code{fishsel_devs_shared_bins}, reducing the number of
#' estimated deviation series. Fleet sharing (\code{"est_shared_f_x"}) is handled
#' in a second pass. The resulting integer map is also stored as
#' \code{$data$map_ln_fishsel_devs} for use in the RTMB objective function.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists.
#' @param fish_sel_devs_spec Character vector of length \code{n_fish_fleets}.
#'   Options:
#'   \describe{
#'     \item{\code{"est_all"}}{Separate deviation series per region × sex × bin.}
#'     \item{\code{"est_shared_r"}}{Shared across regions.}
#'     \item{\code{"est_shared_s"}}{Shared across sexes.}
#'     \item{\code{"est_shared_r_s"}}{Shared across regions and sexes.}
#'     \item{\code{"est_shared_b"}}{Shared across bin groups defined by
#'       \code{fishsel_devs_shared_bins}.}
#'     \item{\code{"est_shared_r_b"}}{Shared across regions and bin groups.}
#'     \item{\code{"est_shared_b_s"}}{Shared across bin groups and sexes.}
#'     \item{\code{"est_shared_r_b_s"}}{Shared across regions, bin groups, and sexes.}
#'     \item{\code{"est_shared_f_x"}}{Copy deviation map from fleet \code{x}.}
#'     \item{\code{"fix"} or \code{"none"}}{All deviations fixed at zero (mapped to \code{NA}).}
#'   }
#'   Bin-sharing options (\code{"est_shared_b"}, etc.) are only valid with
#'   semi-parametric time-varying forms (\code{cont_tv_fish_sel} is in \code{c(3,4,5)}).
#' @param fishsel_devs_shared_bins List of integer vectors defining bin
#'   groupings for bin-sharing options. Each element groups bins that share a
#'   single deviation series, e.g., \code{list(1:5, 6:10, 11:30)}.
#'   Only used when \code{fish_sel_devs_spec} includes \code{"est_shared_b"}.
#' @param bins Number of selectivity bins
#'
#' @return The input \code{input_list} with \code{$map$ln_fishsel_devs} set to a
#'   factor vector and \code{$data$map_ln_fishsel_devs} set to the corresponding
#'   integer array (for use in the RTMB template).
#'
#' @keywords internal
do_fishsel_devs_mapping <- function(input_list, fish_sel_devs_spec, fishsel_devs_shared_bins, bins) {

  # Initialize counter and mapping array for fishery selectivity deviations
  fishsel_devs_counter <- 1
  map_fishsel_devs <- input_list$par$ln_fishsel_devs
  map_fishsel_devs[] <- NA

  for(r in 1:input_list$data$n_regions) {
    for(f in 1:input_list$data$n_fish_fleets) {

      # Validate options
      if(!is.null(fish_sel_devs_spec)) {
        if(!fish_sel_devs_spec[f] %in% c("fix", "none", "est_all", "est_shared_r", "est_shared_s", "est_shared_r_s", "est_shared_b", "est_shared_r_b", "est_shared_r_b_s", "est_shared_b_s") &&
           !stringr::str_detect(fish_sel_devs_spec[f], "est_shared_f_\\d+"))
          stop("fish_sel_devs_spec not correctly specfied. Should be one of these: est_all, est_shared_r, est_shared_r_s, est_shared_s, est_shared_b, est_shared_r_b, est_shared_r_b_s, est_shared_r_s, fix, or est_shared_f_# (where # is fleet number)")
        if(fish_sel_devs_spec[f] %in% c("est_shared_b", "est_shared_r_b", "est_shared_r_b_s", "est_shared_b_s") &&
           !input_list$data$cont_tv_fish_sel[r,f] %in% c(3,4,5)) stop("Sharing age deviations with iid or random walk parametric forms is not supported!")
       }

      # Skip fleet sharing specs in first pass
      if(!is.null(fish_sel_devs_spec)) if(stringr::str_detect(fish_sel_devs_spec[f], "est_shared_f")) next

      for(s in 1:input_list$data$n_sexes) {
        for(y in 1:(length(input_list$data$years) + input_list$data$n_proj_yrs_devs)) {

          # if no time-variation, then fix all parameters for this fleet
          if(input_list$data$cont_tv_fish_sel[r,f] == 0 || (sum(input_list$data$UseCatch[r,,,f]) == 0 && sum(input_list$data$UseCatch_pop[,r,,,f]) == 0) ) {
            map_fishsel_devs[r,y,,s,f] <- NA
          } else {

            # Figure out max number of selectivity parameters for a given region and fleet
            if(unique(input_list$data$fish_sel_model[r,,f]) %in% 2) max_sel_pars <- 1 # exponential
            if(unique(input_list$data$fish_sel_model[r,,f]) %in% c(0,1,3)) max_sel_pars <- 2 # logistic or gamma
            if(unique(input_list$data$fish_sel_model[r,,f]) == 4) max_sel_pars <- 6 # double normal
            if(unique(input_list$data$fish_sel_model[r,,f]) == 5) max_sel_pars <- bins # non-parametric selectivity
            if(unique(input_list$data$fish_sel_model[r,,f]) %in% c(6,7)) max_sel_pars <- 3 # logistic selectivity w/ asmyptote

            # If iid or random walk time-variation for this fleet
            if(input_list$data$cont_tv_fish_sel[r,f] %in% c(1,2)) {

              for(i in 1:max_sel_pars) {
                # Estimating all selectivity deviations across regions, sexes, fleets, and parameter
                if(fish_sel_devs_spec[f] == 'est_all') {
                  map_fishsel_devs[r,y,i,s,f] <- fishsel_devs_counter
                  fishsel_devs_counter <- fishsel_devs_counter + 1
                }

                # Estimating selectivity deviations across sexes, fleets, and parameters, but shared across regions
                if(fish_sel_devs_spec[f] == 'est_shared_r' && r == 1) {
                  map_fishsel_devs[,y,i,s,f] <- fishsel_devs_counter
                  fishsel_devs_counter <- fishsel_devs_counter + 1
                }

                # Estimating selectivity deviations across regions, fleets, and parameters, but shared across sexes
                if(fish_sel_devs_spec[f] == 'est_shared_s' && s == 1) {
                  map_fishsel_devs[r,y,i,,f] <- fishsel_devs_counter
                  fishsel_devs_counter <- fishsel_devs_counter + 1
                }

                # Estimating selectivity deviations across fleets, and parameters, but shared across sexes and regions
                if(fish_sel_devs_spec[f] == 'est_shared_r_s' && r == 1 && s == 1) {
                  map_fishsel_devs[,y,i,,f] <- fishsel_devs_counter
                  fishsel_devs_counter <- fishsel_devs_counter + 1
                }

              } # end i loop
            } # end iid or random walk variation

            # If 3d gmrf for this fleet
            if(input_list$data$cont_tv_fish_sel[r,f] %in% c(3,4,5)) {

              for(i in 1:length(input_list$data$ages)) {
                # Estimating all selectivity deviations across regions, years and bins
                if(fish_sel_devs_spec[f] == 'est_all') {
                  map_fishsel_devs[r,y,i,s,f] <- fishsel_devs_counter
                  fishsel_devs_counter <- fishsel_devs_counter + 1
                }

                # Estimating all selectivity deviations across years and bins, but shared across regions
                if(fish_sel_devs_spec[f] == 'est_shared_r' && r == 1) {
                  map_fishsel_devs[,y,i,s,f] <- fishsel_devs_counter
                  fishsel_devs_counter <- fishsel_devs_counter + 1
                }

                # Estimating all selectivity deviations across years and bins, but shared across sexes
                if(fish_sel_devs_spec[f] == 'est_shared_s' && s == 1) {
                  map_fishsel_devs[r,y,i,,f] <- fishsel_devs_counter
                  fishsel_devs_counter <- fishsel_devs_counter + 1
                }

                # Estimating all selectivity deviations across years and bins, but shared across sexes and regions
                if(fish_sel_devs_spec[f] == 'est_shared_r_s' && s == 1 && r == 1) {
                  map_fishsel_devs[,y,i,,f] <- fishsel_devs_counter
                  fishsel_devs_counter <- fishsel_devs_counter + 1
                }

                if(fish_sel_devs_spec[f] == 'est_shared_b') {
                  for(k in 1:length(fishsel_devs_shared_bins)) {
                    map_fishsel_devs[r,y,fishsel_devs_shared_bins[[k]],s,f] <- fishsel_devs_counter
                    fishsel_devs_counter <- fishsel_devs_counter + 1
                  } # end k loop
                }

                if(fish_sel_devs_spec[f] == 'est_shared_r_b' && r == 1) {
                  for(k in 1:length(fishsel_devs_shared_bins)) {
                    map_fishsel_devs[,y,fishsel_devs_shared_bins[[k]],s,f] <- fishsel_devs_counter
                    fishsel_devs_counter <- fishsel_devs_counter + 1
                  } # end k loop
                }

                if(fish_sel_devs_spec[f] == 'est_shared_b_s' && s == 1) {
                  for(k in 1:length(fishsel_devs_shared_bins)) {
                    map_fishsel_devs[r,y,fishsel_devs_shared_bins[[k]],,f] <- fishsel_devs_counter
                    fishsel_devs_counter <- fishsel_devs_counter + 1
                  } # end k loop
                }

                if(fish_sel_devs_spec[f] == 'est_shared_r_b_s' && s == 1 && r == 1) {
                  for(k in 1:length(fishsel_devs_shared_bins)) {
                    map_fishsel_devs[,y,fishsel_devs_shared_bins[[k]],,f] <- fishsel_devs_counter
                    fishsel_devs_counter <- fishsel_devs_counter + 1
                  } # end k loop
                }

              } # end i loop
            } # end 3d gmrf

          } # end else
        } # end y loop
      } # end s loop

      if(!is.null(fish_sel_devs_spec)) collect_message("fish_sel_devs_spec is specified as: ", fish_sel_devs_spec[f], "for fishery fleet", f, "and region ", r)

    } # end f loop
  } # end r loop

  # Handle fleet sharing after all base mappings are established
  for(f in 1:input_list$data$n_fish_fleets) {
    if(stringr::str_detect(fish_sel_devs_spec[f], "est_shared_f") && !is.null(fish_sel_devs_spec)) {
      # extract fleet sharing index
      flt_shared <- as.numeric(unlist(stringr::str_extract_all(fish_sel_devs_spec[f], "\\d+")))

      # Validate options here
      if(flt_shared > input_list$data$n_fish_fleets || flt_shared < 1) stop("Fleet sharing specification 'est_shared_f", flt_shared, "' for fleet ", f, " references invalid fleet number.")
      if(stringr::str_detect(fish_sel_devs_spec[flt_shared], "est_shared_f")) stop("Fleet ", f, " cannot share with fleet ", flt_shared, " because fleet ", flt_shared, " is self-sharing parameters, which does not make sense.")

      # Copy mapping from reference fleet
      map_fishsel_devs[,,,,f] <- map_fishsel_devs[,,,,flt_shared]
      collect_message("fish_sel_devs_spec is specified as: ", fish_sel_devs_spec[f], " for fishery fleet ", f, " (sharing with fleet ", flt_shared, ")")
    } # end if statement
  } # end f loop

  # input into mapping list
  input_list$map$ln_fishsel_devs <- factor(map_fishsel_devs)
  input_list$data$map_ln_fishsel_devs <- array(as.numeric(input_list$map$ln_fishsel_devs), dim = dim(input_list$par$ln_fishsel_devs))

  return(input_list)
}

#' Map retained fishery selectivity fixed-effect parameters
#'
#' Constructs the factor map for \code{ret_fixed_sel_pars} (e.g., \eqn{a_{50}}, \eqn{k}, \eqn{a_{max}}),
#' controlling whether selectivity shape parameters are estimated independently or shared across
#' regions, sexes, or fleets. Cells with no catch data (\code{UseCatch == 0} and
#' \code{UseCatch_pop == 0}) are automatically mapped to \code{NA}.
#'
#' Non-parametric selectivity (model \code{sel_model == 5}) is supported via
#' \code{ret_sel_nonpar_est_bins}, which defines bin groupings that are treated
#' as individual estimated parameters.
#'
#' Fleet sharing (\code{"est_shared_f_x"}) is applied in a second pass after all
#' base mappings are constructed, copying the full parameter index structure from
#' the reference fleet.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists.
#'
#' @param ret_fixed_sel_pars_spec Character vector of length \code{n_fish_fleets}.
#' Each element specifies the estimation structure for one fleet. Options:
#' \describe{
#'   \item{\code{"est_all"}}{Separate parameters per region × sex × block × bin.}
#'   \item{\code{"est_shared_r"}}{Shared across regions; unique per sex × block × bin.}
#'   \item{\code{"est_shared_s"}}{Shared across sexes; unique per region × block × bin.}
#'   \item{\code{"est_shared_r_s"}}{Shared across regions and sexes; unique per block × bin.}
#'   \item{\code{"est_shared_f_x"}}{Copy parameter mapping from fleet \code{x}. Fleet \code{x}
#'     must not itself use \code{"est_shared_f_y"}.}
#'   \item{\code{"fix"} / \code{"fixed_ret_sel_input"}}{All parameters fixed (mapped to \code{NA}).}
#' }
#'
#' @param bins Number of selectivity bins.
#'
#' @param ret_sel_nonpar_est_bins Optional list specifying bin groupings for
#' non-parametric selectivity. Structure is \code{[[fleet]][[block]]}, where each
#' element is a list of bin index vectors defining grouped parameters.
#'
#' @return Updated \code{input_list} with \code{$map$ret_fixed_sel_pars}
#'   as a factor array.
#'
#' @keywords internal
do_ret_fixed_sel_pars_mapping <- function(input_list, ret_fixed_sel_pars_spec, bins, ret_sel_nonpar_est_bins) {

  # Initialize counter and mapping array for fixed effects retained fishery selectivity
  ret_fixed_sel_pars_counter <- 1
  map_ret_fixed_sel_pars <- input_list$par$ret_fixed_sel_pars
  map_ret_fixed_sel_pars[] <- NA

  for(f in 1:input_list$data$n_fish_fleets) {

    # Validate Options
    if(!ret_fixed_sel_pars_spec[f] %in% c("est_all", "est_shared_r", "est_shared_r_s", "fix", "est_shared_s", 'fix_ret_sel_input') &&
       !stringr::str_detect(ret_fixed_sel_pars_spec[f], "est_shared_f_\\d+"))
      stop("ret_fixed_sel_pars_spec not correctly specfied. Should be one of these: est_all, est_shared_r, est_shared_r_s, est_shared_s, fix, or est_shared_f_# (where # is fleet number)")

    # check fixed selectivity options
    if(input_list$data$use_fixed_ret_sel[f] == 1 && stringr::str_detect(ret_fixed_sel_pars_spec[f], 'est'))
      stop("use_fixed_ret_sel has 1s for a given fleet, but ret_fixed_sel_pars_spec is specified at an est variant.")
    if(input_list$data$use_fixed_ret_sel[f] == 0 && ret_fixed_sel_pars_spec[f] == 'fix_ret_sel_input')
      stop("use_fixed_ret_sel has 0s for a given fleet, but ret_fixed_sel_pars_spec is specified at fix_ret_sel_input.")

    # Skip fleet sharing specs in first pass
    if(stringr::str_detect(ret_fixed_sel_pars_spec[f], "est_shared_f")) next

    for(r in 1:input_list$data$n_regions) {

      # Only add a counter if caatches are avaliable in some years for a given region and fleet combination
      if(sum(input_list$data$UseCatch[r,,,f]) > 0 || sum(input_list$data$UseCatch_pop[,r,,,f]) > 0) {

        # Extract number of retained fishery selectivity blocks
        retsel_blocks_tmp <- unique(as.vector(input_list$data$ret_sel_blocks[r,,f]))

        for(s in 1:input_list$data$n_sexes) {
          for(b in 1:length(retsel_blocks_tmp)) {

            block_years <- which(input_list$data$ret_sel_blocks[r,,f] == retsel_blocks_tmp[b]) # figure out block years
            sel_model_this_block <- unique(input_list$data$ret_sel_model[r, block_years, f]) # get selectivity form for a given block
            if(length(sel_model_this_block) > 1) stop("Block ", retsel_blocks_tmp[b], " for fleet ", f, " region ", r, " has multiple selectivity models assigned to it")

            # determine maximum selectivity parameters
            if(sel_model_this_block == 2) max_sel_pars <- 1 # exponential
            if(sel_model_this_block %in% c(0,1,3)) max_sel_pars <- 2 # logistic or gamma
            if(sel_model_this_block == 4) max_sel_pars <- 6 # double normal
            if(sel_model_this_block %in% c(6,7)) max_sel_pars <- 3 # logistic with an asymptote
            if(sel_model_this_block == 8) { # bicubic spline: flattened bin-node x year-node grid
              n_bin_nodes_this <- unique(input_list$data$ret_sel_bicubic_binnodes[r, block_years, f])
              n_yr_nodes_this <- unique(input_list$data$ret_sel_bicubic_yrnodes[r, block_years, f])
              max_sel_pars <- n_bin_nodes_this * n_yr_nodes_this
            }

            # non-parametric selectivity
            if(sel_model_this_block == 5) {

              if(is.null(ret_sel_nonpar_est_bins)) stop("Non-parametric retention selectivtiy specified, but ret_sel_nonpar_est_bins is NULL. Please specify bins!")

              bin_groups <- ret_sel_nonpar_est_bins[[f]][[b]]
              max_sel_pars <- length(bin_groups)  # number of groups = number of estimated pars

              # validate
              all_bins <- unlist(bin_groups)
              if(any(all_bins < 1) || any(all_bins > bins))
                stop("ret_sel_nonpar_est_bins[[", f, "]][[", b, "]] contains indices outside 1:", bins)
              if(length(all_bins) != length(unique(all_bins)))
                stop("ret_sel_nonpar_est_bins[[", f, "]][[", b, "]] has duplicate bin indices")
            }

            for(i in 1:max_sel_pars) {

              # get non-parametric selectivity bins
              group_bins <- if(sel_model_this_block == 5) bin_groups[[i]] else i

              # Estimate all selectivity fixed effects parameters within the constraints of the defined blocks
              if(ret_fixed_sel_pars_spec[f] == "est_all") {
                for(bi in group_bins) map_ret_fixed_sel_pars[r,bi,b,s,f] <- ret_fixed_sel_pars_counter
                ret_fixed_sel_pars_counter <- ret_fixed_sel_pars_counter + 1
              } # end if

              # Estimating process error parameters shared across regions (but unique for each sex, fleet, parameter)
              if(ret_fixed_sel_pars_spec[f] == 'est_shared_r' && r == 1) {
                for(rr in 1:input_list$data$n_regions) {
                  if(retsel_blocks_tmp[b] %in% input_list$data$ret_sel_blocks[rr,,f]) {
                    for(bi in group_bins) map_ret_fixed_sel_pars[rr, bi, b, s, f] <- ret_fixed_sel_pars_counter
                  } # end if
                } # end rr loop
                ret_fixed_sel_pars_counter <- ret_fixed_sel_pars_counter + 1
              } # end if

              # Estimating process error parameters shared across sexes (but unique for each region, fleet, parameter)
              if(ret_fixed_sel_pars_spec[f] == 'est_shared_s' && s == 1) {
                for(ss in 1:input_list$data$n_sexes) {
                  for(bi in group_bins) map_ret_fixed_sel_pars[r, bi, b, ss, f] <- ret_fixed_sel_pars_counter
                } # end ss loop
                ret_fixed_sel_pars_counter <- ret_fixed_sel_pars_counter + 1
              } # end if

              # Estimating process error parameters shared across regions and sexes (but unique for each fleet, parameter)
              if(ret_fixed_sel_pars_spec[f] == 'est_shared_r_s' && r == 1 && s == 1) {
                for(rr in 1:input_list$data$n_regions) {
                  for(ss in 1:input_list$data$n_sexes) {
                    if(retsel_blocks_tmp[b] %in% input_list$data$ret_sel_blocks[rr,,f]) {
                      for(bi in group_bins) map_ret_fixed_sel_pars[rr, bi, b, ss, f] <- ret_fixed_sel_pars_counter
                    } # end if
                  } # end ss loop
                } #end rr loop
                ret_fixed_sel_pars_counter <- ret_fixed_sel_pars_counter + 1
              } # end if

            } # end i loop
          } # end b loop
        } # end s loop
      } # end if statement
    } # end r loop

    # fix all parameters
    if(ret_fixed_sel_pars_spec[f] %in% c("fix", "fixed_ret_sel_input")) map_ret_fixed_sel_pars[,,,,f] <- NA
    collect_message("ret_fixed_sel_pars_spec is specified as: ", ret_fixed_sel_pars_spec[f], " for fishery fleet ", f)

  } # end f loop

  # Handle fleet sharing after all base mappings are established
  for(f in 1:input_list$data$n_fish_fleets) {
    if(stringr::str_detect(ret_fixed_sel_pars_spec[f], "est_shared_f")) {
      # extract fleet sharing index
      flt_shared <- as.numeric(unlist(stringr::str_extract_all(ret_fixed_sel_pars_spec[f], "\\d+")))

      # Validate options here
      if(flt_shared > input_list$data$n_fish_fleets || flt_shared < 1) stop("Fleet sharing specification 'est_shared_f", flt_shared, "' for fleet ", f, " references invalid fleet number.")
      if(stringr::str_detect(ret_fixed_sel_pars_spec[flt_shared], "est_shared_f")) stop("Fleet ", f, " cannot share with fleet ", flt_shared, " because fleet ", flt_shared, " is self-sharing parameters, which does not make sense.")

      # Copy mapping from reference fleet
      map_ret_fixed_sel_pars[,,,,f] <- map_ret_fixed_sel_pars[,,,,flt_shared]
      collect_message("ret_fixed_sel_pars_spec is specified as: ", ret_fixed_sel_pars_spec[f], " for fishery fleet ", f, " (sharing with fleet ", flt_shared, ")")
    } # end if statement
  } # end f loop

  # input into mapping list
  input_list$map$ret_fixed_sel_pars <- factor(map_ret_fixed_sel_pars)
  return(input_list)
}

#' Map retained fishery selectivity process error hyperparameters
#'
#' Constructs the factor map for \code{retsel_pe_pars}, which contains the
#' variance and correlation hyperparameters governing continuous time-varying
#' selectivity. The set of active parameters depends on the time-variation type
#' (\code{cont_tv_ret_sel}): iid/random-walk forms use up to 2 parameters
#' (log-sigma); 3D GMRF forms use up to 4 (partial correlations for age, year,
#' cohort dimensions plus log-sigma); the 2D AR1 form uses 3 (bin AR1, year AR1,
#' log-sigma). Correlation components can be selectively suppressed via
#' \code{ret_sel_corr_opt_semipar}.
#'
#' Fleet sharing (\code{"est_shared_f_x"}) is handled in a second pass.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists.
#' @param retsel_pe_pars_spec Character vector of length \code{n_fish_fleets}.
#'   Options:
#'   \describe{
#'     \item{\code{"est_all"}}{Separate hyperparameters per region × sex.}
#'     \item{\code{"est_shared_r"}}{Shared across regions; unique per sex.}
#'     \item{\code{"est_shared_s"}}{Shared across sexes; unique per region.}
#'     \item{\code{"est_shared_r_s"}}{Shared across regions and sexes.}
#'     \item{\code{"est_shared_f_x"}}{Copy hyperparameters from fleet \code{x}.}
#'     \item{\code{"fix"} or \code{"none"}}{All parameters fixed (mapped to \code{NA}).}
#'   }
#' @param ret_sel_corr_opt_semipar Character vector of length \code{n_fish_fleets}
#'   specifying which correlation components to suppress for semi-parametric
#'   models. Valid values per fleet: \code{NA} (no suppression),
#'   \code{"corr_zero_y"}, \code{"corr_zero_b"}, \code{"corr_zero_y_b"},
#'   \code{"corr_zero_c"}, \code{"corr_zero_y_c"}, \code{"corr_zero_b_c"},
#'   \code{"corr_zero_y_b_c"}. Cohort options (\code{"corr_zero_c"}, etc.) are
#'   only valid for 3D GMRF forms and will error if applied to the 2D AR1
#'   (\code{cont_tv_ret_sel == 5}).
#'
#' @return The input \code{input_list} with \code{$map$retsel_pe_pars} set to a
#'   factor vector. Index numbering is reset after any correlation suppression to
#'   maintain contiguous integer indices.
#'
#' @keywords internal
do_retsel_pe_pars_mapping <- function(input_list, retsel_pe_pars_spec, ret_sel_corr_opt_semipar, bins) {

  # Initialize counter and mapping array for fishery process errors
  retsel_pe_pars_counter <- 1 # initalize counter
  map_retsel_pe_pars <- input_list$par$retsel_pe_pars # initalize array
  map_retsel_pe_pars[] <- NA

  # Fishery process error parameters
  for(f in 1:input_list$data$n_fish_fleets) {

    # Validate options
    if(!is.null(retsel_pe_pars_spec)) {
      if(!retsel_pe_pars_spec[f] %in% c("fix", "none", "est_all", "est_shared_r", "est_shared_s", "est_shared_r_s") &&
         !stringr::str_detect(retsel_pe_pars_spec[f], "est_shared_f_\\d+"))
        stop("retsel_pe_pars_spec not correctly specfied. Should be one of these: est_all, est_shared_r, est_shared_r_s, est_shared_s, fix, or est_shared_f_# (where # is fleet number)")
    }

    # Skip fleet sharing specs in first pass
    if(!is.null(retsel_pe_pars_spec)) if(stringr::str_detect(retsel_pe_pars_spec[f], "est_shared_f")) next

    for(r in 1:input_list$data$n_regions) {

      # if no time-variation, then fix all parameters for this fleet
      if(input_list$data$cont_tv_ret_sel[r,f] == 0 || (sum(input_list$data$UseCatch[r,,,f]) == 0 && sum(input_list$data$UseCatch_pop[,r,,,f]) == 0)) {
        map_retsel_pe_pars[r,,,f] <- NA
      } else { # if we have time-variation

        # Figure out max number of selectivity parameters for a given region and fleet
        if(unique(input_list$data$ret_sel_model[r,,f]) %in% 2) max_sel_pars <- 1 # exponential
        if(unique(input_list$data$ret_sel_model[r,,f]) %in% c(0,1,3)) max_sel_pars <- 2 # logistic or gamma
        if(unique(input_list$data$ret_sel_model[r,,f]) == 4) max_sel_pars <- 6 # double normal
        if(unique(input_list$data$ret_sel_model[r,,f]) == 5) max_sel_pars <- bins # non-parametric selectivity
        if(unique(input_list$data$ret_sel_model[r,,f]) %in% c(6,7)) max_sel_pars <- 3 # logistic selectivity with asymptote

        for(s in 1:input_list$data$n_sexes) {

          # If iid time-variation or random walk for this fleet
          if(input_list$data$cont_tv_ret_sel[r,f] %in% c(1,2)) {

            for(i in 1:max_sel_pars) {

              # either fixing parameters or not used for a given fleet
              if(retsel_pe_pars_spec[f] %in% c("none", "fix")) map_retsel_pe_pars[r,i,s,f] <- NA

              # Estimating all parameters separately (unique for each region, sex, fleet, parameter)
              if(retsel_pe_pars_spec[f] == "est_all") {
                map_retsel_pe_pars[r,i,s,f] <- retsel_pe_pars_counter
                retsel_pe_pars_counter <- retsel_pe_pars_counter + 1
              } # end est_all

              # Estimating process error parameters shared across regions (but unique for each sex, fleet, parameter)
              if(retsel_pe_pars_spec[f] == 'est_shared_r' && r == 1) {
                map_retsel_pe_pars[,i,s,f] <- retsel_pe_pars_counter
                retsel_pe_pars_counter <- retsel_pe_pars_counter + 1
              }

              # Estimating process error parameters shared across sexes (but unique for each region, fleet, parameter)
              if(retsel_pe_pars_spec[f] == 'est_shared_s' && s == 1) {
                map_retsel_pe_pars[r,i,,f] <- retsel_pe_pars_counter
                retsel_pe_pars_counter <- retsel_pe_pars_counter + 1
              }

              # Estimating process error parameters shared across regions and sexes (but unique for each fleet, parameter)
              if(retsel_pe_pars_spec[f] == 'est_shared_r_s' && r == 1 && s == 1) {
                map_retsel_pe_pars[,i,,f] <- retsel_pe_pars_counter
                retsel_pe_pars_counter <- retsel_pe_pars_counter + 1
              }

            } # end i loop
          } # end iid or random walk variation

          # If 3d gmrf or 2dar1
          if(input_list$data$cont_tv_ret_sel[r,f] %in% c(3,4,5)) {

            # Set up indexing to loop through
            if(input_list$data$cont_tv_ret_sel[r,f] %in% c(3,4)) idx = 1:4 # 3dgmrf (1 = pcorr_age, 2 = pcorr_year, 3= pcorr_cohort, 4 = log_sigma)
            if(input_list$data$cont_tv_ret_sel[r,f] %in% c(5)) idx = c(1,2,4) # 2dar1 (1 = pcorr_bin, 2 = pcorr_year, 4 = log_sigma)
            if(input_list$data$cont_tv_ret_sel[r,f] %in% c(3,4) && input_list$data$ret_selex_type == 1) stop("Cohort-based selectivity deviations are specified, but selectivity is specified as length-based. Please choose another deviation form!")

            for(i in idx) {

              # either fixing parameters or not used for a given fleet
              if(retsel_pe_pars_spec[f] %in% c("none", "fix")) map_retsel_pe_pars[r,i,s,f] <- NA

              # Estimating all process error parameters
              if(retsel_pe_pars_spec[f] == "est_all") {
                map_retsel_pe_pars[r,i,s,f] <- retsel_pe_pars_counter
                retsel_pe_pars_counter <- retsel_pe_pars_counter + 1
              } # end est_all

              # Estimating process error parameters shared across regions (but unique for each sex, fleet, parameter)
              if(retsel_pe_pars_spec[f] == 'est_shared_r' && r == 1) {
                map_retsel_pe_pars[,i,s,f] <- retsel_pe_pars_counter
                retsel_pe_pars_counter <- retsel_pe_pars_counter + 1
              }

              # Estimating process error parameters shared across sexes (but unique for each region, fleet, parameter)
              if(retsel_pe_pars_spec[f] == 'est_shared_s' && s == 1) {
                map_retsel_pe_pars[r,i,,f] <- retsel_pe_pars_counter
                retsel_pe_pars_counter <- retsel_pe_pars_counter + 1
              }

              # Estimating process error parameters shared across regions and sexes (but unique for each fleet, parameter)
              if(retsel_pe_pars_spec[f] == 'est_shared_r_s' && r == 1 && s == 1) {
                map_retsel_pe_pars[,i,,f] <- retsel_pe_pars_counter
                retsel_pe_pars_counter <- retsel_pe_pars_counter + 1
              }

            } # end i loop

            # Options to set correaltions to 0 for 3dgmrf
            if(!is.null(ret_sel_corr_opt_semipar)) {

              opt <- input_list$data$cont_tv_ret_sel[r,f] # get random effects options

              # Validate options
              if(!ret_sel_corr_opt_semipar[f] %in% c(NA, "corr_zero_y", "corr_zero_b", "corr_zero_y_b", "corr_zero_c", "corr_zero_y_c", "corr_zero_b_c", "corr_zero_y_b_c"))
                stop("ret_sel_corr_opt_semipar not correctly specfied. Should be one of these: corr_zero_y, corr_zero_b, corr_zero_y_b, corr_zero_c, corr_zero_y_c, corr_zero_b_c, corr_zero_y_b_c, NA")
              if(opt == 5 && ret_sel_corr_opt_semipar[f] %in% c("corr_zero_c","corr_zero_y_c","corr_zero_b_c","corr_zero_y_b_c"))
                stop("Invalid ret_sel_corr_opt_semipar for 2dar1 (opt=5): cohort correlations are not allowed.")

              if (opt %in% c(3,4,5)) {
                # 2d and 3d options
                if (ret_sel_corr_opt_semipar[f] == "corr_zero_y")    map_retsel_pe_pars[,2,,f]     <- NA
                if (ret_sel_corr_opt_semipar[f] == "corr_zero_b")    map_retsel_pe_pars[,1,,f]     <- NA
                if (ret_sel_corr_opt_semipar[f] == "corr_zero_y_b")  map_retsel_pe_pars[,1:2,,f]   <- NA
              }

              if(opt %in% c(3,4)) {
                # 3d gmrf options only (adds the cohort dimension)
                if (ret_sel_corr_opt_semipar[f] == "corr_zero_c")      map_retsel_pe_pars[,3,,f]   <- NA
                if (ret_sel_corr_opt_semipar[f] == "corr_zero_y_c")    map_retsel_pe_pars[,2:3,,f] <- NA
                if (ret_sel_corr_opt_semipar[f] == "corr_zero_b_c")    map_retsel_pe_pars[,c(1,3),,f] <- NA
                if (ret_sel_corr_opt_semipar[f] == "corr_zero_y_b_c")  map_retsel_pe_pars[,1:3,,f] <- NA
              }

              # Reset numbering for mapping off correlation parameters for clarity
              non_na_positions <- which(!is.na(map_retsel_pe_pars))
              map_retsel_pe_pars[non_na_positions] <- seq_along(non_na_positions)
              collect_message("ret_sel_corr_opt_semipar is specified as: ", ret_sel_corr_opt_semipar[f], "for fishery fleet", f)

            }
          } # end if 3d gmrf marginal or conditional variance

          # fix all parameters
          if(retsel_pe_pars_spec[f] == "fix") map_retsel_pe_pars[r,,s,f] <- NA

        } # end s loop
      } # end else
    } # end r loop

    if(!is.null(retsel_pe_pars_spec)) collect_message("retsel_pe_pars_spec is specified as: ", retsel_pe_pars_spec[f], "for fishery fleet", f)

  } # end f loop

  # Handle fleet sharing after all base mappings are established
  for(f in 1:input_list$data$n_fish_fleets) {
    if(stringr::str_detect(retsel_pe_pars_spec[f], "est_shared_f") && !is.null(retsel_pe_pars_spec)) {
      # extract fleet sharing index
      flt_shared <- as.numeric(unlist(stringr::str_extract_all(retsel_pe_pars_spec[f], "\\d+")))

      # Validate options here
      if(flt_shared > input_list$data$n_fish_fleets || flt_shared < 1) stop("Fleet sharing specification 'est_shared_f", flt_shared, "' for fleet ", f, " references invalid fleet number.")
      if(stringr::str_detect(retsel_pe_pars_spec[flt_shared], "est_shared_f")) stop("Fleet ", f, " cannot share with fleet ", flt_shared, " because fleet ", flt_shared, " is self-sharing parameters, which does not make sense.")

      # Copy mapping from reference fleet
      map_retsel_pe_pars[,,,f] <- map_retsel_pe_pars[,,,flt_shared]
      collect_message("retsel_pe_pars_spec is specified as: ", retsel_pe_pars_spec[f], " for fishery fleet ", f, " (sharing with fleet ", flt_shared, ")")
    } # end if statement
  } # end f loop

  # input into mapping list
  input_list$map$retsel_pe_pars <- factor(map_retsel_pe_pars)

  return(input_list)
}

#' Map retained fishery selectivity deviation parameters
#'
#' Constructs the factor map for \code{ln_retsel_devs}, the annual deviations in
#' continuous time-varying retained fishery selectivity. For iid and random-walk forms,
#' the deviation dimension corresponds to selectivity parameters (up to 6 for
#' double-normal); for semi-parametric forms (3D GMRF, 2D AR1), it corresponds
#' to age or length bins. Cells with no time-variation
#' (\code{cont_tv_ret_sel == 0}) or no catch data are mapped to \code{NA}.
#'
#' Bin-sharing (\code{"est_shared_b"} and related options) groups bins into
#' blocks defined by \code{retsel_devs_shared_bins}, reducing the number of
#' estimated deviation series. Fleet sharing (\code{"est_shared_f_x"}) is handled
#' in a second pass. The resulting integer map is also stored as
#' \code{$data$map_ln_retsel_devs} for use in the RTMB objective function.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists.
#' @param ret_sel_devs_spec Character vector of length \code{n_fish_fleets}.
#'   Options:
#'   \describe{
#'     \item{\code{"est_all"}}{Separate deviation series per region × sex × bin.}
#'     \item{\code{"est_shared_r"}}{Shared across regions.}
#'     \item{\code{"est_shared_s"}}{Shared across sexes.}
#'     \item{\code{"est_shared_r_s"}}{Shared across regions and sexes.}
#'     \item{\code{"est_shared_b"}}{Shared across bin groups defined by
#'       \code{retsel_devs_shared_bins}.}
#'     \item{\code{"est_shared_r_b"}}{Shared across regions and bin groups.}
#'     \item{\code{"est_shared_b_s"}}{Shared across bin groups and sexes.}
#'     \item{\code{"est_shared_r_b_s"}}{Shared across regions, bin groups, and sexes.}
#'     \item{\code{"est_shared_f_x"}}{Copy deviation map from fleet \code{x}.}
#'     \item{\code{"fix"} or \code{"none"}}{All deviations fixed at zero (mapped to \code{NA}).}
#'   }
#'   Bin-sharing options (\code{"est_shared_b"}, etc.) are only valid with
#'   semi-parametric time-varying forms (\code{cont_tv_ret_sel} is in \code{c(3,4,5)}).
#' @param retsel_devs_shared_bins List of integer vectors defining bin
#'   groupings for bin-sharing options. Each element groups bins that share a
#'   single deviation series, e.g., \code{list(1:5, 6:10, 11:30)}.
#'   Only used when \code{ret_sel_devs_spec} includes \code{"est_shared_b"}.
#' @param bins Number of selectivity bins
#'
#' @return The input \code{input_list} with \code{$map$ln_retsel_devs} set to a
#'   factor vector and \code{$data$map_ln_retsel_devs} set to the corresponding
#'   integer array (for use in the RTMB template).
#'
#' @keywords internal
do_retsel_devs_mapping <- function(input_list, ret_sel_devs_spec, retsel_devs_shared_bins, bins) {

  # Initialize counter and mapping array for retained fishery selectivity deviations
  retsel_devs_counter <- 1
  map_retsel_devs <- input_list$par$ln_retsel_devs
  map_retsel_devs[] <- NA

  for(r in 1:input_list$data$n_regions) {
    for(f in 1:input_list$data$n_fish_fleets) {

      # Validate options
      if(!is.null(ret_sel_devs_spec)) {
        if(!ret_sel_devs_spec[f] %in% c("fix", "none", "est_all", "est_shared_r", "est_shared_s", "est_shared_r_s", "est_shared_b", "est_shared_r_b", "est_shared_r_b_s", "est_shared_b_s") &&
           !stringr::str_detect(ret_sel_devs_spec[f], "est_shared_f_\\d+"))
          stop("ret_sel_devs_spec not correctly specfied. Should be one of these: est_all, est_shared_r, est_shared_r_s, est_shared_s, est_shared_b, est_shared_r_b, est_shared_r_b_s, est_shared_r_s, fix, or est_shared_f_# (where # is fleet number)")
        if(ret_sel_devs_spec[f] %in% c("est_shared_b", "est_shared_r_b", "est_shared_r_b_s", "est_shared_b_s") &&
           !input_list$data$cont_tv_ret_sel[r,f] %in% c(3,4,5)) stop("Sharing age deviations with iid or random walk parametric forms is not supported!")
      }

      # Skip fleet sharing specs in first pass
      if(!is.null(ret_sel_devs_spec)) if(stringr::str_detect(ret_sel_devs_spec[f], "est_shared_f")) next

      for(s in 1:input_list$data$n_sexes) {
        for(y in 1:(length(input_list$data$years) + input_list$data$n_proj_yrs_devs)) {

          # if no time-variation, then fix all parameters for this fleet
          if(input_list$data$cont_tv_ret_sel[r,f] == 0 || (sum(input_list$data$UseCatch[r,,,f]) == 0 && sum(input_list$data$UseCatch_pop[,r,,,f]) == 0) ) {
            map_retsel_devs[r,y,,s,f] <- NA
          } else {

            # Figure out max number of selectivity parameters for a given region and fleet
            if(unique(input_list$data$ret_sel_model[r,,f]) %in% 2) max_sel_pars <- 1 # exponential
            if(unique(input_list$data$ret_sel_model[r,,f]) %in% c(0,1,3)) max_sel_pars <- 2 # logistic or gamma
            if(unique(input_list$data$ret_sel_model[r,,f]) == 4) max_sel_pars <- 6 # double normal
            if(unique(input_list$data$ret_sel_model[r,,f]) == 5) max_sel_pars <- bins # non-parametric selectivity
            if(unique(input_list$data$ret_sel_model[r,,f]) %in% c(6,7)) max_sel_pars <- 3 # logistic selectivity with asymptote

            # If iid or random walk time-variation for this fleet
            if(input_list$data$cont_tv_ret_sel[r,f] %in% c(1,2)) {

              for(i in 1:max_sel_pars) {
                # Estimating all selectivity deviations across regions, sexes, fleets, and parameter
                if(ret_sel_devs_spec[f] == 'est_all') {
                  map_retsel_devs[r,y,i,s,f] <- retsel_devs_counter
                  retsel_devs_counter <- retsel_devs_counter + 1
                }

                # Estimating selectivity deviations across sexes, fleets, and parameters, but shared across regions
                if(ret_sel_devs_spec[f] == 'est_shared_r' && r == 1) {
                  map_retsel_devs[,y,i,s,f] <- retsel_devs_counter
                  retsel_devs_counter <- retsel_devs_counter + 1
                }

                # Estimating selectivity deviations across regions, fleets, and parameters, but shared across sexes
                if(ret_sel_devs_spec[f] == 'est_shared_s' && s == 1) {
                  map_retsel_devs[r,y,i,,f] <- retsel_devs_counter
                  retsel_devs_counter <- retsel_devs_counter + 1
                }

                # Estimating selectivity deviations across fleets, and parameters, but shared across sexes and regions
                if(ret_sel_devs_spec[f] == 'est_shared_r_s' && r == 1 && s == 1) {
                  map_retsel_devs[,y,i,,f] <- retsel_devs_counter
                  retsel_devs_counter <- retsel_devs_counter + 1
                }

              } # end i loop
            } # end iid or random walk variation

            # If 3d gmrf for this fleet
            if(input_list$data$cont_tv_ret_sel[r,f] %in% c(3,4,5)) {

              for(i in 1:length(input_list$data$ages)) {
                # Estimating all selectivity deviations across regions, years and bins
                if(ret_sel_devs_spec[f] == 'est_all') {
                  map_retsel_devs[r,y,i,s,f] <- retsel_devs_counter
                  retsel_devs_counter <- retsel_devs_counter + 1
                }

                # Estimating all selectivity deviations across years and bins, but shared across regions
                if(ret_sel_devs_spec[f] == 'est_shared_r' && r == 1) {
                  map_retsel_devs[,y,i,s,f] <- retsel_devs_counter
                  retsel_devs_counter <- retsel_devs_counter + 1
                }

                # Estimating all selectivity deviations across years and bins, but shared across sexes
                if(ret_sel_devs_spec[f] == 'est_shared_s' && s == 1) {
                  map_retsel_devs[r,y,i,,f] <- retsel_devs_counter
                  retsel_devs_counter <- retsel_devs_counter + 1
                }

                # Estimating all selectivity deviations across years and bins, but shared across sexes and regions
                if(ret_sel_devs_spec[f] == 'est_shared_r_s' && s == 1 && r == 1) {
                  map_retsel_devs[,y,i,,f] <- retsel_devs_counter
                  retsel_devs_counter <- retsel_devs_counter + 1
                }

                if(ret_sel_devs_spec[f] == 'est_shared_b') {
                  for(k in 1:length(retsel_devs_shared_bins)) {
                    map_retsel_devs[r,y,retsel_devs_shared_bins[[k]],s,f] <- retsel_devs_counter
                    retsel_devs_counter <- retsel_devs_counter + 1
                  } # end k loop
                }

                if(ret_sel_devs_spec[f] == 'est_shared_r_b' && r == 1) {
                  for(k in 1:length(retsel_devs_shared_bins)) {
                    map_retsel_devs[,y,retsel_devs_shared_bins[[k]],s,f] <- retsel_devs_counter
                    retsel_devs_counter <- retsel_devs_counter + 1
                  } # end k loop
                }

                if(ret_sel_devs_spec[f] == 'est_shared_b_s' && s == 1) {
                  for(k in 1:length(retsel_devs_shared_bins)) {
                    map_retsel_devs[r,y,retsel_devs_shared_bins[[k]],,f] <- retsel_devs_counter
                    retsel_devs_counter <- retsel_devs_counter + 1
                  } # end k loop
                }

                if(ret_sel_devs_spec[f] == 'est_shared_r_b_s' && s == 1 && r == 1) {
                  for(k in 1:length(retsel_devs_shared_bins)) {
                    map_retsel_devs[,y,retsel_devs_shared_bins[[k]],,f] <- retsel_devs_counter
                    retsel_devs_counter <- retsel_devs_counter + 1
                  } # end k loop
                }

              } # end i loop
            } # end 3d gmrf

          } # end else
        } # end y loop
      } # end s loop

      if(!is.null(ret_sel_devs_spec)) collect_message("ret_sel_devs_spec is specified as: ", ret_sel_devs_spec[f], "for fishery fleet", f, "and region ", r)

    } # end f loop
  } # end r loop

  # Handle fleet sharing after all base mappings are established
  for(f in 1:input_list$data$n_fish_fleets) {
    if(stringr::str_detect(ret_sel_devs_spec[f], "est_shared_f") && !is.null(ret_sel_devs_spec)) {
      # extract fleet sharing index
      flt_shared <- as.numeric(unlist(stringr::str_extract_all(ret_sel_devs_spec[f], "\\d+")))

      # Validate options here
      if(flt_shared > input_list$data$n_fish_fleets || flt_shared < 1) stop("Fleet sharing specification 'est_shared_f", flt_shared, "' for fleet ", f, " references invalid fleet number.")
      if(stringr::str_detect(ret_sel_devs_spec[flt_shared], "est_shared_f")) stop("Fleet ", f, " cannot share with fleet ", flt_shared, " because fleet ", flt_shared, " is self-sharing parameters, which does not make sense.")

      # Copy mapping from reference fleet
      map_retsel_devs[,,,,f] <- map_retsel_devs[,,,,flt_shared]
      collect_message("ret_sel_devs_spec is specified as: ", ret_sel_devs_spec[f], " for fishery fleet ", f, " (sharing with fleet ", flt_shared, ")")
    } # end if statement
  } # end f loop

  # input into mapping list
  input_list$map$ln_retsel_devs <- factor(map_retsel_devs)
  input_list$data$map_ln_retsel_devs <- array(as.numeric(input_list$map$ln_retsel_devs), dim = dim(input_list$par$ln_retsel_devs))

  return(input_list)
}

#' Set up retained fishery selectivity
#'
#' Configures all aspects of retained fishery selectivity and catchability for the
#' estimation model, including functional forms, time blocks, continuous time-varying
#' selectivity, process error structure, annual deviations, and fixed or estimated
#' selectivity parameters. Must be called after
#' \code{\link{Setup_Mod_FishIdx_and_Comps}}.
#'
#' Selectivity time-variation via \code{cont_tv_ret_sel} and blocked selectivity via
#' \code{ret_sel_blocks} are mutually exclusive within a fleet. Specifying both for
#' the same fleet will result in an error.
#'
#' @param input_list Named list containing \code{$data}, \code{$par}, \code{$map},
#'   and \code{$verbose} sublists as created by upstream setup functions.
#'
#' @param cont_tv_ret_sel Character vector defining continuous time-varying selectivity
#'   structure per fleet. Format:
#'   \code{"<type>_Fleet_<f>"} where \code{type} is one of:
#'   \describe{
#'     \item{\code{"none"}}{No time variation}
#'     \item{\code{"iid"}}{Independent year effects}
#'     \item{\code{"rw"}}{Random walk process}
#'     \item{\code{"3dmarg"}}{3D marginal structure}
#'     \item{\code{"3dcond"}}{3D conditional structure}
#'     \item{\code{"2dar1"}}{2D AR1 structure}
#'   }
#'   Output is mapped to integer codes internally.
#'
#' @param ret_sel_blocks Character vector defining discrete selectivity blocks.
#'   Format:
#'   \code{"none_Fleet_<f>"} or
#'   \code{"Block_<b>_Year_<start>-<end>_Fleet_<f>"} (or terminal year version).
#'   Values are expanded into a 3D array:
#'   \code{[n_regions × n_years × n_fish_fleets]}.
#'
#' @param ret_sel_model Character vector defining selectivity functional forms.
#'   Format:
#'   \code{"<type>_Fleet_<f>"} or
#'   \code{"<type>_Fleet_<f>_Block_<b>"}.
#'   Supported types:
#'   \describe{
#'     \item{\code{"logist1"}}{Logistic with \eqn{a_{50}} and slope \eqn{k} (2 parameters).}
#'     \item{\code{"logist2"}}{Logistic with \eqn{a_{50}} and \eqn{a_{95}} (2 parameters).}
#'     \item{\code{"gamma"}}{Dome-shaped gamma with \eqn{a_{max}} and \eqn{\delta} (2 parameters).}
#'     \item{\code{"exponential"}}{Exponential with a single power parameter (1 parameter).}
#'     \item{\code{"dbnrml"}}{Double-normal with 6 parameters.}
#'     \item{\code{"nonpar"}}{Non-parametric selectivity defined over discrete age or length bins, where selectivity is estimated as independent parameters (or grouped bins if specified via nonparametric bin mapping). No fixed functional form is imposed.}
#'     \item{\code{"asymplogist1"}}{Logistic selectivity with \eqn{a_{50}} and slope \eqn{k} and asymptotic control (3 parameters).}
#'     \item{\code{"asymplogist2"}}{Logistic selectivity with with \eqn{a_{50}} and \eqn{a_{95}} and asymptotic control (3 parameters).}
#'     \item{\code{"bicubic"}}{Bicubic spline over a bin-node x year-node grid; see \code{ret_sel_model} in \code{\link{Setup_Mod_Fishsel_and_Q}} for full syntax.}
#'   }
#'
#' @param retsel_pe_pars_spec Specification of process error parameters for
#'   time-varying selectivity. Length must equal \code{n_fish_fleets}.
#'
#' @param ret_fixed_sel_pars_spec Specification controlling which fixed
#'   selectivity parameters are estimated vs fixed.
#'
#' @param ret_sel_devs_spec Specification of annual selectivity deviations
#'   structure per fleet.
#'
#' @param ret_sel_corr_opt_semipar Optional correlation structure for semi-parametric
#'   selectivity deviations. Length must equal \code{n_fish_fleets}.
#'
#' @param Use_ret_selex_prior Integer flag (\code{0/1}) indicating whether
#'   selectivity priors are used.
#'
#' @param ret_selex_prior Data frame of priors for selectivity parameters.
#'   Must include columns: \code{region}, \code{fleet}, \code{block}, \code{sex},
#'   \code{par}, \code{mu}, \code{sd}.
#'
#' @param retsel_devs_shared_bins Vector defining shared bins for selectivity
#'   deviation estimation (e.g., age or length grouping structure).
#'
#' @param ret_selex_type Character string indicating selectivity domain:
#'   \code{"age"} or \code{"length"}.
#'
#' @param use_fixed_ret_sel Binary indicator array
#'   \code{[n_pop × n_regions × n_years × n_seas × n_fish_fleets]}.
#'   \code{1} = use fixed selectivity input; \code{0} = estimate.
#'
#' @param ret_sel_input Fixed selectivity input array:
#'   \code{[n_pop × n_regions × n_years × n_seas × (ages or lengths) × n_sexes × n_fish_fleets]}.
#'
#' @param ret_sel_nonpar_est_bins Vector defining estimated bins for
#'   non-parametric selectivity.
#'
#' @param ... Optional starting values for selectivity parameters and deviations.
#'
#' @return The updated \code{input_list} with:
#' \itemize{
#'   \item Parsed selectivity structure arrays in \code{$data}
#'   \item Starting parameter arrays in \code{$par}
#'   \item Factor mapping objects in \code{$map}
#' }
#'
#' @keywords internal
#' @importFrom stringr str_detect
#' @family Model Setup
Setup_Mod_Retsel <- function(input_list,
                             cont_tv_ret_sel,
                             ret_sel_blocks,
                             ret_sel_model,
                             retsel_pe_pars_spec,
                             ret_fixed_sel_pars_spec,
                             ret_sel_devs_spec,
                             ret_sel_corr_opt_semipar,
                             Use_ret_selex_prior,
                             ret_selex_prior,
                             retsel_devs_shared_bins,
                             ret_selex_type,
                             use_fixed_ret_sel,
                             ret_sel_input,
                             ret_sel_nonpar_est_bins,
                             ...
) {

  messages_list <<- character(0) # string to attach to for printing messages
  starting_values <- list(...)
  if(input_list$store_config) input_list$config$Setup_Mod_Fishsel_and_Q <- c(input_list$config$Setup_Mod_Fishsel_and_Q, mget(names(formals()))[-1])


  # Input Validation --------------------------------------------------------
  # Continuous Selectivity Deviations
  if(!is.null(retsel_pe_pars_spec)) if(length(retsel_pe_pars_spec) != input_list$data$n_fish_fleets) stop("retsel_pe_pars_spec is not length n_fish_fleets")
  if(!is.null(ret_sel_devs_spec)) if(length(ret_sel_devs_spec) != input_list$data$n_fish_fleets) stop("ret_sel_devs_spec is not length n_fish_fleets")
  if(!is.null(ret_sel_corr_opt_semipar)) if(length(ret_sel_corr_opt_semipar) != input_list$data$n_fish_fleets) stop("ret_sel_corr_opt_semipar is not length n_fish_fleets")

  # Selectivity Priors
  if(!Use_ret_selex_prior %in% c(0,1)) stop("Values for Use_ret_selex_prior are not valid. They are == 0 (don't use prior), or == 1 (use prior)")
  # Checking selectivity priors
  if(Use_ret_selex_prior == 1) {
    required_cols <- c("region", "fleet", "block", "sex", "par", "mu", "sd")
    missing_cols <- setdiff(required_cols, names(ret_selex_prior))
    if(length(missing_cols) > 0) {
      stop("ret_selex_prior is missing required columns: ", paste(missing_cols, collapse = ", "))
    }
  }
  collect_message("Retained Fishery Selectivity priors are: ", ifelse(Use_ret_selex_prior == 0, "Not Used", "Used"))

  if(any(use_fixed_ret_sel == 1) && is.null(ret_sel_input)) stop("ret_sel_input is NULL, please provide an input array.")
  if(any(use_fixed_ret_sel == 1) && ret_selex_type == 'age') check_data_dimensions(ret_sel_input, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_ages = length(input_list$data$ages), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ret_sel_input_age')
  if(any(use_fixed_ret_sel == 1) && ret_selex_type == 'length') check_data_dimensions(ret_sel_input, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_lens = length(input_list$data$lens), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ret_sel_input_len')

  # Selectivity Options -----------------------------------------------------
  # Age based selectivity
  if(ret_selex_type == 'age') {
    ret_selex_type <- 0
    bins <- length(input_list$data$ages)
    collect_message("Retained Fishery Selectivity is aged-based.")
  } # if age based

  # Length based selectivity
  if(ret_selex_type == 'length') {
    if(input_list$data$fit_lengths == 0) stop("Length composition data are not fit, but retained selectivity is length-based. This is not allowed. Please change to a valid option (either fit lengths or use age-based selectivity).")
    ret_selex_type <- 1
    bins <- length(input_list$data$lens)
    collect_message("Retained Fishery Selectivity is length-based")
  } # if length based


  # Continuous Retained Time-Varying Selectivity Options -----------------------------
  cont_tv_ret_sel_mat <- array(NA, dim = c(input_list$data$n_regions, input_list$data$n_fish_fleets))
  cont_tv_map <- data.frame(type = c("none", "iid", "rw", "3dmarg", "3dcond", "2dar1"), num = c(0,1,2,3,4,5)) # set up values we map to

  for(i in 1:length(cont_tv_ret_sel)) {
    # Extract out components from list
    tmp <- cont_tv_ret_sel[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    cont_tv_type <- tmp_vec[1] # get continuous selex type
    fleet <- as.numeric(tmp_vec[3]) # extract fleet index

    # Validate options
    if(!fleet %in% c(1:input_list$data$n_fish_fleets)) stop("Invalid fleet specified for cont_tv_ret_sel This needs to be specified as timevarytype_Fleet_x")
    if(!cont_tv_type %in% c(cont_tv_map$type)) stop("cont_tv_ret_sel is not correctly specified. This needs to be one of these: none, iid, rw, 3dmarg, 3dcond, 2dar1 (the timevarytypes) and specified as timevarytype_Fleet_x")

    # Input options
    cont_tv_ret_sel_mat[,fleet] <- cont_tv_map$num[which(cont_tv_map$type == cont_tv_type)]
    collect_message("Continuous retained fishery time-varying selectivity specified as: ", cont_tv_type, " for fishery fleet ", fleet)
  }

  if(any(cont_tv_ret_sel_mat > 0) && is.null(retsel_pe_pars_spec) && is.null(ret_sel_devs_spec)) stop("Continuous time-varying selectivity specified, but retsel_pe_pars_spec and/or ret_sel_devs_spec is NULL (i.e., not specified)!")

  # Blocked Retained Time-Varying Selectivity Options --------------------------------
  ret_sel_blocks_arr <- array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets))
  for(i in 1:length(ret_sel_blocks)) {

    # Extract out components from list
    tmp <- ret_sel_blocks[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))

    # Validate options
    if(!tmp_vec[1] %in% c("none", "Block")) stop("Retained Fishery Selectivity Blocks not correctly specified. This should be either none_Fleet_x or Block_x_Year_x-y_Fleet_x")

    # extract out fleets if constant
    if(tmp_vec[1] == "none") {
      fleet <- as.numeric(tmp_vec[3]) # get fleet number
      ret_sel_blocks_arr[,,fleet] <- 1 # input only 1 fishery time block
    }

    if(tmp_vec[1] == "Block") {

      block_val <- as.numeric(tmp_vec[2]) # get block value
      fleet <- as.numeric(tmp_vec[6]) # extract fleet index

      # get year ranges
      if(!str_detect(tmp, "terminal")) { # if not terminal year
        year_range <- as.numeric(unlist(strsplit(tmp_vec[4], "-")))
        years <- year_range[1]:year_range[2] # get sequence of years
      } else { # if terminal year
        year_range <- unlist(strsplit(tmp_vec[4], '-'))[1] # get year range
        years <- as.numeric(year_range):length(input_list$data$years) # get sequence of years
      }

      ret_sel_blocks_arr[,years,fleet] <- block_val
    }

  }

  if(any(is.na(ret_sel_blocks_arr))) stop("Retained Fishery Selectivtiy Blocks are returning an NA. Did you forget to specify the year range of ret_sel_blocks?")
  for(f in 1:input_list$data$n_fish_fleets) collect_message(paste("Retained Fishery Selectivity Time Blocks for fishery", f, "is specified at:", length(unique(ret_sel_blocks_arr[,,f]))))

  # Retained Selectivity Functional Forms --------------------------------------------
  sel_map <- data.frame(sel = c('logist1', "gamma", "exponential", "logist2", "dbnrml", 'nonpar', 'asymplogist1', "asymplogist2", "bicubic"), num = c(0,1,2,3,4,5,6,7,8)) # set up values we can map to
  ret_sel_model_arr <- array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets))
  ret_sel_bicubic_binnodes_arr <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets)) # number of bin nodes, only set where ret_sel_model == 8 (bicubic)
  ret_sel_bicubic_yrnodes_arr <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets)) # number of year nodes, only set where ret_sel_model == 8 (bicubic)
  ret_sel_bicubic_selstyr_arr <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets)) # calendar year the bicubic surface is actually fit from (0 = block's own start year); years before this are edge-held, matching fish_sel_model's SelStyr
  ret_sel_bicubic_nselbins_arr <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets)) # number of bins the bicubic surface is actually fit over (0 = all bins); bins beyond this are edge-held, matching fish_sel_model's NSelBins

  for(i in 1:length(ret_sel_model)) {

    # Extract out retained fishery selectivity components from vector
    tmp_sel_form <- ret_sel_model[i]
    tmp_sel_form_vec <- unlist(strsplit(tmp_sel_form, "_")) # split string
    sel_form <- tmp_sel_form_vec[1] # get selectivity type

    if(sel_form == "bicubic") {
      # bicubic spline: bicubic_Bin_<n_bin_nodes>_Yr_<n_yr_nodes>_Fleet_<f>[_Block_<b>][_SelStyr_<year>][_NSelBins_<n>]
      bin_pos <- which(tmp_sel_form_vec == "Bin")
      yr_pos <- which(tmp_sel_form_vec == "Yr")
      fleet_pos <- which(tmp_sel_form_vec == "Fleet")
      block_pos <- which(tmp_sel_form_vec == "Block")
      selstyr_pos <- which(tmp_sel_form_vec == "SelStyr")
      nselbins_pos <- which(tmp_sel_form_vec == "NSelBins")
      if(length(bin_pos) != 1 || length(yr_pos) != 1 || length(fleet_pos) != 1)
        stop("ret_sel_model 'bicubic' entries must be specified as bicubic_Bin_<n_bin_nodes>_Yr_<n_yr_nodes>_Fleet_<f> or bicubic_Bin_<n_bin_nodes>_Yr_<n_yr_nodes>_Block_<b>_Fleet_<f>, optionally with _SelStyr_<year> and/or _NSelBins_<n>")
      tmp_n_bin_nodes <- suppressWarnings(as.numeric(tmp_sel_form_vec[bin_pos + 1]))
      tmp_n_yr_nodes <- suppressWarnings(as.numeric(tmp_sel_form_vec[yr_pos + 1]))
      tmp_fleet <- suppressWarnings(as.numeric(tmp_sel_form_vec[fleet_pos + 1]))
      tmp_block <- if(length(block_pos) == 1) suppressWarnings(as.numeric(tmp_sel_form_vec[block_pos + 1])) else NULL
      tmp_selstyr <- if(length(selstyr_pos) == 1) suppressWarnings(as.numeric(tmp_sel_form_vec[selstyr_pos + 1])) else 0
      tmp_nselbins <- if(length(nselbins_pos) == 1) suppressWarnings(as.numeric(tmp_sel_form_vec[nselbins_pos + 1])) else 0
      if(is.na(tmp_n_bin_nodes) || tmp_n_bin_nodes < 2) stop("bicubic ret_sel_model requires at least 2 bin nodes (n_bin_nodes >= 2)")
      if(is.na(tmp_n_yr_nodes) || tmp_n_yr_nodes < 1) stop("bicubic ret_sel_model requires at least 1 year node (n_yr_nodes >= 1). Use n_yr_nodes == 1 for a time-invariant bin-only spline.")
      if(length(selstyr_pos) == 1 && (is.na(tmp_selstyr) || !tmp_selstyr %in% input_list$data$years)) stop("bicubic ret_sel_model SelStyr must be a calendar year within the modeled years")
      if(length(nselbins_pos) == 1 && (is.na(tmp_nselbins) || tmp_nselbins < 2 || tmp_nselbins > bins)) stop("bicubic ret_sel_model NSelBins must be an integer between 2 and the total number of bins (ages or lengths)")
    } else {
      # get fleet index
      tmp_fleet <- if(length(tmp_sel_form_vec) == 3) as.numeric(tmp_sel_form_vec[3]) else as.numeric(tmp_sel_form_vec[5]) # fleet index changes if block is included in character vector
      # get block index
      tmp_block <- if(length(tmp_sel_form_vec) == 5) as.numeric(tmp_sel_form_vec[3]) else NULL
    }

    # validate options
    if(!sel_form %in% c(sel_map$sel)) stop("ret_sel_model is not correctly specified. This needs to be one of these: logist1, gamma, exponential, logist2, dbnrml, nonpar, asymplogist1, asymplogist2, bicubic (the seltypes) and specified as seltype_Fleet_x")
    if(!tmp_fleet %in% c(1:input_list$data$n_fish_fleets)) stop("Invalid fleet specified for ret_sel_model This needs to be specified as seltype_Fleet_x or seltype_Fleet_x_Block_x (if blocks are specified to change for a fleet)")

    # Input options
    if(is.null(tmp_block)) {
      ret_sel_model_arr[,,tmp_fleet] <- sel_map$num[which(sel_map$sel == sel_form)] # same selectivity form across blocks
      if(sel_form == "bicubic") {
        ret_sel_bicubic_binnodes_arr[,,tmp_fleet] <- tmp_n_bin_nodes
        ret_sel_bicubic_yrnodes_arr[,,tmp_fleet] <- tmp_n_yr_nodes
        ret_sel_bicubic_selstyr_arr[,,tmp_fleet] <- tmp_selstyr
        ret_sel_bicubic_nselbins_arr[,,tmp_fleet] <- tmp_nselbins
      }
    } else {
      ret_sel_model_arr[,which(ret_sel_blocks_arr[,,tmp_fleet] == tmp_block),tmp_fleet] <- sel_map$num[which(sel_map$sel == sel_form)]
      if(sel_form == "bicubic") {
        ret_sel_bicubic_binnodes_arr[,which(ret_sel_blocks_arr[,,tmp_fleet] == tmp_block),tmp_fleet] <- tmp_n_bin_nodes
        ret_sel_bicubic_yrnodes_arr[,which(ret_sel_blocks_arr[,,tmp_fleet] == tmp_block),tmp_fleet] <- tmp_n_yr_nodes
        ret_sel_bicubic_selstyr_arr[,which(ret_sel_blocks_arr[,,tmp_fleet] == tmp_block),tmp_fleet] <- tmp_selstyr
        ret_sel_bicubic_nselbins_arr[,which(ret_sel_blocks_arr[,,tmp_fleet] == tmp_block),tmp_fleet] <- tmp_nselbins
      }
    }
    rm(tmp_block) # remove tmp block to start next loop
    collect_message("Retained Fishery selectivity functional form specified as:", sel_form, " for fishery fleet ", tmp_fleet)
  }

  # Validate that blocks and continuous time-variation aren't both specified for same fleet
  for(f in 1:input_list$data$n_fish_fleets) {
    has_blocks <- length(unique(ret_sel_blocks_arr[1,,f])) > 1
    has_cont_tv <- cont_tv_ret_sel_mat[1,f] != 0  # 0 = "none"
    if(has_blocks && has_cont_tv) {
      stop("Fleet ", f, " has both selectivity blocks and continuous time-varying selectivity specified. ",
           "These are mutually exclusive - choose one approach to time-variation.")
    }
  }

  # Populate Data List ------------------------------------------------------

  input_list$data$cont_tv_ret_sel <- cont_tv_ret_sel_mat
  input_list$data$ret_sel_blocks <- ret_sel_blocks_arr
  input_list$data$ret_sel_model <- ret_sel_model_arr
  input_list$data$ret_sel_bicubic_binnodes <- ret_sel_bicubic_binnodes_arr
  input_list$data$ret_sel_bicubic_yrnodes <- ret_sel_bicubic_yrnodes_arr
  input_list$data$ret_sel_bicubic_selstyr <- ret_sel_bicubic_selstyr_arr
  input_list$data$ret_sel_bicubic_nselbins <- ret_sel_bicubic_nselbins_arr
  input_list$data$Use_ret_selex_prior <- Use_ret_selex_prior
  input_list$data$ret_selex_prior <- ret_selex_prior
  input_list$data$ret_selex_type <- ret_selex_type
  input_list$data$use_fixed_ret_sel <- use_fixed_ret_sel
  input_list$data$ret_sel_input <- ret_sel_input
  input_list$data$retsel_devs_min_shared_bins <- if(!is.null(retsel_devs_shared_bins)) unlist(lapply(retsel_devs_shared_bins, min)) else 1:length(input_list$data$ages)

  # Populate Parameter List -------------------------------------------------

  # Figure out number of selectivity parameters for a given functional form
  unique_retsel_vals <- unique(as.vector(input_list$data$ret_sel_model))
  sel_pars_vec <- vector() # create empty vector to populate

  for(i in 1:length(unique_retsel_vals)) {
    if(unique_retsel_vals[i] %in% c(2)) sel_pars_vec[i] <- 1 # exponential
    if(unique_retsel_vals[i] %in% c(0,1,3)) sel_pars_vec[i] <- 2 # logistic or gamma
    if(unique_retsel_vals[i] %in% c(4)) sel_pars_vec[i] <- 6 # double normal
    if(unique_retsel_vals[i] == 5) sel_pars_vec[i] <- bins # non-parametric
    if(unique_retsel_vals[i] %in% c(6,7)) sel_pars_vec[i] <- 3 # logistic w/ asymptote parameter
    if(unique_retsel_vals[i] == 8) sel_pars_vec[i] <- max(input_list$data$ret_sel_bicubic_binnodes * input_list$data$ret_sel_bicubic_yrnodes) # bicubic: flattened bin-node x year-node grid
  } # end i loop

  # figure out maximum number of retained fishery selectivity blocks for a given reigon and fleet
  max_retsel_blks <- max(apply(input_list$data$ret_sel_blocks, c(1,3), FUN = function(x) length(unique(x))))
  # maximum number of selectivity parameters across all forms
  max_retsel_pars <- max(sel_pars_vec)
  if("ret_fixed_sel_pars" %in% names(starting_values)) input_list$par$ret_fixed_sel_pars <- starting_values$ret_fixed_sel_pars
  else input_list$par$ret_fixed_sel_pars <- array(0, dim = c(input_list$data$n_regions, max_retsel_pars, max_retsel_blks, input_list$data$n_sexes, input_list$data$n_fish_fleets))

  # Bicubic spline interpolation weight matrices (bin node x year node grid) for retention selectivity,
  # mirroring fish_sel_bicubic_Wbin/Wyr above (see Get_Selex documentation).
  has_bicubic_ret_sel <- any(input_list$data$ret_sel_model == 8)
  max_bin_nodes_bicubic_ret <- if(has_bicubic_ret_sel) max(input_list$data$ret_sel_bicubic_binnodes) else 1
  max_yr_nodes_bicubic_ret <- if(has_bicubic_ret_sel) max(input_list$data$ret_sel_bicubic_yrnodes) else 1
  n_yrs_total_bicubic_ret <- length(input_list$data$years) + input_list$data$n_proj_yrs_devs

  ret_sel_bicubic_Wbin <- array(0, dim = c(input_list$data$n_regions, bins, max_bin_nodes_bicubic_ret, max_retsel_blks, input_list$data$n_fish_fleets))
  ret_sel_bicubic_Wyr <- array(0, dim = c(input_list$data$n_regions, n_yrs_total_bicubic_ret, max_yr_nodes_bicubic_ret, max_retsel_blks, input_list$data$n_fish_fleets))

  if(has_bicubic_ret_sel) {
    for(f in 1:input_list$data$n_fish_fleets) {
      for(r in 1:input_list$data$n_regions) {

        retsel_blocks_tmp <- unique(as.vector(input_list$data$ret_sel_blocks[r,,f]))

        for(b in 1:length(retsel_blocks_tmp)) {

          block_years <- which(input_list$data$ret_sel_blocks[r,,f] == retsel_blocks_tmp[b])
          if(unique(input_list$data$ret_sel_model[r, block_years, f]) != 8) next # only bicubic blocks need weight matrices

          n_bin_nodes_this <- unique(input_list$data$ret_sel_bicubic_binnodes[r, block_years, f])
          n_yr_nodes_this <- unique(input_list$data$ret_sel_bicubic_yrnodes[r, block_years, f])

          # Bin dimension: see fish_sel_bicubic_Wbin construction above for NSelBins truncation/plateau details
          nselbins_this <- unique(input_list$data$ret_sel_bicubic_nselbins[r, block_years, f])
          n_fit_bins <- if(nselbins_this == 0) bins else nselbins_this

          bin_nodes_scaled <- seq(0, 1, length.out = n_bin_nodes_this)
          fit_bin_scaled <- seq(0, 1, length.out = n_fit_bins)
          Wbin_fit <- Get_Natural_Cubic_Spline_Weights(bin_nodes_scaled, fit_bin_scaled)

          Wbin_this <- matrix(0, nrow = bins, ncol = n_bin_nodes_this)
          Wbin_this[1:n_fit_bins, ] <- Wbin_fit
          if(n_fit_bins < bins) Wbin_this[(n_fit_bins + 1):bins, ] <- matrix(Wbin_fit[nrow(Wbin_fit), ], nrow = bins - n_fit_bins, ncol = n_bin_nodes_this, byrow = TRUE)

          ret_sel_bicubic_Wbin[r, , 1:n_bin_nodes_this, b, f] <- Wbin_this

          # Year dimension: see fish_sel_bicubic_Wyr construction above for SelStyr truncation/edge-hold details
          selstyr_this <- unique(input_list$data$ret_sel_bicubic_selstyr[r, block_years, f])
          selstyr_idx <- if(selstyr_this == 0) min(block_years) else which(input_list$data$years == selstyr_this)
          fit_years <- block_years[block_years >= selstyr_idx]
          pre_fit_years <- block_years[block_years < selstyr_idx]

          yr_nodes_scaled <- seq(0, 1, length.out = n_yr_nodes_this)
          fit_yr_scaled <- seq(0, 1, length.out = length(fit_years))
          Wyr_block <- Get_Natural_Cubic_Spline_Weights(yr_nodes_scaled, fit_yr_scaled)

          Wyr_this <- matrix(0, nrow = n_yrs_total_bicubic_ret, ncol = n_yr_nodes_this)
          Wyr_this[fit_years, ] <- Wyr_block
          if(length(pre_fit_years) > 0) Wyr_this[pre_fit_years, ] <- matrix(Wyr_block[1, ], nrow = length(pre_fit_years), ncol = n_yr_nodes_this, byrow = TRUE)
          if(min(block_years) > 1) Wyr_this[1:(min(block_years) - 1), ] <- matrix(Wyr_block[1, ], nrow = min(block_years) - 1, ncol = n_yr_nodes_this, byrow = TRUE)
          if(max(block_years) < n_yrs_total_bicubic_ret) Wyr_this[(max(block_years) + 1):n_yrs_total_bicubic_ret, ] <- matrix(Wyr_block[nrow(Wyr_block), ], nrow = n_yrs_total_bicubic_ret - max(block_years), ncol = n_yr_nodes_this, byrow = TRUE)

          ret_sel_bicubic_Wyr[r, , 1:n_yr_nodes_this, b, f] <- Wyr_this
        } # end b loop
      } # end r loop
    } # end f loop
  } # end if has_bicubic_ret_sel

  input_list$data$ret_sel_bicubic_Wbin <- ret_sel_bicubic_Wbin
  input_list$data$ret_sel_bicubic_Wyr <- ret_sel_bicubic_Wyr

  # Retained Fishery selectivity process error parameters
  if("retsel_pe_pars" %in% names(starting_values)) input_list$par$retsel_pe_pars <- starting_values$retsel_pe_pars
  else input_list$par$retsel_pe_pars <- array(0, dim = c(input_list$data$n_regions, max(max_retsel_pars, 4), input_list$data$n_sexes, input_list$data$n_fish_fleets)) # dimensioned 4 as the max number of pars for process errors (e.g., sigmas), and then just map off if not using

  # Retained Fishery selectivity deviations
  if(input_list$data$ret_selex_type == 0) bins <- length(input_list$data$ages) # age based deviations
  if(input_list$data$ret_selex_type == 1) bins <- length(input_list$data$lens) # length based deviations
  if("ln_retsel_devs" %in% names(starting_values)) input_list$par$ln_retsel_devs <- starting_values$ln_retsel_devs
  else input_list$par$ln_retsel_devs <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years) + input_list$data$n_proj_yrs_devs, bins, input_list$data$n_sexes, input_list$data$n_fish_fleets))

  # Mapping Options ---------------------------------------------------------
  input_list <- do_ret_fixed_sel_pars_mapping(input_list, ret_fixed_sel_pars_spec, bins, ret_sel_nonpar_est_bins)
  input_list <- do_retsel_pe_pars_mapping(input_list, retsel_pe_pars_spec, ret_sel_corr_opt_semipar, bins)
  input_list <- do_retsel_devs_mapping(input_list, ret_sel_devs_spec, retsel_devs_shared_bins, bins)

  return(input_list)
}



#' Set up total and retained fishery selectivity and catchability specifications
#'
#' Configures all aspects of fishery selectivity and catchability for the
#' estimation model: functional forms, time blocks, continuous time-varying
#' structures, process error hyperparameters, annual deviations, and
#' catchability blocks and estimation structure. Must be called after
#' \code{\link{Setup_Mod_FishIdx_and_Comps}}.
#'
#' Selectivity time-variation and blocked selectivity are mutually exclusive
#' within a fleet — specifying both for the same fleet will raise an error.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, \code{$map},
#'   and \code{$verbose} sublists.
#' @param fish_sel_model Character vector specifying the selectivity functional
#'   form for each fleet (and optionally each time block). Each element must
#'   follow one of:
#'   \itemize{
#'     \item \code{"<model>_Fleet_<f>"} — single form for all years of fleet \code{f}.
#'     \item \code{"<model>_Fleet_<f>_Block_<b>"} — form specific to block \code{b}
#'       of fleet \code{f}, as defined in \code{fish_sel_blocks}.
#'   }
#'   Available models:
#'   \describe{
#'     \item{\code{"logist1"}}{Logistic with \eqn{a_{50}} and slope \eqn{k} (2 parameters).}
#'     \item{\code{"logist2"}}{Logistic with \eqn{a_{50}} and \eqn{a_{95}} (2 parameters).}
#'     \item{\code{"gamma"}}{Dome-shaped gamma with \eqn{a_{max}} and \eqn{\delta} (2 parameters).}
#'     \item{\code{"exponential"}}{Exponential with a single power parameter (1 parameter).}
#'     \item{\code{"dbnrml"}}{Double-normal with 6 parameters.}
#'     \item{\code{"nonpar"}}{Non-parametric selectivity defined over discrete age or length bins, where selectivity is estimated as independent parameters (or grouped bins if specified via nonparametric bin mapping). No fixed functional form is imposed.}
#'     \item{\code{"asymplogist1"}}{Logistic selectivity with \eqn{a_{50}} and slope \eqn{k} and asymptotic control (3 parameters).}
#'     \item{\code{"asymplogist2"}}{Logistic selectivity with with \eqn{a_{50}} and \eqn{a_{95}} and asymptotic control (3 parameters).}
#'     \item{\code{"bicubic"}}{Bicubic spline over a bin-node x year-node grid
#'       (see \code{\link{Get_Selex}}, \code{Selex_Model == 8}). Specified as
#'       \code{"bicubic_Bin_<n_bin_nodes>_Yr_<n_yr_nodes>_Fleet_x"} (optionally
#'       with \code{_Block_k}). One generalized form covers a smooth bin x year
#'       surface (\code{n_yr_nodes > 1}), a time-invariant bin-only spline
#'       (\code{n_yr_nodes == 1}), or a bin-only spline re-fit independently
#'       per year-block (\code{n_yr_nodes == 1} within each of several blocks
#'       defined via \code{fish_sel_blocks}). An optional \code{_SelStyr_<year>}
#'       suffix (a calendar year within the block) restricts the actual spline
#'       fit to \code{SelStyr}:block-end; years within the block before
#'       \code{SelStyr} are held constant at the \code{SelStyr} year's fitted
#'       curve, rather than fitting the surface over the whole block. An
#'       optional \code{_NSelBins_<n>} suffix restricts the actual spline fit
#'       to the first \code{n} bins (ages or lengths, per \code{fish_selex_type});
#'       bins beyond \code{n} are held constant at the last fitted bin's
#'       curve.}
#'   }
#'   See the model equations vignette for mathematical definitions.
#' @param cont_tv_fish_sel Character vector of length \code{n_fish_fleets}
#'   specifying continuous time-varying selectivity per fleet. Each element
#'   must be \code{"<type>_Fleet_<f>"}. Valid types:
#'   \describe{
#'     \item{\code{"none"}}{No continuous time-variation (default).}
#'     \item{\code{"iid"}}{IID annual deviations on selectivity parameters.}
#'     \item{\code{"rw"}}{Random walk in selectivity parameters over time.}
#'     \item{\code{"3dmarg"}}{3D GMRF with marginal variance parameterisation.}
#'     \item{\code{"3dcond"}}{3D GMRF with conditional variance parameterisation.}
#'     \item{\code{"2dar1"}}{2D separable AR1 in bin and year dimensions.}
#'   }
#'   If any fleet has \code{cont_tv_fish_sel != "none"}, both
#'   \code{fishsel_pe_pars_spec} and \code{fish_sel_devs_spec} must also be
#'   provided.
#' @param fish_sel_blocks Character vector defining discrete selectivity time
#'   blocks per fleet. Each element follows \code{"Block_<b>_Year_<s>-<e>_Fleet_<f>"}
#'   or \code{"Block_<b>_Year_<s>-terminal_Fleet_<f>"}. Use
#'   \code{"none_Fleet_<f>"} (default) for a single constant block. Blocks must
#'   be non-overlapping and together span all model years for the specified fleet.
#'   Mutually exclusive with \code{cont_tv_fish_sel != "none"} for the same fleet.
#' @param fish_q_blocks Character vector defining catchability time blocks per
#'   fleet, using the same format as \code{fish_sel_blocks}. Default
#'   \code{"none_Fleet_<f>"} gives a single constant block.
#' @param fish_fixed_sel_pars_spec Character vector of length \code{n_fish_fleets}
#'   specifying how fixed-effect selectivity parameters are estimated. See
#'   \code{\link{do_fish_fixed_sel_pars_mapping}} for all options
#'   (\code{"est_all"}, \code{"est_shared_r"}, \code{"est_shared_s"},
#'   \code{"est_shared_r_s"}, \code{"est_shared_f_x"}, \code{"fix"}).
#' @param fish_q_spec Character vector of length \code{n_fish_fleets} specifying
#'   catchability estimation structure. See \code{\link{do_fish_q_mapping}} for
#'   options (\code{"est_all"}, \code{"est_shared_r"}, \code{"fix"}).
#' @param fishsel_pe_pars_spec Character vector of length \code{n_fish_fleets}
#'   specifying the estimation structure for selectivity process error
#'   hyperparameters. Required when any fleet has continuous time-variation.
#'   See \code{\link{do_fishsel_pe_pars_mapping}} for all options.
#' @param fish_sel_devs_spec Character vector of length \code{n_fish_fleets}
#'   specifying the estimation structure for annual selectivity deviations.
#'   Required when any fleet has continuous time-variation. See
#'   \code{\link{do_fishsel_devs_mapping}} for all options including age-sharing
#'   options for semi-parametric forms.
#' @param fishsel_devs_shared_bins List of integer vectors grouping age or length
#'   bins that share a single deviation series. Only used when
#'   \code{fish_sel_devs_spec} contains one of the \code{"est_shared_b"} variants.
#'   Example: \code{list(1:5, 6:10, 11:30)}.
#' @param corr_opt_semipar Character vector of length \code{n_fish_fleets}
#'   controlling which correlation components to suppress in semi-parametric
#'   (3D GMRF or 2D AR1) time-varying selectivity. Set to \code{NA} (default)
#'   for no suppression. See \code{\link{do_fishsel_pe_pars_mapping}} for valid
#'   suppression codes. Cohort-correlation options are invalid for \code{"2dar1"}.
#' @param Use_fish_q_prior Integer flag. \code{1} = apply lognormal priors to
#'   catchability; \code{0} = no priors (default). Requires \code{fish_q_prior}.
#' @param fish_q_prior Data frame of catchability prior hyperparameters. Required
#'   columns: \code{region}, \code{fleet}, \code{block} (block index), \code{mu}
#'   (prior mean on natural scale), \code{sd} (prior SD on log scale). Each row
#'   specifies a \eqn{\text{Normal}(\log(\mu), \sigma)} prior for one catchability
#'   parameter. Only used when \code{Use_fish_q_prior = 1}.
#' @param Use_fish_selex_prior Integer flag. \code{1} = apply lognormal priors to
#'   selectivity parameters; \code{0} = no priors (default). Requires
#'   \code{fish_selex_prior}.
#' @param fish_selex_prior Data frame of selectivity prior hyperparameters.
#'   Required columns: \code{region}, \code{fleet}, \code{block}, \code{sex},
#'   \code{par} (parameter index within the functional form), \code{mu}, \code{sd}.
#'   Only used when \code{Use_fish_selex_prior = 1}.
#' @param ... Optional starting value overrides for selectivity parameters.
#' @param cont_tv_ret_sel Character vector of length \code{n_fish_fleets}
#'   specifying continuous time-varying selectivity per fleet. Each element
#'   must be \code{"<type>_Fleet_<f>"}. Valid types:
#'   \describe{
#'     \item{\code{"none"}}{No continuous time-variation (default).}
#'     \item{\code{"iid"}}{IID annual deviations on selectivity parameters.}
#'     \item{\code{"rw"}}{Random walk in selectivity parameters over time.}
#'     \item{\code{"3dmarg"}}{3D GMRF with marginal variance parameterisation.}
#'     \item{\code{"3dcond"}}{3D GMRF with conditional variance parameterisation.}
#'     \item{\code{"2dar1"}}{2D separable AR1 in bin and year dimensions.}
#'   }
#'   If any fleet has \code{cont_tv_ret_sel != "none"}, both
#'   \code{retsel_pe_pars_spec} and \code{ret_sel_devs_spec} must also be
#'   provided.
#'
#' @param ret_sel_blocks Character vector defining discrete selectivity time
#'   blocks per fleet. Each element follows \code{"Block_<b>_Year_<s>-<e>_Fleet_<f>"}
#'   or \code{"Block_<b>_Year_<s>-terminal_Fleet_<f>"}. Use
#'   \code{"none_Fleet_<f>"} (default) for a single constant block. Blocks must
#'   be non-overlapping and together span all model years for the specified fleet.
#'   Mutually exclusive with \code{cont_tv_ret_sel != "none"} for the same fleet.
#'
#' @param ret_sel_model Character vector specifying the selectivity functional
#'   form for each fleet (and optionally each time block). Each element must
#'   follow one of:
#'   \itemize{
#'     \item \code{"<model>_Fleet_<f>"} — single form for all years of fleet \code{f}.
#'     \item \code{"<model>_Fleet_<f>_Block_<b>"} — form specific to block \code{b}
#'       of fleet \code{f}, as defined in \code{ret_sel_blocks}.
#'   }
#'   Available models:
#'   \describe{
#'     \item{\code{"logist1"}}{Logistic with \eqn{a_{50}} and slope \eqn{k} (2 parameters).}
#'     \item{\code{"logist2"}}{Logistic with \eqn{a_{50}} and \eqn{a_{95}} (2 parameters).}
#'     \item{\code{"gamma"}}{Dome-shaped gamma with \eqn{a_{max}} and \eqn{\delta} (2 parameters).}
#'     \item{\code{"exponential"}}{Exponential with a single power parameter (1 parameter).}
#'     \item{\code{"dbnrml"}}{Double-normal with 6 parameters.}
#'     \item{\code{"nonpar"}}{Non-parametric selectivity defined over discrete age or length bins, where selectivity is estimated as independent parameters (or grouped bins if specified via nonparametric bin mapping). No fixed functional form is imposed.}
#'     \item{\code{"asymplogist1"}}{Logistic selectivity with \eqn{a_{50}} and slope \eqn{k} and asymptotic control (3 parameters).}
#'     \item{\code{"asymplogist2"}}{Logistic selectivity with with \eqn{a_{50}} and \eqn{a_{95}} and asymptotic control (3 parameters).}
#'     \item{\code{"bicubic"}}{Bicubic spline over a bin-node x year-node grid, specified as
#'       \code{"bicubic_Bin_<n_bin_nodes>_Yr_<n_yr_nodes>_Fleet_x"} (optionally with \code{_Block_k},
#'       \code{_SelStyr_<year>}, and/or \code{_NSelBins_<n>}); see \code{fish_sel_model} above for the
#'       full syntax and \code{\link{Get_Selex}} (\code{Selex_Model == 8}) for the underlying math.}
#'   }
#'   See the model equations vignette for mathematical definitions.
#'
#' @param retsel_pe_pars_spec Character vector of length \code{n_fish_fleets}
#'   specifying the estimation structure for selectivity process error
#'   hyperparameters. Required when any fleet has continuous time-variation.
#'   See \code{\link{do_retsel_pe_pars_mapping}} for all options.
#'
#' @param ret_fixed_sel_pars_spec Character vector of length \code{n_fish_fleets}
#'   specifying how fixed-effect selectivity parameters are estimated. See
#'   \code{\link{do_ret_fixed_sel_pars_mapping}} for all options
#'   (\code{"est_all"}, \code{"est_shared_r"}, \code{"est_shared_s"},
#'   \code{"est_shared_r_s"}, \code{"est_shared_f_x"}, \code{"fix"}).
#'
#' @param ret_sel_devs_spec Character vector of length \code{n_fish_fleets}
#'   specifying the estimation structure for annual selectivity deviations.
#'   Required when any fleet has continuous time-variation. See
#'   \code{\link{do_retsel_devs_mapping}} for all options including age-sharing
#'   options for semi-parametric forms.
#'
#' @param retsel_devs_shared_bins List of integer vectors grouping age or length
#'   bins that share a single deviation series. Only used when
#'   \code{ret_sel_devs_spec} contains one of the \code{"est_shared_b"} variants.
#'   Example: \code{list(1:5, 6:10, 11:30)}.
#'
#' @param ret_sel_corr_opt_semipar Character vector of length \code{n_fish_fleets}
#'   controlling which correlation components to suppress in semi-parametric
#'   (3D GMRF or 2D AR1) time-varying selectivity. Set to \code{NA} (default)
#'   for no suppression. See \code{\link{do_retsel_pe_pars_mapping}} for valid
#'   suppression codes. Cohort-correlation options are invalid for \code{"2dar1"}.
#'
#' @param Use_ret_selex_prior Integer flag. \code{1} = apply lognormal priors to
#'   selectivity parameters; \code{0} = no priors (default). Requires
#'   \code{ret_selex_prior}.
#'
#' @param ret_selex_prior Data frame of selectivity prior hyperparameters.
#'   Required columns: \code{region}, \code{fleet}, \code{block}, \code{sex},
#'   \code{par}, \code{mu}, \code{sd}.
#'
#' @param use_fixed_ret_sel Integer vector of length \code{n_fish_fleets}
#'   indicating whether to fix selectivity (\code{1}) or estimate it (\code{0}).
#'
#' @param ret_sel_input Array of fixed selectivity values with dimensions
#'   \code{[n_pop × n_regions × n_years × n_seas × n_bins × n_sexes × n_fish_fleets]}.
#'
#' @param ret_sel_nonpar_est_bins Optional list specifying bin groupings for
#' non-parametric retained selectivity. Structure is \code{[[fleet]][[block]]}, where each
#' element is a list of bin index vectors defining grouped parameters.
#' @param fish_selex_type Character scalar specifying whether selectivity is
#'   age- or length-based. Options:
#'   \describe{
#'     \item{\code{"age"}}{Selectivity is defined over age bins.}
#'     \item{\code{"length"}}{Selectivity is defined over length bins.}
#'   }
#'   Determines the bin dimension used for all fishery selectivity functions,
#'   including parametric, time-varying, and non-parametric forms.

#' @param use_fixed_fish_sel Integer vector of length \code{n_fish_fleets}
#'   indicating whether fishery selectivity is fixed (\code{1}) or estimated
#'   (\code{0}) for each fleet.

#' @param fish_sel_input Array of fixed fishery selectivity values with
#'   dimensions:
#'   \code{[n_pop × n_regions × n_years × n_seas × n_bins × n_sexes × n_fish_fleets]}.
#'   Required when any element of \code{use_fixed_fish_sel == 1}.

#' @param fish_sel_nonpar_est_bins Optional list defining bin groupings for
#' non-parametric fishery selectivity. Structure is \code{[[fleet]][[block]]}, where each
#' element is a list of integer vectors. Each vector defines a group of bins that
#' share a single estimated selectivity parameter. Indices must correspond to the
#' bin dimension defined by \code{fish_selex_type}.

#' @param ret_selex_type Character scalar specifying whether retained selectivity
#'   is age- or length-based. Options:
#'   \describe{
#'     \item{\code{"age"}}{Selectivity is defined over age bins.}
#'     \item{\code{"length"}}{Selectivity is defined over length bins.}
#'   }
#'   Determines the bin dimension used for all retained selectivity functions,
#'   including parametric, time-varying, and non-parametric forms.
#'
#' @return The input \code{input_list} with \code{$data}, \code{$par}, and
#'   \code{$map} updated. Key additions include the parsed integer arrays for
#'   \code{cont_tv_fish_sel}, \code{fish_sel_blocks}, \code{fish_sel_model}, and
#'   \code{fish_q_blocks}; starting value arrays for all four parameter groups;
#'   and their corresponding factor maps.
#'
#'
#' @export Setup_Mod_Fishsel_and_Q
#' @importFrom stringr str_detect
#' @family Model Setup
Setup_Mod_Fishsel_and_Q <- function(input_list,

                                    # Total Selectivity
                                    cont_tv_fish_sel = paste("none_Fleet_", 1:input_list$data$n_fish_fleets, sep = ''),
                                    fish_sel_blocks = paste("none_Fleet_", 1:input_list$data$n_fish_fleets, sep = ''),
                                    fish_sel_model,
                                    Use_fish_q_prior = 0,
                                    fish_q_prior = NA,
                                    fish_q_blocks = paste("none_Fleet_", 1:input_list$data$n_fish_fleets, sep = ''),
                                    fishsel_pe_pars_spec = NULL,
                                    fish_fixed_sel_pars_spec = NULL,
                                    fish_q_spec = NULL,
                                    fish_sel_devs_spec = NULL,
                                    corr_opt_semipar = NULL,
                                    Use_fish_selex_prior = 0,
                                    fish_selex_prior = NULL,
                                    fishsel_devs_shared_bins = NULL,
                                    fish_selex_type = 'age',
                                    use_fixed_fish_sel = rep(0, input_list$data$n_fish_fleets),
                                    fish_sel_input = NULL,
                                    fish_sel_nonpar_est_bins = NULL,

                                    # Retained Selectivity
                                    cont_tv_ret_sel = paste("none_Fleet_", 1:input_list$data$n_fish_fleets, sep = ''),
                                    ret_sel_blocks = paste("none_Fleet_", 1:input_list$data$n_fish_fleets, sep = ''),
                                    ret_sel_model = paste("logist1_Fleet_", 1:input_list$data$n_fish_fleets, sep = ''),
                                    retsel_pe_pars_spec = NULL,
                                    ret_fixed_sel_pars_spec = rep("fix_ret_sel_input", input_list$data$n_fish_fleets),
                                    ret_sel_devs_spec = NULL,
                                    ret_sel_corr_opt_semipar = NULL,
                                    Use_ret_selex_prior = 0,
                                    ret_selex_prior = NULL,
                                    retsel_devs_shared_bins = NULL,
                                    ret_selex_type = 'age',
                                    use_fixed_ret_sel = rep(1, input_list$data$n_fish_fleets),
                                    ret_sel_input = array(1, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, length(input_list$data$ages), input_list$data$n_sexes, input_list$data$n_fish_fleets )),
                                    ret_sel_nonpar_est_bins = NULL,
                                    ...
                                    ) {

  messages_list <<- character(0) # string to attach to for printing messages
  starting_values <- list(...)
  if(input_list$store_config) input_list$config$Setup_Mod_Fishsel_and_Q <- mget(names(formals()))[-1]

  # Input Validation --------------------------------------------------------

  # Selectivity
  # Continuous Selectivity Deviations
  if(!is.null(fishsel_pe_pars_spec)) if(length(fishsel_pe_pars_spec) != input_list$data$n_fish_fleets) stop("fishsel_pe_pars_spec is not length n_fish_fleets")
  if(!is.null(fish_sel_devs_spec)) if(length(fish_sel_devs_spec) != input_list$data$n_fish_fleets) stop("fish_sel_devs_spec is not length n_fish_fleets")
  if(!is.null(corr_opt_semipar)) if(length(corr_opt_semipar) != input_list$data$n_fish_fleets) stop("corr_opt_semipar is not length n_fish_fleets")

  # Catchability Priors
  if(!Use_fish_q_prior %in% c(0,1)) stop("Values for Use_fish_q_prior are not valid. They are == 0 (don't use prior), or == 1 (use prior)")
  # Checking catchability priors
  if(Use_fish_q_prior == 1) {
    required_cols <- c("region", "fleet", "block", "mu", "sd")
    missing_cols <- setdiff(required_cols, names(fish_q_prior))
    if(length(missing_cols) > 0) {
      stop("fish_q_prior is missing required columns: ", paste(missing_cols, collapse = ", "))
    }
  }
  collect_message("Fishery Catchability priors are: ", ifelse(Use_fish_q_prior == 0, "Not Used", "Used"))

  # Selectivity Priors
  if(!Use_fish_selex_prior %in% c(0,1)) stop("Values for Use_fish_selex_prior are not valid. They are == 0 (don't use prior), or == 1 (use prior)")
  # Checking selectivity priors
  if(Use_fish_selex_prior == 1) {
    required_cols <- c("region", "fleet", "block", "sex", "par", "mu", "sd")
    missing_cols <- setdiff(required_cols, names(fish_selex_prior))
    if(length(missing_cols) > 0) {
      stop("fish_selex_prior is missing required columns: ", paste(missing_cols, collapse = ", "))
    }
  }
  collect_message("Fishery Selectivity priors are: ", ifelse(Use_fish_selex_prior == 0, "Not Used", "Used"))

  # Fixed selectivity options
  if(any(use_fixed_fish_sel == 1) && is.null(fish_sel_input)) stop("fish_sel_input is NULL, please provide an input array.")
  if(any(use_fixed_fish_sel == 1) && fish_selex_type == 'age') check_data_dimensions(fish_sel_input, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_ages = length(input_list$data$ages), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'fish_sel_input_age')
  if(any(use_fixed_fish_sel == 1) && fish_selex_type == 'length') check_data_dimensions(fish_sel_input, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_lens = length(input_list$data$lens), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'fish_sel_input_len')

  # Selectivity Options -----------------------------------------------------
  # Age based selectivity
  if(fish_selex_type == 'age') {
    fish_selex_type <- 0
    bins <- length(input_list$data$ages)
    collect_message("Total Fishery Selectivity is aged-based.")
  } # if age based

  # Length based selectivity
  if(fish_selex_type == 'length') {
    if(input_list$data$fit_lengths == 0) stop("Length composition data are not fit, but total selectivity is length-based. This is not allowed. Please change to a valid option (either fit lengths or use age-based selectivity).")
    fish_selex_type <- 1
    bins <- length(input_list$data$lens)
    collect_message("Total Fishery Selectivity is length-based")
  } # if length based

  # Continuous Time-Varying Selectivity Options -----------------------------
  cont_tv_fish_sel_mat <- array(NA, dim = c(input_list$data$n_regions, input_list$data$n_fish_fleets))
  cont_tv_map <- data.frame(type = c("none", "iid", "rw", "3dmarg", "3dcond", "2dar1"), num = c(0,1,2,3,4,5)) # set up values we map to

  for(i in 1:length(cont_tv_fish_sel)) {
    # Extract out components from list
    tmp <- cont_tv_fish_sel[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    cont_tv_type <- tmp_vec[1] # get continuous selex type
    fleet <- as.numeric(tmp_vec[3]) # extract fleet index

    # Validate options
    if(!fleet %in% c(1:input_list$data$n_fish_fleets)) stop("Invalid fleet specified for cont_tv_fish_sel This needs to be specified as timevarytype_Fleet_x")
    if(!cont_tv_type %in% c(cont_tv_map$type)) stop("cont_tv_fish_sel is not correctly specified. This needs to be one of these: none, iid, rw, 3dmarg, 3dcond, 2dar1 (the timevarytypes) and specified as timevarytype_Fleet_x")

    # Input options
    cont_tv_fish_sel_mat[,fleet] <- cont_tv_map$num[which(cont_tv_map$type == cont_tv_type)]
    collect_message("Continuous fishery time-varying selectivity specified as: ", cont_tv_type, " for fishery fleet ", fleet)
  }

  if(any(cont_tv_fish_sel_mat > 0) && is.null(fishsel_pe_pars_spec) && is.null(fish_sel_devs_spec)) stop("Continuous time-varying selectivity specified, but fishsel_pe_pars_spec and/or fish_sel_devs_spec is NULL (i.e., not specified)!")

  # Blocked Time-Varying Selectivity Options --------------------------------
  fish_sel_blocks_arr <- array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets))
  for(i in 1:length(fish_sel_blocks)) {

    # Extract out components from list
    tmp <- fish_sel_blocks[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))

    # Validate options
    if(!tmp_vec[1] %in% c("none", "Block")) stop("Fishery Selectivity Blocks not correctly specified. This should be either none_Fleet_x or Block_x_Year_x-y_Fleet_x")

    # extract out fleets if constant
    if(tmp_vec[1] == "none") {
      fleet <- as.numeric(tmp_vec[3]) # get fleet number
      fish_sel_blocks_arr[,,fleet] <- 1 # input only 1 fishery time block
    }

    if(tmp_vec[1] == "Block") {

      block_val <- as.numeric(tmp_vec[2]) # get block value
      fleet <- as.numeric(tmp_vec[6]) # extract fleet index

      # get year ranges
      if(!str_detect(tmp, "terminal")) { # if not terminal year
        year_range <- as.numeric(unlist(strsplit(tmp_vec[4], "-")))
        years <- year_range[1]:year_range[2] # get sequence of years
      } else { # if terminal year
        year_range <- unlist(strsplit(tmp_vec[4], '-'))[1] # get year range
        years <- as.numeric(year_range):length(input_list$data$years) # get sequence of years
      }

      fish_sel_blocks_arr[,years,fleet] <- block_val
    }

  }

  if(any(is.na(fish_sel_blocks_arr))) stop("Fishery Selectivtiy Blocks are returning an NA. Did you forget to specify the year range of fish_sel_blocks?")
  for(f in 1:input_list$data$n_fish_fleets) collect_message(paste("Fishery Selectivity Time Blocks for fishery", f, "is specified at:", length(unique(fish_sel_blocks_arr[,,f]))))

  # Selectivity Functional Forms --------------------------------------------
  sel_map <- data.frame(sel = c('logist1', "gamma", "exponential", "logist2", "dbnrml", 'nonpar', 'asymplogist1', "asymplogist2", "bicubic"), num = c(0,1,2,3,4,5,6,7,8)) # set up values we can map to
  fish_sel_model_arr <- array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets))
  fish_sel_bicubic_binnodes_arr <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets)) # number of bin nodes, only set where fish_sel_model == 8 (bicubic)
  fish_sel_bicubic_yrnodes_arr <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets)) # number of year nodes, only set where fish_sel_model == 8 (bicubic)
  fish_sel_bicubic_selstyr_arr <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets)) # calendar year the bicubic surface is actually fit from (0 = block's own start year, i.e. no offset); years within the block before this are edge-held at this year's fitted curve
  fish_sel_bicubic_nselbins_arr <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets)) # number of bins (starting from the first) the bicubic surface is actually fit over (0 = all bins, i.e. no truncation); bins beyond this are held flat at the last fitted bin's value

  for(i in 1:length(fish_sel_model)) {

    # Extract out fishery selectivity components from vector
    tmp_sel_form <- fish_sel_model[i]
    tmp_sel_form_vec <- unlist(strsplit(tmp_sel_form, "_")) # split string
    sel_form <- tmp_sel_form_vec[1] # get selectivity type

    if(sel_form == "bicubic") {
      # bicubic spline: bicubic_Bin_<n_bin_nodes>_Yr_<n_yr_nodes>_Fleet_<f>[_Block_<b>][_SelStyr_<year>][_NSelBins_<n>]
      bin_pos <- which(tmp_sel_form_vec == "Bin")
      yr_pos <- which(tmp_sel_form_vec == "Yr")
      fleet_pos <- which(tmp_sel_form_vec == "Fleet")
      block_pos <- which(tmp_sel_form_vec == "Block")
      selstyr_pos <- which(tmp_sel_form_vec == "SelStyr")
      nselbins_pos <- which(tmp_sel_form_vec == "NSelBins")
      if(length(bin_pos) != 1 || length(yr_pos) != 1 || length(fleet_pos) != 1)
        stop("fish_sel_model 'bicubic' entries must be specified as bicubic_Bin_<n_bin_nodes>_Yr_<n_yr_nodes>_Fleet_<f> or bicubic_Bin_<n_bin_nodes>_Yr_<n_yr_nodes>_Block_<b>_Fleet_<f>, optionally with _SelStyr_<year> and/or _NSelBins_<n>")
      tmp_n_bin_nodes <- suppressWarnings(as.numeric(tmp_sel_form_vec[bin_pos + 1]))
      tmp_n_yr_nodes <- suppressWarnings(as.numeric(tmp_sel_form_vec[yr_pos + 1]))
      tmp_fleet <- suppressWarnings(as.numeric(tmp_sel_form_vec[fleet_pos + 1]))
      tmp_block <- if(length(block_pos) == 1) suppressWarnings(as.numeric(tmp_sel_form_vec[block_pos + 1])) else NULL
      tmp_selstyr <- if(length(selstyr_pos) == 1) suppressWarnings(as.numeric(tmp_sel_form_vec[selstyr_pos + 1])) else 0
      tmp_nselbins <- if(length(nselbins_pos) == 1) suppressWarnings(as.numeric(tmp_sel_form_vec[nselbins_pos + 1])) else 0
      if(is.na(tmp_n_bin_nodes) || tmp_n_bin_nodes < 2) stop("bicubic fish_sel_model requires at least 2 bin nodes (n_bin_nodes >= 2)")
      if(is.na(tmp_n_yr_nodes) || tmp_n_yr_nodes < 1) stop("bicubic fish_sel_model requires at least 1 year node (n_yr_nodes >= 1). Use n_yr_nodes == 1 for a time-invariant bin-only spline.")
      if(length(selstyr_pos) == 1 && (is.na(tmp_selstyr) || !tmp_selstyr %in% input_list$data$years)) stop("bicubic fish_sel_model SelStyr must be a calendar year within the modeled years")
      if(length(nselbins_pos) == 1 && (is.na(tmp_nselbins) || tmp_nselbins < 2 || tmp_nselbins > bins)) stop("bicubic fish_sel_model NSelBins must be an integer between 2 and the total number of bins (ages or lengths)")
    } else {
      # get fleet index
      tmp_fleet <- if(length(tmp_sel_form_vec) == 3) as.numeric(tmp_sel_form_vec[3]) else as.numeric(tmp_sel_form_vec[5]) # fleet index changes if block is included in character vector
      # get block index
      tmp_block <- if(length(tmp_sel_form_vec) == 5) as.numeric(tmp_sel_form_vec[3]) else NULL
    }

    # validate options
    if(!sel_form %in% c(sel_map$sel)) stop("fish_sel_model is not correctly specified. This needs to be one of these: logist1, gamma, exponential, logist2, dbnrml, nonpar, asymplogist1, asymplogist2, bicubic (the seltypes) and specified as seltype_Fleet_x")
    if(!tmp_fleet %in% c(1:input_list$data$n_fish_fleets)) stop("Invalid fleet specified for fish_sel_model This needs to be specified as seltype_Fleet_x or seltype_Fleet_x_Block_x (if blocks are specified to change for a fleet)")

    # Input options
    if(is.null(tmp_block)) {
      fish_sel_model_arr[,,tmp_fleet] <- sel_map$num[which(sel_map$sel == sel_form)] # same selectivity form across blocks
      if(sel_form == "bicubic") {
        fish_sel_bicubic_binnodes_arr[,,tmp_fleet] <- tmp_n_bin_nodes
        fish_sel_bicubic_yrnodes_arr[,,tmp_fleet] <- tmp_n_yr_nodes
        fish_sel_bicubic_selstyr_arr[,,tmp_fleet] <- tmp_selstyr
        fish_sel_bicubic_nselbins_arr[,,tmp_fleet] <- tmp_nselbins
      }
    } else {
      fish_sel_model_arr[,which(fish_sel_blocks_arr[,,tmp_fleet] == tmp_block),tmp_fleet] <- sel_map$num[which(sel_map$sel == sel_form)]
      if(sel_form == "bicubic") {
        fish_sel_bicubic_binnodes_arr[,which(fish_sel_blocks_arr[,,tmp_fleet] == tmp_block),tmp_fleet] <- tmp_n_bin_nodes
        fish_sel_bicubic_yrnodes_arr[,which(fish_sel_blocks_arr[,,tmp_fleet] == tmp_block),tmp_fleet] <- tmp_n_yr_nodes
        fish_sel_bicubic_selstyr_arr[,which(fish_sel_blocks_arr[,,tmp_fleet] == tmp_block),tmp_fleet] <- tmp_selstyr
        fish_sel_bicubic_nselbins_arr[,which(fish_sel_blocks_arr[,,tmp_fleet] == tmp_block),tmp_fleet] <- tmp_nselbins
      }
    }
    rm(tmp_block) # remove tmp block to start next loop
    collect_message("Fishery selectivity functional form specified as:", sel_form, " for fishery fleet ", tmp_fleet)
  }

  # Validate that blocks and continuous time-variation aren't both specified for same fleet
  for(f in 1:input_list$data$n_fish_fleets) {
    has_blocks <- length(unique(fish_sel_blocks_arr[1,,f])) > 1
    has_cont_tv <- cont_tv_fish_sel_mat[1,f] != 0  # 0 = "none"
    if(has_blocks && has_cont_tv) {
      stop("Fleet ", f, " has both selectivity blocks and continuous time-varying selectivity specified. ",
           "These are mutually exclusive - choose one approach to time-variation.")
    }
  }

  # Blocked Catchability Options --------------------------------------------
  fish_q_blocks_arr <- array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets))
  for(i in 1:length(fish_q_blocks)) {
    # Extract out components from list
    tmp <- fish_q_blocks[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))

    # Validate options
    if(!tmp_vec[1] %in% c("none", "Block")) stop("Fishery Catchability Blocks not correctly specified. This should be either none_Fleet_x or Block_x_Year_x-y_Fleet_x")

    # extract out fleets if constant
    if(tmp_vec[1] == "none") {
      fleet <- as.numeric(tmp_vec[3]) # get fleet number
      fish_q_blocks_arr[,,fleet] <- 1 # input only 1 fishery catchability time block
    }

    if(tmp_vec[1] == "Block") {

      block_val <- as.numeric(tmp_vec[2]) # get block value
      fleet <- as.numeric(tmp_vec[6]) # get fleet number

      # get year ranges
      if(!str_detect(tmp, "terminal")) { # if not terminal year
        year_range <- as.numeric(unlist(strsplit(tmp_vec[4], "-")))
        years <- year_range[1]:year_range[2] # get sequence of years
      } else { # if terminal year
        year_range <- unlist(strsplit(tmp_vec[4], '-'))[1] # get year range
        years <- as.numeric(year_range):length(input_list$data$years) # get sequence of years
      }

      fish_q_blocks_arr[,years,fleet] <- block_val # input catchability time block
    }
  }

  if(any(is.na(fish_q_blocks))) stop("Fishery Catchability Blocks are returning an NA. Did you forget to specify the year range of fish_q_blocks?")
  for(f in 1:input_list$data$n_fish_fleets) collect_message(paste("Fishery Catchability Time Blocks for fishery", f, "is specified at:", length(unique(fish_q_blocks_arr[,,f]))))

  # Populate Data List ------------------------------------------------------

  input_list$data$cont_tv_fish_sel <- cont_tv_fish_sel_mat
  input_list$data$fish_sel_blocks <- fish_sel_blocks_arr
  input_list$data$fish_sel_model <- fish_sel_model_arr
  input_list$data$fish_sel_bicubic_binnodes <- fish_sel_bicubic_binnodes_arr
  input_list$data$fish_sel_bicubic_yrnodes <- fish_sel_bicubic_yrnodes_arr
  input_list$data$fish_sel_bicubic_selstyr <- fish_sel_bicubic_selstyr_arr
  input_list$data$fish_sel_bicubic_nselbins <- fish_sel_bicubic_nselbins_arr
  input_list$data$fish_q_blocks <- fish_q_blocks_arr
  input_list$data$fish_q_prior <- fish_q_prior
  input_list$data$Use_fish_q_prior <- Use_fish_q_prior
  input_list$data$Use_fish_selex_prior <- Use_fish_selex_prior
  input_list$data$fish_selex_prior <- fish_selex_prior
  input_list$data$fish_selex_type <- fish_selex_type
  input_list$data$use_fixed_fish_sel <- use_fixed_fish_sel
  input_list$data$fish_sel_input <- fish_sel_input
  input_list$data$fishsel_devs_min_shared_bins <- if(!is.null(fishsel_devs_shared_bins)) unlist(lapply(fishsel_devs_shared_bins, min)) else 1:length(input_list$data$ages)

  # Populate Parameter List -------------------------------------------------

  # Figure out number of selectivity parameters for a given functional form
  unique_fishsel_vals <- unique(as.vector(input_list$data$fish_sel_model))
  sel_pars_vec <- vector() # create empty vector to populate

  for(i in 1:length(unique_fishsel_vals)) {
    if(unique_fishsel_vals[i] %in% c(2)) sel_pars_vec[i] <- 1 # exponential
    if(unique_fishsel_vals[i] %in% c(0,1,3)) sel_pars_vec[i] <- 2 # logistic or gamma
    if(unique_fishsel_vals[i] %in% c(4)) sel_pars_vec[i] <- 6 # double normal
    if(unique_fishsel_vals[i] == 5) sel_pars_vec[i] <- bins # non-parametric selex
    if(unique_fishsel_vals[i] %in% c(6,7)) sel_pars_vec[i] <- 3 # logistic selex w/ asymptote
    if(unique_fishsel_vals[i] == 8) sel_pars_vec[i] <- max(input_list$data$fish_sel_bicubic_binnodes * input_list$data$fish_sel_bicubic_yrnodes) # bicubic: flattened bin-node x year-node grid
  } # end i loop

  # figure out maximum number of fishery selectivity blocks for a given reigon and fleet
  max_fishsel_blks <- max(apply(input_list$data$fish_sel_blocks, c(1,3), FUN = function(x) length(unique(x))))

  # Bicubic spline interpolation weight matrices (bin node x year node grid), built here so they can be
  # threaded through SPoRC_rtmb.R alongside the flattened node parameters (see Get_Selex, Selex_Model == 8).
  # Padded with zeros to a common width across regions/blocks/fleets; padding is harmless because unused
  # (zero-weight) columns/rows never contribute to the resulting selectivity (see Get_Selex documentation).
  has_bicubic_fish_sel <- any(input_list$data$fish_sel_model == 8)
  max_bin_nodes_bicubic <- if(has_bicubic_fish_sel) max(input_list$data$fish_sel_bicubic_binnodes) else 1
  max_yr_nodes_bicubic <- if(has_bicubic_fish_sel) max(input_list$data$fish_sel_bicubic_yrnodes) else 1
  n_yrs_total_bicubic <- length(input_list$data$years) + input_list$data$n_proj_yrs_devs

  fish_sel_bicubic_Wbin <- array(0, dim = c(input_list$data$n_regions, bins, max_bin_nodes_bicubic, max_fishsel_blks, input_list$data$n_fish_fleets))
  fish_sel_bicubic_Wyr <- array(0, dim = c(input_list$data$n_regions, n_yrs_total_bicubic, max_yr_nodes_bicubic, max_fishsel_blks, input_list$data$n_fish_fleets))

  if(has_bicubic_fish_sel) {
    for(f in 1:input_list$data$n_fish_fleets) {
      for(r in 1:input_list$data$n_regions) {

        fishsel_blocks_tmp <- unique(as.vector(input_list$data$fish_sel_blocks[r,,f]))

        for(b in 1:length(fishsel_blocks_tmp)) {

          block_years <- which(input_list$data$fish_sel_blocks[r,,f] == fishsel_blocks_tmp[b])
          if(unique(input_list$data$fish_sel_model[r, block_years, f]) != 8) next # only bicubic blocks need weight matrices

          n_bin_nodes_this <- unique(input_list$data$fish_sel_bicubic_binnodes[r, block_years, f])
          n_yr_nodes_this <- unique(input_list$data$fish_sel_bicubic_yrnodes[r, block_years, f])

          # Age/length dimension: nodes evenly spaced over [0,1]. By default (NSelBins unset, i.e. 0)
          # the spline is evaluated over all bins, as before. When NSelBins is set , the spline surface is only actually fit over the first NSelBins bins;
          # bins beyond that are edge-held at the last fitted bin's weights ("plateau")
          nselbins_this <- unique(input_list$data$fish_sel_bicubic_nselbins[r, block_years, f])
          n_fit_bins <- if(nselbins_this == 0) bins else nselbins_this

          bin_nodes_scaled <- seq(0, 1, length.out = n_bin_nodes_this)
          fit_bin_scaled <- seq(0, 1, length.out = n_fit_bins)
          Wbin_fit <- Get_Natural_Cubic_Spline_Weights(bin_nodes_scaled, fit_bin_scaled)

          Wbin_this <- matrix(0, nrow = bins, ncol = n_bin_nodes_this)
          Wbin_this[1:n_fit_bins, ] <- Wbin_fit
          if(n_fit_bins < bins) Wbin_this[(n_fit_bins + 1):bins, ] <- matrix(Wbin_fit[nrow(Wbin_fit), ], nrow = bins - n_fit_bins, ncol = n_bin_nodes_this, byrow = TRUE)

          fish_sel_bicubic_Wbin[r, , 1:n_bin_nodes_this, b, f] <- Wbin_this

          # Year dimension: nodes evenly spaced over the block's own contiguous fit range. By default
          # (SelStyr unset, i.e. 0) the fit range is the whole block, as before. When SelStyr is set, only years from SelStyr through the block's end are
          # actually spline-fit; years within the block before SelStyr are edge-held at the SelStyr
          # row's weights ("previous years are filled"). Years outside the block entirely (before it,
          # after it, and any projection years, since projections reuse the terminal modeled year's
          # block) hold the boundary node weights constant, which for a spline evaluated exactly at
          # its first/last node reduces to full weight on that node.
          selstyr_this <- unique(input_list$data$fish_sel_bicubic_selstyr[r, block_years, f])
          selstyr_idx <- if(selstyr_this == 0) min(block_years) else which(input_list$data$years == selstyr_this)
          fit_years <- block_years[block_years >= selstyr_idx]
          pre_fit_years <- block_years[block_years < selstyr_idx]

          yr_nodes_scaled <- seq(0, 1, length.out = n_yr_nodes_this)
          fit_yr_scaled <- seq(0, 1, length.out = length(fit_years))
          Wyr_block <- Get_Natural_Cubic_Spline_Weights(yr_nodes_scaled, fit_yr_scaled)

          Wyr_this <- matrix(0, nrow = n_yrs_total_bicubic, ncol = n_yr_nodes_this)
          Wyr_this[fit_years, ] <- Wyr_block
          if(length(pre_fit_years) > 0) Wyr_this[pre_fit_years, ] <- matrix(Wyr_block[1, ], nrow = length(pre_fit_years), ncol = n_yr_nodes_this, byrow = TRUE)
          if(min(block_years) > 1) Wyr_this[1:(min(block_years) - 1), ] <- matrix(Wyr_block[1, ], nrow = min(block_years) - 1, ncol = n_yr_nodes_this, byrow = TRUE)
          if(max(block_years) < n_yrs_total_bicubic) Wyr_this[(max(block_years) + 1):n_yrs_total_bicubic, ] <- matrix(Wyr_block[nrow(Wyr_block), ], nrow = n_yrs_total_bicubic - max(block_years), ncol = n_yr_nodes_this, byrow = TRUE)

          fish_sel_bicubic_Wyr[r, , 1:n_yr_nodes_this, b, f] <- Wyr_this

        } # end b loop
      } # end r loop
    } # end f loop
  } # end if has_bicubic_fish_sel

  input_list$data$fish_sel_bicubic_Wbin <- fish_sel_bicubic_Wbin
  input_list$data$fish_sel_bicubic_Wyr <- fish_sel_bicubic_Wyr

  # maximum number of selectivity parameters across all forms
  max_fishsel_pars <- max(sel_pars_vec)
  if("fish_fixed_sel_pars" %in% names(starting_values)) input_list$par$fish_fixed_sel_pars <- starting_values$fish_fixed_sel_pars
  else input_list$par$fish_fixed_sel_pars <- array(0, dim = c(input_list$data$n_regions, max_fishsel_pars, max_fishsel_blks, input_list$data$n_sexes, input_list$data$n_fish_fleets))

  # Fishery catchability
  max_fishq_blks <- max(apply(input_list$data$fish_q_blocks, c(1,3), FUN = function(x) length(unique(x)))) # figure out maximum number of fishery catchability blocks for a given reigon and fleet
  if("ln_fish_q" %in% names(starting_values)) input_list$par$ln_fish_q <- starting_values$ln_fish_q
  else input_list$par$ln_fish_q <- array(0, dim = c(input_list$data$n_regions, max_fishq_blks, input_list$data$n_fish_fleets))

  # Fishery selectivity process error parameters
  if("fishsel_pe_pars" %in% names(starting_values)) input_list$par$fishsel_pe_pars <- starting_values$fishsel_pe_pars
  else input_list$par$fishsel_pe_pars <- array(0, dim = c(input_list$data$n_regions, max(max_fishsel_pars, 4), input_list$data$n_sexes, input_list$data$n_fish_fleets)) # dimensioned 4 as the max number of pars for process errors (e.g., sigmas), and then just map off if not using

  # Fishery selectivity deviations
  if("ln_fishsel_devs" %in% names(starting_values)) input_list$par$ln_fishsel_devs <- starting_values$ln_fishsel_devs
  else input_list$par$ln_fishsel_devs <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years) + input_list$data$n_proj_yrs_devs, bins, input_list$data$n_sexes, input_list$data$n_fish_fleets))


  # Mapping Options ---------------------------------------------------------
  input_list <- do_fish_fixed_sel_pars_mapping(input_list, fish_fixed_sel_pars_spec, bins, fish_sel_nonpar_est_bins)
  input_list <- do_fish_q_mapping(input_list, fish_q_spec)
  input_list <- do_fishsel_pe_pars_mapping(input_list, fishsel_pe_pars_spec, corr_opt_semipar, bins)
  input_list <- do_fishsel_devs_mapping(input_list, fish_sel_devs_spec, fishsel_devs_shared_bins, bins)


  # Retained Selectivity ---------------------------------------------
  input_list <- Setup_Mod_Retsel(input_list,
                                 cont_tv_ret_sel = cont_tv_ret_sel,
                                 ret_sel_blocks = ret_sel_blocks,
                                 ret_sel_model = ret_sel_model,
                                 retsel_pe_pars_spec = retsel_pe_pars_spec,
                                 ret_fixed_sel_pars_spec = ret_fixed_sel_pars_spec,
                                 ret_sel_devs_spec = ret_sel_devs_spec,
                                 ret_sel_corr_opt_semipar = ret_sel_corr_opt_semipar,
                                 Use_ret_selex_prior = Use_ret_selex_prior,
                                 ret_selex_prior = ret_selex_prior,
                                 retsel_devs_shared_bins = retsel_devs_shared_bins,
                                 ret_selex_type = ret_selex_type,
                                 use_fixed_ret_sel = use_fixed_ret_sel,
                                 ret_sel_input = ret_sel_input,
                                 ret_sel_nonpar_est_bins = ret_sel_nonpar_est_bins,
                                 ...)

  # Print Messages ----------------------------------------------------------
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)
}



