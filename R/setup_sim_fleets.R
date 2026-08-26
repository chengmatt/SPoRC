# Stage 1 of 3: model setup
#
# Operating model fleet inputs. Setup_Sim_Fishing and Setup_Sim_Survey each
# populate sim_list with a whole fleet's worth of settings at once: mortality or
# catchability, selectivity, timing, index type and the composition likelihood
# choices with their overdispersion and correlation parameters. They live
# together because each one spans every topic that the setup_fishery_*.R and
# setup_survey_*.R files split apart on the estimation side.

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
#' @param FishIdx_LikeType Character or numeric vector, length `n_fish_fleets`.
#'   Error structure each fleet's index is drawn under: \code{"lognormal"} (0),
#'   \code{"normal"} (1), or \code{"mvn"} (2), matching the estimation model's
#'   \code{FishIdx_LikeType}. An mvn fleet draws from \code{FishIdx_Cov} through
#'   a common-factor decomposition (see \code{\link{cov_to_factor}}) instead of
#'   \code{ObsFishIdx_SE}, and its population-specific stream stays lognormal.
#'   Default: lognormal for every fleet.
#' @param FishIdx_Cov List with one element per fishery fleet holding the fixed
#'   covariance over that fleet's fitted index observations, ordered by scanning
#'   \code{UseFishIdx} in array order (region fastest, then year, then season).
#'   Required for mvn fleets. Default: \code{NULL}.
#' @param UseFishIdx Numeric array \code{[n_regions x n_yrs x n_seas x n_fish_fleets]}
#'   of fit flags from the estimation model, used to position each simulated cell
#'   in the covariance. Its year dimension may be shorter than the simulation,
#'   in which case later years draw with the mean factor scale and loading.
#'   Required for mvn fleets. Default: \code{NULL}.
#'
#' @param t_fish Numeric array \code{[n_regions x n_seas x n_fish_fleets]} giving
#'   the fishery index timing, the fraction of the season elapsed when the index
#'   is observed. Numbers at age are decayed by \code{exp(-t_fish * ZAA)} before
#'   the index is formed, matching \code{t_srv} for surveys and the estimation
#'   model's own \code{t_fish}. Defaults to \code{0} (start of season).
#' @param comp_fishage_like Numeric vector. Likelihood for age composition
#'   (0 = Multinomial, 1 = Dirichlet-Multinomial, 2-4 = Logistic-Normal variants),
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
#'   (0 = Multinomial, 1 = Dirichlet-Multinomial, 2-4 = Logistic-Normal variants),
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
#'   (0 = Multinomial, 1 = Dirichlet-Multinomial, 2-4 = Logistic-Normal variants),
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
#'   (0 = Multinomial, 1 = Dirichlet-Multinomial, 2-4 = Logistic-Normal variants),
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
#'   (0 = Multinomial, 1 = Dirichlet-Multinomial, 2-4 = Logistic-Normal variants, 999 = none),
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
#'   (0 = Multinomial, 1 = Dirichlet-Multinomial, 2-4 = Logistic-Normal variants, 999 = none),
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
#'   (0 = Multinomial, 1 = Dirichlet-Multinomial, 2-4 = Logistic-Normal variants, 999 = none),
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
#'   (0 = Multinomial, 1 = Dirichlet-Multinomial, 2-4 = Logistic-Normal variants, 999 = none),
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
                              ln_sigmaCAA = array(log(0.2), dim = c(sim_list$n_ages, sim_list$n_fish_fleets)),
                              ln_sigmaDAA = array(log(0.2), dim = c(sim_list$n_ages, sim_list$n_fish_fleets)),
                              ln_sigmaFishIdxAA = array(log(0.2), dim = c(sim_list$n_ages, sim_list$n_fish_fleets)),
                              UseCatchAA = array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_fish_fleets)),
                              UseDiscardAA = array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_fish_fleets)),
                              UseFishIdxAA = array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_fish_fleets)),
                              use_catch_aa = rep(0, sim_list$n_fish_fleets),
                              use_discard_aa = rep(0, sim_list$n_fish_fleets),
                              use_fish_idx_aa = rep(0, sim_list$n_fish_fleets),
                              catch_units = array(1, dim = c(sim_list$n_fish_fleets)),
                              init_F_val = array(0, dim = c(sim_list$n_regions, sim_list$n_seas, sim_list$n_fish_fleets)),
                              Fmort_input = array(0.1, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets, sim_list$n_sims)),
                              fish_sel_input,
                              fish_q_input = array(1, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_fish_fleets, sim_list$n_sims)),
                              ObsFishIdx_SE = array(0.2, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets)),
                              ObsFishIdx_pop_SE = array(0.2, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets)),
                              fish_idx_type = array(1, dim = c(sim_list$n_regions, sim_list$n_fish_fleets)),
                              FishIdx_LikeType = rep(0, sim_list$n_fish_fleets),
                              FishIdx_Cov = NULL,
                              UseFishIdx = NULL,
                              t_fish = array(0, dim = c(sim_list$n_regions, sim_list$n_seas, sim_list$n_fish_fleets)),

                              # Conditional age-at-length. Off unless an ISS array is supplied; the
                              # thetas default to log(1) and are only read by the Dirichlet-multinomial
                              comp_fish_caal_like = rep(999, sim_list$n_fish_fleets),
                              ISS_Fish_caal = NULL,
                              ln_Fish_caal_theta = NULL,
                              ln_Fish_caal_theta_agg = NULL,
                              Fish_caal_Type = array(999, dim = c(sim_list$n_yrs, sim_list$n_fish_fleets)),

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
  FishIdx_LikeType <- convert_to_numeric(FishIdx_LikeType, list(lognormal = 0, normal = 1, mvn = 2))
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
  check_sim_dimensions(FishIdx_LikeType, n_fish_fleets = sim_list$n_fish_fleets, what = "FishIdx_LikeType")

  # Multivariate normal index fleets draw from the supplied covariance rather than
  # the SE array, so the covariance is validated and factor-decomposed once here.
  fish_idx_mvn <- NULL
  if(any(FishIdx_LikeType == 2)) {
    if(is.null(UseFishIdx)) stop("UseFishIdx must be supplied when any FishIdx_LikeType is mvn, to position each observation in the covariance.")
    if(length(dim(UseFishIdx)) != 4 || any(dim(UseFishIdx)[c(1,3,4)] != c(sim_list$n_regions, sim_list$n_seas, sim_list$n_fish_fleets)) || dim(UseFishIdx)[2] > sim_list$n_yrs)
      stop("UseFishIdx must be an n_regions x (at most n_yrs) x n_seas x n_fish_fleets array.")
    fish_idx_mvn <- build_idx_factor(FishIdx_Cov, FishIdx_LikeType, UseFishIdx, sim_list$n_fish_fleets, "FishIdx_Cov")
  }

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
  sim_list$ln_sigmaC_pop <- ln_sigmaC_pop
  sim_list$ln_sigmaCAA <- ln_sigmaCAA
  sim_list$ln_sigmaDAA <- ln_sigmaDAA
  sim_list$ln_sigmaFishIdxAA <- ln_sigmaFishIdxAA
  sim_list$UseCatchAA <- UseCatchAA
  sim_list$UseDiscardAA <- UseDiscardAA
  sim_list$UseFishIdxAA <- UseFishIdxAA
  sim_list$use_catch_aa <- use_catch_aa
  sim_list$use_discard_aa <- use_discard_aa
  sim_list$use_fish_idx_aa <- use_fish_idx_aa # observation sd for pop-specific catch
  sim_list$init_F <- init_F_val # initial F value
  sim_list$fish_sel <- fish_sel_input # fishery selectivity
  sim_list$fish_q <- fish_q_input # fishery catchability
  sim_list$ObsFishIdx_SE <- ObsFishIdx_SE # fishery index SE
  sim_list$ObsFishIdx_pop_SE <- ObsFishIdx_pop_SE # fishery index SE pop-specific
  sim_list$fish_idx_type <- fish_idx_type # fishery index type
  sim_list$FishIdx_LikeType <- FishIdx_LikeType # fishery index error structure
  if(!is.null(fish_idx_mvn)) {
    sim_list$fish_idx_mvn <- fish_idx_mvn # factor parameters for mvn index fleets
    sim_list$fish_idx_u <- matrix(NA_real_, sim_list$n_fish_fleets, sim_list$n_sims) # shared factor draw, filled per fleet and replicate
  }
  sim_list$t_fish <- t_fish # fishery index timing within the season

  # Fishery age compositions
  sim_list$comp_fishage_like <- comp_fishage_like
  sim_list$ISS_FishAgeComps <- ISS_FishAgeComps
  sim_list$ln_FishAge_theta <- ln_FishAge_theta
  sim_list$ln_FishAge_theta_agg <- ln_FishAge_theta_agg
  sim_list$FishAge_corr_pars_agg <- FishAge_corr_pars_agg
  sim_list$FishAge_corr_pars <- FishAge_corr_pars
  sim_list$FishAgeComps_Type <- FishAgeComps_Type

  # Fishery conditional age-at-length
  comp_fish_caal_like <- convert_to_numeric(comp_fish_caal_like, list(Multinomial = 0, `Dirichlet-Multinomial` = 1, none = 999))
  Fish_caal_Type <- convert_to_numeric(Fish_caal_Type, list(agg = 0, spltRspltS = 1, spltRjntS = 2, none = 999))
  if(is.null(ln_Fish_caal_theta)) ln_Fish_caal_theta <- array(log(1), dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets))
  if(is.null(ln_Fish_caal_theta_agg)) ln_Fish_caal_theta_agg <- rep(log(1), sim_list$n_fish_fleets)
  sim_list$do_fish_caal <- !is.null(ISS_Fish_caal) && any(comp_fish_caal_like != 999)
  if(sim_list$do_fish_caal) {
    if(is.null(sim_list$n_lens)) stop("ISS_Fish_caal was supplied, but the simulation has no length bins (n_lens is NULL)")
    if(length(dim(ISS_Fish_caal)) != 7 || !all(dim(ISS_Fish_caal) == c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_lens, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims)))
      stop("Dimensions of ISS_Fish_caal are not correct. Should be n_regions, n_years, n_seas, n_lens, n_sexes, n_fish_fleets, and n_sims")
  }
  sim_list$comp_fish_caal_like <- comp_fish_caal_like
  sim_list$ISS_Fish_caal <- ISS_Fish_caal
  sim_list$ln_Fish_caal_theta <- ln_Fish_caal_theta
  sim_list$ln_Fish_caal_theta_agg <- ln_Fish_caal_theta_agg
  sim_list$Fish_caal_Type <- Fish_caal_Type

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

#' Set up survey parameterisation for the operating model simulation
#'
#' Populates \code{sim_list} with all survey-related inputs needed by the
#' operating model: catchability, selectivity, survey timing, index type,
#' and age/length composition likelihood settings including overdispersion
#' and correlation parameters. Must be called after \code{\link{Setup_Sim_Dim}}.
#'
#' @param sim_list Simulation list returned by \code{\link{Setup_Sim_Dim}}.
#' @param srv_sel_input Survey selectivity array
#'   \code{[n_pop x n_regions x n_yrs x n_seas × n_ages × n_sexes × n_srv_fleets × n_sims]}.
#'   No default; must be provided.
#' @param srv_q_input Survey catchability array
#'   \code{[n_regions × n_yrs × n_srv_fleets × n_sims]}. Default: 1 for all cells.
#' @param ObsSrvIdx_SE Lognormal observation error SD for survey index,
#'   array \code{[n_regions × n_yrs × n_seas × n_srv_fleets]}. Default: 0.2.
#' @param ObsSrvIdx_pop_SE As above, but for population-specific indices,
#'   array \code{[n_pop × n_regions × n_yrs × n_seas × n_srv_fleets]}.
#' @param t_srv Survey timing as fraction of year or season, array
#'   \code{[n_regions × n_seas × n_srv_fleets]}. Default: 1.
#' @param srv_idx_type Integer vector \code{[n_srv_fleets]} specifying survey
#'   index type. Default: all 1 (biomass). Options: 0/“abd” (abundance),
#'   1/“biom” (biomass).
#' @param SrvIdx_LikeType Character or numeric vector, length `n_srv_fleets`.
#'   Error structure each fleet's index is drawn under: \code{"lognormal"} (0),
#'   \code{"normal"} (1), or \code{"mvn"} (2), matching the estimation model's
#'   \code{SrvIdx_LikeType}. An mvn fleet draws from \code{SrvIdx_Cov} through
#'   a common-factor decomposition (see \code{\link{cov_to_factor}}) instead of
#'   \code{ObsSrvIdx_SE}, and its population-specific stream stays lognormal.
#'   Default: lognormal for every fleet.
#' @param SrvIdx_Cov List with one element per survey fleet holding the fixed
#'   covariance over that fleet's fitted index observations, ordered by scanning
#'   \code{UseSrvIdx} in array order (region fastest, then year, then season).
#'   Required for mvn fleets. Default: \code{NULL}.
#' @param UseSrvIdx Numeric array \code{[n_regions x n_yrs x n_seas x n_srv_fleets]}
#'   of fit flags from the estimation model, used to position each simulated cell
#'   in the covariance. Its year dimension may be shorter than the simulation,
#'   in which case later years draw with the mean factor scale and loading.
#'   Required for mvn fleets. Default: \code{NULL}.
#' @param comp_srvage_like Integer or character vector \code{[n_srv_fleets]}
#'   specifying likelihood for survey age compositions. Default: all 0
#'   (multinomial). Options: 0/“Multinomial”, 1/“Dirichlet-Multinomial”,
#'   2/“iid-Logistic-Normal”, 3/“1d-Logistic-Normal”, 4/“2d-Logistic-Normal”.
#' @param ISS_SrvAgeComps Array \code{[n_regions × n_yrs × n_seas × n_sexes × n_srv_fleets × n_sims]}
#'   of sample sizes or overdispersion for survey age compositions. Default: 100.
#' @param ln_SrvAge_theta Log-scale overdispersion array
#'   \code{[n_regions × n_sexes × n_srv_fleets]}. Used for likelihoods 1-4.
#'   Default: log(1).
#' @param ln_SrvAge_theta_agg Log-scale overdispersion for aggregated survey
#'   age compositions, vector \code{[n_srv_fleets]}. Default: log(1).
#' @param SrvAge_corr_pars Correlation parameters array
#'   \code{[n_regions × n_sexes × n_srv_fleets × 2]} (age AR1, sex). Only for
#'   likelihoods 3-4. Default: 0.01.
#' @param SrvAge_corr_pars_agg Vector \code{[n_srv_fleets]} for aggregated
#'   survey age correlations. Only for likelihood 3. Default: 0.01.
#' @param SrvAgeComps_Type Array \code{[n_yrs × n_srv_fleets]} specifying
#'   composition structure. Default: 2 (split by region, joint sexes).
#'   Options: 0/“agg”, 1/“spltRspltS”, 2/“spltRjntS”, 999/“none”.
#' @param comp_srvlen_like Integer or character vector \code{[n_srv_fleets]}
#'   specifying likelihood for survey length compositions. Default: all 0.
#' @param ISS_SrvLenComps Array \code{[n_regions × n_yrs × n_seas × n_sexes × n_srv_fleets × n_sims]}
#'   of sample sizes or overdispersion for survey length compositions. Default: 100.
#' @param ln_SrvLen_theta Log-scale overdispersion array
#'   \code{[n_regions × n_sexes × n_srv_fleets]}. Default: log(1).
#' @param ln_SrvLen_theta_agg Vector \code{[n_srv_fleets]} for aggregated
#'   length composition overdispersion. Default: log(1).
#' @param SrvLen_corr_pars Array \code{[n_regions × n_sexes × n_srv_fleets × 2]}
#'   correlation parameters for length comps. Default: 0.01.
#' @param SrvLen_corr_pars_agg Vector \code{[n_srv_fleets]} for aggregated
#'   length composition correlations. Default: 0.01.
#' @param SrvLenComps_Type Array \code{[n_yrs × n_srv_fleets]} specifying
#'   length composition structure. Default: 2.
#'
#' @param comp_srvage_pop_like Integer or character vector \code{[n_srv_fleets]}
#'   specifying likelihood for population-specific survey age compositions. Default: all 0.
#' @param ISS_SrvAgeComps_pop Array \code{[n_pop × n_regions × n_yrs × n_seas × n_sexes × n_srv_fleets × n_sims]}
#'   of population-specific sample sizes or overdispersion. Default: 100.
#' @param ln_SrvAge_pop_theta Log-scale overdispersion array
#'   \code{[n_pop × n_regions × n_sexes × n_srv_fleets]}. Default: log(1).
#' @param ln_SrvAge_pop_theta_agg Array \code{[n_pop × n_srv_fleets]} for
#'   aggregated population-specific overdispersion. Default: log(1).
#' @param SrvAge_pop_corr_pars Array \code{[n_pop × n_regions × n_sexes × n_srv_fleets × 2]}
#'   correlation parameters (age AR1, sex) for population-specific age compositions. Default: 0.01.
#' @param SrvAge_pop_corr_pars_agg Array \code{[n_pop × n_srv_fleets]} for
#'   aggregated population-specific age correlations. Default: 0.01.
#' @param SrvAgeComps_pop_Type Array \code{[n_yrs × n_srv_fleets]} specifying
#'   population-specific age composition structure. Default: 2.
#' @param comp_srvlen_pop_like Integer or character vector \code{[n_srv_fleets]}
#'   specifying likelihood for population-specific survey length compositions. Default: all 0.
#' @param ISS_SrvLenComps_pop Array \code{[n_pop × n_regions × n_yrs × n_seas × n_sexes × n_srv_fleets × n_sims]}
#'   of population-specific sample sizes or overdispersion. Default: 100.
#' @param ln_SrvLen_pop_theta Array \code{[n_pop × n_regions × n_sexes × n_srv_fleets]}
#'   log-scale overdispersion for population-specific lengths. Default: log(1).
#' @param ln_SrvLen_pop_theta_agg Array \code{[n_pop × n_srv_fleets]} for
#'   aggregated population-specific length overdispersion. Default: log(1).
#' @param SrvLen_pop_corr_pars Array \code{[n_pop × n_regions × n_sexes × n_srv_fleets × 2]}
#'   correlation parameters for population-specific length comps. Default: 0.01.
#' @param SrvLen_pop_corr_pars_agg Array \code{[n_pop × n_srv_fleets]} for
#'   aggregated population-specific length correlations. Default: 0.01.
#' @param SrvLenComps_pop_Type Array \code{[n_yrs × n_srv_fleets]} specifying
#'   population-specific length composition structure. Default: 2.
#'
#' @return The input \code{sim_list} with survey-related fields appended:
#'   \code{$srv_sel}, \code{$srv_q}, \code{$ObsSrvIdx_SE}, \code{$ObsSrvIdx_pop_SE},
#'   \code{$t_srv}, \code{$srv_idx_type}, \code{$comp_srvage_like}, \code{$ISS_SrvAgeComps},
#'   \code{$ln_SrvAge_theta}, \code{$ln_SrvAge_theta_agg}, \code{$SrvAge_corr_pars_agg},
#'   \code{$SrvAge_corr_pars}, \code{$SrvAgeComps_Type}, \code{$comp_srvlen_like},
#'   \code{$ISS_SrvLenComps}, \code{$ln_SrvLen_theta}, \code{$ln_SrvLen_theta_agg},
#'   \code{$SrvLen_corr_pars_agg}, \code{$SrvLen_corr_pars}, \code{$SrvLenComps_Type},
#'   \code{$comp_srvage_pop_like}, \code{$ISS_SrvAgeComps_pop}, \code{$ln_SrvAge_pop_theta},
#'   \code{$ln_SrvAge_pop_theta_agg}, \code{$SrvAge_pop_corr_pars_agg}, \code{$SrvAge_pop_corr_pars},
#'   \code{$SrvAgeComps_pop_Type}, \code{$comp_srvlen_pop_like}, \code{$ISS_SrvLenComps_pop},
#'   \code{$ln_SrvLen_pop_theta}, \code{$ln_SrvLen_pop_theta_agg}, \code{$SrvLen_pop_corr_pars_agg},
#'   \code{$SrvLen_pop_corr_pars}, \code{$SrvLenComps_pop_Type}. Character-coded
#'   inputs are converted to integer equivalents before storage.
#'
#' @export Setup_Sim_Survey
#' @family Simulation Setup
Setup_Sim_Survey <- function(sim_list,
                             srv_sel_input,
                             ObsSrvIdx_SE = array(0.2, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas,  sim_list$n_srv_fleets)),
                             ln_sigmaSrvIdxAA = array(log(0.2), dim = c(sim_list$n_ages, sim_list$n_srv_fleets)),
                             UseSrvIdxAA = array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_srv_fleets)),
                             use_srv_idx_aa = rep(0, sim_list$n_srv_fleets),
                             ObsSrvIdx_pop_SE = array(0.2, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas,  sim_list$n_srv_fleets)),
                             srv_q_input = array(1, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_srv_fleets, sim_list$n_sims)),
                             t_srv = array(1, dim = c(sim_list$n_regions, sim_list$n_seas, sim_list$n_srv_fleets)),
                             srv_idx_type = array(1, dim = c(sim_list$n_srv_fleets)),
                             SrvIdx_LikeType = rep(0, sim_list$n_srv_fleets),
                             SrvIdx_Cov = NULL,
                             UseSrvIdx = NULL,
                             comp_srv_caal_like = rep(999, sim_list$n_srv_fleets),
                             ISS_Srv_caal = NULL,
                             ln_Srv_caal_theta = NULL,
                             ln_Srv_caal_theta_agg = NULL,
                             Srv_caal_Type = array(999, dim = c(sim_list$n_yrs, sim_list$n_srv_fleets)),
                             comp_srvage_like = rep(0, sim_list$n_srv_fleets),
                             ISS_SrvAgeComps = array(100, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims)),
                             ln_SrvAge_theta = array(log(1), dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_srv_fleets)),
                             ln_SrvAge_theta_agg = rep(log(1), sim_list$n_srv_fleets),
                             SrvAge_corr_pars_agg = rep(0.01, sim_list$n_srv_fleets),
                             SrvAge_corr_pars = array(0.01, dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_srv_fleets, 2)),
                             SrvAgeComps_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_srv_fleets)),
                             comp_srvlen_like = rep(0, sim_list$n_srv_fleets),
                             ISS_SrvLenComps = array(100, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims)),
                             ln_SrvLen_theta = array(log(1), dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_srv_fleets)),
                             ln_SrvLen_theta_agg = rep(log(1), sim_list$n_srv_fleets),
                             SrvLen_corr_pars_agg = rep(0.01, sim_list$n_srv_fleets),
                             SrvLen_corr_pars = array(0.01, dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_srv_fleets, 2)),
                             SrvLenComps_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_srv_fleets)),
                             comp_srvage_pop_like = rep(0, sim_list$n_srv_fleets),
                             ISS_SrvAgeComps_pop = array(100, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims)),
                             ln_SrvAge_pop_theta = array(log(1), dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_sexes, sim_list$n_srv_fleets)),
                             ln_SrvAge_pop_theta_agg = array(log(1), dim = c(sim_list$n_pop, sim_list$n_srv_fleets)),
                             SrvAge_pop_corr_pars = array(0.01, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_sexes, sim_list$n_srv_fleets, 2)),
                             SrvAge_pop_corr_pars_agg = array(0.01, dim = c(sim_list$n_pop, sim_list$n_srv_fleets)),
                             SrvAgeComps_pop_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_srv_fleets)),
                             comp_srvlen_pop_like = rep(0, sim_list$n_srv_fleets),
                             ISS_SrvLenComps_pop = array(100, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims)),
                             ln_SrvLen_pop_theta = array(log(1), dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_sexes, sim_list$n_srv_fleets)),
                             ln_SrvLen_pop_theta_agg = array(log(1), dim = c(sim_list$n_pop, sim_list$n_srv_fleets)),
                             SrvLen_pop_corr_pars = array(0.01, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_sexes, sim_list$n_srv_fleets, 2)),
                             SrvLen_pop_corr_pars_agg = array(0.01, dim = c(sim_list$n_pop, sim_list$n_srv_fleets)),
                             SrvLenComps_pop_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_srv_fleets))
                             ) {

  # Convert character inputs to numeric codes
  srv_idx_type <- convert_to_numeric(srv_idx_type, list(abd = 0, biom = 1, recdev = 2))
  SrvIdx_LikeType <- convert_to_numeric(SrvIdx_LikeType, list(lognormal = 0, normal = 1, mvn = 2))
  comp_srvage_like <- convert_to_numeric(comp_srvage_like, list(Multinomial = 0,  `Dirichlet-Multinomial` = 1, `iid-Logistic-Normal` = 2, `1d-Logistic-Normal` = 3, `2d-Logistic-Normal` = 4))
  comp_srvlen_like <- convert_to_numeric(comp_srvlen_like, list(Multinomial = 0, `Dirichlet-Multinomial` = 1, `iid-Logistic-Normal` = 2, `1d-Logistic-Normal` = 3, `2d-Logistic-Normal` = 4))
  SrvAgeComps_Type <- convert_to_numeric(SrvAgeComps_Type,  list(agg = 0, spltRspltS = 1, spltRjntS = 2, none = 999))
  SrvLenComps_Type <- convert_to_numeric(SrvLenComps_Type,  list(agg = 0, spltRspltS = 1, spltRjntS = 2, none = 999))
  SrvAgeComps_pop_Type <- convert_to_numeric(SrvAgeComps_pop_Type,  list(agg = 0, spltRspltS = 1, spltRjntS = 2, none = 999))
  SrvLenComps_pop_Type <- convert_to_numeric(SrvLenComps_pop_Type,  list(agg = 0, spltRspltS = 1, spltRjntS = 2, none = 999))
  comp_srvage_pop_like <- convert_to_numeric(comp_srvage_pop_like, list( Multinomial = 0, `Dirichlet-Multinomial` = 1, `iid-Logistic-Normal` = 2, `1d-Logistic-Normal` = 3, `2d-Logistic-Normal` = 4))
  comp_srvlen_pop_like <- convert_to_numeric(comp_srvlen_pop_like, list( Multinomial = 0, `Dirichlet-Multinomial` = 1, `iid-Logistic-Normal` = 2, `1d-Logistic-Normal` = 3, `2d-Logistic-Normal` = 4))

  # Validate dimensions of all input parameters
  check_sim_dimensions(srv_sel_input, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs,
                       n_ages = sim_list$n_ages, n_sexes = sim_list$n_sexes, n_pop = sim_list$n_pop, n_seas = sim_list$n_seas,
                       n_srv_fleets = sim_list$n_srv_fleets, n_sims = sim_list$n_sims, what = "srv_sel_input")
  check_sim_dimensions(srv_q_input, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs,
                       n_srv_fleets = sim_list$n_srv_fleets, n_sims = sim_list$n_sims, what = "srv_q_input")
  check_sim_dimensions(ObsSrvIdx_SE, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs, n_seas = sim_list$n_seas,
                       n_srv_fleets = sim_list$n_srv_fleets, what = "ObsSrvIdx_SE")
  check_sim_dimensions(ObsSrvIdx_pop_SE, n_pop = sim_list$n_pop, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs, n_seas = sim_list$n_seas,
                       n_srv_fleets = sim_list$n_srv_fleets, what = "ObsSrvIdx_pop_SE")
  check_sim_dimensions(t_srv, n_regions = sim_list$n_regions, n_seas = sim_list$n_seas, n_srv_fleets = sim_list$n_srv_fleets, what = "t_srv")
  check_sim_dimensions(SrvIdx_LikeType, n_srv_fleets = sim_list$n_srv_fleets, what = "SrvIdx_LikeType")

  # Multivariate normal index fleets draw from the supplied covariance rather than
  # the SE array, so the covariance is validated and factor-decomposed once here.
  srv_idx_mvn <- NULL
  if(any(SrvIdx_LikeType == 2)) {
    if(is.null(UseSrvIdx)) stop("UseSrvIdx must be supplied when any SrvIdx_LikeType is mvn, to position each observation in the covariance.")
    if(length(dim(UseSrvIdx)) != 4 || any(dim(UseSrvIdx)[c(1,3,4)] != c(sim_list$n_regions, sim_list$n_seas, sim_list$n_srv_fleets)) || dim(UseSrvIdx)[2] > sim_list$n_yrs)
      stop("UseSrvIdx must be an n_regions x (at most n_yrs) x n_seas x n_srv_fleets array.")
    srv_idx_mvn <- build_idx_factor(SrvIdx_Cov, SrvIdx_LikeType, UseSrvIdx, sim_list$n_srv_fleets, "SrvIdx_Cov")
  }

  # Validate survey age composition parameters
  check_sim_dimensions(comp_srvage_like, n_srv_fleets = sim_list$n_srv_fleets, what = "comp_srvage_like")
  check_sim_dimensions(ISS_SrvAgeComps, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs, n_seas = sim_list$n_seas,
                       n_sexes = sim_list$n_sexes, n_srv_fleets = sim_list$n_srv_fleets,
                       n_sims = sim_list$n_sims, what = "ISS_SrvAgeComps")
  check_sim_dimensions(ln_SrvAge_theta, n_regions = sim_list$n_regions, n_sexes = sim_list$n_sexes,
                       n_srv_fleets = sim_list$n_srv_fleets, what = "ln_SrvAge_theta")
  check_sim_dimensions(ln_SrvAge_theta_agg, n_srv_fleets = sim_list$n_srv_fleets, what = "ln_SrvAge_theta_agg")
  check_sim_dimensions(SrvAge_corr_pars_agg, n_srv_fleets = sim_list$n_srv_fleets, what = "SrvAge_corr_pars_agg")
  check_sim_dimensions(SrvAge_corr_pars, n_regions = sim_list$n_regions, n_sexes = sim_list$n_sexes,
                       n_srv_fleets = sim_list$n_srv_fleets, what = "SrvAge_corr_pars")
  check_sim_dimensions(SrvAgeComps_Type, n_years = sim_list$n_yrs, n_srv_fleets = sim_list$n_srv_fleets,
                       what = "SrvAgeComps_Type")
  check_sim_dimensions(comp_srvage_pop_like, n_srv_fleets = sim_list$n_srv_fleets, what = "comp_srvage_pop_like")
  check_sim_dimensions(ISS_SrvAgeComps_pop, n_pop = sim_list$n_pop, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs, n_seas = sim_list$n_seas, n_sexes = sim_list$n_sexes, n_srv_fleets = sim_list$n_srv_fleets, n_sims = sim_list$n_sims, what = "ISS_SrvAgeComps_pop")
  check_sim_dimensions(ln_SrvAge_pop_theta, n_pop = sim_list$n_pop, n_regions = sim_list$n_regions, n_sexes = sim_list$n_sexes, n_srv_fleets = sim_list$n_srv_fleets, what = "ln_SrvAge_pop_theta")
  check_sim_dimensions(ln_SrvAge_pop_theta_agg, n_pop = sim_list$n_pop, n_srv_fleets = sim_list$n_srv_fleets, what = "ln_SrvAge_pop_theta_agg")
  check_sim_dimensions(SrvAge_pop_corr_pars_agg, n_pop = sim_list$n_pop, n_srv_fleets = sim_list$n_srv_fleets, what = "SrvAge_pop_corr_pars_agg")
  check_sim_dimensions(SrvAge_pop_corr_pars, n_pop = sim_list$n_pop, n_regions = sim_list$n_regions, n_sexes = sim_list$n_sexes, n_srv_fleets = sim_list$n_srv_fleets, what = "SrvAge_pop_corr_pars")
  check_sim_dimensions(SrvAgeComps_pop_Type, n_years = sim_list$n_yrs, n_srv_fleets = sim_list$n_srv_fleets, what = "SrvAgeComps_pop_Type")


  # Validate survey length composition parameters
  check_sim_dimensions(comp_srvlen_like, n_srv_fleets = sim_list$n_srv_fleets, what = "comp_srvlen_like")
  check_sim_dimensions(ISS_SrvLenComps, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs, n_seas = sim_list$n_seas,
                       n_sexes = sim_list$n_sexes, n_srv_fleets = sim_list$n_srv_fleets,
                       n_sims = sim_list$n_sims, what = "ISS_SrvLenComps")
  check_sim_dimensions(ln_SrvLen_theta, n_regions = sim_list$n_regions, n_sexes = sim_list$n_sexes,
                       n_srv_fleets = sim_list$n_srv_fleets, what = "ln_SrvLen_theta")
  check_sim_dimensions(ln_SrvLen_theta_agg, n_srv_fleets = sim_list$n_srv_fleets, what = "ln_SrvLen_theta_agg")
  check_sim_dimensions(SrvLen_corr_pars_agg, n_srv_fleets = sim_list$n_srv_fleets, what = "SrvLen_corr_pars_agg")
  check_sim_dimensions(SrvLen_corr_pars, n_regions = sim_list$n_regions, n_sexes = sim_list$n_sexes,
                       n_srv_fleets = sim_list$n_srv_fleets, what = "SrvLen_corr_pars")
  check_sim_dimensions(SrvLenComps_Type, n_years = sim_list$n_yrs, n_srv_fleets = sim_list$n_srv_fleets,
                       what = "SrvLenComps_Type")
  check_sim_dimensions(comp_srvlen_pop_like, n_srv_fleets = sim_list$n_srv_fleets, what = "comp_srvlen_pop_like")
  check_sim_dimensions(ISS_SrvLenComps_pop, n_pop = sim_list$n_pop, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs, n_seas = sim_list$n_seas, n_sexes = sim_list$n_sexes, n_srv_fleets = sim_list$n_srv_fleets, n_sims = sim_list$n_sims, what = "ISS_SrvLenComps_pop")
  check_sim_dimensions(ln_SrvLen_pop_theta, n_pop = sim_list$n_pop, n_regions = sim_list$n_regions, n_sexes = sim_list$n_sexes, n_srv_fleets = sim_list$n_srv_fleets, what = "ln_SrvLen_pop_theta")
  check_sim_dimensions(ln_SrvLen_pop_theta_agg, n_pop = sim_list$n_pop, n_srv_fleets = sim_list$n_srv_fleets, what = "ln_SrvLen_pop_theta_agg")
  check_sim_dimensions(SrvLen_pop_corr_pars_agg, n_pop = sim_list$n_pop, n_srv_fleets = sim_list$n_srv_fleets, what = "SrvLen_pop_corr_pars_agg")
  check_sim_dimensions(SrvLen_pop_corr_pars, n_pop = sim_list$n_pop, n_regions = sim_list$n_regions, n_sexes = sim_list$n_sexes, n_srv_fleets = sim_list$n_srv_fleets, what = "SrvLen_pop_corr_pars")
  check_sim_dimensions(SrvLenComps_pop_Type, n_years = sim_list$n_yrs, n_srv_fleets = sim_list$n_srv_fleets, what = "SrvLenComps_pop_Type")

  # output into list
  sim_list$srv_sel <- srv_sel_input
  sim_list$srv_q <- srv_q_input
  sim_list$ObsSrvIdx_SE <- ObsSrvIdx_SE
  sim_list$ln_sigmaSrvIdxAA <- ln_sigmaSrvIdxAA
  sim_list$UseSrvIdxAA <- UseSrvIdxAA
  sim_list$use_srv_idx_aa <- use_srv_idx_aa
  sim_list$ObsSrvIdx_pop_SE <- ObsSrvIdx_pop_SE
  sim_list$t_srv <- t_srv
  sim_list$srv_idx_type <- srv_idx_type
  sim_list$SrvIdx_LikeType <- SrvIdx_LikeType # survey index error structure
  if(!is.null(srv_idx_mvn)) {
    sim_list$srv_idx_mvn <- srv_idx_mvn # factor parameters for mvn index fleets
    sim_list$srv_idx_u <- matrix(NA_real_, sim_list$n_srv_fleets, sim_list$n_sims) # shared factor draw, filled per fleet and replicate
  }

  # Survey age compositions
  sim_list$comp_srvage_like <- comp_srvage_like
  sim_list$ISS_SrvAgeComps <- ISS_SrvAgeComps
  sim_list$ln_SrvAge_theta <- ln_SrvAge_theta
  sim_list$ln_SrvAge_theta_agg <- ln_SrvAge_theta_agg
  sim_list$SrvAge_corr_pars_agg <- SrvAge_corr_pars_agg
  sim_list$SrvAge_corr_pars <- SrvAge_corr_pars
  sim_list$SrvAgeComps_Type <- SrvAgeComps_Type

  # Survey conditional age-at-length
  comp_srv_caal_like <- convert_to_numeric(comp_srv_caal_like, list(Multinomial = 0, `Dirichlet-Multinomial` = 1, none = 999))
  Srv_caal_Type <- convert_to_numeric(Srv_caal_Type, list(agg = 0, spltRspltS = 1, spltRjntS = 2, none = 999))
  if(is.null(ln_Srv_caal_theta)) ln_Srv_caal_theta <- array(log(1), dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_srv_fleets))
  if(is.null(ln_Srv_caal_theta_agg)) ln_Srv_caal_theta_agg <- rep(log(1), sim_list$n_srv_fleets)
  sim_list$do_srv_caal <- !is.null(ISS_Srv_caal) && any(comp_srv_caal_like != 999)
  if(sim_list$do_srv_caal) {
    if(is.null(sim_list$n_lens)) stop("ISS_Srv_caal was supplied, but the simulation has no length bins (n_lens is NULL)")
    if(length(dim(ISS_Srv_caal)) != 7 || !all(dim(ISS_Srv_caal) == c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_lens, sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims)))
      stop("Dimensions of ISS_Srv_caal are not correct. Should be n_regions, n_years, n_seas, n_lens, n_sexes, n_srv_fleets, and n_sims")
  }
  sim_list$comp_srv_caal_like <- comp_srv_caal_like
  sim_list$ISS_Srv_caal <- ISS_Srv_caal
  sim_list$ln_Srv_caal_theta <- ln_Srv_caal_theta
  sim_list$ln_Srv_caal_theta_agg <- ln_Srv_caal_theta_agg
  sim_list$Srv_caal_Type <- Srv_caal_Type

  # Survey length compositions
  sim_list$comp_srvlen_like <- comp_srvlen_like
  sim_list$ISS_SrvLenComps <- ISS_SrvLenComps
  sim_list$ln_SrvLen_theta <- ln_SrvLen_theta
  sim_list$ln_SrvLen_theta_agg <- ln_SrvLen_theta_agg
  sim_list$SrvLen_corr_pars_agg <- SrvLen_corr_pars_agg
  sim_list$SrvLen_corr_pars <- SrvLen_corr_pars
  sim_list$SrvLenComps_Type <- SrvLenComps_Type

  # Population-specific stuff
  sim_list$comp_srvage_pop_like <- comp_srvage_pop_like
  sim_list$ISS_SrvAgeComps_pop <- ISS_SrvAgeComps_pop
  sim_list$ln_SrvAge_pop_theta <- ln_SrvAge_pop_theta
  sim_list$ln_SrvAge_pop_theta_agg <- ln_SrvAge_pop_theta_agg
  sim_list$SrvAge_pop_corr_pars <- SrvAge_pop_corr_pars
  sim_list$SrvAge_pop_corr_pars_agg <- SrvAge_pop_corr_pars_agg
  sim_list$SrvAgeComps_pop_Type <- SrvAgeComps_pop_Type

  sim_list$comp_srvlen_pop_like <- comp_srvlen_pop_like
  sim_list$ISS_SrvLenComps_pop <- ISS_SrvLenComps_pop
  sim_list$ln_SrvLen_pop_theta <- ln_SrvLen_pop_theta
  sim_list$ln_SrvLen_pop_theta_agg <- ln_SrvLen_pop_theta_agg
  sim_list$SrvLen_pop_corr_pars <- SrvLen_pop_corr_pars
  sim_list$SrvLen_pop_corr_pars_agg <- SrvLen_pop_corr_pars_agg
  sim_list$SrvLenComps_pop_Type <- SrvLenComps_pop_Type

  return(sim_list)

} # end function
