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
#' @param comp_srvage_like Integer or character vector \code{[n_srv_fleets]}
#'   specifying likelihood for survey age compositions. Default: all 0
#'   (multinomial). Options: 0/“Multinomial”, 1/“Dirichlet-Multinomial”,
#'   2/“iid-Logistic-Normal”, 3/“1d-Logistic-Normal”, 4/“2d-Logistic-Normal”.
#' @param ISS_SrvAgeComps Array \code{[n_regions × n_yrs × n_seas × n_sexes × n_srv_fleets × n_sims]}
#'   of sample sizes or overdispersion for survey age compositions. Default: 100.
#' @param ln_SrvAge_theta Log-scale overdispersion array
#'   \code{[n_regions × n_sexes × n_srv_fleets]}. Used for likelihoods 1–4.
#'   Default: log(1).
#' @param ln_SrvAge_theta_agg Log-scale overdispersion for aggregated survey
#'   age compositions, vector \code{[n_srv_fleets]}. Default: log(1).
#' @param SrvAge_corr_pars Correlation parameters array
#'   \code{[n_regions × n_sexes × n_srv_fleets × 2]} (age AR1, sex). Only for
#'   likelihoods 3–4. Default: 0.01.
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
                             ObsSrvIdx_pop_SE = array(0.2, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas,  sim_list$n_srv_fleets)),
                             srv_q_input = array(1, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_srv_fleets, sim_list$n_sims)),
                             t_srv = array(1, dim = c(sim_list$n_regions, sim_list$n_seas, sim_list$n_srv_fleets)),
                             srv_idx_type = array(1, dim = c(sim_list$n_srv_fleets)),
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
  srv_idx_type <- convert_to_numeric(srv_idx_type, list(abd = 0, biom = 1))
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
  sim_list$ObsSrvIdx_pop_SE <- ObsSrvIdx_pop_SE
  sim_list$t_srv <- t_srv
  sim_list$srv_idx_type <- srv_idx_type

  # Survey age compositions
  sim_list$comp_srvage_like <- comp_srvage_like
  sim_list$ISS_SrvAgeComps <- ISS_SrvAgeComps
  sim_list$ln_SrvAge_theta <- ln_SrvAge_theta
  sim_list$ln_SrvAge_theta_agg <- ln_SrvAge_theta_agg
  sim_list$SrvAge_corr_pars_agg <- SrvAge_corr_pars_agg
  sim_list$SrvAge_corr_pars <- SrvAge_corr_pars
  sim_list$SrvAgeComps_Type <- SrvAgeComps_Type

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

#' Set up observed survey indices and composition data
#'
#' Ingests observed survey index, age composition, and length composition data
#' (both pooled and population-specific) into \code{input_list$data},
#' initialises overdispersion and correlation starting values in
#' \code{input_list$par}, and constructs parameter maps via
#' \code{\link{do_comp_theta_mapping}} and \code{\link{do_comp_corr_pars_mapping}}
#' (called with \code{comp_prefix = "SrvAge"}/\code{"SrvLen"} and
#' \code{fleet_field = "n_srv_fleets"}). When \code{ISS_SrvAgeComps},
#' \code{ISS_SrvLenComps}, \code{ISS_SrvAgeComps_pop}, or
#' \code{ISS_SrvLenComps_pop} is \code{NULL}, input sample sizes are derived
#' automatically by summing observed composition counts across the appropriate
#' dimensions each year. Must be called after \code{\link{Setup_Mod_Dim}} and
#' before model compilation.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, \code{$map},
#'   and \code{$verbose} sublists, as returned by upstream setup functions.
#' @param ObsSrvIdx Observed survey index array
#'   \code{[n_regions × n_years × n_seas × n_srv_fleets]}.
#' @param ObsSrvIdx_SE Lognormal standard errors for \code{ObsSrvIdx}, same
#'   dimensions \code{[n_regions × n_years × n_seas × n_srv_fleets]}.
#' @param UseSrvIdx Binary indicator array
#'   \code{[n_regions × n_years × n_seas × n_srv_fleets]}. \code{1} = include
#'   in likelihood; \code{0} = exclude.
#' @param srv_idx_type Character vector \code{[n_srv_fleets]} specifying the
#'   index type per fleet. One of \code{"biom"} (biomass), \code{"abd"}
#'   (abundance), or \code{"none"} (no index for that fleet). Converted to
#'   integer codes (\code{1}, \code{0}, \code{999}) before storage.
#' @param ObsSrvIdx_pop Observed population-specific survey index array
#'   \code{[n_pop × n_regions × n_years × n_seas × n_srv_fleets]}.
#' @param ObsSrvIdx_pop_SE Lognormal standard errors for \code{ObsSrvIdx_pop},
#'   same dimensions \code{[n_pop × n_regions × n_years × n_seas × n_srv_fleets]}.
#' @param UseSrvIdx_pop Binary indicator array
#'   \code{[n_pop × n_regions × n_years × n_seas × n_srv_fleets]}. \code{1} =
#'   include population-specific index in likelihood; \code{0} = exclude.
#'   Default: all zeros.
#' @param ObsSrvAgeComps Observed survey age compositions, array
#'   \code{[n_regions × n_years × n_seas × n_ages × n_sexes × n_srv_fleets]}.
#'   Values may be counts or proportions on a comparable scale.
#' @param UseSrvAgeComps Binary indicator array
#'   \code{[n_regions × n_years × n_seas × n_srv_fleets]}. \code{1} = fit age
#'   compositions; \code{0} = exclude.
#' @param ISS_SrvAgeComps Input sample sizes for survey age compositions, array
#'   \code{[n_regions × n_years × n_seas × n_sexes × n_srv_fleets]}, or
#'   \code{NULL} to derive automatically by summing \code{ObsSrvAgeComps}
#'   across the age dimension each year, respecting \code{SrvAgeComps_Type}.
#' @param ObsSrvLenComps Observed survey length compositions, array
#'   \code{[n_regions × n_years × n_seas × n_lens × n_sexes × n_srv_fleets]}.
#'   Only validated when \code{input_list$data$fit_lengths = 1} in \code{$data}.
#' @param UseSrvLenComps Binary indicator array
#'   \code{[n_regions × n_years × n_seas × n_srv_fleets]}. \code{1} = fit
#'   length compositions; \code{0} = exclude.
#' @param ISS_SrvLenComps Input sample sizes for survey length compositions,
#'   same structure as \code{ISS_SrvAgeComps}, or \code{NULL} for automatic
#'   derivation from \code{ObsSrvLenComps}.
#' @param SrvAgeComps_LikeType Character vector \code{[n_srv_fleets]}
#'   specifying the likelihood for survey age compositions. One of
#'   \code{"none"}, \code{"Multinomial"}, \code{"Dirichlet-Multinomial"},
#'   \code{"iid-Logistic-Normal"}, \code{"1d-Logistic-Normal"},
#'   \code{"2d-Logistic-Normal"}. Converted to integer codes
#'   (\code{999}, \code{0}–\code{4}) before storage.
#' @param SrvLenComps_LikeType Character vector \code{[n_srv_fleets]}
#'   specifying the likelihood for survey length compositions. Same options
#'   as \code{SrvAgeComps_LikeType}.
#' @param SrvAgeComps_Type Character vector defining the survey age composition
#'   structure per fleet and year range. Each element follows the format
#'   \code{"<type>_Year_<start>-<end>_Fleet_<fleet>"}. Use \code{"terminal"}
#'   in place of the end year to extend to the final model year. Valid types:
#'   \describe{
#'     \item{\code{"agg"}}{Aggregated across regions and sexes. Not compatible
#'       with \code{"2d-Logistic-Normal"}.}
#'     \item{\code{"spltRspltS"}}{Split by region and sex.}
#'     \item{\code{"spltRjntS"}}{Split by region, joint across sexes.}
#'     \item{\code{"none"}}{No composition data used.}
#'   }
#'   Parsed into a \code{[n_years × n_srv_fleets]} integer matrix before
#'   storage. An error is raised if any cell remains \code{NA} after parsing,
#'   indicating an incomplete year range specification.
#' @param SrvLenComps_Type Character vector defining the survey length
#'   composition structure. Same format and options as \code{SrvAgeComps_Type}.
#' @param ObsSrvAgeComps_pop Observed population-specific survey age
#'   composition array
#'   \code{[n_pop × n_regions × n_years × n_seas × n_ages × n_sexes × n_srv_fleets]}.
#'   Required when any element of \code{UseSrvAgeComps_pop} is \code{1}.
#' @param UseSrvAgeComps_pop Binary indicator array
#'   \code{[n_pop × n_regions × n_years × n_seas × n_srv_fleets]}.
#'   \code{1} = fit population-specific age compositions; \code{0} = exclude.
#'   Default: all zeros.
#' @param ISS_SrvAgeComps_pop Input sample size array for population-specific
#'   survey age compositions
#'   \code{[n_pop × n_regions × n_years × n_seas × n_sexes × n_srv_fleets]}.
#'   If \code{NULL} (default), computed automatically by summing
#'   \code{ObsSrvAgeComps_pop} within each population–year–fleet–season–region
#'   cell according to \code{SrvAgeComps_pop_Type}.
#' @param ObsSrvLenComps_pop Observed population-specific survey length
#'   composition array
#'   \code{[n_pop × n_regions × n_years × n_seas × n_lens × n_sexes × n_srv_fleets]}.
#'   Required when \code{input_list$data$fit_lengths == 1} and any element of
#'   \code{UseSrvLenComps_pop} is \code{1}.
#' @param UseSrvLenComps_pop Binary indicator array
#'   \code{[n_pop × n_regions × n_years × n_seas × n_srv_fleets]}.
#'   \code{1} = fit population-specific length compositions; \code{0} = exclude.
#'   Default: all zeros.
#' @param ISS_SrvLenComps_pop Input sample size array for population-specific
#'   survey length compositions
#'   \code{[n_pop × n_regions × n_years × n_seas × n_sexes × n_srv_fleets]}.
#'   If \code{NULL} (default), derived automatically from
#'   \code{ObsSrvLenComps_pop}.
#' @param SrvAgeComps_pop_LikeType Character vector of length
#'   \code{n_srv_fleets} specifying the likelihood for population-specific
#'   survey age compositions. Same options as \code{SrvAgeComps_LikeType}.
#'   Default: \code{"none"} for all fleets.
#' @param SrvLenComps_pop_LikeType Character vector of length
#'   \code{n_srv_fleets} specifying the likelihood for population-specific
#'   survey length compositions. Same options as \code{SrvLenComps_LikeType}.
#'   Default: \code{"none"} for all fleets.
#' @param SrvAgeComps_pop_Type Character vector defining the composition
#'   structure for population-specific survey age compositions. Same format and
#'   options as \code{SrvAgeComps_Type}. Default: \code{"none"} for all fleets
#'   across all years.
#' @param SrvLenComps_pop_Type Character vector defining the composition
#'   structure for population-specific survey length compositions. Same format
#'   and options as \code{SrvLenComps_Type}. Default: \code{"none"} for all
#'   fleets across all years.
#' @param ... Optional named starting values for overdispersion and correlation
#'   parameters.
#'
#' @return The input \code{input_list} with survey data stored in \code{$data}
#'   (\code{ObsSrvIdx}, \code{ObsSrvIdx_SE}, \code{UseSrvIdx},
#'   \code{ObsSrvIdx_pop}, \code{ObsSrvIdx_pop_SE}, \code{UseSrvIdx_pop},
#'   \code{ObsSrvAgeComps}, \code{UseSrvAgeComps}, \code{ISS_SrvAgeComps},
#'   \code{ObsSrvLenComps}, \code{UseSrvLenComps}, \code{ISS_SrvLenComps},
#'   \code{ObsSrvAgeComps_pop}, \code{UseSrvAgeComps_pop},
#'   \code{ISS_SrvAgeComps_pop}, \code{ObsSrvLenComps_pop},
#'   \code{UseSrvLenComps_pop}, \code{ISS_SrvLenComps_pop},
#'   \code{SrvAgeComps_LikeType}, \code{SrvLenComps_LikeType},
#'   \code{SrvAgeComps_pop_LikeType}, \code{SrvLenComps_pop_LikeType},
#'   \code{SrvAgeComps_Type}, \code{SrvLenComps_Type},
#'   \code{SrvAgeComps_pop_Type}, \code{SrvLenComps_pop_Type},
#'   \code{srv_idx_type}); overdispersion and correlation starting values in
#'   \code{$par}; and factor maps in \code{$map} for all pooled and
#'   population-specific overdispersion and correlation parameter arrays.
#'
#' @export Setup_Mod_SrvIdx_and_Comps
#' @importFrom stringr str_detect
#' @family Model Setup
Setup_Mod_SrvIdx_and_Comps <- function(input_list,
                                       ObsSrvIdx,
                                       ObsSrvIdx_SE,
                                       UseSrvIdx,
                                       ObsSrvIdx_pop = NULL,
                                       ObsSrvIdx_pop_SE = NULL,
                                       UseSrvIdx_pop = array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_srv_fleets)),
                                       srv_idx_type,
                                       ObsSrvAgeComps,
                                       UseSrvAgeComps,
                                       ObsSrvLenComps,
                                       UseSrvLenComps,
                                       ISS_SrvAgeComps = NULL,
                                       ISS_SrvLenComps = NULL,
                                       SrvAgeComps_LikeType,
                                       SrvLenComps_LikeType,
                                       SrvAgeComps_Type,
                                       SrvLenComps_Type,
                                       ObsSrvAgeComps_pop = NULL,
                                       UseSrvAgeComps_pop = array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_srv_fleets)),
                                       ISS_SrvAgeComps_pop = NULL,
                                       ObsSrvLenComps_pop = NULL,
                                       UseSrvLenComps_pop = array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_srv_fleets)),
                                       ISS_SrvLenComps_pop = NULL,
                                       SrvAgeComps_pop_LikeType = rep("none", input_list$data$n_srv_fleets),
                                       SrvLenComps_pop_LikeType = rep("none", input_list$data$n_srv_fleets),
                                       SrvAgeComps_pop_Type = paste("none_Year_1-terminal_Fleet_", 1:input_list$data$n_srv_fleets, sep = ''),
                                       SrvLenComps_pop_Type = paste("none_Year_1-terminal_Fleet_", 1:input_list$data$n_srv_fleets, sep = ''),
                                       ...
                                       ) {

  messages_list <<- character(0) # string to attach to for printing messages
  starting_values <- list(...)
  if(input_list$store_config) input_list$config$Setup_Mod_SrvIdx_and_Comps <- mget(names(formals()))[-1]

  # Input Validation --------------------------------------------------------

  # Survey Indices
  check_data_dimensions(ObsSrvIdx, n_regions = input_list$data$n_regions, n_seas = input_list$data$n_seas, n_years = length(input_list$data$years), n_srv_fleets = input_list$data$n_srv_fleets, what = 'ObsSrvIdx')
  check_data_dimensions(ObsSrvIdx_SE, n_regions = input_list$data$n_regions, n_seas = input_list$data$n_seas, n_years = length(input_list$data$years), n_srv_fleets = input_list$data$n_srv_fleets, what = 'ObsSrvIdx_SE')
  check_data_dimensions(UseSrvIdx, n_regions = input_list$data$n_regions, n_seas = input_list$data$n_seas, n_years = length(input_list$data$years), n_srv_fleets = input_list$data$n_srv_fleets, what = 'UseSrvIdx')
  if(any(UseSrvIdx_pop == 1)) {
    check_data_dimensions(ObsSrvIdx_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_srv_fleets = input_list$data$n_srv_fleets, what = 'ObsSrvIdx_pop')
    check_data_dimensions(ObsSrvIdx_pop_SE, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_srv_fleets = input_list$data$n_srv_fleets, what = 'ObsSrvIdx_pop_SE')
    check_data_dimensions(UseSrvIdx_pop, n_pop = input_list$data$n_pop,  n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_srv_fleets = input_list$data$n_srv_fleets, what = 'UseSrvIdx_pop')
  }
  if(!all(srv_idx_type %in% c("biom", "abd", "none"))) stop("Invalid specification for srv_idx_type. Should be either abd, biom, or none")

  # Survey compositions
  check_data_dimensions(ObsSrvAgeComps, n_regions = input_list$data$n_regions, n_seas = input_list$data$n_seas, n_years = length(input_list$data$years), n_sexes = input_list$data$n_sexes, n_srv_fleets = input_list$data$n_srv_fleets, what = 'ObsSrvAgeComps')
  check_data_dimensions(UseSrvAgeComps, n_regions = input_list$data$n_regions, n_seas = input_list$data$n_seas, n_years = length(input_list$data$years), n_srv_fleets = input_list$data$n_srv_fleets, what = 'UseSrvAgeComps')
  check_data_dimensions(UseSrvLenComps, n_regions = input_list$data$n_regions, n_seas = input_list$data$n_seas, n_years = length(input_list$data$years), n_srv_fleets = input_list$data$n_srv_fleets, what = 'UseSrvLenComps')
  if(input_list$data$fit_lengths == 1) check_data_dimensions(ObsSrvLenComps, n_seas = input_list$data$n_seas, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_lens = length(input_list$data$lens), n_sexes = input_list$data$n_sexes, n_srv_fleets = input_list$data$n_srv_fleets, what = 'ObsSrvLenComps')
  if(!is.null(ISS_SrvAgeComps)) check_data_dimensions(ISS_SrvAgeComps, n_seas = input_list$data$n_seas, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_sexes = input_list$data$n_sexes, n_srv_fleets = input_list$data$n_srv_fleets, what = 'ISS_SrvAgeComps')
  if(!is.null(ISS_SrvLenComps)) check_data_dimensions(ISS_SrvLenComps, n_seas = input_list$data$n_seas, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_sexes = input_list$data$n_sexes, n_srv_fleets = input_list$data$n_srv_fleets, what = 'ISS_SrvLenComps')
  check_data_dimensions(SrvAgeComps_LikeType, n_srv_fleets = input_list$data$n_srv_fleets, what = 'SrvAgeComps_LikeType')
  check_data_dimensions(SrvLenComps_LikeType, n_srv_fleets = input_list$data$n_srv_fleets, what = 'SrvLenComps_LikeType')
  if(!all(SrvAgeComps_LikeType %in% c("none", "Multinomial", "Dirichlet-Multinomial", "iid-Logistic-Normal", "1d-Logistic-Normal", "2d-Logistic-Normal")))
    stop("Invalid specification for SrvAgeComps_LikeType Should be either none, Multinomial, Dirichlet-Multinomial, iid-Logistic-Normal, 1d-Logistic-Normal, 2d-Logistic-Normal")
  if(!all(SrvLenComps_LikeType %in% c("none", "Multinomial", "Dirichlet-Multinomial", "iid-Logistic-Normal", "1d-Logistic-Normal", "2d-Logistic-Normal")))
    stop("Invalid specification for SrvLenComps_LikeType Should be either none, Multinomial, Dirichlet-Multinomial, iid-Logistic-Normal, 1d-Logistic-Normal, 2d-Logistic-Normal")

  # Survey compositions (population-specific)
  if(any(UseSrvAgeComps_pop == 1)) check_data_dimensions(ObsSrvAgeComps_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_sexes = input_list$data$n_sexes, n_srv_fleets = input_list$data$n_srv_fleets, what = 'ObsSrvAgeComps_pop')
  check_data_dimensions(UseSrvAgeComps_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_srv_fleets = input_list$data$n_srv_fleets, what = 'UseSrvAgeComps_pop')
  check_data_dimensions(UseSrvLenComps_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_srv_fleets = input_list$data$n_srv_fleets, what = 'UseSrvLenComps_pop')
  if(input_list$data$fit_lengths == 1 && any(UseSrvLenComps_pop == 1)) check_data_dimensions(ObsSrvLenComps_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_seas = input_list$data$n_seas, n_years = length(input_list$data$years), n_lens = length(input_list$data$lens), n_sexes = input_list$data$n_sexes, n_srv_fleets = input_list$data$n_srv_fleets, what = 'ObsSrvLenComps_pop')
  if(!is.null(ISS_SrvAgeComps_pop)) check_data_dimensions(ISS_SrvAgeComps_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_seas = input_list$data$n_seas, n_years = length(input_list$data$years), n_sexes = input_list$data$n_sexes, n_srv_fleets = input_list$data$n_srv_fleets, what = 'ISS_SrvAgeComps_pop')
  if(!is.null(ISS_SrvLenComps_pop)) check_data_dimensions(ISS_SrvLenComps_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_seas = input_list$data$n_seas, n_years = length(input_list$data$years), n_sexes = input_list$data$n_sexes, n_srv_fleets = input_list$data$n_srv_fleets, what = 'ISS_SrvLenComps_pop')
  check_data_dimensions(SrvAgeComps_pop_LikeType, n_srv_fleets = input_list$data$n_srv_fleets, what = 'SrvAgeComps_pop_LikeType')
  check_data_dimensions(SrvLenComps_pop_LikeType, n_srv_fleets = input_list$data$n_srv_fleets, what = 'SrvLenComps_pop_LikeType')
  if(!all(SrvAgeComps_pop_LikeType %in% c("none", "Multinomial", "Dirichlet-Multinomial", "iid-Logistic-Normal", "1d-Logistic-Normal", "2d-Logistic-Normal")))
    stop("Invalid specification for SrvAgeComps_pop_LikeType Should be either none, Multinomial, Dirichlet-Multinomial, iid-Logistic-Normal, 1d-Logistic-Normal, 2d-Logistic-Normal")
  if(!all(SrvLenComps_pop_LikeType %in% c("none", "Multinomial", "Dirichlet-Multinomial", "iid-Logistic-Normal", "1d-Logistic-Normal", "2d-Logistic-Normal")))
    stop("Invalid specification for SrvLenComps_pop_LikeType Should be either none, Multinomial, Dirichlet-Multinomial, iid-Logistic-Normal, 1d-Logistic-Normal, 2d-Logistic-Normal")

  # checking to make sure defaults are not applied
  if(any(UseSrvAgeComps_pop == 1)) {
    if(is.null(ObsSrvAgeComps_pop)) stop("ObsSrvAgeComps_pop is NULL, but UseSrvAgeComps_pop contains 1s!")
    if(any(str_detect(SrvAgeComps_pop_LikeType, "none"))) warning("SrvAgeComps_pop_LikeType has nones, but UseSrvAgeComps_pop contains 1s! Please verify!")
    if(any(str_detect(SrvAgeComps_pop_Type, "none"))) warning("SrvAgeComps_pop_Type has nones, but UseSrvAgeComps_pop contains 1s! Please verify!")
  }

  if(any(UseSrvLenComps_pop == 1)) {
    if(is.null(ObsSrvLenComps_pop)) stop("ObsSrvLenComps_pop is NULL, but UseSrvAgeComps_pop contains 1s!")
    if(any(str_detect(SrvLenComps_pop_LikeType, "none"))) warning("SrvLenComps_pop_LikeType has nones, but UseSrvAgeComps_pop contains 1s! Please verify!")
    if(any(str_detect(SrvLenComps_pop_Type, "none"))) warning("SrvLenComps_pop_Type has nones, but UseSrvAgeComps_pop contains 1s! Please verify!")
  }

  # Survey Index Options ----------------------------------------------------

  srv_idx_type_vals <- array(NA, dim = c( input_list$data$n_srv_fleets))
  for(f in 1:input_list$data$n_srv_fleets) {
    if(srv_idx_type[f] == 'biom') srv_idx_type_vals[f] <- 1 # biomass
    if(srv_idx_type[f] == 'abd') srv_idx_type_vals[f] <- 0 # abundance
    if(srv_idx_type[f] == 'none') srv_idx_type_vals[f] <- 999 # none
    collect_message(paste("Survey Index", "for survey fleet", f, "specified as:" , srv_idx_type[f]))
  } # end f loop


  # Survey Age Composition Options ------------------------------------------

  comp_srvage_like_vals <- vector()
  for(f in 1:input_list$data$n_srv_fleets) {
    if(SrvAgeComps_LikeType[f] == 'none') comp_srvage_like_vals <- c(comp_srvage_like_vals, 999)
    if(SrvAgeComps_LikeType[f] == "Multinomial") comp_srvage_like_vals <- c(comp_srvage_like_vals, 0)
    if(SrvAgeComps_LikeType[f] == "Dirichlet-Multinomial") comp_srvage_like_vals <- c(comp_srvage_like_vals, 1)
    if(SrvAgeComps_LikeType[f] == "iid-Logistic-Normal") comp_srvage_like_vals <- c(comp_srvage_like_vals, 2)
    if(SrvAgeComps_LikeType[f] == "1d-Logistic-Normal") comp_srvage_like_vals <- c(comp_srvage_like_vals, 3)
    if(SrvAgeComps_LikeType[f] == "2d-Logistic-Normal") comp_srvage_like_vals <- c(comp_srvage_like_vals, 4)
    collect_message(paste("Survey Age Composition Likelihoods", "for survey fleet", f, "specified as:" , SrvAgeComps_LikeType[f]))
  } # end f loop

  # Specifying composition type
  SrvAgeComps_Type_Mat <- array(NA, dim = c(length(input_list$data$years), input_list$data$n_srv_fleets))
  for(i in 1:length(SrvAgeComps_Type)) {

    # Extract out components from list
    tmp <- SrvAgeComps_Type[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    comps_type_tmp <- tmp_vec[1] # get composition type
    fleet <- as.numeric(tmp_vec[5]) # extract fleet index

    # Checking character string
    if(!comps_type_tmp %in% c("agg", "spltRspltS", "spltRjntS", 'none')) stop("SrvAgeComps_Type not specified correctly. Must be one of: agg, spltRspltS, spltRjntS, none")
    if(!fleet %in% c(1:input_list$data$n_srv_fleets)) stop("Invalid fleet specified for SrvAgeComps_Type This needs to be specified as CompType_Year_x-y_Fleet_x")

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
      if(comp_srvage_like_vals[fleet] == 4) stop("Age composition likelihood specified as 2d-Logistic-Normal, but composition type is aggregated. This is not valid.")
      comps_type_val <- 0
    }
    if(comps_type_tmp == "spltRspltS") comps_type_val <- 1
    if(comps_type_tmp == "spltRjntS") comps_type_val <- 2
    if(comps_type_tmp == "none") comps_type_val <- 999

    # input into matrix
    SrvAgeComps_Type_Mat[years,fleet] <- comps_type_val
  } # end i

  if(any(is.na(SrvAgeComps_Type_Mat))) stop("SrvAgeComps_Type_Mat is returning an NA. Did you update the year range of SrvAgeComps_Type_Mat?")

  # Specifying composition likelihood for population-specific data
  comp_srvage_pop_like_vals <- vector()
  for(f in 1:input_list$data$n_srv_fleets) {
    if(SrvAgeComps_pop_LikeType[f] == 'none') comp_srvage_pop_like_vals <- c(comp_srvage_pop_like_vals, 999)
    if(SrvAgeComps_pop_LikeType[f] == "Multinomial") comp_srvage_pop_like_vals <- c(comp_srvage_pop_like_vals, 0)
    if(SrvAgeComps_pop_LikeType[f] == "Dirichlet-Multinomial") comp_srvage_pop_like_vals <- c(comp_srvage_pop_like_vals, 1)
    if(SrvAgeComps_pop_LikeType[f] == "iid-Logistic-Normal") comp_srvage_pop_like_vals <- c(comp_srvage_pop_like_vals, 2)
    if(SrvAgeComps_pop_LikeType[f] == "1d-Logistic-Normal") comp_srvage_pop_like_vals <- c(comp_srvage_pop_like_vals, 3)
    if(SrvAgeComps_pop_LikeType[f] == "2d-Logistic-Normal") comp_srvage_pop_like_vals <- c(comp_srvage_pop_like_vals, 4)
    collect_message(paste("Population Survey Age Composition Likelihoods", "for survey fleet", f, "specified as:" , SrvAgeComps_pop_LikeType[f]))
  } # end f loop

  # Specifying composition type
  SrvAgeComps_pop_Type_Mat <- array(NA, dim = c(length(input_list$data$years), input_list$data$n_srv_fleets))
  for(i in 1:length(SrvAgeComps_pop_Type)) {

    # Extract out components from list
    tmp <- SrvAgeComps_pop_Type[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    comps_type_tmp <- tmp_vec[1] # get composition type
    fleet <- as.numeric(tmp_vec[5]) # extract fleet index

    # Checking character string
    if(!comps_type_tmp %in% c("agg", "spltRspltS", "spltRjntS", 'none')) stop("SrvAgeComps_pop_Type not specified correctly. Must be one of: agg, spltRspltS, spltRjntS, none")
    if(!fleet %in% c(1:input_list$data$n_srv_fleets)) stop("Invalid fleet specified for SrvAgeComps_pop_Type. This needs to be specified as CompType_Year_x-y_Fleet_x")

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
      if(comp_srvage_pop_like_vals[fleet] == 4) stop("Population Age composition likelihood specified as 2d-Logistic-Normal, but composition type is aggregated. This is not valid.")
      comps_type_val <- 0
    }
    if(comps_type_tmp == "spltRspltS") comps_type_val <- 1
    if(comps_type_tmp == "spltRjntS") comps_type_val <- 2
    if(comps_type_tmp == "none") comps_type_val <- 999

    # input into matrix
    SrvAgeComps_pop_Type_Mat[years,fleet] <- comps_type_val
  } # end i

  if(any(is.na(SrvAgeComps_pop_Type_Mat))) stop("SrvAgeComps_pop_Type is returning an NA. Did you update the year range of SrvAgeComps_pop_Type?")

  # Survey Length Composition Options ---------------------------------------

  comp_srvlen_like_vals <- vector()
  for(f in 1:input_list$data$n_srv_fleets) {
    if(SrvLenComps_LikeType[f] == 'none') comp_srvlen_like_vals <- c(comp_srvlen_like_vals, 999)
    if(SrvLenComps_LikeType[f] == "Multinomial") comp_srvlen_like_vals <- c(comp_srvlen_like_vals, 0)
    if(SrvLenComps_LikeType[f] == "Dirichlet-Multinomial") comp_srvlen_like_vals <- c(comp_srvlen_like_vals, 1)
    if(SrvLenComps_LikeType[f] == "iid-Logistic-Normal") comp_srvlen_like_vals <- c(comp_srvlen_like_vals, 2)
    if(SrvLenComps_LikeType[f] == "1d-Logistic-Normal") comp_srvlen_like_vals <- c(comp_srvlen_like_vals, 3)
    if(SrvLenComps_LikeType[f] == "2d-Logistic-Normal") comp_srvlen_like_vals <- c(comp_srvlen_like_vals, 4)
    collect_message(paste("Survey Length Composition Likelihoods", "for survey fleet", f, "specified as:" , SrvLenComps_LikeType[f]))
  } # end f loop

  SrvLenComps_Type_Mat <- array(NA, dim = c(length(input_list$data$years), input_list$data$n_srv_fleets))
  for(i in 1:length(SrvLenComps_Type)) {

    # Extract out components from list
    tmp <- SrvLenComps_Type[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    comps_type_tmp <- tmp_vec[1] # get composition type
    fleet <- as.numeric(tmp_vec[5]) # extract fleet index

    # Checking character string
    if(!comps_type_tmp %in% c("agg", "spltRspltS", "spltRjntS", 'none')) stop("SrvLenComps_Type not specified correctly. Must be one of: agg, spltRspltS, spltRjntS, none")
    if(!fleet %in% c(1:input_list$data$n_srv_fleets)) stop("Invalid fleet specified for SrvLenComps_Type This needs to be specified as CompType_Year_x-y_Fleet_x")

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
      if(comp_srvlen_like_vals[fleet] == 4) stop("Length composition likelihood specified as 2d-Logistic-Normal, but composition type is aggregated. This is not valid.")
      comps_type_val <- 0
    }
    if(comps_type_tmp == "spltRspltS") comps_type_val <- 1
    if(comps_type_tmp == "spltRjntS") comps_type_val <- 2
    if(comps_type_tmp == "none") comps_type_val <- 999

    # input into matrix
    SrvLenComps_Type_Mat[years,fleet] <- comps_type_val
  } # end i

  if(any(is.na(SrvLenComps_Type_Mat))) stop("SrvLenComps_Type_Mat is returning an NA. Did you update the year range of SrvLenComps_Type_Mat?")

  # Specifying composition likelihood for population-specific data
  comp_srvlen_pop_like_vals <- vector()
  for(f in 1:input_list$data$n_srv_fleets) {
    if(SrvLenComps_pop_LikeType[f] == 'none') comp_srvlen_pop_like_vals <- c(comp_srvlen_pop_like_vals, 999)
    if(SrvLenComps_pop_LikeType[f] == "Multinomial") comp_srvlen_pop_like_vals <- c(comp_srvlen_pop_like_vals, 0)
    if(SrvLenComps_pop_LikeType[f] == "Dirichlet-Multinomial") comp_srvlen_pop_like_vals <- c(comp_srvlen_pop_like_vals, 1)
    if(SrvLenComps_pop_LikeType[f] == "iid-Logistic-Normal") comp_srvlen_pop_like_vals <- c(comp_srvlen_pop_like_vals, 2)
    if(SrvLenComps_pop_LikeType[f] == "1d-Logistic-Normal") comp_srvlen_pop_like_vals <- c(comp_srvlen_pop_like_vals, 3)
    if(SrvLenComps_pop_LikeType[f] == "2d-Logistic-Normal") comp_srvlen_pop_like_vals <- c(comp_srvlen_pop_like_vals, 4)
    collect_message(paste("Population Survey Length Composition Likelihoods", "for survey fleet", f, "specified as:" , SrvLenComps_pop_LikeType[f]))
  } # end f loop

  # Specifying composition type
  SrvLenComps_pop_Type_Mat <- array(NA, dim = c(length(input_list$data$years), input_list$data$n_srv_fleets))
  for(i in 1:length(SrvLenComps_pop_Type)) {

    # Extract out components from list
    tmp <- SrvLenComps_pop_Type[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    comps_type_tmp <- tmp_vec[1] # get composition type
    fleet <- as.numeric(tmp_vec[5]) # extract fleet index

    # Checking character string
    if(!comps_type_tmp %in% c("agg", "spltRspltS", "spltRjntS", 'none')) stop("SrvLenComps_pop_Type not specified correctly. Must be one of: agg, spltRspltS, spltRjntS, none")
    if(!fleet %in% c(1:input_list$data$n_srv_fleets)) stop("Invalid fleet specified for SrvLenComps_pop_Type. This needs to be specified as CompType_Year_x-y_Fleet_x")

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
      if(comp_srvlen_pop_like_vals[fleet] == 4) stop("Population Len composition likelihood specified as 2d-Logistic-Normal, but composition type is aggregated. This is not valid.")
      comps_type_val <- 0
    }
    if(comps_type_tmp == "spltRspltS") comps_type_val <- 1
    if(comps_type_tmp == "spltRjntS") comps_type_val <- 2
    if(comps_type_tmp == "none") comps_type_val <- 999

    # input into matrix
    SrvLenComps_pop_Type_Mat[years,fleet] <- comps_type_val
  } # end i

  if(any(is.na(SrvLenComps_pop_Type_Mat))) stop("SrvLenComps_pop_Type is returning an NA. Did you update the year range of SrvLenComps_pop_Type?")

  # ISS Munging -------------------------------------------------------------

  # Survey Ages
  if(is.null(ISS_SrvAgeComps)) {
    collect_message("No ISS is specified for SrvAgeComps. ISS weighting is calculated by summing up values from ObsSrvAgeComps each year")
    ISS_SrvAgeComps <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_srv_fleets))
    for(y in 1:length(input_list$data$years)) {
      for(f in 1:input_list$data$n_srv_fleets) {
        for(seas in 1:input_list$data$n_seas) {
          # if aggregated across sexes and regions (0)
          if(SrvAgeComps_Type_Mat[y,f] == 0) ISS_SrvAgeComps[1,y,seas,1,f] <- sum(ObsSrvAgeComps[,y,seas,,,f])
          # if split by region and sex
          if(SrvAgeComps_Type_Mat[y,f] == 1) ISS_SrvAgeComps[,y,seas,,f] <- apply(ObsSrvAgeComps[,y,seas,,,f, drop = FALSE], c(1,4), sum)
          # if split by region, joint by sex
          if(SrvAgeComps_Type_Mat[y,f] == 2) ISS_SrvAgeComps[,y,seas,1,f] <- apply(ObsSrvAgeComps[,y,seas,,,f, drop = FALSE], 1, sum)
        } # end seas loop
      } # end f loop
    } # end y loop
  }

  # Survey Lengths
  if(is.null(ISS_SrvLenComps)) {
    collect_message("No ISS is specified for SrvLenComps. ISS weighting is calculated by summing up values from ObsSrvLenComps each year")
    ISS_SrvLenComps <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_srv_fleets))
    for(y in 1:length(input_list$data$years)) {
      for(f in 1:input_list$data$n_srv_fleets) {
        for(seas in 1:input_list$data$n_seas) {
          # if aggregated across sexes and regions (0)
          if(SrvLenComps_Type_Mat[y,f] == 0) ISS_SrvLenComps[1,y,seas,1,f] <- sum(ObsSrvLenComps[,y,seas,,,f])
          # if split by region and sex
          if(SrvLenComps_Type_Mat[y,f] == 1) ISS_SrvLenComps[,y,seas,,f] <- apply(ObsSrvLenComps[,y,seas,,,f, drop = FALSE], c(1,4), sum)
          # if split by region, joint by sex
          if(SrvLenComps_Type_Mat[y,f] == 2) ISS_SrvLenComps[,y,seas,1,f] <- apply(ObsSrvLenComps[,y,seas,,,f, drop = FALSE], 1, sum)
        } # end seas loop
      } # end f loop
    } # end y loop
  }

  # Srvery Ages
  if(is.null(ISS_SrvAgeComps_pop)) {
    collect_message("No ISS is specified for pop_SrvAgeComps. ISS weighting is calculated by summing up values from ObsSrvAgeComps_pop each year")
    ISS_SrvAgeComps_pop <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_srv_fleets))
    for(p in 1:input_list$data$n_pop) {
      for(y in 1:length(input_list$data$years)) {
        for(f in 1:input_list$data$n_srv_fleets) {
          for(seas in 1:input_list$data$n_seas) {
            # if aggregated across sexes and regions (0) or joint across sexes
            if(SrvAgeComps_pop_Type_Mat[y,f] == 0) ISS_SrvAgeComps_pop[p,1,y,seas,1,f] <- sum(ObsSrvAgeComps_pop[p,,y,seas,,,f])
            # if split by region and sex
            if(SrvAgeComps_pop_Type_Mat[y,f] == 1) ISS_SrvAgeComps_pop[p,,y,seas,,f] <- apply(ObsSrvAgeComps_pop[p,,y,seas,,,f, drop = FALSE], c(2,5), sum)
            # if split by region, joint by sex
            if(SrvAgeComps_pop_Type_Mat[y,f] == 2) ISS_SrvAgeComps_pop[p,,y,seas,1,f] <- apply(ObsSrvAgeComps_pop[p,,y,seas,,,f, drop = FALSE], 2, sum)
          } # end seas loop
        } # end f loop
      } # end y loop
    } # end p loop
  }

  # Srvery Lengths
  if(is.null(ISS_SrvLenComps_pop)) {
    collect_message("No ISS is specified for pop_SrvLenComps. ISS weighting is calculated by summing up values from ObsSrvLenComps_pop each year")
    ISS_SrvLenComps_pop <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_srv_fleets))
    for(p in 1:input_list$data$n_pop) {
      for(y in 1:length(input_list$data$years)) {
        for(f in 1:input_list$data$n_srv_fleets) {
          for(seas in 1:input_list$data$n_seas) {
            # if aggregated across sexes and regions (0)
            if(SrvLenComps_pop_Type_Mat[y,f] == 0) ISS_SrvLenComps_pop[p,1,y,seas,1,f] <- sum(ObsSrvLenComps_pop[p,,y,seas,,,f])
            # if split by region and sex
            if(SrvLenComps_pop_Type_Mat[y,f] == 1) ISS_SrvLenComps_pop[p,,y,seas,,f] <- apply(ObsSrvLenComps_pop[p,,y,seas,,,f, drop = FALSE], c(2,5), sum)
            # if split by region, joint by sex
            if(SrvLenComps_pop_Type_Mat[y,f] == 2) ISS_SrvLenComps_pop[p,,y,seas,1,f] <- apply(ObsSrvLenComps_pop[p,,y,seas,,,f, drop = FALSE], 2, sum)
          } # end seas loop
        } # end f loop
      } # end y loop
    } # end p loop
  }


  # Populate Data List ------------------------------------------------------

  input_list$data$ISS_SrvAgeComps <- ISS_SrvAgeComps
  input_list$data$ISS_SrvLenComps <- ISS_SrvLenComps
  input_list$data$ISS_SrvAgeComps_pop <- ISS_SrvAgeComps_pop
  input_list$data$ISS_SrvLenComps_pop <- ISS_SrvLenComps_pop
  input_list$data$ObsSrvIdx <- ObsSrvIdx
  input_list$data$ObsSrvIdx_SE <- ObsSrvIdx_SE
  input_list$data$UseSrvIdx <- UseSrvIdx
  input_list$data$ObsSrvIdx_pop <- ObsSrvIdx_pop
  input_list$data$ObsSrvIdx_pop_SE <- ObsSrvIdx_pop_SE
  input_list$data$UseSrvIdx_pop <- UseSrvIdx_pop
  input_list$data$ObsSrvAgeComps <- ObsSrvAgeComps
  input_list$data$UseSrvAgeComps <- UseSrvAgeComps
  input_list$data$ObsSrvLenComps <- ObsSrvLenComps
  input_list$data$UseSrvLenComps <- UseSrvLenComps
  input_list$data$ObsSrvAgeComps_pop <- ObsSrvAgeComps_pop
  input_list$data$UseSrvAgeComps_pop <- UseSrvAgeComps_pop
  input_list$data$ObsSrvLenComps_pop <- ObsSrvLenComps_pop
  input_list$data$UseSrvLenComps_pop <- UseSrvLenComps_pop
  input_list$data$SrvAgeComps_LikeType <- comp_srvage_like_vals
  input_list$data$SrvLenComps_LikeType <- comp_srvlen_like_vals
  input_list$data$SrvAgeComps_pop_LikeType <- comp_srvage_pop_like_vals
  input_list$data$SrvLenComps_pop_LikeType <- comp_srvlen_pop_like_vals
  input_list$data$SrvAgeComps_Type <- SrvAgeComps_Type_Mat
  input_list$data$SrvLenComps_Type <- SrvLenComps_Type_Mat
  input_list$data$srv_idx_type <- srv_idx_type_vals
  input_list$data$SrvAgeComps_pop_Type <- SrvAgeComps_pop_Type_Mat
  input_list$data$SrvLenComps_pop_Type <- SrvLenComps_pop_Type_Mat

  # Populate Parameter List -------------------------------------------------

  # Dispersion parameters for the survey age comps
  if("ln_SrvAge_theta" %in% names(starting_values)) input_list$par$ln_SrvAge_theta <- starting_values$ln_SrvAge_theta
  else input_list$par$ln_SrvAge_theta <- array(0, dim = c(input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_srv_fleets))

  # logistic normal correlation parameters for survey age comps
  if("SrvAge_corr_pars" %in% names(starting_values)) input_list$par$SrvAge_corr_pars <- starting_values$SrvAge_corr_pars
  else input_list$par$SrvAge_corr_pars <- array(0.01, dim = c(input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_srv_fleets, 2))

  # aggregated
  if("ln_SrvAge_theta_agg" %in% names(starting_values)) input_list$par$ln_SrvAge_theta_agg <- starting_values$ln_SrvAge_theta_agg
  else input_list$par$ln_SrvAge_theta_agg <- array(0, dim = c(input_list$data$n_srv_fleets))

  # aggregated correlation parameters
  if("SrvAge_corr_pars_agg" %in% names(starting_values)) input_list$par$SrvAge_corr_pars_agg <- starting_values$SrvAge_corr_pars_agg
  else input_list$par$SrvAge_corr_pars_agg <- array(0.01, dim = c(input_list$data$n_srv_fleets))

  # Dispersion parameters for survey length comps
  if("ln_SrvLen_theta" %in% names(starting_values)) input_list$par$ln_SrvLen_theta <- starting_values$ln_SrvLen_theta
  else input_list$par$ln_SrvLen_theta <- array(0, dim = c(input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_srv_fleets))

  # logistic normal correlation parameters for survey length comps
  if("SrvLen_corr_pars" %in% names(starting_values)) input_list$par$SrvLen_corr_pars <- starting_values$SrvLen_corr_pars
  else input_list$par$SrvLen_corr_pars <- array(0.01, dim = c(input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_srv_fleets, 2))

  # aggregated
  if("ln_SrvLen_theta_agg" %in% names(starting_values)) input_list$par$ln_SrvLen_theta_agg <- starting_values$ln_SrvLen_theta_agg
  else input_list$par$ln_SrvLen_theta_agg <- array(0, dim = c(input_list$data$n_srv_fleets))

  if("SrvLen_corr_pars_agg" %in% names(starting_values)) input_list$par$SrvLen_corr_pars_agg <- starting_values$SrvLen_corr_pars_agg
  else input_list$par$SrvLen_corr_pars_agg <- array(0.01, dim = c(input_list$data$n_srv_fleets))

  # Dispersion parameters for the population survey age comps
  if("ln_SrvAge_pop_theta" %in% names(starting_values)) input_list$par$ln_SrvAge_pop_theta <- starting_values$ln_SrvAge_pop_theta
  else input_list$par$ln_SrvAge_pop_theta <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_srv_fleets))

  # logistic normal correlation parameters for population survey age comps
  if("SrvAge_pop_corr_pars" %in% names(starting_values)) input_list$par$SrvAge_pop_corr_pars <- starting_values$SrvAge_pop_corr_pars
  else input_list$par$SrvAge_pop_corr_pars <- array(0.01, dim = c(input_list$data$n_pop,input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_srv_fleets, 2))

  # aggregated population pars
  if("ln_SrvAge_pop_theta_agg" %in% names(starting_values)) input_list$par$ln_SrvAge_pop_theta_agg <- starting_values$ln_SrvAge_pop_theta_agg
  else input_list$par$ln_SrvAge_pop_theta_agg <- array(0, dim = c(input_list$data$n_pop,input_list$data$n_srv_fleets))

  # aggregated population correlation parameters
  if("SrvAge_pop_corr_pars_agg" %in% names(starting_values)) input_list$par$SrvAge_pop_corr_pars_agg <- starting_values$SrvAge_pop_corr_pars_agg
  else input_list$par$SrvAge_pop_corr_pars_agg <- array(0.01, dim = c(input_list$data$n_pop,input_list$data$n_srv_fleets))

  # Dispersion parameters for population survey length comps
  if("ln_SrvLen_pop_theta" %in% names(starting_values)) input_list$par$ln_SrvLen_pop_theta <- starting_values$ln_SrvLen_pop_theta
  else input_list$par$ln_SrvLen_pop_theta <- array(0, dim = c(input_list$data$n_pop,input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_srv_fleets))

  # logistic normal correlation parameters for population survey length comps
  if("SrvLen_pop_corr_pars" %in% names(starting_values)) input_list$par$SrvLen_pop_corr_pars <- starting_values$SrvLen_pop_corr_pars
  else input_list$par$SrvLen_pop_corr_pars <- array(0.01, dim = c(input_list$data$n_pop,input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_srv_fleets, 2))

  # aggregated population pars
  if("ln_SrvLen_pop_theta_agg" %in% names(starting_values)) input_list$par$ln_SrvLen_pop_theta_agg <- starting_values$ln_SrvLen_pop_theta_agg
  else input_list$par$ln_SrvLen_pop_theta_agg <- array(0, dim = c(input_list$data$n_pop,input_list$data$n_srv_fleets))

  if("SrvLen_pop_corr_pars_agg" %in% names(starting_values)) input_list$par$SrvLen_pop_corr_pars_agg <- starting_values$SrvLen_pop_corr_pars_agg
  else input_list$par$SrvLen_pop_corr_pars_agg <- array(0.01, dim = c(input_list$data$n_pop, input_list$data$n_srv_fleets))


  # Mapping Options ---------------------------------------------------------

  input_list <- do_comp_theta_mapping(input_list, comp_prefix = "SrvAge", fleet_field = "n_srv_fleets")
  input_list <- do_comp_theta_mapping(input_list, comp_prefix = "SrvLen", fleet_field = "n_srv_fleets")
  input_list <- do_comp_corr_pars_mapping(input_list, comp_prefix = "SrvAge", fleet_field = "n_srv_fleets")
  input_list <- do_comp_corr_pars_mapping(input_list, comp_prefix = "SrvLen", fleet_field = "n_srv_fleets")

  input_list <- do_comp_theta_mapping(input_list, comp_prefix = "SrvAge", has_pop = TRUE, fleet_field = "n_srv_fleets")
  input_list <- do_comp_theta_mapping(input_list, comp_prefix = "SrvLen", has_pop = TRUE, fleet_field = "n_srv_fleets")
  input_list <- do_comp_corr_pars_mapping(input_list, comp_prefix = "SrvAge", has_pop = TRUE, fleet_field = "n_srv_fleets")
  input_list <- do_comp_corr_pars_mapping(input_list, comp_prefix = "SrvLen", has_pop = TRUE, fleet_field = "n_srv_fleets")

  # Print Messages ----------------------------------------------------------
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)
}

#' Set up survey selectivity and catchability specifications
#'
#' Configures all survey selectivity and catchability components of the
#' estimation model: continuous and blocked time-varying selectivity,
#' selectivity functional forms, catchability blocks and optional
#' environmental covariate effects, process error and deviation mapping, and
#' selectivity/catchability priors. Delegates parameter mapping to four
#' internal helpers (\code{\link{do_fixed_sel_pars_mapping}},
#' \code{\link{do_q_mapping}}, \code{\link{do_sel_pe_pars_mapping}},
#' \code{\link{do_sel_devs_mapping}}). Must be called after
#' \code{\link{Setup_Mod_SrvIdx_and_Comps}} and before model compilation.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, \code{$map},
#'   and \code{$verbose} sublists, as returned by upstream setup functions.
#'   \code{$data$srv_selex_type} must already be set by
#'   \code{\link{Setup_Mod_Biologicals}}.
#'
#' @param srv_sel_model Character vector specifying the selectivity functional
#'   form per fleet, and optionally per time block. Each element follows one
#'   of:
#'   \describe{
#'     \item{\code{"<model>_Fleet_x"}}{Single form applied across all years for
#'       fleet \code{x}.}
#'     \item{\code{"<model>_Fleet_x_Block_k"}}{Form applied only to block
#'       \code{k} for fleet \code{x}, as defined in \code{srv_sel_blocks}.
#'       Required when multiple blocks are defined for a fleet.}
#'   }
#'   Available models (see the model equations vignette for parameterisations):
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
#'       defined via \code{srv_sel_blocks}). An optional \code{_SelStyr_<year>}
#'       suffix (a calendar year within the block) restricts the actual spline
#'       fit to \code{SelStyr}:block-end; years within the block before
#'       \code{SelStyr} are held constant at the \code{SelStyr} year's fitted
#'       curve, rather than fitting the surface over the whole block. An
#'       optional \code{_NSelBins_<n>} suffix restricts the actual spline fit
#'       to the first \code{n} bins (ages or lengths, per \code{srv_selex_type});
#'       bins beyond \code{n} are held constant at the last fitted bin's
#'       curve.}
#'   }
#'   No default; must be provided.
#' @param srv_fixed_sel_pars_spec Character vector \code{[n_srv_fleets]}.
#'   Sharing structure for fixed-effect selectivity parameters. See
#'   \code{\link{do_fixed_sel_pars_mapping}} for full option descriptions.
#'   No default; must be provided.
#' @param srv_sel_blocks Character vector defining discrete time blocks for
#'   survey selectivity. Each element follows \code{"Block_k_Year_a-b_Fleet_x"}
#'   or \code{"none_Fleet_x"} (constant selectivity). Use \code{"terminal"} in
#'   place of the end year to extend to the final model year. Parsed into an
#'   internal array \code{[n_regions × n_years × n_srv_fleets]}. Blocked and
#'   continuous time-varying selectivity are mutually exclusive for a given
#'   fleet. Default: \code{"none_Fleet_x"} for each fleet.
#'
#' @param cont_tv_srv_sel Character vector defining the continuous
#'   time-variation form per fleet. Each element follows
#'   \code{"<type>_Fleet_x"}. Options:
#'   \describe{
#'     \item{\code{"none"}}{No continuous time variation (default).}
#'     \item{\code{"iid"}}{IID deviations across years.}
#'     \item{\code{"rw"}}{Random walk in time.}
#'     \item{\code{"3dmarg"}}{3D marginal GMRF (age × year × cohort).}
#'     \item{\code{"3dcond"}}{3D conditional GMRF.}
#'     \item{\code{"2dar1"}}{2D AR1 (bin × year).}
#'   }
#'   When any fleet uses a non-\code{"none"} type, both
#'   \code{srvsel_pe_pars_spec} and \code{srv_sel_devs_spec} must be
#'   specified. Default: \code{"none_Fleet_x"} for each fleet.
#' @param srvsel_pe_pars_spec Character vector \code{[n_srv_fleets]} or
#'   \code{NULL}. Sharing structure for process error hyperparameters. See
#'   \code{\link{do_sel_pe_pars_mapping}} for full option descriptions.
#'   Default \code{NULL}.
#' @param srv_sel_devs_spec Character vector \code{[n_srv_fleets]} or
#'   \code{NULL}. Sharing structure for selectivity deviation time series.
#'   See \code{\link{do_sel_devs_mapping}} for full option descriptions.
#'   Default \code{NULL}.
#' @param srvsel_devs_shared_bins List of integer vectors defining bin groups
#'   for age/length-sharing of deviations under semi-parametric forms (e.g.,
#'   \code{list(1:5, 6:10, 11:30)}). Required when \code{srv_sel_devs_spec}
#'   includes any \code{"est_shared_b"} variant. Default \code{NULL}.
#' @param corr_opt_semipar Character vector \code{[n_srv_fleets]} or
#'   \code{NULL}. Specifies correlation components to suppress for 3D GMRF or
#'   2D AR1 forms. See \code{\link{do_sel_pe_pars_mapping}} for valid
#'   values. Default \code{NULL}.
#'
#' @param srv_q_blocks Character vector defining discrete time blocks for
#'   survey catchability. Same format as \code{srv_sel_blocks}:
#'   \code{"Block_k_Year_a-b_Fleet_x"} or \code{"none_Fleet_x"}. Parsed into
#'   an array \code{[n_regions × n_years × n_srv_fleets]}. Default:
#'   \code{"none_Fleet_x"} for each fleet.
#' @param srv_q_spec Character vector \code{[n_srv_fleets]} or \code{NULL}.
#'   Sharing structure for catchability. See \code{\link{do_q_mapping}}
#'   for full option descriptions. Default \code{NULL}.
#' @param Use_srv_q_prior Integer (0/1). Whether log-normal priors are applied
#'   to survey catchability parameters. Default \code{0}.
#' @param srv_q_prior Data frame of catchability prior specifications. Required
#'   columns: \code{region}, \code{fleet}, \code{block}, \code{mu} (prior mean
#'   on natural scale), \code{sd} (prior SD on log scale). Each row specifies
#'   a \eqn{\log\text{N}(\log(\mu), \text{sd})} prior. Ignored when
#'   \code{Use_srv_q_prior = 0}. Default \code{NA}.
#' @param srv_q_formula Named list of R formulas specifying environmental
#'   covariate relationships for catchability per region–fleet combination.
#'   Names follow the convention \code{"Region_r_Fleet_f"}. Covariates must
#'   be present in \code{srv_q_cov_dat}. If \code{NULL}, no covariate effects
#'   are included. Default \code{NULL}.
#' @param srv_q_cov_dat Named list of numeric vectors (length = \code{n_years})
#'   containing covariate time series referenced in \code{srv_q_formula}.
#'   All vectors must be the same length and contain no missing values; set
#'   values to \code{0} for years when the survey is not active. If
#'   \code{NULL}, covariate effects are excluded. Default \code{NULL}.
#'
#' @param t_srv Survey timing fraction within a given year (annual models) or
#'   season (seasonal models), array
#'   \code{[n_regions × n_seas × n_srv_fleets]}. Default: \code{1}
#'   (end of period).
#' @param Use_srv_selex_prior Integer (0/1). Whether log-normal priors are
#'   applied to survey selectivity parameters. Default \code{0}.
#' @param srv_selex_prior Data frame of selectivity prior specifications.
#'   Required columns: \code{region}, \code{fleet}, \code{block}, \code{sex},
#'   \code{par}, \code{mu}, \code{sd}. Ignored when
#'   \code{Use_srv_selex_prior = 0}. Default \code{NULL}.
#'
#' @param ... Optional named starting values for selectivity and catchability
#'   parameters.
#' @param srv_selex_type Character. Whether survey selectivity type is 'age' or 'length' based. Default: \code{age}.
#' @param use_fixed_srv_sel Integer vector of length \code{n_srv_fleets}
#'   indicating whether survey selectivity is fixed (\code{1}) or estimated
#'   (\code{0}) for each survey index.
#'
#' @param srv_sel_input Array of fixed survey selectivity values used when
#'   \code{use_fixed_srv_sel == 1}. Dimensions:
#'   \code{[n_pop × n_regions × n_years × n_seas × n_bins × n_sexes × n_srv_fleets]}.
#'   Required whenever any survey has fixed selectivity specified.
#'
#' @param srv_sel_nonpar_est_bins Optional list defining bin groupings for
#'   non-parametric survey selectivity. Structure is
#'   \code{[[survey]][[block]]}, where each element is a list of integer vectors.
#'   Each vector defines a group of bins that share a single estimated
#'   selectivity parameter. Indices must correspond to the bin dimension
#'   defined by the survey selectivity type (age or length).
#'
#' @return The input \code{input_list} with selectivity and catchability
#'   configuration stored in \code{$data} (\code{cont_tv_srv_sel}, \code{srv_sel_blocks},
#'   \code{srv_sel_model}, \code{srv_q_blocks}, \code{srv_q_prior},
#'   \code{Use_srv_q_prior}, \code{do_srv_q_cov}, \code{srv_q_cov},
#'   \code{Use_srv_selex_prior}, \code{srv_selex_prior}, \code{t_srv});
#'   starting values in \code{$par} for \code{srv_fixed_sel_pars},
#'   \code{ln_srv_q}, \code{srvsel_pe_pars}, \code{ln_srvsel_devs}, and
#'   \code{srv_q_coeff}; and factor maps in \code{$map} for all five
#'   parameter arrays plus \code{srv_q_coeff}.
#'
#' @export Setup_Mod_Srvsel_and_Q
#' @importFrom stringr str_detect
#' @family Model Setup
Setup_Mod_Srvsel_and_Q <- function(input_list,
                                   cont_tv_srv_sel = paste("none_Fleet_", 1:input_list$data$n_srv_fleets, sep = ''),
                                   srv_sel_blocks = paste("none_Fleet_", 1:input_list$data$n_srv_fleets, sep = ''),
                                   srv_sel_model,
                                   Use_srv_q_prior = 0,
                                   srv_q_prior = NA,
                                   srv_q_blocks = paste("none_Fleet_", 1:input_list$data$n_srv_fleets, sep = ''),
                                   srvsel_pe_pars_spec = NULL,
                                   srv_fixed_sel_pars_spec,
                                   srv_q_spec = NULL,
                                   srv_sel_devs_spec = NULL,
                                   corr_opt_semipar = NULL,
                                   srv_q_formula = NULL,
                                   srv_q_cov_dat = NULL,
                                   Use_srv_selex_prior = 0,
                                   srv_selex_prior = NULL,
                                   t_srv = array(1, dim = c(input_list$data$n_regions, input_list$data$n_seas, input_list$data$n_srv_fleets)),
                                   srvsel_devs_shared_bins = NULL,
                                   srv_selex_type = 'age',
                                   use_fixed_srv_sel = rep(0, input_list$data$n_srv_fleets),
                                   srv_sel_input = NULL,
                                   srv_sel_nonpar_est_bins = NULL,
                                   ...
                                   ) {

  messages_list <<- character(0) # string to attach to for printing messages
  starting_values <- list(...)
  if(input_list$store_config) input_list$config$Setup_Mod_Srvsel_and_Q <- mget(names(formals()))[-1]

  # Input Validation --------------------------------------------------------

  # Selectivity
  # Continuous Selectivity Deviations
  if(!is.null(srvsel_pe_pars_spec)) if(length(srvsel_pe_pars_spec) != input_list$data$n_srv_fleets) stop("srvsel_pe_pars_spec is not length n_srv_fleets")
  if(!is.null(srv_sel_devs_spec)) if(length(srv_sel_devs_spec) != input_list$data$n_srv_fleets) stop("srv_sel_devs_spec is not length n_srv_fleets")
  if(!is.null(corr_opt_semipar)) if(length(corr_opt_semipar) != input_list$data$n_srv_fleets) stop("corr_opt_semipar is not length n_srv_fleets")

  # Catchability Priors
  if(!Use_srv_q_prior %in% c(0,1)) stop("Values for Use_srv_q_prior are not valid. They are == 0 (don't use prior), or == 1 (use prior)")
  # Checking catchability priors
  if(Use_srv_q_prior == 1) {
    required_cols <- c("region", "fleet", "block", "mu", "sd")
    missing_cols <- setdiff(required_cols, names(srv_q_prior))
    if(length(missing_cols) > 0) {
      stop("srv_q_prior is missing required columns: ", paste(missing_cols, collapse = ", "))
    }
  }
  collect_message("Survey Catchability priors are: ", ifelse(Use_srv_q_prior == 0, "Not Used", "Used"))

  # Selectivity Priors
  if(!Use_srv_selex_prior %in% c(0,1)) stop("Values for Use_srv_selex_prior are not valid. They are == 0 (don't use prior), or == 1 (use prior)")
  # Checking selectivity priors
  if(Use_srv_selex_prior == 1) {
    required_cols <- c("region", "fleet", "block", "sex", "par", "mu", "sd")
    missing_cols <- setdiff(required_cols, names(srv_selex_prior))
    if(length(missing_cols) > 0) {
      stop("srv_selex_prior is missing required columns: ", paste(missing_cols, collapse = ", "))
    }
  }
  collect_message("Survey Selectivity priors are: ", ifelse(Use_srv_selex_prior == 0, "Not Used", "Used"))

  if(any(use_fixed_srv_sel == 1) && is.null(srv_sel_input)) stop("srv_sel_input is NULL, please provide an input array.")
  if(any(use_fixed_srv_sel == 1) && srv_selex_type == 'age') check_data_dimensions(srv_sel_input, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_ages = length(input_list$data$ages), n_sexes = input_list$data$n_sexes, n_srv_fleets = input_list$data$n_srv_fleets, what = 'srv_sel_input_age')
  if(any(use_fixed_srv_sel == 1) && srv_selex_type == 'length') check_data_dimensions(srv_sel_input, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_lens = length(input_list$data$lens), n_sexes = input_list$data$n_sexes, n_srv_fleets = input_list$data$n_srv_fleets, what = 'srv_sel_input_len')

  # Selectivity Options -----------------------------------------------------
  # Age based selectivity
  if(srv_selex_type == 'age') {
    srv_selex_type <- 0
    bins <- length(input_list$data$ages)
    collect_message("Survey Selectivity is aged-based.")
  } # if age based

  # Length based selectivity
  if(srv_selex_type == 'length') {
    if(input_list$data$fit_lengths == 0) stop("Length composition data are not fit, but survey selectivity is length-based. This is not allowed. Please change to a valid option (either fit lengths or use age-based selectivity).")
    srv_selex_type <- 1
    bins <- length(input_list$data$lens)
    collect_message("Survey Selectivity is length-based")
  } # if length based

  # Continuous Time-Varying Selectivity Options -----------------------------
  # define for continuous time-varying selectivity
  cont_tv_srv_sel_mat <- array(NA, dim = c(input_list$data$n_regions, input_list$data$n_srv_fleets))
  cont_tv_map <- data.frame(type = c("none", "iid", "rw", "3dmarg", "3dcond", "2dar1"), num = c(0,1,2,3,4,5)) # set up values we map to

  for(i in 1:length(cont_tv_srv_sel)) {
    # Extract out components from list
    tmp <- cont_tv_srv_sel[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    cont_tv_type <- tmp_vec[1] # get continuous selex type
    fleet <- as.numeric(tmp_vec[3]) # extract fleet index

    # Validate options
    if(!fleet %in% c(1:input_list$data$n_srv_fleets)) stop("Invalid fleet specified for cont_tv_srv_sel This needs to be specified as timevarytype_Fleet_x")
    if(!cont_tv_type %in% c(cont_tv_map$type)) stop("cont_tv_srv_sel is not correctly specified. This needs to be one of these: none, iid, rw, 3dmarg, 3dcond, 2dar1 (the timevarytypes) and specified as timevarytype_Fleet_x")

    # Input options
    cont_tv_srv_sel_mat[,fleet] <- cont_tv_map$num[which(cont_tv_map$type == cont_tv_type)]
    collect_message("Continuous survey time-varying selectivity specified as: ", cont_tv_type, " for survey fleet ", fleet)
  }

  if(any(cont_tv_srv_sel_mat > 0) && is.null(srvsel_pe_pars_spec) && is.null(srv_sel_devs_spec)) stop("Continuous time-varying selectivity specified, but srvsel_pe_pars_spec and/or srv_sel_devs_spec is NULL (i.e., not specified)!")

  # Blocked Time-Varying Selectivity Options --------------------------------
  srv_sel_blocks_arr <- array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_srv_fleets))
  for(i in 1:length(srv_sel_blocks)) {

    # Extract out components from list
    tmp <- srv_sel_blocks[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))

    if(!tmp_vec[1] %in% c("none", "Block")) stop("Survey Selectivity Blocks not correctly specified. This should be either none_Fleet_x or Block_x_Year_x-y_Fleet_x")

    # extract out fleets if constant
    if(tmp_vec[1] == "none") {
      fleet <- as.numeric(tmp_vec[3]) # get fleet number
      srv_sel_blocks_arr[,,fleet] <- 1 # input only 1 survey time block
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

      srv_sel_blocks_arr[,years,fleet] <- block_val
    }
  }

  if(any(is.na(srv_sel_blocks_arr))) stop("Survey Selectivity Blocks are returning an NA. Did you update the year range of srv_sel_blocks?")
  for(f in 1:input_list$data$n_srv_fleets) collect_message(paste("Survey Selectivity Time Blocks for survey", f, "is specified at:", length(unique(srv_sel_blocks_arr[,,f]))))

  # Selectivity Functional Forms --------------------------------------------
  sel_map <- data.frame(sel = c('logist1', "gamma", "exponential", "logist2", "dbnrml", 'nonpar', 'asymplogist1', "asymplogist2", "bicubic"), num = c(0,1,2,3,4,5,6,7,8)) # set up values we can map to
  srv_sel_model_arr <- array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_srv_fleets))
  srv_sel_bicubic_binnodes_arr <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_srv_fleets)) # number of bin nodes, only set where srv_sel_model == 8 (bicubic)
  srv_sel_bicubic_yrnodes_arr <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_srv_fleets)) # number of year nodes, only set where srv_sel_model == 8 (bicubic)
  srv_sel_bicubic_selstyr_arr <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_srv_fleets)) # calendar year the bicubic surface is actually fit from (0 = block's own start year, i.e. no offset); years within the block before this are edge-held at this year's fitted curve
  srv_sel_bicubic_nselbins_arr <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_srv_fleets)) # number of bins (starting from the first) the bicubic surface is actually fit over (0 = all bins, i.e. no truncation); bins beyond this are held flat at the last fitted bin's value

  for(i in 1:length(srv_sel_model)) {

    # Extract out survey selectivity components from vector
    tmp_sel_form <- srv_sel_model[i]
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
        stop("srv_sel_model 'bicubic' entries must be specified as bicubic_Bin_<n_bin_nodes>_Yr_<n_yr_nodes>_Fleet_<f> or bicubic_Bin_<n_bin_nodes>_Yr_<n_yr_nodes>_Block_<b>_Fleet_<f>, optionally with _SelStyr_<year> and/or _NSelBins_<n>")
      tmp_n_bin_nodes <- suppressWarnings(as.numeric(tmp_sel_form_vec[bin_pos + 1]))
      tmp_n_yr_nodes <- suppressWarnings(as.numeric(tmp_sel_form_vec[yr_pos + 1]))
      tmp_fleet <- suppressWarnings(as.numeric(tmp_sel_form_vec[fleet_pos + 1]))
      tmp_block <- if(length(block_pos) == 1) suppressWarnings(as.numeric(tmp_sel_form_vec[block_pos + 1])) else NULL
      tmp_selstyr <- if(length(selstyr_pos) == 1) suppressWarnings(as.numeric(tmp_sel_form_vec[selstyr_pos + 1])) else 0
      tmp_nselbins <- if(length(nselbins_pos) == 1) suppressWarnings(as.numeric(tmp_sel_form_vec[nselbins_pos + 1])) else 0
      if(is.na(tmp_n_bin_nodes) || tmp_n_bin_nodes < 2) stop("bicubic srv_sel_model requires at least 2 bin nodes (n_bin_nodes >= 2)")
      if(is.na(tmp_n_yr_nodes) || tmp_n_yr_nodes < 1) stop("bicubic srv_sel_model requires at least 1 year node (n_yr_nodes >= 1). Use n_yr_nodes == 1 for a time-invariant bin-only spline.")
      if(length(selstyr_pos) == 1 && (is.na(tmp_selstyr) || !tmp_selstyr %in% input_list$data$years)) stop("bicubic srv_sel_model SelStyr must be a calendar year within the modeled years")
      if(length(nselbins_pos) == 1 && (is.na(tmp_nselbins) || tmp_nselbins < 2 || tmp_nselbins > bins)) stop("bicubic srv_sel_model NSelBins must be an integer between 2 and the total number of bins (ages or lengths)")
    } else {
      # get fleet index
      tmp_fleet <- if(length(tmp_sel_form_vec) == 3) as.numeric(tmp_sel_form_vec[3]) else as.numeric(tmp_sel_form_vec[5]) # fleet index changes if block is included in character vector
      # get block index
      tmp_block <- if(length(tmp_sel_form_vec) == 5) as.numeric(tmp_sel_form_vec[3]) else NULL
    }

    # validate options
    if(!sel_form %in% c(sel_map$sel)) stop("srv_sel_model is not correctly specified. This needs to be one of these: logist1, gamma, exponential, logist2, dbnrml, nonpar, asymplogist1, asymplogist2, bicubic (the seltypes) and specified as seltype_Fleet_x")
    if(!tmp_fleet %in% c(1:input_list$data$n_srv_fleets)) stop("Invalid fleet specified for srv_sel_model This needs to be specified as seltype_Fleet_x or seltype_Fleet_x_Block_x (if blocks are specified to change for a fleet)")

    # Input options
    if(is.null(tmp_block)) {
      srv_sel_model_arr[,,tmp_fleet] <- sel_map$num[which(sel_map$sel == sel_form)] # same selectivity form across blocks
      if(sel_form == "bicubic") {
        srv_sel_bicubic_binnodes_arr[,,tmp_fleet] <- tmp_n_bin_nodes
        srv_sel_bicubic_yrnodes_arr[,,tmp_fleet] <- tmp_n_yr_nodes
        srv_sel_bicubic_selstyr_arr[,,tmp_fleet] <- tmp_selstyr
        srv_sel_bicubic_nselbins_arr[,,tmp_fleet] <- tmp_nselbins
      }
    } else {
      srv_sel_model_arr[,which(srv_sel_blocks_arr[,,tmp_fleet] == tmp_block),tmp_fleet] <- sel_map$num[which(sel_map$sel == sel_form)]
      if(sel_form == "bicubic") {
        srv_sel_bicubic_binnodes_arr[,which(srv_sel_blocks_arr[,,tmp_fleet] == tmp_block),tmp_fleet] <- tmp_n_bin_nodes
        srv_sel_bicubic_yrnodes_arr[,which(srv_sel_blocks_arr[,,tmp_fleet] == tmp_block),tmp_fleet] <- tmp_n_yr_nodes
        srv_sel_bicubic_selstyr_arr[,which(srv_sel_blocks_arr[,,tmp_fleet] == tmp_block),tmp_fleet] <- tmp_selstyr
        srv_sel_bicubic_nselbins_arr[,which(srv_sel_blocks_arr[,,tmp_fleet] == tmp_block),tmp_fleet] <- tmp_nselbins
      }
    }
    rm(tmp_block) # remove tmp block to start next loop
    collect_message("Survey selectivity functional form specified as:", sel_form, " for survey fleet ", tmp_fleet)
  }

  # Validate that blocks and continuous time-variation aren't both specified for same fleet
  for(f in 1:input_list$data$n_srv_fleets) {
    has_blocks <- length(unique(srv_sel_blocks_arr[1,,f])) > 1
    has_cont_tv <- cont_tv_srv_sel_mat[1,f] != 0  # 0 = "none"
    if(has_blocks && has_cont_tv) {
      stop("Fleet ", f, " has both selectivity blocks and continuous time-varying selectivity specified. ",
           "These are mutually exclusive - choose one approach to time-variation.")
    }
  }

  # Blocked Catchability Options --------------------------------------------
  srv_q_blocks_arr <- array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_srv_fleets))
  for(i in 1:length(srv_q_blocks)) {

    # Extract out components from list
    tmp <- srv_q_blocks[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))

    # Vakudate option
    if(!tmp_vec[1] %in% c("none", "Block")) stop("Survey Catchability Blocks not correctly specified. This should be either none_Fleet_x or Block_x_Year_x-y_Fleet_x")

    # extract out fleets if constant
    if(tmp_vec[1] == "none") {
      fleet <- as.numeric(tmp_vec[3]) # get fleet number
      srv_q_blocks_arr[,,fleet] <- 1 # input only 1 survey catchability time block
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
      srv_q_blocks_arr[,years,fleet] <- block_val # input catchability time block
    }
  }
  if(any(is.na(srv_q_blocks_arr))) stop("Survey Catchability Blocks are returning an NA. Did you update the year range of srv_q_blocks?")
  for(f in 1:input_list$data$n_srv_fleets) collect_message(paste("Survey Catchability Time Blocks for survey", f, "is specified at:", length(unique(srv_q_blocks_arr[,,f]))))

  # Covariate Catchability Options ------------------------------------------
  if(!is.null(srv_q_cov_dat) && !is.null(srv_q_formula)) collect_message("Using covariates to predict survey catchability")

  # Figure out the total number of regression coefficients that could be estimated
  if(!is.null(srv_q_cov_dat) && !is.null(srv_q_formula)) {
    n_srv_q_cov <- max(sapply(names(srv_q_formula), function(key) { # sapply to extract names from formula
      tmp_formula <- srv_q_formula[[key]] # get formula
      var_names <- all.vars(tmp_formula) # get var names
      tmp_dat <- data.frame(srv_q_cov_dat[var_names]) # make dataframe
      ncol(model.matrix(tmp_formula, data = tmp_dat)) # figure out number of columns for formula (max number of coefficients to estimate)
    }))
  } else {
    do_srv_q_cov <- 0 # Indicator for whether covariates are included into survey catchability
    n_srv_q_cov <- 1 # dummy to initialize the array
  }

  # Catchability covariate containers
  srv_q_cov <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_srv_fleets, n_srv_q_cov)) # environmental time series
  srv_q_coeff <- array(0, dim = c(input_list$data$n_regions, input_list$data$n_srv_fleets, n_srv_q_cov)) # coefficients to be estimated
  map_srv_q_coeff <- array(NA, dim = c(input_list$data$n_regions, input_list$data$n_srv_fleets, n_srv_q_cov)) # coefficients to be mapped off

  # Loop through to map stuff off and populate containers
  if(!is.null(srv_q_cov_dat) && !is.null(srv_q_formula)) {

    # Validate covariate length
    cov_lengths <- lengths(srv_q_cov_dat)

    # Validate options
    # Check all covariates are the same length
    if (length(unique(cov_lengths)) != 1) stop("All covariates in 'srv_q_cov_dat' must have the same length. If some years are missing data, either impute some value, or set at 0 (if it is not used in the calculation).")
    # Check that length matches the model year structure
    if (unique(cov_lengths) != length(input_list$data$years)) stop(paste0("Covariate length mismatch: expected ",  length(input_list$data$years),  " years but got ", unique(cov_lengths),  "."))

    do_srv_q_cov <- 1 # Indicator for whether covariates are included into survey catchability
    coeff_counter <- 0 # setup counter for mapping

    for(r in 1:input_list$data$n_regions) {
      for(f in 1:input_list$data$n_srv_fleets) {

        # Get key to index
        key <- paste(paste("Region", r, sep = "_"), "_Fleet_", f, sep = "")
        # get temporary formula
        tmp_formula <- srv_q_formula[[key]]
        # extract variable names
        var_names <- all.vars(tmp_formula)
        if(length(var_names) == 0) next # skip if no variables
        # get environmental covariates from environmental data list, based on model formula
        tmp_dat <- data.frame(srv_q_cov_dat[var_names])
        # Generate design matrix
        tmp_design_mat <- model.matrix(tmp_formula, data = tmp_dat)
        # store covariate effects into container
        srv_q_cov[r,,f,1:ncol(tmp_design_mat)] <- tmp_design_mat

        # setup mapping - assign unique counter values for each coefficient
        for(i in 1:ncol(tmp_design_mat)) {
          coeff_counter <- coeff_counter + 1
          map_srv_q_coeff[r,f,i] <- coeff_counter
        } # end i loop

      } # end sf loop
    } # end r loop
  } # if using covariates

  # Populate Data List ------------------------------------------------------
  input_list$data$cont_tv_srv_sel <- cont_tv_srv_sel_mat
  input_list$data$srv_sel_blocks <- srv_sel_blocks_arr
  input_list$data$srv_sel_model <- srv_sel_model_arr
  input_list$data$srv_sel_bicubic_binnodes <- srv_sel_bicubic_binnodes_arr
  input_list$data$srv_sel_bicubic_yrnodes <- srv_sel_bicubic_yrnodes_arr
  input_list$data$srv_sel_bicubic_selstyr <- srv_sel_bicubic_selstyr_arr
  input_list$data$srv_sel_bicubic_nselbins <- srv_sel_bicubic_nselbins_arr
  input_list$data$srv_q_blocks <- srv_q_blocks_arr
  input_list$data$srv_q_prior <- srv_q_prior
  input_list$data$Use_srv_q_prior <- Use_srv_q_prior
  input_list$data$do_srv_q_cov <- do_srv_q_cov
  input_list$data$srv_q_cov <- srv_q_cov
  input_list$data$Use_srv_selex_prior <- Use_srv_selex_prior
  input_list$data$srv_selex_prior <- srv_selex_prior
  input_list$data$t_srv <- t_srv
  input_list$data$srv_selex_type <- srv_selex_type
  input_list$data$use_fixed_srv_sel <- use_fixed_srv_sel
  input_list$data$srv_sel_input <- srv_sel_input
  input_list$data$srvsel_devs_min_shared_bins <- if(!is.null(srvsel_devs_shared_bins)) unlist(lapply(srvsel_devs_shared_bins, min)) else 1:length(input_list$data$ages)

  # Populate Parameter List -------------------------------------------------
  # Figure out number of selectivity parameters for a given functional form
  unique_srvsel_vals <- unique(as.vector(input_list$data$srv_sel_model))
  sel_pars_vec <- vector() # create empty vector to populate

  for(i in 1:length(unique_srvsel_vals)) {
    if(unique_srvsel_vals[i] %in% c(2)) sel_pars_vec[i] <- 1 # exponential
    if(unique_srvsel_vals[i] %in% c(0,1,3)) sel_pars_vec[i] <- 2 # logistic or gamma
    if(unique_srvsel_vals[i] == 4) sel_pars_vec[i] <- 6 # double normal
    if(unique_srvsel_vals[i] == 5) sel_pars_vec[i] <- bins # non-parametric selex
    if(unique_srvsel_vals[i] %in% c(6,7)) sel_pars_vec[i] <- 3 # logistic selex w/ asymptote
    if(unique_srvsel_vals[i] == 8) sel_pars_vec[i] <- max(input_list$data$srv_sel_bicubic_binnodes * input_list$data$srv_sel_bicubic_yrnodes) # bicubic: flattened bin-node x year-node grid
  } # end i loop

  max_srvsel_blks <- max(apply(input_list$data$srv_sel_blocks, c(1,3), FUN = function(x) length(unique(x)))) # figure out maximum number of survey selectivity blocks for a given reigon and fleet

  # Bicubic spline interpolation weight matrices (bin node x year node grid), built here so they can be
  # threaded through SPoRC_rtmb.R alongside the flattened node parameters (see Get_Selex, Selex_Model == 8).
  # Padded with zeros to a common width across regions/blocks/fleets; padding is harmless because unused
  # (zero-weight) columns/rows never contribute to the resulting selectivity (see Get_Selex documentation).
  has_bicubic_srv_sel <- any(input_list$data$srv_sel_model == 8)
  max_bin_nodes_bicubic <- if(has_bicubic_srv_sel) max(input_list$data$srv_sel_bicubic_binnodes) else 1
  max_yr_nodes_bicubic <- if(has_bicubic_srv_sel) max(input_list$data$srv_sel_bicubic_yrnodes) else 1
  n_yrs_total_bicubic <- length(input_list$data$years) + input_list$data$n_proj_yrs_devs

  srv_sel_bicubic_Wbin <- array(0, dim = c(input_list$data$n_regions, bins, max_bin_nodes_bicubic, max_srvsel_blks, input_list$data$n_srv_fleets))
  srv_sel_bicubic_Wyr <- array(0, dim = c(input_list$data$n_regions, n_yrs_total_bicubic, max_yr_nodes_bicubic, max_srvsel_blks, input_list$data$n_srv_fleets))

  if(has_bicubic_srv_sel) {
    for(f in 1:input_list$data$n_srv_fleets) {
      for(r in 1:input_list$data$n_regions) {

        srvsel_blocks_tmp <- unique(as.vector(input_list$data$srv_sel_blocks[r,,f]))

        for(b in 1:length(srvsel_blocks_tmp)) {

          block_years <- which(input_list$data$srv_sel_blocks[r,,f] == srvsel_blocks_tmp[b])
          if(unique(input_list$data$srv_sel_model[r, block_years, f]) != 8) next # only bicubic blocks need weight matrices

          n_bin_nodes_this <- unique(input_list$data$srv_sel_bicubic_binnodes[r, block_years, f])
          n_yr_nodes_this <- unique(input_list$data$srv_sel_bicubic_yrnodes[r, block_years, f])

          # Bin dimension: nodes evenly spaced over [0,1]. By default (NSelBins unset, i.e. 0) the
          # spline is evaluated over all bins, as before. When NSelBins is set, the spline surface is only actually fit over the first NSelBins bins;
          # bins beyond that are edge-held at the last fitted bin's weights ("plateau").
          nselbins_this <- unique(input_list$data$srv_sel_bicubic_nselbins[r, block_years, f])
          n_fit_bins <- if(nselbins_this == 0) bins else nselbins_this

          bin_nodes_scaled <- seq(0, 1, length.out = n_bin_nodes_this)
          fit_bin_scaled <- seq(0, 1, length.out = n_fit_bins)
          Wbin_fit <- Get_Natural_Cubic_Spline_Weights(bin_nodes_scaled, fit_bin_scaled)

          Wbin_this <- matrix(0, nrow = bins, ncol = n_bin_nodes_this)
          Wbin_this[1:n_fit_bins, ] <- Wbin_fit
          if(n_fit_bins < bins) Wbin_this[(n_fit_bins + 1):bins, ] <- matrix(Wbin_fit[nrow(Wbin_fit), ], nrow = bins - n_fit_bins, ncol = n_bin_nodes_this, byrow = TRUE)

          srv_sel_bicubic_Wbin[r, , 1:n_bin_nodes_this, b, f] <- Wbin_this

          # Year dimension: nodes evenly spaced over the block's own contiguous fit range. By default
          # (SelStyr unset, i.e. 0) the fit range is the whole block, as before. When SelStyr is set , only years from SelStyr through the block's end are
          # actually spline-fit; years within the block before SelStyr are edge-held at the SelStyr
          # row's weights ("previous years are filled"). Years outside the block entirely (before it,
          # after it, and any projection years, since projections reuse the terminal modeled year's
          # block) hold the boundary node weights constant, which for a spline evaluated exactly at
          # its first/last node reduces to full weight on that node.
          selstyr_this <- unique(input_list$data$srv_sel_bicubic_selstyr[r, block_years, f])
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

          srv_sel_bicubic_Wyr[r, , 1:n_yr_nodes_this, b, f] <- Wyr_this

        } # end b loop
      } # end r loop
    } # end f loop
  } # end if has_bicubic_srv_sel

  input_list$data$srv_sel_bicubic_Wbin <- srv_sel_bicubic_Wbin
  input_list$data$srv_sel_bicubic_Wyr <- srv_sel_bicubic_Wyr

  max_srvsel_pars <- max(sel_pars_vec) # maximum number of selectivity parameters across all forms
  if("srv_fixed_sel_pars" %in% names(starting_values)) input_list$par$srv_fixed_sel_pars <- starting_values$srv_fixed_sel_pars
  else input_list$par$srv_fixed_sel_pars <- array(0, dim = c(input_list$data$n_regions, max_srvsel_pars, max_srvsel_blks, input_list$data$n_sexes, input_list$data$n_srv_fleets))

  # Survey catchability
  max_srvq_blks <- max(apply(input_list$data$srv_q_blocks, c(1,3), FUN = function(x) length(unique(x)))) # figure out maximum number of survey catchability blocks for a given reigon and fleet
  if("ln_srv_q" %in% names(starting_values)) input_list$par$ln_srv_q <- starting_values$ln_srv_q
  else input_list$par$ln_srv_q <- array(0, dim = c(input_list$data$n_regions, max_srvq_blks, input_list$data$n_srv_fleets))

  # Survey selectivity process error parameters
  if("srvsel_pe_pars" %in% names(starting_values)) input_list$par$srvsel_pe_pars <- starting_values$srvsel_pe_pars
  else input_list$par$srvsel_pe_pars <- array(0, dim = c(input_list$data$n_regions, max(max_srvsel_pars, 4), input_list$data$n_sexes, input_list$data$n_srv_fleets)) # dimensioned 4 as the max number of pars for process errors (e.g., sigmas), and then just map off if not using

  # Survey selectivity deviations
  if("ln_srvsel_devs" %in% names(starting_values)) input_list$par$ln_srvsel_devs <- starting_values$ln_srvsel_devs
  else input_list$par$ln_srvsel_devs <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years) + input_list$data$n_proj_yrs_devs, bins, input_list$data$n_sexes, input_list$data$n_srv_fleets))

  # Survey catchability covariate effects
  if("srv_q_coeff" %in% names(starting_values)) input_list$par$srv_q_coeff <- starting_values$srv_q_coeff
  else input_list$par$srv_q_coeff <- srv_q_coeff # input parameter array

  # Mapping Options ---------------------------------------------------------
  input_list$map$srv_q_coeff <- factor(map_srv_q_coeff) # set up mapping for catchability covariate
  input_list <- do_fixed_sel_pars_mapping(input_list, srv_fixed_sel_pars_spec, bins, srv_sel_nonpar_est_bins,
                                          prefix = "srv", fleet_field = "n_srv_fleets", use_field = "SrvIdx", fleet_label = "survey fleet")
  input_list <- do_q_mapping(input_list, srv_q_spec, prefix = "srv", fleet_field = "n_srv_fleets", fleet_label = "survey fleet")
  input_list <- do_sel_pe_pars_mapping(input_list, srvsel_pe_pars_spec, corr_opt_semipar, bins,
                                       prefix = "srv", fleet_field = "n_srv_fleets", use_field = "SrvIdx", fleet_label = "survey fleet")
  input_list <- do_sel_devs_mapping(input_list, srv_sel_devs_spec, srvsel_devs_shared_bins, bins,
                                    prefix = "srv", fleet_field = "n_srv_fleets", use_field = "SrvIdx", fleet_label = "survey fleet")


  # Print Messages ----------------------------------------------------------
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)
}
