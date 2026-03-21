#' Set up survey parameterisation for the operating model simulation
#'
#' Populates \code{sim_list} with all survey-related inputs needed by the
#' operating model: catchability, selectivity, survey timing, index type,
#' and age/length composition likelihood settings including overdispersion
#' and correlation parameters. Must be called after \code{\link{Setup_Sim_Dim}}.
#'
#' @param sim_list Simulation list returned by \code{\link{Setup_Sim_Dim}}.
#' @param srv_sel_input Survey selectivity array
#'   \code{[n_regions × n_yrs × n_ages × n_sexes × n_srv_fleets × n_sims]}.
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
#' @param pop_comp_srvage_like Integer or character vector \code{[n_srv_fleets]}
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
#' @param pop_SrvAgeComps_Type Array \code{[n_yrs × n_srv_fleets]} specifying
#'   population-specific age composition structure. Default: 2.
#' @param pop_comp_srvlen_like Integer or character vector \code{[n_srv_fleets]}
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
#' @param pop_SrvLenComps_Type Array \code{[n_yrs × n_srv_fleets]} specifying
#'   population-specific length composition structure. Default: 2.
#'
#' @return The input \code{sim_list} with survey-related fields appended:
#'   \code{$srv_sel}, \code{$srv_q}, \code{$ObsSrvIdx_SE}, \code{$ObsSrvIdx_pop_SE},
#'   \code{$t_srv}, \code{$srv_idx_type}, \code{$comp_srvage_like}, \code{$ISS_SrvAgeComps},
#'   \code{$ln_SrvAge_theta}, \code{$ln_SrvAge_theta_agg}, \code{$SrvAge_corr_pars_agg},
#'   \code{$SrvAge_corr_pars}, \code{$SrvAgeComps_Type}, \code{$comp_srvlen_like},
#'   \code{$ISS_SrvLenComps}, \code{$ln_SrvLen_theta}, \code{$ln_SrvLen_theta_agg},
#'   \code{$SrvLen_corr_pars_agg}, \code{$SrvLen_corr_pars}, \code{$SrvLenComps_Type},
#'   \code{$pop_comp_srvage_like}, \code{$ISS_SrvAgeComps_pop}, \code{$ln_SrvAge_pop_theta},
#'   \code{$ln_SrvAge_pop_theta_agg}, \code{$SrvAge_pop_corr_pars_agg}, \code{$SrvAge_pop_corr_pars},
#'   \code{$pop_SrvAgeComps_Type}, \code{$pop_comp_srvlen_like}, \code{$ISS_SrvLenComps_pop},
#'   \code{$ln_SrvLen_pop_theta}, \code{$ln_SrvLen_pop_theta_agg}, \code{$SrvLen_pop_corr_pars_agg},
#'   \code{$SrvLen_pop_corr_pars}, \code{$pop_SrvLenComps_Type}. Character-coded
#'   inputs are converted to integer equivalents before storage.
#'
#' @export Setup_Sim_Survey
#' @family Simulation Setup
Setup_Sim_Survey <- function(ObsSrvIdx_SE = array(0.2, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas,  sim_list$n_srv_fleets)),
                             ObsSrvIdx_pop_SE = array(0.2, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas,  sim_list$n_srv_fleets)),
                             sim_list,
                             srv_sel_input,
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
                             pop_comp_srvage_like = rep(0, sim_list$n_srv_fleets),
                             ISS_SrvAgeComps_pop = array(100, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims)),
                             ln_SrvAge_pop_theta = array(log(1), dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_sexes, sim_list$n_srv_fleets)),
                             ln_SrvAge_pop_theta_agg = array(log(1), dim = c(sim_list$n_pop, sim_list$n_srv_fleets)),
                             SrvAge_pop_corr_pars = array(0.01, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_sexes, sim_list$n_srv_fleets, 2)),
                             SrvAge_pop_corr_pars_agg = array(0.01, dim = c(sim_list$n_pop, sim_list$n_srv_fleets)),
                             pop_SrvAgeComps_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_srv_fleets)),
                             pop_comp_srvlen_like = rep(0, sim_list$n_srv_fleets),
                             ISS_SrvLenComps_pop = array(100, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims)),
                             ln_SrvLen_pop_theta = array(log(1), dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_sexes, sim_list$n_srv_fleets)),
                             ln_SrvLen_pop_theta_agg = array(log(1), dim = c(sim_list$n_pop, sim_list$n_srv_fleets)),
                             SrvLen_pop_corr_pars = array(0.01, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_sexes, sim_list$n_srv_fleets, 2)),
                             SrvLen_pop_corr_pars_agg = array(0.01, dim = c(sim_list$n_pop, sim_list$n_srv_fleets)),
                             pop_SrvLenComps_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_srv_fleets))
                             ) {

  # Convert character inputs to numeric codes
  srv_idx_type <- convert_to_numeric(srv_idx_type, list(abd = 0, biom = 1))
  comp_srvage_like <- convert_to_numeric(comp_srvage_like, list(Multinomial = 0,  `Dirichlet-Multinomial` = 1, `iid-Logistic-Normal` = 2, `1d-Logistic-Normal` = 3, `2d-Logistic-Normal` = 4))
  comp_srvlen_like <- convert_to_numeric(comp_srvlen_like, list(Multinomial = 0, `Dirichlet-Multinomial` = 1, `iid-Logistic-Normal` = 2, `1d-Logistic-Normal` = 3, `2d-Logistic-Normal` = 4))
  SrvAgeComps_Type <- convert_to_numeric(SrvAgeComps_Type,  list(agg = 0, spltRspltS = 1, spltRjntS = 2, none = 999))
  SrvLenComps_Type <- convert_to_numeric(SrvLenComps_Type,  list(agg = 0, spltRspltS = 1, spltRjntS = 2, none = 999))
  pop_SrvAgeComps_Type <- convert_to_numeric(pop_SrvAgeComps_Type,  list(agg = 0, spltRspltS = 1, spltRjntS = 2, none = 999))
  pop_SrvLenComps_Type <- convert_to_numeric(pop_SrvLenComps_Type,  list(agg = 0, spltRspltS = 1, spltRjntS = 2, none = 999))

  # Validate dimensions of all input parameters
  check_sim_dimensions(srv_sel_input, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs,
                       n_ages = sim_list$n_ages, n_sexes = sim_list$n_sexes,
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
  check_sim_dimensions(pop_comp_srvage_like, n_srv_fleets = sim_list$n_srv_fleets, what = "pop_comp_srvage_like")
  check_sim_dimensions(ISS_SrvAgeComps_pop, n_pop = sim_list$n_pop, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs, n_seas = sim_list$n_seas, n_sexes = sim_list$n_sexes, n_srv_fleets = sim_list$n_srv_fleets, n_sims = sim_list$n_sims, what = "ISS_SrvAgeComps_pop")
  check_sim_dimensions(ln_SrvAge_pop_theta, n_pop = sim_list$n_pop, n_regions = sim_list$n_regions, n_sexes = sim_list$n_sexes, n_srv_fleets = sim_list$n_srv_fleets, what = "ln_SrvAge_pop_theta")
  check_sim_dimensions(ln_SrvAge_pop_theta_agg, n_pop = sim_list$n_pop, n_srv_fleets = sim_list$n_srv_fleets, what = "ln_SrvAge_pop_theta_agg")
  check_sim_dimensions(SrvAge_pop_corr_pars_agg, n_pop = sim_list$n_pop, n_srv_fleets = sim_list$n_srv_fleets, what = "SrvAge_pop_corr_pars_agg")
  check_sim_dimensions(SrvAge_pop_corr_pars, n_pop = sim_list$n_pop, n_regions = sim_list$n_regions, n_sexes = sim_list$n_sexes, n_srv_fleets = sim_list$n_srv_fleets, what = "SrvAge_pop_corr_pars")
  check_sim_dimensions(pop_SrvAgeComps_Type, n_years = sim_list$n_yrs, n_srv_fleets = sim_list$n_srv_fleets, what = "pop_SrvAgeComps_Type")


  # Validate suvey length composition parameters
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
  check_sim_dimensions(pop_comp_srvlen_like, n_srv_fleets = sim_list$n_srv_fleets, what = "pop_comp_srvlen_like")
  check_sim_dimensions(ISS_SrvLenComps_pop, n_pop = sim_list$n_pop, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs, n_seas = sim_list$n_seas, n_sexes = sim_list$n_sexes, n_srv_fleets = sim_list$n_srv_fleets, n_sims = sim_list$n_sims, what = "ISS_SrvLenComps_pop")
  check_sim_dimensions(ln_SrvLen_pop_theta, n_pop = sim_list$n_pop, n_regions = sim_list$n_regions, n_sexes = sim_list$n_sexes, n_srv_fleets = sim_list$n_srv_fleets, what = "ln_SrvLen_pop_theta")
  check_sim_dimensions(ln_SrvLen_pop_theta_agg, n_pop = sim_list$n_pop, n_srv_fleets = sim_list$n_srv_fleets, what = "ln_SrvLen_pop_theta_agg")
  check_sim_dimensions(SrvLen_pop_corr_pars_agg, n_pop = sim_list$n_pop, n_srv_fleets = sim_list$n_srv_fleets, what = "SrvLen_pop_corr_pars_agg")
  check_sim_dimensions(SrvLen_pop_corr_pars, n_pop = sim_list$n_pop, n_regions = sim_list$n_regions, n_sexes = sim_list$n_sexes, n_srv_fleets = sim_list$n_srv_fleets, what = "SrvLen_pop_corr_pars")
  check_sim_dimensions(pop_SrvLenComps_Type, n_years = sim_list$n_yrs, n_srv_fleets = sim_list$n_srv_fleets, what = "pop_SrvLenComps_Type")

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
  sim_list$pop_comp_srvage_like <- pop_comp_srvage_like
  sim_list$ISS_SrvAgeComps_pop <- ISS_SrvAgeComps_pop
  sim_list$ln_SrvAge_pop_theta <- ln_SrvAge_pop_theta
  sim_list$ln_SrvAge_pop_theta_agg <- ln_SrvAge_pop_theta_agg
  sim_list$SrvAge_pop_corr_pars <- SrvAge_pop_corr_pars
  sim_list$SrvAge_pop_corr_pars_agg <- SrvAge_pop_corr_pars_agg
  sim_list$pop_SrvAgeComps_Type <- pop_SrvAgeComps_Type

  sim_list$pop_comp_srvlen_like <- pop_comp_srvlen_like
  sim_list$ISS_SrvLenComps_pop <- ISS_SrvLenComps_pop
  sim_list$ln_SrvLen_pop_theta <- ln_SrvLen_pop_theta
  sim_list$ln_SrvLen_pop_theta_agg <- ln_SrvLen_pop_theta_agg
  sim_list$SrvLen_pop_corr_pars <- SrvLen_pop_corr_pars
  sim_list$SrvLen_pop_corr_pars_agg <- SrvLen_pop_corr_pars_agg
  sim_list$pop_SrvLenComps_Type <- pop_SrvLenComps_Type

  return(sim_list)

} # end function

#' Map survey age composition overdispersion parameters
#'
#' Internal helper called by \code{\link{Setup_Mod_SrvIdx_and_Comps}} to
#' construct the TMB/RTMB factor maps for \code{ln_SrvAge_theta}
#' \code{[n_regions × n_sexes × n_srv_fleets]} and
#' \code{ln_SrvAge_theta_agg} \code{[n_srv_fleets]}, the log-scale
#' overdispersion parameters for survey age composition likelihoods.
#'
#' Parameters are activated only when the fleet's likelihood type is not
#' multinomial (\code{SrvAgeComps_LikeType != 0}) and at least one year of
#' age composition data is used (\code{UseSrvAgeComps > 0}). Within active
#' fleets, the composition type (\code{SrvAgeComps_Type}) determines which
#' parameter array is populated: aggregated compositions use
#' \code{ln_SrvAge_theta_agg}; split or joint-by-sex compositions use
#' \code{ln_SrvAge_theta}. Fleets using multinomial likelihoods or with no
#' active age composition data have all parameters mapped to \code{NA}.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists. Requires \code{$data$n_srv_fleets}, \code{$data$n_regions},
#'   \code{$data$n_sexes}, \code{$data$SrvAgeComps_Type},
#'   \code{$data$SrvAgeComps_LikeType}, and \code{$data$UseSrvAgeComps}.
#'
#' @return The input \code{input_list} with \code{$map$ln_SrvAge_theta} and
#'   \code{$map$ln_SrvAge_theta_agg} set to factor vectors. Active parameters
#'   receive sequential integer indices; inactive parameters are \code{NA}.
#'
#' @seealso \code{\link{do_SrvLen_theta_mapping}} for the analogous length
#'   composition overdispersion mapping;
#'   \code{\link{do_SrvAge_corr_pars_mapping}} for the associated correlation
#'   parameter mapping.
#'
#' @keywords internal
do_SrvAge_theta_mapping <- function(input_list) {

  # setup counters
  counter_srvage_agg <- 1
  counter_srvage <- 1

  # initialize array to set up mapping
  map_SrvAge_theta <- input_list$par$ln_SrvAge_theta
  map_SrvAge_theta_agg <- input_list$par$ln_SrvAge_theta_agg
  map_SrvAge_theta[] <- NA
  map_SrvAge_theta_agg[] <- NA

  for(f in 1:input_list$data$n_srv_fleets) {

    # get unique survey comp types
    srvage_comp_type <- unique(input_list$data$SrvAgeComps_Type[,f])

    # If aggregated (ages)
    if(any(srvage_comp_type == 0) && input_list$data$SrvAgeComps_LikeType[f] != 0) {
      map_SrvAge_theta_agg[f] <- counter_srvage_agg
      counter_srvage_agg <- counter_srvage_agg + 1 # aggregated
    }

    # Loop through to make sure mapping stuff off correctly
    for(r in 1:input_list$data$n_regions) {
      for(s in 1:input_list$data$n_sexes) {

        # if split by sex and region
        if(any(srvage_comp_type == 1) && input_list$data$SrvAgeComps_LikeType[f] != 0) {
          map_SrvAge_theta[r,s,f] <- counter_srvage
          counter_srvage <- counter_srvage + 1 # split by sex and region
        }

        # joint by sex, split by region
        if(any(srvage_comp_type == 2) && input_list$data$SrvAgeComps_LikeType[f] != 0 && s == 1) {
          map_SrvAge_theta[r,1,f] <- counter_srvage
          counter_srvage <- counter_srvage + 1 # joint by sex, split by region
        }

      } # end s loop
    } # end r loop

    # If we are using a multinomial or there aren't any age comps for a given fleet
    if(input_list$data$SrvAgeComps_LikeType[f] == 0 || sum(input_list$data$UseSrvAgeComps[,,,f]) == 0) {
      map_SrvAge_theta[,,f] <- NA
      map_SrvAge_theta_agg[f] <- NA
    }

  } # end f loop

  # Input into mapping list
  input_list$map$ln_SrvAge_theta <- factor(map_SrvAge_theta)
  input_list$map$ln_SrvAge_theta_agg <- factor(map_SrvAge_theta_agg)

  return(input_list)
}

#' Map population survey age composition overdispersion parameters
#'
#' Constructs factor maps for \code{ln_SrvAge_pop_theta} (fleet- region- and
#' sex-specific overdispersion) and \code{ln_SrvAge_pop_theta_agg} (aggregated
#' overdispersion) based on the composition type and likelihood specified in
#' \code{$data$pop_SrvAgeComps_Type} and \code{$data$pop_SrvAgeComps_LikeType}.
#' Parameters are mapped to \code{NA} for fleets using multinomial likelihoods
#' (\code{LikeType == 0}) or with no observed age compositions.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists. Requires \code{pop_SrvAgeComps_Type}, \code{pop_SrvAgeComps_LikeType},
#'   and \code{UseSrvAgeComps_pop} to be set by
#'   \code{\link{Setup_Mod_SrvIdx_and_Comps}}.
#'
#' @return The input \code{input_list} with \code{$map$ln_SrvAge_pop_theta} and
#'   \code{$map$ln_SrvAge_pop_theta_agg} set to factor vectors. Active parameters
#'   receive sequential integer indices; inactive parameters are \code{NA}.
#'
#' @keywords internal
do_SrvAge_pop_theta_mapping <- function(input_list) {

  # setup counters
  counter_srvage_agg <- 1
  counter_srvage <- 1

  # initialize array to set up mapping
  map_SrvAge_pop_theta <- input_list$par$ln_SrvAge_pop_theta
  map_SrvAge_pop_theta_agg <- input_list$par$ln_SrvAge_pop_theta_agg
  map_SrvAge_pop_theta[] <- NA
  map_SrvAge_pop_theta_agg[] <- NA

  for(p in 1:input_list$data$n_pop) {
    for(f in 1:input_list$data$n_srv_fleets) {

      # get unique survey comp types
      srvage_comp_type <- unique(input_list$data$pop_SrvAgeComps_Type[,f])

      # If aggregated (ages)
      if(any(srvage_comp_type == 0) && input_list$data$pop_SrvAgeComps_LikeType[f] != 0) {
        map_SrvAge_pop_theta_agg[p,f] <- counter_srvage_agg
        counter_srvage_agg <- counter_srvage_agg + 1 # aggregated
      }

      # Loop through to make sure mapping stuff off correctly
      for(r in 1:input_list$data$n_regions) {
        for(s in 1:input_list$data$n_sexes) {

          # if split by sex and region
          if(any(srvage_comp_type == 1) && input_list$data$pop_SrvAgeComps_LikeType[f] != 0) {
            map_SrvAge_pop_theta[p,r,s,f] <- counter_srvage
            counter_srvage <- counter_srvage + 1 # split by sex and region
          }

          # joint by sex, split by region
          if(any(srvage_comp_type == 2) && input_list$data$pop_SrvAgeComps_LikeType[f] != 0 && s == 1) {
            map_SrvAge_pop_theta[p,r,1,f] <- counter_srvage
            counter_srvage <- counter_srvage + 1 # joint by sex, split by region
          }

        } # end s loop
      } # end r loop

      # If we are using a multinomial or there aren't any age comps for a given fleet
      if(input_list$data$pop_SrvAgeComps_LikeType[f] == 0 || sum(input_list$data$UseSrvAgeComps_pop[p,,,,f]) == 0) {
        map_SrvAge_pop_theta[p,,,f] <- NA
        map_SrvAge_pop_theta_agg[p,f] <- NA
      }

    } # end f loop
  } # end p loop

  # Input into mapping list
  input_list$map$ln_SrvAge_pop_theta <- factor(map_SrvAge_pop_theta)
  input_list$map$ln_SrvAge_pop_theta_agg <- factor(map_SrvAge_pop_theta_agg)

  return(input_list)
}


#' Map survey length composition overdispersion parameters
#'
#' Internal helper called by \code{\link{Setup_Mod_SrvIdx_and_Comps}} to
#' construct the TMB/RTMB factor maps for \code{ln_SrvLen_theta}
#' \code{[n_regions × n_sexes × n_srv_fleets]} and
#' \code{ln_SrvLen_theta_agg} \code{[n_srv_fleets]}, the log-scale
#' overdispersion parameters for survey length composition likelihoods.
#' Follows identical activation logic to \code{\link{do_SrvAge_theta_mapping}}
#' but operates on \code{SrvLenComps_LikeType}, \code{SrvLenComps_Type}, and
#' \code{UseSrvLenComps}.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists. Requires \code{$data$n_srv_fleets}, \code{$data$n_regions},
#'   \code{$data$n_sexes}, \code{$data$SrvLenComps_Type},
#'   \code{$data$SrvLenComps_LikeType}, and \code{$data$UseSrvLenComps}.
#'
#' @return The input \code{input_list} with \code{$map$ln_SrvLen_theta} and
#'   \code{$map$ln_SrvLen_theta_agg} set to factor vectors. Active parameters
#'   receive sequential integer indices; inactive parameters are \code{NA}.
#'
#' @seealso \code{\link{do_SrvAge_theta_mapping}} for the analogous age
#'   composition overdispersion mapping;
#'   \code{\link{do_SrvLen_corr_pars_mapping}} for the associated correlation
#'   parameter mapping.
#'
#' @keywords internal
do_SrvLen_theta_mapping <- function(input_list) {

  # setup counters
  counter_srvlen_agg <- 1
  counter_srvlen <- 1

  # initialize array to set up mapping
  map_SrvLen_theta <- input_list$par$ln_SrvLen_theta
  map_SrvLen_theta_agg <- input_list$par$ln_SrvLen_theta_agg
  map_SrvLen_theta[] <- NA
  map_SrvLen_theta_agg[] <- NA

  for(f in 1:input_list$data$n_srv_fleets) {

    # get unique survey comp types
    srvlen_comp_type <- unique(input_list$data$SrvLenComps_Type[,f])

    # If aggregated (ages)
    if(any(srvlen_comp_type == 0) && input_list$data$SrvLenComps_LikeType[f] != 0) {
      map_SrvLen_theta_agg[f] <- counter_srvlen_agg
      counter_srvlen_agg <- counter_srvlen_agg + 1 # aggregated
    }

    # Loop through to make sure mapping stuff off correctly
    for(r in 1:input_list$data$n_regions) {
      for(s in 1:input_list$data$n_sexes) {

        # if split by sex and region
        if(any(srvlen_comp_type == 1) && input_list$data$SrvLenComps_LikeType[f] != 0) {
          map_SrvLen_theta[r,s,f] <- counter_srvlen
          counter_srvlen <- counter_srvlen + 1 # split by sex and region
        }

        # joint by sex, split by region
        if(any(srvlen_comp_type == 2) && input_list$data$SrvLenComps_LikeType[f] != 0 && s == 1) {
          map_SrvLen_theta[r,1,f] <- counter_srvlen
          counter_srvlen <- counter_srvlen + 1 # joint by sex, split by region
        }

      } # end s loop
    } # end r loop

    # If we are using a multinomial or there aren't any lenght comps for a given fleet
    if(input_list$data$SrvLenComps_LikeType[f] == 0 || sum(input_list$data$UseSrvLenComps[,,,f]) == 0) {
      map_SrvLen_theta[,,f] <- NA
      map_SrvLen_theta_agg[f] <- NA
    }

  } # end f loop

  # Input into mapping list
  input_list$map$ln_SrvLen_theta <- factor(map_SrvLen_theta)
  input_list$map$ln_SrvLen_theta_agg <- factor(map_SrvLen_theta_agg)

  return(input_list)
}

#' Map population survey length composition overdispersion parameters
#'
#' Constructs factor maps for \code{ln_SrvLen_pop_theta} (fleet- region- and
#' sex-specific overdispersion) and \code{ln_SrvLen_pop_theta_agg} (aggregated
#' overdispersion) based on the composition type and likelihood specified in
#' \code{$data$pop_SrvLenComps_Type} and \code{$data$pop_SrvLenComps_LikeType}.
#' Parameters are mapped to \code{NA} for fleets using multinomial likelihoods
#' (\code{LikeType == 0}) or with no observed len compositions.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists. Requires \code{pop_SrvLenComps_Type}, \code{pop_SrvLenComps_LikeType},
#'   and \code{UseSrvLenComps_pop} to be set by
#'   \code{\link{Setup_Mod_SrvIdx_and_Comps}}.
#'
#' @return The input \code{input_list} with \code{$map$ln_SrvLen_pop_theta} and
#'   \code{$map$ln_SrvLen_pop_theta_agg} set to factor vectors. Active parameters
#'   receive sequential integer indices; inactive parameters are \code{NA}.
#'
#' @keywords internal
do_SrvLen_pop_theta_mapping <- function(input_list) {

  # setup counters
  counter_srvlen_agg <- 1
  counter_srvlen <- 1

  # initialize array to set up mapping
  map_SrvLen_pop_theta <- input_list$par$ln_SrvLen_pop_theta
  map_SrvLen_pop_theta_agg <- input_list$par$ln_SrvLen_pop_theta_agg
  map_SrvLen_pop_theta[] <- NA
  map_SrvLen_pop_theta_agg[] <- NA

  for(p in 1:input_list$data$n_pop) {
    for(f in 1:input_list$data$n_srv_fleets) {

      # get unique survey comp types
      srvlen_comp_type <- unique(input_list$data$pop_SrvLenComps_Type[,f])

      # If aggregated (lens)
      if(any(srvlen_comp_type == 0) && input_list$data$pop_SrvLenComps_LikeType[f] != 0) {
        map_SrvLen_pop_theta_agg[p,f] <- counter_srvlen_agg
        counter_srvlen_agg <- counter_srvlen_agg + 1 # aggregated
      }

      # Loop through to make sure mapping stuff off correctly
      for(r in 1:input_list$data$n_regions) {
        for(s in 1:input_list$data$n_sexes) {

          # if split by sex and region
          if(any(srvlen_comp_type == 1) && input_list$data$pop_SrvLenComps_LikeType[f] != 0) {
            map_SrvLen_pop_theta[p,r,s,f] <- counter_srvlen
            counter_srvlen <- counter_srvlen + 1 # split by sex and region
          }

          # joint by sex, split by region
          if(any(srvlen_comp_type == 2) && input_list$data$pop_SrvLenComps_LikeType[f] != 0 && s == 1) {
            map_SrvLen_pop_theta[p,r,1,f] <- counter_srvlen
            counter_srvlen <- counter_srvlen + 1 # joint by sex, split by region
          }

        } # end s loop
      } # end r loop

      # If we are using a multinomial or there aren't any len comps for a given fleet
      if(input_list$data$pop_SrvLenComps_LikeType[f] == 0 || sum(input_list$data$UseSrvLenComps_pop[p,,,,f]) == 0) {
        map_SrvLen_pop_theta[p,,,f] <- NA
        map_SrvLen_pop_theta_agg[p,f] <- NA
      }

    } # end f loop
  } # end p loop

  # Input into mapping list
  input_list$map$ln_SrvLen_pop_theta <- factor(map_SrvLen_pop_theta)
  input_list$map$ln_SrvLen_pop_theta_agg <- factor(map_SrvLen_pop_theta_agg)

  return(input_list)
}

#' Map survey age composition correlation parameters
#'
#' Internal helper called by \code{\link{Setup_Mod_SrvIdx_and_Comps}} to
#' construct the TMB/RTMB factor maps for \code{SrvAge_corr_pars}
#' \code{[n_regions × n_sexes × n_srv_fleets × 2]} and
#' \code{SrvAge_corr_pars_agg} \code{[n_srv_fleets]}, the correlation
#' parameters for logistic-normal survey age composition likelihoods.
#'
#' The trailing dimension of \code{SrvAge_corr_pars} distinguishes the age
#' AR1 correlation (index 1, used by likelihoods 3 and 4) from the sex
#' correlation (index 2, used only by likelihood 4). Parameters are activated
#' conditionally on the fleet's likelihood type and composition type:
#' aggregated compositions (\code{type = 0}) activate
#' \code{SrvAge_corr_pars_agg} only for likelihood 3; split or joint-by-sex
#' compositions activate \code{SrvAge_corr_pars} with index 2 skipped when
#' \code{n_sexes = 1}. Fleets using multinomial likelihoods or with no active
#' age composition data have all parameters mapped to \code{NA}.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists. Requires \code{$data$n_srv_fleets}, \code{$data$n_regions},
#'   \code{$data$n_sexes}, \code{$data$SrvAgeComps_Type},
#'   \code{$data$SrvAgeComps_LikeType}, and \code{$data$UseSrvAgeComps}.
#'
#' @return The input \code{input_list} with \code{$map$SrvAge_corr_pars} and
#'   \code{$map$SrvAge_corr_pars_agg} set to factor vectors. Active parameters
#'   receive sequential integer indices; inactive parameters are \code{NA}.
#'
#' @seealso \code{\link{do_SrvAge_theta_mapping}} for the overdispersion
#'   parameter mapping; \code{\link{do_SrvLen_corr_pars_mapping}} for the
#'   analogous length composition correlation mapping.
#'
#' @keywords internal
do_SrvAge_corr_pars_mapping <- function(input_list) {

  # setup counters
  counter_srvage_corr <- 1
  counter_srvage_corr_agg <- 1

  # initialize array to set up mapping
  map_SrvAge_corr_pars <- input_list$par$SrvAge_corr_pars
  map_SrvAge_corr_pars_agg <- input_list$par$SrvAge_corr_pars_agg
  map_SrvAge_corr_pars[] <- NA
  map_SrvAge_corr_pars_agg[] <- NA

  for(f in 1:input_list$data$n_srv_fleets) {

    # No overdispersion parameters estimated
    if(input_list$data$SrvAgeComps_LikeType[f] == 0 || sum(input_list$data$UseSrvAgeComps[,,,f]) == 0) {
      map_SrvAge_corr_pars[,,f,] <- NA
      map_SrvAge_corr_pars_agg[f] <- NA
      next # skip if none
    }

    # get unique survey comp types
    srvage_comp_type <- unique(input_list$data$SrvAgeComps_Type[,f])

    # Aggregated Correlation Parameters
    if(any(srvage_comp_type == 0) && input_list$data$SrvAgeComps_LikeType[f] != 0) {
      if(input_list$data$SrvAgeComps_LikeType[f] == 3) {
        map_SrvAge_corr_pars_agg[f] <- counter_srvage_corr_agg
        counter_srvage_corr_agg <- counter_srvage_corr_agg + 1 # aggregated
      }
    }

    # Loop through to make sure mapping stuff off correctly
    for(r in 1:input_list$data$n_regions) {
      for(s in 1:input_list$data$n_sexes) {

        # Split by region and sex
        if(any(srvage_comp_type == 1) && input_list$data$SrvAgeComps_LikeType[f] != 0) {
          if(input_list$data$SrvAgeComps_LikeType[f] == 3) {
            map_SrvAge_corr_pars[r,s,f,1] <- counter_srvage_corr
            counter_srvage_corr <- counter_srvage_corr + 1
          }
        }

        # Joint by sex, split by region
        if(any(srvage_comp_type == 2) && input_list$data$SrvAgeComps_LikeType[f] != 0 && s == 1) {

          # 1dar1 correlation
          if(input_list$data$SrvAgeComps_LikeType[f] == 3) {
            map_SrvAge_corr_pars[r,1,f,1] <- counter_srvage_corr
            counter_srvage_corr <- counter_srvage_corr + 1
          }

          # 2dar1 correlation
          if(input_list$data$SrvAgeComps_LikeType[f] == 4) {
            for(i in 1:2) {
              if(i == 2 && input_list$data$n_sexes == 1) next # skip if we only have 1 sex
              map_SrvAge_corr_pars[r,1,f,i] <- counter_srvage_corr
              counter_srvage_corr <- counter_srvage_corr + 1
            } # end i
          } # end if

        }
      } # end s loop
    } # end r loop

  } # end f loop

  # Input into mapping list
  input_list$map$SrvAge_corr_pars_agg <- factor(map_SrvAge_corr_pars_agg)
  input_list$map$SrvAge_corr_pars <- factor(map_SrvAge_corr_pars)

  return(input_list)
}

#' Map population survey age composition correlation parameters
#'
#' Constructs factor maps for \code{SrvAge_pop_corr_pars} (region- and sex-specific
#' AR1 and sex correlation parameters) and \code{SrvAge_pop_corr_pars_agg}
#' (aggregated correlation parameters) for 1D and 2D logistic-normal age
#' composition likelihoods.
#'
#' Parameters are activated only when
#' \code{pop_SrvAgeComps_LikeType} is in \code{c(3, 4)}.
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
#'   \code{map\$SrvAge_pop_corr_pars} and
#'   \code{map\$SrvAge_pop_corr_pars_agg} set to factor vectors.
#'
#' @keywords internal
do_SrvAge_pop_corr_pars_mapping <- function(input_list) {

  # setup counters
  counter_srvage_corr <- 1
  counter_srvage_corr_agg <- 1

  # initialize array to set up mapping
  map_SrvAge_pop_corr_pars <- input_list$par$SrvAge_pop_corr_pars
  map_SrvAge_pop_corr_pars_agg <- input_list$par$SrvAge_pop_corr_pars_agg
  map_SrvAge_pop_corr_pars[] <- NA
  map_SrvAge_pop_corr_pars_agg[] <- NA

  for(p in 1:input_list$data$n_pop) {
    for(f in 1:input_list$data$n_srv_fleets) {

      # No overdispersion parameters estimated
      if(input_list$data$pop_SrvAgeComps_LikeType[f] == 0 || sum(input_list$data$UseSrvAgeComps_pop[p,,,,f]) == 0) {
        map_SrvAge_pop_corr_pars[p,,,f,] <- NA
        map_SrvAge_pop_corr_pars_agg[p,f] <- NA
        next # skip if none
      }

      # get unique survey comp types
      srvage_comp_type <- unique(input_list$data$SrvAgeComps_Type[,f])

      # Aggregated Correlation Parameters
      if(any(srvage_comp_type == 0) && input_list$data$pop_SrvAgeComps_LikeType[f] != 0) {
        if(input_list$data$pop_SrvAgeComps_LikeType[f] == 3) {
          map_SrvAge_pop_corr_pars_agg[p,f] <- counter_srvage_corr_agg
          counter_srvage_corr_agg <- counter_srvage_corr_agg + 1 # aggregated
        }
      }

      # Loop through to make sure mapping stuff off correctly
      for(r in 1:input_list$data$n_regions) {
        for(s in 1:input_list$data$n_sexes) {

          # Split by region and sex
          if(any(srvage_comp_type == 1) && input_list$data$pop_SrvAgeComps_LikeType[f] != 0) {
            if(input_list$data$pop_SrvAgeComps_LikeType[f] == 3) {
              map_SrvAge_pop_corr_pars[p,r,s,f,1] <- counter_srvage_corr
              counter_srvage_corr <- counter_srvage_corr + 1
            }
          }

          # Joint by sex, split by region
          if(any(srvage_comp_type == 2) && input_list$data$pop_SrvAgeComps_LikeType[f] != 0 && s == 1) {

            # 1dar1 correlation
            if(input_list$data$pop_SrvAgeComps_LikeType[f] == 3) {
              map_SrvAge_pop_corr_pars[p,r,1,f,1] <- counter_srvage_corr
              counter_srvage_corr <- counter_srvage_corr + 1
            }

            # 2dar1 correlation
            if(input_list$data$pop_SrvAgeComps_LikeType[f] == 4) {
              for(i in 1:2) {
                if(i == 2 && input_list$data$n_sexes == 1) next # skip if we only have 1 sex
                map_SrvAge_pop_corr_pars[p,r,1,f,i] <- counter_srvage_corr
                counter_srvage_corr <- counter_srvage_corr + 1
              } # end i
            } # end if

          }
        } # end s loop
      } # end r loop

    } # end f loop
  } # end p loop

  # Input into mapping list
  input_list$map$SrvAge_pop_corr_pars_agg <- factor(map_SrvAge_pop_corr_pars_agg)
  input_list$map$SrvAge_pop_corr_pars <- factor(map_SrvAge_pop_corr_pars)

  return(input_list)
}

#' Map survey length composition correlation parameters
#'
#' Internal helper called by \code{\link{Setup_Mod_SrvIdx_and_Comps}} to
#' construct the TMB/RTMB factor maps for \code{SrvLen_corr_pars}
#' \code{[n_regions × n_sexes × n_srv_fleets × 2]} and
#' \code{SrvLen_corr_pars_agg} \code{[n_srv_fleets]}, the correlation
#' parameters for logistic-normal survey length composition likelihoods.
#' Follows identical activation logic to
#' \code{\link{do_SrvAge_corr_pars_mapping}} but operates on
#' \code{SrvLenComps_LikeType}, \code{SrvLenComps_Type}, and
#' \code{UseSrvLenComps}.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists. Requires \code{$data$n_srv_fleets}, \code{$data$n_regions},
#'   \code{$data$n_sexes}, \code{$data$SrvLenComps_Type},
#'   \code{$data$SrvLenComps_LikeType}, and \code{$data$UseSrvLenComps}.
#'
#' @return The input \code{input_list} with \code{$map$SrvLen_corr_pars} and
#'   \code{$map$SrvLen_corr_pars_agg} set to factor vectors. Active parameters
#'   receive sequential integer indices; inactive parameters are \code{NA}.
#'
#' @seealso \code{\link{do_SrvAge_corr_pars_mapping}} for the analogous age
#'   composition correlation mapping; \code{\link{do_SrvLen_theta_mapping}}
#'   for the overdispersion parameter mapping.
#'
#' @keywords internal
do_SrvLen_corr_pars_mapping <- function(input_list) {

  # setup counters
  counter_srvlen_corr <- 1
  counter_srvlen_corr_agg <- 1

  # initialize array to set up mapping
  map_SrvLen_corr_pars <- input_list$par$SrvLen_corr_pars
  map_SrvLen_corr_pars_agg <- input_list$par$SrvLen_corr_pars_agg
  map_SrvLen_corr_pars[] <- NA
  map_SrvLen_corr_pars_agg[] <- NA

  for(f in 1:input_list$data$n_srv_fleets) {

    # No overdispersion parameters estimated
    if(input_list$data$SrvLenComps_LikeType[f] == 0 || sum(input_list$data$UseSrvLenComps[,,,f]) == 0) {
      map_SrvLen_corr_pars[,,f,] <- NA
      map_SrvLen_corr_pars_agg[f] <- NA
      next # skip if none should be estimated
    }

    # get unique survey comp types
    srvlen_comp_type <- unique(input_list$data$SrvLenComps_Type[,f])

    # Aggregated Correlation Parameters
    if(any(srvlen_comp_type == 0) && input_list$data$SrvLenComps_LikeType[f] != 0) {
      if(input_list$data$SrvLenComps_LikeType[f] == 3) {
        map_SrvLen_corr_pars_agg[f] <- counter_srvlen_corr_agg
        counter_srvlen_corr_agg <- counter_srvlen_corr_agg + 1 # aggregated
      }
    }

    # Loop through to make sure mapping stuff off correctly
    for(r in 1:input_list$data$n_regions) {
      for(s in 1:input_list$data$n_sexes) {

        # Split by region and sex
        if(any(srvlen_comp_type == 1) && input_list$data$SrvLenComps_LikeType[f] != 0) {
          if(input_list$data$SrvLenComps_LikeType[f] == 3) {
            map_SrvLen_corr_pars[r,s,f,1] <- counter_srvlen_corr
            counter_srvlen_corr <- counter_srvlen_corr + 1
          }
        }

        # Joint by sex, split by region
        if(any(srvlen_comp_type == 2) && input_list$data$SrvLenComps_LikeType[f] != 0 && s == 1) {

          # 1dar1 correlation
          if(input_list$data$SrvLenComps_LikeType[f] == 3) {
            map_SrvLen_corr_pars[r,1,f,1] <- counter_srvlen_corr
            counter_srvlen_corr <- counter_srvlen_corr + 1
          }

          # 2dar1 correlation
          if(input_list$data$SrvLenComps_LikeType[f] == 4) {
            for(i in 1:2) {
              if(i == 2 && input_list$data$n_sexes == 1) next # skip if we only have 1 sex
              map_SrvLen_corr_pars[r,1,f,i] <- counter_srvlen_corr
              counter_srvlen_corr <- counter_srvlen_corr + 1
            } # end i
          } # end if

        }
      } # end s loop
    } # end r loop

  } # end f loop

  # Input into mapping list
  input_list$map$SrvLen_corr_pars_agg <- factor(map_SrvLen_corr_pars_agg)
  input_list$map$SrvLen_corr_pars <- factor(map_SrvLen_corr_pars)

  return(input_list)
}

#' Map population survey length composition correlation parameters
#'
#' Constructs factor maps for \code{SrvLen_pop_corr_pars} (region- and sex-specific
#' AR1 and sex correlation parameters) and \code{SrvLen_pop_corr_pars_agg}
#' (aggregated correlation parameters) for 1D and 2D logistic-normal len
#' composition likelihoods.
#'
#' Parameters are activated only when
#' \code{pop_SrvLenComps_LikeType} is in \code{c(3, 4)}.
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
#'   \code{map\$SrvLen_pop_corr_pars} and
#'   \code{map\$SrvLen_pop_corr_pars_agg} set to factor vectors.
#'
#' @keywords internal
do_SrvLen_pop_corr_pars_mapping <- function(input_list) {

  # setup counters
  counter_srvlen_corr <- 1
  counter_srvlen_corr_agg <- 1

  # initialize array to set up mapping
  map_SrvLen_pop_corr_pars <- input_list$par$SrvLen_pop_corr_pars
  map_SrvLen_pop_corr_pars_agg <- input_list$par$SrvLen_pop_corr_pars_agg
  map_SrvLen_pop_corr_pars[] <- NA
  map_SrvLen_pop_corr_pars_agg[] <- NA

  for(p in 1:input_list$data$n_pop) {
    for(f in 1:input_list$data$n_srv_fleets) {

      # No overdispersion parameters estimated
      if(input_list$data$pop_SrvLenComps_LikeType[f] == 0 || sum(input_list$data$UseSrvLenComps_pop[p,,,,f]) == 0) {
        map_SrvLen_pop_corr_pars[p,,,f,] <- NA
        map_SrvLen_pop_corr_pars_agg[p,f] <- NA
        next # skip if none
      }

      # get unique survey comp types
      srvlen_comp_type <- unique(input_list$data$SrvLenComps_Type[,f])

      # Aggregated Correlation Parameters
      if(any(srvlen_comp_type == 0) && input_list$data$pop_SrvLenComps_LikeType[f] != 0) {
        if(input_list$data$pop_SrvLenComps_LikeType[f] == 3) {
          map_SrvLen_pop_corr_pars_agg[p,f] <- counter_srvlen_corr_agg
          counter_srvlen_corr_agg <- counter_srvlen_corr_agg + 1 # aggregated
        }
      }

      # Loop through to make sure mapping stuff off correctly
      for(r in 1:input_list$data$n_regions) {
        for(s in 1:input_list$data$n_sexes) {

          # Split by region and sex
          if(any(srvlen_comp_type == 1) && input_list$data$pop_SrvLenComps_LikeType[f] != 0) {
            if(input_list$data$pop_SrvLenComps_LikeType[f] == 3) {
              map_SrvLen_pop_corr_pars[p,r,s,f,1] <- counter_srvlen_corr
              counter_srvlen_corr <- counter_srvlen_corr + 1
            }
          }

          # Joint by sex, split by region
          if(any(srvlen_comp_type == 2) && input_list$data$pop_SrvLenComps_LikeType[f] != 0 && s == 1) {

            # 1dar1 correlation
            if(input_list$data$pop_SrvLenComps_LikeType[f] == 3) {
              map_SrvLen_pop_corr_pars[p,r,1,f,1] <- counter_srvlen_corr
              counter_srvlen_corr <- counter_srvlen_corr + 1
            }

            # 2dar1 correlation
            if(input_list$data$pop_SrvLenComps_LikeType[f] == 4) {
              for(i in 1:2) {
                if(i == 2 && input_list$data$n_sexes == 1) next # skip if we only have 1 sex
                map_SrvLen_pop_corr_pars[p,r,1,f,i] <- counter_srvlen_corr
                counter_srvlen_corr <- counter_srvlen_corr + 1
              } # end i
            } # end if

          }
        } # end s loop
      } # end r loop

    } # end f loop
  } # end p loop

  # Input into mapping list
  input_list$map$SrvLen_pop_corr_pars_agg <- factor(map_SrvLen_pop_corr_pars_agg)
  input_list$map$SrvLen_pop_corr_pars <- factor(map_SrvLen_pop_corr_pars)

  return(input_list)
}

#' Set up observed survey indices and composition data
#'
#' Ingests observed survey index, age composition, and length composition data
#' (both pooled and population-specific) into \code{input_list$data},
#' initialises overdispersion and correlation starting values in
#' \code{input_list$par}, and constructs parameter maps via
#' \code{\link{do_SrvAge_theta_mapping}}, \code{\link{do_SrvLen_theta_mapping}},
#' \code{\link{do_SrvAge_corr_pars_mapping}}, and
#' \code{\link{do_SrvLen_corr_pars_mapping}}. When \code{ISS_SrvAgeComps},
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
#'   Only validated when \code{fit_lengths = 1} in \code{$data}.
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
#'   cell according to \code{pop_SrvAgeComps_Type}.
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
#' @param pop_SrvAgeComps_LikeType Character vector of length
#'   \code{n_srv_fleets} specifying the likelihood for population-specific
#'   survey age compositions. Same options as \code{SrvAgeComps_LikeType}.
#'   Default: \code{"none"} for all fleets.
#' @param pop_SrvLenComps_LikeType Character vector of length
#'   \code{n_srv_fleets} specifying the likelihood for population-specific
#'   survey length compositions. Same options as \code{SrvLenComps_LikeType}.
#'   Default: \code{"none"} for all fleets.
#' @param pop_SrvAgeComps_Type Character vector defining the composition
#'   structure for population-specific survey age compositions. Same format and
#'   options as \code{SrvAgeComps_Type}. Default: \code{"none"} for all fleets
#'   across all years.
#' @param pop_SrvLenComps_Type Character vector defining the composition
#'   structure for population-specific survey length compositions. Same format
#'   and options as \code{SrvLenComps_Type}. Default: \code{"none"} for all
#'   fleets across all years.
#' @param ... Optional named starting values for overdispersion and correlation
#'   parameters. Supported names and default dimensions:
#'   \describe{
#'     \item{\code{ln_SrvAge_theta}}{\code{[n_regions × n_sexes × n_srv_fleets]}.
#'       Default: \code{0}.}
#'     \item{\code{ln_SrvAge_theta_agg}}{\code{[n_srv_fleets]}.
#'       Default: \code{0}.}
#'     \item{\code{SrvAge_corr_pars}}{\code{[n_regions × n_sexes × n_srv_fleets × 2]}.
#'       Default: \code{0.01}.}
#'     \item{\code{SrvAge_corr_pars_agg}}{\code{[n_srv_fleets]}.
#'       Default: \code{0.01}.}
#'     \item{\code{ln_SrvLen_theta}}{\code{[n_regions × n_sexes × n_srv_fleets]}.
#'       Default: \code{0}.}
#'     \item{\code{ln_SrvLen_theta_agg}}{\code{[n_srv_fleets]}.
#'       Default: \code{0}.}
#'     \item{\code{SrvLen_corr_pars}}{\code{[n_regions × n_sexes × n_srv_fleets × 2]}.
#'       Default: \code{0.01}.}
#'     \item{\code{SrvLen_corr_pars_agg}}{\code{[n_srv_fleets]}.
#'       Default: \code{0.01}.}
#'     \item{\code{ln_SrvAge_pop_theta}}{\code{[n_pop × n_regions × n_sexes × n_srv_fleets]}.
#'       Default: \code{0}.}
#'     \item{\code{ln_SrvAge_pop_theta_agg}}{\code{[n_pop × n_srv_fleets]}.
#'       Default: \code{0}.}
#'     \item{\code{SrvAge_pop_corr_pars}}{\code{[n_pop × n_regions × n_sexes × n_srv_fleets × 2]}.
#'       Default: \code{0.01}.}
#'     \item{\code{SrvAge_pop_corr_pars_agg}}{\code{[n_pop × n_srv_fleets]}.
#'       Default: \code{0.01}.}
#'     \item{\code{ln_SrvLen_pop_theta}}{\code{[n_pop × n_regions × n_sexes × n_srv_fleets]}.
#'       Default: \code{0}.}
#'     \item{\code{ln_SrvLen_pop_theta_agg}}{\code{[n_pop × n_srv_fleets]}.
#'       Default: \code{0}.}
#'     \item{\code{SrvLen_pop_corr_pars}}{\code{[n_pop × n_regions × n_sexes × n_srv_fleets × 2]}.
#'       Default: \code{0.01}.}
#'     \item{\code{SrvLen_pop_corr_pars_agg}}{\code{[n_pop × n_srv_fleets]}.
#'       Default: \code{0.01}.}
#'   }
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
#'   \code{pop_SrvAgeComps_LikeType}, \code{pop_SrvLenComps_LikeType},
#'   \code{SrvAgeComps_Type}, \code{SrvLenComps_Type},
#'   \code{pop_SrvAgeComps_Type}, \code{pop_SrvLenComps_Type},
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
                                       pop_SrvAgeComps_LikeType = rep("none", input_list$data$n_srv_fleets),
                                       pop_SrvLenComps_LikeType = rep("none", input_list$data$n_srv_fleets),
                                       pop_SrvAgeComps_Type = paste("none_Year_1-terminal_Fleet_", 1:input_list$data$n_srv_fleets, sep = ''),
                                       pop_SrvLenComps_Type = paste("none_Year_1-terminal_Fleet_", 1:input_list$data$n_srv_fleets, sep = ''),
                                       ...
                                       ) {

  messages_list <<- character(0) # string to attach to for printing messages
  starting_values <- list(...)


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
  check_data_dimensions(pop_SrvAgeComps_LikeType, n_srv_fleets = input_list$data$n_srv_fleets, what = 'pop_SrvAgeComps_LikeType')
  check_data_dimensions(pop_SrvLenComps_LikeType, n_srv_fleets = input_list$data$n_srv_fleets, what = 'pop_SrvLenComps_LikeType')
  if(!all(pop_SrvAgeComps_LikeType %in% c("none", "Multinomial", "Dirichlet-Multinomial", "iid-Logistic-Normal", "1d-Logistic-Normal", "2d-Logistic-Normal")))
    stop("Invalid specification for pop_SrvAgeComps_LikeType Should be either none, Multinomial, Dirichlet-Multinomial, iid-Logistic-Normal, 1d-Logistic-Normal, 2d-Logistic-Normal")
  if(!all(pop_SrvLenComps_LikeType %in% c("none", "Multinomial", "Dirichlet-Multinomial", "iid-Logistic-Normal", "1d-Logistic-Normal", "2d-Logistic-Normal")))
    stop("Invalid specification for pop_SrvLenComps_LikeType Should be either none, Multinomial, Dirichlet-Multinomial, iid-Logistic-Normal, 1d-Logistic-Normal, 2d-Logistic-Normal")

  # checking to make sure defaults are not applied
  if(any(UseSrvAgeComps_pop == 1)) {
    if(is.null(ObsSrvAgeComps_pop)) stop("ObsSrvAgeComps_pop is NULL, but UseSrvAgeComps_pop contains 1s!")
    if(any(str_detect(pop_SrvAgeComps_LikeType, "none"))) warning("pop_SrvAgeComps_LikeType has nones, but UseSrvAgeComps_pop contains 1s! Please verify!")
    if(any(str_detect(pop_SrvAgeComps_Type, "none"))) warning("pop_SrvAgeComps_Type has nones, but UseSrvAgeComps_pop contains 1s! Please verify!")
  }

  if(any(UseSrvLenComps_pop == 1)) {
    if(is.null(ObsSrvLenComps_pop)) stop("ObsSrvLenComps_pop is NULL, but UseSrvAgeComps_pop contains 1s!")
    if(any(str_detect(pop_SrvLenComps_LikeType, "none"))) warning("pop_SrvLenComps_LikeType has nones, but UseSrvAgeComps_pop contains 1s! Please verify!")
    if(any(str_detect(pop_SrvLenComps_Type, "none"))) warning("pop_SrvLenComps_Type has nones, but UseSrvAgeComps_pop contains 1s! Please verify!")
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
  pop_comp_srvage_like_vals <- vector()
  for(f in 1:input_list$data$n_srv_fleets) {
    if(pop_SrvAgeComps_LikeType[f] == 'none') pop_comp_srvage_like_vals <- c(pop_comp_srvage_like_vals, 999)
    if(pop_SrvAgeComps_LikeType[f] == "Multinomial") pop_comp_srvage_like_vals <- c(pop_comp_srvage_like_vals, 0)
    if(pop_SrvAgeComps_LikeType[f] == "Dirichlet-Multinomial") pop_comp_srvage_like_vals <- c(pop_comp_srvage_like_vals, 1)
    if(pop_SrvAgeComps_LikeType[f] == "iid-Logistic-Normal") pop_comp_srvage_like_vals <- c(pop_comp_srvage_like_vals, 2)
    if(pop_SrvAgeComps_LikeType[f] == "1d-Logistic-Normal") pop_comp_srvage_like_vals <- c(pop_comp_srvage_like_vals, 3)
    if(pop_SrvAgeComps_LikeType[f] == "2d-Logistic-Normal") pop_comp_srvage_like_vals <- c(pop_comp_srvage_like_vals, 4)
    collect_message(paste("Population Survey Age Composition Likelihoods", "for survey fleet", f, "specified as:" , pop_SrvAgeComps_LikeType[f]))
  } # end f loop

  # Specifying composition type
  pop_SrvAgeComps_Type_Mat <- array(NA, dim = c(length(input_list$data$years), input_list$data$n_srv_fleets))
  for(i in 1:length(pop_SrvAgeComps_Type)) {

    # Extract out components from list
    tmp <- pop_SrvAgeComps_Type[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    comps_type_tmp <- tmp_vec[1] # get composition type
    fleet <- as.numeric(tmp_vec[5]) # extract fleet index

    # Checking character string
    if(!comps_type_tmp %in% c("agg", "spltRspltS", "spltRjntS", 'none')) stop("pop_SrvAgeComps_Type not specified correctly. Must be one of: agg, spltRspltS, spltRjntS, none")
    if(!fleet %in% c(1:input_list$data$n_srv_fleets)) stop("Invalid fleet specified for pop_SrvAgeComps_Type. This needs to be specified as CompType_Year_x-y_Fleet_x")

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
      if(pop_comp_srvage_like_vals[fleet] == 4) stop("Population Age composition likelihood specified as 2d-Logistic-Normal, but composition type is aggregated. This is not valid.")
      comps_type_val <- 0
    }
    if(comps_type_tmp == "spltRspltS") comps_type_val <- 1
    if(comps_type_tmp == "spltRjntS") comps_type_val <- 2
    if(comps_type_tmp == "none") comps_type_val <- 999

    # input into matrix
    pop_SrvAgeComps_Type_Mat[years,fleet] <- comps_type_val
  } # end i

  if(any(is.na(pop_SrvAgeComps_Type_Mat))) stop("pop_SrvAgeComps_Type is returning an NA. Did you update the year range of pop_SrvAgeComps_Type?")

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
  pop_comp_srvlen_like_vals <- vector()
  for(f in 1:input_list$data$n_srv_fleets) {
    if(pop_SrvLenComps_LikeType[f] == 'none') pop_comp_srvlen_like_vals <- c(pop_comp_srvlen_like_vals, 999)
    if(pop_SrvLenComps_LikeType[f] == "Multinomial") pop_comp_srvlen_like_vals <- c(pop_comp_srvlen_like_vals, 0)
    if(pop_SrvLenComps_LikeType[f] == "Dirichlet-Multinomial") pop_comp_srvlen_like_vals <- c(pop_comp_srvlen_like_vals, 1)
    if(pop_SrvLenComps_LikeType[f] == "iid-Logistic-Normal") pop_comp_srvlen_like_vals <- c(pop_comp_srvlen_like_vals, 2)
    if(pop_SrvLenComps_LikeType[f] == "1d-Logistic-Normal") pop_comp_srvlen_like_vals <- c(pop_comp_srvlen_like_vals, 3)
    if(pop_SrvLenComps_LikeType[f] == "2d-Logistic-Normal") pop_comp_srvlen_like_vals <- c(pop_comp_srvlen_like_vals, 4)
    collect_message(paste("Population Survey Length Composition Likelihoods", "for survey fleet", f, "specified as:" , pop_SrvLenComps_LikeType[f]))
  } # end f loop

  # Specifying composition type
  pop_SrvLenComps_Type_Mat <- array(NA, dim = c(length(input_list$data$years), input_list$data$n_srv_fleets))
  for(i in 1:length(pop_SrvLenComps_Type)) {

    # Extract out components from list
    tmp <- pop_SrvLenComps_Type[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    comps_type_tmp <- tmp_vec[1] # get composition type
    fleet <- as.numeric(tmp_vec[5]) # extract fleet index

    # Checking character string
    if(!comps_type_tmp %in% c("agg", "spltRspltS", "spltRjntS", 'none')) stop("pop_SrvLenComps_Type not specified correctly. Must be one of: agg, spltRspltS, spltRjntS, none")
    if(!fleet %in% c(1:input_list$data$n_srv_fleets)) stop("Invalid fleet specified for pop_SrvLenComps_Type. This needs to be specified as CompType_Year_x-y_Fleet_x")

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
      if(pop_comp_srvlen_like_vals[fleet] == 4) stop("Population Len composition likelihood specified as 2d-Logistic-Normal, but composition type is aggregated. This is not valid.")
      comps_type_val <- 0
    }
    if(comps_type_tmp == "spltRspltS") comps_type_val <- 1
    if(comps_type_tmp == "spltRjntS") comps_type_val <- 2
    if(comps_type_tmp == "none") comps_type_val <- 999

    # input into matrix
    pop_SrvLenComps_Type_Mat[years,fleet] <- comps_type_val
  } # end i

  if(any(is.na(pop_SrvLenComps_Type_Mat))) stop("pop_SrvLenComps_Type is returning an NA. Did you update the year range of pop_SrvLenComps_Type?")

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
            if(pop_SrvAgeComps_Type_Mat[y,f] == 0) ISS_SrvAgeComps_pop[p,1,y,seas,1,f] <- sum(ObsSrvAgeComps_pop[p,,y,seas,,,f])
            # if split by region and sex
            if(pop_SrvAgeComps_Type_Mat[y,f] == 1) ISS_SrvAgeComps_pop[p,,y,seas,,f] <- apply(ObsSrvAgeComps_pop[p,,y,seas,,,f, drop = FALSE], c(2,5), sum)
            # if split by region, joint by sex
            if(pop_SrvAgeComps_Type_Mat[y,f] == 2) ISS_SrvAgeComps_pop[p,,y,seas,1,f] <- apply(ObsSrvAgeComps_pop[p,,y,seas,,,f, drop = FALSE], 2, sum)
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
            if(pop_SrvLenComps_Type_Mat[y,f] == 0) ISS_SrvLenComps_pop[p,1,y,seas,1,f] <- sum(ObsSrvLenComps_pop[p,,y,seas,,,f])
            # if split by region and sex
            if(pop_SrvLenComps_Type_Mat[y,f] == 1) ISS_SrvLenComps_pop[p,,y,seas,,f] <- apply(ObsSrvLenComps_pop[p,,y,seas,,,f, drop = FALSE], c(2,5), sum)
            # if split by region, joint by sex
            if(pop_SrvLenComps_Type_Mat[y,f] == 2) ISS_SrvLenComps_pop[p,,y,seas,1,f] <- apply(ObsSrvLenComps_pop[p,,y,seas,,,f, drop = FALSE], 2, sum)
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
  input_list$data$pop_SrvAgeComps_LikeType <- pop_comp_srvage_like_vals
  input_list$data$pop_SrvLenComps_LikeType <- pop_comp_srvlen_like_vals
  input_list$data$SrvAgeComps_Type <- SrvAgeComps_Type_Mat
  input_list$data$SrvLenComps_Type <- SrvLenComps_Type_Mat
  input_list$data$srv_idx_type <- srv_idx_type_vals
  input_list$data$pop_SrvAgeComps_Type <- pop_SrvAgeComps_Type_Mat
  input_list$data$pop_SrvLenComps_Type <- pop_SrvLenComps_Type_Mat

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

  input_list <- do_SrvAge_theta_mapping(input_list)
  input_list <- do_SrvLen_theta_mapping(input_list)
  input_list <- do_SrvAge_corr_pars_mapping(input_list)
  input_list <- do_SrvLen_corr_pars_mapping(input_list)

  input_list <- do_SrvAge_pop_theta_mapping(input_list)
  input_list <- do_SrvLen_pop_theta_mapping(input_list)
  input_list <- do_SrvAge_pop_corr_pars_mapping(input_list)
  input_list <- do_SrvLen_pop_corr_pars_mapping(input_list)

  # Print Messages ----------------------------------------------------------
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)
}


#' Map survey selectivity fixed-effect parameters
#'
#' Internal helper called by \code{\link{Setup_Mod_Srvsel_and_Q}} to construct
#' the TMB/RTMB factor map for \code{ln_srv_fixed_sel_pars}
#' \code{[n_regions × max_sel_pars × max_sel_blocks × n_sexes × n_srv_fleets]},
#' the log-scale fixed-effect parameters of the survey selectivity functional
#' forms (e.g., \eqn{a_{50}}, \eqn{k}, \eqn{a_{max}}). The number of active
#' parameters per block depends on the selectivity model: exponential = 1,
#' logistic/gamma = 2, double-normal = 6.
#'
#' Parameters are only activated for region–fleet combinations where survey
#' index data are used (\code{sum(UseSrvIdx[r,,,f]) > 0}). Fleet-sharing
#' (\code{"est_shared_f_x"}) is handled in a second pass after all base
#' mappings are established, copying the map from the reference fleet.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists. Requires \code{$data$n_srv_fleets}, \code{$data$n_regions},
#'   \code{$data$n_sexes}, \code{$data$srv_sel_blocks},
#'   \code{$data$srv_sel_model}, and \code{$data$UseSrvIdx}.
#' @param srv_fixed_sel_pars_spec Character vector \code{[n_srv_fleets]}
#'   specifying the sharing structure for fixed-effect selectivity parameters.
#'   One of:
#'   \describe{
#'     \item{\code{"est_all"}}{Separate parameters for each region, sex, block,
#'       and parameter index.}
#'     \item{\code{"est_shared_r"}}{Shared across regions; separate by sex,
#'       block, and parameter index.}
#'     \item{\code{"est_shared_s"}}{Shared across sexes; separate by region,
#'       block, and parameter index.}
#'     \item{\code{"est_shared_r_s"}}{Shared across both regions and sexes;
#'       separate by block and parameter index.}
#'     \item{\code{"est_shared_f_x"}}{Copy the map from fleet \code{x}.
#'       Fleet \code{x} must not itself use \code{"est_shared_f_x"}.}
#'     \item{\code{"fix"}}{All parameters fixed at starting values
#'       (mapped to \code{NA}).}
#'   }
#'
#' @return The input \code{input_list} with \code{$map$ln_srv_fixed_sel_pars}
#'   set to a factor vector. Active parameters receive sequential integer
#'   indices; inactive region–fleet combinations and fixed parameters are
#'   \code{NA}.
#'
#'
#' @keywords internal
do_srv_fixed_sel_pars_mapping <- function(input_list, srv_fixed_sel_pars_spec) {

  # Initialize counter and mapping array for fixed effects survey selectivity
  srv_fixed_sel_pars_counter <- 1
  map_srv_fixed_sel_pars <- input_list$par$ln_srv_fixed_sel_pars
  map_srv_fixed_sel_pars[] <- NA

  for(f in 1:input_list$data$n_srv_fleets) {

    # Validate Options
    if(!srv_fixed_sel_pars_spec[f] %in% c("est_all", "est_shared_r", "est_shared_r_s", "fix", "est_shared_s") &&
       !stringr::str_detect(srv_fixed_sel_pars_spec[f], "est_shared_f_\\d+"))
      stop("srv_fixed_sel_pars_spec not correctly specfied. Should be one of these: est_all, est_shared_r, est_shared_r_s, est_shared_s, fix, or est_shared_f_# (where # is fleet number)")

    # Skip fleet sharing specs in first pass
    if(stringr::str_detect(srv_fixed_sel_pars_spec[f], "est_shared_f")) next

    for(r in 1:input_list$data$n_regions) {

      # Only add a counter if caatches are avaliable in some years for a given region and fleet combination
      if(sum(input_list$data$UseSrvIdx[r,,,f]) > 0) {

        # Extract number of survey selectivity blocks
        srvsel_blocks_tmp <- unique(as.vector(input_list$data$srv_sel_blocks[r,,f]))

        for(s in 1:input_list$data$n_sexes) {
          for(b in 1:length(srvsel_blocks_tmp)) {

            block_years <- which(input_list$data$srv_sel_blocks[r,,f] == srvsel_blocks_tmp[b]) # figure out block years
            sel_model_this_block <- unique(input_list$data$srv_sel_model[r, block_years, f]) # get selectivity form for a given block
            if(length(sel_model_this_block) > 1) stop("Block ", srvsel_blocks_tmp[b], " for fleet ", f, " region ", r, " has multiple selectivity models assigned to it")

            # determine maximum selectivity parameters
            if(sel_model_this_block == 2) max_sel_pars <- 1 # exponential
            if(sel_model_this_block %in% c(0,1,3)) max_sel_pars <- 2 # logistic or gamma
            if(sel_model_this_block == 4) max_sel_pars <- 6 # double normal

            for(i in 1:max_sel_pars) {

              # Estimate all selectivity fixed effects parameters within the constraints of the defined blocks
              if(srv_fixed_sel_pars_spec[f] == "est_all") {
                map_srv_fixed_sel_pars[r,i,b,s,f] <- srv_fixed_sel_pars_counter
                srv_fixed_sel_pars_counter <- srv_fixed_sel_pars_counter + 1
              } # end if

              # Estimating process error parameters shared across regions (but unique for each sex, fleet, parameter)
              if(srv_fixed_sel_pars_spec[f] == 'est_shared_r' && r == 1) {
                for(rr in 1:input_list$data$n_regions) {
                  if(srvsel_blocks_tmp[b] %in% input_list$data$srv_sel_blocks[rr,,f]) {
                    map_srv_fixed_sel_pars[rr, i, b, s, f] <- srv_fixed_sel_pars_counter
                  } # end if
                } # end rr loop
                srv_fixed_sel_pars_counter <- srv_fixed_sel_pars_counter + 1
              } # end if

              # Estimating process error parameters shared across sexes (but unique for each region, fleet, parameter)
              if(srv_fixed_sel_pars_spec[f] == 'est_shared_s' && s == 1) {
                for(ss in 1:input_list$data$n_sexes) {
                  map_srv_fixed_sel_pars[r, i, b, ss, f] <- srv_fixed_sel_pars_counter
                } # end ss loop
                srv_fixed_sel_pars_counter <- srv_fixed_sel_pars_counter + 1
              } # end if

              # Estimating process error parameters shared across regions and sexes (but unique for each fleet, parameter)
              if(srv_fixed_sel_pars_spec[f] == 'est_shared_r_s' && r == 1 && s == 1) {
                for(rr in 1:input_list$data$n_regions) {
                  for(ss in 1:input_list$data$n_sexes) {
                    if(srvsel_blocks_tmp[b] %in% input_list$data$srv_sel_blocks[rr,,f]) {
                      map_srv_fixed_sel_pars[rr, i, b, ss, f] <- srv_fixed_sel_pars_counter
                    } # end if
                  } # end ss loop
                } #end rr loop
                srv_fixed_sel_pars_counter <- srv_fixed_sel_pars_counter + 1
              } # end if

            } # end i loop
          } # end b loop
        } # end s loop
      } # end if statement
    } # end r loop

    # fix all parameters
    if(srv_fixed_sel_pars_spec[f] == "fix") map_srv_fixed_sel_pars[,,,,f] <- NA
    collect_message("srv_fixed_sel_pars_spec is specified as: ", srv_fixed_sel_pars_spec[f], " for survey fleet ", f)

  } # end f loop

  # Handle fleet sharing after all base mappings are established
  for(f in 1:input_list$data$n_srv_fleets) {
    if(stringr::str_detect(srv_fixed_sel_pars_spec[f], "est_shared_f")) {
      # extract fleet sharing index
      flt_shared <- as.numeric(unlist(stringr::str_extract_all(srv_fixed_sel_pars_spec[f], "\\d+")))

      # Validate options here
      if(flt_shared > input_list$data$n_srv_fleets || flt_shared < 1) stop("Fleet sharing specification 'est_shared_f", flt_shared, "' for fleet ", f, " references invalid fleet number.")
      if(stringr::str_detect(srv_fixed_sel_pars_spec[flt_shared], "est_shared_f")) stop("Fleet ", f, " cannot share with fleet ", flt_shared, " because fleet ", flt_shared, " is self-sharing parameters, which does not make sense.")

      # Copy mapping from reference fleet
      map_srv_fixed_sel_pars[,,,,f] <- map_srv_fixed_sel_pars[,,,,flt_shared]
      collect_message("srv_fixed_sel_pars_spec is specified as: ", srv_fixed_sel_pars_spec[f], " for survey fleet ", f, " (sharing with fleet ", flt_shared, ")")
    } # end if statement
  } # end f loop

  # input into mapping list
  input_list$map$ln_srv_fixed_sel_pars <- factor(map_srv_fixed_sel_pars)

  return(input_list)
}


#' Map survey catchability parameters
#'
#' Internal helper called by \code{\link{Setup_Mod_Srvsel_and_Q}} to construct
#' the TMB/RTMB factor map for \code{ln_srv_q}
#' \code{[n_regions × max_q_blocks × n_srv_fleets]}, the log-scale survey
#' catchability parameters. Parameters are only activated for region–fleet
#' combinations where survey index data are used
#' (\code{sum(UseSrvIdx[r,,,f]) > 0}); unused combinations are mapped to
#' \code{NA}.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists. Requires \code{$data$n_srv_fleets}, \code{$data$n_regions},
#'   \code{$data$srv_q_blocks}, and \code{$data$UseSrvIdx}.
#' @param srv_q_spec Character vector \code{[n_srv_fleets]} specifying the
#'   sharing structure for catchability. One of:
#'   \describe{
#'     \item{\code{"est_all"}}{Separate catchability per region and block.}
#'     \item{\code{"est_shared_r"}}{Single catchability per block shared across
#'       all regions. Block membership is checked per region before assignment.}
#'     \item{\code{"fix"}}{All catchability parameters fixed at starting values
#'       (mapped to \code{NA}).}
#'   }
#'
#' @return The input \code{input_list} with \code{$map$ln_srv_q} set to a
#'   factor vector. Active parameters receive sequential integer indices;
#'   unused region–fleet combinations and fixed parameters are \code{NA}.
#'
#'
#' @keywords internal
do_srv_q_mapping <- function(input_list, srv_q_spec) {

  # Initialize counter and mapping array for survey catchability
  srv_q_counter <- 1
  map_srv_q <- input_list$par$ln_srv_q
  map_srv_q[] <- NA

  for(f in 1:input_list$data$n_srv_fleets) {

    # Validate options
    if(!is.null(srv_q_spec)) {
      if(!srv_q_spec[f] %in% c("est_all", "est_shared_r", "fix"))
        stop("srv_q_spec not correctly specfied. Should be one of these: est_all, est_shared_r, fix")
    }

    for(r in 1:input_list$data$n_regions) {

      if(sum(input_list$data$UseSrvIdx[r,,,f]) == 0) {
        map_srv_q[r,,f] <- NA # fix parameters if we are not using survey indices for these fleets and regions
      } else {

        # Extract number of survey catchability blocks
        srvq_blocks_tmp <- unique(as.vector(input_list$data$srv_q_blocks[r,,f]))

        for(b in 1:length(srvq_blocks_tmp)) {

          # Estimate for all regions
          if(srv_q_spec[f] == 'est_all') {
            map_srv_q[r,b,f] <- srv_q_counter
            srv_q_counter <- srv_q_counter + 1
          } # end if

          # Estimate but share q across regions
          if(srv_q_spec[f] == 'est_shared_r' && r == 1) {
            for(rr in 1:input_list$data$n_regions) {
              if(srvq_blocks_tmp[b] %in% input_list$data$srv_q_blocks[rr,,f]) {
                map_srv_q[rr, b, f] <- srv_q_counter
              } # end if
            } # end rr loop
            srv_q_counter <- srv_q_counter + 1
          } # end if

        } # end b loop
      } # end else loop
    } # end r loop

    # fix all parameters
    if(srv_q_spec[f] == 'fix') map_srv_q[,,f] <- NA
    collect_message("srv_q_spec is specified as: ", srv_q_spec[f], " for survey fleet ", f)
  } # end f loop

  # input into mapping list
  input_list$map$ln_srv_q <- factor(map_srv_q)

  return(input_list)
}

#' Map survey selectivity process error (variance/correlation) parameters
#'
#' Internal helper called by \code{\link{Setup_Mod_Srvsel_and_Q}} to construct
#' the TMB/RTMB factor map for \code{srvsel_pe_pars}
#' \code{[n_regions × max(max_sel_pars, 4) × n_sexes × n_srv_fleets]}, the
#' hyperparameters governing the variance and correlation structure of
#' continuous time-varying survey selectivity. The second dimension is padded
#' to at least 4 to accommodate the maximum number of process error parameters
#' across all semi-parametric forms:
#' \describe{
#'   \item{iid or random walk (\code{cont_tv_srv_sel} 1–2)}{Up to
#'     \code{max_sel_pars} variance parameters (1, 2, or 6 depending on
#'     selectivity model).}
#'   \item{3D GMRF (\code{cont_tv_srv_sel} 3–4)}{4 parameters: age partial
#'     correlation (index 1), year partial correlation (index 2), cohort
#'     partial correlation (index 3), log-sigma (index 4).}
#'   \item{2D AR1 (\code{cont_tv_srv_sel} 5)}{3 parameters at indices 1, 2,
#'     and 4 (bin, year, log-sigma); index 3 is always \code{NA}.}
#' }
#'
#' When \code{corr_opt_semipar} is non-\code{NULL}, specified correlation
#' indices are set to \code{NA} and the remaining active parameters are
#' re-indexed sequentially. Fleet-sharing (\code{"est_shared_f_x"}) is
#' handled in a second pass. Parameters are automatically fixed for fleets
#' with no continuous time-variation (\code{cont_tv_srv_sel = 0}) or no
#' active survey index data.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists. Requires \code{$data$n_srv_fleets}, \code{$data$n_regions},
#'   \code{$data$n_sexes}, \code{$data$cont_tv_srv_sel},
#'   \code{$data$srv_sel_model}, \code{$data$UseSrvIdx}, and
#'   \code{$data$Selex_Type}.
#' @param srvsel_pe_pars_spec Character vector \code{[n_srv_fleets]} specifying
#'   the sharing structure. Same options as
#'   \code{\link{do_srvsel_devs_mapping}}: \code{"est_all"},
#'   \code{"est_shared_r"}, \code{"est_shared_s"}, \code{"est_shared_r_s"},
#'   \code{"est_shared_f_x"}, \code{"fix"}, \code{"none"}.
#' @param corr_opt_semipar Character vector \code{[n_srv_fleets]} or
#'   \code{NULL}. Specifies which correlation components to suppress for
#'   3D GMRF or 2D AR1 likelihoods by setting the corresponding parameter
#'   indices to \code{NA}. Valid values: \code{NA}, \code{"corr_zero_y"},
#'   \code{"corr_zero_b"}, \code{"corr_zero_y_b"}, \code{"corr_zero_c"}
#'   (3D GMRF only), \code{"corr_zero_y_c"} (3D GMRF only),
#'   \code{"corr_zero_b_c"} (3D GMRF only), \code{"corr_zero_y_b_c"}
#'   (3D GMRF only). After suppression, remaining non-\code{NA} parameters
#'   are re-indexed sequentially.
#'
#' @return The input \code{input_list} with \code{$map$srvsel_pe_pars} set to
#'   a factor vector. Active parameters receive sequential integer indices;
#'   fixed, inactive, and correlation-suppressed parameters are \code{NA}.
#'
#'
#' @keywords internal
do_srvsel_pe_pars_mapping <- function(input_list, srvsel_pe_pars_spec, corr_opt_semipar) {

  # Initialize counter and mapping array for survey process errors
  srvsel_pe_pars_counter <- 1 # initalize counter
  map_srvsel_pe_pars <- input_list$par$srvsel_pe_pars # initalize array
  map_srvsel_pe_pars[] <- NA

  # Survey process error parameters
  for(f in 1:input_list$data$n_srv_fleets) {

    # Validate options
    if(!is.null(srvsel_pe_pars_spec)) {
      if(!srvsel_pe_pars_spec[f] %in% c("fix", "none", "est_all", "est_shared_r", "est_shared_s", "est_shared_r_s") &&
         !stringr::str_detect(srvsel_pe_pars_spec[f], "est_shared_f_\\d+"))
        stop("srvsel_pe_pars_spec not correctly specfied. Should be one of these: est_all, est_shared_r, est_shared_r_s, est_shared_s, fix, or est_shared_f_# (where # is fleet number)")
    }

    # Skip fleet sharing specs in first pass
    if(!is.null(srvsel_pe_pars_spec)) if(stringr::str_detect(srvsel_pe_pars_spec[f], "est_shared_f")) next

    for(r in 1:input_list$data$n_regions) {

      # if no time-variation, then fix all parameters for this fleet
      if(input_list$data$cont_tv_srv_sel[r,f] == 0 || sum(input_list$data$UseSrvIdx[r,,,f]) == 0) {
        map_srvsel_pe_pars[r,,,f] <- NA
      } else { # if we have time-variation

        # Figure out max number of selectivity parameters for a given region and fleet
        if(unique(input_list$data$srv_sel_model[r,,f]) %in% 2) max_sel_pars <- 1 # exponential
        if(unique(input_list$data$srv_sel_model[r,,f]) %in% c(0,1,3)) max_sel_pars <- 2 # logistic or gamma
        if(unique(input_list$data$srv_sel_model[r,,f]) == 4) max_sel_pars <- 6 # double normal

        for(s in 1:input_list$data$n_sexes) {

          # If iid time-variation or random walk for this fleet
          if(input_list$data$cont_tv_srv_sel[r,f] %in% c(1,2)) {

            for(i in 1:max_sel_pars) {

              # either fixing parameters or not used for a given fleet
              if(srvsel_pe_pars_spec[f] %in% c("none", "fix")) map_srvsel_pe_pars[r,i,s,f] <- NA

              # Estimating all parameters separately (unique for each region, sex, fleet, parameter)
              if(srvsel_pe_pars_spec[f] == "est_all") {
                map_srvsel_pe_pars[r,i,s,f] <- srvsel_pe_pars_counter
                srvsel_pe_pars_counter <- srvsel_pe_pars_counter + 1
              } # end est_all

              # Estimating process error parameters shared across regions (but unique for each sex, fleet, parameter)
              if(srvsel_pe_pars_spec[f] == 'est_shared_r' && r == 1) {
                map_srvsel_pe_pars[,i,s,f] <- srvsel_pe_pars_counter
                srvsel_pe_pars_counter <- srvsel_pe_pars_counter + 1
              }

              # Estimating process error parameters shared across sexes (but unique for each region, fleet, parameter)
              if(srvsel_pe_pars_spec[f] == 'est_shared_s' && s == 1) {
                map_srvsel_pe_pars[r,i,,f] <- srvsel_pe_pars_counter
                srvsel_pe_pars_counter <- srvsel_pe_pars_counter + 1
              }

              # Estimating process error parameters shared across regions and sexes (but unique for each fleet, parameter)
              if(srvsel_pe_pars_spec[f] == 'est_shared_r_s' && r == 1 && s == 1) {
                map_srvsel_pe_pars[,i,,f] <- srvsel_pe_pars_counter
                srvsel_pe_pars_counter <- srvsel_pe_pars_counter + 1
              }

            } # end i loop
          } # end iid or random walk variation

          # If 3d gmrf or 2dar1
          if(input_list$data$cont_tv_srv_sel[r,f] %in% c(3,4,5)) {

            # Set up indexing to loop through
            if(input_list$data$cont_tv_srv_sel[r,f] %in% c(3,4)) idx = 1:4 # 3dgmrf (1 = pcorr_age, 2 = pcorr_year, 3= pcorr_cohort, 4 = log_sigma)
            if(input_list$data$cont_tv_srv_sel[r,f] %in% c(5)) idx = c(1,2,4) # 2dar1 (1 = pcorr_bin, 2 = pcorr_year, 4 = log_sigma)
            if(input_list$data$cont_tv_srv_sel[r,f] %in% c(3,4) && input_list$data$Selex_Type == 1) stop("Cohort-based selectivity deviations are specified, but selectivity is specified as length-based. Please choose another deviation form!")

            for(i in idx) {

              # either fixing parameters or not used for a given fleet
              if(srvsel_pe_pars_spec[f] %in% c("none", "fix")) map_srvsel_pe_pars[r,i,s,f] <- NA

              # Estimating all process error parameters
              if(srvsel_pe_pars_spec[f] == "est_all") {
                map_srvsel_pe_pars[r,i,s,f] <- srvsel_pe_pars_counter
                srvsel_pe_pars_counter <- srvsel_pe_pars_counter + 1
              } # end est_all

              # Estimating process error parameters shared across regions (but unique for each sex, fleet, parameter)
              if(srvsel_pe_pars_spec[f] == 'est_shared_r' && r == 1) {
                map_srvsel_pe_pars[,i,s,f] <- srvsel_pe_pars_counter
                srvsel_pe_pars_counter <- srvsel_pe_pars_counter + 1
              }

              # Estimating process error parameters shared across sexes (but unique for each region, fleet, parameter)
              if(srvsel_pe_pars_spec[f] == 'est_shared_s' && s == 1) {
                map_srvsel_pe_pars[r,i,,f] <- srvsel_pe_pars_counter
                srvsel_pe_pars_counter <- srvsel_pe_pars_counter + 1
              }

              # Estimating process error parameters shared across regions and sexes (but unique for each fleet, parameter)
              if(srvsel_pe_pars_spec[f] == 'est_shared_r_s' && r == 1 && s == 1) {
                map_srvsel_pe_pars[,i,,f] <- srvsel_pe_pars_counter
                srvsel_pe_pars_counter <- srvsel_pe_pars_counter + 1
              }

            } # end i loop

            # Options to set correaltions to 0 for 3dgmrf
            if(!is.null(corr_opt_semipar)) {

              opt <- input_list$data$cont_tv_srv_sel[r,f] # get random effects options

              # Validate options
              if(!corr_opt_semipar[f] %in% c(NA, "corr_zero_y", "corr_zero_b", "corr_zero_y_b", "corr_zero_c", "corr_zero_y_c", "corr_zero_b_c", "corr_zero_y_b_c"))
                stop("corr_opt_semipar not correctly specfied. Should be one of these: corr_zero_y, corr_zero_b, corr_zero_y_b, corr_zero_c, corr_zero_y_c, corr_zero_b_c, corr_zero_y_b_c, NA")
              if(opt == 5 && corr_opt_semipar[f] %in% c("corr_zero_c","corr_zero_y_c","corr_zero_b_c","corr_zero_y_b_c"))
                stop("Invalid corr_opt_semipar for 2dar1 (opt=5): cohort correlations are not allowed.")

              if (opt %in% c(3,4,5)) {
                # 2d and 3d options
                if (corr_opt_semipar[f] == "corr_zero_y")    map_srvsel_pe_pars[,2,,f]     <- NA
                if (corr_opt_semipar[f] == "corr_zero_b")    map_srvsel_pe_pars[,1,,f]     <- NA
                if (corr_opt_semipar[f] == "corr_zero_y_b")  map_srvsel_pe_pars[,1:2,,f]   <- NA
              }

              if(opt %in% c(3,4)) {
                # 3d gmrf options only (adds the cohort dimension)
                if (corr_opt_semipar[f] == "corr_zero_c")      map_srvsel_pe_pars[,3,,f]   <- NA
                if (corr_opt_semipar[f] == "corr_zero_y_c")    map_srvsel_pe_pars[,2:3,,f] <- NA
                if (corr_opt_semipar[f] == "corr_zero_b_c")    map_srvsel_pe_pars[,c(1,3),,f] <- NA
                if (corr_opt_semipar[f] == "corr_zero_y_b_c")  map_srvsel_pe_pars[,1:3,,f] <- NA
              }

              # Reset numbering for mapping off correlation parameters for clarity
              non_na_positions <- which(!is.na(map_srvsel_pe_pars))
              map_srvsel_pe_pars[non_na_positions] <- seq_along(non_na_positions)
              collect_message("corr_opt_semipar is specified as: ", corr_opt_semipar[f], "for survey fleet", f)

            }
          } # end if 3d gmrf marginal or conditional variance

          # fix all parameters
          if(srvsel_pe_pars_spec[f] == "fix") map_srvsel_pe_pars[r,,s,f] <- NA

        } # end s loop
      } # end else
    } # end r loop

    if(!is.null(srvsel_pe_pars_spec)) collect_message("srvsel_pe_pars_spec is specified as: ", srvsel_pe_pars_spec[f], "for survey fleet", f)

  } # end f loop

  # Handle fleet sharing after all base mappings are established
  for(f in 1:input_list$data$n_srv_fleets) {
    if(stringr::str_detect(srvsel_pe_pars_spec[f], "est_shared_f") && !is.null(srvsel_pe_pars_spec)) {
      # extract fleet sharing index
      flt_shared <- as.numeric(unlist(stringr::str_extract_all(srvsel_pe_pars_spec[f], "\\d+")))

      # Validate options here
      if(flt_shared > input_list$data$n_srv_fleets || flt_shared < 1) stop("Fleet sharing specification 'est_shared_f", flt_shared, "' for fleet ", f, " references invalid fleet number.")
      if(stringr::str_detect(srvsel_pe_pars_spec[flt_shared], "est_shared_f")) stop("Fleet ", f, " cannot share with fleet ", flt_shared, " because fleet ", flt_shared, " is self-sharing parameters, which does not make sense.")

      # Copy mapping from reference fleet
      map_srvsel_pe_pars[,,,f] <- map_srvsel_pe_pars[,,,flt_shared]
      collect_message("srvsel_pe_pars_spec is specified as: ", srvsel_pe_pars_spec[f], " for survey fleet ", f, " (sharing with fleet ", flt_shared, ")")
    } # end if statement
  } # end f loop

  # input into mapping list
  input_list$map$srvsel_pe_pars <- factor(map_srvsel_pe_pars)

  return(input_list)
}

#' Map survey selectivity deviation parameters
#'
#' Internal helper called by \code{\link{Setup_Mod_Srvsel_and_Q}} to construct
#' the TMB/RTMB factor map for \code{ln_srvsel_devs}
#' \code{[n_regions × (n_years + n_proj_yrs_devs) × n_bins × n_sexes × n_srv_fleets]},
#' where \code{n_bins} is \code{n_ages} for age-based selectivity or
#' \code{n_lens} for length-based selectivity. For iid and random walk forms,
#' the bin dimension indexes selectivity parameters (1, 2, or 6 depending on
#' functional form); for semi-parametric forms (3D GMRF, 2D AR1), it indexes
#' age or length bins directly.
#'
#' Age-sharing options (\code{"est_shared_a"} and variants) are only valid for
#' semi-parametric forms and require \code{srvsel_devs_shared_ages} to define
#' which bin groups share a common deviation. Fleet-sharing
#' (\code{"est_shared_f_x"}) is handled in a second pass. Parameters are
#' automatically fixed for fleets with no continuous time-variation or no
#' active survey index data. The integer-valued map is also stored in
#' \code{$data$map_ln_srvsel_devs} for use within the TMB objective function.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists. Requires \code{$data$n_srv_fleets}, \code{$data$n_regions},
#'   \code{$data$n_sexes}, \code{$data$cont_tv_srv_sel},
#'   \code{$data$srv_sel_model}, \code{$data$UseSrvIdx},
#'   \code{$data$n_proj_yrs_devs}, and \code{$data$Selex_Type}.
#' @param srv_sel_devs_spec Character vector \code{[n_srv_fleets]} specifying
#'   the sharing structure for selectivity deviations. One of:
#'   \describe{
#'     \item{\code{"est_all"}}{Separate deviation time series per region, sex,
#'       bin, and fleet.}
#'     \item{\code{"est_shared_r"}}{Shared across regions; separate by sex
#'       and bin.}
#'     \item{\code{"est_shared_s"}}{Shared across sexes; separate by region
#'       and bin.}
#'     \item{\code{"est_shared_r_s"}}{Shared across regions and sexes.}
#'     \item{\code{"est_shared_a"}}{Shared across bin groups defined by
#'       \code{srvsel_devs_shared_ages}. Semi-parametric forms only.}
#'     \item{\code{"est_shared_r_a"}}{Shared across regions and bin groups.
#'       Semi-parametric forms only.}
#'     \item{\code{"est_shared_a_s"}}{Shared across bin groups and sexes.
#'       Semi-parametric forms only.}
#'     \item{\code{"est_shared_r_a_s"}}{Shared across regions, bin groups,
#'       and sexes. Semi-parametric forms only.}
#'     \item{\code{"est_shared_f_x"}}{Copy the map from fleet \code{x}.
#'       Fleet \code{x} must not itself use \code{"est_shared_f_x"}.}
#'     \item{\code{"fix"}/\code{"none"}}{All deviations fixed at zero
#'       (mapped to \code{NA}).}
#'   }
#' @param srvsel_devs_shared_ages List of integer vectors defining bin groups
#'   for age/length sharing. Each element specifies the bin indices within one
#'   group (e.g., \code{list(1:5, 6:10, 11:30)}). Required when
#'   \code{srv_sel_devs_spec} includes \code{"est_shared_a"} or its
#'   variants; ignored otherwise.
#'
#' @return The input \code{input_list} with \code{$map$ln_srvsel_devs} set to
#'   a factor vector and \code{$data$map_ln_srvsel_devs} set to the equivalent
#'   integer array \code{[n_regions × (n_years + n_proj_yrs_devs) × n_bins ×
#'   n_sexes × n_srv_fleets]}. Active parameters receive sequential integer
#'   indices; inactive parameters are \code{NA}.
#'
#'
#' @keywords internal
do_srvsel_devs_mapping <- function(input_list, srv_sel_devs_spec, srvsel_devs_shared_ages) {

  # Initialize counter and mapping array for survey selectivity deviations
  srvsel_devs_counter <- 1
  map_srvsel_devs <- input_list$par$ln_srvsel_devs
  map_srvsel_devs[] <- NA

  for(r in 1:input_list$data$n_regions) {
    for(f in 1:input_list$data$n_srv_fleets) {

      # Validate options
      if(!is.null(srv_sel_devs_spec)) {
        if(!srv_sel_devs_spec[f] %in% c("fix", "none", "est_all", "est_shared_r", "est_shared_s", "est_shared_r_s", "est_shared_a", "est_shared_r_a", "est_shared_r_a_s", "est_shared_a_s") &&
           !stringr::str_detect(srv_sel_devs_spec[f], "est_shared_f_\\d+"))
          stop("srv_sel_devs_spec not correctly specfied. Should be one of these: est_all, est_shared_r, est_shared_r_s, est_shared_s, est_shared_a, est_shared_r_a, est_shared_r_a_s, est_shared_r_s, fix, or est_shared_f_# (where # is fleet number)")
        if(srv_sel_devs_spec[f] %in% c("est_shared_a", "est_shared_r_a", "est_shared_r_a_s", "est_shared_a_s") &&
           !input_list$data$cont_tv_srv_sel[r,f] %in% c(3,4,5)) stop("Sharing age deviations with iid or random walk parametric forms is not supported!")
      }

      # Skip fleet sharing specs in first pass
      if(!is.null(srv_sel_devs_spec)) if(stringr::str_detect(srv_sel_devs_spec[f], "est_shared_f")) next

      for(s in 1:input_list$data$n_sexes) {
        for(y in 1:(length(input_list$data$years) + input_list$data$n_proj_yrs_devs)) {

          # if no time-variation, then fix all parameters for this fleet
          if(input_list$data$cont_tv_srv_sel[r,f] == 0 || sum(input_list$data$UseSrvIdx[r,,,f]) == 0) {
            map_srvsel_devs[r,y,,s,f] <- NA
          } else {

            # Figure out max number of selectivity parameters for a given region and fleet
            if(unique(input_list$data$srv_sel_model[r,,f]) %in% 2) max_sel_pars <- 1 # exponential
            if(unique(input_list$data$srv_sel_model[r,,f]) %in% c(0,1,3)) max_sel_pars <- 2 # logistic or gamma
            if(unique(input_list$data$srv_sel_model[r,,f]) == 4) max_sel_pars <- 6 # double normal

            # If iid or random walk time-variation for this fleet
            if(input_list$data$cont_tv_srv_sel[r,f] %in% c(1,2)) {

              for(i in 1:max_sel_pars) {
                # Estimating all selectivity deviations across regions, sexes, fleets, and parameter
                if(srv_sel_devs_spec[f] == 'est_all') {
                  map_srvsel_devs[r,y,i,s,f] <- srvsel_devs_counter
                  srvsel_devs_counter <- srvsel_devs_counter + 1
                }

                # Estimating selectivity deviations across sexes, fleets, and parameters, but shared across regions
                if(srv_sel_devs_spec[f] == 'est_shared_r' && r == 1) {
                  map_srvsel_devs[,y,i,s,f] <- srvsel_devs_counter
                  srvsel_devs_counter <- srvsel_devs_counter + 1
                }

                # Estimating selectivity deviations across regions, fleets, and parameters, but shared across sexes
                if(srv_sel_devs_spec[f] == 'est_shared_s' && s == 1) {
                  map_srvsel_devs[r,y,i,,f] <- srvsel_devs_counter
                  srvsel_devs_counter <- srvsel_devs_counter + 1
                }

                # Estimating selectivity deviations across fleets, and parameters, but shared across sexes and regions
                if(srv_sel_devs_spec[f] == 'est_shared_r_s' && r == 1 && s == 1) {
                  map_srvsel_devs[,y,i,,f] <- srvsel_devs_counter
                  srvsel_devs_counter <- srvsel_devs_counter + 1
                }

              } # end i loop
            } # end iid or random walk variation

            # If 3d gmrf for this fleet
            if(input_list$data$cont_tv_srv_sel[r,f] %in% c(3,4,5)) {

              for(i in 1:length(input_list$data$ages)) {
                # Estimating all selectivity deviations across regions, years and bins
                if(srv_sel_devs_spec[f] == 'est_all') {
                  map_srvsel_devs[r,y,i,s,f] <- srvsel_devs_counter
                  srvsel_devs_counter <- srvsel_devs_counter + 1
                }

                # Estimating all selectivity deviations across years and bins, but shared across regions
                if(srv_sel_devs_spec[f] == 'est_shared_r' && r == 1) {
                  map_srvsel_devs[,y,i,s,f] <- srvsel_devs_counter
                  srvsel_devs_counter <- srvsel_devs_counter + 1
                }

                # Estimating all selectivity deviations across years and bins, but shared across sexes
                if(srv_sel_devs_spec[f] == 'est_shared_s' && s == 1) {
                  map_srvsel_devs[r,y,i,,f] <- srvsel_devs_counter
                  srvsel_devs_counter <- srvsel_devs_counter + 1
                }

                # Estimating all selectivity deviations across years and bins, but shared across sexes and regions
                if(srv_sel_devs_spec[f] == 'est_shared_r_s' && s == 1 && r == 1) {
                  map_srvsel_devs[,y,i,,f] <- srvsel_devs_counter
                  srvsel_devs_counter <- srvsel_devs_counter + 1
                }

                if(srv_sel_devs_spec[f] == 'est_shared_a') {
                  for(k in 1:length(srvsel_devs_shared_ages)) {
                    map_srvsel_devs[r,y,srvsel_devs_shared_ages[[k]],s,f] <- srvsel_devs_counter
                    srvsel_devs_counter <- srvsel_devs_counter + 1
                  } # end k loop
                }

                if(srv_sel_devs_spec[f] == 'est_shared_r_a' && r == 1) {
                  for(k in 1:length(srvsel_devs_shared_ages)) {
                    map_srvsel_devs[,y,srvsel_devs_shared_ages[[k]],s,f] <- srvsel_devs_counter
                    srvsel_devs_counter <- srvsel_devs_counter + 1
                  } # end k loop
                }

                if(srv_sel_devs_spec[f] == 'est_shared_a_s' && s == 1) {
                  for(k in 1:length(srvsel_devs_shared_ages)) {
                    map_srvsel_devs[r,y,srvsel_devs_shared_ages[[k]],,f] <- srvsel_devs_counter
                    srvsel_devs_counter <- srvsel_devs_counter + 1
                  } # end k loop
                }

                if(srv_sel_devs_spec[f] == 'est_shared_r_a_s' && s == 1 && r == 1) {
                  for(k in 1:length(srvsel_devs_shared_ages)) {
                    map_srvsel_devs[,y,srvsel_devs_shared_ages[[k]],,f] <- srvsel_devs_counter
                    srvsel_devs_counter <- srvsel_devs_counter + 1
                  } # end k loop
                }

              } # end i loop
            } # end 3d gmrf

          } # end else
        } # end y loop
      } # end s loop

      if(!is.null(srv_sel_devs_spec)) collect_message("srv_sel_devs_spec is specified as: ", srv_sel_devs_spec[f], "for survey fleet", f, "and region ", r)

    } # end f loop
  } # end r loop

  # Handle fleet sharing after all base mappings are established
  for(f in 1:input_list$data$n_srv_fleets) {
    if(stringr::str_detect(srv_sel_devs_spec[f], "est_shared_f") && !is.null(srv_sel_devs_spec)) {
      # extract fleet sharing index
      flt_shared <- as.numeric(unlist(stringr::str_extract_all(srv_sel_devs_spec[f], "\\d+")))

      # Validate options here
      if(flt_shared > input_list$data$n_srv_fleets || flt_shared < 1) stop("Fleet sharing specification 'est_shared_f", flt_shared, "' for fleet ", f, " references invalid fleet number.")
      if(stringr::str_detect(srv_sel_devs_spec[flt_shared], "est_shared_f")) stop("Fleet ", f, " cannot share with fleet ", flt_shared, " because fleet ", flt_shared, " is self-sharing parameters, which does not make sense.")

      # Copy mapping from reference fleet
      map_srvsel_devs[,,,,f] <- map_srvsel_devs[,,,,flt_shared]
      collect_message("srv_sel_devs_spec is specified as: ", srv_sel_devs_spec[f], " for survey fleet ", f, " (sharing with fleet ", flt_shared, ")")
    } # end if statement
  } # end f loop

  # input into mapping list
  input_list$map$ln_srvsel_devs <- factor(map_srvsel_devs)
  input_list$data$map_ln_srvsel_devs <- array(as.numeric(input_list$map$ln_srvsel_devs), dim = dim(input_list$par$ln_srvsel_devs))

  return(input_list)
}


#' Set up survey selectivity and catchability specifications
#'
#' Configures all survey selectivity and catchability components of the
#' estimation model: continuous and blocked time-varying selectivity,
#' selectivity functional forms, catchability blocks and optional
#' environmental covariate effects, process error and deviation mapping, and
#' selectivity/catchability priors. Delegates parameter mapping to four
#' internal helpers (\code{\link{do_srv_fixed_sel_pars_mapping}},
#' \code{\link{do_srv_q_mapping}}, \code{\link{do_srvsel_pe_pars_mapping}},
#' \code{\link{do_srvsel_devs_mapping}}). Must be called after
#' \code{\link{Setup_Mod_SrvIdx_and_Comps}} and before model compilation.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, \code{$map},
#'   and \code{$verbose} sublists, as returned by upstream setup functions.
#'   \code{$data$Selex_Type} must already be set by
#'   \code{\link{Setup_Mod_Biologicals}}.
#'
#' @section Survey Selectivity Model:
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
#'     \item{\code{"logist1"}}{Logistic: \eqn{a_{50}} and \eqn{k}. 2 params.}
#'     \item{\code{"logist2"}}{Logistic: \eqn{a_{50}} and \eqn{a_{95}}.
#'       2 params.}
#'     \item{\code{"gamma"}}{Dome-shaped gamma: \eqn{a_{max}} and
#'       \eqn{\delta}. 2 params.}
#'     \item{\code{"exponential"}}{Exponential with power parameter. 1 param.}
#'     \item{\code{"dbnrml"}}{Double-normal. 6 params.}
#'   }
#'   No default; must be provided.
#' @param srv_fixed_sel_pars_spec Character vector \code{[n_srv_fleets]}.
#'   Sharing structure for fixed-effect selectivity parameters. See
#'   \code{\link{do_srv_fixed_sel_pars_mapping}} for full option descriptions.
#'   No default; must be provided.
#' @param srv_sel_blocks Character vector defining discrete time blocks for
#'   survey selectivity. Each element follows \code{"Block_k_Year_a-b_Fleet_x"}
#'   or \code{"none_Fleet_x"} (constant selectivity). Use \code{"terminal"} in
#'   place of the end year to extend to the final model year. Parsed into an
#'   internal array \code{[n_regions × n_years × n_srv_fleets]}. Blocked and
#'   continuous time-varying selectivity are mutually exclusive for a given
#'   fleet. Default: \code{"none_Fleet_x"} for each fleet.
#'
#' @section Continuous Time-Varying Selectivity:
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
#' @param cont_tv_srv_sel_penalty Logical. Whether to apply penalties to
#'   continuous time-varying selectivity deviations. Default \code{TRUE}.
#' @param srvsel_pe_pars_spec Character vector \code{[n_srv_fleets]} or
#'   \code{NULL}. Sharing structure for process error hyperparameters. See
#'   \code{\link{do_srvsel_pe_pars_mapping}} for full option descriptions.
#'   Default \code{NULL}.
#' @param srv_sel_devs_spec Character vector \code{[n_srv_fleets]} or
#'   \code{NULL}. Sharing structure for selectivity deviation time series.
#'   See \code{\link{do_srvsel_devs_mapping}} for full option descriptions.
#'   Default \code{NULL}.
#' @param srvsel_devs_shared_ages List of integer vectors defining bin groups
#'   for age/length-sharing of deviations under semi-parametric forms (e.g.,
#'   \code{list(1:5, 6:10, 11:30)}). Required when \code{srv_sel_devs_spec}
#'   includes any \code{"est_shared_a"} variant. Default \code{NULL}.
#' @param corr_opt_semipar Character vector \code{[n_srv_fleets]} or
#'   \code{NULL}. Specifies correlation components to suppress for 3D GMRF or
#'   2D AR1 forms. See \code{\link{do_srvsel_pe_pars_mapping}} for valid
#'   values. Default \code{NULL}.
#'
#' @section Survey Catchability:
#'
#' @param srv_q_blocks Character vector defining discrete time blocks for
#'   survey catchability. Same format as \code{srv_sel_blocks}:
#'   \code{"Block_k_Year_a-b_Fleet_x"} or \code{"none_Fleet_x"}. Parsed into
#'   an array \code{[n_regions × n_years × n_srv_fleets]}. Default:
#'   \code{"none_Fleet_x"} for each fleet.
#' @param srv_q_spec Character vector \code{[n_srv_fleets]} or \code{NULL}.
#'   Sharing structure for catchability. See \code{\link{do_srv_q_mapping}}
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
#' @section Survey Timing and Priors:
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
#'   parameters. Supported names and default dimensions:
#'   \code{ln_srv_fixed_sel_pars}
#'     \code{[n_regions × max_sel_pars × max_sel_blocks × n_sexes × n_srv_fleets]},
#'   \code{ln_srv_q}
#'     \code{[n_regions × max_q_blocks × n_srv_fleets]},
#'   \code{srvsel_pe_pars}
#'     \code{[n_regions × max(max_sel_pars, 4) × n_sexes × n_srv_fleets]},
#'   \code{ln_srvsel_devs}
#'     \code{[n_regions × (n_years + n_proj_yrs_devs) × n_bins × n_sexes × n_srv_fleets]},
#'   \code{srv_q_coeff}
#'     \code{[n_regions × n_srv_fleets × n_covariates]}.
#'
#' @return The input \code{input_list} with selectivity and catchability
#'   configuration stored in \code{$data} (\code{cont_tv_srv_sel},
#'   \code{cont_tv_srv_sel_penalty}, \code{srv_sel_blocks},
#'   \code{srv_sel_model}, \code{srv_q_blocks}, \code{srv_q_prior},
#'   \code{Use_srv_q_prior}, \code{do_srv_q_cov}, \code{srv_q_cov},
#'   \code{Use_srv_selex_prior}, \code{srv_selex_prior}, \code{t_srv});
#'   starting values in \code{$par} for \code{ln_srv_fixed_sel_pars},
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
                                   cont_tv_srv_sel_penalty = TRUE,
                                   srvsel_devs_shared_ages = NULL,
                                   ...
                                   ) {

  messages_list <<- character(0) # string to attach to for printing messages
  starting_values <- list(...)

  # Input Validation --------------------------------------------------------

  # Selectivity
  # Selectivity Type
  if(is.null(input_list$data$Selex_Type)) stop("Selectivity type (age or length-based) has not been specified yet! Make sure to first specify biological inputs with Setup_Mod_Biologicals.")

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
  sel_map <- data.frame(sel = c('logist1', "gamma", "exponential", "logist2", "dbnrml"), num = c(0,1,2,3,4)) # set up values we can map to
  srv_sel_model_arr <- array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_srv_fleets))
  for(i in 1:length(srv_sel_model)) {

    # Extract out survey selectivity components from vector
    tmp_sel_form <- srv_sel_model[i]
    tmp_sel_form_vec <- unlist(strsplit(tmp_sel_form, "_")) # split string
    sel_form <- tmp_sel_form_vec[1] # get selectivity type

    # get fleet index
    tmp_fleet <- if(length(tmp_sel_form_vec) == 3) as.numeric(tmp_sel_form_vec[3]) else as.numeric(tmp_sel_form_vec[5]) # fleet index changes if block is included in character vector
    # get block index
    tmp_block <- if(length(tmp_sel_form_vec) == 5) as.numeric(tmp_sel_form_vec[3]) else NULL

    # validate options
    if(!sel_form %in% c(sel_map$sel)) stop("srv_sel_model is not correctly specified. This needs to be one of these: logist1, gamma, exponential, logist2, dbnrml (the seltypes) and specified as seltype_Fleet_x")
    if(!tmp_fleet %in% c(1:input_list$data$n_srv_fleets)) stop("Invalid fleet specified for srv_sel_model This needs to be specified as seltype_Fleet_x or seltype_Fleet_x_Block_x (if blocks are specified to change for a fleet)")

    # Input options
    if(is.null(tmp_block)) srv_sel_model_arr[,,tmp_fleet] <- sel_map$num[which(sel_map$sel == sel_form)] # same selectivity form across blocks
    else srv_sel_model_arr[,which(srv_sel_blocks_arr[,,tmp_fleet] == tmp_block),tmp_fleet] <- sel_map$num[which(sel_map$sel == sel_form)]
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
  input_list$data$cont_tv_srv_sel_penalty <- cont_tv_srv_sel_penalty
  input_list$data$srv_sel_blocks <- srv_sel_blocks_arr
  input_list$data$srv_sel_model <- srv_sel_model_arr
  input_list$data$srv_q_blocks <- srv_q_blocks_arr
  input_list$data$srv_q_prior <- srv_q_prior
  input_list$data$Use_srv_q_prior <- Use_srv_q_prior
  input_list$data$do_srv_q_cov <- do_srv_q_cov
  input_list$data$srv_q_cov <- srv_q_cov
  input_list$data$Use_srv_selex_prior <- Use_srv_selex_prior
  input_list$data$srv_selex_prior <- srv_selex_prior
  input_list$data$t_srv <- t_srv


  # Populate Parameter List -------------------------------------------------
  # Figure out number of selectivity parameters for a given functional form
  unique_srvsel_vals <- unique(as.vector(input_list$data$srv_sel_model))
  sel_pars_vec <- vector() # create empty vector to populate

  for(i in 1:length(unique_srvsel_vals)) {
    if(unique_srvsel_vals[i] %in% c(2)) sel_pars_vec[i] <- 1 # exponential
    if(unique_srvsel_vals[i] %in% c(0,1,3)) sel_pars_vec[i] <- 2 # logistic or gamma
    if(unique_srvsel_vals[i] == 4) sel_pars_vec[i] <- 6 # double normal
  } # end i loop

  max_srvsel_blks <- max(apply(input_list$data$srv_sel_blocks, c(1,3), FUN = function(x) length(unique(x)))) # figure out maximum number of survey selectivity blocks for a given reigon and fleet
  max_srvsel_pars <- max(sel_pars_vec) # maximum number of selectivity parameters across all forms
  if("ln_srv_fixed_sel_pars" %in% names(starting_values)) input_list$par$ln_srv_fixed_sel_pars <- starting_values$ln_srv_fixed_sel_pars
  else input_list$par$ln_srv_fixed_sel_pars <- array(0, dim = c(input_list$data$n_regions, max_srvsel_pars, max_srvsel_blks, input_list$data$n_sexes, input_list$data$n_srv_fleets))

  # Survey catchability
  max_srvq_blks <- max(apply(input_list$data$srv_q_blocks, c(1,3), FUN = function(x) length(unique(x)))) # figure out maximum number of survey catchability blocks for a given reigon and fleet
  if("ln_srv_q" %in% names(starting_values)) input_list$par$ln_srv_q <- starting_values$ln_srv_q
  else input_list$par$ln_srv_q <- array(0, dim = c(input_list$data$n_regions, max_srvq_blks, input_list$data$n_srv_fleets))

  # Survey selectivity process error parameters
  if("srvsel_pe_pars" %in% names(starting_values)) input_list$par$srvsel_pe_pars <- starting_values$srvsel_pe_pars
  else input_list$par$srvsel_pe_pars <- array(0, dim = c(input_list$data$n_regions, max(max_srvsel_pars, 4), input_list$data$n_sexes, input_list$data$n_srv_fleets)) # dimensioned 4 as the max number of pars for process errors (e.g., sigmas), and then just map off if not using

  # Survey selectivity deviations
  if(input_list$data$Selex_Type == 0) bins <- length(input_list$data$ages) # age based deviations
  if(input_list$data$Selex_Type == 1) bins <- length(input_list$data$lens) # length based deviations
  if("ln_srvsel_devs" %in% names(starting_values)) input_list$par$ln_srvsel_devs <- starting_values$ln_srvsel_devs
  else input_list$par$ln_srvsel_devs <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years) + input_list$data$n_proj_yrs_devs, bins, input_list$data$n_sexes, input_list$data$n_srv_fleets))

  # Survey catchability covariate effects
  if("srv_q_coeff" %in% names(starting_values)) input_list$par$srv_q_coeff <- starting_values$srv_q_coeff
  else input_list$par$srv_q_coeff <- srv_q_coeff # input parameter array

  # Mapping Options ---------------------------------------------------------
  input_list$map$srv_q_coeff <- factor(map_srv_q_coeff) # set up mapping for catchability covariate
  input_list <- do_srv_fixed_sel_pars_mapping(input_list, srv_fixed_sel_pars_spec)
  input_list <- do_srv_q_mapping(input_list, srv_q_spec)
  input_list <- do_srvsel_pe_pars_mapping(input_list, srvsel_pe_pars_spec, corr_opt_semipar)
  input_list <- do_srvsel_devs_mapping(input_list, srv_sel_devs_spec, srvsel_devs_shared_ages)


  # Print Messages ----------------------------------------------------------
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)
}
