# Operating model
#
# Self testing: simulate data from known parameters, fit the estimation model
# back to it, and compare. simulation_data_to_SPoRC is the bridge that turns
# operating model output into the data and par lists fit_model expects.

#' Run a simulation self-test of a fitted RTMB estimation model
#'
#' Validates model performance by: (1) generating \code{n_sims} new datasets
#' from the fitted model parameters using \code{\link{Simulate_Pop_Static}},
#' (2) re-fitting the estimation model to each simulated dataset, and (3)
#' storing user-specified report quantities for comparison against true values.
#' Supports sequential or parallel execution via
#' \code{future}/\code{future.apply}. Likelihood weights from the original fit
#' are propagated into the simulation (e.g., ISS scaled by
#' \code{Wt_FishAgeComps}; \code{ObsSrvIdx_SE} divided by
#' \code{sqrt(Wt_SrvIdx)}); all weights are reset to 1 when re-fitting.
#' Failed replicates are silently stored as \code{NA}.
#'
#' @param data Named list of model data from a fitted RTMB object
#'   (\code{$data}).
#' @param parameters Named list of fitted parameter values (\code{$par} or
#'   equivalent).
#' @param mapping Named list of parameter factor maps (\code{$map}).
#' @param random Character vector of random effect names passed to RTMB.
#' @param rep Named list of report values from the fitted model
#'   (\code{obj$rep}).
#' @param sd_rep \code{sdreport} object from the fitted model, used to
#'   extract optimised parameter values in list format via
#'   \code{get_optim_param_list}.
#' @param n_sims Integer. Number of simulation replicates.
#' @param newton_loops Integer. Number of Newton refinement steps applied
#'   during re-fitting. Default \code{3}.
#' @param do_sdrep Logical. Whether to compute \code{sdreport} for each
#'   fitted replicate. Results stored as \code{$sd_rep} in the output list;
#'   failed \code{sdreport} calls stored as \code{NA}. Default \code{FALSE}.
#' @param do_par Logical. Whether to run replicates in parallel via
#'   \code{future::multisession}. Default \code{FALSE}.
#' @param n_cores Integer. Number of parallel workers. If \code{NULL}
#'   (default), \code{parallel::detectCores() - 1} is used.
#' @param output_path Character string. Path to save the simulated dataset
#'   RDS file. Passed to \code{\link{Simulate_Pop_Static}}. Default
#'   \code{NULL}.
#' @param what Character vector. Names of report elements (keys of
#'   \code{rep}) to extract and store from each replicate. An error is raised
#'   if any name is not found in \code{rep}. Default \code{c("SSB", "Rec")}.
#' @param sim_recruitment Character. How the operating model generates
#'   recruitment. \code{"input"} (default, and the historical behaviour) feeds
#'   the estimated recruitment series in as \code{Rec_input}, so every simulated
#'   replicate carries the same recruitment and \code{rec_model} has no effect on
#'   the data. That conditions away recruitment and tests everything downstream
#'   of it, but it cannot test the stock-recruit relationship itself, because
#'   steepness is then informed only by its penalty. \code{"model"} withholds the
#'   input so recruitment is generated from the fitted curve under
#'   \code{rec_model}, which is what to use when the point of the test is whether
#'   steepness and \code{R0} are recoverable.
#'
#' @return Named list with one element per entry in \code{what}, each an
#'   array with the last dimension indexing simulation replicates (via
#'   \code{simplify2array}). If \code{do_sdrep = TRUE}, an additional element
#'   \code{"sd_rep"} contains a list of \code{sdreport} objects (or \code{NA}
#'   for failed replicates).
#'
#'
#' @export
#' @family Simulation Setup
#'
#' @examples
#' \dontrun{
#' res <- simulation_self_test(
#'   data = obj$data, parameters = obj$par, mapping = obj$map,
#'   random = obj$random, rep = obj$rep, sd_rep = obj$sd_rep,
#'   n_sims = 100, what = c("SSB", "Rec", "Fmort")
#' )
#' str(res$SSB)
#' }
simulation_self_test <- function(data,
                                 parameters,
                                 mapping,
                                 random,
                                 rep,
                                 sd_rep,
                                 n_sims,
                                 newton_loops = 3,
                                 do_sdrep = FALSE,
                                 do_par = FALSE,
                                 n_cores = NULL,
                                 output_path = NULL,
                                 what = c('SSB', 'Rec'),
                                 sim_recruitment = c("input", "model")
                                 ) {

  sim_recruitment <- match.arg(sim_recruitment)

  missing_names <- setdiff(what, names(rep))
  if(length(missing_names) > 0)  stop(paste("The following elements in 'what' are not found in rep:",  paste(missing_names, collapse = ", ")))
  optim_parameters_list <- get_optim_param_list(parameters, mapping, sd_rep, random) # get optimized parameters in original list format

  # Likelihood weights are converted back into simulation standard deviations as
  # sd / sqrt(wt). A weight of zero means the datum is excluded from the fit, not
  # that it carries infinite error, so dividing by sqrt(0) gives Inf and
  # rnorm(1, 0, Inf) gives NaN, which then propagates through every replicate and
  # is only caught as a silent failure at the end. Excluded cells keep their
  # nominal sd instead; they are never fit, so the value does not matter.
  deweight <- function(sd, wt) {
    w <- wt
    w[!is.finite(w) | w <= 0] <- 1
    sd / sqrt(w)
  }

  # Modify any data weights that are NA to 0
  if(any(is.na(data$Wt_Catch))) data$Wt_Catch[is.na(data$Wt_Catch)] <- 0
  if(any(is.na(data$Wt_Catch_pop))) data$Wt_Catch_pop[is.na(data$Wt_Catch_pop)] <- 0
  if(any(is.na(data$Wt_FishAgeComps))) data$Wt_FishAgeComps[is.na(data$Wt_FishAgeComps)] <- 0
  if(any(is.na(data$Wt_FishLenComps))) data$Wt_FishLenComps[is.na(data$Wt_FishLenComps)] <- 0
  if(any(is.na(data$Wt_FishAgeComps_discard))) data$Wt_FishAgeComps_discard[is.na(data$Wt_FishAgeComps_discard)] <- 0
  if(any(is.na(data$Wt_FishLenComps_discard))) data$Wt_FishLenComps_discard[is.na(data$Wt_FishLenComps_discard)] <- 0
  if(any(is.na(data$Wt_FishIdx))) data$Wt_FishIdx[is.na(data$Wt_FishIdx)] <- 0
  if(any(is.na(data$Wt_FishAgeComps_pop))) data$Wt_FishAgeComps_pop[is.na(data$Wt_FishAgeComps_pop)] <- 0
  if(any(is.na(data$Wt_FishLenComps_pop))) data$Wt_FishLenComps_pop[is.na(data$Wt_FishLenComps_pop)] <- 0
  if(any(is.na(data$Wt_FishAgeComps_discard_pop))) data$Wt_FishAgeComps_pop[is.na(data$Wt_FishAgeComps_discard_pop)] <- 0
  if(any(is.na(data$Wt_FishLenComps_discard_pop))) data$Wt_FishLenComps_pop[is.na(data$Wt_FishLenComps_discard_pop)] <- 0
  if(any(is.na(data$Wt_FishIdx_pop))) data$Wt_FishIdx_pop[is.na(data$Wt_FishIdx_pop)] <- 0
  if(any(is.na(data$Wt_SrvAgeComps))) data$Wt_SrvAgeComps[is.na(data$Wt_SrvAgeComps)] <- 0
  if(any(is.na(data$Wt_SrvLenComps))) data$Wt_SrvLenComps[is.na(data$Wt_SrvLenComps)] <- 0
  if(any(is.na(data$Wt_SrvIdx))) data$Wt_SrvIdx[is.na(data$Wt_SrvIdx)] <- 0
  if(any(is.na(data$Wt_SrvAgeComps_pop))) data$Wt_SrvAgeComps_pop[is.na(data$Wt_SrvAgeComps_pop)] <- 0
  if(any(is.na(data$Wt_SrvLenComps_pop))) data$Wt_SrvLenComps_pop[is.na(data$Wt_SrvLenComps_pop)] <- 0
  if(any(is.na(data$Wt_SrvIdx_pop))) data$Wt_SrvIdx_pop[is.na(data$Wt_SrvIdx_pop)] <- 0
  if(any(is.na(data$Wt_Tagging))) data$Wt_Tagging[is.na(data$Wt_Tagging)] <- 0

  # Setup Model Dimensions --------------------------------------------------
  sim_list <- Setup_Sim_Dim(n_sims = n_sims, # number of simulations
                            n_yrs = length(data$years), # number of years
                            n_regions = data$n_regions,  # number of regions
                            n_ages = length(data$ages), # number of ages
                            # Use fishery or survey observed ages depending on what is availiable
                            n_obs_ages = if(any(data$UseFishAgeComps == 1)) {
                              dim(data$ObsFishAgeComps)[4]
                            } else if(any(data$UseFishAgeComps_pop == 1)) {
                              dim(data$ObsFishAgeComps_pop)[5]
                            } else if(any(data$UseSrvAgeComps == 1)) {
                              dim(data$ObsSrvAgeComps)[4]
                            } else if(any(data$UseSrvAgeComps_pop == 1)) {
                              dim(data$ObsSrvAgeComps_pop)[5]
                            } else if(!is.null(data$UseFish_caal) && any(data$UseFish_caal == 1)) {
                              dim(data$ObsFish_caal)[5] # conditional age-at-length carries the observed ages
                            } else if(!is.null(data$UseSrv_caal) && any(data$UseSrv_caal == 1)) {
                              dim(data$ObsSrv_caal)[5]
                            } else {
                              length(data$ages) # no age observations at all, so the containers take the model ages
                            },
                            n_lens = length(data$lens), # number of lengths
                            n_sexes = data$n_sexes, # number of sexes
                            n_fish_fleets = data$n_fish_fleets, # number of fishery fleets
                            n_srv_fleets = data$n_srv_fleets, # number of survey fleets
                            # Seasonal stuff
                            n_seas = data$n_seas,
                            seasdur = data$seasdur,
                            # Population stuff
                            n_pop = data$n_pop,
                            natal_region = data$natal_region
  )

  # Setup Simulation Containers ---------------------------------------------
  sim_list <- Setup_Sim_Containers(sim_list)

  # Setup Fishing Processes -------------------------------------------------

  # Region-specific sigmaC
  ln_sigmaC <- array(NA, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets)) # setup sigmaC container
  # Loop through to populate ln_sigmaC with associated weights
  for(r in 1:sim_list$n_regions) for(f in 1:sim_list$n_fish_fleets) {
    if(!is.vector(data$Wt_Catch)) ln_sigmaC[r,,,f] <- log(deweight(exp(optim_parameters_list$ln_sigmaC[r,,,f]), data$Wt_Catch[r,,,f]))
    else ln_sigmaC[r,,,f] <- log(deweight(exp(optim_parameters_list$ln_sigmaC[r,,,f]), data$Wt_Catch))
  }

  # Population-specific sigmaC
  ln_sigmaC_pop <- array(NA, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets)) # setup sigmaC container
  # Loop through to populate ln_sigmaC with associated weights
  for(p in 1:sim_list$n_pop) for(r in 1:sim_list$n_regions) for(f in 1:sim_list$n_fish_fleets) {
    if(!is.vector(data$Wt_Catch_pop)) ln_sigmaC_pop[p,r,,,f] <- log(deweight(exp(optim_parameters_list$ln_sigmaC_pop[p,r,,,f]), data$Wt_Catch_pop[p,r,,,f]))
    else ln_sigmaC_pop[p,r,,,f] <- log(deweight(exp(optim_parameters_list$ln_sigmaC_pop[p,r,,,f]), data$Wt_Catch_pop))
  }

  # Region-specific sigmaD
  ln_sigmaD <- array(NA, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets)) # setup sigmaD container
  # Loop through to populate ln_sigmaD with associated weights
  for(r in 1:sim_list$n_regions) for(f in 1:sim_list$n_fish_fleets) {
    if(!is.vector(data$Wt_Discard)) ln_sigmaD[r,,,f] <- log(deweight(exp(optim_parameters_list$ln_sigmaD[r,,,f]), data$Wt_Discard[r,,,f]))
    else ln_sigmaD[r,,,f] <- log(deweight(exp(optim_parameters_list$ln_sigmaD[r,,,f]), data$Wt_Discard))
  }

  # Population-specific sigmaD
  ln_sigmaD_pop <- array(NA, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets)) # setup sigmaD container
  # Loop through to populate ln_sigmaD with associated weights
  for(p in 1:sim_list$n_pop) for(r in 1:sim_list$n_regions) for(f in 1:sim_list$n_fish_fleets) {
    if(!is.vector(data$Wt_Discard_pop)) ln_sigmaD_pop[p,r,,,f] <- log(deweight(exp(optim_parameters_list$ln_sigmaD_pop[p,r,,,f]), data$Wt_Discard_pop[p,r,,,f]))
    else ln_sigmaD_pop[p,r,,,f] <- log(deweight(exp(optim_parameters_list$ln_sigmaD_pop[p,r,,,f]), data$Wt_Discard_pop))
  }

  # setup fishery simulation processes
  sim_list <- Setup_Sim_Fishing(sim_list = sim_list,
                                ln_sigmaC = ln_sigmaC,
                                ln_sigmaC_pop = ln_sigmaC_pop,
                                ln_sigmaD = ln_sigmaD,
                                ln_sigmaD_pop = ln_sigmaD_pop,
                                catch_units = data$catch_units,
                                discard_units = data$discard_units,
                                Fmort_input = replicate(n = sim_list$n_sims, rep$Fmort[,1:length(data$years),,,drop = FALSE]),
                                dmr_input = replicate(n = sim_list$n_sims, rep$dmr[,1:length(data$years),,,drop = FALSE]),
                                fish_sel_input = replicate(n = sim_list$n_sims, rep$fish_sel[,,1:length(data$years),,,,,drop = FALSE]),
                                ret_sel_input = replicate(n = sim_list$n_sims, rep$ret_sel[,,1:length(data$years),,,,,drop = FALSE]),
                                fish_q_input = replicate(n = sim_list$n_sims, rep$fish_q[,1:length(data$years),,drop = FALSE]),
                                ObsFishIdx_SE = deweight(if(is.null(rep$FishIdx_SD)) data$ObsFishIdx_SE else rep$FishIdx_SD,
                             data$Wt_FishIdx),
                                ObsFishIdx_pop_SE = if(any(data$UseFishIdx_pop == 1)) {
                                  deweight(data$ObsFishIdx_pop_SE, data$Wt_FishIdx_pop)
                                } else {
                                  array(0.2, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets))
                                },
                                fish_idx_type = data$fish_idx_type,
                                FishIdx_LikeType = if(is.null(data$FishIdx_LikeType)) rep(0, data$n_fish_fleets) else data$FishIdx_LikeType,
                                FishIdx_Cov = data$FishIdx_Cov,
                                UseFishIdx = data$UseFishIdx,
                                init_F_val = rep$init_F,

                                # fishery age composition specifications
                                comp_fishage_like = data$FishAgeComps_LikeType,
                                FishAgeComps_Type = data$FishAgeComps_Type,
                                ISS_FishAgeComps = replicate(sim_list$n_sims, data$ISS_FishAgeComps[,,,,,drop = F] * data$Wt_FishAgeComps),
                                ln_FishAge_theta = optim_parameters_list$ln_FishAge_theta[,,,drop = F],
                                ln_FishAge_theta_agg = optim_parameters_list$ln_FishAge_theta_agg,
                                FishAge_corr_pars_agg = optim_parameters_list$FishAge_corr_pars_agg,
                                FishAge_corr_pars = optim_parameters_list$FishAge_corr_pars[,,,,drop = F],

                                # fishery length composition specifications
                                comp_fishlen_like = data$FishLenComps_LikeType,
                                FishLenComps_Type = data$FishLenComps_Type,
                                ISS_FishLenComps = replicate(sim_list$n_sims, data$ISS_FishLenComps[,,,,,drop = F] * data$Wt_FishLenComps),
                                ln_FishLen_theta = optim_parameters_list$ln_FishLen_theta[,,,drop = F],
                                ln_FishLen_theta_agg = optim_parameters_list$ln_FishLen_theta_agg,
                                FishLen_corr_pars_agg = optim_parameters_list$FishLen_corr_pars_agg,
                                FishLen_corr_pars = optim_parameters_list$FishLen_corr_pars[,,,,drop = F],

                                # population-specific age composition specifications
                                comp_fishage_pop_like = data$pop_FishAgeComps_LikeType,
                                FishAgeComps_pop_Type = data$FishAgeComps_pop_Type,
                                ISS_FishAgeComps_pop = if(any(data$UseFishAgeComps_pop == 1)) {
                                  replicate(sim_list$n_sims, data$ISS_FishAgeComps_pop[,,,,,,drop = F] * data$Wt_FishAgeComps_pop)
                                } else {
                                  array(100, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims))
                                },
                                ln_FishAge_pop_theta = optim_parameters_list$ln_FishAge_pop_theta[,,,,drop = F],
                                ln_FishAge_pop_theta_agg = optim_parameters_list$ln_FishAge_pop_theta_agg,
                                FishAge_pop_corr_pars_agg = optim_parameters_list$FishAge_pop_corr_pars_agg,
                                FishAge_pop_corr_pars = optim_parameters_list$FishAge_pop_corr_pars[,,,,,drop = F],

                                # population-specific length composition specifications
                                comp_fishlen_pop_like = data$FishLenComps_pop_LikeType,
                                FishLenComps_pop_Type = data$FishLenComps_pop_Type,
                                ISS_FishLenComps_pop = if(any(data$UseFishLenComps_pop == 1)) {
                                  replicate(sim_list$n_sims, data$ISS_FishLenComps_pop[,,,,,,drop = F] * data$Wt_FishLenComps_pop)
                                } else {
                                  array(100, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims))
                                },
                                ln_FishLen_pop_theta = optim_parameters_list$ln_FishLen_pop_theta[,,,,drop = F],
                                ln_FishLen_pop_theta_agg = optim_parameters_list$ln_FishLen_pop_theta_agg,
                                FishLen_pop_corr_pars_agg = optim_parameters_list$FishLen_pop_corr_pars_agg,
                                FishLen_pop_corr_pars = optim_parameters_list$FishLen_pop_corr_pars[,,,,,drop = F],

                                # discarded fishery age composition specifications
                                comp_fishage_discard_like = data$FishAgeComps_discard_LikeType,
                                FishAgeComps_discard_Type = data$FishAgeComps_discard_Type,
                                ISS_FishAgeComps_discard = replicate(sim_list$n_sims, data$ISS_FishAgeComps_discard[,,,,,drop = F] * data$Wt_FishAgeComps_discard),
                                ln_FishAge_discard_theta = optim_parameters_list$ln_FishAge_discard_theta[,,,drop = F],
                                ln_FishAge_discard_theta_agg = optim_parameters_list$ln_FishAge_discard_theta_agg,
                                FishAge_discard_corr_pars_agg = optim_parameters_list$FishAge_discard_corr_pars_agg,
                                FishAge_discard_corr_pars = optim_parameters_list$FishAge_discard_corr_pars[,,,,drop = F],

                                # discarded fishery length composition specifications
                                comp_fishlen_discard_like = data$FishLenComps_discard_LikeType,
                                FishLenComps_discard_Type = data$FishLenComps_discard_Type,
                                ISS_FishLenComps_discard = replicate(sim_list$n_sims, data$ISS_FishLenComps_discard[,,,,,drop = F] * data$Wt_FishLenComps_discard),
                                ln_FishLen_discard_theta = optim_parameters_list$ln_FishLen_discard_theta[,,,drop = F],
                                ln_FishLen_discard_theta_agg = optim_parameters_list$ln_FishLen_discard_theta_agg,
                                FishLen_discard_corr_pars_agg = optim_parameters_list$FishLen_discard_corr_pars_agg,
                                FishLen_discard_corr_pars = optim_parameters_list$FishLen_discard_corr_pars[,,,,drop = F],

                                # discarded population-specific age composition specifications
                                comp_fishage_discard_pop_like = data$FishAgeComps_discard_pop_LikeType,
                                FishAgeComps_discard_pop_Type = data$FishAgeComps_discard_pop_Type,
                                ISS_FishAgeComps_discard_pop = if(any(data$UseFishAgeComps_discard_pop == 1)) {
                                  replicate(sim_list$n_sims, data$ISS_FishAgeComps_discard_pop[,,,,,,drop = F] * data$Wt_FishAgeComps_discard_pop)
                                } else {
                                  array(100, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims))
                                },
                                ln_FishAge_discard_pop_theta = optim_parameters_list$ln_FishAge_discard_pop_theta[,,,,drop = F],
                                ln_FishAge_discard_pop_theta_agg = optim_parameters_list$ln_FishAge_discard_pop_theta_agg,
                                FishAge_discard_pop_corr_pars_agg = optim_parameters_list$FishAge_discard_pop_corr_pars_agg,
                                FishAge_discard_pop_corr_pars = optim_parameters_list$FishAge_discard_pop_corr_pars[,,,,,drop = F],

                                # discarded population-specific length composition specifications
                                comp_fishlen_discard_pop_like = data$FishLenComps_discard_pop_LikeType,
                                FishLenComps_discard_pop_Type = data$FishLenComps_discard_pop_Type,
                                ISS_FishLenComps_discard_pop = if(any(data$UseFishLenComps_discard_pop == 1)) {
                                  replicate(sim_list$n_sims, data$ISS_FishLenComps_discard_pop[,,,,,,drop = F] * data$Wt_FishLenComps_discard_pop)
                                } else {
                                  array(100, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims))
                                },
                                ln_FishLen_discard_pop_theta = optim_parameters_list$ln_FishLen_discard_pop_theta[,,,,drop = F],
                                ln_FishLen_discard_pop_theta_agg = optim_parameters_list$ln_FishLen_discard_pop_theta_agg,
                                FishLen_discard_pop_corr_pars_agg = optim_parameters_list$FishLen_discard_pop_corr_pars_agg,
                                FishLen_discard_pop_corr_pars = optim_parameters_list$FishLen_discard_pop_corr_pars[,,,,,drop = F],

                                # conditional age-at-length specifications; absent on models built
                                # before the stream existed, which the defaults leave off
                                comp_fish_caal_like = if(is.null(data$Fish_caal_LikeType)) rep(999, sim_list$n_fish_fleets) else data$Fish_caal_LikeType,
                                Fish_caal_Type = if(is.null(data$Fish_caal_Type)) array(999, dim = c(sim_list$n_yrs, sim_list$n_fish_fleets)) else data$Fish_caal_Type,
                                ISS_Fish_caal = if(!is.null(data$UseFish_caal) && any(data$UseFish_caal == 1)) {
                                  replicate(sim_list$n_sims, data$ISS_Fish_caal[,,,,,,drop = F] * data$Wt_Fish_caal)
                                } else NULL,
                                ln_Fish_caal_theta = optim_parameters_list$ln_Fish_caal_theta,
                                ln_Fish_caal_theta_agg = optim_parameters_list$ln_Fish_caal_theta_agg
  )

  # Setup Survey Processes --------------------------------------------------
  sim_list <- Setup_Sim_Survey(
    sim_list = sim_list,
    srv_sel_input = replicate(n = sim_list$n_sims, rep$srv_sel[,,1:length(data$years),,,,,drop = FALSE]),
    srv_q_input = replicate(n = sim_list$n_sims, rep$srv_q[,1:length(data$years),,drop = FALSE]),
    ObsSrvIdx_SE = deweight(if(is.null(rep$SrvIdx_SD)) data$ObsSrvIdx_SE else rep$SrvIdx_SD, data$Wt_SrvIdx),
    ObsSrvIdx_pop_SE = if(any(data$UseSrvIdx_pop == 1)) {
      deweight(data$ObsSrvIdx_pop_SE, data$Wt_SrvIdx_pop)
    } else {
      array(0.2, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_srv_fleets))
    },
    srv_idx_type = data$srv_idx_type,
    SrvIdx_LikeType = if(is.null(data$SrvIdx_LikeType)) rep(0, data$n_srv_fleets) else data$SrvIdx_LikeType,
    SrvIdx_Cov = data$SrvIdx_Cov,
    UseSrvIdx = data$UseSrvIdx,
    t_srv = data$t_srv,

    # survey age composition specifications
    comp_srvage_like = data$SrvAgeComps_LikeType,
    SrvAgeComps_Type = data$SrvAgeComps_Type,
    ISS_SrvAgeComps = replicate(sim_list$n_sims, data$ISS_SrvAgeComps[,,,,,drop = F] * data$Wt_SrvAgeComps),
    ln_SrvAge_theta = optim_parameters_list$ln_SrvAge_theta[,,,drop = F],
    ln_SrvAge_theta_agg = optim_parameters_list$ln_SrvAge_theta_agg,
    SrvAge_corr_pars_agg = optim_parameters_list$SrvAge_corr_pars_agg,
    SrvAge_corr_pars = optim_parameters_list$SrvAge_corr_pars[,,,,drop = F],

    # survey length composition specifications
    comp_srvlen_like = data$SrvLenComps_LikeType,
    SrvLenComps_Type = data$SrvLenComps_Type,
    ISS_SrvLenComps = replicate(sim_list$n_sims, data$ISS_SrvLenComps[,,,,,drop = F] * data$Wt_SrvLenComps),
    ln_SrvLen_theta = optim_parameters_list$ln_SrvLen_theta[,,,drop = F],
    ln_SrvLen_theta_agg = optim_parameters_list$ln_SrvLen_theta_agg,
    SrvLen_corr_pars_agg = optim_parameters_list$SrvLen_corr_pars_agg,
    SrvLen_corr_pars = optim_parameters_list$SrvLen_corr_pars[,,,,drop = F],

    # population-specific age composition specifications
    comp_srvage_pop_like = data$SrvAgeComps_pop_LikeType,
    SrvAgeComps_pop_Type = data$SrvAgeComps_pop_Type,
    ISS_SrvAgeComps_pop = if(any(data$UseSrvAgeComps_pop == 1)) {
      replicate(sim_list$n_sims, data$ISS_SrvAgeComps_pop[,,,,,,drop = F] * data$Wt_SrvAgeComps_pop)
    } else {
      array(100, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims))
    },
    ln_SrvAge_pop_theta = optim_parameters_list$ln_SrvAge_pop_theta[,,,,drop = F],
    ln_SrvAge_pop_theta_agg = optim_parameters_list$ln_SrvAge_pop_theta_agg,
    SrvAge_pop_corr_pars_agg = optim_parameters_list$SrvAge_pop_corr_pars_agg,
    SrvAge_pop_corr_pars = optim_parameters_list$SrvAge_pop_corr_pars[,,,,,drop = F],

    # population-specific length composition specifications
    comp_srvlen_pop_like = data$SrvLenComps_pop_LikeType,
    SrvLenComps_pop_Type = data$SrvLenComps_pop_Type,
    ISS_SrvLenComps_pop = if(any(data$UseSrvLenComps_pop == 1)) {
      replicate(sim_list$n_sims, data$ISS_SrvLenComps_pop[,,,,,,drop = F] * data$Wt_SrvLenComps_pop)
    } else {
      array(100, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims))
    },
    ln_SrvLen_pop_theta = optim_parameters_list$ln_SrvLen_pop_theta[,,,,drop = F],
    ln_SrvLen_pop_theta_agg = optim_parameters_list$ln_SrvLen_pop_theta_agg,
    SrvLen_pop_corr_pars_agg = optim_parameters_list$SrvLen_pop_corr_pars_agg,
    SrvLen_pop_corr_pars = optim_parameters_list$SrvLen_pop_corr_pars[,,,,,drop = F],

    # conditional age-at-length specifications
    comp_srv_caal_like = if(is.null(data$Srv_caal_LikeType)) rep(999, sim_list$n_srv_fleets) else data$Srv_caal_LikeType,
    Srv_caal_Type = if(is.null(data$Srv_caal_Type)) array(999, dim = c(sim_list$n_yrs, sim_list$n_srv_fleets)) else data$Srv_caal_Type,
    ISS_Srv_caal = if(!is.null(data$UseSrv_caal) && any(data$UseSrv_caal == 1)) {
      replicate(sim_list$n_sims, data$ISS_Srv_caal[,,,,,,drop = F] * data$Wt_Srv_caal)
    } else NULL,
    ln_Srv_caal_theta = optim_parameters_list$ln_Srv_caal_theta,
    ln_Srv_caal_theta_agg = optim_parameters_list$ln_Srv_caal_theta_agg
  )

  # Setup Biological Dynamics -----------------------------------------------
  sim_list <- Setup_Sim_Biologicals(
    sim_list = sim_list, # simualtion list
    natmort_input = replicate(n = sim_list$n_sims, rep$natmort[,,1:length(data$years),,,drop = FALSE]), # natuyral mortality
    # derived by the growth module when present, otherwise the data the model was given
    WAA_input = replicate(n = sim_list$n_sims, (if(is.null(rep$WAA)) data$WAA else rep$WAA)[,,1:length(data$years),,,,drop = FALSE]), # weight at age
    WAA_fish_input = replicate(n = sim_list$n_sims, (if(is.null(rep$WAA_fish)) data$WAA_fish else rep$WAA_fish)[,,1:length(data$years),,,,,drop = FALSE]), # fishery weight at age
    WAA_srv_input = replicate(n = sim_list$n_sims, (if(is.null(rep$WAA_srv)) data$WAA_srv else rep$WAA_srv)[,,1:length(data$years),,,,,drop = FALSE]), # survey weight at age
    MatAA_input = replicate(n = sim_list$n_sims, data$MatAA[,,1:length(data$years),,,,drop = FALSE]), # maturity at age
    AgeingError_input = replicate(n = sim_list$n_sims, data$AgeingError[1:length(data$years),,,drop = FALSE]), # ageing error
    # fleet-specific ageing error, absent from data lists written before it existed,
    # in which case the operating model falls back on the shared matrix
    AgeingError_fish_input = if(is.null(data$AgeingError_fish)) NULL else replicate(n = sim_list$n_sims, data$AgeingError_fish[1:length(data$years),,,,drop = FALSE]),
    AgeingError_srv_input = if(is.null(data$AgeingError_srv)) NULL else replicate(n = sim_list$n_sims, data$AgeingError_srv[1:length(data$years),,,,drop = FALSE]),
    SizeAgeTrans_input = if(data$fit_lengths == 0 || is.null(data$SizeAgeTrans) || all(is.na(data$SizeAgeTrans))) NULL else replicate(n = sim_list$n_sims, data$SizeAgeTrans[,,1:length(data$years),,,,,drop = FALSE]),
    # keys per fleet from the growth module, each at its fleet's own timing
    SizeAgeTrans_fish_input = if(is.null(rep$SizeAgeTrans_fish)) NULL else replicate(n = sim_list$n_sims, rep$SizeAgeTrans_fish[,,1:length(data$years),,,,,,drop = FALSE]),
    SizeAgeTrans_srv_input = if(is.null(rep$SizeAgeTrans_srv)) NULL else replicate(n = sim_list$n_sims, rep$SizeAgeTrans_srv[,,1:length(data$years),,,,,,drop = FALSE]) # size age transition matrix, derived by the growth module when present
  )

  # Movement
  sim_list$Movement <- replicate(n = sim_list$n_sims, rep$Movement[,,,1:length(data$years),,,,drop = FALSE])
  sim_list$sgl_seas_spawning_movement <- replicate(n = sim_list$n_sims, rep$sgl_seas_spawning_movement[,,,1:length(data$years),,,drop = FALSE])
  # Movement / mortality sequencing; absent for models built before this option existed
  sim_list$move_timing <- if(is.null(data$move_timing)) 0 else data$move_timing
  # The instantaneous rate matrix only exists for an estimated CTMC, and is only needed
  # for continuous movement
  sim_list$Mrate <- if(sim_list$move_timing == 2)
    replicate(n = sim_list$n_sims, rep$Mrate[,,,1:length(data$years),,,,drop = FALSE]) else NULL

  # Setup Recruitment Processes ---------------------------------------------
  sim_list <- Setup_Sim_Rec(
    sim_list = sim_list,
    spawn_seas = data$spawn_seas, # spawning season
    do_recruits_move = data$do_recruits_move, # whether recruits move
    t_spawn = data$t_spawn, # spawn timing
    init_age_strc = data$init_age_strc, # initilaizing age structure
    h_input = replicate(n = sim_list$n_sims, array(rep$h_trans, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs))), # steepness
    R0_input = {
      tmp = array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_sims))
      for(p in 1:sim_list$n_pop) for(r in 1:sim_list$n_regions) tmp[p,r,,] = rep$R0[p] * rep$rec_region_prop[p,r]
      tmp
    },
    rinit_input = {
      tmp = array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_sims))
      for(p in 1:sim_list$n_pop) for(r in 1:sim_list$n_regions) tmp[p,r,] = rep$rinit[p] * rep$rec_region_prop[p,r]
      tmp
    },
    use_rinit = data$use_rinit,
    sexratio_input = replicate(n = sim_list$n_sims, expr = rep$sexratio[,,1:length(data$years),,drop = FALSE]), # sex ratio
    # Rescaling by the recruitment weight is only an identity for a single scalar,
    # and a deviation-specific weight has no equivalent sigma. Both recruitment
    # and the initial age deviations are supplied directly below, so ln_sigmaR is
    # never drawn from here and the unscaled value is passed through instead.
    ln_sigmaR = if(length(data$Wt_Rec) == 1) optim_parameters_list$ln_sigmaR / sqrt(data$Wt_Rec) else optim_parameters_list$ln_sigmaR, # ln_sigmaR
    # Supplying Rec_input for every year makes use_rec_input TRUE throughout
    # sim_population(), so recruitment_opt is never consulted and rec_model has no
    # effect on the simulated data. That is the right conditioning for testing
    # everything downstream of recruitment, but it cannot test the stock-recruit
    # relationship itself: steepness is then informed only by its penalty.
    # sim_recruitment = "model" withholds the input so the curve generates
    # recruitment and the self-test has to recover it. Setup_Sim_Rec only assigns
    # sim_list$Rec_input when it is non-NULL, so NULL is what turns this off.
    Rec_input = if(sim_recruitment == "model") NULL else replicate(n = sim_list$n_sims, expr = rep$Rec[,,1:length(data$years),drop = FALSE]), # recruitment time series
    ln_InitDevs_input = replicate(sim_list$n_sims, optim_parameters_list$ln_InitDevs),  # init devs
    stray_rate_input = replicate(sim_list$n_sims, data$stray_rate[,1:length(data$years), drop = FALSE]),
    rec_seas_prop_input = array(
      rep(rep$rec_seas_prop, times = sim_list$n_sims),
      dim = c(data$n_pop, data$n_seas, sim_list$n_sims)), # seasonal recruitment apportionment

    # Not needed; already specified in Rec_input and ln_InitDevs_input
    recruitment_opt = data$rec_model,
    rec_dd = data$rec_dd,
    init_dd = data$rec_dd,
    rec_lag = data$rec_lag
  )

  # Setup Tagging -----------------------------------------------------------
  if(!is.na(sum(data$conv_tagged_fish))) n_tags_rel_input <- apply(data$conv_tagged_fish, 1, sum) else n_tags_rel_input <- NA
  if(exists("conv_tag_release_indicator", data)) conv_tag_release_indicator <- data$conv_tag_release_indicator  else conv_tag_release_indicator <- NA
  conv_tag_fish_reporting_input <- if(!is.null(rep$conv_tag_fish_reporting)) replicate(n = sim_list$n_sims, rep$conv_tag_fish_reporting) else NULL

  sim_list <- Setup_Sim_Tagging(
    sim_list = sim_list, # simulation list
    conv_tag_max_liberty = data$conv_tag_max_liberty, # maximum tag liberty
    conv_tag_t_tagging = data$conv_tag_t_tagging, # time of tagging
    n_tags_rel_input = n_tags_rel_input * data$Wt_Tagging,  # number of tags to release per event
    conv_tag_release_indicator = conv_tag_release_indicator,  # tag release indicator
    ln_init_conv_tag_mort = optim_parameters_list$ln_init_conv_tag_mort,  # inital tagging mortality
    ln_conv_tag_shed = optim_parameters_list$ln_conv_tag_shed, # chronic tag shedding
    conv_tag_fish_reporting_input = conv_tag_fish_reporting_input, # tag reporting rates
    use_conv_fish_tagging = data$use_conv_fish_tagging, # whether or not tagging is used / simulated
    conv_fish_tag_like = data$conv_fish_tag_like, # tag likelihood
    ln_conv_fish_tag_theta = parameters$ln_conv_fish_tag_theta, # tag overdispersion
    conv_tag_release_platform = data$conv_tag_release_platform,  # tag release platform
    conv_fish_tag_attr = data$conv_fish_tag_attr # tag attributes
  )


  # Run Simulation ----------------------------------------------------------

  # storage
  store_res_list <- vector("list", length(what) + 1) # get list
  names(store_res_list) <- c(what, "sd_rep") # name list
  for(j in 1:length(what)) store_res_list[[j]] <- vector("list", n_sims) # stick in n_sims lists into storage

  sim_obj <- Simulate_Pop_Static(sim_list = sim_list, output_path = output_path) # get simulated datasets

  if(do_par == FALSE) {

    for(i in 1:n_sims) {

      tryCatch({

        # set up data stuff
        tmp_data <- data
        tmp_pars <- parameters
        tmp_data$ObsFishIdx <- array(sim_obj$ObsFishIdx[,,,,i], dim = dim(tmp_data$ObsFishIdx))
        tmp_data$ObsSrvIdx <- array(sim_obj$ObsSrvIdx[,,,,i], dim = dim(tmp_data$ObsSrvIdx))
        tmp_data$ObsCatch <- array(sim_obj$ObsCatch[,,,,i], dim = dim(tmp_data$ObsCatch))
        tmp_data$ObsFishAgeComps <- array(sim_obj$ObsFishAgeComps[,,,,,,i], dim = dim(tmp_data$ObsFishAgeComps))
        tmp_data$ObsSrvAgeComps  <- array(sim_obj$ObsSrvAgeComps[,,,,,,i], dim = dim(tmp_data$ObsSrvAgeComps))
        if(tmp_data$fit_lengths != 0) {
          tmp_data$ObsFishLenComps <- array(sim_obj$ObsFishLenComps[,,,,,,i], dim = dim(tmp_data$ObsFishLenComps))
          tmp_data$ObsSrvLenComps  <- array(sim_obj$ObsSrvLenComps[,,,,,,i], dim = dim(tmp_data$ObsSrvLenComps))
        }

        # population-specific observations
        if(any(tmp_data$UseFishIdx_pop == 1)) {
          tmp_data$ObsFishIdx_pop <- array(sim_obj$ObsFishIdx_pop[,,,,,i], dim = dim(tmp_data$ObsFishIdx_pop))
        }
        if(any(tmp_data$UseSrvIdx_pop == 1)) {
          tmp_data$ObsSrvIdx_pop <- array(sim_obj$ObsSrvIdx_pop[,,,,,i], dim = dim(tmp_data$ObsSrvIdx_pop))
        }
        if(any(tmp_data$UseFishAgeComps_pop == 1)) {
          tmp_data$ObsFishAgeComps_pop <- array(sim_obj$ObsFishAgeComps_pop[,,,,,,,i], dim = dim(tmp_data$ObsFishAgeComps_pop))
        }
        if(any(tmp_data$UseSrvAgeComps_pop == 1)) {
          tmp_data$ObsSrvAgeComps_pop <- array(sim_obj$ObsSrvAgeComps_pop[,,,,,,,i], dim = dim(tmp_data$ObsSrvAgeComps_pop))
        }
        if(tmp_data$fit_lengths != 0) {
          if(any(tmp_data$UseFishLenComps_pop == 1)) {
            tmp_data$ObsFishLenComps_pop <- array(sim_obj$ObsFishLenComps_pop[,,,,,,,i], dim = dim(tmp_data$ObsFishLenComps_pop))
          }
          if(any(tmp_data$UseSrvLenComps_pop == 1)) {
            tmp_data$ObsSrvLenComps_pop <- array(sim_obj$ObsSrvLenComps_pop[,,,,,,,i], dim = dim(tmp_data$ObsSrvLenComps_pop))
          }
        }
        if(any(tmp_data$UseCatch_pop == 1)) {
          tmp_data$ObsCatch_pop <- array(sim_obj$ObsCatch_pop[,,,,,i], dim = dim(tmp_data$ObsCatch_pop))
        }

        # set up discarding stuff
        if(any(tmp_data$UseDiscard == 1)) tmp_data$ObsDiscard <- array(sim_obj$ObsDiscard[,,,,i], dim = dim(tmp_data$ObsDiscard))
        if(any(tmp_data$UseDiscard_pop == 1)) tmp_data$ObsDiscard_pop <- array(sim_obj$ObsDiscard_pop[,,,,i], dim = dim(tmp_data$ObsDiscard_pop))
        if(any(tmp_data$UseFishAgeComps_discard == 1)) tmp_data$ObsFishAgeComps_discard <- array(sim_obj$ObsFishAgeComps_discard[,,,,,,i], dim = dim(tmp_data$ObsFishAgeComps_discard))
        if(tmp_data$fit_lengths != 0) if(any(tmp_data$UseFishLenComps_discard == 1)) tmp_data$ObsFishLenComps_discard <- array(sim_obj$ObsFishLenComps_discard[,,,,,,i], dim = dim(tmp_data$ObsFishLenComps_discard))
        if(any(tmp_data$UseFishAgeComps_discard_pop == 1)) tmp_data$ObsFishAgeComps_discard_pop <- array(sim_obj$ObsFishAgeComps_discard_pop[,,,,,,i], dim = dim(tmp_data$ObsFishAgeComps_discard_pop))
        if(tmp_data$fit_lengths != 0) if(any(tmp_data$UseFishLenComps_discard_pop == 1)) tmp_data$ObsFishLenComps_discard_pop <- array(sim_obj$ObsFishLenComps_discard_pop[,,,,,,i], dim = dim(tmp_data$ObsFishLenComps_discard_pop))

        # setup tagging data stuff if tagging is done
        if(any(tmp_data$use_conv_fish_tagging == 1)) {
          tmp_data$conv_tagged_fish <- array(sim_obj$conv_tagged_fish[,,,,i], dim = dim(tmp_data$conv_tagged_fish))
          tmp_data$obs_conv_tag_fish_recap <- array(sim_obj$obs_conv_tag_fish_recap[,,,,,,,,i], dim = dim(tmp_data$obs_conv_tag_fish_recap))
          tmp_data$conv_tag_release_indicator <- sim_obj$conv_tag_release_indicator
        }

        # reset weights
        tmp_data$Wt_Rec[] <- 1
        tmp_data$Wt_D <- 1
        tmp_data$Wt_Tagging <- 1
        tmp_data$Wt_Catch[] <- 1
        tmp_data$Wt_Discard[] <- 1
        tmp_data$Wt_FishAgeComps[] <- 1
        tmp_data$Wt_FishAgeComps_discard[] <- 1
        tmp_data$Wt_FishIdx[] <- 1
        tmp_data$Wt_FishLenComps[] <- 1
        tmp_data$Wt_FishLenComps_discard[] <- 1
        tmp_data$Wt_SrvAgeComps[] <- 1
        tmp_data$Wt_SrvIdx[] <- 1
        tmp_data$Wt_SrvLenComps[] <- 1
        tmp_data$Wt_Catch_pop[] <- 1
        tmp_data$Wt_Discard_pop[] <- 1
        tmp_data$Wt_FishIdx_pop[] <- 1
        tmp_data$Wt_SrvIdx_pop[] <- 1
        tmp_data$Wt_FishAgeComps_pop[] <- 1
        tmp_data$Wt_FishAgeComps_discard_pop[] <- 1
        tmp_data$Wt_SrvAgeComps_pop[] <- 1
        tmp_data$Wt_FishLenComps_pop[] <- 1
        tmp_data$Wt_FishLenComps_discard_pop[] <- 1
        tmp_data$Wt_SrvLenComps_pop[] <- 1

        # input simulated uncertainty
        tmp_pars$ln_sigmaC[] <- sim_list$ln_sigmaC
        tmp_pars$ln_sigmaC_pop[] <- sim_list$ln_sigmaC_pop
        tmp_pars$ln_sigmaD[] <- sim_list$ln_sigmaD
        tmp_pars$ln_sigmaD_pop[] <- sim_list$ln_sigmaD_pop
        tmp_data$ObsFishIdx_SE[] <- sim_list$ObsFishIdx_SE
        tmp_data$ObsSrvIdx_SE[] <- sim_list$ObsSrvIdx_SE
        if(!is.null(tmp_data$ObsCatchAA) && !is.null(sim_list$ObsCatchAA))
          tmp_data$ObsCatchAA[] <- sim_list$ObsCatchAA[,,,,,i]
        if(!is.null(tmp_data$ObsSrvIdxAA) && !is.null(sim_list$ObsSrvIdxAA))
          tmp_data$ObsSrvIdxAA[] <- sim_list$ObsSrvIdxAA[,,,,,i]
        if(!is.null(tmp_pars$ln_sigmaCAA)) tmp_pars$ln_sigmaCAA[] <- parameters$ln_sigmaCAA
        if(!is.null(tmp_pars$ln_sigmaSrvIdxAA)) tmp_pars$ln_sigmaSrvIdxAA[] <- parameters$ln_sigmaSrvIdxAA
        if(!is.null(tmp_pars$ln_sigmaFishIdx)) tmp_pars$ln_sigmaFishIdx[] <- parameters$ln_sigmaFishIdx
        if(!is.null(tmp_pars$ln_sigmaSrvIdx)) tmp_pars$ln_sigmaSrvIdx[] <- parameters$ln_sigmaSrvIdx
        tmp_data$ISS_FishAgeComps[] <- sim_list$ISS_FishAgeComps[,,,,,i]
        tmp_data$ISS_FishLenComps[] <- sim_list$ISS_FishLenComps[,,,,,i]
        if(any(tmp_data$UseFishAgeComps_discard == 1)) tmp_data$ISS_FishAgeComps_discard[] <- sim_list$ISS_FishAgeComps_discard[,,,,,i]
        if(any(tmp_data$UseFishLenComps_discard == 1)) tmp_data$ISS_FishLenComps_discard[] <- sim_list$ISS_FishLenComps_discard[,,,,,i]
        tmp_data$ISS_SrvAgeComps[] <- sim_list$ISS_SrvAgeComps[,,,,,i]
        tmp_data$ISS_SrvLenComps[] <- sim_list$ISS_SrvLenComps[,,,,,i]
        if(any(tmp_data$UseFishIdx_pop == 1)) tmp_data$ObsFishIdx_pop_SE[] <- sim_list$ObsFishIdx_pop_SE
        if(any(tmp_data$UseSrvIdx_pop == 1)) tmp_data$ObsSrvIdx_pop_SE[] <- sim_list$ObsSrvIdx_pop_SE
        if(any(tmp_data$UseFishAgeComps_pop == 1)) tmp_data$ISS_FishAgeComps_pop[] <- sim_list$ISS_FishAgeComps_pop[,,,,,,i]
        if(any(tmp_data$UseFishLenComps_pop == 1)) tmp_data$ISS_FishLenComps_pop[] <- sim_list$ISS_FishLenComps_pop[,,,,,,i]
        if(any(tmp_data$UseSrvAgeComps_pop == 1)) tmp_data$ISS_SrvAgeComps_pop[] <- sim_list$ISS_SrvAgeComps_pop[,,,,,,i]
        if(any(tmp_data$UseSrvLenComps_pop == 1)) tmp_data$ISS_SrvLenComps_pop[] <- sim_list$ISS_SrvLenComps_pop[,,,,,,i]
        if(any(tmp_data$UseFishAgeComps_discard_pop == 1)) tmp_data$ISS_FishAgeComps_discard_pop[] <- sim_list$ISS_FishAgeComps_discard_pop[,,,,,i]
        if(any(tmp_data$UseFishLenComps_discard_pop == 1)) tmp_data$ISS_FishLenComps_discard_pop[] <- sim_list$ISS_FishLenComps_discard_pop[,,,,,i]

        # conditional age-at-length observations
        if(!is.null(tmp_data$UseFish_caal) && any(tmp_data$UseFish_caal == 1)) {
          tmp_data$ObsFish_caal <- array(sim_obj$ObsFish_caal[,,,,,,,i], dim = dim(tmp_data$ObsFish_caal))
          tmp_data$ISS_Fish_caal[] <- sim_list$ISS_Fish_caal[,,,,,,i]
          tmp_data$Wt_Fish_caal[] <- 1
        }
        if(!is.null(tmp_data$UseSrv_caal) && any(tmp_data$UseSrv_caal == 1)) {
          tmp_data$ObsSrv_caal <- array(sim_obj$ObsSrv_caal[,,,,,,,i], dim = dim(tmp_data$ObsSrv_caal))
          tmp_data$ISS_Srv_caal[] <- sim_list$ISS_Srv_caal[,,,,,,i]
          tmp_data$Wt_Srv_caal[] <- 1
        }

        # This replicate's observations are not the ones setup reconciled against,
        # so a bin restriction may have emptied a block that was full before
        tmp_data <- resync_fitted_blocks(tmp_data)

        # Fit model
        obj <- fit_model(
          data = tmp_data,
          parameters = tmp_pars,
          mapping = mapping,
          random = random,
          newton_loops = newton_loops,
          silent = TRUE
        )

        # Populate results into store list
        for(j in 1:length(what)) store_res_list[[j]][[i]] <- obj$rep[[what[j]]]

        if(do_sdrep == TRUE) {
          tryCatch({
            obj$sd_rep <- RTMB::sdreport(obj)
            store_res_list[[length(what) + 1]][[i]] <- obj$sd_rep # input sd report
          }, error = function(e) {
            store_res_list[[length(what) + 1]][[i]] <- NA
          })
        }

      }, error = function(e) {
        # Skip failed simulations, saying why
        warning(sprintf("simulation %d failed: %s", i, conditionMessage(e)), call. = FALSE)
        for(j in 1:length(what)) store_res_list[[j]][[i]] <- NA
        if(do_sdrep == TRUE) store_res_list[[length(what) + 1]][[i]] <- NA
      })

    } # end i loop

    # Convert result lists to array
    for(j in 1:length(what)) store_res_list[[j]] <- simplify2array(store_res_list[[j]])

  } # not doing parallelization

  if(do_par == TRUE) {

    # initialize cores
    if(is.null(n_cores)) n_cores <- parallel::detectCores() - 1
    options(future.globals.maxSize = 5e3 * 1024^2)  # increase parrlalel size
    future::plan(future::multisession, workers = n_cores) # set up cores
    progressr::with_progress({
      p <- progressr::progressor(along = 1:n_sims) # progress bar

      sim_results <- future.apply::future_lapply(1:n_sims, function(i) {

        tryCatch({

          # set up data stuff
          tmp_data <- data
          tmp_pars <- parameters
          tmp_data$ObsFishIdx <- array(sim_obj$ObsFishIdx[,,,,i], dim = dim(tmp_data$ObsFishIdx))
          tmp_data$ObsSrvIdx <- array(sim_obj$ObsSrvIdx[,,,,i], dim = dim(tmp_data$ObsSrvIdx))
          tmp_data$ObsCatch <- array(sim_obj$ObsCatch[,,,,i], dim = dim(tmp_data$ObsCatch))
          tmp_data$ObsFishAgeComps <- array(sim_obj$ObsFishAgeComps[,,,,,,i], dim = dim(tmp_data$ObsFishAgeComps))
          tmp_data$ObsSrvAgeComps  <- array(sim_obj$ObsSrvAgeComps[,,,,,,i], dim = dim(tmp_data$ObsSrvAgeComps))
          if(tmp_data$fit_lengths != 0) {
            tmp_data$ObsFishLenComps <- array(sim_obj$ObsFishLenComps[,,,,,,i], dim = dim(tmp_data$ObsFishLenComps))
            tmp_data$ObsSrvLenComps  <- array(sim_obj$ObsSrvLenComps[,,,,,,i], dim = dim(tmp_data$ObsSrvLenComps))
          }

          # population-specific observations
          if(any(tmp_data$UseFishIdx_pop == 1)) {
            tmp_data$ObsFishIdx_pop <- array(sim_obj$ObsFishIdx_pop[,,,,,i], dim = dim(tmp_data$ObsFishIdx_pop))
          }
          if(any(tmp_data$UseSrvIdx_pop == 1)) {
            tmp_data$ObsSrvIdx_pop <- array(sim_obj$ObsSrvIdx_pop[,,,,,i], dim = dim(tmp_data$ObsSrvIdx_pop))
          }
          if(any(tmp_data$UseFishAgeComps_pop == 1)) {
            tmp_data$ObsFishAgeComps_pop <- array(sim_obj$ObsFishAgeComps_pop[,,,,,,,i], dim = dim(tmp_data$ObsFishAgeComps_pop))
          }
          if(any(tmp_data$UseSrvAgeComps_pop == 1)) {
            tmp_data$ObsSrvAgeComps_pop <- array(sim_obj$ObsSrvAgeComps_pop[,,,,,,,i], dim = dim(tmp_data$ObsSrvAgeComps_pop))
          }
          if(tmp_data$fit_lengths != 0) {
            if(any(tmp_data$UseFishLenComps_pop == 1)) {
              tmp_data$ObsFishLenComps_pop <- array(sim_obj$ObsFishLenComps_pop[,,,,,,,i], dim = dim(tmp_data$ObsFishLenComps_pop))
            }
            if(any(tmp_data$UseSrvLenComps_pop == 1)) {
              tmp_data$ObsSrvLenComps_pop <- array(sim_obj$ObsSrvLenComps_pop[,,,,,,,i], dim = dim(tmp_data$ObsSrvLenComps_pop))
            }
          }
          if(any(tmp_data$UseCatch_pop == 1)) {
            tmp_data$ObsCatch_pop <- array(sim_obj$ObsCatch_pop[,,,,,i], dim = dim(tmp_data$ObsCatch_pop))
          }

          # set up discarding stuff
          if(any(tmp_data$UseDiscard == 1)) tmp_data$ObsDiscard <- array(sim_obj$ObsDiscard[,,,,i], dim = dim(tmp_data$ObsDiscard))
          if(any(tmp_data$UseDiscard_pop == 1)) tmp_data$ObsDiscard_pop <- array(sim_obj$ObsDiscard_pop[,,,,i], dim = dim(tmp_data$ObsDiscard_pop))
          if(any(tmp_data$UseFishAgeComps_discard == 1)) tmp_data$ObsFishAgeComps_discard <- array(sim_obj$ObsFishAgeComps_discard[,,,,,,i], dim = dim(tmp_data$ObsFishAgeComps_discard))
          if(tmp_data$fit_lengths != 0) if(any(tmp_data$UseFishLenComps_discard == 1)) tmp_data$ObsFishLenComps_discard <- array(sim_obj$ObsFishLenComps_discard[,,,,,,i], dim = dim(tmp_data$ObsFishLenComps_discard))
          if(any(tmp_data$UseFishAgeComps_discard_pop == 1)) tmp_data$ObsFishAgeComps_discard_pop <- array(sim_obj$ObsFishAgeComps_discard_pop[,,,,,,i], dim = dim(tmp_data$ObsFishAgeComps_discard_pop))
          if(tmp_data$fit_lengths != 0) if(any(tmp_data$UseFishLenComps_discard_pop == 1)) tmp_data$ObsFishLenComps_discard_pop <- array(sim_obj$ObsFishLenComps_discard_pop[,,,,,,i], dim = dim(tmp_data$ObsFishLenComps_discard_pop))

          # setup tagging data stuff if tagging is done
          if(any(tmp_data$use_conv_fish_tagging == 1)) {
            tmp_data$conv_tagged_fish <- array(sim_obj$conv_tagged_fish[,,,,i], dim = dim(tmp_data$conv_tagged_fish))
            tmp_data$obs_conv_tag_fish_recap <- array(sim_obj$obs_conv_tag_fish_recap[,,,,,,,,i], dim = dim(tmp_data$obs_conv_tag_fish_recap))
            tmp_data$conv_tag_release_indicator <- sim_obj$conv_tag_release_indicator
          }

          # reset weights
          tmp_data$Wt_Rec[] <- 1
          tmp_data$Wt_D <- 1
          tmp_data$Wt_Tagging <- 1
          tmp_data$Wt_Catch[] <- 1
          tmp_data$Wt_Discard[] <- 1
          tmp_data$Wt_FishAgeComps[] <- 1
          tmp_data$Wt_FishAgeComps_discard[] <- 1
          tmp_data$Wt_FishIdx[] <- 1
          tmp_data$Wt_FishLenComps[] <- 1
          tmp_data$Wt_FishLenComps_discard[] <- 1
          tmp_data$Wt_SrvAgeComps[] <- 1
          tmp_data$Wt_SrvIdx[] <- 1
          tmp_data$Wt_SrvLenComps[] <- 1
          tmp_data$Wt_Catch_pop[] <- 1
          tmp_data$Wt_Discard_pop[] <- 1
          tmp_data$Wt_FishIdx_pop[] <- 1
          tmp_data$Wt_SrvIdx_pop[] <- 1
          tmp_data$Wt_FishAgeComps_pop[] <- 1
          tmp_data$Wt_FishAgeComps_discard_pop[] <- 1
          tmp_data$Wt_SrvAgeComps_pop[] <- 1
          tmp_data$Wt_FishLenComps_pop[] <- 1
          tmp_data$Wt_FishLenComps_discard_pop[] <- 1
          tmp_data$Wt_SrvLenComps_pop[] <- 1

          # input simulated uncertainty
          tmp_pars$ln_sigmaC[] <- sim_list$ln_sigmaC
          tmp_pars$ln_sigmaC_pop[] <- sim_list$ln_sigmaC_pop
          tmp_pars$ln_sigmaD[] <- sim_list$ln_sigmaD
          tmp_pars$ln_sigmaD_pop[] <- sim_list$ln_sigmaD_pop
          tmp_data$ObsFishIdx_SE[] <- sim_list$ObsFishIdx_SE
          tmp_data$ObsSrvIdx_SE[] <- sim_list$ObsSrvIdx_SE
          if(!is.null(tmp_data$ObsCatchAA) && !is.null(sim_list$ObsCatchAA))
            tmp_data$ObsCatchAA[] <- sim_list$ObsCatchAA[,,,,,i]
          if(!is.null(tmp_data$ObsSrvIdxAA) && !is.null(sim_list$ObsSrvIdxAA))
            tmp_data$ObsSrvIdxAA[] <- sim_list$ObsSrvIdxAA[,,,,,i]
          if(!is.null(tmp_pars$ln_sigmaCAA)) tmp_pars$ln_sigmaCAA[] <- parameters$ln_sigmaCAA
          if(!is.null(tmp_pars$ln_sigmaSrvIdxAA)) tmp_pars$ln_sigmaSrvIdxAA[] <- parameters$ln_sigmaSrvIdxAA
          if(!is.null(tmp_pars$ln_sigmaFishIdx)) tmp_pars$ln_sigmaFishIdx[] <- parameters$ln_sigmaFishIdx
          if(!is.null(tmp_pars$ln_sigmaSrvIdx)) tmp_pars$ln_sigmaSrvIdx[] <- parameters$ln_sigmaSrvIdx
          tmp_data$ISS_FishAgeComps[] <- sim_list$ISS_FishAgeComps[,,,,,i]
          tmp_data$ISS_FishLenComps[] <- sim_list$ISS_FishLenComps[,,,,,i]
          if(any(tmp_data$UseFishAgeComps_discard == 1)) tmp_data$ISS_FishAgeComps_discard[] <- sim_list$ISS_FishAgeComps_discard[,,,,,i]
          if(any(tmp_data$UseFishLenComps_discard == 1)) tmp_data$ISS_FishLenComps_discard[] <- sim_list$ISS_FishLenComps_discard[,,,,,i]
          tmp_data$ISS_SrvAgeComps[] <- sim_list$ISS_SrvAgeComps[,,,,,i]
          tmp_data$ISS_SrvLenComps[] <- sim_list$ISS_SrvLenComps[,,,,,i]
          if(any(tmp_data$UseFishIdx_pop == 1)) tmp_data$ObsFishIdx_pop_SE[] <- sim_list$ObsFishIdx_pop_SE
          if(any(tmp_data$UseSrvIdx_pop == 1)) tmp_data$ObsSrvIdx_pop_SE[] <- sim_list$ObsSrvIdx_pop_SE
          if(any(tmp_data$UseFishAgeComps_pop == 1)) tmp_data$ISS_FishAgeComps_pop[] <- sim_list$ISS_FishAgeComps_pop[,,,,,,i]
          if(any(tmp_data$UseFishLenComps_pop == 1)) tmp_data$ISS_FishLenComps_pop[] <- sim_list$ISS_FishLenComps_pop[,,,,,,i]
          if(any(tmp_data$UseSrvAgeComps_pop == 1)) tmp_data$ISS_SrvAgeComps_pop[] <- sim_list$ISS_SrvAgeComps_pop[,,,,,,i]
          if(any(tmp_data$UseSrvLenComps_pop == 1)) tmp_data$ISS_SrvLenComps_pop[] <- sim_list$ISS_SrvLenComps_pop[,,,,,,i]
          if(any(tmp_data$UseFishAgeComps_discard_pop == 1)) tmp_data$ISS_FishAgeComps_discard_pop[] <- sim_list$ISS_FishAgeComps_discard_pop[,,,,,i]
          if(any(tmp_data$UseFishLenComps_discard_pop == 1)) tmp_data$ISS_FishLenComps_discard_pop[] <- sim_list$ISS_FishLenComps_discard_pop[,,,,,i]

          # conditional age-at-length observations
          if(!is.null(tmp_data$UseFish_caal) && any(tmp_data$UseFish_caal == 1)) {
            tmp_data$ObsFish_caal <- array(sim_obj$ObsFish_caal[,,,,,,,i], dim = dim(tmp_data$ObsFish_caal))
            tmp_data$ISS_Fish_caal[] <- sim_list$ISS_Fish_caal[,,,,,,i]
            tmp_data$Wt_Fish_caal[] <- 1
          }
          if(!is.null(tmp_data$UseSrv_caal) && any(tmp_data$UseSrv_caal == 1)) {
            tmp_data$ObsSrv_caal <- array(sim_obj$ObsSrv_caal[,,,,,,,i], dim = dim(tmp_data$ObsSrv_caal))
            tmp_data$ISS_Srv_caal[] <- sim_list$ISS_Srv_caal[,,,,,,i]
            tmp_data$Wt_Srv_caal[] <- 1
          }


          # see the note at the single-model path above
          tmp_data <- resync_fitted_blocks(tmp_data)

          # Fit model
          obj <- fit_model(
            data = tmp_data,
            parameters = tmp_pars,
            mapping = mapping,
            random = random,
            newton_loops = newton_loops,
            silent = TRUE
          )

          # Extract what we need and return
          result <- list()
          for(j in 1:length(what)) result[[what[j]]] <- obj$rep[[what[j]]]

          if(do_sdrep == TRUE) {
            tryCatch({
              obj$sd_rep <- RTMB::sdreport(obj) # get sdreport
              result[[length(what) + 1]] <- obj$sd_rep # input sd report
            }, error = function(e) {
              result[[length(what) + 1]] <- NA
            })
          }

          p() # update progress
          return(result)

        }, error = function(e) {
          # Skip failed simulations, saying why
          warning(sprintf("simulation %d failed: %s", i, conditionMessage(e)), call. = FALSE)
          result <- list()
          for(j in 1:length(what)) result[[what[j]]] <- NA
          if(do_sdrep == TRUE) result[[length(what) + 1]] <- NA

          p() # update progress
          return(result)
        })

      }, future.seed = TRUE)

      future::plan(future::sequential)  # Reset
    })

    # Populate results from parallel run
    for(i in 1:n_sims) for(j in 1:length(what)) store_res_list[[j]][[i]] <- sim_results[[i]][[what[j]]]
    if(do_sdrep == TRUE) for(i in 1:n_sims) store_res_list[[length(what) + 1]][[i]] <- sim_results[[i]][[length(what) + 1]]
    for(j in 1:length(what)) store_res_list[[j]] <- simplify2array(store_res_list[[j]])  # Convert lists to array
  }

  return(store_res_list)

}

#' Extract simulation outputs into SPoRC estimation model format
#'
#' Subsets and reshapes biological, tagging, fishery, and survey arrays from a
#' simulation environment or output list to cover years \code{1:y} and
#' simulation replicate \code{sim}, producing a named list ready for direct
#' use in \code{\link{Setup_Mod_Biologicals}}, \code{\link{Setup_Mod_Catch_and_F}},
#' \code{\link{Setup_Mod_SrvIdx_and_Comps}}, and \code{\link{Setup_Mod_Tagging}}.
#' Binary \code{Use*} indicator arrays are derived automatically from the
#' extracted observation arrays (non-NA, positive values set to 1).
#'
#' Population-specific arrays (\code{ObsCatch_pop}, \code{ObsFishIdx_pop},
#' \code{ObsFishAgeComps_pop}, \code{ObsFishLenComps_pop}, \code{ObsSrvIdx_pop},
#' \code{ObsSrvAgeComps_pop}, \code{ObsSrvLenComps_pop}) are extracted when
#' present in \code{sim_env}; corresponding \code{Use*_PopSpec} flags are
#' derived automatically.
#'
#' Input sample sizes for age and length compositions are extracted for both
#' aggregate and population-specific data streams, including retained
#' (\code{ISS_FishAgeComps}, \code{ISS_FishLenComps}, \code{ISS_FishAgeComps_pop},
#' \code{ISS_FishLenComps_pop}, \code{ISS_SrvAgeComps}, \code{ISS_SrvLenComps},
#' \code{ISS_SrvAgeComps_pop}, \code{ISS_SrvLenComps_pop}) and discard
#' (\code{ISS_FishAgeComps_discard}, \code{ISS_FishLenComps_discard},
#' \code{ISS_FishAgeComps_discard_pop}, \code{ISS_FishLenComps_discard_pop})
#' data streams.
#'
#' Length composition outputs (\code{ObsFishLenComps}, \code{ObsSrvLenComps},
#' and their population-specific and discard variants) and \code{SizeAgeTrans}
#' are \code{NULL} when no size-age transition matrix is present in
#' \code{sim_env}. Tagging outputs are \code{NULL} when
#' \code{use_conv_fish_tagging = 0}; otherwise, only cohorts with release
#' years in \code{1:y} are retained.
#'
#' @param sim_env Simulation environment or list (e.g., output from
#'   \code{\link{Simulate_Pop_Static}} or a \code{\link{Setup_sim_env}}
#'   environment) containing all operating model arrays.
#' @param y Integer. Last year to include; years \code{1:y} are retained.
#' @param sim Integer. Simulation replicate index to extract.
#'
#' @return Named list with the following elements (all arrays have \code{y}
#'   in the year dimension unless noted):
#'   \code{WAA} \code{[n_pop x n_regions x y x n_seas x n_ages x n_sexes]},
#'   \code{WAA_fish} \code{[... x n_fish_fleets]},
#'   \code{WAA_srv} \code{[... x n_srv_fleets]},
#'   \code{MatAA}, \code{SizeAgeTrans} (or \code{NULL}),
#'   \code{AgeingError} \code{[y x n_obs_ages x n_ages]},
#'   \code{use_conv_fish_tagging},
#'   \code{conv_tag_release_indicator}, \code{obs_conv_tag_fish_recap},
#'   \code{conv_tagged_fish}, \code{conv_tagged_fish_attr},
#'   \code{n_tag_cohorts} (all \code{NULL} when tagging inactive),
#'   \code{ObsCatch}, \code{ln_sigmaC}, \code{UseCatch},
#'   \code{ObsCatch_pop}, \code{ln_sigmaC_pop}, \code{UseCatch_pop},
#'   \code{ObsDiscard}, \code{ln_sigmaD}, \code{UseDiscard},
#'   \code{ObsDiscard_pop}, \code{ln_sigmaD_pop}, \code{UseDiscard_pop},
#'   \code{ObsFishIdx}, \code{ObsFishIdx_SE}, \code{UseFishIdx},
#'   \code{ObsFishIdx_pop}, \code{ObsFishIdx_pop_SE}, \code{UseFishIdx_pop},
#'   \code{ObsFishAgeComps}, \code{ISS_FishAgeComps}, \code{UseFishAgeComps},
#'   \code{ObsFishAgeComps_pop}, \code{ISS_FishAgeComps_pop},
#'   \code{UseFishAgeComps_pop},
#'   \code{ObsFishLenComps} (or \code{NULL}), \code{ISS_FishLenComps} (or \code{NULL}),
#'   \code{UseFishLenComps},
#'   \code{ObsFishLenComps_pop} (or \code{NULL}), \code{ISS_FishLenComps_pop} (or \code{NULL}),
#'   \code{UseFishLenComps_pop},
#'   \code{ObsFishAgeComps_discard}, \code{ISS_FishAgeComps_discard}, \code{UseFishAgeComps_discard},
#'   \code{ObsFishAgeComps_discard_pop}, \code{ISS_FishAgeComps_discard_pop},
#'   \code{UseFishAgeComps_discard_pop},
#'   \code{ObsFishLenComps_discard} (or \code{NULL}), \code{ISS_FishLenComps_discard} (or \code{NULL}),
#'   \code{UseFishLenComps_discard},
#'   \code{ObsFishLenComps_discard_pop} (or \code{NULL}), \code{ISS_FishLenComps_discard_pop} (or \code{NULL}),
#'   \code{UseFishLenComps_discard_pop},
#'   \code{ObsSrvIdx}, \code{ObsSrvIdx_SE}, \code{UseSrvIdx},
#'   \code{ObsSrvIdx_pop}, \code{ObsSrvIdx_pop_SE}, \code{UseSrvIdx_pop},
#'   \code{ObsSrvAgeComps}, \code{ISS_SrvAgeComps}, \code{UseSrvAgeComps},
#'   \code{ObsSrvAgeComps_pop}, \code{ISS_SrvAgeComps_pop},
#'   \code{UseSrvAgeComps_pop},
#'   \code{ObsSrvLenComps} (or \code{NULL}), \code{ISS_SrvLenComps} (or \code{NULL}),
#'   \code{UseSrvLenComps},
#'   \code{ObsSrvLenComps_pop} (or \code{NULL}), \code{ISS_SrvLenComps_pop} (or \code{NULL}),
#'   \code{UseSrvLenComps_pop}.
#'
#' @seealso \code{\link{Setup_Mod_Biologicals}}, \code{\link{Setup_Mod_Catch_and_F}},
#'   \code{\link{Setup_Mod_SrvIdx_and_Comps}}, \code{\link{Setup_Mod_Tagging}},
#'   \code{\link{Simulate_Pop_Static}}, \code{\link{Setup_sim_env}}
#'
#' @export simulation_data_to_SPoRC
#' @family Simulation Utilities
simulation_data_to_SPoRC <- function(sim_env,
                                     y,
                                     sim) {

  # Biologicals
  WAA <- array(sim_env$WAA[,,1:y,,,,sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_ages, sim_env$n_sexes))
  WAA_fish <- array(sim_env$WAA_fish[,,1:y,,,,,sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_ages, sim_env$n_sexes, sim_env$n_fish_fleets))
  WAA_srv <- array(sim_env$WAA_srv[,,1:y,,,,,sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_ages, sim_env$n_sexes, sim_env$n_srv_fleets))
  MatAA <- array(sim_env$MatAA[,,1:y,,,,sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_ages, sim_env$n_sexes))
  SizeAgeTrans <- if(!is.null(sim_env$SizeAgeTrans)) {
    array(sim_env$SizeAgeTrans[,,1:y,,,,,sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_lens, sim_env$n_ages, sim_env$n_sexes))
  } else NULL
  AgeingError <- array(sim_env$AgeingError[1:y,,,sim, drop = FALSE],
                       dim = c(length(1:y), dim(sim_env$AgeingError)[2], dim(sim_env$AgeingError)[3]))
  # the fleet-specific matrices the estimation model reads, peeled to the same years
  AgeingError_fish <- if(is.null(sim_env$AgeingError_fish)) NULL else {
    array(sim_env$AgeingError_fish[1:y,,,,sim, drop = FALSE], dim = c(length(1:y), dim(sim_env$AgeingError_fish)[2:4]))
  }
  AgeingError_srv <- if(is.null(sim_env$AgeingError_srv)) NULL else {
    array(sim_env$AgeingError_srv[1:y,,,,sim, drop = FALSE], dim = c(length(1:y), dim(sim_env$AgeingError_srv)[2:4]))
  }

  # Tagging
  if(sim_env$use_conv_fish_tagging == 1) {
    keep_tag_cohorts <- which(sim_env$conv_tag_release_indicator[,2] %in% 1:y)
    conv_tag_release_indicator <- sim_env$conv_tag_release_indicator[keep_tag_cohorts,,drop = FALSE]
    obs_conv_tag_fish_recap <- array(sim_env$obs_conv_tag_fish_recap[,,keep_tag_cohorts,,,,,,sim], dim = dim(sim_env$obs_conv_tag_fish_recap)[-length(dim(sim_env$obs_conv_tag_fish_recap))])
    conv_tagged_fish <- array(sim_env$conv_tagged_fish[keep_tag_cohorts,,,,sim], dim = c(dim(sim_env$conv_tagged_fish)[-length(dim(sim_env$conv_tagged_fish))]))
    conv_tagged_fish_attr <- array(sim_env$conv_tagged_fish_attr[keep_tag_cohorts,,,,sim], dim = c(dim(sim_env$conv_tagged_fish_attr)[-length(dim(sim_env$conv_tagged_fish_attr))]))
    n_tag_cohorts <- nrow(conv_tag_release_indicator)
  } else {
    conv_tag_release_indicator = obs_conv_tag_fish_recap = conv_tagged_fish = conv_tagged_fish_attr = n_tag_cohorts = NULL
  }

  # Fishery Landed Catches
  ObsCatch <- array(sim_env$ObsCatch[,1:y,,,sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))
  ln_sigmaC <- array(sim_env$ln_sigmaC[,1:y,,, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))
  UseCatch <- array(0, dim = dim(ObsCatch))
  UseCatch[!is.na(ObsCatch) & ObsCatch > 0] <- 1

  # Population-specific catches
  ObsCatch_pop <- array(sim_env$ObsCatch_pop[,,1:y,,,sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))
  UseCatch_pop <- array(0, dim = dim(ObsCatch_pop))
  UseCatch_pop[!is.na(ObsCatch_pop) & ObsCatch_pop > 0] <- 1
  ln_sigmaC_pop <- array(sim_env$ln_sigmaC_pop[,,1:y,,, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))

  # Discards
  ObsDiscard <- array(sim_env$ObsDiscard[,1:y,,,sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))
  ln_sigmaD <- array(sim_env$ln_sigmaD[,1:y,,, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))
  UseDiscard <- array(0, dim = dim(ObsDiscard))
  UseDiscard[!is.na(ObsDiscard) & ObsDiscard > 0] <- 1

  # Population-specific discards
  ObsDiscard_pop <- array(sim_env$ObsDiscard_pop[,,1:y,,,sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))
  ln_sigmaD_pop <- array(sim_env$ln_sigmaD_pop[,,1:y,,, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))
  UseDiscard_pop <- array(0, dim = dim(ObsDiscard_pop))
  UseDiscard_pop[!is.na(ObsDiscard_pop) & ObsDiscard_pop > 0] <- 1

  # Fishery Indices
  ObsFishIdx <- array(sim_env$ObsFishIdx[,1:y,,,sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))
  ObsFishIdx_SE <- array(sim_env$ObsFishIdx_SE[,1:y,,, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))
  UseFishIdx <- array(0, dim = dim(ObsFishIdx))
  UseFishIdx[!is.na(ObsFishIdx) & ObsFishIdx > 0] <- 1

  # Population-specific fishery indices
  ObsFishIdx_pop <- array(sim_env$ObsFishIdx_pop[,,1:y,,,sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))
  UseFishIdx_pop <- array(0, dim = dim(ObsFishIdx_pop))
  UseFishIdx_pop[!is.na(ObsFishIdx_pop) & ObsFishIdx_pop > 0] <- 1
  ObsFishIdx_pop_SE <- array(sim_env$ObsFishIdx_pop_SE[,,1:y,,, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))

  # Retained Fishery Compositions
  ObsFishAgeComps <- array(sim_env$ObsFishAgeComps[,1:y,,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, dim(sim_env$AgeingError)[3], sim_env$n_sexes, sim_env$n_fish_fleets))
  ISS_FishAgeComps <- array(sim_env$ISS_FishAgeComps[,1:y,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_sexes, sim_env$n_fish_fleets))
  UseFishAgeComps <- apply(ObsFishAgeComps, c(1,2,3,6), sum)
  UseFishAgeComps[!is.na(UseFishAgeComps) & UseFishAgeComps > 0] <- 1

  ObsFishLenComps <- if(!is.null(sim_env$n_lens)) {
    array(sim_env$ObsFishLenComps[,1:y,,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_lens, sim_env$n_sexes, sim_env$n_fish_fleets))
  } else NULL
  ISS_FishLenComps <- if(!is.null(sim_env$n_lens)) {
    array(sim_env$ISS_FishLenComps[,1:y,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_sexes, sim_env$n_fish_fleets))
  } else NULL
  if(!is.null(sim_env$n_lens)) {
    UseFishLenComps <- apply(ObsFishLenComps, c(1,2,3,6), sum)
    UseFishLenComps[!is.na(UseFishLenComps) & UseFishLenComps > 0] <- 1
  } else UseFishLenComps <- array(0, dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))

  # Population-specific retained fishery compositions
  ObsFishAgeComps_pop <- array(sim_env$ObsFishAgeComps_pop[,,1:y,,,,,sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, dim(sim_env$AgeingError)[3], sim_env$n_sexes, sim_env$n_fish_fleets))
  ISS_FishAgeComps_pop <- array(sim_env$ISS_FishAgeComps_pop[,,1:y,,,, sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_sexes, sim_env$n_fish_fleets))
  UseFishAgeComps_pop <- apply(ObsFishAgeComps_pop, c(1,2,3,4,7), sum)
  UseFishAgeComps_pop[!is.na(UseFishAgeComps_pop) & UseFishAgeComps_pop > 0] <- 1

  ObsFishLenComps_pop <- if(!is.null(sim_env$n_lens)) {
    array(sim_env$ObsFishLenComps_pop[,,1:y,,,,,sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_lens, sim_env$n_sexes, sim_env$n_fish_fleets))
  } else NULL
  ISS_FishLenComps_pop <- if(!is.null(sim_env$n_lens)) {
    array(sim_env$ISS_FishLenComps_pop[,,1:y,,,, sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_sexes, sim_env$n_fish_fleets))
  } else NULL
  if(!is.null(sim_env$n_lens)) {
    UseFishLenComps_pop <- apply(ObsFishLenComps_pop, c(1,2,3,4,7), sum)
    UseFishLenComps_pop[!is.na(UseFishLenComps_pop) & UseFishLenComps_pop > 0] <- 1
  } else UseFishLenComps_pop <- array(0, dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))

  # Discarded Fishery Compositions
  ObsFishAgeComps_discard <- array(sim_env$ObsFishAgeComps_discard[,1:y,,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, dim(sim_env$AgeingError)[3], sim_env$n_sexes, sim_env$n_fish_fleets))
  ISS_FishAgeComps_discard <- array(sim_env$ISS_FishAgeComps_discard[,1:y,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_sexes, sim_env$n_fish_fleets))
  UseFishAgeComps_discard <- apply(ObsFishAgeComps_discard, c(1,2,3,6), sum)
  UseFishAgeComps_discard[!is.na(UseFishAgeComps_discard) & UseFishAgeComps_discard > 0] <- 1

  ObsFishLenComps_discard <- if(!is.null(sim_env$n_lens)) {
    array(sim_env$ObsFishLenComps_discard[,1:y,,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_lens, sim_env$n_sexes, sim_env$n_fish_fleets))
  } else NULL
  ISS_FishLenComps_discard <- if(!is.null(sim_env$n_lens)) {
    array(sim_env$ISS_FishLenComps_discard[,1:y,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_sexes, sim_env$n_fish_fleets))
  } else NULL
  if(!is.null(sim_env$n_lens)) {
    UseFishLenComps_discard <- apply(ObsFishLenComps_discard, c(1,2,3,6), sum)
    UseFishLenComps_discard[!is.na(UseFishLenComps_discard) & UseFishLenComps_discard > 0] <- 1
  } else UseFishLenComps_discard <- array(0, dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))

  # Population-specific discarded fishery compositions
  ObsFishAgeComps_discard_pop <- array(sim_env$ObsFishAgeComps_discard_pop[,,1:y,,,,,sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, dim(sim_env$AgeingError)[3], sim_env$n_sexes, sim_env$n_fish_fleets))
  ISS_FishAgeComps_discard_pop <- array(sim_env$ISS_FishAgeComps_discard_pop[,,1:y,,,, sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_sexes, sim_env$n_fish_fleets))
  UseFishAgeComps_discard_pop <- apply(ObsFishAgeComps_discard_pop, c(1,2,3,4,7), sum)
  UseFishAgeComps_discard_pop[!is.na(UseFishAgeComps_discard_pop) & UseFishAgeComps_discard_pop > 0] <- 1

  ObsFishLenComps_discard_pop <- if(!is.null(sim_env$n_lens)) {
    array(sim_env$ObsFishLenComps_discard_pop[,,1:y,,,,,sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_lens, sim_env$n_sexes, sim_env$n_fish_fleets))
  } else NULL
  ISS_FishLenComps_discard_pop <- if(!is.null(sim_env$n_lens)) {
    array(sim_env$ISS_FishLenComps_discard_pop[,,1:y,,,, sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_sexes, sim_env$n_fish_fleets))
  } else NULL
  if(!is.null(sim_env$n_lens)) {
    UseFishLenComps_discard_pop <- apply(ObsFishLenComps_discard_pop, c(1,2,3,4,7), sum)
    UseFishLenComps_discard_pop[!is.na(UseFishLenComps_discard_pop) & UseFishLenComps_discard_pop > 0] <- 1
  } else UseFishLenComps_discard_pop <- array(0, dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))

  # Survey Indices
  ObsSrvIdx <- array(sim_env$ObsSrvIdx[,1:y,,,sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_srv_fleets))
  ObsSrvIdx_SE <- array(sim_env$ObsSrvIdx_SE[,1:y,,, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_srv_fleets))
  UseSrvIdx <- array(0, dim = dim(ObsSrvIdx))
  UseSrvIdx[!is.na(ObsSrvIdx) & ObsSrvIdx > 0] <- 1

  # Population-specific survey indices
  ObsSrvIdx_pop <- array(sim_env$ObsSrvIdx_pop[,,1:y,,,sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_srv_fleets))
  UseSrvIdx_pop <- array(0, dim = dim(ObsSrvIdx_pop))
  UseSrvIdx_pop[!is.na(ObsSrvIdx_pop) & ObsSrvIdx_pop > 0] <- 1
  ObsSrvIdx_pop_SE <- array(sim_env$ObsSrvIdx_pop_SE[,,1:y,,, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_srv_fleets))

  # Survey Compositions
  ObsSrvAgeComps <- array(sim_env$ObsSrvAgeComps[,1:y,,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, dim(sim_env$AgeingError)[3], sim_env$n_sexes, sim_env$n_srv_fleets))
  ISS_SrvAgeComps <- array(sim_env$ISS_SrvAgeComps[,1:y,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_sexes, sim_env$n_srv_fleets))
  UseSrvAgeComps <- apply(ObsSrvAgeComps, c(1,2,3,6), sum)
  UseSrvAgeComps[!is.na(UseSrvAgeComps) & UseSrvAgeComps > 0] <- 1

  ObsSrvLenComps <- if(!is.null(sim_env$n_lens)) {
    array(sim_env$ObsSrvLenComps[,1:y,,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_lens, sim_env$n_sexes, sim_env$n_srv_fleets))
  } else NULL
  ISS_SrvLenComps <- if(!is.null(sim_env$n_lens)) {
    array(sim_env$ISS_SrvLenComps[,1:y,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_sexes, sim_env$n_srv_fleets))
  } else NULL
  if(!is.null(sim_env$n_lens)) {
    UseSrvLenComps <- apply(ObsSrvLenComps, c(1,2,3,6), sum)
    UseSrvLenComps[!is.na(UseSrvLenComps) & UseSrvLenComps > 0] <- 1
  } else UseSrvLenComps <- array(0, dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_srv_fleets))

  # Population-specific survey compositions
  ObsSrvAgeComps_pop <- array(sim_env$ObsSrvAgeComps_pop[,,1:y,,,,,sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, dim(sim_env$AgeingError)[3], sim_env$n_sexes, sim_env$n_srv_fleets))
  ISS_SrvAgeComps_pop <- array(sim_env$ISS_SrvAgeComps_pop[,,1:y,,,, sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_sexes, sim_env$n_srv_fleets))
  UseSrvAgeComps_pop <- apply(ObsSrvAgeComps_pop, c(1,2,3,4,7), sum)
  UseSrvAgeComps_pop[!is.na(UseSrvAgeComps_pop) & UseSrvAgeComps_pop > 0] <- 1

  ObsSrvLenComps_pop <- if(!is.null(sim_env$n_lens)) {
    array(sim_env$ObsSrvLenComps_pop[,,1:y,,,,,sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_lens, sim_env$n_sexes, sim_env$n_srv_fleets))
  } else NULL
  ISS_SrvLenComps_pop <- if(!is.null(sim_env$n_lens)) {
    array(sim_env$ISS_SrvLenComps_pop[,,1:y,,,, sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_sexes, sim_env$n_srv_fleets))
  } else NULL
  if(!is.null(sim_env$n_lens)) {
    UseSrvLenComps_pop <- apply(ObsSrvLenComps_pop, c(1,2,3,4,7), sum)
    UseSrvLenComps_pop[!is.na(UseSrvLenComps_pop) & UseSrvLenComps_pop > 0] <- 1
  } else UseSrvLenComps_pop <- array(0, dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_srv_fleets))

  # Conditional age-at-length, present only when the simulation drew it. The use
  # flag marks length bins that received at least one aged fish.
  n_obs_ages <- dim(sim_env$AgeingError)[3]
  if(isTRUE(sim_env$do_fish_caal)) {
    ObsFish_caal <- array(sim_env$ObsFish_caal[,1:y,,,,,,sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_lens, n_obs_ages, sim_env$n_sexes, sim_env$n_fish_fleets))
    ISS_Fish_caal <- array(sim_env$ISS_Fish_caal[,1:y,,,,,sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_lens, sim_env$n_sexes, sim_env$n_fish_fleets))
    UseFish_caal <- apply(ObsFish_caal, c(1,2,3,4,7), sum)
    UseFish_caal[] <- as.numeric(!is.na(UseFish_caal) & UseFish_caal > 0)
  } else ObsFish_caal <- ISS_Fish_caal <- UseFish_caal <- NULL
  if(isTRUE(sim_env$do_srv_caal)) {
    ObsSrv_caal <- array(sim_env$ObsSrv_caal[,1:y,,,,,,sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_lens, n_obs_ages, sim_env$n_sexes, sim_env$n_srv_fleets))
    ISS_Srv_caal <- array(sim_env$ISS_Srv_caal[,1:y,,,,,sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_lens, sim_env$n_sexes, sim_env$n_srv_fleets))
    UseSrv_caal <- apply(ObsSrv_caal, c(1,2,3,4,7), sum)
    UseSrv_caal[] <- as.numeric(!is.na(UseSrv_caal) & UseSrv_caal > 0)
  } else ObsSrv_caal <- ISS_Srv_caal <- UseSrv_caal <- NULL

  # Return
  return(list(
    # Biologicals
    WAA = WAA,
    WAA_fish = WAA_fish,
    WAA_srv = WAA_srv,
    MatAA = MatAA,
    SizeAgeTrans = SizeAgeTrans,
    AgeingError = AgeingError,
    AgeingError_fish = AgeingError_fish,
    AgeingError_srv = AgeingError_srv,

    # Tagging
    use_conv_fish_tagging = sim_env$use_conv_fish_tagging,
    conv_tag_release_indicator = conv_tag_release_indicator,
    obs_conv_tag_fish_recap = obs_conv_tag_fish_recap,
    conv_tagged_fish = conv_tagged_fish,
    conv_tagged_fish_attr = conv_tagged_fish_attr,
    n_tag_cohorts = n_tag_cohorts,

    # Aggregated catches
    ObsCatch = ObsCatch,
    ln_sigmaC = ln_sigmaC,
    UseCatch = UseCatch,
    ObsDiscard = ObsDiscard,
    ln_sigmaD = ln_sigmaD,
    UseDiscard = UseDiscard,

    # Population-specific catches
    ObsCatch_pop = ObsCatch_pop,
    ln_sigmaC_pop = ln_sigmaC_pop,
    UseCatch_pop = UseCatch_pop,
    ObsDiscard_pop = ObsDiscard_pop,
    ln_sigmaD_pop = ln_sigmaD_pop,
    UseDiscard_pop = UseDiscard_pop,

    # Aggregated fishery indices
    ObsFishIdx = ObsFishIdx,
    ObsFishIdx_SE = ObsFishIdx_SE,
    UseFishIdx = UseFishIdx,

    # Population-specific fishery indices
    ObsFishIdx_pop = ObsFishIdx_pop,
    ObsFishIdx_pop_SE = ObsFishIdx_pop_SE,
    UseFishIdx_pop = UseFishIdx_pop,

    # Aggregated retained fishery compositions
    ObsFishAgeComps = ObsFishAgeComps,
    ISS_FishAgeComps = ISS_FishAgeComps,
    UseFishAgeComps = UseFishAgeComps,
    ObsFishLenComps = ObsFishLenComps,
    ISS_FishLenComps = ISS_FishLenComps,
    UseFishLenComps = UseFishLenComps,

    # Population-specific retained fishery compositions
    ObsFishAgeComps_pop = ObsFishAgeComps_pop,
    ISS_FishAgeComps_pop = ISS_FishAgeComps_pop,
    UseFishAgeComps_pop = UseFishAgeComps_pop,
    ObsFishLenComps_pop = ObsFishLenComps_pop,
    ISS_FishLenComps_pop = ISS_FishLenComps_pop,
    UseFishLenComps_pop = UseFishLenComps_pop,

    # Aggregated discarded fishery compositions
    ObsFishAgeComps_discard = ObsFishAgeComps_discard,
    ISS_FishAgeComps_discard = ISS_FishAgeComps_discard,
    UseFishAgeComps_discard = UseFishAgeComps_discard,
    ObsFishLenComps_discard = ObsFishLenComps_discard,
    ISS_FishLenComps_discard = ISS_FishLenComps_discard,
    UseFishLenComps_discard = UseFishLenComps_discard,

    # Population-specific discarded fishery compositions
    ObsFishAgeComps_discard_pop = ObsFishAgeComps_discard_pop,
    ISS_FishAgeComps_discard_pop = ISS_FishAgeComps_discard_pop,
    UseFishAgeComps_discard_pop = UseFishAgeComps_discard_pop,
    ObsFishLenComps_discard_pop = ObsFishLenComps_discard_pop,
    ISS_FishLenComps_discard_pop = ISS_FishLenComps_discard_pop,
    UseFishLenComps_discard_pop = UseFishLenComps_discard_pop,

    # Aggregated survey indices
    ObsSrvIdx = ObsSrvIdx,
    ObsSrvIdx_SE = ObsSrvIdx_SE,
    UseSrvIdx = UseSrvIdx,

    # Population-specific survey indices
    ObsSrvIdx_pop = ObsSrvIdx_pop,
    ObsSrvIdx_pop_SE = ObsSrvIdx_pop_SE,
    UseSrvIdx_pop = UseSrvIdx_pop,

    # Aggregated survey compositions
    ObsSrvAgeComps = ObsSrvAgeComps,
    ISS_SrvAgeComps = ISS_SrvAgeComps,
    UseSrvAgeComps = UseSrvAgeComps,
    ObsSrvLenComps = ObsSrvLenComps,
    ISS_SrvLenComps = ISS_SrvLenComps,
    UseSrvLenComps = UseSrvLenComps,

    # Population-specific survey compositions
    ObsSrvAgeComps_pop = ObsSrvAgeComps_pop,
    ISS_SrvAgeComps_pop = ISS_SrvAgeComps_pop,
    UseSrvAgeComps_pop = UseSrvAgeComps_pop,
    ObsSrvLenComps_pop = ObsSrvLenComps_pop,
    ISS_SrvLenComps_pop = ISS_SrvLenComps_pop,
    UseSrvLenComps_pop = UseSrvLenComps_pop,

    # Conditional age-at-length
    ObsFish_caal = ObsFish_caal,
    ISS_Fish_caal = ISS_Fish_caal,
    UseFish_caal = UseFish_caal,
    ObsSrv_caal = ObsSrv_caal,
    ISS_Srv_caal = ISS_Srv_caal,
    UseSrv_caal = UseSrv_caal
  ))

}
