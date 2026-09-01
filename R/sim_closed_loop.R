# Operating model
#
# Closed loop simulation: condition an operating model on a fitted assessment,
# then run the assessment and a harvest control rule forward against it. The
# catch_to_F_* helpers invert a catch target to a fishing mortality, which is what
# a control rule needs to hand back to the operating model.

#' Construct and Condition Closed-Loop Simulation Inputs
#'
#' Initializes and conditions a simulation object used for closed-loop
#' population projections. The function reconstructs fitted model quantities,
#' extends time-varying processes into projection years, and prepares all
#' simulation components required by the operating model.
#'
#' Simulation inputs are generated for biological dynamics, fishing,
#' surveys, recruitment, tagging, and movement processes. Historical years
#' are conditioned on outputs from the fitted assessment model, while
#' projection years are initialized using extension rules or user-supplied
#' overrides.
#'
#' Users may replace any internally generated component by supplying a
#' correctly named object via \code{...}.
#'
#' @param closed_loop_yrs Integer. Number of projection years added beyond
#'   the fitted data period.
#' @param n_sims Integer. Number of stochastic simulation replicates.
#' @param data List. Data object used to fit the assessment model.
#' @param parameters List. Parameter vector from the fitted model.
#' @param mapping List. Parameter mapping object used during estimation.
#' @param sd_rep List. Standard deviation report from the fitted model.
#' @param rep List. Model report object produced by the fitted model.
#' @param random Character vector of estimated random effects.
#'
#' @param FishIdx_SE_fill Character or numeric specifying how pooled fishery
#'   index standard errors are extended into projection years.
#' @param SrvIdx_SE_fill Character or numeric specifying how pooled survey
#'   index standard errors are extended into projection years.
#' @param FishIdx_SE_pop_fill Character or numeric specifying how
#'   population-specific fishery index standard errors are extended into
#'   projection years. Default \code{"mean"}.
#' @param SrvIdx_SE_pop_fill Character or numeric specifying how
#'   population-specific survey index standard errors are extended into
#'   projection years. Default \code{"mean"}.
#' @param ISS_FishAgeComps_fill Character or numeric specifying how pooled
#'   fishery age-composition input sample sizes are extended into projection
#'   years.
#' @param ISS_FishLenComps_fill Same behavior as \code{ISS_FishAgeComps_fill}
#'   for pooled fishery length compositions.
#' @param ISS_SrvAgeComps_fill Character or numeric specifying how pooled
#'   survey age-composition input sample sizes are extended into projection
#'   years.
#' @param ISS_SrvLenComps_fill Same behavior as \code{ISS_SrvAgeComps_fill}
#'   for pooled survey length compositions.
#' @param ISS_FishAgeComps_pop_fill Character or numeric specifying how
#'   population-specific fishery age-composition input sample sizes are
#'   extended into projection years. Default \code{"mean"}.
#' @param ISS_FishLenComps_pop_fill Same behavior as
#'   \code{ISS_FishAgeComps_pop_fill} for population-specific fishery length
#'   compositions. Default \code{"mean"}.
#' @param ISS_SrvAgeComps_pop_fill Character or numeric specifying how
#'   population-specific survey age-composition input sample sizes are
#'   extended into projection years. Default \code{"mean"}.
#' @param ISS_SrvLenComps_pop_fill Same behavior as
#'   \code{ISS_SrvAgeComps_pop_fill} for population-specific survey length
#'   compositions. Default \code{"mean"}.
#' @param ISS_FishAgeComps_discard_fill Same behavior as \code{ISS_FishAgeComps_fill}
#'   for pooled fishery length compositions.
#' @param ISS_FishLenComps_discard_fill Same behavior as \code{ISS_FishAgeComps_fill}
#'   for pooled fishery length compositions.
#' @param ISS_FishAgeComps_discard_pop_fill Same behavior as
#'   \code{ISS_FishAgeComps_pop_fill} for population-specific fishery length
#'   compositions. Default \code{"mean"}.
#' @param ISS_FishLenComps_discard_pop_fill Same behavior as
#'   \code{ISS_FishLenComps_pop_fill} for population-specific fishery length
#'   compositions. Default \code{"mean"}.
#'
#' Extension rules for all \code{*_fill} arguments may be:
#'
#' \describe{
#'   \item{\code{"zeros"}}{Fill projection years with zeros.}
#'   \item{\code{"last"}}{Repeat the final observed year.}
#'   \item{\code{"mean"}}{Use the mean of the historical series.}
#'   \item{\code{"F_pattern"}}{Scale fishery composition sample sizes
#'     proportionally to the simulated fishing mortality pattern (pooled
#'     fishery ISS only).}
#'   \item{numeric}{Use a supplied constant or array.}
#' }
#'
#' If a fill argument receives an array directly, it is interpreted as the
#' fully specified input and stored accordingly; the fill rule is ignored.
#'
#' @param ... Optional named simulation inputs that override internally
#'   generated components. Any argument expected by
#'   \code{\link{Setup_Sim_Fishing}}, \code{\link{Setup_Sim_Survey}},
#'   \code{\link{Setup_Sim_Biologicals}}, \code{\link{Setup_Sim_Rec}}, or
#'   \code{\link{Setup_Sim_Tagging}} may be provided.
#'
#'   Common overrides include:
#'
#'   \strong{Fishing processes}
#'   \itemize{
#'     \item \code{Fmort_input}
#'     \item \code{fish_sel_input}
#'     \item \code{fish_q_input}
#'   }
#'
#'   \strong{Survey processes}
#'   \itemize{
#'     \item \code{srv_sel_input}
#'     \item \code{srv_q_input}
#'   }
#'
#'   \strong{Biological processes}
#'   \itemize{
#'     \item \code{WAA_input}
#'     \item \code{MatAA_input}
#'     \item \code{natmort_input}
#'   }
#'
#'   \strong{Recruitment processes}
#'   \itemize{
#'     \item \code{R0_input}
#'     \item \code{rinit_input}
#'     \item \code{h_input}
#'     \item \code{Rec_input}
#'   }
#'
#'   \strong{Tagging processes}
#'   \itemize{
#'     \item \code{conv_tag_fish_reporting_input}
#'   }
#'
#'   \strong{Movement processes}
#'   \itemize{
#'     \item \code{Movement}
#'   }
#'
#'   Supplied objects must have dimensions consistent with model structure,
#'   the total number of years (\code{length(data$years) + closed_loop_yrs}),
#'   and the number of simulations (\code{n_sims}).
#'
#' @details
#' Simulation years consist of two periods:
#'
#' \itemize{
#'   \item \strong{Conditioning years}: historical years corresponding to
#'     the fitted assessment model.
#'   \item \strong{Projection years}: future years simulated under
#'     closed-loop management.
#' }
#'
#' During conditioning years, model processes are reconstructed from the
#' fitted model report objects. For projection years, quantities are extended
#' using the specified fill rules or user-supplied inputs.
#'
#' By default:
#'
#' \itemize{
#'   \item Biological inputs, selectivity, and catchability are extended
#'     using values from the final estimated year.
#'   \item Fishing mortality is initialized to zero in projection years.
#'   \item Recruitment is simulated forward when not fully specified by
#'     \code{Rec_input}.
#'   \item Population-specific data streams (\code{ObsFishIdx_pop_SE},
#'     \code{ObsSrvIdx_pop_SE}, and all \code{*_pop} ISS arrays) fall back
#'     to uninformative defaults when the corresponding \code{Use*_pop}
#'     flags contain no ones.
#' }
#'
#' Closed-loop feedback begins in the first projection year and allows
#' management actions (e.g., fishing mortality adjustments) to update
#' dynamically during the simulation.
#'
#' @return
#' A fully initialized \code{sim_list} object containing:
#'
#' \itemize{
#'   \item model dimensions and simulation containers
#'   \item biological process inputs
#'   \item pooled and population-specific fishing and survey processes
#'   \item recruitment dynamics
#'   \item tagging processes
#'   \item spatial movement structures
#'   \item replicated arrays across \code{n_sims}
#' }
#'
#' @export condition_closed_loop_simulations
#' @family Closed Loop Simulations
condition_closed_loop_simulations <- function(closed_loop_yrs,
                                              n_sims,
                                              data,
                                              parameters,
                                              mapping,
                                              sd_rep,
                                              rep,
                                              random = random,
                                              FishIdx_SE_fill = "mean",
                                              SrvIdx_SE_fill = "mean",
                                              FishIdx_SE_pop_fill = "mean",
                                              SrvIdx_SE_pop_fill = "mean",
                                              ISS_FishAgeComps_fill = "mean",
                                              ISS_FishLenComps_fill = "mean",
                                              ISS_FishAgeComps_discard_fill = "mean",
                                              ISS_FishLenComps_discard_fill = "mean",
                                              ISS_SrvAgeComps_fill = "mean",
                                              ISS_SrvLenComps_fill = "mean",
                                              ISS_FishAgeComps_pop_fill = "mean",
                                              ISS_FishLenComps_pop_fill = "mean",
                                              ISS_FishAgeComps_discard_pop_fill = "mean",
                                              ISS_FishLenComps_discard_pop_fill = "mean",
                                              ISS_SrvAgeComps_pop_fill = "mean",
                                              ISS_SrvLenComps_pop_fill = "mean",
                                              ...
                                              ) {

  # Additional user inputs as desired
  args <- list(...)

  # Detect partial matching: if *_fill received an array, it was meant as data
  if(is.array(ISS_FishAgeComps_fill)) {
    args$ISS_FishAgeComps <- ISS_FishAgeComps_fill
    ISS_FishAgeComps_fill <- "placeholder"
  }
  if(is.array(ISS_FishLenComps_fill)) {
    args$ISS_FishLenComps <- ISS_FishLenComps_fill
    ISS_FishLenComps_fill <- "placeholder"
  }
  if(is.array(ISS_FishAgeComps_discard_fill)) {
    args$ISS_FishAgeComps_discard <- ISS_FishAgeComps_discard_fill
    ISS_FishAgeComps_discard_fill <- "placeholder"
  }
  if(is.array(ISS_FishLenComps_discard_fill)) {
    args$ISS_FishLenComps_discard <- ISS_FishLenComps_discard_fill
    ISS_FishLenComps_discard_fill <- "placeholder"
  }
  if(is.array(ISS_SrvAgeComps_fill)) {
    args$ISS_SrvAgeComps <- ISS_SrvAgeComps_fill
    ISS_SrvAgeComps_fill <- "placeholder"
  }
  if(is.array(ISS_SrvLenComps_fill)) {
    args$ISS_SrvLenComps <- ISS_SrvLenComps_fill
    ISS_SrvLenComps_fill <- "placeholder"
  }
  if(is.array(FishIdx_SE_fill)) {
    args$ObsFishIdx_SE <- FishIdx_SE_fill
    FishIdx_SE_fill <- "placeholder"
  }
  if(is.array(SrvIdx_SE_fill)) {
    args$ObsSrvIdx_SE <- SrvIdx_SE_fill
    SrvIdx_SE_fill <- "placeholder"
  }
  if(is.array(ISS_FishAgeComps_pop_fill)) {
    args$ISS_FishAgeComps_pop <- ISS_FishAgeComps_pop_fill
    ISS_FishAgeComps_pop_fill <- "placeholder"
  }
  if(is.array(ISS_FishLenComps_pop_fill)) {
    args$ISS_FishLenComps_pop <- ISS_FishLenComps_pop_fill
    ISS_FishLenComps_pop_fill <- "placeholder"
  }
  if(is.array(ISS_FishAgeComps_discard_pop_fill)) {
    args$ISS_FishAgeComps_discard_pop <- ISS_FishAgeComps_discard_pop_fill
    ISS_FishAgeComps_discard_pop_fill <- "placeholder"
  }
  if(is.array(ISS_FishLenComps_discard_pop_fill)) {
    args$ISS_FishLenComps_discard_pop <- ISS_FishLenComps_discard_pop_fill
    ISS_FishLenComps_discard_pop_fill <- "placeholder"
  }
  if(is.array(ISS_SrvAgeComps_pop_fill)) {
    args$ISS_SrvAgeComps_pop <- ISS_SrvAgeComps_pop_fill
    ISS_SrvAgeComps_pop_fill <- "placeholder"
  }
  if(is.array(ISS_SrvLenComps_pop_fill)) {
    args$ISS_SrvLenComps_pop <- ISS_SrvLenComps_pop_fill
    ISS_SrvLenComps_pop_fill <- "placeholder"
  }
  if(is.array(FishIdx_SE_pop_fill)) {
    args$ObsFishIdx_pop_SE <- FishIdx_SE_pop_fill
    FishIdx_SE_pop_fill <- "placeholder"
  }
  if(is.array(SrvIdx_SE_pop_fill)) {
    args$ObsSrvIdx_pop_SE <- SrvIdx_SE_pop_fill
    SrvIdx_SE_pop_fill <- "placeholder"
  }

  optim_parameters_list <- get_optim_param_list(parameters, mapping, sd_rep, random) # get optimized parameters in original list format

  # Setup Model Dimensions --------------------------------------------------
  sim_list <- Setup_Sim_Dim(n_sims = n_sims,
                            n_yrs = length(data$years) + closed_loop_yrs,
                            n_regions = data$n_regions,
                            n_ages = length(data$ages),
                            n_obs_ages = if(any(data$UseFishAgeComps == 1)) {
                              dim(data$ObsFishAgeComps)[4]
                            } else if(any(data$UseFishAgeComps_pop == 1)) {
                              dim(data$ObsFishAgeComps_pop)[5]
                            } else if(any(data$UseSrvAgeComps == 1)) {
                              dim(data$ObsSrvAgeComps)[4]
                            } else if(any(data$UseSrvAgeComps_pop == 1)) {
                              dim(data$ObsSrvAgeComps_pop)[5]
                            },
                            n_lens = length(data$lens),
                            n_sexes = data$n_sexes,
                            n_fish_fleets = data$n_fish_fleets,
                            n_srv_fleets = data$n_srv_fleets,
                            feedback_start_yr = length(data$years),
                            n_seas = data$n_seas,
                            seasdur = data$seasdur,
                            n_pop = data$n_pop,
                            natal_region = data$natal_region,
                            run_feedback = TRUE
  )

  # Setup Simulation Containers ---------------------------------------------
  sim_list <- Setup_Sim_Containers(sim_list) # set up simulation containers to use

  # Setup Fishing Processes -------------------------------------------------
  # Catch uncertainty
  ln_sigmaC <- if(!"ln_sigmaC" %in% names(args)) {
    tmp <- array(NA, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets))
    for(r in 1:sim_list$n_regions) for(f in 1:sim_list$n_fish_fleets) {
      if(!is.vector(data$Wt_Catch)) {
        tmp[r,,,f] <- mean(log(exp(optim_parameters_list$ln_sigmaC[r,,,f]) / sqrt(data$Wt_Catch[r,,,f])))
      } else {
        tmp[r,,,f] <- mean(log(exp(optim_parameters_list$ln_sigmaC[r,,,f]) / sqrt(data$Wt_Catch)))
      }
    }
    tmp
  } else args$ln_sigmaC

  ln_sigmaC_pop <- if(!"ln_sigmaC_pop" %in% names(args)) {
    tmp <- array(NA, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets))
    for(p in 1:sim_list$n_pop) for(r in 1:sim_list$n_regions) for(f in 1:sim_list$n_fish_fleets) {
      if(!is.vector(data$Wt_Catch_pop)) {
        tmp[p,r,,,f] <- mean(log(exp(optim_parameters_list$ln_sigmaC_pop[p,r,,,f]) / sqrt(data$Wt_Catch_pop[p,r,,,f])))
      } else {
        tmp[p,r,,,f] <- mean(log(exp(optim_parameters_list$ln_sigmaC_pop[p,r,,,f]) / sqrt(data$Wt_Catch_pop)))
      }
    }
    tmp
  } else args$ln_sigmaC_pop

  # Catch uncertainty
  ln_sigmaD <- if(!"ln_sigmaD" %in% names(args)) {
    tmp <- array(NA, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets))
    for(r in 1:sim_list$n_regions) for(f in 1:sim_list$n_fish_fleets) {
      if(!is.vector(data$Wt_Discard)) {
        tmp[r,,,f] <- mean(log(exp(optim_parameters_list$ln_sigmaD[r,,,f]) / sqrt(data$Wt_Discard[r,,,f])))
      } else {
        tmp[r,,,f] <- mean(log(exp(optim_parameters_list$ln_sigmaD[r,,,f]) / sqrt(data$Wt_Discard)))
      }
    }
    tmp
  } else args$ln_sigmaD

  ln_sigmaD_pop <- if(!"ln_sigmaD_pop" %in% names(args)) {
    tmp <- array(NA, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets))
    for(p in 1:sim_list$n_pop) for(r in 1:sim_list$n_regions) for(f in 1:sim_list$n_fish_fleets) {
      if(!is.vector(data$Wt_Discard_pop)) {
        tmp[p,r,,,f] <- mean(log(exp(optim_parameters_list$ln_sigmaD_pop[p,r,,,f]) / sqrt(data$Wt_Discard_pop[p,r,,,f])))
      } else {
        tmp[p,r,,,f] <- mean(log(exp(optim_parameters_list$ln_sigmaD_pop[p,r,,,f]) / sqrt(data$Wt_Discard_pop)))
      }
    }
    tmp
  } else args$ln_sigmaD_pop

  # Fishery selectivity
  fish_sel_input <- if(!"fish_sel_input" %in% names(args)) {
    extend_years(replicate(n = sim_list$n_sims, rep$fish_sel[,,1:length(data$years),,,,,drop = FALSE]), n_years = closed_loop_yrs, 3, fill = 'last')
  } else args$fish_sel_input
  # Retained selectivity
  ret_sel_input <- if(!"ret_sel_input" %in% names(args)) {
    extend_years(replicate(n = sim_list$n_sims, rep$ret_sel[,,1:length(data$years),,,,,drop = FALSE]), n_years = closed_loop_yrs, 3, fill = 'last')
  } else args$ret_sel_input
  # Fishery catchability
  fish_q_input <- if(!"fish_q_input" %in% names(args)) {
    extend_years(replicate(n = sim_list$n_sims, rep$fish_q[,1:length(data$years),,drop = FALSE]), n_years = closed_loop_yrs, 2, fill = 'last')
  } else args$fish_q_input
  # Fishery index uncertainty
  ObsFishIdx_SE <- if(!"ObsFishIdx_SE" %in% names(args)) {
    extend_years(arr = data$ObsFishIdx_SE / sqrt(data$Wt_FishIdx), n_years = closed_loop_yrs, 2, fill = FishIdx_SE_fill)
  } else args$ObsFishIdx_SE

  # Fishery age compositions
  comp_fishage_like <- if(!"comp_fishage_like" %in% names(args)) data$FishAgeComps_LikeType else args$comp_fishage_like
  FishAgeComps_Type <- if(!"FishAgeComps_Type" %in% names(args)) extend_years(data$FishAgeComps_Type, closed_loop_yrs, 1, 'last') else args$FishAgeComps_Type
  ISS_FishAgeComps <- if(!"ISS_FishAgeComps" %in% names(args)) {
    extend_years(replicate(sim_list$n_sims, data$ISS_FishAgeComps[,,,,,drop = FALSE] * data$Wt_FishAgeComps), closed_loop_yrs, 2, fill = ISS_FishAgeComps_fill)
  } else args$ISS_FishAgeComps
  ln_FishAge_theta <- if(!"ln_FishAge_theta" %in% names(args)) optim_parameters_list$ln_FishAge_theta[,,,drop = FALSE] else args$ln_FishAge_theta
  ln_FishAge_theta_agg <- if(!"ln_FishAge_theta_agg" %in% names(args)) optim_parameters_list$ln_FishAge_theta_agg else args$ln_FishAge_theta_agg
  FishAge_corr_pars_agg <- if(!"FishAge_corr_pars_agg" %in% names(args)) optim_parameters_list$FishAge_corr_pars_agg else args$FishAge_corr_pars_agg
  FishAge_corr_pars <- if(!"FishAge_corr_pars" %in% names(args)) optim_parameters_list$FishAge_corr_pars[,,,,drop = FALSE] else args$FishAge_corr_pars

  # Fishery length compositions
  comp_fishlen_like <- if(!"comp_fishlen_like" %in% names(args)) data$FishLenComps_LikeType else args$comp_fishlen_like
  FishLenComps_Type <- if(!"FishLenComps_Type" %in% names(args)) extend_years(data$FishLenComps_Type, closed_loop_yrs, 1, 'last') else args$FishLenComps_Type
  ISS_FishLenComps <- if(!"ISS_FishLenComps" %in% names(args)) {
    extend_years(replicate(sim_list$n_sims, data$ISS_FishLenComps[,,,,,drop = FALSE] * data$Wt_FishLenComps), closed_loop_yrs, 2, fill = ISS_FishLenComps_fill)
  } else args$ISS_FishLenComps
  ln_FishLen_theta <- if(!"ln_FishLen_theta" %in% names(args)) optim_parameters_list$ln_FishLen_theta[,,,drop = FALSE] else args$ln_FishLen_theta
  ln_FishLen_theta_agg <- if(!"ln_FishLen_theta_agg" %in% names(args)) optim_parameters_list$ln_FishLen_theta_agg else args$ln_FishLen_theta_agg
  FishLen_corr_pars_agg <- if(!"FishLen_corr_pars_agg" %in% names(args)) optim_parameters_list$FishLen_corr_pars_agg else args$FishLen_corr_pars_agg
  FishLen_corr_pars <- if(!"FishLen_corr_pars" %in% names(args)) optim_parameters_list$FishLen_corr_pars[,,,,drop = FALSE] else args$FishLen_corr_pars

  # Discard Fishery age compositions
  comp_fishage_discard_like <- if(!"comp_fishage_discard_like" %in% names(args)) data$FishAgeComps_discard_LikeType else args$comp_fishage_discard_like
  FishAgeComps_discard_Type <- if(!"FishAgeComps_discard_Type" %in% names(args)) extend_years(data$FishAgeComps_discard_Type, closed_loop_yrs, 1, 'last') else args$FishAgeComps_discard_Type
  ISS_FishAgeComps_discard <- if(!"ISS_FishAgeComps_discard" %in% names(args)) {
    extend_years(replicate(sim_list$n_sims, data$ISS_FishAgeComps_discard[,,,,,drop = FALSE] * data$Wt_FishAgeComps_discard), closed_loop_yrs, 2, fill = ISS_FishAgeComps_discard_fill)
  } else args$ISS_FishAgeComps_discard_fill
  ln_FishAge_discard_theta <- if(!"ln_FishAge_discard_theta" %in% names(args)) optim_parameters_list$ln_FishAge_discard_theta[,,,drop = FALSE] else args$ln_FishAge_discard_theta
  ln_FishAge_discard_theta_agg <- if(!"ln_FishAge_discard_theta_agg" %in% names(args)) optim_parameters_list$ln_FishAge_discard_theta_agg else args$ln_FishAge_discard_theta_agg
  FishAge_discard_corr_pars_agg <- if(!"FishAge_discard_corr_pars_agg" %in% names(args)) optim_parameters_list$FishAge_discard_corr_pars_agg else args$FishAge_discard_corr_pars_agg
  FishAge_discard_corr_pars <- if(!"FishAge_discard_corr_pars" %in% names(args)) optim_parameters_list$FishAge_discard_corr_pars[,,,,drop = FALSE] else args$FishAge_discard_corr_pars

  # Discard Fishery length compositions
  comp_fishlen_discard_like <- if(!"comp_fishlen_discard_like" %in% names(args)) data$FishLenComps_discard_LikeType else args$comp_fishlen_discard_like
  FishLenComps_discard_Type <- if(!"FishLenComps_discard_Type" %in% names(args)) extend_years(data$FishLenComps_discard_Type, closed_loop_yrs, 1, 'last') else args$FishLenComps_discard_Type
  ISS_FishLenComps_discard <- if(!"ISS_FishLenComps_discard" %in% names(args)) {
    extend_years(replicate(sim_list$n_sims, data$ISS_FishLenComps_discard[,,,,,drop = FALSE] * data$Wt_FishLenComps_discard), closed_loop_yrs, 2, fill = ISS_FishLenComps_discard_fill)
  } else args$ISS_FishLenComps_discard_fill
  ln_FishLen_discard_theta <- if(!"ln_FishLen_discard_theta" %in% names(args)) optim_parameters_list$ln_FishLen_discard_theta[,,,drop = FALSE] else args$ln_FishLen_discard_theta
  ln_FishLen_discard_theta_agg <- if(!"ln_FishLen_discard_theta_agg" %in% names(args)) optim_parameters_list$ln_FishLen_discard_theta_agg else args$ln_FishLen_discard_theta_agg
  FishLen_discard_corr_pars_agg <- if(!"FishLen_discard_corr_pars_agg" %in% names(args)) optim_parameters_list$FishLen_discard_corr_pars_agg else args$FishLen_discard_corr_pars_agg
  FishLen_discard_corr_pars <- if(!"FishLen_discard_corr_pars" %in% names(args)) optim_parameters_list$FishLen_discard_corr_pars[,,,,drop = FALSE] else args$FishLen_discard_corr_pars

  # Population-specific fishery index SE
  ObsFishIdx_pop_SE <- if(!"ObsFishIdx_pop_SE" %in% names(args)) {
    if(any(data$UseFishIdx_pop == 1)) {
      extend_years(data$ObsFishIdx_pop_SE / sqrt(data$Wt_FishIdx_pop), closed_loop_yrs, 3, fill = FishIdx_SE_pop_fill)
    } else {
      array(0.2, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets))
    }
  } else args$ObsFishIdx_pop_SE

  # Population-specific fishery age compositions
  comp_fishage_pop_like <- if(!"comp_fishage_pop_like" %in% names(args)) data$FishAgeComps_pop_LikeType else args$comp_fishage_pop_like
  FishAgeComps_pop_Type <- if(!"FishAgeComps_pop_Type" %in% names(args)) extend_years(data$FishAgeComps_pop_Type, closed_loop_yrs, 1, 'last') else args$FishAgeComps_pop_Type
  ISS_FishAgeComps_pop <- if(!"ISS_FishAgeComps_pop" %in% names(args)) extend_years(replicate(sim_list$n_sims, data$ISS_FishAgeComps_pop[,,,,,,drop = FALSE] * data$Wt_FishAgeComps_pop), closed_loop_yrs, 3, fill = ISS_FishAgeComps_pop_fill) else args$ISS_FishAgeComps_pop
  ln_FishAge_pop_theta <- if(!"ln_FishAge_pop_theta" %in% names(args)) optim_parameters_list$ln_FishAge_pop_theta[,,,,drop = FALSE] else args$ln_FishAge_pop_theta
  ln_FishAge_pop_theta_agg <- if(!"ln_FishAge_pop_theta_agg" %in% names(args)) optim_parameters_list$ln_FishAge_pop_theta_agg else args$ln_FishAge_pop_theta_agg
  FishAge_pop_corr_pars_agg <- if(!"FishAge_pop_corr_pars_agg" %in% names(args)) optim_parameters_list$FishAge_pop_corr_pars_agg else args$FishAge_pop_corr_pars_agg
  FishAge_pop_corr_pars <- if(!"FishAge_pop_corr_pars" %in% names(args)) optim_parameters_list$FishAge_pop_corr_pars[,,,,,drop = FALSE] else args$FishAge_pop_corr_pars

  # Population-specific fishery length compositions
  comp_fishlen_pop_like <- if(!"comp_fishlen_pop_like" %in% names(args)) data$FishLenComps_pop_LikeType else args$comp_fishlen_pop_like
  FishLenComps_pop_Type <- if(!"FishLenComps_pop_Type" %in% names(args)) extend_years(data$FishLenComps_pop_Type, closed_loop_yrs, 1, 'last') else args$FishLenComps_pop_Type
  ISS_FishLenComps_pop <- if(!"ISS_FishLenComps_pop" %in% names(args)) extend_years(replicate(sim_list$n_sims, data$ISS_FishLenComps_pop[,,,,,,drop = FALSE] * data$Wt_FishLenComps_pop), closed_loop_yrs, 3, fill = ISS_FishLenComps_pop_fill) else args$ISS_FishLenComps_pop
  ln_FishLen_pop_theta <- if(!"ln_FishLen_pop_theta" %in% names(args)) optim_parameters_list$ln_FishLen_pop_theta[,,,,drop = FALSE] else args$ln_FishLen_pop_theta
  ln_FishLen_pop_theta_agg <- if(!"ln_FishLen_pop_theta_agg" %in% names(args)) optim_parameters_list$ln_FishLen_pop_theta_agg else args$ln_FishLen_pop_theta_agg
  FishLen_pop_corr_pars_agg <- if(!"FishLen_pop_corr_pars_agg" %in% names(args)) optim_parameters_list$FishLen_pop_corr_pars_agg else args$FishLen_pop_corr_pars_agg
  FishLen_pop_corr_pars <- if(!"FishLen_pop_corr_pars" %in% names(args)) optim_parameters_list$FishLen_pop_corr_pars[,,,,,drop = FALSE] else args$FishLen_pop_corr_pars

  # Discarded Population-specific fishery age compositions
  comp_fishage_discard_pop_like <- if(!"comp_fishage_discard_pop_like" %in% names(args)) data$FishAgeComps_pop_LikeType else args$comp_fishage_discard_pop_like
  FishAgeComps_discard_pop_Type <- if(!"FishAgeComps_discard_pop_Type" %in% names(args)) extend_years(data$FishAgeComps_discard_pop_Type, closed_loop_yrs, 1, 'last') else args$FishAgeComps_discard_pop_Type
  ISS_FishAgeComps_discard_pop <- if(!"ISS_FishAgeComps_discard_pop" %in% names(args)) extend_years(replicate(sim_list$n_sims, data$ISS_FishAgeComps_discard_pop[,,,,,,drop = FALSE] * data$Wt_FishAgeComps_discard_pop), closed_loop_yrs, 3, fill = ISS_FishAgeComps_discard_pop_fill) else args$ISS_FishAgeComps_discard_pop
  ln_FishAge_discard_pop_theta <- if(!"ln_FishAge_discard_pop_theta" %in% names(args)) optim_parameters_list$ln_FishAge_discard_pop_theta[,,,,drop = FALSE] else args$ln_FishAge_discard_pop_theta
  ln_FishAge_discard_pop_theta_agg <- if(!"ln_FishAge_discard_pop_theta_agg" %in% names(args)) optim_parameters_list$ln_FishAge_discard_pop_theta_agg else args$ln_FishAge_discard_pop_theta_agg
  FishAge_discard_pop_corr_pars_agg <- if(!"FishAge_discard_pop_corr_pars_agg" %in% names(args)) optim_parameters_list$FishAge_discard_pop_corr_pars_agg else args$FishAge_discard_pop_corr_pars_agg
  FishAge_discard_pop_corr_pars <- if(!"FishAge_discard_pop_corr_pars" %in% names(args)) optim_parameters_list$FishAge_discard_pop_corr_pars[,,,,,drop = FALSE] else args$FishAge_discard_pop_corr_pars

  # Discarded Population-specific fishery length compositions
  comp_fishlen_discard_pop_like <- if(!"comp_fishlen_discard_pop_like" %in% names(args)) data$FishLenComps_pop_LikeType else args$comp_fishlen_discard_pop_like
  FishLenComps_discard_pop_Type <- if(!"FishLenComps_discard_pop_Type" %in% names(args)) extend_years(data$FishLenComps_discard_pop_Type, closed_loop_yrs, 1, 'last') else args$FishLenComps_discard_pop_Type
  ISS_FishLenComps_discard_pop <- if(!"ISS_FishLenComps_discard_pop" %in% names(args)) extend_years(replicate(sim_list$n_sims, data$ISS_FishLenComps_discard_pop[,,,,,,drop = FALSE] * data$Wt_FishLenComps_discard_pop), closed_loop_yrs, 3, fill = ISS_FishLenComps_discard_pop_fill) else args$ISS_FishLenComps_discard_pop
  ln_FishLen_discard_pop_theta <- if(!"ln_FishLen_discard_pop_theta" %in% names(args)) optim_parameters_list$ln_FishLen_discard_pop_theta[,,,,drop = FALSE] else args$ln_FishLen_discard_pop_theta
  ln_FishLen_discard_pop_theta_agg <- if(!"ln_FishLen_discard_pop_theta_agg" %in% names(args)) optim_parameters_list$ln_FishLen_discard_pop_theta_agg else args$ln_FishLen_discard_pop_theta_agg
  FishLen_discard_pop_corr_pars_agg <- if(!"FishLen_discard_pop_corr_pars_agg" %in% names(args)) optim_parameters_list$FishLen_discard_pop_corr_pars_agg else args$FishLen_discard_pop_corr_pars_agg
  FishLen_discard_pop_corr_pars <- if(!"FishLen_discard_pop_corr_pars" %in% names(args)) optim_parameters_list$FishLen_discard_pop_corr_pars[,,,,,drop = FALSE] else args$FishLen_discard_pop_corr_pars

  # setup fishery simulation processes
  sim_list <- Setup_Sim_Fishing(
    sim_list = sim_list, # update simulate list
    ln_sigmaC = ln_sigmaC,
    # Age-disaggregated observation error, carried through unweighted: these are
    # keyed by age and fleet and are not scaled by a likelihood weight the way
    # the aggregated sigmas are.
    ln_sigmaCAA = optim_parameters_list$ln_sigmaCAA,
    ln_sigmaDAA = optim_parameters_list$ln_sigmaDAA,
    # the at-age streams carry year on their second margin and no simulation
    # margin, and the closed-loop years keep observing whatever the terminal year
    # observed, the same rule the compositions and selectivities follow
    UseCatchAA = extend_years(data$UseCatchAA, closed_loop_yrs, 2, fill = 'last'),
    UseDiscardAA = extend_years(data$UseDiscardAA, closed_loop_yrs, 2, fill = 'last'),
    use_catch_aa = data$use_catch_aa,
    use_discard_aa = data$use_discard_aa,
    ObsCatchAA_SE = extend_years(data$ObsCatchAA_SE, closed_loop_yrs, 2, fill = 'last'),
    ObsDiscardAA_SE = extend_years(data$ObsDiscardAA_SE, closed_loop_yrs, 2, fill = 'last'),
    CatchAA_Type = extend_years(data$CatchAA_Type, closed_loop_yrs, 1, fill = 'last'),
    DiscardAA_Type = extend_years(data$DiscardAA_Type, closed_loop_yrs, 1, fill = 'last'),
    CatchAA_LikeType = data$CatchAA_LikeType, DiscardAA_LikeType = data$DiscardAA_LikeType,
    CatchAA_sigma_form = data$CatchAA_sigma_form, DiscardAA_sigma_form = data$DiscardAA_sigma_form,
    ln_sigmaC_pop = ln_sigmaC_pop,
    Fmort_input = extend_years(replicate(n = sim_list$n_sims, rep$Fmort[,1:length(data$years),,,drop = FALSE]), n_years = closed_loop_yrs, 2, fill = 'zeros'),
    ln_sigmaD = ln_sigmaD,
    ln_sigmaD_pop = ln_sigmaD_pop,
    dmr_input = extend_years(replicate(n = sim_list$n_sims, rep$dmr[,1:length(data$years),,,drop = FALSE]), n_years = closed_loop_yrs, 2, fill = 'zeros'),
    fish_sel_input = fish_sel_input,
    ret_sel_input = ret_sel_input,
    fish_q_input = fish_q_input,
    ObsFishIdx_SE = ObsFishIdx_SE,
    ObsFishIdx_pop_SE = ObsFishIdx_pop_SE,
    fish_idx_type = data$fish_idx_type,
    init_F_val = rep$init_F,
    catch_units = data$catch_units,
    discard_units = data$discard_units,

    # fishery age composition specifications
    comp_fishage_like = comp_fishage_like,
    FishAgeComps_Type = FishAgeComps_Type,
    ISS_FishAgeComps = ISS_FishAgeComps,
    ln_FishAge_theta = ln_FishAge_theta ,
    ln_FishAge_theta_agg = ln_FishAge_theta_agg,
    FishAge_corr_pars_agg = FishAge_corr_pars_agg,
    FishAge_corr_pars = FishAge_corr_pars,

    # fishery length composition specifications
    comp_fishlen_like = comp_fishlen_like,
    FishLenComps_Type = FishLenComps_Type,
    ISS_FishLenComps =ISS_FishLenComps,
    ln_FishLen_theta = ln_FishLen_theta,
    ln_FishLen_theta_agg = ln_FishLen_theta_agg,
    FishLen_corr_pars_agg = FishLen_corr_pars_agg,
    FishLen_corr_pars = FishLen_corr_pars,

    # population-specific age composition specifications
    comp_fishage_pop_like = comp_fishage_pop_like,
    FishAgeComps_pop_Type = FishAgeComps_pop_Type,
    ISS_FishAgeComps_pop = ISS_FishAgeComps_pop,
    ln_FishAge_pop_theta = ln_FishAge_pop_theta,
    ln_FishAge_pop_theta_agg = ln_FishAge_pop_theta_agg,
    FishAge_pop_corr_pars_agg = FishAge_pop_corr_pars_agg,
    FishAge_pop_corr_pars = FishAge_pop_corr_pars,

    # population-specific length composition specifications
    comp_fishlen_pop_like = comp_fishlen_pop_like,
    FishLenComps_pop_Type = FishLenComps_pop_Type,
    ISS_FishLenComps_pop = ISS_FishLenComps_pop,
    ln_FishLen_pop_theta = ln_FishLen_pop_theta,
    ln_FishLen_pop_theta_agg = ln_FishLen_pop_theta_agg,
    FishLen_pop_corr_pars_agg = FishLen_pop_corr_pars_agg,
    FishLen_pop_corr_pars = FishLen_pop_corr_pars,

    # discarded fishery age composition specifications
    comp_fishage_discard_like = comp_fishage_discard_like,
    FishAgeComps_discard_Type = FishAgeComps_discard_Type,
    ISS_FishAgeComps_discard = ISS_FishAgeComps_discard,
    ln_FishAge_discard_theta = ln_FishAge_discard_theta ,
    ln_FishAge_discard_theta_agg = ln_FishAge_discard_theta_agg,
    FishAge_discard_corr_pars_agg = FishAge_discard_corr_pars_agg,
    FishAge_discard_corr_pars = FishAge_discard_corr_pars,

    # discarded fishery length composition specifications
    comp_fishlen_discard_like = comp_fishlen_discard_like,
    FishLenComps_discard_Type = FishLenComps_discard_Type,
    ISS_FishLenComps_discard =ISS_FishLenComps_discard,
    ln_FishLen_discard_theta = ln_FishLen_discard_theta,
    ln_FishLen_discard_theta_agg = ln_FishLen_discard_theta_agg,
    FishLen_discard_corr_pars_agg = FishLen_discard_corr_pars_agg,
    FishLen_discard_corr_pars = FishLen_discard_corr_pars,

    # discarded population-specific age composition specifications
    comp_fishage_discard_pop_like = comp_fishage_discard_pop_like,
    FishAgeComps_discard_pop_Type = FishAgeComps_discard_pop_Type,
    ISS_FishAgeComps_discard_pop = ISS_FishAgeComps_discard_pop,
    ln_FishAge_discard_pop_theta = ln_FishAge_discard_pop_theta,
    ln_FishAge_discard_pop_theta_agg = ln_FishAge_discard_pop_theta_agg,
    FishAge_discard_pop_corr_pars_agg = FishAge_discard_pop_corr_pars_agg,
    FishAge_discard_pop_corr_pars = FishAge_discard_pop_corr_pars,

    # discarded population-specific length composition specifications
    comp_fishlen_discard_pop_like = comp_fishlen_discard_pop_like,
    FishLenComps_discard_pop_Type = FishLenComps_discard_pop_Type,
    ISS_FishLenComps_discard_pop = ISS_FishLenComps_discard_pop,
    ln_FishLen_discard_pop_theta = ln_FishLen_discard_pop_theta,
    ln_FishLen_discard_pop_theta_agg = ln_FishLen_discard_pop_theta_agg,
    FishLen_discard_pop_corr_pars_agg = FishLen_discard_pop_corr_pars_agg,
    FishLen_discard_pop_corr_pars = FishLen_discard_pop_corr_pars
  )

  # add in ISS F pattern into simulation list
  if(ISS_FishAgeComps_fill == 'F_pattern') sim_list$ISS_FishAgeComps_fill <- "F_pattern"
  if(ISS_FishLenComps_fill == 'F_pattern') sim_list$ISS_FishLenComps_fill <- "F_pattern"
  if(ISS_FishAgeComps_pop_fill == 'F_pattern') sim_list$ISS_FishAgeComps_pop_fill <- "F_pattern"
  if(ISS_FishLenComps_pop_fill == 'F_pattern') sim_list$ISS_FishLenComps_pop_fill <- "F_pattern"
  if(ISS_FishAgeComps_discard_fill == 'F_pattern') sim_list$ISS_FishAgeComps_discard_fill <- "F_pattern"
  if(ISS_FishLenComps_discard_fill == 'F_pattern') sim_list$ISS_FishLenComps_discard_fill <- "F_pattern"
  if(ISS_FishAgeComps_discard_pop_fill == 'F_pattern') sim_list$ISS_FishAgeComps_discard_pop_fill <- "F_pattern"
  if(ISS_FishLenComps_discard_pop_fill == 'F_pattern') sim_list$ISS_FishLenComps_discard_pop_fill <- "F_pattern"

  # Setup Survey Processes --------------------------------------------------
  # Survey selectivity
  srv_sel_input <- if(!"srv_sel_input" %in% names(args)) {
    extend_years(replicate(n = sim_list$n_sims, rep$srv_sel[,,1:length(data$years),,,,,drop = FALSE]), closed_loop_yrs, 3, 'last')
  } else args$srv_sel_input
  # Survey catchability / q
  srv_q_input <- if(!"srv_q_input" %in% names(args)) {
    extend_years(replicate(n = sim_list$n_sims, rep$srv_q[,1:length(data$years),,drop = FALSE]), closed_loop_yrs, 2, 'last')
  } else args$srv_q_input
  # Survey index uncertainty
  ObsSrvIdx_SE <- if(!"ObsSrvIdx_SE" %in% names(args)) {
    extend_years(arr = data$ObsSrvIdx_SE / sqrt(data$Wt_SrvIdx), n_years = closed_loop_yrs, 2, fill = SrvIdx_SE_fill)
  } else args$ObsSrvIdx_SE

  # Survey age compositions
  comp_srvage_like <- if(!"comp_srvage_like" %in% names(args)) data$SrvAgeComps_LikeType else args$comp_srvage_like
  SrvAgeComps_Type <- if(!"SrvAgeComps_Type" %in% names(args)) extend_years(data$SrvAgeComps_Type, closed_loop_yrs, 1, 'last') else args$SrvAgeComps_Type
  ISS_SrvAgeComps <- if(!"ISS_SrvAgeComps" %in% names(args)) {
    extend_years(replicate(sim_list$n_sims, data$ISS_SrvAgeComps[,,,,,drop = FALSE] * data$Wt_SrvAgeComps), closed_loop_yrs, 2, fill = ISS_SrvAgeComps_fill)
  } else args$ISS_SrvAgeComps
  ln_SrvAge_theta <- if(!"ln_SrvAge_theta" %in% names(args)) optim_parameters_list$ln_SrvAge_theta[,,,drop = FALSE] else args$ln_SrvAge_theta
  ln_SrvAge_theta_agg <- if(!"ln_SrvAge_theta_agg" %in% names(args)) optim_parameters_list$ln_SrvAge_theta_agg else args$ln_SrvAge_theta_agg
  SrvAge_corr_pars_agg <- if(!"SrvAge_corr_pars_agg" %in% names(args)) optim_parameters_list$SrvAge_corr_pars_agg else args$SrvAge_corr_pars_agg
  SrvAge_corr_pars <- if(!"SrvAge_corr_pars" %in% names(args)) optim_parameters_list$SrvAge_corr_pars[,,,,drop = FALSE] else args$SrvAge_corr_pars

  # Survey length compositions
  comp_srvlen_like <- if(!"comp_srvlen_like" %in% names(args)) data$SrvLenComps_LikeType else args$comp_srvlen_like
  SrvLenComps_Type <- if(!"SrvLenComps_Type" %in% names(args)) extend_years(data$SrvLenComps_Type, closed_loop_yrs, 1, 'last') else args$SrvLenComps_Type
  ISS_SrvLenComps <- if(!"ISS_SrvLenComps" %in% names(args)) {
    extend_years(replicate(sim_list$n_sims, data$ISS_SrvLenComps[,,,,,drop = FALSE] * data$Wt_SrvLenComps), closed_loop_yrs, 2, fill = ISS_SrvLenComps_fill)
  } else args$ISS_SrvLenComps
  ln_SrvLen_theta <- if(!"ln_SrvLen_theta" %in% names(args)) optim_parameters_list$ln_SrvLen_theta[,,,drop = FALSE] else args$ln_SrvLen_theta
  ln_SrvLen_theta_agg <- if(!"ln_SrvLen_theta_agg" %in% names(args)) optim_parameters_list$ln_SrvLen_theta_agg else args$ln_SrvLen_theta_agg
  SrvLen_corr_pars_agg <- if(!"SrvLen_corr_pars_agg" %in% names(args)) optim_parameters_list$SrvLen_corr_pars_agg else args$SrvLen_corr_pars_agg
  SrvLen_corr_pars <- if(!"SrvLen_corr_pars" %in% names(args)) optim_parameters_list$SrvLen_corr_pars[,,,,drop = FALSE] else args$SrvLen_corr_pars

  # Population-specific survey index SE
  ObsSrvIdx_pop_SE <- if(!"ObsSrvIdx_pop_SE" %in% names(args)) {
    if(any(data$UseSrvIdx_pop == 1)) {
      extend_years(data$ObsSrvIdx_pop_SE / sqrt(data$Wt_SrvIdx_pop), closed_loop_yrs, 3, fill = SrvIdx_SE_pop_fill)
    } else {
      array(0.2, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_srv_fleets))
    }
  } else args$ObsSrvIdx_pop_SE

  # Population-specific survey age compositions
  comp_srvage_pop_like <- if(!"comp_srvage_pop_like" %in% names(args)) data$SrvAgeComps_pop_LikeType else args$comp_srvage_pop_like
  SrvAgeComps_pop_Type <- if(!"SrvAgeComps_pop_Type" %in% names(args)) extend_years(data$SrvAgeComps_pop_Type, closed_loop_yrs, 1, 'last') else args$SrvAgeComps_pop_Type
  ISS_SrvAgeComps_pop <- if(!"ISS_SrvAgeComps_pop" %in% names(args)) {
    if(any(data$UseSrvAgeComps_pop == 1)) {
      extend_years(replicate(sim_list$n_sims, data$ISS_SrvAgeComps_pop[,,,,,,drop = FALSE] * data$Wt_SrvAgeComps_pop), closed_loop_yrs, 3, fill = ISS_SrvAgeComps_pop_fill)
    } else {
      array(100, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims))
    }
  } else args$ISS_SrvAgeComps_pop
  ln_SrvAge_pop_theta <- if(!"ln_SrvAge_pop_theta" %in% names(args)) optim_parameters_list$ln_SrvAge_pop_theta[,,,,drop = FALSE] else args$ln_SrvAge_pop_theta
  ln_SrvAge_pop_theta_agg <- if(!"ln_SrvAge_pop_theta_agg" %in% names(args)) optim_parameters_list$ln_SrvAge_pop_theta_agg else args$ln_SrvAge_pop_theta_agg
  SrvAge_pop_corr_pars_agg <- if(!"SrvAge_pop_corr_pars_agg" %in% names(args)) optim_parameters_list$SrvAge_pop_corr_pars_agg else args$SrvAge_pop_corr_pars_agg
  SrvAge_pop_corr_pars <- if(!"SrvAge_pop_corr_pars" %in% names(args)) optim_parameters_list$SrvAge_pop_corr_pars[,,,,,drop = FALSE] else args$SrvAge_pop_corr_pars

  # Population-specific survey length compositions
  comp_srvlen_pop_like <- if(!"comp_srvlen_pop_like" %in% names(args)) data$SrvLenComps_pop_LikeType else args$comp_srvlen_pop_like
  SrvLenComps_pop_Type <- if(!"SrvLenComps_pop_Type" %in% names(args)) extend_years(data$SrvLenComps_pop_Type, closed_loop_yrs, 1, 'last') else args$SrvLenComps_pop_Type
  ISS_SrvLenComps_pop <- if(!"ISS_SrvLenComps_pop" %in% names(args)) {
    if(any(data$UseSrvLenComps_pop == 1)) {
      extend_years(replicate(sim_list$n_sims, data$ISS_SrvLenComps_pop[,,,,,,drop = FALSE] * data$Wt_SrvLenComps_pop), closed_loop_yrs, 3, fill = ISS_SrvLenComps_pop_fill)
    } else {
      array(100, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims))
    }
  } else args$ISS_SrvLenComps_pop
  ln_SrvLen_pop_theta <- if(!"ln_SrvLen_pop_theta" %in% names(args)) optim_parameters_list$ln_SrvLen_pop_theta[,,,,drop = FALSE] else args$ln_SrvLen_pop_theta
  ln_SrvLen_pop_theta_agg <- if(!"ln_SrvLen_pop_theta_agg" %in% names(args)) optim_parameters_list$ln_SrvLen_pop_theta_agg else args$ln_SrvLen_pop_theta_agg
  SrvLen_pop_corr_pars_agg <- if(!"SrvLen_pop_corr_pars_agg" %in% names(args)) optim_parameters_list$SrvLen_pop_corr_pars_agg else args$SrvLen_pop_corr_pars_agg
  SrvLen_pop_corr_pars <- if(!"SrvLen_pop_corr_pars" %in% names(args)) optim_parameters_list$SrvLen_pop_corr_pars[,,,,,drop = FALSE] else args$SrvLen_pop_corr_pars

  # setup survey simulation processes
  sim_list <- Setup_Sim_Survey(
    sim_list = sim_list,
    srv_sel_input = srv_sel_input,
    srv_q_input = srv_q_input,
    ObsSrvIdx_SE = ObsSrvIdx_SE,
    ln_sigmaSrvIdxAA = optim_parameters_list$ln_sigmaSrvIdxAA,
    UseSrvIdxAA = extend_years(data$UseSrvIdxAA, closed_loop_yrs, 2, fill = 'last'),
    use_srv_idx_aa = data$use_srv_idx_aa,
    ObsSrvIdxAA_SE = extend_years(data$ObsSrvIdxAA_SE, closed_loop_yrs, 2, fill = 'last'),
    SrvIdxAA_Type = extend_years(data$SrvIdxAA_Type, closed_loop_yrs, 1, fill = 'last'),
    SrvIdxAA_LikeType = data$SrvIdxAA_LikeType, SrvIdxAA_sigma_form = data$SrvIdxAA_sigma_form,
    ObsSrvIdx_pop_SE = ObsSrvIdx_pop_SE,
    srv_idx_type = data$srv_idx_type,
    t_srv = data$t_srv,

    # Survey age composition specifications
    comp_srvage_like = comp_srvage_like,
    SrvAgeComps_Type = SrvAgeComps_Type,
    ISS_SrvAgeComps = ISS_SrvAgeComps,
    ln_SrvAge_theta = ln_SrvAge_theta,
    ln_SrvAge_theta_agg = ln_SrvAge_theta_agg,
    SrvAge_corr_pars_agg = SrvAge_corr_pars_agg,
    SrvAge_corr_pars = SrvAge_corr_pars,

    # Survey length composition specifications
    comp_srvlen_like = comp_srvlen_like,
    SrvLenComps_Type = SrvLenComps_Type,
    ISS_SrvLenComps = ISS_SrvLenComps,
    ln_SrvLen_theta = ln_SrvLen_theta,
    ln_SrvLen_theta_agg = ln_SrvLen_theta_agg,
    SrvLen_corr_pars_agg = SrvLen_corr_pars_agg,
    SrvLen_corr_pars = SrvLen_corr_pars,

    # population-specific age composition specifications
    comp_srvage_pop_like = comp_srvage_pop_like,
    SrvAgeComps_pop_Type = SrvAgeComps_pop_Type,
    ISS_SrvAgeComps_pop = ISS_SrvAgeComps_pop,
    ln_SrvAge_pop_theta = ln_SrvAge_pop_theta,
    ln_SrvAge_pop_theta_agg = ln_SrvAge_pop_theta_agg,
    SrvAge_pop_corr_pars_agg = SrvAge_pop_corr_pars_agg,
    SrvAge_pop_corr_pars = SrvAge_pop_corr_pars,

    # population-specific length composition specifications
    comp_srvlen_pop_like = comp_srvlen_pop_like,
    SrvLenComps_pop_Type = SrvLenComps_pop_Type,
    ISS_SrvLenComps_pop = ISS_SrvLenComps_pop,
    ln_SrvLen_pop_theta = ln_SrvLen_pop_theta,
    ln_SrvLen_pop_theta_agg = ln_SrvLen_pop_theta_agg,
    SrvLen_pop_corr_pars_agg = SrvLen_pop_corr_pars_agg,
    SrvLen_pop_corr_pars = SrvLen_pop_corr_pars
  )

  # Setup Biological Dynamics -----------------------------------------------
  natmort_input <- if(!"natmort_input" %in% names(args)) {
    extend_years(replicate(n = sim_list$n_sims, rep$natmort[,,1:length(data$years),,,drop = FALSE]), closed_loop_yrs, 3, 'last')
  } else args$natmort_input
  # biologicals the growth module derives come from the report, the rest from the data
  WAA_input <- if(!"WAA_input" %in% names(args)) {
    extend_years(replicate(n = sim_list$n_sims, (if(is.null(rep$WAA)) data$WAA else rep$WAA)[,,1:length(data$years),,,,drop = FALSE]), closed_loop_yrs, 3, 'last')
  } else args$WAA_input
  WAA_fish_input <- if(!"WAA_fish_input" %in% names(args)) {
    extend_years(replicate(n = sim_list$n_sims, (if(is.null(rep$WAA_fish)) data$WAA_fish else rep$WAA_fish)[,,1:length(data$years),,,,,drop = FALSE]), closed_loop_yrs, 3, 'last')
  } else args$WAA_fish_input
  WAA_srv_input <- if(!"WAA_srv_input" %in% names(args)) {
    extend_years(replicate(n = sim_list$n_sims, (if(is.null(rep$WAA_srv)) data$WAA_srv else rep$WAA_srv)[,,1:length(data$years),,,,,drop = FALSE]), closed_loop_yrs, 3, 'last')
  } else args$WAA_srv_input
  MatAA_input <- if(!"MatAA_input" %in% names(args)) {
    extend_years(replicate(n = sim_list$n_sims, data$MatAA[,,1:length(data$years),,,,drop = FALSE]), closed_loop_yrs, 3, 'last')
  } else args$MatAA_input
  AgeingError_input <- if(!"AgeingError_input" %in% names(args)) {
    extend_years(replicate(n = sim_list$n_sims, data$AgeingError[1:length(data$years),,,drop = FALSE]), closed_loop_yrs, 1, 'last')
  } else args$AgeingError_input
  AgeingError_fish_input <- if(!"AgeingError_fish_input" %in% names(args)) {
    if(is.null(data$AgeingError_fish)) NULL else extend_years(replicate(n = sim_list$n_sims, data$AgeingError_fish[1:length(data$years),,,,drop = FALSE]), closed_loop_yrs, 1, 'last')
  } else args$AgeingError_fish_input
  AgeingError_srv_input <- if(!"AgeingError_srv_input" %in% names(args)) {
    if(is.null(data$AgeingError_srv)) NULL else extend_years(replicate(n = sim_list$n_sims, data$AgeingError_srv[1:length(data$years),,,,drop = FALSE]), closed_loop_yrs, 1, 'last')
  } else args$AgeingError_srv_input
  SizeAgeTrans_input <- if(!"SizeAgeTrans_input" %in% names(args)) {
    if(data$fit_lengths == 0) array(NA, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_lens, sim_list$n_ages, sim_list$n_sexes))
    if(data$fit_lengths == 1 && !is.null(data$SizeAgeTrans) && !all(is.na(data$SizeAgeTrans))) extend_years(replicate(n = sim_list$n_sims, data$SizeAgeTrans[,,1:length(data$years),,,,,drop = FALSE]), closed_loop_yrs, 3, 'last')
  } else args$SizeAgeTrans_input

  # setup biologicals
  sim_list <- Setup_Sim_Biologicals(
    sim_list = sim_list,
    natmort_input = natmort_input,
    WAA_input = WAA_input,
    WAA_fish_input = WAA_fish_input,
    WAA_srv_input = WAA_srv_input,
    MatAA_input = MatAA_input,
    AgeingError_input = AgeingError_input,
    AgeingError_fish_input = AgeingError_fish_input,
    AgeingError_srv_input = AgeingError_srv_input,
    SizeAgeTrans_input = SizeAgeTrans_input,
    # keys per fleet from the growth module, each at its fleet's own timing
    SizeAgeTrans_fish_input = if(is.null(rep$SizeAgeTrans_fish)) NULL else extend_years(replicate(n = sim_list$n_sims, rep$SizeAgeTrans_fish[,,1:length(data$years),,,,,,drop = FALSE]), closed_loop_yrs, 3, 'last'),
    SizeAgeTrans_srv_input = if(is.null(rep$SizeAgeTrans_srv)) NULL else extend_years(replicate(n = sim_list$n_sims, rep$SizeAgeTrans_srv[,,1:length(data$years),,,,,,drop = FALSE]), closed_loop_yrs, 3, 'last')
  )

  # Setup Recruitment Processes ---------------------------------------------
  h_input <- if(!"h_input" %in% names(args)) {
    replicate(n = sim_list$n_sims, array(rep$h_trans, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs)))
  } else args$h_input
  R0_input <- if(!"R0_input" %in% names(args)) {
    R0_r = array(0, dim = c(sim_list$n_pop, sim_list$n_regions)) # container
    for(p in 1:sim_list$n_pop) R0_r[p,] = rep$R0 [p] * rep$rec_region_prop[p,]
    replicate(n = sim_list$n_sims, expr = array(R0_r, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs)))
  } else args$R0_input
  rinit_input <- if(!"rinit_input" %in% names(args)) {
    rinit_r = array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_sims)) # container
    for(p in 1:sim_list$n_pop) for(s in 1:sim_list$n_sims) rinit_r[p,,s] <- rep$rinit[p] * rep$rec_region_prop[p,]
    rinit_r
  } else args$rinit_input
  sexratio_input <- if(!"sexratio_input" %in% names(args)) {
    extend_years(replicate(n = sim_list$n_sims, expr = rep$sexratio[,,1:length(data$years),,drop = FALSE]), closed_loop_yrs, 3, 'last')
  } else args$sexratio_input
  ln_sigmaR <- if(!"ln_sigmaR" %in% names(args)) optim_parameters_list$ln_sigmaR else args$ln_sigmaR
  stray_rate_input <- if(!"stray_rate_input" %in% names(args)) {
    extend_years(replicate(n = sim_list$n_sims, expr = data$stray_rate[,1:length(data$years),drop = FALSE]), closed_loop_yrs, 2, 'last')
  } else args$stray_rate_input
  Rec_input <- if(!"Rec_input" %in% names(args)) {
    replicate(n = sim_list$n_sims, expr = rep$Rec[,,1:length(data$years),drop = FALSE])
  } else args$Rec_input
  ln_InitDevs_input <- if(!"ln_InitDevs_input" %in% names(args)) {
    replicate(sim_list$n_sims, optim_parameters_list$ln_InitDevs)
  } else args$ln_InitDevs_input

  # recruitment options
  rec_dd <- if(!"rec_dd" %in% names(args)) data$rec_dd else args$rec_dd
  init_dd <- if(!"init_dd" %in% names(args)) data$rec_dd else args$init_dd
  rec_lag <- if(!"rec_lag" %in% names(args)) data$rec_lag else args$rec_lag
  recruitment_opt <- if(!"recruitment_opt" %in% names(args)) data$rec_model else args$recruitment_opt
  rec_seas_prop_input <- if(!"rec_seas_prop_input" %in% names(args)) array(replicate(sim_list$n_sims, rep$rec_seas_prop), dim = c(dim(rep$rec_seas_prop), sim_list$n_sims))
  else args$rec_seas_prop_input

  # setup recruitment simulation
  sim_list <- Setup_Sim_Rec(
    sim_list = sim_list,
    spawn_seas = data$spawn_seas,
    do_recruits_move = data$do_recruits_move,
    t_spawn = data$t_spawn,
    init_age_strc = data$init_age_strc,
    h_input = h_input,
    R0_input = R0_input,
    rinit_input = rinit_input,
    use_rinit = data$use_rinit,
    sexratio_input = sexratio_input,
    ln_sigmaR = ln_sigmaR,
    Rec_input = Rec_input,
    ln_InitDevs_input = ln_InitDevs_input,
    recruitment_opt = recruitment_opt,
    rec_dd = rec_dd,
    init_dd = init_dd,
    rec_lag = rec_lag,
    stray_rate_input = stray_rate_input,
    rec_seas_prop_input = rec_seas_prop_input
  )

  # Setup Tagging -----------------------------------------------------------
  n_tags_rel_input <- if(!"n_tags_rel_input" %in% names(args)) {
    if(!is.na(sum(data$conv_tagged_fish))) apply(data$conv_tagged_fish, 1, sum) else NA
  } else args$n_tags_rel_input
  n_tags <- if(!"n_tags" %in% names(args)) NULL else args$n_tags
  conv_tag_release_indicator <- if(!"conv_tag_release_indicator" %in% names(args)) {
    if(exists("conv_tag_release_indicator", data)) data$conv_tag_release_indicator else NA
  } else args$conv_tag_release_indicator
  ln_init_conv_tag_mort <- if(!"ln_init_conv_tag_mort" %in% names(args)) optim_parameters_list$ln_init_conv_tag_mort else args$ln_init_conv_tag_mort
  ln_conv_tag_shed <- if(!"ln_conv_tag_shed" %in% names(args)) optim_parameters_list$ln_conv_tag_shed else args$ln_conv_tag_shed
  conv_tag_fish_reporting_input <- if(!"conv_tag_fish_reporting_input" %in% names(args)) {
    if(is.null(rep$conv_tag_fish_reporting)) NULL else extend_years(replicate(n = sim_list$n_sims, rep$conv_tag_fish_reporting), closed_loop_yrs, 2, 'last')
  } else args$conv_tag_fish_reporting_input
  use_conv_fish_tagging <- if(!"use_conv_fish_tagging" %in% names(args)) data$use_conv_fish_tagging else args$use_conv_fish_tagging
  tag_selex <- if(!"tag_selex" %in% names(args)) data$tag_selex else args$tag_selex
  tag_natmort <- if(!"tag_natmort" %in% names(args)) data$tag_natmort else args$tag_natmort
  conv_fish_tag_like <- if(!"conv_fish_tag_like" %in% names(args)) data$conv_fish_tag_like else args$conv_fish_tag_like
  ln_conv_fish_tag_theta <- if(!"ln_conv_fish_tag_theta" %in% names(args)) parameters$ln_conv_fish_tag_theta else args$ln_conv_fish_tag_theta

  # setup tagging simulation
  if(!is.null(n_tags)) sim_list$n_tags_rel_input <- NULL # set release input to NULL if n_tags is specified.
  sim_list <- Setup_Sim_Tagging(
    sim_list = sim_list,
    conv_tag_max_liberty = data$conv_tag_max_liberty,
    conv_tag_t_tagging = data$conv_tag_t_tagging,
    n_tags = n_tags,
    n_tags_rel_input = n_tags_rel_input * data$Wt_Tagging,
    conv_tag_release_indicator = conv_tag_release_indicator,
    conv_tag_release_platform = data$conv_tag_release_platform,
    ln_init_conv_tag_mort = ln_init_conv_tag_mort,
    ln_conv_tag_shed = ln_conv_tag_shed,
    conv_fish_tag_attr = data$conv_fish_tag_attr,
    conv_tag_fish_reporting_input = conv_tag_fish_reporting_input,
    use_conv_fish_tagging = use_conv_fish_tagging,
    conv_fish_tag_like = conv_fish_tag_like,
    ln_conv_fish_tag_theta = ln_conv_fish_tag_theta
  )

  # Movement ----------------------------------------------------------------
  Movement <- if(!"Movement" %in% names(args)) extend_years(replicate(n = sim_list$n_sims, rep$Movement[,,,1:length(data$years),,,,drop = FALSE]), closed_loop_yrs, 4, 'last') else args$Movement
  sim_list$Movement <- Movement
  sgl_seas_spawning_movement <- if(!"sgl_seas_spawning_movement" %in% names(args)) extend_years(replicate(n = sim_list$n_sims, data$sgl_seas_spawning_movement[,,,1:length(data$years),,,drop = FALSE]), closed_loop_yrs, 4, 'last') else args$sgl_seas_spawning_movement
  sim_list$sgl_seas_spawning_movement <- sgl_seas_spawning_movement
  sim_list$move_timing <- if(!"move_timing" %in% names(args)) {
    if(is.null(data$move_timing)) 0 else data$move_timing
  } else args$move_timing
  sim_list$expm_nsub <- if(!"move_expm_nsub" %in% names(args)) {
    if(is.null(data$move_expm_nsub)) 0 else data$move_expm_nsub
  } else args$move_expm_nsub
  sim_list$Mrate <- if("Mrate" %in% names(args)) args$Mrate else {
    if(sim_list$move_timing != 2 || is.null(rep$Mrate)) NULL
    else extend_years(replicate(n = sim_list$n_sims, rep$Mrate[,,,1:length(data$years),,,,drop = FALSE]), closed_loop_yrs, 4, 'last')
  }

  # State-space numbers at age ----------------------------------------------
  sim_list$NAA_re <- if("NAA_re" %in% names(args)) args$NAA_re else {
    if(is.null(data$NAA_re)) 0 else data$NAA_re
  }
  state_on <- isTRUE(sim_list$NAA_re > 0)

  # Extend sigma blocks for years etc
  sim_list$sigmaNAA <- if("sigmaNAA" %in% names(args)) args$sigmaNAA else {
    if(!state_on || is.null(optim_parameters_list$ln_sigmaNAA)) 0 else {
      blk <- data$naa_sigma_blocks[,,1:length(data$years),,,drop = FALSE]
      extend_years(array(exp(optim_parameters_list$ln_sigmaNAA)[as.vector(blk)], dim = dim(blk)),
                   closed_loop_yrs, 3, 'last')
    }
  }

  # The simulator takes one correlation per margin, so shared parameters come across as they are
  # and per-cell ones are averaged before transforming.
  sim_list$naa_rho <- if("naa_rho" %in% names(args)) args$naa_rho else {
    pe <- optim_parameters_list$NAA_pe_pars
    if(!state_on || is.null(pe)) c(age = 0, year = 0, cohort = 0)
    else c(age = rho_trans(mean(pe[,,1,])), year = rho_trans(mean(pe[,,2,])),
           cohort = rho_trans(mean(pe[,,3,])))
  }

  for(nm in c("NAA_re_pop", "NAA_re_region", "NAA_re_sex")) {
    sim_list[[nm]] <- if(nm %in% names(args)) args[[nm]] else {
      if(!state_on || is.null(data[[nm]])) 0 else data[[nm]]
    }
  } # end nm loop

  # Unconstrained parameters become the correlations themselves, in the lower-triangle order the
  # simulator's unstructured factor reads them in.
  corr_from <- function(pars, n) {
    if(is.null(pars) || n < 2) return(0)
    C <- build_us_corr(as.vector(pars), n)
    C[lower.tri(C)]
  }
  sim_list$naa_pop_corr <- if("naa_pop_corr" %in% names(args)) args$naa_pop_corr else {
    if(!state_on || !isTRUE(sim_list$NAA_re_pop == 1)) 0
    else corr_from(optim_parameters_list$NAA_pop_corr_pars, sim_list$n_pop)
  }
  sim_list$naa_region_corr <- if("naa_region_corr" %in% names(args)) args$naa_region_corr else {
    if(!state_on || !isTRUE(sim_list$NAA_re_region == 1)) 0
    else corr_from(optim_parameters_list$NAA_region_corr_pars, sim_list$n_regions)
  }
  sim_list$naa_sex_corr <- if("naa_sex_corr" %in% names(args)) args$naa_sex_corr else {
    if(!state_on || !isTRUE(sim_list$NAA_re_sex == 1)) 0
    else corr_from(optim_parameters_list$NAA_sex_corr_pars, sim_list$n_sexes)
  }

  # The active ages carry over unchanged, and the active years run on through the projection: the
  # operating model should keep generating process error in the future, not stop at the data.
  sim_list$naa_re_ages <- if("naa_re_ages" %in% names(args)) args$naa_re_ages else {
    if(!state_on || is.null(data$naa_re_ages)) integer(0) else data$naa_re_ages
  }
  sim_list$naa_re_yrs <- if("naa_re_yrs" %in% names(args)) args$naa_re_yrs else {
    if(!state_on || is.null(data$naa_re_yrs) || !length(data$naa_re_yrs)) integer(0)
    else seq(min(data$naa_re_yrs), sim_list$n_yrs)
  }

  return(sim_list)
}

#' Get Closed Loop Reference Points
#'
#' Computes fishery and biological reference points either using "true" simulated values
#' from the operating model or using assessment-derived data and report objects. Supports
#' single-region and multi-region reference points.
#'
#' @param use_true_values Logical. If TRUE, uses values from the simulation environment
#'   (`sim_env`) for calculating reference points. If FALSE, uses `asmt_data` and `asmt_rep`.
#' @param asmt_data Optional list. Assessment data object (from RTMB) if not using true values.
#' @param asmt_rep Optional list. Assessment report object (from RTMB) if not using true values.
#' @param y Integer. Number of years to include in calculations (usually the last year of the assessment or simulation).
#' @param sim Integer. Index of the simulation replicate in `sim_env`.
#' @param reference_points_opt List. Options for reference point calculations:
#'   \describe{
#'     \item{n_avg_yrs}{Number of years to average over demographic rates. Default is 1.}
#'     \item{SPR_x}{Target SPR fraction for reference point calculations. Default is 0.4.}
#'     \item{calc_rec_st_yr}{Year to start calculating mean recruitment. Default is 1.}
#'     \item{rec_age}{Age at recruitment. Default is 1.}
#'     \item{type}{Reference point type: "single_region" or "multi_region". Default is "single_region".}
#'     \item{what}{Method for reference point calculation. Options include "SPR", "MSY",
#'       "independent_SPR", "independent_MSY", "global_SPR", "global_MSY". Default is "SPR".}
#'     \item{is_discard_fleet}{Integer vector \code{[n_fish_fleets]}. Indicator for
#'       fleets whose catch should be excluded from landed yield in reference point
#'       calculations (0 = landing fleet, 1 = discard-only fleet). These fleets still
#'       contribute to total fishing mortality and population dynamics. Default is
#'       \code{NULL}, which treats all fleets as landing fleets.}
#'   }
#' @param sim_env Simulation environment
#' @param n_proj_yrs Number of projection years
#' @param t_spawn Spawn timing within a given season / year, default uses sim_env true values
#'
#' @return A list with elements:
#'   \describe{
#'     \item{f_ref_pt}{Array of fishing reference points by region and projection year.}
#'     \item{b_ref_pt}{Array of biological reference points by region and projection year.}
#'     \item{virgin_b_ref_pt}{Array of unfished biological reference points by region and projection year.}
#'     \item{pop_b_ref_pt}{Array of population-level biological reference points by population, region, and projection year.}
#'     \item{virgin_pop_b_ref_pt}{Array of unfished population-level biological reference points by population, region, and projection year.}
#'   }
#'
#' @export get_closed_loop_reference_points
#' @family Closed Loop Simulations
get_closed_loop_reference_points <- function(use_true_values,
                                             sim_env,
                                             asmt_data = NULL,
                                             asmt_rep = NULL,
                                             t_spawn = sim_env$t_spawn,
                                             y,
                                             sim,
                                             reference_points_opt = list(
                                               n_avg_yrs = 1,
                                               SPR_x = 0.4,
                                               calc_rec_st_yr = 1,
                                               rec_age = 1,
                                               type = 'single_region',
                                               what = "SPR",
                                               is_discard_fleet = NULL
                                             ),
                                             n_proj_yrs
                                             ) {

  if(use_true_values) {

    # Build data and report objects to feed into reference points
    data_obj <- list(
      ages = 1:sim_env$n_ages,
      years = 1:y,
      n_pop = sim_env$n_pop,
      natal_region = sim_env$natal_region,
      n_seas = sim_env$n_seas,
      seasdur = sim_env$seasdur,
      spawn_seas = sim_env$spawn_seas,
      n_fish_fleets = sim_env$n_fish_fleets,
      n_regions = sim_env$n_regions,
      WAA = array(sim_env$WAA[,,1:y, ,, , sim], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_ages, sim_env$n_sexes)),
      MatAA = array(sim_env$MatAA[,, 1:y,, , , sim], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_ages, sim_env$n_sexes)),
      do_recruits_move = sim_env$do_recruits_move,
      move_timing = if(is.null(sim_env$move_timing)) 0 else sim_env$move_timing,
      move_expm_nsub = if(is.null(sim_env$expm_nsub)) 0 else sim_env$expm_nsub
    )

    # Build rep list if not doing assessment (using truth)
    rep_obj <- list(
      Fmort = array(sim_env$Fmort[, 1:y,, , sim], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets)),
      dmr = array(sim_env$dmr[, 1:y,, , sim], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets)),
      fish_sel = array(sim_env$fish_sel[,, 1:y,, , , , sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_ages, sim_env$n_sexes, sim_env$n_fish_fleets)),
      ret_sel = array(sim_env$ret_sel[,, 1:y,, , , , sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_ages, sim_env$n_sexes, sim_env$n_fish_fleets)),
      natmort = array(sim_env$natmort[,, 1:y, , , sim], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_ages, sim_env$n_sexes)),
      h_trans = array(sim_env$h[,, y, sim], dim = c(sim_env$n_pop, sim_env$n_regions)),
      R0 = apply(sim_env$R0[,, y, sim, drop = FALSE], 1, sum),
      stray_rate = array(sim_env$stray_rate[,1:y,sim], dim = c(sim_env$n_pop, length(1:y))),
      rec_seas_prop = array(sim_env$rec_seas_prop[,,sim], dim = c(sim_env$n_pop, sim_env$n_seas)),
      # Get_Reference_Points averages the generator on the rate scale under continuous
      # movement (mean(expm(Q)) != expm(mean(Q))), so Mrate must reach rep as well
      Mrate = if(is.null(sim_env$Mrate)) NULL else
        array(sim_env$Mrate[,,, 1:y,,,, sim], dim = c(sim_env$n_pop, sim_env$n_regions, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_ages, sim_env$n_sexes)),
      rec_region_prop = {
        R0_slice <- sim_env$R0[,, y, sim, drop = FALSE]
        row_sums <- rowSums(R0_slice)
        R0_prop <- array(R0_slice / row_sums, dim = c(sim_env$n_pop, sim_env$n_regions))
      },
      Rec = array(sim_env$Rec[, , 1:y, sim], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y))),
      Movement = array(sim_env$Movement[, , , 1:y, , , , sim],  dim = c(sim_env$n_pop, sim_env$n_regions, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_ages, sim_env$n_sexes)),
      sgl_seas_spawning_movement = array(sim_env$sgl_seas_spawning_movement[, , , 1:y, , , sim],  dim = c(sim_env$n_pop, sim_env$n_regions, sim_env$n_regions, length(1:y), sim_env$n_ages, sim_env$n_sexes)),
      move_timing = if(is.null(sim_env$move_timing)) 0 else sim_env$move_timing,
      Mrate = if(is.null(sim_env$Mrate)) NULL else
        array(sim_env$Mrate[, , , 1:y, , , , sim], dim = c(sim_env$n_pop, sim_env$n_regions, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_ages, sim_env$n_sexes))
    )

    # get sex ratio
    tmp_sex_ratio_f <- if(sim_env$n_sexes == 1) array(0.5, dim = c(sim_env$n_pop, sim_env$n_regions)) else sim_env$sexratio[,,y,1,sim]

  } else {
    data_obj <- asmt_data
    rep_obj <- asmt_rep
    tmp_sex_ratio_f <- if(data_obj$n_sexes == 1) array(0.5, dim = c(sim_env$n_pop, sim_env$n_regions)) else rep_obj$sexratio[,,y,1]
  }

  # dealing with whether there are any discard fleets
  reference_points_opt$is_discard_fleet <- if(!is.null(reference_points_opt$is_discard_fleet)) {
    reference_points_opt$is_discard_fleet
  } else {
    rep(0, data_obj$n_fish_fleets)
  }

  # get reference points based on true values
  reference_points <- Get_Reference_Points(data = data_obj,
                                           rep = rep_obj,
                                           SPR_x = reference_points_opt$SPR_x,
                                           t_spawn = t_spawn,
                                           sex_ratio_f = tmp_sex_ratio_f,
                                           calc_rec_st_yr = reference_points_opt$calc_rec_st_yr,
                                           rec_age = reference_points_opt$rec_age,
                                           type = reference_points_opt$type,
                                           what = reference_points_opt$what,
                                           n_avg_yrs = reference_points_opt$n_avg_yrs,
                                           is_discard_fleet = reference_points_opt$is_discard_fleet
  )

  # extract fishery and biological reference points
  f_ref_pt <- array(reference_points$f_ref_pt, dim = c(data_obj$n_regions, n_proj_yrs)) # fishery reference points
  b_ref_pt <- array(reference_points$b_ref_pt, dim = c(data_obj$n_pop, data_obj$n_regions, n_proj_yrs)) # biological reference points
  virgin_b_ref_pt <- array(reference_points$virgin_b_ref_pt, dim = c(data_obj$n_pop, data_obj$n_regions, n_proj_yrs)) # biological reference points
  pop_b_ref_pt <- array(reference_points$pop_b_ref_pt, dim = c(data_obj$n_pop, data_obj$n_regions, n_proj_yrs)) # biological reference points
  virgin_pop_b_ref_pt <- array(reference_points$virgin_pop_b_ref_pt, dim = c(data_obj$n_pop, data_obj$n_regions, n_proj_yrs)) # biological reference points

  return(list(f_ref_pt = f_ref_pt, b_ref_pt = b_ref_pt, virgin_b_ref_pt = virgin_b_ref_pt, pop_b_ref_pt = pop_b_ref_pt, virgin_pop_b_ref_pt = virgin_pop_b_ref_pt))
}



#' Convert a target catch to fishing mortality for a single fleet via bisection
#'
#' Uses interval bisection to find the scalar fishing mortality \eqn{F} that
#' produces a predicted catch (Baranov catch equation, summed over ages and
#' sexes in biomass) equal to \code{catch}. Retained and discarded catch
#' components are handled separately via retention selectivity and discard
#' mortality. Intended for closed-loop MSE harvest control rules where a TAC
#' must be translated into an \eqn{F} for the operating model.
#'
#' @param f_guess Numeric. Initial \eqn{F} guess (not used directly by the
#'   bisection algorithm but retained for API consistency).
#' @param catch Numeric. Target catch in biomass units.
#' @param NAA Numeric matrix \code{[n_ages x n_sexes]}. Numbers-at-age at
#'   the start of the time step.
#' @param WAA Numeric matrix \code{[n_ages x n_sexes]}. Weight-at-age.
#' @param natmort Numeric matrix \code{[n_ages x n_sexes]}. Instantaneous
#'   natural mortality rate.
#' @param fish_sel Numeric matrix \code{[n_ages x n_sexes]}. Fishery
#'   selectivity (scaled to a maximum of 1).
#' @param ret_sel Numeric matrix \code{[n_ages x n_sexes]}. Retention
#'   selectivity, i.e. the proportion of selected catch that is retained
#'   (vs. discarded). Default is an array of 1s matching \code{dim(fish_sel)}
#'   (full retention).
#' @param dmr Numeric scalar. Discard mortality rate, i.e. the fraction of
#'   discarded fish that die. Default \code{0} (all discards survive).
#' @param n.iter Integer. Number of bisection iterations. Default \code{20};
#'   approximately \eqn{\log_2((ub - lb) / \epsilon)} iterations are required
#'   for tolerance \eqn{\epsilon}.
#' @param lb Numeric. Lower bound of the \eqn{F} search interval.
#'   Default \code{0}.
#' @param ub Numeric. Upper bound of the \eqn{F} search interval.
#'   Default \code{2}.
#'
#' @return Scalar numeric. The \eqn{F} value at the final bisection midpoint
#'   that most closely produces \code{catch}.
#'
#' @export catch_to_F_singlefleet
#' @family Closed Loop Simulations
catch_to_F_singlefleet <- function(f_guess,
                                   catch,
                                   NAA,
                                   WAA,
                                   natmort,
                                   fish_sel,
                                   ret_sel = {
                                     # fill dimensions of fish_sel with 1s
                                     tmp = fish_sel[]
                                     tmp[] <- 1
                                     tmp
                                   },
                                   dmr = 0,
                                   n.iter = 20,
                                   lb = 0,
                                   ub = 2) {

  range <- vector(length=2) # F range
  range[1] <- lb # Lower bound
  range[2] <- ub # Upper bound

  for(i in 1:n.iter) {

    # Get midpoint of range
    midpoint <- mean(range)

    # Caclulate baranov's
    ret_FAA <- (midpoint * fish_sel * ret_sel)
    disc_FAA <- (midpoint * fish_sel * (1 - ret_sel) * dmr)
    ZAA <- ret_FAA + natmort + disc_FAA
    pred_catch <- sum((ret_FAA / ZAA * NAA * (1 - exp(-ZAA))) * WAA)

    if(pred_catch < catch) {
      range[1] <- midpoint
      range[2] <- range[2]
    }else {
      range[1] <- range[1]
      range[2] <- midpoint
    }

  } # end i loop

  return(midpoint)
}

#' Convert target catches to fishing mortality rates for multiple fleets
#'
#' Solves for the vector of fleet-specific fishing mortality rates
#' \eqn{\mathbf{F} = (F_1, \ldots, F_k)} that simultaneously satisfy the
#' Baranov catch equations for all fleets, given a vector of target catches.
#' Total mortality at age accounts for retained and discard mortality
#' contributions from all fleets:
#' \eqn{Z_a = M_a + \sum_f (F_f \, s_{a,f} \, r_{a,f} + F_f \, s_{a,f} \, (1 - r_{a,f}) \, d_f)},
#' where \eqn{r_{a,f}} is retention selectivity and \eqn{d_f} is fleet-specific
#' discard mortality. The system of equations is solved via
#' \code{\link[nleqslv]{nleqslv}}. Intended for closed-loop MSE harvest
#' control rules with multiple interacting fishery fleets.
#'
#' @param target_catch Numeric vector \code{[n_fleets]}. Target catch in
#'   biomass units for each fleet.
#' @param NAA Numeric matrix \code{[n_ages x n_sexes]}. Numbers-at-age.
#' @param WAA Numeric matrix \code{[n_ages x n_sexes]}. Weight-at-age.
#' @param natmort Numeric matrix \code{[n_ages x n_sexes]}. Instantaneous
#'   natural mortality rate.
#' @param fish_sel Numeric array \code{[n_ages x n_sexes x n_fleets]}.
#'   Fishery selectivity for each fleet.
#' @param ret_sel Numeric array \code{[n_ages x n_sexes x n_fleets]}.
#'   Retention selectivity for each fleet, i.e. the proportion of selected
#'   catch that is retained. Default is an array of 1s matching
#'   \code{dim(fish_sel)} (full retention for all fleets).
#' @param dmr Numeric vector \code{[n_fleets]}. Fleet-specific discard
#'   mortality rates, i.e. the fraction of discarded fish that die.
#'   Default \code{rep(0, length(target_catch))} (all discards survive).
#' @param f_init Numeric scalar or vector \code{[n_fleets]}. Starting values
#'   for the \eqn{F} solver. If a scalar is supplied it is recycled across
#'   all fleets. Default \code{0.05}.
#' @param control Named list of control parameters passed to
#'   \code{\link[nleqslv]{nleqslv}}. Default \code{list(btol = 1e-6)}.
#'
#' @return Numeric vector \code{[n_fleets]} of solved fishing mortality rates,
#'   one per fleet.
#'
#' @export catch_to_F_multifleet
#' @family Closed Loop Simulations
catch_to_F_multifleet <- function(target_catch,
                                  NAA,
                                  WAA,
                                  natmort,
                                  fish_sel,
                                  ret_sel = array(1, dim = dim(fish_sel)),
                                  dmr = rep(0, length(target_catch)),
                                  f_init = 0.05,
                                  control = list(btol = 1e-6)) {

  n_fleets <- length(target_catch)

  # Expand f_init if scalar
  if(length(f_init) == 1) f_init <- rep(f_init, n_fleets)

  # Function to minimize: difference between predicted and target catch for all fleets
  catch_diff <- function(f_vec) {
    pred_catches <- numeric(n_fleets)

    for(f in 1:n_fleets) {

      # retained F-at-age for this fleet
      ret_FAA <- f_vec[f] * fish_sel[, , f] * ret_sel[,, f]
      disc_FAA <- f_vec[f] * fish_sel[, , f] * (1 - ret_sel[,, f]) * dmr[f]

      # Total Z includes F from ALL fleets
      ZAA_total <- natmort
      for(ff in 1:n_fleets) {
        ZAA_total <- ZAA_total + (f_vec[ff] * fish_sel[, , ff] * ret_sel[,, ff]) +  # retained F
          (f_vec[ff] * fish_sel[, , ff] * (1 - ret_sel[,, ff]) * dmr[ff]) # discarded F
      }

      # Predicted catch for this fleet (Baranov catch equation)
      pred_catches[f] <- sum((ret_FAA / ZAA_total * NAA * (1 - exp(-ZAA_total))) * WAA)
    }

    return(pred_catches - target_catch)  # Difference from target
  }

  # Solve for F vector
  result <- nleqslv::nleqslv(f_init, catch_diff, control = control)

  return(result$x)
}
