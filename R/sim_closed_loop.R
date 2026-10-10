# Operating model
#
# Closed loop simulation: condition an operating model on a fitted assessment, then run the assessment
# and a control rule forward against it. catch_to_F_om inverts catch advice to the operating model's F.

#' Construct and Condition Closed-Loop Simulation Inputs
#'
#' Builds the simulation object a closed-loop projection runs on. Historical years
#' are conditioned on the fitted model's report, and projection years are extended
#' by the \code{*_fill} rules or by inputs supplied through \code{...}, which
#' replace any internally generated component.
#'
#' @param closed_loop_yrs Integer projection years beyond the fitted data period.
#' @param n_sims Integer number of stochastic replicates.
#' @param data List. Data object used to fit the assessment model.
#' @param parameters List. Parameter vector from the fitted model.
#' @param mapping List. Parameter mapping used during estimation.
#' @param sd_rep List. Standard deviation report from the fitted model.
#' @param rep List. Report object from the fitted model.
#' @param random Character vector of estimated random effects.
#' @param FishIdx_SE_fill,SrvIdx_SE_fill,FishIdx_SE_pop_fill,SrvIdx_SE_pop_fill How
#'   the pooled and population-specific index standard errors are extended into the
#'   projection years. The population ones default to \code{"mean"}.
#' @param ISS_FishAgeComps_fill,ISS_FishLenComps_fill,ISS_SrvAgeComps_fill,ISS_SrvLenComps_fill,ISS_FishAgeComps_discard_fill,ISS_FishLenComps_discard_fill
#'   How the pooled input sample sizes are extended into the projection years.
#' @param ISS_FishAgeComps_pop_fill,ISS_FishLenComps_pop_fill,ISS_SrvAgeComps_pop_fill,ISS_SrvLenComps_pop_fill,ISS_FishAgeComps_discard_pop_fill,ISS_FishLenComps_discard_pop_fill
#'   The population-specific counterparts, each defaulting to \code{"mean"}.
#'
#'   Every \code{*_fill} argument takes \code{"zeros"}, \code{"last"} (repeat the
#'   final observed year), \code{"mean"}, \code{"F_pattern"} (scale the sample
#'   sizes with the simulated fishing mortality pattern; fishery input sample
#'   sizes only, pooled and population-specific, retained and discarded), or a
#'   constant. An array passed instead is taken as the fully specified input and
#'   the fill rule is ignored.
#' @param ... Optional named simulation inputs overriding what is generated
#'   internally. Any argument of \code{\link{Setup_Sim_Fishing}},
#'   \code{\link{Setup_Sim_Survey}}, \code{\link{Setup_Sim_Biologicals}},
#'   \code{\link{Setup_Sim_Rec}} or \code{\link{Setup_Sim_Tagging}} may be given;
#'   the common ones are \code{Fmort_input}, \code{fish_sel_input},
#'   \code{fish_q_input}, \code{srv_sel_input}, \code{srv_q_input},
#'   \code{WAA_input}, \code{MatAA_input}, \code{natmort_input}, \code{R0_input},
#'   \code{rinit_input}, \code{h_input}, \code{Rec_input},
#'   \code{conv_tag_fish_reporting_input} and \code{Movement}. Dimensions must
#'   match the model structure, \code{length(data$years) + closed_loop_yrs} years
#'   and \code{n_sims} simulations. \code{bias_correct_pe} and
#'   \code{bias_correct_oe} of \code{\link{Setup_Sim_Dim}} may be given as well;
#'   without them the operating model takes the fit's, as it does
#'   \code{sigmaR_switch}.
#'
#' @details
#' The conditioning years are the fitted model's own, reconstructed from its report
#' objects; the projection years are simulated under closed-loop management and
#' extended by the fill rules or the supplied inputs. By default the biological
#' inputs, selectivity and catchability are extended at the final estimated year's
#' values, fishing mortality starts at zero in the projection years, recruitment is
#' simulated forward unless \code{Rec_input} specifies it, and the
#' population-specific data sources fall back to uninformative defaults when their
#' \code{Use*_pop} flags hold no ones. Feedback begins in the first projection
#' year.
#'
#' @return A \code{sim_list} holding the model dimensions and simulation
#'   containers, the biological inputs, the pooled and population-specific fishing
#'   and survey processes, recruitment, tagging and movement, each replicated
#'   across \code{n_sims}.
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

  # old reports have no season dim
  rep$natmort <- expand_natmort_seasons(rep$natmort, data$n_seas)

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
                            } else {
                              dim(data$AgeingError)[3] # otherwise the ageing error's observed ages, which the at-age data sit on
                            },
                            n_lens = length(data$lens),
                            n_obs_lens = if(is.null(data$LenBinMap)) length(data$lens) else ncol(data$LenBinMap), # length bins the comps are recorded on
                            n_caal_lens = length(caal_row_lens(data)), # age-at-length rows, each covering the bins CAAL_LenBinMap gives it
                            n_sexes = data$n_sexes,
                            n_fish_fleets = data$n_fish_fleets,
                            n_srv_fleets = data$n_srv_fleets,
                            feedback_start_yr = length(data$years),
                            n_seas = data$n_seas,
                            seasdur = data$seasdur,
                            n_pop = data$n_pop,
                            natal_region = data$natal_region,
                            run_feedback = TRUE,
                            bias_correct_pe = if("bias_correct_pe" %in% names(args)) args$bias_correct_pe else if(is.null(data$bias_correct_pe)) "rec" else data$bias_correct_pe,
                            bias_correct_oe = if("bias_correct_oe" %in% names(args)) args$bias_correct_oe else if(is.null(data$bias_correct_oe)) 0 else data$bias_correct_oe
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
    extend_years(replicate(n = sim_list$n_sims, rep$fish_sel[,,seq_along(data$years),,,,,drop = FALSE]), n_years = closed_loop_yrs, 3, fill = 'last')
  } else args$fish_sel_input
  # Retained selectivity
  ret_sel_input <- if(!"ret_sel_input" %in% names(args)) {
    extend_years(replicate(n = sim_list$n_sims, rep$ret_sel[,,seq_along(data$years),,,,,drop = FALSE]), n_years = closed_loop_yrs, 3, fill = 'last')
  } else args$ret_sel_input
  # Catchability: the reported value is the block mean times the fit's deviation, and the operating
  # model takes the mean as its level, with the deviations read back over the conditioning years
  fit_yrs <- seq_along(data$years)
  fish_q_fit <- split_reported_q(rep$fish_q[,fit_yrs,,drop = FALSE], optim_parameters_list$ln_fish_q,
                                 data$fish_q_blocks[,fit_yrs,,drop = FALSE], data$fish_q_type)
  srv_q_fit <- split_reported_q(rep$srv_q[,fit_yrs,,drop = FALSE], optim_parameters_list$ln_srv_q,
                                data$srv_q_blocks[,fit_yrs,,drop = FALSE], data$srv_q_type)

  # Fishery catchability
  fish_q_input <- if(!"fish_q_input" %in% names(args)) {
    extend_years(replicate(n = sim_list$n_sims, fish_q_fit$q_mean), n_years = closed_loop_yrs, 2, fill = 'last')
  } else args$fish_q_input

  # determine how to draw index se... index with no form = rescale / deweight it, otherwise use the raw se
  idx_draw_se <- function(se, wt, form) if(is.null(form) || form == 0) se / sqrt(wt) else se

  # Fishery index uncertainty
  ObsFishIdx_SE <- if(!"ObsFishIdx_SE" %in% names(args)) {
    extend_years(
      arr = idx_draw_se(data$ObsFishIdx_SE, data$Wt_FishIdx, data$sigmaFishIdx_form),
      n_years = closed_loop_yrs,
      2,
      fill = FishIdx_SE_fill
    )
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
  } else args$ISS_FishAgeComps_discard
  ln_FishAge_discard_theta <- if(!"ln_FishAge_discard_theta" %in% names(args)) optim_parameters_list$ln_FishAge_discard_theta[,,,drop = FALSE] else args$ln_FishAge_discard_theta
  ln_FishAge_discard_theta_agg <- if(!"ln_FishAge_discard_theta_agg" %in% names(args)) optim_parameters_list$ln_FishAge_discard_theta_agg else args$ln_FishAge_discard_theta_agg
  FishAge_discard_corr_pars_agg <- if(!"FishAge_discard_corr_pars_agg" %in% names(args)) optim_parameters_list$FishAge_discard_corr_pars_agg else args$FishAge_discard_corr_pars_agg
  FishAge_discard_corr_pars <- if(!"FishAge_discard_corr_pars" %in% names(args)) optim_parameters_list$FishAge_discard_corr_pars[,,,,drop = FALSE] else args$FishAge_discard_corr_pars

  # Discard Fishery length compositions
  comp_fishlen_discard_like <- if(!"comp_fishlen_discard_like" %in% names(args)) data$FishLenComps_discard_LikeType else args$comp_fishlen_discard_like
  FishLenComps_discard_Type <- if(!"FishLenComps_discard_Type" %in% names(args)) extend_years(data$FishLenComps_discard_Type, closed_loop_yrs, 1, 'last') else args$FishLenComps_discard_Type
  ISS_FishLenComps_discard <- if(!"ISS_FishLenComps_discard" %in% names(args)) {
    extend_years(replicate(sim_list$n_sims, data$ISS_FishLenComps_discard[,,,,,drop = FALSE] * data$Wt_FishLenComps_discard), closed_loop_yrs, 2, fill = ISS_FishLenComps_discard_fill)
  } else args$ISS_FishLenComps_discard
  ln_FishLen_discard_theta <- if(!"ln_FishLen_discard_theta" %in% names(args)) optim_parameters_list$ln_FishLen_discard_theta[,,,drop = FALSE] else args$ln_FishLen_discard_theta
  ln_FishLen_discard_theta_agg <- if(!"ln_FishLen_discard_theta_agg" %in% names(args)) optim_parameters_list$ln_FishLen_discard_theta_agg else args$ln_FishLen_discard_theta_agg
  FishLen_discard_corr_pars_agg <- if(!"FishLen_discard_corr_pars_agg" %in% names(args)) optim_parameters_list$FishLen_discard_corr_pars_agg else args$FishLen_discard_corr_pars_agg
  FishLen_discard_corr_pars <- if(!"FishLen_discard_corr_pars" %in% names(args)) optim_parameters_list$FishLen_discard_corr_pars[,,,,drop = FALSE] else args$FishLen_discard_corr_pars

  # Population-specific fishery index SE
  ObsFishIdx_pop_SE <- if(!"ObsFishIdx_pop_SE" %in% names(args)) {
    if(any(data$UseFishIdx_pop == 1)) {
      extend_years(idx_draw_se(data$ObsFishIdx_pop_SE, data$Wt_FishIdx_pop, data$sigmaFishIdx_pop_form), closed_loop_yrs, 3, fill = FishIdx_SE_pop_fill)
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
  comp_fishage_discard_pop_like <- if(!"comp_fishage_discard_pop_like" %in% names(args)) data$FishAgeComps_discard_pop_LikeType else args$comp_fishage_discard_pop_like
  FishAgeComps_discard_pop_Type <- if(!"FishAgeComps_discard_pop_Type" %in% names(args)) extend_years(data$FishAgeComps_discard_pop_Type, closed_loop_yrs, 1, 'last') else args$FishAgeComps_discard_pop_Type
  ISS_FishAgeComps_discard_pop <- if(!"ISS_FishAgeComps_discard_pop" %in% names(args)) extend_years(replicate(sim_list$n_sims, data$ISS_FishAgeComps_discard_pop[,,,,,,drop = FALSE] * data$Wt_FishAgeComps_discard_pop), closed_loop_yrs, 3, fill = ISS_FishAgeComps_discard_pop_fill) else args$ISS_FishAgeComps_discard_pop
  ln_FishAge_discard_pop_theta <- if(!"ln_FishAge_discard_pop_theta" %in% names(args)) optim_parameters_list$ln_FishAge_discard_pop_theta[,,,,drop = FALSE] else args$ln_FishAge_discard_pop_theta
  ln_FishAge_discard_pop_theta_agg <- if(!"ln_FishAge_discard_pop_theta_agg" %in% names(args)) optim_parameters_list$ln_FishAge_discard_pop_theta_agg else args$ln_FishAge_discard_pop_theta_agg
  FishAge_discard_pop_corr_pars_agg <- if(!"FishAge_discard_pop_corr_pars_agg" %in% names(args)) optim_parameters_list$FishAge_discard_pop_corr_pars_agg else args$FishAge_discard_pop_corr_pars_agg
  FishAge_discard_pop_corr_pars <- if(!"FishAge_discard_pop_corr_pars" %in% names(args)) optim_parameters_list$FishAge_discard_pop_corr_pars[,,,,,drop = FALSE] else args$FishAge_discard_pop_corr_pars

  # Discarded Population-specific fishery length compositions
  comp_fishlen_discard_pop_like <- if(!"comp_fishlen_discard_pop_like" %in% names(args)) data$FishLenComps_discard_pop_LikeType else args$comp_fishlen_discard_pop_like
  FishLenComps_discard_pop_Type <- if(!"FishLenComps_discard_pop_Type" %in% names(args)) extend_years(data$FishLenComps_discard_pop_Type, closed_loop_yrs, 1, 'last') else args$FishLenComps_discard_pop_Type
  ISS_FishLenComps_discard_pop <- if(!"ISS_FishLenComps_discard_pop" %in% names(args)) extend_years(replicate(sim_list$n_sims, data$ISS_FishLenComps_discard_pop[,,,,,,drop = FALSE] * data$Wt_FishLenComps_discard_pop), closed_loop_yrs, 3, fill = ISS_FishLenComps_discard_pop_fill) else args$ISS_FishLenComps_discard_pop
  ln_FishLen_discard_pop_theta <- if(!"ln_FishLen_discard_pop_theta" %in% names(args)) optim_parameters_list$ln_FishLen_discard_pop_theta[,,,,drop = FALSE] else args$ln_FishLen_discard_pop_theta
  ln_FishLen_discard_pop_theta_agg <- if(!"ln_FishLen_discard_pop_theta_agg" %in% names(args)) optim_parameters_list$ln_FishLen_discard_pop_theta_agg else args$ln_FishLen_discard_pop_theta_agg
  FishLen_discard_pop_corr_pars_agg <- if(!"FishLen_discard_pop_corr_pars_agg" %in% names(args)) optim_parameters_list$FishLen_discard_pop_corr_pars_agg else args$FishLen_discard_pop_corr_pars_agg
  FishLen_discard_pop_corr_pars <- if(!"FishLen_discard_pop_corr_pars" %in% names(args)) optim_parameters_list$FishLen_discard_pop_corr_pars[,,,,,drop = FALSE] else args$FishLen_discard_pop_corr_pars

  n_obs_om <- sim_list$n_obs_ages # observed ages the operating model draws on
  catch_aa_used <- any(data$UseCatchAA == 1) # whether the fit observes catch at age
  discard_aa_used <- any(data$UseDiscardAA == 1) # discards at age
  srv_idx_aa_used <- any(data$UseSrvIdxAA == 1) # survey index at age
  catch_aa_pop_used <- any(data$UseCatchAA_pop == 1) # population-specific catch at age
  discard_aa_pop_used <- any(data$UseDiscardAA_pop == 1) # population-specific discards at age
  srv_idx_aa_pop_used <- any(data$UseSrvIdxAA_pop == 1) # population-specific survey index at age

  # setup fishery simulation processes
  sim_list <- Setup_Sim_Fishing(
    sim_list = sim_list, # update simulate list
    ln_sigmaC = ln_sigmaC,
    # age-disaggregated observation error passes through unweighted: it is keyed by age and fleet
    # and is not scaled by a likelihood weight the way the aggregated sigmas are
    ln_sigmaCAA = unused_at_age_on_obs_ages(optim_parameters_list$ln_sigmaCAA, catch_aa_used, 1, n_obs_om, log(0.5)),
    ln_sigmaDAA = unused_at_age_on_obs_ages(optim_parameters_list$ln_sigmaDAA, discard_aa_used, 1, n_obs_om, log(0.5)),
    # the at-age data sources have year on their second dim and no simulation dim, and closed-loop
    # years keep observing whatever the terminal year observed
    UseCatchAA = unused_at_age_on_obs_ages(extend_years(data$UseCatchAA, closed_loop_yrs, 2, fill = 'last'), catch_aa_used, 4, n_obs_om),
    UseDiscardAA = unused_at_age_on_obs_ages(extend_years(data$UseDiscardAA, closed_loop_yrs, 2, fill = 'last'), discard_aa_used, 4, n_obs_om),
    use_catch_aa = data$use_catch_aa,
    use_discard_aa = data$use_discard_aa,
    # the fit's correlation of each fleet's at-age residuals, under the estimation model's names
    AgeObsCorr_catch = data$AgeObsCorr_catch,
    AgeObsCorr_discard = data$AgeObsCorr_discard,
    trans_rho_catch = if(catch_aa_used) optim_parameters_list$trans_rho_catch,
    trans_rho_catch_year = if(catch_aa_used) optim_parameters_list$trans_rho_catch_year,
    trans_rho_catch_us = if(catch_aa_used) optim_parameters_list$trans_rho_catch_us,
    trans_rho_discard = if(discard_aa_used) optim_parameters_list$trans_rho_discard,
    trans_rho_discard_year = if(discard_aa_used) optim_parameters_list$trans_rho_discard_year,
    trans_rho_discard_us = if(discard_aa_used) optim_parameters_list$trans_rho_discard_us,
    # the population-specific at-age data sources and the year totals, as the fit reads them
    CatchAA_seas_Type = data$CatchAA_seas_Type,
    DiscardAA_seas_Type = data$DiscardAA_seas_Type,
    UseCatchAA_pop = if(catch_aa_pop_used) extend_years(data$UseCatchAA_pop, closed_loop_yrs, 3, fill = 'last'),
    UseDiscardAA_pop = if(discard_aa_pop_used) extend_years(data$UseDiscardAA_pop, closed_loop_yrs, 3, fill = 'last'),
    ln_sigmaCAA_pop = if(catch_aa_pop_used) optim_parameters_list$ln_sigmaCAA_pop,
    ln_sigmaDAA_pop = if(discard_aa_pop_used) optim_parameters_list$ln_sigmaDAA_pop,
    ObsCatchAA_pop_SE = if(catch_aa_pop_used) extend_years(data$ObsCatchAA_pop_SE, closed_loop_yrs, 3, fill = 'last'),
    ObsDiscardAA_pop_SE = if(discard_aa_pop_used) extend_years(data$ObsDiscardAA_pop_SE, closed_loop_yrs, 3, fill = 'last'),
    CatchAA_pop_Type = data$CatchAA_pop_Type,
    DiscardAA_pop_Type = data$DiscardAA_pop_Type,
    CatchAA_pop_LikeType = data$CatchAA_pop_LikeType,
    DiscardAA_pop_LikeType = data$DiscardAA_pop_LikeType,
    CatchAA_pop_sigma_form = data$CatchAA_pop_sigma_form,
    DiscardAA_pop_sigma_form = data$DiscardAA_pop_sigma_form,
    CatchAA_pop_seas_Type = data$CatchAA_pop_seas_Type,
    DiscardAA_pop_seas_Type = data$DiscardAA_pop_seas_Type,
    AgeObsCorr_catch_pop = data$AgeObsCorr_catch_pop,
    AgeObsCorr_discard_pop = data$AgeObsCorr_discard_pop,
    trans_rho_catch_pop = if(catch_aa_pop_used) optim_parameters_list$trans_rho_catch_pop,
    trans_rho_catch_pop_year = if(catch_aa_pop_used) optim_parameters_list$trans_rho_catch_pop_year,
    trans_rho_catch_pop_us = if(catch_aa_pop_used) optim_parameters_list$trans_rho_catch_pop_us,
    trans_rho_discard_pop = if(discard_aa_pop_used) optim_parameters_list$trans_rho_discard_pop,
    trans_rho_discard_pop_year = if(discard_aa_pop_used) optim_parameters_list$trans_rho_discard_pop_year,
    trans_rho_discard_pop_us = if(discard_aa_pop_used) optim_parameters_list$trans_rho_discard_pop_us,
    ObsCatchAA_SE = unused_at_age_on_obs_ages(extend_years(data$ObsCatchAA_SE, closed_loop_yrs, 2, fill = 'last'), catch_aa_used, 4, n_obs_om),
    ObsDiscardAA_SE = unused_at_age_on_obs_ages(extend_years(data$ObsDiscardAA_SE, closed_loop_yrs, 2, fill = 'last'), discard_aa_used, 4, n_obs_om),
    CatchAA_Type = extend_years(data$CatchAA_Type, closed_loop_yrs, 1, fill = 'last'),
    DiscardAA_Type = extend_years(data$DiscardAA_Type, closed_loop_yrs, 1, fill = 'last'),
    CatchAA_LikeType = data$CatchAA_LikeType, DiscardAA_LikeType = data$DiscardAA_LikeType,
    # data sources the fitted model reports once a year stay annual in the operating model
    Catch_seas_Type = data$Catch_seas_Type,
    Catch_pop_seas_Type = data$Catch_pop_seas_Type,
    FishIdx_seas_Type = data$FishIdx_seas_Type,
    FishIdx_pop_seas_Type = data$FishIdx_pop_seas_Type,
    FishAgeComps_seas_Type = data$FishAgeComps_seas_Type,
    Discard_seas_Type = data$Discard_seas_Type,
    Discard_pop_seas_Type = data$Discard_pop_seas_Type,
    FishLenComps_seas_Type = data$FishLenComps_seas_Type,
    FishAgeComps_pop_seas_Type = data$FishAgeComps_pop_seas_Type,
    FishLenComps_pop_seas_Type = data$FishLenComps_pop_seas_Type,
    FishAgeComps_discard_seas_Type = data$FishAgeComps_discard_seas_Type,
    FishLenComps_discard_seas_Type = data$FishLenComps_discard_seas_Type,
    FishAgeComps_discard_pop_seas_Type = data$FishAgeComps_discard_pop_seas_Type,
    FishLenComps_discard_pop_seas_Type = data$FishLenComps_discard_pop_seas_Type,
    seas_agg_slot = seas_agg_slot_list(data, length(data$years) + closed_loop_yrs, "fish"), # the season each year total sits in, the last fitted year's after the fit
    CatchAA_sigma_form = data$CatchAA_sigma_form, DiscardAA_sigma_form = data$DiscardAA_sigma_form,
    ln_sigmaC_pop = ln_sigmaC_pop,
    Fmort_input = extend_years(replicate(n = sim_list$n_sims, rep$Fmort[,seq_along(data$years),,,drop = FALSE]), n_years = closed_loop_yrs, 2, fill = 'zeros'),
    ln_sigmaD = ln_sigmaD,
    ln_sigmaD_pop = ln_sigmaD_pop,
    dmr_input = extend_years(replicate(n = sim_list$n_sims, rep$dmr[,seq_along(data$years),,,drop = FALSE]), n_years = closed_loop_yrs, 2, fill = 'zeros'),
    fish_sel_input = fish_sel_input,
    ret_sel_input = ret_sel_input,
    # length comps selected at length read the fit's selectivity at length, kept at its last year after the data
    FishLenComps_sel = if(is.null(data$fish_len_comp_sel)) rep("age", data$n_fish_fleets) else ifelse(data$fish_len_comp_sel == 1, "length", "age"),
    fish_sel_l_input = if(any(data$fish_len_comp_sel == 1)) extend_years(replicate(n = sim_list$n_sims, rep$fish_sel_l[,seq_along(data$years),,,,drop = FALSE]), closed_loop_yrs, 2, 'last') else NULL,
    ret_sel_l_input = if(any(data$fish_len_comp_sel == 1) && isTRUE(data$ret_selex_type == 1)) extend_years(replicate(n = sim_list$n_sims, rep$ret_sel_l[,seq_along(data$years),,,,drop = FALSE]), closed_loop_yrs, 2, 'last') else NULL,
    fish_q_input = fish_q_input,
    ObsFishIdx_SE = ObsFishIdx_SE,
    sigmaFishIdx_form = if(is.null(data$sigmaFishIdx_form)) 0 else data$sigmaFishIdx_form, # any estimated part of the index sd, drawn on top of the reported errors
    ln_sigmaFishIdx = optim_parameters_list$ln_sigmaFishIdx,
    sigmaFishIdx_pop_form = if(is.null(data$sigmaFishIdx_pop_form)) 0 else data$sigmaFishIdx_pop_form,
    ln_sigmaFishIdx_pop = optim_parameters_list$ln_sigmaFishIdx_pop,
    ObsFishIdx_pop_SE = ObsFishIdx_pop_SE,
    fish_idx_type = data$fish_idx_type,
    fish_idx_ages = data$fish_idx_ages, # ages in each index total
    t_fish = if(is.null(data$t_fish)) array(0, dim = c(data$n_regions, data$n_seas, data$n_fish_fleets)) else data$t_fish, # fishery index timing
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
    ISS_FishLenComps = ISS_FishLenComps,
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
    ISS_FishLenComps_discard = ISS_FishLenComps_discard,
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
    extend_years(replicate(n = sim_list$n_sims, rep$srv_sel[,,seq_along(data$years),,,,,drop = FALSE]), closed_loop_yrs, 3, 'last')
  } else args$srv_sel_input
  # Survey catchability / q
  srv_q_input <- if(!"srv_q_input" %in% names(args)) {
    extend_years(replicate(n = sim_list$n_sims, srv_q_fit$q_mean), closed_loop_yrs, 2, 'last')
  } else args$srv_q_input
  # Survey index uncertainty
  ObsSrvIdx_SE <- if(!"ObsSrvIdx_SE" %in% names(args)) {
    extend_years(
      arr = idx_draw_se(data$ObsSrvIdx_SE, data$Wt_SrvIdx, data$sigmaSrvIdx_form),
      n_years = closed_loop_yrs,
      2,
      fill = SrvIdx_SE_fill
    )
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
      extend_years(idx_draw_se(data$ObsSrvIdx_pop_SE, data$Wt_SrvIdx_pop, data$sigmaSrvIdx_pop_form), closed_loop_yrs, 3, fill = SrvIdx_SE_pop_fill)
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
    # length comps selected at length read the fit's selectivity at length, kept at its last year after the data
    SrvLenComps_sel = if(is.null(data$srv_len_comp_sel)) rep("age", data$n_srv_fleets) else ifelse(data$srv_len_comp_sel == 1, "length", "age"),
    srv_sel_l_input = if(any(data$srv_len_comp_sel == 1)) extend_years(replicate(n = sim_list$n_sims, rep$srv_sel_l[,seq_along(data$years),,,,drop = FALSE]), closed_loop_yrs, 2, 'last') else NULL,
    srv_q_input = srv_q_input,
    ObsSrvIdx_SE = ObsSrvIdx_SE,
    sigmaSrvIdx_form = if(is.null(data$sigmaSrvIdx_form)) 0 else data$sigmaSrvIdx_form, # any estimated part of the index sd, drawn on top of the reported errors
    ln_sigmaSrvIdx = optim_parameters_list$ln_sigmaSrvIdx,
    sigmaSrvIdx_pop_form = if(is.null(data$sigmaSrvIdx_pop_form)) 0 else data$sigmaSrvIdx_pop_form,
    ln_sigmaSrvIdx_pop = optim_parameters_list$ln_sigmaSrvIdx_pop,
    ln_sigmaSrvIdxAA = unused_at_age_on_obs_ages(optim_parameters_list$ln_sigmaSrvIdxAA, srv_idx_aa_used, 1, n_obs_om, log(0.5)),
    UseSrvIdxAA = unused_at_age_on_obs_ages(extend_years(data$UseSrvIdxAA, closed_loop_yrs, 2, fill = 'last'), srv_idx_aa_used, 4, n_obs_om),
    use_srv_idx_aa = data$use_srv_idx_aa,
    AgeObsCorr_srv_idx = data$AgeObsCorr_srv_idx,
    trans_rho_srv_idx = if(srv_idx_aa_used) optim_parameters_list$trans_rho_srv_idx,
    trans_rho_srv_idx_year = if(srv_idx_aa_used) optim_parameters_list$trans_rho_srv_idx_year,
    trans_rho_srv_idx_us = if(srv_idx_aa_used) optim_parameters_list$trans_rho_srv_idx_us,
    SrvIdxAA_seas_Type = data$SrvIdxAA_seas_Type,
    UseSrvIdxAA_pop = if(srv_idx_aa_pop_used) extend_years(data$UseSrvIdxAA_pop, closed_loop_yrs, 3, fill = 'last'),
    ln_sigmaSrvIdxAA_pop = if(srv_idx_aa_pop_used) optim_parameters_list$ln_sigmaSrvIdxAA_pop,
    ObsSrvIdxAA_pop_SE = if(srv_idx_aa_pop_used) extend_years(data$ObsSrvIdxAA_pop_SE, closed_loop_yrs, 3, fill = 'last'),
    SrvIdxAA_pop_Type = data$SrvIdxAA_pop_Type,
    SrvIdxAA_pop_LikeType = data$SrvIdxAA_pop_LikeType,
    SrvIdxAA_pop_sigma_form = data$SrvIdxAA_pop_sigma_form,
    SrvIdxAA_pop_seas_Type = data$SrvIdxAA_pop_seas_Type,
    AgeObsCorr_srv_idx_pop = data$AgeObsCorr_srv_idx_pop,
    trans_rho_srv_idx_pop = if(srv_idx_aa_pop_used) optim_parameters_list$trans_rho_srv_idx_pop,
    trans_rho_srv_idx_pop_year = if(srv_idx_aa_pop_used) optim_parameters_list$trans_rho_srv_idx_pop_year,
    trans_rho_srv_idx_pop_us = if(srv_idx_aa_pop_used) optim_parameters_list$trans_rho_srv_idx_pop_us,
    ObsSrvIdxAA_SE = unused_at_age_on_obs_ages(extend_years(data$ObsSrvIdxAA_SE, closed_loop_yrs, 2, fill = 'last'), srv_idx_aa_used, 4, n_obs_om),
    SrvIdxAA_Type = extend_years(data$SrvIdxAA_Type, closed_loop_yrs, 1, fill = 'last'),
    SrvIdxAA_LikeType = data$SrvIdxAA_LikeType, SrvIdxAA_sigma_form = data$SrvIdxAA_sigma_form,
    # data sources the fitted model reports once a year stay annual in the operating model
    SrvIdx_seas_Type = data$SrvIdx_seas_Type,
    SrvIdx_pop_seas_Type = data$SrvIdx_pop_seas_Type,
    SrvAgeComps_seas_Type = data$SrvAgeComps_seas_Type,
    SrvLenComps_seas_Type = data$SrvLenComps_seas_Type,
    SrvAgeComps_pop_seas_Type = data$SrvAgeComps_pop_seas_Type,
    SrvLenComps_pop_seas_Type = data$SrvLenComps_pop_seas_Type,
    seas_agg_slot = seas_agg_slot_list(data, length(data$years) + closed_loop_yrs, "srv"), # the season each year total sits in, the last fitted year's after the fit
    ObsSrvIdx_pop_SE = ObsSrvIdx_pop_SE,
    srv_idx_type = data$srv_idx_type,
    srv_idx_ages = data$srv_idx_ages, # ages in each index total
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
    extend_years(replicate(n = sim_list$n_sims, truncate_years(rep$natmort, length(data$years))), closed_loop_yrs, 3, 'last')
  } else args$natmort_input
  # biologicals the growth module derives come from the report, the rest from the data
  WAA_input <- if(!"WAA_input" %in% names(args)) {
    extend_years(replicate(n = sim_list$n_sims, (if(is.null(rep$WAA)) data$WAA else rep$WAA)[,,seq_along(data$years),,,,drop = FALSE]), closed_loop_yrs, 3, 'last')
  } else args$WAA_input
  WAA_fish_input <- if(!"WAA_fish_input" %in% names(args)) {
    extend_years(replicate(n = sim_list$n_sims, (if(is.null(rep$WAA_fish)) data$WAA_fish else rep$WAA_fish)[,,seq_along(data$years),,,,,drop = FALSE]), closed_loop_yrs, 3, 'last')
  } else args$WAA_fish_input
  WAA_srv_input <- if(!"WAA_srv_input" %in% names(args)) {
    extend_years(replicate(n = sim_list$n_sims, (if(is.null(rep$WAA_srv)) data$WAA_srv else rep$WAA_srv)[,,seq_along(data$years),,,,,drop = FALSE]), closed_loop_yrs, 3, 'last')
  } else args$WAA_srv_input
  MatAA_input <- if(!"MatAA_input" %in% names(args)) {
    extend_years(replicate(n = sim_list$n_sims, data$MatAA[,,seq_along(data$years),,,,drop = FALSE]), closed_loop_yrs, 3, 'last')
  } else args$MatAA_input
  AgeingError_input <- if(!"AgeingError_input" %in% names(args)) {
    extend_years(replicate(n = sim_list$n_sims, data$AgeingError[seq_along(data$years),,,drop = FALSE]), closed_loop_yrs, 1, 'last')
  } else args$AgeingError_input
  AgeingError_fish_input <- if(!"AgeingError_fish_input" %in% names(args)) {
    if(is.null(data$AgeingError_fish)) NULL else extend_years(replicate(n = sim_list$n_sims, data$AgeingError_fish[seq_along(data$years),,,,drop = FALSE]), closed_loop_yrs, 1, 'last')
  } else args$AgeingError_fish_input
  AgeingError_srv_input <- if(!"AgeingError_srv_input" %in% names(args)) {
    if(is.null(data$AgeingError_srv)) NULL else extend_years(replicate(n = sim_list$n_sims, data$AgeingError_srv[seq_along(data$years),,,,drop = FALSE]), closed_loop_yrs, 1, 'last')
  } else args$AgeingError_srv_input
  SizeAgeTrans_input <- if(!"SizeAgeTrans_input" %in% names(args)) {
    if(data$fit_lengths == 0) array(NA, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_lens, sim_list$n_ages, sim_list$n_sexes))
    if(data$fit_lengths == 1 && !is.null(data$SizeAgeTrans) && !all(is.na(data$SizeAgeTrans))) extend_years(replicate(n = sim_list$n_sims, data$SizeAgeTrans[,,seq_along(data$years),,,,,drop = FALSE]), closed_loop_yrs, 3, 'last')
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
    LenBinMap_input = data$LenBinMap, # model length bins onto the recorded ones, as the fit maps them
    CAAL_LenBinMap_input = data$CAAL_LenBinMap, # model length bins each age-at-length row covers
    SizeAgeTrans_input = SizeAgeTrans_input,
    # keys per fleet from the growth module, each at its fleet's own timing
    SizeAgeTrans_fish_input = if(is.null(rep$SizeAgeTrans_fish)) NULL else extend_years(replicate(n = sim_list$n_sims, rep$SizeAgeTrans_fish[,,seq_along(data$years),,,,,,drop = FALSE]), closed_loop_yrs, 3, 'last'),
    SizeAgeTrans_srv_input = if(is.null(rep$SizeAgeTrans_srv)) NULL else extend_years(replicate(n = sim_list$n_sims, rep$SizeAgeTrans_srv[,,seq_along(data$years),,,,,,drop = FALSE]), closed_loop_yrs, 3, 'last')
  )

  # growth deviations stuff; conditional during conditioning years, new draws after
  sim_list <- Setup_Sim_Growth_RE(sim_list, data, optim_parameters_list, rep = rep)

  # Setup Recruitment Processes ---------------------------------------------
  h_input <- if(!"h_input" %in% names(args)) {
    replicate(n = sim_list$n_sims, array(rep$h_trans, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs)))
  } else args$h_input
  warn_R0_ref_block_om(data, "condition_closed_loop_simulations")
  R0_input <- if(!"R0_input" %in% names(args)) {
    # R0 can have time blocks, so the conditioned operating model holds the terminal year's value
    # forward through the projection. identical to the year-by-year value in an unblocked model
    n_hist <- length(data$years)
    R0_by_yr <- if(is.null(rep$R0_yr)) matrix(rep$R0, sim_list$n_pop, n_hist) else rep$R0_yr[, 1:n_hist, drop = FALSE]
    R0_r = array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs)) # container
    for(p in 1:sim_list$n_pop) for(r in 1:sim_list$n_regions)
      R0_r[p,r,] = c(R0_by_yr[p,], rep(R0_by_yr[p, n_hist], max(0, sim_list$n_yrs - n_hist)))[1:sim_list$n_yrs] * rep$rec_region_prop[p,r]
    replicate(n = sim_list$n_sims, expr = R0_r)
  } else args$R0_input
  rinit_input <- if(!"rinit_input" %in% names(args)) {
    rinit_r = array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_sims)) # container
    for(p in 1:sim_list$n_pop) for(s in 1:sim_list$n_sims) rinit_r[p,,s] <- rep$rinit[p] * rep$rec_region_prop[p,]
    rinit_r
  } else args$rinit_input
  sexratio_input <- if(!"sexratio_input" %in% names(args)) {
    extend_years(replicate(n = sim_list$n_sims, expr = rep$sexratio[,,seq_along(data$years),,drop = FALSE]), closed_loop_yrs, 3, 'last')
  } else args$sexratio_input
  ln_sigmaR <- if(!"ln_sigmaR" %in% names(args)) optim_parameters_list$ln_sigmaR else args$ln_sigmaR
  stray_rate_input <- if(!"stray_rate_input" %in% names(args)) {
    extend_years(replicate(n = sim_list$n_sims, expr = data$stray_rate[,seq_along(data$years),drop = FALSE]), closed_loop_yrs, 2, 'last')
  } else args$stray_rate_input
  Rec_input <- if(!"Rec_input" %in% names(args)) {
    replicate(n = sim_list$n_sims, expr = rep$Rec[,,seq_along(data$years),drop = FALSE])
  } else args$Rec_input
  ln_InitDevs_input <- if(!"ln_InitDevs_input" %in% names(args)) {
    replicate(sim_list$n_sims, optim_parameters_list$ln_InitDevs)
  } else args$ln_InitDevs_input

  # recruitment options
  rec_dd <- if(!"rec_dd" %in% names(args)) data$rec_dd else args$rec_dd
  init_dd <- if(!"init_dd" %in% names(args)) data$rec_dd else args$init_dd
  rec_lag <- if(!"rec_lag" %in% names(args)) data$rec_lag else args$rec_lag
  recruitment_opt <- if(!"recruitment_opt" %in% names(args)) data$rec_model else args$recruitment_opt
  RecDevs_model <- if(!"RecDevs_model" %in% names(args)) c("iid", "rw", "ar1")[if(is.null(data$RecDevs_model)) 1 else data$RecDevs_model] else args$RecDevs_model
  RecDevs_rho <- if(!"RecDevs_rho" %in% names(args)) {
    if(is.null(optim_parameters_list$RecDevs_rho)) array(0, dim = c(sim_list$n_pop, sim_list$n_regions))
    else rho_trans(array(optim_parameters_list$RecDevs_rho, dim = c(sim_list$n_pop, sim_list$n_regions)))
  } else args$RecDevs_rho
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
    sigmaR_switch = if("sigmaR_switch" %in% names(args)) args$sigmaR_switch else if(is.null(data$sigmaR_switch)) 1 else data$sigmaR_switch, # first year on the late sigmaR
    Rec_input = Rec_input,
    ln_InitDevs_input = ln_InitDevs_input,
    recruitment_opt = recruitment_opt,
    rec_dd = rec_dd,
    init_dd = init_dd,
    rec_lag = rec_lag,
    SR_ref_yr = if(is.null(data$SR_ref_yr)) 1 else data$SR_ref_yr, # match the fit's per-recruit reference year
    RecDevs_model = RecDevs_model, # recruitment process error the projection draws under
    RecDevs_rho = RecDevs_rho, # ar1 correlation, natural scale
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
  conv_fish_tag_like <- if(!"conv_fish_tag_like" %in% names(args)) data$conv_fish_tag_like else args$conv_fish_tag_like
  ln_conv_fish_tag_theta <- if(!"ln_conv_fish_tag_theta" %in% names(args)) optim_parameters_list$ln_conv_fish_tag_theta else args$ln_conv_fish_tag_theta

  # setup tagging simulation
  if(!is.null(n_tags)) sim_list$n_tags_rel_input <- NULL # set release input to NULL if n_tags is specified.
  sim_list <- Setup_Sim_Tagging(
    sim_list = sim_list,
    conv_tag_max_liberty = data$conv_tag_max_liberty,
    conv_tag_t_tagging = data$conv_tag_t_tagging,
    n_tags = n_tags,
    n_tags_rel_input = n_tags_rel_input * data$Wt_Tagging,
    conv_tag_release_indicator = conv_tag_release_indicator,
    conv_tag_release_platform = if(is.null(data$conv_tag_release_platform)) default_tag_release_platform(conv_tag_release_indicator) else data$conv_tag_release_platform,
    ln_init_conv_tag_mort = ln_init_conv_tag_mort,
    ln_conv_tag_shed = ln_conv_tag_shed,
    conv_fish_tag_attr = data$conv_fish_tag_attr,
    conv_tag_fish_reporting_input = conv_tag_fish_reporting_input,
    use_conv_fish_tagging = use_conv_fish_tagging,
    conv_fish_tag_like = conv_fish_tag_like,
    ln_conv_fish_tag_theta = ln_conv_fish_tag_theta,
    conv_tag_pop_pool = data$conv_tag_pop_pool, # recaptures the fit pools into one count
    conv_tag_age_pool = data$conv_tag_age_pool,
    conv_tag_sex_pool = data$conv_tag_sex_pool,
    conv_tagged_fish_input = if("conv_tagged_fish_input" %in% names(args)) args$conv_tagged_fish_input
                             else if(!anyNA(data$conv_tagged_fish)) data$conv_tagged_fish * data$Wt_Tagging # the fit's releases by age
  )

  # Movement ----------------------------------------------------------------
  Movement <- if(!"Movement" %in% names(args)) extend_years(replicate(n = sim_list$n_sims, rep$Movement[,,,seq_along(data$years),,,,drop = FALSE]), closed_loop_yrs, 4, 'last') else args$Movement
  sim_list$Movement <- Movement
  sgl_seas_spawning_movement <- if(!"sgl_seas_spawning_movement" %in% names(args)) extend_years(replicate(n = sim_list$n_sims, data$sgl_seas_spawning_movement[,,,seq_along(data$years),,,drop = FALSE]), closed_loop_yrs, 4, 'last') else args$sgl_seas_spawning_movement
  sim_list$sgl_seas_spawning_movement <- sgl_seas_spawning_movement
  sim_list$move_timing <- if(!"move_timing" %in% names(args)) {
    if(is.null(data$move_timing)) 0 else data$move_timing
  } else args$move_timing
  sim_list$expm_nsub <- if(!"move_expm_nsub" %in% names(args)) {
    if(is.null(data$move_expm_nsub)) 0 else data$move_expm_nsub
  } else args$move_expm_nsub
  sim_list$Mrate <- if("Mrate" %in% names(args)) args$Mrate else {
    if(sim_list$move_timing != 2 || is.null(rep$Mrate)) NULL
    else extend_years(replicate(n = sim_list$n_sims, rep$Mrate[,,,seq_along(data$years),,,,drop = FALSE]), closed_loop_yrs, 4, 'last')
  }

  # movement deviations; conditional during conditioning years, new draws after
  sim_list <- Setup_Sim_Movement(sim_list, data, optim_parameters_list)

  # selectivity deviations; the fit's over the conditioning years, and then drawn on through the projection. F is
  # the control rule's past the fit, so its deviations and discard mortality's stay at the fit
  if(!any(c("fish_sel_input", "ret_sel_input", "srv_sel_input") %in% names(args))) {
    sim_list <- Setup_Sim_Fleet_Devs(sim_list, data, optim_parameters_list, random = random)
    if(!is.null(sim_list$fleet_devs_on)) sim_list$fleet_devs_on[c("F", "dmr")] <- FALSE
  }

  # State-space numbers at age ----------------------------------------------
  sim_list$NAA_re <- if("NAA_re" %in% names(args)) args$NAA_re else {
    if(is.null(data$NAA_re)) 0 else data$NAA_re
  }
  state_on <- isTRUE(sim_list$NAA_re > 0)

  # the fit's process error in naa, the sd kept at the terminal year through the projection years
  naa_process <- if(state_on && !is.null(optim_parameters_list$ln_sigmaNAA)) naa_process_from_fit(data, optim_parameters_list)
  if(!is.null(naa_process)) naa_process$sigmaNAA <- extend_years(naa_process$sigmaNAA, closed_loop_yrs, 3, 'last')
  for(opt_name in c("sigmaNAA", "naa_rho", "NAA_re_pop", "NAA_re_region", "NAA_re_sex", "NAA_re_season",
                    "naa_pop_corr", "naa_region_corr", "naa_sex_corr")) {
    sim_list[[opt_name]] <- if(opt_name %in% names(args)) args[[opt_name]]
                            else if(!is.null(naa_process)) naa_process[[opt_name]]
                            else if(opt_name == "naa_rho") c(age = 0, year = 0, cohort = 0)
                            else 0
  } # end opt_name loop

  # The active ages are reused unchanged, and the active years run on through the projection: the
  # operating model should keep generating process error in the future, not stop at the data.
  sim_list$naa_re_ages <- if("naa_re_ages" %in% names(args)) args$naa_re_ages else {
    if(!state_on || is.null(data$naa_re_ages)) integer(0) else data$naa_re_ages
  }
  sim_list$naa_re_yrs <- if("naa_re_yrs" %in% names(args)) args$naa_re_yrs else {
    if(!state_on || is.null(data$naa_re_yrs) || !length(data$naa_re_yrs)) integer(0)
    else seq(min(data$naa_re_yrs), sim_list$n_yrs)
  }

  # the active seasons are reused unchanged; season one alone is the annual state
  sim_list$naa_re_seas <- if("naa_re_seas" %in% names(args)) args$naa_re_seas else {
    if(!state_on) 1 else if(is.null(data$naa_re_seas)) 1 else data$naa_re_seas
  }
  sim_list$naa_season_corr <- if("naa_season_corr" %in% names(args)) args$naa_season_corr
                              else if(!is.null(naa_process)) naa_process$naa_season_corr
                              else 0

  # number of conditioning years
  sim_list$n_cond_yrs <- if("n_cond_yrs" %in% names(args)) args$n_cond_yrs else length(data$years)

  # the fit's own deviations over the conditioning years and zero after, so a projection walk starts
  # from the fit's last one. a supplied catchability is taken as it stands and gets none
  if(!"fish_q_input" %in% names(args))
    sim_list$ln_fish_q_devs <- extend_years(replicate(n = sim_list$n_sims, fish_q_fit$devs), closed_loop_yrs, 2)
  if(!"srv_q_input" %in% names(args))
    sim_list$ln_srv_q_devs <- extend_years(replicate(n = sim_list$n_sims, srv_q_fit$devs), closed_loop_yrs, 2)

  # conditioning on the numbers-at-age process error
  state_estimates_known <- state_on && !is.null(sd_rep$par.random) && "ln_NAA" %in% names(sd_rep$par.random)
  if(state_on && !state_estimates_known && !("naa_eta_input" %in% names(args)))
    warning("The fit has a state on the numbers at age, but sd_rep carries no random effects for it, so its ",
            "innovations cannot be read and the conditioning years will draw their own. Pass sd_rep with ",
            "par.random, or naa_eta_input directly, for the conditioning years to reproduce the fit.", call. = FALSE)

  # derive NAA PE to use in conditioning period
  if(state_on && state_estimates_known && !("naa_eta_input" %in% names(args)) && !is.null(rep$NAA_pred)) {
    ln_NAA <- optim_parameters_list$ln_NAA
    n_cond <- dim(ln_NAA)[3]
    pred <- rep$NAA_pred[,,seq_len(n_cond),,,,drop = FALSE]
    fit_eta <- array(0, dim = dim(ln_NAA))
    fit_eta[pred > 0] <- ln_NAA[pred > 0] - log(pred[pred > 0])
    eta <- array(0, dim = c(dim(ln_NAA), sim_list$n_sims))
    for(i in seq_len(sim_list$n_sims)) eta[,,,,,,i] <- fit_eta
    sim_list$naa_eta_input <- eta
  } else if("naa_eta_input" %in% names(args)) sim_list$naa_eta_input <- args$naa_eta_input

  # setup dsem - conditional on historical fits
  if(!is.null(data$dsem_model)) sim_list <- Setup_Sim_DSEM(sim_list, data, optim_parameters_list, rep = rep, condition_on_fit = TRUE)

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
      natmort = array(sim_env$natmort[,, 1:y, , , , sim], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_ages, sim_env$n_sexes)),
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
        array(R0_slice / row_sums, dim = c(sim_env$n_pop, sim_env$n_regions))
      },
      Rec = array(sim_env$Rec[, , 1:y, sim], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y))),
      Movement = array(sim_env$Movement[, , , 1:y, , , , sim],  dim = c(sim_env$n_pop, sim_env$n_regions, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_ages, sim_env$n_sexes)),
      sgl_seas_spawning_movement = array(sim_env$sgl_seas_spawning_movement[, , , 1:y, , , sim],  dim = c(sim_env$n_pop, sim_env$n_regions, sim_env$n_regions, length(1:y), sim_env$n_ages, sim_env$n_sexes)),
      move_timing = if(is.null(sim_env$move_timing)) 0 else sim_env$move_timing
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

  return(list(
    f_ref_pt = f_ref_pt,
    b_ref_pt = b_ref_pt,
    virgin_b_ref_pt = virgin_b_ref_pt,
    pop_b_ref_pt = pop_b_ref_pt,
    virgin_pop_b_ref_pt = virgin_pop_b_ref_pt
  ))
}



# Catch Advice To Fishing Mortality --------------------------------------------
#
# Call order, outermost first, once per closed loop year:
#
#   catch_to_F_om()            sweeps the seasons in order, earlier seasons settled first
#    +- solve_om_season_F()    tapes one season's catch with RTMB for solve_log_catch_newton()
#        +- om_season_catch()      what catch does a trial F give this season?
#
# Regions are solved one at a time, except where one region's F changes another's catch: under
# move_timing 2, and in the spawning season under rec_lag 0, where recruits come from every region.

#' Operating Model Catch Over One Season At A Trial F
#'
#' Runs one season of the operating model from the numbers at its start, with
#' the same mortality, movement and catch equation as \code{apply_pop_dy} and
#' \code{generate_fishery_catch_comp_idx}. Only the rec_lag 0 recruitment step
#' writes into \code{sim_env}, values the annual cycle sets again before reading.
#'
#' @param N Numeric array \code{[n_pop, n_regions, n_ages, n_sexes]}. Numbers at
#'   age at the start of the season, with recruits already known added.
#' @param F_seas Numeric matrix \code{[n_regions, n_fish_fleets]}. Trial F.
#' @param y,seas,sim Year, season and replicate.
#' @param sim_env Simulation environment from \code{\link{Setup_sim_env}}.
#' @param target_units Numeric \code{[n_fish_fleets]}. 1 sums catch in biomass
#'   through \code{WAA_fish}, 0 in numbers.
#' @param spawn_recruits Logical. Whether this year's recruits arrive in this
#'   season from its own spawning biomass (\code{rec_lag = 0} at \code{spawn_seas}).
#'
#' @return Named list with \code{catch} \code{[n_regions, n_fish_fleets]} in
#'   \code{target_units}, \code{N_end}, the survivors at the end of the season
#'   before ageing, shaped like \code{N}, and \code{Rec_y} \code{[n_pop, n_regions]},
#'   this year's recruitment when \code{spawn_recruits}, \code{NULL} otherwise.
#' @keywords internal
om_season_catch <- function(N, F_seas, y, seas, sim, sim_env, target_units, spawn_recruits = FALSE) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  # get dimensions
  n_pop <- sim_env$n_pop
  n_regions <- sim_env$n_regions
  n_ages <- sim_env$n_ages
  n_sexes <- sim_env$n_sexes
  n_fish_fleets <- sim_env$n_fish_fleets
  move_timing <- if(is.null(sim_env$move_timing)) 0 else sim_env$move_timing
  expm_nsub <- if(is.null(sim_env$expm_nsub)) 0 else sim_env$expm_nsub
  seasdur <- sim_env$seasdur[seas]

  # total mortality, with every fleet fishing the region adding to it
  ZAA <- array(0, dim = c(n_pop, n_regions, n_ages, n_sexes))
  for(p in 1:n_pop) {
    for(r in 1:n_regions) {

      ZAA[p,r,,] <- sim_env$natmort[p,r,y,seas,,,sim] * seasdur # natural mortality is a rate per year
      for(f in 1:n_fish_fleets) {
        sel <- sim_env$fish_sel[p,r,y,seas,,,f,sim] # total selectivity
        ret <- sim_env$ret_sel[p,r,y,seas,,,f,sim] # proportion of selected fish retained
        dmr <- sim_env$dmr[r,y,seas,f,sim] # proportion of discards that die
        ZAA[p,r,,] <- ZAA[p,r,,] + F_seas[r,f] * sel * (ret + (1 - ret) * dmr) # retained plus dead discards
      } # end f loop

    } # end r loop
  } # end p loop

  # under move_timing 0 fish move at the start of the season and are caught where they land
  N_fished <- N
  if(n_regions > 1 && move_timing == 0) {
    first_age <- if(sim_env$do_recruits_move == 1) 1 else 2 # recruits stay put unless they move from birth
    for(p in 1:n_pop) {
      for(a in first_age:n_ages) {
        for(s in 1:n_sexes) {
          N_fished[p,,a,s] <- as.vector(t(N[p,,a,s]) %*% sim_env$Movement[p,,,y,seas,a,s,sim])
        } # end s loop
      } # end a loop
    } # end p loop
  } # end if movement at the start of the season

  # under rec_lag 0 this year's recruits arrive now, from the spawning biomass this season's F leaves
  Rec_y <- NULL
  if(spawn_recruits) {

    spawn_biom <- compute_biom_y_sim(y, seas, sim, sim_env,
                                     NAA_s = array(N_fished, dim = c(n_pop, n_regions, 1, 1, n_ages, n_sexes)), # survivors only, as apply_pop_dy
                                     ZAA_s = array(ZAA, dim = c(n_pop, n_regions, 1, 1, n_ages, n_sexes)))
    SSB_vals <- array(sim_env$SSB[,,,sim], dim = c(n_pop, n_regions, sim_env$n_yrs))
    SSB_vals[,,y] <- spawn_biom$SSB_y # this year's spawning biomass at the trial F
    det_rec <- array(sim_det_recruitment(y, sim, sim_env, SSB_vals = SSB_vals), dim = c(n_pop, n_regions))

    # deviations drawn at the end of last year, or the conditioned recruitment where Rec_input covers the year
    use_rec_input <- !is.null(sim_env$Rec_input) && y <= dim(sim_env$Rec_input)[3]
    if(use_rec_input) Rec_y <- array(sim_env$Rec_input[,,y,sim], dim = c(n_pop, n_regions))
    else Rec_y <- det_rec * exp(array(sim_env$ln_RecDevs[,,y,sim], dim = c(n_pop, n_regions)))

    for(p in 1:n_pop) {
      for(r in 1:n_regions) {
        for(s in 1:n_sexes) {
          N[p,r,1,s] <- Rec_y[p,r] * sim_env$rec_seas_prop[p,seas,sim] * sim_env$sexratio[p,r,y,s,sim]
          N_fished[p,r,1,s] <- N[p,r,1,s]
        } # end s loop
      } # end r loop
    } # end p loop

    # recruits that move from birth missed the move at the start of the season
    if(n_regions > 1 && move_timing == 0 && sim_env$do_recruits_move == 1) {
      for(p in 1:n_pop) {
        for(s in 1:n_sexes) N_fished[p,,1,s] <- as.vector(t(N_fished[p,,1,s]) %*% sim_env$Movement[p,,,y,seas,1,s,sim])
      } # end p loop
    }

  } # end if recruits from this year's spawning

  # numbers removable over the season. under move_timing 2 fish move while dying, so integrate over the season
  if(move_timing == 2) {
    Avail <- array(0, dim = c(n_pop, n_regions, n_ages, n_sexes))
    for(p in 1:n_pop) {
      for(a in 1:n_ages) {
        for(s in 1:n_sexes) {
          Avail[p,,a,s] <- integrate_seas_abundance(N[p,,a,s],
                                                    ZAA[p,,a,s],
                                                    sim_env$Mrate[p,,,y,seas,a,s,sim],
                                                    seasdur,
                                                    expm_nsub = expm_nsub)
        } # end s loop
      } # end a loop
    } # end p loop
  } else {
    Avail <- N_fished * (1 - exp(-ZAA)) / ZAA # Baranov
  }

  # retained catch by region and fleet, summed over populations, ages and sexes
  catch <- array(0, dim = c(n_regions, n_fish_fleets))
  for(r in 1:n_regions) {
    for(f in 1:n_fish_fleets) {
      for(p in 1:n_pop) {

        CAA <- F_seas[r,f] * sim_env$fish_sel[p,r,y,seas,,,f,sim] * sim_env$ret_sel[p,r,y,seas,,,f,sim] * Avail[p,r,,] # retained catch at age
        wt <- if(target_units[f] == 1) sim_env$WAA_fish[p,r,y,seas,,,f,sim] else 1 # biomass, or numbers
        catch[r,f] <- catch[r,f] + sum(CAA * wt)

      } # end p loop
    } # end f loop
  } # end r loop

  # survivors at the end of the season
  if(move_timing == 0 || n_regions == 1) {
    N_end <- N_fished * exp(-ZAA)
  } else {
    N_end <- array(0, dim = c(n_pop, n_regions, n_ages, n_sexes))
    for(p in 1:n_pop) {
      for(a in 1:n_ages) {

        moves <- (sim_env$do_recruits_move == 1 || a > 1) # recruits stay put unless they move from birth
        for(s in 1:n_sexes) {
          Mv <- if(moves) sim_env$Movement[p,,,y,seas,a,s,sim] else diag(n_regions)
          Qv <- if(moves) sim_env$Mrate[p,,,y,seas,a,s,sim] else matrix(0, n_regions, n_regions)
          N_end[p,,a,s] <- advance_seas(N[p,,a,s],
                                        Mv,
                                        ZAA[p,,a,s],
                                        Qv,
                                        seasdur,
                                        move_timing,
                                        expm_nsub = expm_nsub)
        } # end s loop

      } # end a loop
    } # end p loop
  } # end if movement at the end of the season or continuous

  return(list(catch = catch, N_end = N_end, Rec_y = Rec_y))

} # end function


#' Solve One Season's F For A Block Of Regions
#'
#' Tapes log catch in the cells being solved as a function of their log F with
#' \code{RTMB::MakeTape} and solves them jointly with \code{solve_log_catch_newton},
#' since catch in a cell falls as other fleets in the same region fish harder.
#'
#' @param F_seas Numeric matrix \code{[n_regions, n_fish_fleets]}. This season's
#'   F so far; cells outside \code{free} are kept fixed.
#' @param free Integer vector. Cells to solve, counted down regions first.
#' @param target_seas Numeric matrix \code{[n_regions, n_fish_fleets]}.
#' @param N,y,seas,sim,sim_env,target_units,spawn_recruits Passed to \code{om_season_catch}.
#' @param catch_f_max,catch_tol,catch_max_iter Solver settings.
#'
#' @return Named list with \code{F_seas}, solved in the \code{free} cells, and
#'   \code{capped}, the cells left at \code{catch_f_max} because the target is
#'   not reachable there.
#' @keywords internal
solve_om_season_F <- function(F_seas, free, target_seas, N, y, seas, sim, sim_env, target_units, spawn_recruits,
                              catch_f_max, catch_tol, catch_max_iter) {

  n_regions <- sim_env$n_regions

  # log catch in the free cells as a function of their log F, with F everywhere else kept fixed
  log_catch_tape <- RTMB::MakeTape(function(theta) {
    "[<-" <- RTMB::ADoverload("[<-")
    F_trial <- F_seas
    F_trial[free] <- exp(theta)
    log(om_season_catch(N, F_trial, y, seas, sim, sim_env, target_units, spawn_recruits)$catch[free])
  }, rep(log(0.1), length(free)))

  # start every cell at F = 0.1
  newton <- solve_log_catch_newton(log_catch_tape = log_catch_tape,
                                   theta_start = rep(log(0.1), length(free)),
                                   log_target = log(target_seas[free]),
                                   theta_max = log(catch_f_max),
                                   catch_tol = catch_tol,
                                   catch_max_iter = catch_max_iter)
  F_seas[free] <- exp(newton$theta)
  capped <- free[newton$capped_pos]

  if(length(capped) > 0) {
    cell_names <- paste0("region ", (capped - 1) %% n_regions + 1, " fleet ", (capped - 1) %/% n_regions + 1)
    warning(paste0("Catch target in year ", y, ", season ", seas, " is not reachable for ",
                   paste(cell_names, collapse = ", "), " at the F bound catch_f_max = ",
                   catch_f_max, ". F is capped there and the target is undershot."))
  }

  return(list(F_seas = F_seas, capped = capped))

} # end function


#' Convert Catch Advice To Fishing Mortality In The Operating Model
#'
#' Finds the F by region, season and fleet that makes the operating model's
#' retained catch in year \code{y} equal \code{target}, replaying the year one
#' season at a time with the operating model's own dynamics. Fleets fishing the
#' same region share total mortality, what earlier seasons caught is gone before
#' later ones are fished, fish move between regions as the operating model moves
#' them, and under \code{rec_lag = 0} this year's recruits come from the
#' spawning biomass the trial F leaves.
#'
#' @param target Numeric array \code{[n_regions, n_seas, n_fish_fleets]}. Retained
#'   catch for year \code{y}, in \code{target_units}. Zero means no fishing.
#' @param y Integer. Year being fished. Its start of year numbers at age must
#'   exist, so call after \code{run_annual_cycle(y - 1, ...)}.
#' @param sim Integer. Simulation replicate.
#' @param sim_env Simulation environment from \code{\link{Setup_sim_env}}. Under
#'   cohort growth, year \code{y}'s growth is formed here, from the same numbers
#'   the annual cycle forms it from at the start of the year.
#' @param target_units Units of \code{target} by fleet, \code{"biom"} (default,
#'   through \code{WAA_fish}) or \code{"abd"} for numbers, or 1 and 0. One value
#'   is used for every fleet. Independent of the fleet's \code{catch_units}.
#' @param catch_f_max,catch_tol,catch_max_iter Upper bound on F, relative catch
#'   tolerance and iterations per solve, as in \code{\link{Do_Population_Projection}}.
#'
#' @return Named list with \code{Fmort} \code{[n_regions, n_seas, n_fish_fleets]},
#'   to write into \code{sim_env$Fmort[, y, , , sim]}, and \code{resid}, the
#'   relative catch miss in each cell.
#'
#' @export catch_to_F_om
#' @family Closed Loop Simulations
catch_to_F_om <- function(target,
                          y,
                          sim,
                          sim_env,
                          target_units = "biom",
                          catch_f_max = 5,
                          catch_tol = 1e-6,
                          catch_max_iter = 100) {

  # get dimensions
  n_pop <- sim_env$n_pop
  n_regions <- sim_env$n_regions
  n_seas <- sim_env$n_seas
  n_ages <- sim_env$n_ages
  n_sexes <- sim_env$n_sexes
  n_fish_fleets <- sim_env$n_fish_fleets
  move_timing <- if(is.null(sim_env$move_timing)) 0 else sim_env$move_timing
  rec_lag <- sim_env$rec_lag
  spawn_seas <- sim_env$spawn_seas

  target <- array(target, dim = c(n_regions, n_seas, n_fish_fleets))
  target_units <- convert_to_numeric(target_units, list(abd = 0, biom = 1)) # same codes as catch_units
  if(length(target_units) == 1) target_units <- rep(target_units, n_fish_fleets) # one setting for every fleet
  Fmort = resid = array(0, dim = c(n_regions, n_seas, n_fish_fleets))

  # cohort growth for year y comes from its start of year numbers, which exist already. the annual
  # cycle forms it again from the same numbers, so forming it early changes nothing
  if(!is.null(sim_env$growth_state) && y >= sim_growth_args(sim_env)$growth_cohort_styr) advance_sim_growth_year(y, sim, sim_env)

  # start of year numbers at age, with this year's recruits already in under rec_lag > 0
  N <- array(sim_env$NAA[,,y,1,,,sim], dim = c(n_pop, n_regions, n_ages, n_sexes))
  Rec_y <- array(sim_env$Rec[,,y,sim], dim = c(n_pop, n_regions)) # replaced at spawn_seas under rec_lag 0

  for(seas in 1:n_seas) {

    # recruits entering this season after their first, as apply_pop_dy adds them
    if(if(rec_lag != 0) seas > 1 else seas > spawn_seas) {
      for(p in 1:n_pop) {
        for(r in 1:n_regions) {
          for(s in 1:n_sexes) {
            N[p,r,1,s] <- N[p,r,1,s] + Rec_y[p,r] * sim_env$rec_seas_prop[p,seas,sim] * sim_env$sexratio[p,r,y,s,sim]
          } # end s loop
        } # end r loop
      } # end p loop
    } # end if seasonal recruitment
    spawn_recruits <- rec_lag == 0 && seas == spawn_seas # recruits from this season's own spawning

    # regions solved one at a time, or together when fish move between them while being caught, or when
    # this season's recruits come from every region's spawning biomass
    if(n_regions > 1 && (move_timing == 2 || spawn_recruits)) region_blocks <- list(1:n_regions) else region_blocks <- as.list(1:n_regions)

    # earlier seasons are settled, so only this season's F is solved
    target_seas <- array(target[,seas,], dim = c(n_regions, n_fish_fleets))
    F_seas <- array(0, dim = c(n_regions, n_fish_fleets))
    capped <- c()
    for(block in region_blocks) {

      in_block <- array(FALSE, dim = c(n_regions, n_fish_fleets))
      in_block[block,] <- TRUE
      free <- which(in_block & target_seas > 0) # a zero target is F = 0, not something to solve
      if(length(free) == 0) next

      block_solve <- solve_om_season_F(F_seas = F_seas,
                                       free = free,
                                       target_seas = target_seas,
                                       N = N,
                                       y = y,
                                       seas = seas,
                                       sim = sim,
                                       sim_env = sim_env,
                                       target_units = target_units,
                                       spawn_recruits = spawn_recruits,
                                       catch_f_max = catch_f_max,
                                       catch_tol = catch_tol,
                                       catch_max_iter = catch_max_iter)
      F_seas <- block_solve$F_seas
      capped <- c(capped, block_solve$capped)

    } # end block loop
    Fmort[,seas,] <- F_seas

    # realized catch at the solved F, and the survivors going into next season
    seas_out <- om_season_catch(N, F_seas, y, seas, sim, sim_env, target_units, spawn_recruits)
    if(spawn_recruits) Rec_y <- seas_out$Rec_y # this year's recruitment at the solved F
    resid_seas <- array(0, dim = c(n_regions, n_fish_fleets))
    pos <- target_seas > 0
    resid_seas[pos] <- (seas_out$catch[pos] - target_seas[pos]) / target_seas[pos]
    resid[,seas,] <- resid_seas

    missed <- setdiff(which(pos & abs(resid_seas) > max(catch_tol, 1e-4)), capped)
    if(length(missed) > 0) {
      warning(paste0("Catch target in year ", y, ", season ", seas, " did not converge. Largest relative catch error is ",
                     signif(max(abs(resid_seas[missed])), 3), "."))
    }

    N <- seas_out$N_end
    # state-space numbers at age at the within-year boundary, drawn ahead of time in year 1
    if(seas < n_seas && isTRUE(sim_env$NAA_re > 0)) {
      N <- N * exp(array(sim_env$naa_eta[,,y,seas + 1,,], dim = c(n_pop, n_regions, n_ages, n_sexes)))
    }

  } # end seas loop

  return(list(Fmort = Fmort, resid = resid))

} # end function
