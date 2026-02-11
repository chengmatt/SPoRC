#' Set up simulation list for closed-loop projections
#'
#' This function creates and initializes a simulation list for closed-loop
#' projections of population dynamics. All components of the simulation list
#' must match the expected names used internally by the setup functions.
#' Users can provide **custom definitions** for any component by passing them
#' through `...` using the correct name (e.g., `WAA_input`, `fish_sel_input`).
#'
#' @param closed_loop_yrs Integer. Number of years to project in the closed loop.
#' @param n_sims Integer. Number of simulation replicates.
#' @param data List. Observed data and configuration for the population.
#' @param parameters List. Parameter values for the model.
#' @param mapping List. Mapping of parameters for optimization.
#' @param sd_rep List. Standard deviation reports from fitted model.
#' @param ... Optional named arguments for custom inputs.
#'   Each name must correspond to a component expected by the simulation setup functions, and be dimensioned appropriately:
#'   `Setup_Sim_Fishing()`, `Setup_Sim_Survey()`, `Setup_Sim_Biologicals()`,
#'   `Setup_Sim_Rec()`, and `Setup_Sim_Tagging()`. Examples include:
#'
#'   - **Fishing**: `fish_sel_input`, `ln_sigmaC`, `Fmort_input`, `fish_q_input`, etc.
#'   - **Survey**: `srv_sel_input`, `srv_q_input`, `ObsSrvIdx_SE`, etc.
#'   - **Biologicals**: `WAA_input`, `MatAA_input`, `natmort_input`, `AgeingError_input`, `SizeAgeTrans_input`
#'   - **Recruitment inputs**:
#'       - `R0_input`, `h_input`, `sexratio_input`, `ln_InitDevs_input`, `Rec_input`
#'       - `Rec_input`:
#'           - If shorter than the number of projection years, new recruitment deviates will be simulated
#'             based on `recruitment_opt`, `R0_input`, and `h_input` (supports changing regimes across years).
#'           - If you want fixed recruitment for all projection years, provide a `Rec_input` array that spans all years.
#'   - **Tagging**: `Tag_Reporting_input`, `ln_Init_Tag_Mort`, `ln_Tag_Shed`, `tag_selex`, `tag_natmort`, `n_tags`, `n_tags_rel_input`
#'   - **Movement**: `Movement` (must match the expected dimensions and be named exactly `Movement`)
#'
#'   The values must have the correct dimensions expected by each component.
#'   If a component is not provided, default behavior will extend the last year (or zeros for fishing mortality, which can filled in subsequently)
#' @param random Character vector of random effects estimated
#' @param rep List. Report from fitted model.
#' @param FishIdx_SE_fill Character or numeric. Fill method for fishery index
#'   standard errors when extending to simulation years. Options are:
#'   - `"zeros"`: fill with zeros
#'   - `"last"`: repeat last non-NA slice
#'   - `"mean"`: fill with the mean of the observed series
#'   - Numeric: constant scalar or array value
#' @param SrvIdx_SE_fill Character or numeric. Fill method for survey index
#'   standard errors. Same options as `FishIdx_SE_fill`.
#' @param ISS_FishAgeComps_fill Character or numeric. Fill method for fishery
#'   age composition input sample sizes. Options are:
#'   - `"zeros"`, `"last"`, `"mean"` (as above)
#'   - `"F_pattern"`: fill based on fishing mortality pattern in the closed-loop simulation
#'   - Numeric: constant scalar or array value
#' @param ISS_FishLenComps_fill Character or numeric. Fill method for fishery
#'   length composition input sample sizes. Same options as `ISS_FishAgeComps_fill`.
#' @param ISS_SrvAgeComps_fill Character or numeric. Fill method for survey
#'   age composition input sample sizes. Options are `"zeros"`, `"last"`, `"mean"`,
#'   or a numeric constant.
#' @param ISS_SrvLenComps_fill Character or numeric. Fill method for survey
#'   length composition input sample sizes. Same options as `ISS_SrvAgeComps_fill`.
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
                                              ISS_FishAgeComps_fill = "mean",
                                              ISS_FishLenComps_fill = "mean",
                                              ISS_SrvAgeComps_fill = "mean",
                                              ISS_SrvLenComps_fill = "mean",
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

  optim_parameters_list <- get_optim_param_list(parameters, mapping, sd_rep, random) # get optimized parameters in original list format

  # Setup Model Dimensions --------------------------------------------------
  sim_list <- Setup_Sim_Dim(n_sims = n_sims, # number of simulations
                            n_yrs = length(data$years) + closed_loop_yrs, # number of years
                            n_regions = data$n_regions,  # number of regions
                            n_ages = length(data$ages), # number of ages
                            n_obs_ages = dim(data$ObsFishAgeComps)[3], # number of observed ages
                            n_lens = length(data$lens), # number of lengths
                            n_sexes = data$n_sexes, # number of sexes
                            n_fish_fleets = data$n_fish_fleets, # number of fishery fleets
                            n_srv_fleets = data$n_srv_fleets, # number of survey fleets
                            feedback_start_yr = length(data$years), # when to start closed loop feedback
                            n_seas = data$n_seas, # number of seasons
                            seasdur = data$seasdur, # seasonal duration
                            run_feedback = TRUE # whether or not to run feedback (closed loop)
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

  # Fishery selectivity
  fish_sel_input <- if(!"fish_sel_input" %in% names(args)) {
    extend_years(replicate(n = sim_list$n_sims, rep$fish_sel[,1:length(data$years),,,,drop = FALSE]), n_years = closed_loop_yrs, 2, fill = 'last')
  } else args$fish_sel_input
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

  # setup fishery simulation processes
  sim_list <- Setup_Sim_Fishing(
    sim_list = sim_list, # update simulate list
    ln_sigmaC = ln_sigmaC,
    Fmort_input = extend_years(replicate(n = sim_list$n_sims, rep$Fmort[,1:length(data$years),,drop = FALSE]), n_years = closed_loop_yrs, 2, fill = 'zeros'),
    fish_sel_input = fish_sel_input,
    fish_q_input = fish_q_input,
    ObsFishIdx_SE = ObsFishIdx_SE,
    fish_idx_type = data$fish_idx_type,
    init_F_val = rep$init_F,

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
    FishLen_corr_pars = FishLen_corr_pars
  )

  # add in ISS F pattern into simulation list
  if(ISS_FishAgeComps_fill == 'F_pattern') sim_list$ISS_FishAgeComps_fill <- "F_pattern"
  if(ISS_FishLenComps_fill == 'F_pattern') sim_list$ISS_FishLenComps_fill <- "F_pattern"

  # Setup Survey Processes --------------------------------------------------
  # Survey selectivity
  srv_sel_input <- if(!"srv_sel_input" %in% names(args)) {
    extend_years(replicate(n = sim_list$n_sims, rep$srv_sel[,1:length(data$years),,,,drop = FALSE]), closed_loop_yrs, 2, 'last')
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

  # setup survey simulation processes
  sim_list <- Setup_Sim_Survey(
    sim_list = sim_list,
    srv_sel_input = srv_sel_input,
    srv_q_input = srv_q_input,
    ObsSrvIdx_SE = ObsSrvIdx_SE,
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
    SrvLen_corr_pars = SrvLen_corr_pars
  )

  # Setup Biological Dynamics -----------------------------------------------
  natmort_input <- if(!"natmort_input" %in% names(args)) {
    extend_years(replicate(n = sim_list$n_sims, rep$natmort[,1:length(data$years),,,drop = FALSE]), closed_loop_yrs, 2, 'last')
  } else args$natmort_input
  WAA_input <- if(!"WAA_input" %in% names(args)) {
    extend_years(replicate(n = sim_list$n_sims, data$WAA[,1:length(data$years),,,,drop = FALSE]), closed_loop_yrs, 2, 'last')
  } else args$WAA_input
  WAA_fish_input <- if(!"WAA_fish_input" %in% names(args)) {
    extend_years(replicate(n = sim_list$n_sims, data$WAA_fish[,1:length(data$years),,,,,drop = FALSE]), closed_loop_yrs, 2, 'last')
  } else args$WAA_fish_input
  WAA_srv_input <- if(!"WAA_srv_input" %in% names(args)) {
    extend_years(replicate(n = sim_list$n_sims, data$WAA_srv[,1:length(data$years),,,,,drop = FALSE]), closed_loop_yrs, 2, 'last')
  } else args$WAA_srv_input
  MatAA_input <- if(!"MatAA_input" %in% names(args)) {
    extend_years(replicate(n = sim_list$n_sims, data$MatAA[,1:length(data$years),,,,drop = FALSE]), closed_loop_yrs, 2, 'last')
  } else args$MatAA_input
  AgeingError_input <- if(!"AgeingError_input" %in% names(args)) {
    extend_years(replicate(n = sim_list$n_sims, data$AgeingError[1:length(data$years),,,drop = FALSE]), closed_loop_yrs, 1, 'last')
  } else args$AgeingError_input
  SizeAgeTrans_input <- if(!"SizeAgeTrans_input" %in% names(args)) {
    extend_years(replicate(n = sim_list$n_sims, data$SizeAgeTrans[,1:length(data$years),,,,,drop = FALSE]), closed_loop_yrs, 2, 'last')
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
    SizeAgeTrans_input = SizeAgeTrans_input
  )

  # Setup Recruitment Processes ---------------------------------------------
  h_input <- if(!"h_input" %in% names(args)) {
    replicate(n = sim_list$n_sims, array(rep$h_trans, dim = c(sim_list$n_regions, sim_list$n_yrs)))
  } else args$h_input
  R0_input <- if(!"R0_input" %in% names(args)) {
    replicate(n = sim_list$n_sims, expr = array(rep$R0 * rep$Rec_trans_prop, dim = c(sim_list$n_regions, sim_list$n_yrs)))
  } else args$R0_input
  sexratio_input <- if(!"sexratio_input" %in% names(args)) {
    extend_years(replicate(n = sim_list$n_sims, expr = rep$sexratio[,1:length(data$years),,drop = FALSE]), closed_loop_yrs, 2, 'last')
  } else args$sexratio_input
  ln_sigmaR <- if(!"ln_sigmaR" %in% names(args)) optim_parameters_list$ln_sigmaR else args$ln_sigmaR
  Rec_input <- if(!"Rec_input" %in% names(args)) {
    replicate(n = sim_list$n_sims, expr = rep$Rec[,1:length(data$years),drop = FALSE])
  } else args$Rec_input
  ln_InitDevs_input <- if(!"ln_InitDevs_input" %in% names(args)) {
    replicate(sim_list$n_sims, optim_parameters_list$ln_InitDevs)
  } else args$ln_InitDevs_input

  # recruitment options
  rec_dd <- if(!"rec_dd" %in% names(args)) "global" else args$rec_dd
  init_dd <- if(!"init_dd" %in% names(args)) "global" else args$init_dd
  rec_lag <- if(!"rec_lag" %in% names(args)) 1 else args$rec_lag
  recruitment_opt <- if(!"recruitment_opt" %in% names(args)) "bh_rec" else args$recruitment_opt

  # setup recruitment simulation
  sim_list <- Setup_Sim_Rec(
    sim_list = sim_list,
    spawn_eas = data$spawn_seas,
    do_recruits_move = data$do_recruits_move,
    t_spawn = data$t_spawn,
    init_age_strc = data$init_age_strc,
    h_input = h_input,
    R0_input = R0_input,
    sexratio_input = sexratio_input,
    ln_sigmaR = ln_sigmaR,
    Rec_input = Rec_input,
    ln_InitDevs_input = ln_InitDevs_input,
    recruitment_opt = recruitment_opt,
    rec_dd = rec_dd,
    init_dd = init_dd,
    rec_lag = rec_lag
  )

  # Setup Tagging -----------------------------------------------------------
  n_tags_rel_input <- if(!"n_tags_rel_input" %in% names(args)) {
    if(!is.na(sum(data$Tagged_Fish))) apply(data$Tagged_Fish, 1, sum) else NA
  } else args$n_tags_rel_input
  n_tags <- if(!"n_tags" %in% names(args)) NULL else args$n_tags
  tag_release_indicator <- if(!"tag_release_indicator" %in% names(args)) {
    if(exists("tag_release_indicator", data)) data$tag_release_indicator else NA
  } else args$tag_release_indicator
  ln_Init_Tag_Mort <- if(!"ln_Init_Tag_Mort" %in% names(args)) optim_parameters_list$ln_Init_Tag_Mort else args$ln_Init_Tag_Mort
  ln_Tag_Shed <- if(!"ln_Tag_Shed" %in% names(args)) optim_parameters_list$ln_Tag_Shed else args$ln_Tag_Shed
  Tag_Reporting_input <- if(!"Tag_Reporting_input" %in% names(args)) {
    if(is.null(rep$Tag_Reporting)) NULL else extend_years(replicate(n = sim_list$n_sims, rep$Tag_Reporting), closed_loop_yrs, 2, 'last')
  } else args$Tag_Reporting_input
  UseTagging <- if(!"UseTagging" %in% names(args)) data$UseTagging else args$UseTagging
  tag_selex <- if(!"tag_selex" %in% names(args)) data$tag_selex else args$tag_selex
  tag_natmort <- if(!"tag_natmort" %in% names(args)) data$tag_natmort else args$tag_natmort
  tag_like <- if(!"tag_like" %in% names(args)) data$Tag_LikeType else args$tag_like
  ln_tag_theta <- if(!"ln_tag_theta" %in% names(args)) parameters$ln_tag_theta else args$ln_tag_theta

  # setup tagging simulation
  if(!is.null(n_tags)) sim_list$n_tags_rel_input <- NULL # set release input to NULL if n_tags is specified.
  sim_list <- Setup_Sim_Tagging(
    sim_list = sim_list,
    max_liberty = data$max_tag_liberty,
    t_tagging = data$t_tagging,
    n_tags = n_tags,
    n_tags_rel_input = n_tags_rel_input * data$Wt_Tagging,
    tag_release_indicator = tag_release_indicator,
    ln_Init_Tag_Mort = ln_Init_Tag_Mort,
    ln_Tag_Shed = ln_Tag_Shed,
    Tag_Reporting_input = Tag_Reporting_input,
    UseTagging = UseTagging,
    tag_selex = tag_selex,
    tag_natmort = tag_natmort,
    tag_like = tag_like,
    ln_tag_theta = ln_tag_theta
  )

  # Movement ----------------------------------------------------------------
  Movement <- if(!"Movement" %in% names(args)) extend_years(replicate(n = sim_list$n_sims, rep$Movement[,,1:length(data$years),,,,drop = FALSE]), closed_loop_yrs, 3, 'last') else args$Movement
  sim_list$Movement <- Movement

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
#'     \item{what}{Method for reference point calculation. Options include "SPR", "BH_MSY", "independent_SPR", "independent_BH_MSY", "global_SPR", "global_BH_MSY". Default is "SPR".}
#'   }
#' @param sim_env Simulation environment
#' @param n_proj_yrs Number of projection years
#' @param t_spawn Spawn timing wihtin a given season / year, default uses sim_env true values
#'
#' @return A list with elements:
#'   \describe{
#'     \item{f_ref_pt}{Array of fishing reference points by region and projection year.}
#'     \item{b_ref_pt}{Array of biological reference points by region and projection year.}
#'     \item{virgin_b_ref_pt}{Array of unfished biological reference points by region and projection year.}
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
                                               what = "SPR"
                                             ),
                                             n_proj_yrs
                                             ) {

  if(use_true_values) {
    # Build data and report objects to feed into reference points
    data_obj <- list(
      ages = 1:sim_env$n_ages,
      years = 1:y,
      n_seas = sim_env$n_seas,
      seasdur = sim_env$seasdur,
      spawn_seas = sim_env$spawn_seas,
      n_fish_fleets = sim_env$n_fish_fleets,
      n_regions = sim_env$n_regions,
      WAA = array(sim_env$WAA[, 1:y, ,, , sim], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_ages, sim_env$n_sexes)),
      MatAA = array(sim_env$MatAA[, 1:y,, , , sim], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_ages, sim_env$n_sexes)),
      do_recruits_move = sim_env$do_recruits_move
    )

    # Build rep list if not doing assessment (using truth)
    rep_obj <- list(
      Fmort = array(sim_env$Fmort[, 1:y,, , sim], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets)),
      fish_sel = array(sim_env$fish_sel[, 1:y, , , , sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_ages, sim_env$n_sexes, sim_env$n_fish_fleets)),
      natmort = array(sim_env$natmort[, 1:y, , , sim], dim = c(sim_env$n_regions, length(1:y), sim_env$n_ages, sim_env$n_sexes)),
      h_trans = sim_env$h[, y, sim],
      R0 = sum(sim_env$R0[, y, sim]),
      Rec_trans_prop = sim_env$R0[, y, sim] / sum(sim_env$R0[, y, sim]),
      Rec = array(sim_env$Rec[, 1:y, sim], dim = c(sim_env$n_regions, length(1:y))),
      Movement = array(sim_env$Movement[, , 1:y, , , , sim],  dim = c(sim_env$n_regions, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_ages, sim_env$n_sexes))
    )

    # get sex ratio
    tmp_sex_ratio_f <- if(sim_env$n_sexes == 1) rep(0.5, sim_env$n_regions) else sim_env$sexratio[,y,1,sim]

  } else {
    data_obj <- asmt_data
    rep_obj <- asmt_rep
    tmp_sex_ratio_f <- if(data_obj$n_sexes == 1) rep(0.5, sim_env$n_regions) else rep_obj$sexratio[,y,1]
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
                                           n_avg_yrs = reference_points_opt$n_avg_yrs
  )

  # extract fishery and biological reference points
  f_ref_pt <- array(reference_points$f_ref_pt, dim = c(data_obj$n_regions, n_proj_yrs)) # fishery reference points
  b_ref_pt <- array(reference_points$b_ref_pt, dim = c(data_obj$n_regions, n_proj_yrs)) # biological reference points
  virgin_b_ref_pt <- array(reference_points$virgin_b_ref_pt, dim = c(data_obj$n_regions, n_proj_yrs)) # biological reference points

  return(list(f_ref_pt = f_ref_pt, b_ref_pt = b_ref_pt, virgin_b_ref_pt = virgin_b_ref_pt))
}


