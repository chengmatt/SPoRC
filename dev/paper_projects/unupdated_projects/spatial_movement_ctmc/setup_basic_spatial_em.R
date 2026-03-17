# Function to setup basic spatial simualtion EM
setup_basic_spatial_sim_em <- function(sim_obj, sim) {

  # Extract simulation data for current year and replicate
  sim_data <- simulation_data_to_SPoRC(sim_env = sim_obj, y = sim_obj$n_years, sim = sim)

  # Setup model dimensions
  input_list <- Setup_Mod_Dim(
    years = 1:sim_obj$n_years,
    ages = 1:sim_obj$n_ages,
    lens = sim_obj$n_lens,
    n_regions = sim_obj$n_regions,
    n_sexes = sim_obj$n_sexes,
    n_fish_fleets = sim_obj$n_fish_fleets,
    n_srv_fleets = sim_obj$n_srv_fleets,
    n_seas = sim_obj$n_seas,
    verbose = F
  )

  # Recruitment setup
  input_list <- Setup_Mod_Rec(
    input_list = input_list,
    do_rec_bias_ramp = 0, # not doing bias ramp
    sigmaR_switch = 1, # when to switch from early to late sigmaR (switch in first year)
    ln_sigmaR = rep(log(1) , 2), # 2 values for early and late sigma
    rec_model = "mean_rec",
    sigmaR_spec = "fix", # fix early sigmaR and late sigmaR
    init_age_strc = 2, # scalar geometric series to derive initial age structure
    equil_init_age_strc = 2, # estimating all intial age deviations
    ln_global_R0 = log(20)
  )

  # Biological setup
  input_list <- Setup_Mod_Biologicals(
    input_list = input_list,
    # Data inputs
    WAA = sim_data$WAA,
    MatAA = sim_data$MatAA,
    WAA_fish = sim_data$WAA_fish,
    WAA_srv = sim_data$WAA_srv,
    fit_lengths = 0, # not fitting lengths
    AgeingError = sim_data$AgeingError,
    M_spec = "fix",     # fixing natural mortality
    Fixed_natmort = array(0.2, dim = c(input_list$data$n_regions, length(input_list$data$years),
                                       length(input_list$data$ages), input_list$data$n_sexes))
  )

  # Tagging setup
  input_list <- Setup_Mod_Tagging(input_list = input_list,
                                  UseTagging = 1,
                                  tag_release_indicator = sim_data$tag_release_indicator,
                                  max_tag_liberty = 15,
                                  Tagged_Fish = sim_data$Tagged_Fish,
                                  Obs_Tag_Recap = sim_data$Obs_Tag_Recap,
                                  Tag_LikeType = "Multinomial_Release",
                                  t_tagging = 0.5,
                                  mixing_period = 1,
                                  tag_selex = "SexSp_AllFleet",
                                  tag_natmort = "AgeSp_SexSp",
                                  ln_Init_Tag_Mort = log(1e-100),
                                  ln_Tag_Shed = log(1e-100),
                                  TagRep_spec = "est_shared_r",
                                  Init_Tag_Mort_spec = "fix",
                                  Tag_Shed_spec = "fix",
                                  move_age_tag_pool = as.list(1:15),
                                  move_sex_tag_pool = as.list(1)
  )

  # Fishery catch & fishing mortality
  input_list <- Setup_Mod_Catch_and_F(
    input_list = input_list,
    # Data inputs
    ObsCatch = sim_data$ObsCatch,
    UseCatch = sim_data$UseCatch,
    # Model options
    Use_F_pen = 1,
    sigmaC_spec = "fix",
    # Fixing sigma C and F
    ln_sigmaC = sim_data$ln_sigmaC,
    ln_sigmaF = array(log(3), dim = c(input_list$data$n_regions, input_list$data$n_seas, input_list$data$n_fish_fleets))
  )

  # Survey selectivity and catchability
  input_list <- Setup_Mod_FishIdx_and_Comps(
    input_list = input_list,
    # Data inputs
    ObsFishIdx = sim_data$ObsFishIdx,
    ObsFishIdx_SE = sim_data$ObsFishIdx_SE,
    UseFishIdx = array(0, dim = c(sim_obj$n_regions, sim_obj$n_years, input_list$data$n_seas,  sim_obj$n_fish_fleets)),
    ObsFishAgeComps = sim_data$ObsFishAgeComps,
    ObsFishLenComps = sim_data$ObsFishLenComps,
    UseFishAgeComps = sim_data$UseFishAgeComps,
    UseFishLenComps = sim_data$UseFishLenComps,
    ISS_FishAgeComps = sim_data$ISS_FishAgeComps,
    ISS_FishLenComps = sim_data$ISS_FishLenComps,
    # Model options
    fish_idx_type = c("none"),
    FishAgeComps_LikeType = c("Multinomial"),
    FishLenComps_LikeType = c("none"),
    FishAgeComps_Type = c("spltRjntS_Year_1-terminal_Fleet_1"),
    FishLenComps_Type = c("none_Year_1-terminal_Fleet_1")
  )

  # Survey indices and compositions
  input_list <- Setup_Mod_SrvIdx_and_Comps(
    input_list = input_list,
    # Data inputs
    ObsSrvIdx = sim_data$ObsSrvIdx,
    ObsSrvIdx_SE = sim_data$ObsSrvIdx_SE,
    UseSrvIdx = sim_data$UseSrvIdx,
    ObsSrvAgeComps = sim_data$ObsSrvAgeComps,
    ObsSrvLenComps = sim_data$ObsSrvLenComps,
    UseSrvAgeComps = sim_data$UseSrvAgeComps,
    UseSrvLenComps = sim_data$UseSrvLenComps,
    ISS_SrvAgeComps = sim_data$ISS_SrvAgeComps,
    ISS_SrvLenComps = sim_data$ISS_SrvLenComps,
    # Model options
    srv_idx_type = c("biom"),
    SrvAgeComps_LikeType = c("Multinomial"),
    SrvLenComps_LikeType = c("none"),
    SrvAgeComps_Type = c("spltRjntS_Year_1-terminal_Fleet_1"),
    SrvLenComps_Type = c("none_Year_1-terminal_Fleet_1")
  )

  # Fishery selectivity and catchability
  input_list <- Setup_Mod_Fishsel_and_Q(
    input_list = input_list,
    # Model options
    fish_sel_model = c("logist1_Fleet_1"), # fishery selex model (NOTE: ASSUMES DOMED)
    fish_fixed_sel_pars_spec = c("est_shared_r"), # whether to estiamte all fixed effects for fishery selectivity
    fish_q_spec = "fix" # estimate fishery q
  )

  # Survey selectivity and catchability
  input_list <- Setup_Mod_Srvsel_and_Q(
    input_list = input_list,
    # Model options
    srv_sel_model = c("logist1_Fleet_1"), # survey selectivity form
    srv_fixed_sel_pars_spec = c("est_shared_r"), # whether to estimate all fixed effects for survey selectivity
    srv_q_spec = c("est_shared_r")  # whether to estiamte all fixed effects for survey catchability
  )

  # Data weighting
  input_list <- Setup_Mod_Weighting(
    input_list = input_list,
    Wt_Catch = 1,
    Wt_FishIdx = 1,
    Wt_SrvIdx = 1,
    Wt_Rec = 1,
    Wt_F = 1,
    Wt_Tagging = 1,
    Wt_FishAgeComps = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas,
                                       input_list$data$n_sexes, input_list$data$n_fish_fleets)),
    Wt_FishLenComps = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas,
                                       input_list$data$n_sexes, input_list$data$n_fish_fleets)),
    Wt_SrvAgeComps = array(1, dim = c(input_list$data$n_regions,length(input_list$data$years), input_list$data$n_seas,
                                      input_list$data$n_sexes, input_list$data$n_srv_fleets)),
    Wt_SrvLenComps = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas,
                                      input_list$data$n_sexes, input_list$data$n_srv_fleets))
  )

  return(input_list)
}

# Function to run one single EM
fit_single_model <- function(row_idx, sim_obj, model_df, ctmc_data, adjacency,
                             diffusion_formula, preference_formula) {

  # Extract parameters for this run
  sim_id <- param_grid$sim_id[row_idx]
  model_id <- param_grid$model_id[row_idx]

  tryCatch({
    # Make input list
    input_list <- setup_basic_spatial_sim_em(sim_obj = sim_obj, sim = sim_id)

    # Setup movement
    input_list <- Setup_Mod_Movement(
      input_list = input_list,
      do_recruits_move = 0,
      ctmc_move_dat = ctmc_data,
      adjacency_mat = adjacency,
      area_r = rep(1, 3), # equal areas

      # formulas are set within a given trial
      diffusion_formula = diffusion_formula,
      preference_formula = preference_formula,

      # arguments that are varied
      move_type = model_df$move_type[model_id],
      ctmc_diffusion_bounds = model_df$ctmc_diffusion_bounds[model_id],
      Movement_ageblk_spec = model_df$ageblock[[model_id]],
      Movement_yearblk_spec = model_df$yearblock[[model_id]]
    )

    # Extract stuff
    data <- input_list$data
    pars <- input_list$par
    map <- input_list$map

    # mess around with starting values
    if(model_df$move_type[model_id] == 1) {
      pars$log_move_diffusion_pars <- rep(log(2), length(pars$log_move_diffusion_pars))
      pars$move_preference_pars <- rnorm(length(pars$move_preference_pars), 0, 0.1)
    }

    # Fit model
    mod <- fit_model(
      data,
      pars,
      map,
      random = NULL,
      newton_loops = 3,
      silent = F,
      do_optim = T
    )

    # Save sdreport
    mod$sd_rep <- sdreport(mod)

    # plot(sim_obj$Movement[3,3,1,-1,1,sim_id], ylim = c(0,1))
    # lines(mod$rep$Movement[3,3,1,-1,1])
    #
    # plot(sim_obj$Movement[3,3,,2,1,sim_id])
    # lines(mod$rep$Movement[3,3,,2,1])
    #
    # image(mod$rep$Movement[2,2,,-1,1])
    # image(sim_obj$Movement[2,2,,-1,1,sim_id])
    #
    # plot(mod$rep$SSB[2,])
    # lines(sim_obj$SSB[2,,sim_id])

    # Return results with metadata
    return(list(
      sim_id = sim_id,
      model_id = model_id,
      move_type = model_df$move_type[model_id],
      ctmc_diffusion_bounds = model_df$ctmc_diffusion_bounds[model_id],
      model = mod,
      pdHess = mod$sd_rep$pdHess,
      objective = mod$optim$objective,
      marg_AIC = marg_AIC(mod$optim)
    ))

  }, error = function(e) {
    # Return error information if model fails
    return(list(
      sim_id = sim_id,
      model_id = model_id,
      move_type = model_df$move_type[model_id],
      ctmc_diffusion_bounds = model_df$ctmc_diffusion_bounds[model_id],
      model = NULL,
      error = as.character(e)
    ))
  })
}
