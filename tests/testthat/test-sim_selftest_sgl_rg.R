library(SPoRC)
library(testthat)
library(reshape2)
library(dplyr)

test_that("Simulation self-test produces approximately unbiased SSB results", {

  # Setup Operating Model ---------------------------------------------------
  ### Setup Model Dimensions --------------------------------------------------
  sim_list <- Setup_Sim_Dim(n_sims = 50, # number of simulations
                            n_yrs = 30, # number of years
                            n_regions = 1,  # number of regions
                            n_ages = 10, # number of ages
                            n_lens = NULL, # number of lengths
                            n_sexes = 1, # number of sexes
                            n_fish_fleets = 1, # number of fishery fleets
                            n_srv_fleets = 1, # number of survey fleets
                            n_pop = 1 # number of populations
  )

  ### Setup Simulation Containers ---------------------------------------------
  sim_list <- Setup_Sim_Containers(sim_list)

  ### Setup Fishing Processes -------------------------------------------------
  sim_list <- Setup_Sim_Fishing(sim_list = sim_list, # update simulate list
                                # Logistic selectivity
                                fish_sel_input = replicate(
                                  n = sim_list$n_sims,
                                  array(rep(1 / (1 + exp(-3 * ((1:sim_list$n_ages) - 5))), each = sim_list$n_yrs),
                                        dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages,
                                                sim_list$n_sexes, sim_list$n_fish_fleets))
                                )
  )

  ### Setup Survey Process ----------------------------------------------------
  sim_list <- Setup_Sim_Survey(
    sim_list = sim_list,
    # Logistic selectivity
    srv_sel_input = replicate(
      n = sim_list$n_sims,
      array(rep(1 / (1 + exp(-1 * ((1:sim_list$n_ages) - 3))), each = sim_list$n_yrs),
            dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages,
                    sim_list$n_sexes, sim_list$n_srv))
    )
  )

  ### Setup Biological Dynamics -----------------------------------------------
  sim_list <- suppressWarnings(
    Setup_Sim_Biologicals(
      sim_list = sim_list, # simualtion list
      natmort_input = replicate(n = sim_list$n_sims, array(0.3, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs,
                                                                        sim_list$n_ages, sim_list$n_sexes))), # natural mortality
      WAA_input = replicate(n = sim_list$n_sims, array(rep(5 / (1 + exp(-3 * ((1:sim_list$n_ages) - 3))), each = sim_list$n_yrs),
                                                       dim = c(sim_list$n_pop,sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes))), # weight at age
      WAA_fish_input = replicate(n = sim_list$n_sims, array(rep(5 / (1 + exp(-3 * ((1:sim_list$n_ages) - 3))), each = sim_list$n_yrs),
                                                            dim = c(sim_list$n_pop,  sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes, sim_list$n_fish_fleets))), # fishery weight at age
      WAA_srv_input = replicate(n = sim_list$n_sims, array(rep(5 / (1 + exp(-3 * ((1:sim_list$n_ages) - 3))), each = sim_list$n_yrs),
                                                           dim = c(sim_list$n_pop,  sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes, sim_list$n_srv_fleets))), # survey weight at age
      MatAA_input = replicate(n = sim_list$n_sims, array(rep(1 / (1 + exp(-3 * ((1:sim_list$n_ages) - 3))), each = sim_list$n_yrs),
                                                         dim = c(sim_list$n_pop,  sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes))) # maturity at age
    )
  )

  ### Setup Tagging and Movement -----------------------------------------------------------
  sim_list <- Setup_Sim_Tagging(
    sim_list = sim_list, # simulation list
    use_conv_fish_tagging = 0
  )

  # No Movement
  sim_list$Movement <- array(1, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes, sim_list$n_sims))

  ### Setup Recruitment Processes ---------------------------------------------
  sim_list <- Setup_Sim_Rec(
    sim_list = sim_list,
    R0_input = replicate(n = sim_list$n_sims, expr = array(5, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs))), # R0
    rinit_input = array(2, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_sims)), # rinit
    use_rinit = 1,
    ln_sigmaR = array(log(1), dim = c(2, sim_list$n_pop, sim_list$n_region)),
    recruitment_opt = 'mean_rec',
    init_age_strc = 1
  )

  ## Simulate Data -----------------------------------------------------------
  set.seed(777)
  sim_obj <- Simulate_Pop_Static(sim_list = sim_list, output_path = NULL) # get simulated datasets

  # Define Estimation Model -------------------------------------------------
  setup_em <- function(sim_obj, sim) {

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
      n_pop = sim_obj$n_pop,
      natal_region = sim_obj$natal_region,
      verbose = F
    )

    # Recruitment setup
    input_list <- Setup_Mod_Rec(
      input_list = input_list,
      do_rec_bias_ramp = 0, # not doing bias ramp
      sigmaR_switch = 1, # when to switch from early to late sigmaR (switch in first year)
      ln_sigmaR = array(log(1), c(2, input_list$data$n_pop, input_list$data$n_regions)), # 2 values for early and late sigma
      rec_model = "mean_rec",
      use_rinit = 1,
      sigmaR_spec = "fix", # fix early sigmaR and late sigmaR
      init_age_strc = 1, # scalar geometric series to derive initial age structure
      equil_init_age_strc = 2, # estimating all intial age deviations
      ln_global_R0 = log(5),
      ln_rinit = log(2)
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
      Fixed_natmort = array(0.3, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years),
                                         length(input_list$data$ages), input_list$data$n_sexes))
    )

    # Movement and tagging
    input_list <- Setup_Mod_Tagging(input_list = input_list, use_conv_fish_tagging = 0)
    input_list <- Setup_Mod_Movement(
      input_list = input_list,
      use_fixed_movement = 1,
      Fixed_Movement = NA,
      do_recruits_move = 0
    )

    # Fishery catch & fishing mortality
    suppressWarnings(
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
        ln_sigmaF = array(log(1), dim = c(input_list$data$n_regions, input_list$data$n_seas, input_list$data$n_fish_fleets))
      )
    )

    # Survey selectivity and catchability
    input_list <- Setup_Mod_FishIdx_and_Comps(
      input_list = input_list,
      # Data inputs
      ObsFishIdx = sim_data$ObsFishIdx,
      ObsFishIdx_SE = sim_data$ObsFishIdx_SE,
      UseFishIdx = sim_data$UseFishIdx,
      ObsFishAgeComps = sim_data$ObsFishAgeComps,
      ObsFishLenComps = sim_data$ObsFishLenComps,
      UseFishAgeComps = sim_data$UseFishAgeComps,
      UseFishLenComps = sim_data$UseFishLenComps,
      ISS_FishAgeComps = sim_data$ISS_FishAgeComps,
      ISS_FishLenComps = sim_data$ISS_FishLenComps,
      # Model options
      fish_idx_type = c("biom"),
      FishAgeComps_LikeType = c("Multinomial"),
      FishLenComps_LikeType = c("none"),
      FishAgeComps_Type = c("agg_Year_1-terminal_Fleet_1"),
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
      SrvAgeComps_Type = c("agg_Year_1-terminal_Fleet_1"),
      SrvLenComps_Type = c("none_Year_1-terminal_Fleet_1")
    )


    # Fishery selectivity and catchability
    input_list <- Setup_Mod_Fishsel_and_Q(
      input_list = input_list,
      # Model options
      fish_sel_model = c("logist2_Fleet_1"), # fishery selex model
      fish_fixed_sel_pars_spec = c("est_all"), # whether to estiamte all fixed effects for fishery selectivity
      fish_q_spec = "est_all" # estimate fishery q
    )

    # Survey selectivity and catchability
    input_list <- Setup_Mod_Srvsel_and_Q(
      input_list = input_list,
      # Model options
      srv_sel_model = c("logist2_Fleet_1"), # survey selectivity form
      srv_fixed_sel_pars_spec = c("est_all"), # whether to estimate all fixed effects for survey selectivity
      srv_q_spec = c("est_all")  # whether to estiamte all fixed effects for survey catchability
    )

    # Data weighting
    input_list <- Setup_Mod_Weighting(
      input_list = input_list,
      Wt_Catch = 1,
      Wt_FishIdx = 1,
      Wt_SrvIdx = 1,
      Wt_Rec = 1,
      Wt_F = 1,
      Wt_Tagging = 0,
      Wt_FishAgeComps = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas,
                                         input_list$data$n_sexes, input_list$data$n_fish_fleets)),
      Wt_FishLenComps = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas,
                                         input_list$data$n_sexes, input_list$data$n_fish_fleets)),
      Wt_SrvAgeComps = array(1, dim = c(input_list$data$n_regions,length(input_list$data$years), input_list$data$n_seas,
                                        input_list$data$n_sexes, input_list$data$n_srv_fleets)),
      Wt_SrvLenComps = array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas,
                                        input_list$data$n_sexes, input_list$data$n_srv_fleets))
    )

    return(input_list)

  }

  ssb_results <- array(NA, dim = c(sim_list$n_yrs, sim_list$n_sims)) # storage container
  rinit_results <- vector()
  r0_results <- vector()

  for(i in 1:sim_obj$n_sims) {

    input_list <- setup_em(sim_obj, sim = i) # setup EM

    # fit model
    model <- fit_model(input_list$data,
                       input_list$par,
                       input_list$map,
                       random = NULL,
                       silent = T
    )

    # the EM configuration is identical across simulations, so one check covers it
    if(i == 1) expect_jnLL_decomposes(model)

    ssb_results[,i] <- as.vector(model$rep$SSB) # save results
    rinit_results[i] <- model$rep$rinit
    r0_results[i] <- model$rep$R0

  } # end i loop

  (ssb_results[,1] - sim_obj$SSB[,,,1]) / sim_obj$SSB[,,,1]

  # Process SSB results
  ssb_df_res <- reshape2::melt(ssb_results) %>%
    rename(Year = Var1, Sim = Var2, Est = value) %>%
    dplyr::left_join(reshape2::melt(sim_obj$SSB) %>%
                       dplyr::rename(Pop = Var1, Region = Var2, Year = Var3, Sim = Var4, True = value),
                     by = c("Year", "Sim")) %>%
    dplyr::mutate(RE = (Est - True) / True)

  # check to see if relative error is within 2%
  expect_equal(median(ssb_df_res$RE), 0, tolerance = 0.02)
  expect_equal(median(((rinit_results - 2) / 2) ), 0, tolerance = 0.05)
  expect_equal(median(((r0_results - 5) / 5) ), 0, tolerance = 0.05)

})

