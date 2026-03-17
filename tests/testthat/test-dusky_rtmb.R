library(SPoRC)
library(testthat)
data("sgl_rg_dusky_data")

test_that("Dusky RTMB model produces expected results", {

  # Setup Model -------------------------------------------------------------

  input_list <- Setup_Mod_Dim(
    # Number of populations
    n_pop = sgl_rg_dusky_data$n_pop,
    natal_region = sgl_rg_dusky_data$natal_region,
    years = sgl_rg_dusky_data$years,
    # vector of years
    ages = sgl_rg_dusky_data$mod_ages,
    # vector of ages
    lens = sgl_rg_dusky_data$lens,
    # number of lengths
    n_regions = sgl_rg_dusky_data$n_regions,
    # number of regions
    n_sexes = sgl_rg_dusky_data$n_sexes,
    # number of sexes
    n_fish_fleets = sgl_rg_dusky_data$n_fish_fleets,
    # number of fishery fleets
    n_srv_fleets = sgl_rg_dusky_data$n_srv_fleets, # number of survey fleets
    n_seas = sgl_rg_dusky_data$n_seas,
    verbose = TRUE # whether to output messages
  )

  # Setup recruitment stuff (using defaults for other stuff)
  input_list <- Setup_Mod_Rec(
    input_list = input_list,

    # Model options
    # Doing bias ramp, but basically setting it so that no lognormal bias correction happens (as in the dusky model)
    do_rec_bias_ramp = 1,
    bias_year = rep(length(sgl_rg_dusky_data$years), 4),
    # do bias ramp (0 == don't do bias ramp, 1 == do bias ramp)
    sigmaR_switch = 1,
    # when to switch from early to late sigmaR (switch in first year)
    ln_sigmaR = array(-0.1068576, dim = c(2, input_list$data$n_pop, input_list$data$n_regions)),
    # Starting values for early and late sigmaR
    rec_model = "mean_rec",
    sigmaR_spec = "fix",
    # fix early sigmaR and late sigmaR
    # recruitment sex ratio
    init_age_strc = 1, # geometric series to derive initial age structure
    ln_global_R0 = log(2.7), # starting value for mean_rec
    t_spawn = sgl_rg_dusky_data$spwn_month
  )

  input_list <- Setup_Mod_Biologicals(
    input_list = input_list,

    # Data inputs
    WAA = sgl_rg_dusky_data$waa_arr,
    MatAA = sgl_rg_dusky_data$mataa_arr,

    # Model options
    # fit lengths
    fit_lengths = 1,
    SizeAgeTrans = sgl_rg_dusky_data$sizeage,
    AgeingError = sgl_rg_dusky_data$age_error_matrix,
    M_spec = "fix",
    # fixing natural mortality
    Fixed_natmort = sgl_rg_dusky_data$fix_natmort,
    addtocomp = 0.00001
  )

  # Setup movement stuff (using defaults for other stuff)
  input_list <- Setup_Mod_Movement(
    input_list = input_list,
    use_fixed_movement = 1,
    Fixed_Movement = NA,
    do_recruits_move = 0
  )

  input_list <- Setup_Mod_Tagging(input_list = input_list, use_conv_fish_tagging = 0)


  input_list <- Setup_Mod_Catch_and_F(
    input_list = input_list,

    # Data inputs
    ObsCatch = sgl_rg_dusky_data$ObsCatch,
    UseCatch = sgl_rg_dusky_data$UseCatch,

    # Model options
    Use_F_pen = 1,
    # whether to use f penalty, == 0 don't use, == 1 use
    sigmaC_spec = "fix",

    # Fixing sigma C and F
    ln_sigmaC = array(log(sqrt(1/2)), dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets)),
    ln_sigmaF = array(log(sqrt(1/2)), dim = c(input_list$data$n_regions, input_list$data$n_seas, input_list$data$n_fish_fleets))
  )

  input_list <- Setup_Mod_FishIdx_and_Comps(
    input_list = input_list,

    # data inputs
    ObsFishIdx = sgl_rg_dusky_data$ObsFishIdx, # fishery index
    ObsFishIdx_SE = sgl_rg_dusky_data$ObsFishIdx_SE, # standard errors
    UseFishIdx = sgl_rg_dusky_data$UseFishIdx, # whether fishery indices are used
    ObsFishAgeComps = sgl_rg_dusky_data$ObsFishAgeComps, # observed fishery ages
    UseFishAgeComps = sgl_rg_dusky_data$UseFishAgeComps, # whether fishery ages are used
    ISS_FishAgeComps = sgl_rg_dusky_data$ISS_FishAgeComps, # input sample size for fishery ages
    ObsFishLenComps = sgl_rg_dusky_data$ObsFishLenComps, # observed fishery lengths
    UseFishLenComps = sgl_rg_dusky_data$UseFishLenComps, # whether fishery lengths are used
    ISS_FishLenComps = sgl_rg_dusky_data$ISS_FishLenComps, # input sample size for fishery lengths

    # Model options
    fish_idx_type = c("none"),
    # indices for fishery
    FishAgeComps_LikeType = c("Multinomial"),
    # age comp likelihoods for fishery fleet
    FishLenComps_LikeType = c("Multinomial"),
    # length comp likelihoods for fishery
    FishAgeComps_Type = c("agg_Year_1-terminal_Fleet_1"),
    # age comp structure for fishery
    FishLenComps_Type = c("agg_Year_1-terminal_Fleet_1")
    # length comp structure for fishery
  )

  # Setup survey indices and compositions
  input_list <- Setup_Mod_SrvIdx_and_Comps(
    input_list = input_list,

    # data inputs
    ObsSrvIdx = sgl_rg_dusky_data$ObsSrvIdx, # observed survey index
    ObsSrvIdx_SE = sgl_rg_dusky_data$ObsSrvIdx_SE / sgl_rg_dusky_data$ObsSrvIdx, # lognormal SD
    UseSrvIdx = sgl_rg_dusky_data$UseSrvIdx, # whether survey indices are used
    ObsSrvAgeComps = sgl_rg_dusky_data$ObsSrvAgeComps, # observed survey ages
    ISS_SrvAgeComps = sgl_rg_dusky_data$ISS_SrvAgeComps, # input sample size for survey ages
    UseSrvAgeComps = sgl_rg_dusky_data$UseSrvAgeComps, # whether survey ages are used
    ObsSrvLenComps = sgl_rg_dusky_data$ObsSrvLenComps, # observed survey lengths
    UseSrvLenComps = sgl_rg_dusky_data$UseSrvLenComps, # whether survey lengths are used
    ISS_SrvLenComps = sgl_rg_dusky_data$ISS_SrvLenComps, # input sample size for survey lengths

    # Model options
    srv_idx_type = c("biom"),
    # abundance and biomass for survey fleet 1
    SrvAgeComps_LikeType = c("Multinomial"),
    # survey age composition likelihood for survey fleet 1
    SrvLenComps_LikeType = c("Multinomial"),
    #  survey length composition likelihood for survey fleet 1
    SrvAgeComps_Type = c(
      "agg_Year_1-terminal_Fleet_1"
    ),
    # survey age comp type

    SrvLenComps_Type = c(
      "agg_Year_1-terminal_Fleet_1"
    )
    # survey length comp type
  )


  input_list <- Setup_Mod_Fishsel_and_Q(

    input_list = input_list,

    # Model options
    # fishery selectivity, whether continuous time-varying
    cont_tv_fish_sel = c("none_Fleet_1"),
    # fishery selectivity blocks
    fish_sel_blocks = c("none_Fleet_1"),
    # fishery selectivity form
    fish_sel_model = c("logist2_Fleet_1"),
    # fishery catchability blocks
    fish_q_blocks = c("none_Fleet_1"),
    # whether to estiamte all fixed effects for fishery selectivity
    fish_fixed_sel_pars_spec = c("est_all"),
    # whether to estiamte all fixed effects for fishery catchability
    fish_q_spec = c("fix")
  )

  # Setup survey selectivity and catchability
  # Set up prior for survey catchability
  srv_q_prior <- data.frame(
    region = 1,
    block = 1,
    fleet = 1,
    mu = 1,
    sd = 0.447213595
  )

  input_list <- Setup_Mod_Srvsel_and_Q(
    input_list = input_list,

    # Model options
    # survey selectivity, whether continuous time-varying
    cont_tv_srv_sel = c("none_Fleet_1"),
    # survey selectivity blocks
    srv_sel_blocks = c("none_Fleet_1"),
    # survey selectivity form
    srv_sel_model = c("logist2_Fleet_1"),
    # survey catchability blocks
    srv_q_blocks = c("none_Fleet_1"),
    # whether to estiamte all fixed effects for survey selectivity
    srv_fixed_sel_pars_spec = c("est_all"),
    # whether to estiamte all fixed effects for survey catchability
    srv_q_spec = c("est_all"),
    Use_srv_q_prior = 1,
    # Use catchability prior
    srv_q_prior = srv_q_prior,
    # survey timing
    t_srv = array(0, dim = c(input_list$data$n_regions,
                             input_list$data$n_seas,
                             input_list$data$n_srv_fleets))
  )

  # catch weigthing for duskies
  Wt_Catch <- array(0, dim = c(sgl_rg_dusky_data$n_regions, length(sgl_rg_dusky_data$years), sgl_rg_dusky_data$n_seas, sgl_rg_dusky_data$n_fish_fleets))
  Wt_Catch[,which(sgl_rg_dusky_data$years %in% 1977:1991),1,] <- 2
  Wt_Catch[,-which(sgl_rg_dusky_data$years %in% 1977:1991),1,] <- 50

  input_list <- Setup_Mod_Weighting(
    input_list = input_list,
    Wt_Catch = Wt_Catch,
    Wt_FishIdx = 1,
    Wt_SrvIdx = 1.66,
    Wt_Rec = 1,
    Wt_F = 2,
    Wt_Tagging = 0,
    Wt_FishAgeComps = array(1, dim = c(input_list$data$n_regions,
                                       length(input_list$data$years),
                                       input_list$data$n_seas,
                                       input_list$data$n_sexes,
                                       input_list$data$n_fish_fleets)),
    Wt_FishLenComps = array(1, dim = c(input_list$data$n_regions,
                                       length(input_list$data$years),
                                       input_list$data$n_seas,
                                       input_list$data$n_sexes,
                                       input_list$data$n_fish_fleets)),
    Wt_SrvAgeComps = array(1, dim = c(input_list$data$n_regions,
                                      length(input_list$data$years),
                                      input_list$data$n_seas,
                                      input_list$data$n_sexes,
                                      input_list$data$n_srv_fleets)),
    Wt_SrvLenComps = array(0, dim = c(input_list$data$n_regions,
                                      length(input_list$data$years),
                                      input_list$data$n_seas,
                                      input_list$data$n_sexes,
                                      input_list$data$n_srv_fleets))
  )

  data <- input_list$data
  parameters <- input_list$par
  mapping <- input_list$map

  # Fit model
  dusky_rtmb_model <- fit_model(data,
                                parameters,
                                mapping,
                                random = NULL,
                                newton_loops = 3,
                                silent = TRUE
  )

  dusky_rtmb_model$sdrep <- RTMB::sdreport(dusky_rtmb_model) # get standard error report

  ssb_expected_vec <- c(
    12797.54, 12250.75, 11909.20, 11599.08, 11192.37,
    10790.43, 10507.33, 10338.19, 10565.86, 11306.05,
    12302.58, 13455.75, 14217.22, 14785.66, 15393.80,
    15934.92, 15768.58, 16105.93, 17027.92, 18489.64,
    20460.80, 22348.89, 23729.31, 24165.39, 24860.83,
    25933.28, 26978.12, 28242.44, 29758.03, 31503.33,
    33118.19, 34220.20, 35018.81, 35787.45, 36175.29,
    36467.88, 35801.61, 35360.99, 35004.50, 34987.86,
    35038.66, 35694.16, 36410.03, 37326.83, 38208.01,
    38414.28, 38341.49, 37409.68
  )

  rec_expected_vec <- c(
    1.719130, 1.821452, 2.325670, 5.221616,
    6.055137, 6.385826, 3.798459, 4.855644,
    3.303878, 3.093645, 2.270344, 8.978767,
    5.840166, 17.847508, 12.477896, 10.773610,
    3.022077, 7.685572, 5.743625, 16.754576,
    3.156288, 9.194834, 18.978940, 2.694959,
    12.000002, 14.767476, 6.897940, 9.725435,
    9.246636, 4.339117, 4.744178, 5.182619,
    6.320576, 7.264755, 11.971134, 9.814706,
    7.959109, 17.536878, 4.838579, 7.376432,
    6.538344, 6.301670, 2.189419, 2.404738,
    2.061323, 2.393618, 2.529818, 2.838106
  )

  expect_equal(dusky_rtmb_model$rep$SSB[1,1,], ssb_expected_vec, tolerance = 1e-2)
  expect_equal(dusky_rtmb_model$rep$Rec[1,1,], rec_expected_vec, tolerance = 1e-2)
  expect_true(dusky_rtmb_model$sdrep$pdHess)

})
