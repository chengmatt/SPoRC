# Regression test. Expected SSB and recruitment come from a previously validated SPoRC fit, not from
# hand-derived values, so a mismatch means something moved a fitted result. See tests/README.md.

library(SPoRC)
library(testthat)
data("sgl_rg_ebswp_data")

test_that("Single-region EBS Pollock RTMB model produces expected results", {

  n_srv <- sgl_rg_ebswp_data$n_srv_fleets

  ## Initialize model dimensions and data list----
  input_list <- Setup_Mod_Dim(
    years = sgl_rg_ebswp_data$years,
    # vector of years
    ages = sgl_rg_ebswp_data$ages,
    # vector of ages
    lens = NA,
    # number of lengths
    n_regions = 1,
    # number of regions
    n_sexes = 1,
    # number of sexes
    n_fish_fleets = 1,
    # number of fishery fleets
    n_srv_fleets = n_srv, # number of survey fleets
    # number of seasons
    n_seas = sgl_rg_ebswp_data$n_seas,
    # Populaiton stuff
    n_pop = sgl_rg_ebswp_data$n_pop,
    natal_region = sgl_rg_ebswp_data$natal_region,
    verbose = FALSE
  )

  inv_steepness <- function(s) qlogis((s - 0.2) / 0.8)

  # Setup recruitment stuff (using defaults for other stuff)
  input_list <- Setup_Mod_Rec(
    input_list = input_list,

    # Model options
    do_rec_bias_ramp = 0,
    # do bias ramp (0 == don't do bias ramp, 1 == do bias ramp)
    sigmaR_switch = 1,
    # when to switch from early to late sigmaR (switch in first year)
    ln_sigmaR = array(log(1), dim = c(2, input_list$data$n_pop, input_list$data$n_regions)),
    # Starting values for early and late sigmaR
    rec_model = "bh_rec",
    # recruitment model
    steepness_h = array(inv_steepness(0.623013), dim = c(input_list$data$n_pop, input_list$data$n_regions)),
    h_spec = "fix",
    # fixing steepness
    sigmaR_spec = "fix",
    # fix early sigmaR and late sigmaR
    init_age_strc = 1,
    ln_global_R0 = 10,
    t_spawn = 0.25,
    equil_init_age_strc = 2
    # starting value for r0
  )

  # Setup a fixed natural mortality array for use
  fix_natmort <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years), length(input_list$data$ages), 1))
  fix_natmort[,,,1,] <- 0.9 # age 1 M
  fix_natmort[,,,2,] <- 0.45 # age 2 M
  fix_natmort[,,,-c(1,2),] <- 0.3 # age 3+ M

  suppressWarnings(
    input_list <- Setup_Mod_Biologicals(
      input_list = input_list,

      # Data inputs
      WAA = sgl_rg_ebswp_data$WAA,
      # the assessment has a separate weight at age matrix for the fishery
      # and for each survey index
      WAA_fish = sgl_rg_ebswp_data$WAA_fish,
      WAA_srv = sgl_rg_ebswp_data$WAA_srv,
      MatAA = sgl_rg_ebswp_data$MatAA,

      # Model options
      # mean and sd for M prior
      fit_lengths = 0,
      # don't fit length compositions
      M_spec = "fix",
      # fixing natural mortality
      Fixed_natmort = fix_natmort
    )
  )

  # Setup movement stuff (using defaults for other stuff)
  input_list <- Setup_Mod_Movement(
    input_list = input_list,
    use_fixed_movement = 1,
    Fixed_Movement = NA,
    do_recruits_move = 0
  )

  suppressWarnings(
    input_list <- Setup_Mod_Catch_and_F(
      input_list = input_list,

      # Data inputs
      ObsCatch = sgl_rg_ebswp_data$ObsCatch,
      UseCatch = sgl_rg_ebswp_data$UseCatch,

      # Model options
      Use_F_pen = 1,
      # whether to use f penalty, == 0 don't use, == 1 use
      sigmaC_spec = "fix",
      # fixing catch standard deviation
      ln_sigmaC = array(log(0.05), dim = c(1, length(input_list$data$years), input_list$data$n_seas, 1))
      # starting / fixed value for catch standard deviation
    )
  )

  input_list <- Setup_Mod_FishIdx_and_Comps(
    input_list = input_list,
    # data inputs
    ObsFishIdx = sgl_rg_ebswp_data$ObsFishIdx,
    ObsFishIdx_SE = sgl_rg_ebswp_data$ObsFishIdx_SE,
    UseFishIdx = sgl_rg_ebswp_data$UseFishIdx,
    ObsFishAgeComps = sgl_rg_ebswp_data$ObsFishAgeComps,
    UseFishAgeComps = sgl_rg_ebswp_data$UseFishAgeComps,
    ISS_FishAgeComps = sgl_rg_ebswp_data$ISS_FishAgeComps,
    ObsFishLenComps = array(NA_real_, dim = c(1, length(input_list$data$years), input_list$data$n_seas, length(input_list$data$lens), 1, 1)),
    UseFishLenComps = array(0, dim = c(1, length(input_list$data$years), input_list$data$n_seas, 1)),
    ISS_FishLenComps = array(0, dim = c(1, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes, 1)),

    # Model options
    fish_idx_type = c("biom"),
    # indices for fishery
    FishAgeComps_LikeType = c("Multinomial"),
    # age comp likelihoods for fishery fleet
    FishLenComps_LikeType = c("none"),
    # length comp likelihoods for fishery
    FishAgeComps_Type = c("agg_Year_1-terminal_Fleet_1"),
    # age comp structure for fishery
    FishLenComps_Type = c("none_Year_1-terminal_Fleet_1")
    # length comp structure for fishery
  )

  # Setup survey indices and compositions
  input_list <- Setup_Mod_SrvIdx_and_Comps(
    input_list = input_list,

    # data inputs
    ObsSrvIdx = sgl_rg_ebswp_data$ObsSrvIdx,
    ObsSrvIdx_SE = sgl_rg_ebswp_data$ObsSrvIdx_SE,
    UseSrvIdx = sgl_rg_ebswp_data$UseSrvIdx,
    ObsSrvAgeComps = sgl_rg_ebswp_data$ObsSrvAgeComps,
    ISS_SrvAgeComps = sgl_rg_ebswp_data$ISS_SrvAgeComps,
    UseSrvAgeComps = sgl_rg_ebswp_data$UseSrvAgeComps,
    ObsSrvLenComps = array(NA_real_, dim = c(1, length(input_list$data$years), input_list$data$n_seas, length(input_list$data$lens), 1, n_srv)),
    UseSrvLenComps = array(0, dim = c(1, length(input_list$data$years), input_list$data$n_seas, n_srv)),
    ISS_SrvLenComps = array(0, dim = c(1, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes, n_srv)),

    # Model options
    srv_idx_type = c("biom", "biom", "biom", "abd"),
    # biomass for survey fleets 1 to 3, abundance for survey fleet 4
    srv_idx_ages = list(NULL, NULL, NULL, 1),
    # fleet 4 is the acoustic survey's age 1 abundance, so it sees age 1 only
    SrvIdx_LikeType = c("mvn", "lognormal", "normal", "lognormal"),
    # index likelihood for survey fleet 1, 2, 3, and 4
    SrvIdx_Cov = list(sgl_rg_ebswp_data$SrvIdx_Cov, NULL, NULL, NULL),
    # the bottom trawl index is fit with a full covariance matrix
    SrvAgeComps_LikeType = c("Multinomial", "Multinomial", "none", "none"),
    # survey age composition likelihood for survey fleet 1, 2, 3, and 4
    SrvAgeComps_bins = list(NULL, 2:15, NULL, NULL),
    # the acoustic compositions are normalized over ages 2-15 only
    SrvLenComps_LikeType = rep("none", n_srv),
    #  survey length composition likelihood for survey fleet 1, 2, 3, and 4
    SrvAgeComps_Type = c(
      "agg_Year_1-terminal_Fleet_1",
      "agg_Year_1-terminal_Fleet_2",
      "none_Year_1-terminal_Fleet_3",
      "none_Year_1-terminal_Fleet_4"
    ),
    # survey age comp type

    SrvLenComps_Type = paste0("none_Year_1-terminal_Fleet_", 1:n_srv),
    # survey length comp type
    t_srv = array(c(0.5, 0.5, 0, 0.5), dim = c(1, 1, n_srv))
    # fraction of the year elapsed when each survey occurs
  )


  # Setup fishery selectivity and catchability
  input_list <- Setup_Mod_Fishsel_and_Q(

    input_list = input_list,

    # Model options
    # fishery selectivity, whether continuous time-varying
    cont_tv_fish_sel = c("2dar1_Fleet_1"),
    fishsel_pe_pars_spec = "fix", # doing penalized likelihood for selex devs
    fish_sel_devs_spec = "est_all", # estimating all sel devs
    corr_opt_semipar = "corr_zero_y_b", # making sure 2d correaltions are 0, collapses to a simple iid case
    # fishery selectivity blocks
    fish_sel_blocks = c("none_Fleet_1"),
    # fishery selectivity form
    fish_sel_model = c("logist1_Fleet_1"),
    # fishery catchability blocks
    fish_q_blocks = c("none_Fleet_1"),
    # whether to estiamte all fixed effects for fishery selectivity
    fish_fixed_sel_pars_spec = c("est_all"),
    # whether to estiamte all fixed effects for fishery catchability
    fish_q_spec = c("est_all")
  )

  # Setup survey selectivity and catchability
  input_list <- Setup_Mod_Srvsel_and_Q(
    input_list = input_list,

    # Model options
    # survey selectivity, whether continuous time-varying
    cont_tv_srv_sel = c("iid_Fleet_1", "2dar1_Fleet_2", "2dar1_Fleet_3", "none_Fleet_4"),
    srvsel_pe_pars_spec = rep("fix", n_srv), # penalize survey selex devs
    srv_sel_devs_spec = c("est_all", "est_all", "est_shared_f_2", "fix"), # estimating all srv selex devs
    corr_opt_semipar = c(NA, "corr_zero_y_b", "corr_zero_y_b", NA), # setting corelations at 0, so 2dar1 collapses to simple iid semi-parametric devs

    # survey selectivity blocks
    srv_sel_blocks = paste0("none_Fleet_", 1:n_srv),
    # survey selectivity form
    srv_sel_model = paste0("logist1_Fleet_", 1:n_srv),
    # survey catchability blocks
    srv_q_blocks = paste0("none_Fleet_", 1:n_srv),
    # the vessel of opportunity index shares the acoustic survey's selectivity, and fleet 4 sees
    # age 1 only, where selectivity is absorbed into catchability and so is fixed
    srv_fixed_sel_pars_spec = c("est_all", "est_all", "est_shared_f_2", "fix"),
    # whether to estiamte all fixed effects for survey catchability
    srv_q_spec = rep("est_all", n_srv),
    t_srv = array(c(0.5, 0.5, 0, 0.5), dim = c(1, 1, n_srv))
  )

  # Setup tagging stuff
  input_list <- Setup_Mod_Tagging(input_list = input_list, use_conv_fish_tagging = 0)

  input_list <- Setup_Mod_Weighting(
    input_list = input_list,
    Wt_Catch = 1,
    Wt_FishIdx = 1,
    Wt_SrvIdx = 1,
    Wt_Rec = 1,
    Wt_F = 1,
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
    Wt_SrvLenComps = array(1, dim = c(input_list$data$n_regions,
                                      length(input_list$data$years),
                                      input_list$data$n_seas,
                                      input_list$data$n_sexes,
                                      input_list$data$n_srv_fleets))
  )


  # extract out lists updated with helper functions
  data <- input_list$data
  parameters <- input_list$par
  mapping <- input_list$map

  # selex sigma to fix at, given penalized likelihood
  parameters$fishsel_pe_pars[,4,,] <- log(0.075) # fishery selex variance
  parameters$srvsel_pe_pars[,1:2,,1] <- log(0.075) # survey BTS - a50 and delta variance
  parameters$srvsel_pe_pars[,4,,2] <- log(0.15) # survey ATS and ato variance


  # Fit model
  ebswp_rtmb_model <- fit_model(data,
                                parameters,
                                mapping,
                                random = NULL,
                                newton_loops = 3,
                                silent = TRUE
  )

  ebswp_rtmb_model$sdrep <- RTMB::sdreport(ebswp_rtmb_model)

  ssb_expected_vec <- c(
    694.8326, 699.3691, 704.0999, 745.4423, 803.8686,
    962.8408, 1131.7134, 1282.9910, 1190.9568, 934.7983,
    533.1818, 506.8241, 615.7872, 799.4256, 698.6487,
    701.6111, 958.1317, 1594.5831, 2291.7074, 3143.5176,
    3259.7215, 3699.5436, 3583.1245, 3659.2346, 3535.4308,
    3122.3689, 2683.8432, 2176.3735, 2097.6557, 2808.7953,
    3273.5585, 3413.6688, 3247.6993, 3059.7111, 2578.7963,
    2740.1746, 2708.9987, 2815.2148, 2657.8062, 2716.5948,
    3041.1099, 2714.9118, 2532.4614, 2133.4825, 1695.9670,
    1860.6397, 1908.5460, 2255.6411, 2597.5750, 2857.4472,
    2777.7454, 2633.4563, 2923.9240, 3345.0265, 3107.4955,
    2912.7880, 2138.6085, 2389.1280, 3534.8520, 3703.9488,
    4077.2457
  )

  rec_expected_vec <- c(
    3812.4279, 15060.2989, 10703.1408, 24403.7352, 22478.2716,
    24800.2775, 19296.6973, 10564.1322, 9839.5068, 28857.0392,
    16665.3268, 13833.9796, 13383.3619, 14306.0665, 24266.0858,
    56356.3490, 27787.6534, 31848.6662, 16023.2883, 46858.9081,
    13793.0410, 34400.7138, 14534.7859, 7385.0682, 5554.9825,
    11524.9723, 52018.8598, 27263.0309, 20981.5547, 39633.1342,
    13441.3224, 9896.8885, 23446.5729, 32115.7603, 15444.8400,
    16581.6138, 25140.5460, 35326.6209, 24129.8781, 15674.5988,
    7305.0766, 5083.3264, 13334.4987, 28543.3140, 12636.8873,
    46936.8527, 21532.4394, 14016.3195, 12484.1604, 47318.1110,
    48218.6427, 17951.9588, 6079.9893, 6087.6974, 13348.3361,
    90544.7358, 26971.2456, 19199.0511, 17507.4517, 30594.0324,
    39440.3696
  )

  expect_equal(ebswp_rtmb_model$rep$SSB[1,1,], ssb_expected_vec, tolerance = 1e-2)
  expect_equal(ebswp_rtmb_model$rep$Rec[1,1,], rec_expected_vec, tolerance = 1e-2)
  expect_true(ebswp_rtmb_model$sdrep$pdHess)
  expect_jnLL_decomposes(ebswp_rtmb_model)


})
