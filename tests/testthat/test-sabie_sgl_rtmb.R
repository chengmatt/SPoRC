library(SPoRC)
library(testthat)
data("sgl_rg_sable_data")

test_that("Single-region Sablefish RTMB model produces expected results", {

  # Initialize model dimensions and data list
  input_list <- Setup_Mod_Dim(years = sgl_rg_sable_data$years, # vector of years
                              ages = sgl_rg_sable_data$ages, # vector of ages
                              lens = sgl_rg_sable_data$lens, # number of lengths
                              n_regions = sgl_rg_sable_data$n_regions, # number of regions
                              n_sexes = sgl_rg_sable_data$n_sexes, # number of sexes == 1, female, == 2 male
                              n_fish_fleets = sgl_rg_sable_data$n_fish_fleets, # number of fishery fleet == 1, fixed gear, == 2 trawl gear
                              n_srv_fleets = sgl_rg_sable_data$n_srv_fleets, # number of survey fleets
                              n_seas = sgl_rg_sable_data$n_seas, # number of seasons
                              verbose = FALSE
  )

  # Setup recruitment stuff (using defaults for other stuff)
  input_list <- Setup_Mod_Rec(input_list = input_list, # input data list from above
                              # Model options
                              do_rec_bias_ramp = 1, # do bias ramp (0 == don't do bias ramp, 1 == do bias ramp)
                              # breakpoints for bias ramp (1 == no bias ramp - 1960 - 1980, 2 == ascending limb of bias ramp - 1980 - 1990,
                              # 3 == full bias correction - 1990 - 2022, == 4 no bias correction - terminal year of recruitment estimate)
                              bias_year = c(length(1960:1979), length(1960:1989), (length(1960:2023) - 5), length(1960:2024) - 2) + 1,
                              sigmaR_switch = as.integer(length(1960:1975)), # when to switch from early to late sigmaR
                              dont_est_recdev_last = 1, # don't estimate last recruitment deviate
                              ln_sigmaR = array(log(c(0.4, 1.2)), dim = c(2, input_list$data$n_pop, input_list$data$n_regions)),
                              rec_model = "mean_rec", # recruitment model
                              sigmaR_spec = "fix_early_est_late", # fix early sigmaR, estiamte late sigmaR
                              InitDevs_spec = NULL, # estimate all initial deviations
                              RecDevs_spec = NULL, # stiamte all recruitment deivations
                              init_age_strc = 1,
                              init_F_prop = 0.1
  )

  # Specificying natural mortality fixed array
  fixed_natmort <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years), length(input_list$data$ages), input_list$data$n_sexes))
  fixed_natmort[,,,,1] <- 0.1134156 # fix female M
  fixed_natmort[,,,,2] <- 0.1052175 # fix male M

  input_list <- Setup_Mod_Biologicals(input_list = input_list,
                                      # Data inputs
                                      WAA = sgl_rg_sable_data$WAA, # weight-at-age
                                      MatAA = sgl_rg_sable_data$MatAA, # maturity at age
                                      AgeingError = as.matrix(sgl_rg_sable_data$age_error), # ageing error
                                      SizeAgeTrans = sgl_rg_sable_data$SizeAgeTrans, # size age transition matrix
                                      # Model options
                                      Use_M_prior = 0, # use natural mortality prior
                                      fit_lengths = 1, # fitting length compositions
                                      M_spec = "fix",
                                      Fixed_natmort = fixed_natmort)

  # Setup movement stuff (using defaults for other stuff)
  input_list <- Setup_Mod_Movement(input_list = input_list,
                                   use_fixed_movement = 1,
                                   Fixed_Movement = NA,
                                   do_recruits_move = 0
  )

  # Setup catch and fishing mortality stuff
  suppressWarnings(
    input_list <- Setup_Mod_Catch_and_F(input_list = input_list,
                                        # Data inputs
                                        ObsCatch = sgl_rg_sable_data$ObsCatch,
                                        UseCatch = sgl_rg_sable_data$UseCatch,
                                        # Model options
                                        Use_F_pen = 1, # whether to use f penalty, == 0 don't use, == 1 use
                                        sigmaC_spec = 'fix'
    )
  )

  # Setup fishery indices and compositions
  input_list <- Setup_Mod_FishIdx_and_Comps(input_list = input_list,
                                            # data inputs
                                            ObsFishIdx = sgl_rg_sable_data$ObsFishIdx,
                                            ObsFishIdx_SE = sgl_rg_sable_data$ObsFishIdx_SE,
                                            UseFishIdx =  sgl_rg_sable_data$UseFishIdx,
                                            ObsFishAgeComps = sgl_rg_sable_data$ObsFishAgeComps,
                                            UseFishAgeComps = sgl_rg_sable_data$UseFishAgeComps,
                                            ISS_FishAgeComps = sgl_rg_sable_data$ISS_FishAgeComps,
                                            ObsFishLenComps = sgl_rg_sable_data$ObsFishLenComps,
                                            UseFishLenComps = sgl_rg_sable_data$UseFishLenComps,
                                            ISS_FishLenComps = sgl_rg_sable_data$ISS_FishLenComps,

                                            # Model options
                                            fish_idx_type = c("biom", "none"), # biomass indices for fishery fleet 1 and 2
                                            FishAgeComps_LikeType = c("Multinomial", "none"), # age comp likelihoods for fishery fleet 1 and 2
                                            FishLenComps_LikeType = c("Multinomial", "Multinomial"), # length comp likelihoods for fishery fleet 1 and 2
                                            FishAgeComps_Type =  c("agg_Year_1-terminal_Fleet_1",
                                                                   "none_Year_1-terminal_Fleet_2"), # age comp structure for fishery fleet 1 and 2

                                            FishLenComps_Type =  c("spltRspltS_Year_1-terminal_Fleet_1",
                                                                   "spltRspltS_Year_1-terminal_Fleet_2")  # length comp structure for fishery fleet 1 and 2
  )

  # Setup survey indices and compositions
  input_list <- Setup_Mod_SrvIdx_and_Comps(input_list = input_list,
                                           # data inputs
                                           ObsSrvIdx = sgl_rg_sable_data$ObsSrvIdx,
                                           ObsSrvIdx_SE = sgl_rg_sable_data$ObsSrvIdx_SE,
                                           UseSrvIdx =  sgl_rg_sable_data$UseSrvIdx,
                                           ObsSrvAgeComps = sgl_rg_sable_data$ObsSrvAgeComps,
                                           ISS_SrvAgeComps = sgl_rg_sable_data$ISS_SrvAgeComps,
                                           UseSrvAgeComps = sgl_rg_sable_data$UseSrvAgeComps,
                                           ObsSrvLenComps = sgl_rg_sable_data$ObsSrvLenComps,
                                           UseSrvLenComps = sgl_rg_sable_data$UseSrvLenComps,
                                           ISS_SrvLenComps = sgl_rg_sable_data$ISS_SrvLenComps,

                                           # Model options
                                           srv_idx_type = c("abd", "biom", "abd"), # abundance and biomass for survey fleet 1, 2, and 3
                                           SrvAgeComps_LikeType = c("Multinomial", "none", "Multinomial"), # survey age composition likelihood for survey fleet 1, 2, and 3
                                           SrvLenComps_LikeType = c("Multinomial", "Multinomial", "Multinomial"), #  survey length composition likelihood for survey fleet 1, 2, and 3
                                           SrvAgeComps_Type = c("agg_Year_1-terminal_Fleet_1",
                                                                "none_Year_1-terminal_Fleet_2",
                                                                "agg_Year_1-terminal_Fleet_3"), # survey age comp type

                                           SrvLenComps_Type = c("spltRspltS_Year_1-terminal_Fleet_1",
                                                                "spltRspltS_Year_1-terminal_Fleet_2",
                                                                "spltRspltS_Year_1-terminal_Fleet_3") # survey length comp type
  )

  # Setup fishery selectivity and catchability
  input_list <- Setup_Mod_Fishsel_and_Q(input_list = input_list,
                                        # Model options
                                        # fishery selectivity, whether continuous time-varying
                                        cont_tv_fish_sel = c("none_Fleet_1", "none_Fleet_2"),

                                        # fishery selectivity blocks
                                        fish_sel_blocks = c("Block_1_Year_1-35_Fleet_1", # block 1, fishery ll selex
                                                            "Block_2_Year_36-56_Fleet_1", # block 2 fishery ll selex
                                                            "Block_3_Year_57-terminal_Fleet_1",  # block 3 fishery ll selex
                                                            "none_Fleet_2"), # no blocks for trawl fishery

                                        # fishery selectivity form
                                        fish_sel_model = c("logist1_Fleet_1", "gamma_Fleet_2"),

                                        # fishery catchability blocks
                                        fish_q_blocks = c("Block_1_Year_1-35_Fleet_1", # block 1, fishery ll selex
                                                          "Block_2_Year_36-56_Fleet_1", # block 2 fishery ll selex
                                                          "Block_3_Year_57-terminal_Fleet_1",  # block 3 fishery ll selex
                                                          "none_Fleet_2"), # no blocks for trawl fishery

                                        # whether to estiamte all fixed effects for fishery selectivity
                                        fish_fixed_sel_pars_spec = c("est_all", "est_all"),

                                        # whether to estiamte all fixed effects for fishery catchability
                                        fish_q_spec = c("est_all", "fix") # estiamte fishery q for fleet 1, not for fleet 2
  )

  # mapping for fishery selectivity
  # sharing delta across sexes from early domestic fishery (first time block)
  # also fixing parameters so that no time block for trawl fishery
  input_list$map$fish_fixed_sel_pars <- factor(c(1:7, 2, 8:11, rep(12:13,3), rep(c(14,13),3)))

  # Setup survey selectivity and catchability
  input_list <- Setup_Mod_Srvsel_and_Q(input_list = input_list,

                                       # Model options
                                       # survey selectivity, whether continuous time-varying
                                       cont_tv_srv_sel = c("none_Fleet_1", "none_Fleet_2", "none_Fleet_3"),

                                       # survey selectivity blocks
                                       srv_sel_blocks = c("Block_1_Year_1-56_Fleet_1",  # block 1 for domestic ll survey
                                                          "Block_2_Year_57-terminal_Fleet_1", # block 2 for domestic ll survey
                                                          "none_Fleet_2", "none_Fleet_3"), # no blocks for trawl and jp survey

                                       # survey selectivity form
                                       srv_sel_model = c("logist1_Fleet_1", "exponential_Fleet_2", "logist1_Fleet_3"),

                                       # survey catchability blocks
                                       srv_q_blocks = c("none_Fleet_1", "none_Fleet_2", "none_Fleet_3"),

                                       # whether to estiamte all fixed effects for survey selectivity
                                       srv_fixed_sel_pars_spec = c("est_all", "est_all", "est_all"),

                                       # whether to estiamte all fixed effects for survey catchability
                                       srv_q_spec = c("est_all", "est_all", "est_all"),
                                       t_srv = array(0.5, dim = c(input_list$data$n_regions, input_list$data$n_seas, input_list$data$n_srv_fleets)),

  )

  # ll survey, share delta female (index 2) across time blocks and to the coop jp ll survey delta
  # ll survey, share delta male (index 5) across time blocks and to the coop jp ll survey delta
  # coop jp survey does not estimate parameters and shares deltas with longline survey
  # single time block with trawl survey and only one parameter hence, only one parameter estimated across blocks (indices 7 and 8)
  input_list$map$srv_fixed_sel_pars <- factor(c(1:3, 2, 4:6, 5,rep(7,4), rep(8, 4), rep(c(NA,2), 2), rep(c(NA, 5), 2)))

  # Coop JP Survey (Logistic) Single time block (these estimates are fixed!)
  input_list$par$srv_fixed_sel_pars[1,,,1,3] <- c(0.980660760456, 0.9287775)
  input_list$par$srv_fixed_sel_pars[1,,,2,3] <- c(1.22224502478, 0.8831787)

  # Setup tagging stuff
  input_list <- Setup_Mod_Tagging(input_list = input_list,
                                  use_conv_fish_tagging = c(0,0)
  )

  # set up data weighting stuff
  Wt_FishAgeComps <- array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas,
                                       input_list$data$n_sexes, input_list$data$n_fish_fleets)) # weights for fishery age comps
  Wt_FishAgeComps[1,,,1,1] <- 0.826107286513784 # Weight for fixed gear age comps

  # Fishery length comps
  Wt_FishLenComps <- array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas,
                                       input_list$data$n_sexes, input_list$data$n_fish_fleets)) # weights for fishery age comps
  Wt_FishLenComps[1,,,1,1] <- 4.1837057381917 # Weight for fixed gear len comps females
  Wt_FishLenComps[1,,,2,1] <- 4.26969350917589 # Weight for fixed gear len comps males
  Wt_FishLenComps[1,,,1,2] <- 0.316485920691651 # Weight for trawl gear len comps females
  Wt_FishLenComps[1,,,2,2] <- 0.229396580680981 # Weight for trawl gear len comps males

  # survey age comps
  Wt_SrvAgeComps <- array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas,
                                      input_list$data$n_sexes, input_list$data$n_srv_fleets)) # weights for survey age comps
  Wt_SrvAgeComps[1,,,1,1] <- 3.79224544725927 # Weight for domestic survey ll gear age comps
  Wt_SrvAgeComps[1,,,1,3] <- 1.31681114024037 # Weight for coop jp survey ll gear age comps

  # Survey length comps
  Wt_SrvLenComps <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas,
                                     input_list$data$n_sexes, input_list$data$n_srv_fleets)) # weights for survey age comps
  Wt_SrvLenComps[1,,,1,1] <- 1.43792019016567 # Weight for domestic ll survey len comps females
  Wt_SrvLenComps[1,,,2,1] <- 1.07053763450712 # Weight for domestic ll survey len comps males
  Wt_SrvLenComps[1,,,1,2] <- 0.670883273592302 # Weight for domestic trawl survey len comps females
  Wt_SrvLenComps[1,,,2,2] <- 0.465207132450763 # Weight for domestic trawl survey len comps males
  Wt_SrvLenComps[1,,,1,3] <- 1.27772810174693 # Weight for coop jp ll survey len comps females
  Wt_SrvLenComps[1,,,2,3] <- 0.857519546948587 # Weight for coop jp ll survey len comps males

  input_list <- Setup_Mod_Weighting(input_list = input_list,
                                    Wt_Catch = 50,
                                    Wt_FishIdx = 0.448,
                                    Wt_SrvIdx = 0.448,
                                    Wt_Rec = 1.5,
                                    Wt_F = 0.1,
                                    Wt_Tagging = 0,
                                    Wt_FishAgeComps = Wt_FishAgeComps,
                                    Wt_FishLenComps = Wt_FishLenComps,
                                    Wt_SrvAgeComps = Wt_SrvAgeComps,
                                    Wt_SrvLenComps = Wt_SrvLenComps
  )


  data <- input_list$data
  parameters <- input_list$par
  mapping <- input_list$map

  data$ObsSrvIdx_SE <- data$ObsSrvIdx_SE / data$ObsSrvIdx
  data$ObsFishIdx_SE <- data$ObsFishIdx_SE / data$ObsFishIdx
  parameters$ln_sigmaC[] <- log(sqrt(1/2))
  parameters$ln_sigmaF[] <- log(sqrt(1/2))

  sabie_rtmb_model <- fit_model(data,
                                parameters,
                                mapping,
                                random = NULL,
                                newton_loops = 5,
                                silent = TRUE)

  # Save model results
  sabie_rtmb_model$sd_rep <- RTMB::sdreport(sabie_rtmb_model)

  ssb_expected_vec <- c(
    279.3366, 278.944, 272.1888, 260.9827, 256.9878,
    258.8751, 261.7645, 262.6049, 262.6017, 257.0804,
    248.0993, 236.6576, 221.8412, 201.3895, 187.0074,
    172.6539, 160.2203, 146.5831, 137.7933, 135.2556,
    132.5095, 132.3235, 135.2722, 144.382, 159.5045,
    177.4083, 194.8352, 203.3201, 203.3629, 195.2643,
    182.8816, 169.138, 155.4292, 142.4806, 129.3439,
    118.4992, 110.9616, 105.731, 101.7776, 98.29818,
    95.18359, 92.15907, 91.78957, 93.02909, 95.09047,
    97.5949, 100.6822, 103.3642, 103.8453, 102.4421,
    99.69182, 96.2965, 92.24577, 88.37464, 85.49064,
    84.00741, 82.82673, 82.17328, 82.36075, 85.59484,
    93.89017, 109.6345, 132.3238, 159.8183, 190.1768
  )

  rec_expected_vec <- c(
    28.689986, 30.502832, 32.360868, 33.787049,
    34.202214, 33.736239, 31.821223, 29.274228,
    26.716539, 24.531354, 22.665874, 21.099214,
    19.895769, 19.128259, 19.128786, 9.711081,
    11.140808, 11.136428, 14.105271, 82.809297,
    44.950098, 22.939653, 67.838359, 39.545937,
    13.879768, 15.280121, 25.049890, 9.655648,
    6.736593, 7.540328, 12.276368, 23.223443,
    7.428620, 25.202457, 6.822210, 7.489681,
    12.448327, 20.433272, 9.405259, 36.805104,
    16.090334, 15.672504, 43.113469, 13.167349,
    9.871526, 11.965683, 7.476270, 10.498650,
    9.772347, 15.781236, 21.284253, 10.060608,
    11.389419, 4.465652, 7.521326, 14.132114,
    49.407918, 21.402132, 95.187127, 93.448527,
    42.593962, 81.488045, 25.248837, 28.146785,
    25.904198
  )

  expect_equal(sabie_rtmb_model$rep$SSB[1,1,], ssb_expected_vec, tolerance = 1e-2)
  expect_equal(sabie_rtmb_model$rep$Rec[1,1,], rec_expected_vec, tolerance = 1e-2)
  expect_true(sabie_rtmb_model$sd_rep$pdHess)

})

