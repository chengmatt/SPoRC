library(SPoRC)
library(testthat)
data("three_rg_sable_data")

test_that("Three-region Sablefish RTMB model produces expected results", {

  # Initialize model dimensions and data list
  input_list <- Setup_Mod_Dim(years = 1:length(three_rg_sable_data$years),
                              # vector of years (1 - 62)
                              ages = 1:length(three_rg_sable_data$ages),
                              # vector of ages (1 - 30)
                              lens = three_rg_sable_data$lens,
                              # number of lengths (41 - 99)
                              n_regions = three_rg_sable_data$n_regions,
                              # number of regions (5)
                              n_sexes = three_rg_sable_data$n_sexes,
                              # number of sexes (2)
                              n_fish_fleets = three_rg_sable_data$n_fish_fleets,
                              # number of fishery fleet (2)
                              n_srv_fleets = three_rg_sable_data$n_srv_fleets,
                              # number of survey fleets (2)
                              n_seas = three_rg_sable_data$n_seas,
                              # number of seasons
                              n_pop = three_rg_sable_data$n_pop,
                              natal_region = three_rg_sable_data$natal_region,
                              # population stuff
                              verbose = FALSE
  )

  # Setup recruitment stuff (using defaults for other stuff)
  input_list <- Setup_Mod_Rec(input_list = input_list, # input data list from above
                              do_rec_bias_ramp = 0, # not using bias ramp
                              sigmaR_switch = 16, # switch to using late sigma in year 16
                              dont_est_recdev_last = 1, # don't estimate last rec dev
                              # Model options
                              rec_model = "mean_rec", # recruitment model
                              sigmaR_spec = "fix", # fixing
                              InitDevs_spec = "est_shared_r",
                              # initial deviations are shared across regions,
                              # but recruitment deviations are region specific
                              ln_sigmaR = array(log(c(0.4, 1.2)), dim = c(2, input_list$data$n_pop, input_list$data$n_regions)),
                              # values to fix sigmaR at, or starting values
                              ln_global_R0 = log(20),
                              # starting value for global R0
                              rec_region_prop_pars = array(c(0.2, 0.2, 0.2, 0.2), dim = c(input_list$data$n_pop, input_list$data$n_regions - 1))
                              # starting value for R0 proportions in multinomial logit space
  )

  # Setup biological stuff (using defaults for other stuff)
  input_list <- Setup_Mod_Biologicals(input_list = input_list,
                                      WAA = three_rg_sable_data$WAA, # weight at age
                                      MatAA = three_rg_sable_data$MatAA, # maturity at age
                                      AgeingError = three_rg_sable_data$AgeingError,
                                      # ageing error matrix
                                      fit_lengths = 1, # fitting lengths
                                      SizeAgeTrans = three_rg_sable_data$SizeAgeTrans,
                                      # size age transition matrix
                                      M_spec = "fix", # fix natural mortality
                                      Fixed_natmort = array(0.104884, dim = c(three_rg_sable_data$n_pop,
                                                                              three_rg_sable_data$n_regions,
                                                                              length(three_rg_sable_data$years),
                                                                              length(three_rg_sable_data$ages),
                                                                              three_rg_sable_data$n_sexes))
                                      # values to fix natural mortality at
  )

  # setting up movement parameterization
  Movement_prior <- expand.grid(
    pop = 1, # populations
    region_from = 1:3, # regions
    year = 1, # penalize first year since no blocks
    seas = 1,
    age = c(6,7,16), # age blocks
    sex = 1, # sex
    alpha = I(list(rep(3, 3))) # prior alpha to each row
  )

  input_list <- Setup_Mod_Movement(input_list = input_list,
                                   # Model options
                                   Movement_ageblk_spec = list(c(1:6), c(7:15), c(16:30)),
                                   # estimating movement in 3 age blocks
                                   # (ages 1-6, ages 7-15, ages 16-30)
                                   Movement_popblk_spec = 'constant', # population-invariant movement
                                   Movement_yearblk_spec = "constant", # time-invariant movement
                                   Movement_sexblk_spec = "constant", # sex-invariant movement
                                   Movement_seasblk_spec = 'constant',
                                   do_recruits_move = 0, # recruits do not move
                                   use_fixed_movement = 0, # estimating movement
                                   Use_Movement_Prior = 1, # priors used for movement
                                   Movement_prior = Movement_prior
                                   # vague prior to penalize movement away from the extremes
  )

  # setting up tagging parameterization

  # setup tagging priors
  tag_prior <- data.frame(
    region = 1,
    block = c(1,2),
    fleet = 1,
    mu = NA, # no mean, since symmetric beta
    sd = 5, # sd = 5
    type = 0 # symmetric beta
  )

  input_list <- Setup_Mod_Tagging(input_list = input_list,
                                  use_conv_fish_tagging = c(1,0), # using tagging data for fixed gear
                                  conv_tag_max_liberty = 15, # maximum number of years to track a cohort

                                  # Data Inputs
                                  conv_tag_release_indicator = three_rg_sable_data$conv_tag_release_indicator,
                                  # tag release indicator (first col = tag region,
                                  # second col = tag year),
                                  # total number of rows = number of tagged cohorts
                                  conv_tagged_fish = three_rg_sable_data$conv_tagged_fish, # Released fish
                                  # dimensioned by total number of tagged cohorts, (implicitly
                                  # tracks the release year and region), pop, age, and sex
                                  obs_conv_tag_fish_recap = three_rg_sable_data$obs_conv_tag_fish_recap,
                                  # dimensioned by max tag liberty, tagged cohorts, pop, regions,
                                  # ages, and sexes

                                  # Model options
                                  conv_fish_tag_like = "NegBin", # Negative Binomial
                                  conv_tag_mixing_period = 2, # Don't fit tagging until release year + 1
                                  conv_tag_t_tagging = 0.5, # tagging happens midway through the year,
                                  # movement does not occur within that year
                                  use_conv_tag_fishrep_prior = 1, # tag reporting rate priors are used
                                  conv_tag_fishrep_prior = tag_prior,
                                  conv_tag_age_pool = list(c(1:6), c(7:15), c(16:30)), # whether or
                                  # not to pool tagging data when fitting (for computational cost)
                                  conv_tag_sex_pool = list(c(1:2)), # whether or not to pool
                                  # sex-specific data when fitting
                                  init_conv_tag_mort_spec = "fix", # fixing initial tag mortality
                                  conv_tag_shed_spec = "fix", # fixing chronic shedding
                                  conv_tagrep_spec = "est_shared_r_f", # tag reporting rates are
                                  # not region specific
                                  # Time blocks for tag reporting rates
                                  conv_tag_fish_reporting_blocks = c(
                                    apply(expand.grid(1:input_list$data$n_regions, 1:input_list$data$n_fish_fleets), 1, function(x)
                                      paste0("Block_1_Year_1-35_Region_", x[1], "_Fleet_", x[2])),
                                    apply(expand.grid(1:input_list$data$n_regions, 1:input_list$data$n_fish_fleets), 1, function(x)
                                      paste0("Block_2_Year_36-terminal_Region_", x[1], "_Fleet_", x[2]))
                                  ),
                                  conv_fish_tag_attr = 'p_a_s',

                                  # Specify starting values or fixing values
                                  ln_init_conv_tag_mort = log(0.1), # fixing initial tag mortality
                                  ln_conv_tag_shed = log(0.02),  # fixing tag shedding
                                  ln_conv_fish_tag_theta = log(0.5),
                                  # starting value for tagging overdispersion
                                  conv_tag_fish_reporting_pars = array(log(0.2 / (1-0.2)), dim = c(input_list$data$n_regions, 2, input_list$data$n_fish_fleets))
                                  # starting values for tag reporting pars

  )

  input_list$par$conv_tag_fish_reporting_pars
  input_list$par$ln_conv_tag_shed
  input_list$par$ln_conv_fish_tag_theta
  input_list$par$ln_init_conv_tag_mort

  input_list$map$conv_tag_fish_reporting_pars
  input_list$par$ln_conv_tag_shed
  input_list$par$ln_conv_fish_tag_theta
  input_list$par$ln_init_conv_tag_mort


  # setting up catch data
  input_list <- Setup_Mod_Catch_and_F(input_list = input_list,
                                      # Data inputs
                                      ObsCatch = three_rg_sable_data$ObsCatch,
                                      UseCatch = three_rg_sable_data$UseCatch,
                                      # Model options
                                      Use_F_pen = 1,
                                      # whether to use f penalty, == 0 don't use, == 1 use
                                      sigmaC_spec = 'fix',
                                      ln_sigmaC =
                                        array(log(0.05), dim = c(input_list$data$n_regions,
                                                                 length(input_list$data$years),
                                                                 input_list$data$n_seas,
                                                                 input_list$data$n_fish_fleets)),
                                      # fixing catch sd at small value
                                      ln_F_mean = array(-2, dim = c(input_list$data$n_regions,
                                                                    input_list$data$n_seas,
                                                                    input_list$data$n_fish_fleets))
                                      # some starting values for fishing mortality
  )

  # Fishery Indices and Compositions
  input_list <- Setup_Mod_FishIdx_and_Comps(input_list = input_list,
                                            # data inputs
                                            ObsFishIdx = three_rg_sable_data$ObsFishIdx,
                                            ObsFishIdx_SE = three_rg_sable_data$ObsFishIdx_SE,
                                            UseFishIdx =  three_rg_sable_data$UseFishIdx,
                                            ObsFishAgeComps = three_rg_sable_data$ObsFishAgeComps,
                                            UseFishAgeComps = three_rg_sable_data$UseFishAgeComps,
                                            ISS_FishAgeComps = three_rg_sable_data$ISS_FishAgeComps,
                                            ObsFishLenComps = three_rg_sable_data$ObsFishLenComps,
                                            UseFishLenComps = three_rg_sable_data$UseFishLenComps,
                                            ISS_FishLenComps = three_rg_sable_data$ISS_FishLenComps,

                                            # Model options
                                            fish_idx_type = c("none", "none"),
                                            # fishery indices not used
                                            FishAgeComps_LikeType =
                                              c("Multinomial", "none"),
                                            # age comp likelihoods for fishery fleet 1 and 2
                                            FishLenComps_LikeType =
                                              c("Multinomial", "Multinomial"),
                                            # length comp likelihoods for fishery fleet 1 and 2
                                            FishAgeComps_Type =
                                              c("spltRjntS_Year_1-terminal_Fleet_1",
                                                "none_Year_1-terminal_Fleet_2"),
                                            # age comp structure for fishery fleet 1 and 2
                                            FishLenComps_Type =
                                              c("spltRjntS_Year_1-terminal_Fleet_1",
                                                "spltRjntS_Year_1-terminal_Fleet_2")
                                            # length comp structure for fishery fleet 1 and 2
  )

  # Survey Indices and Compositions
  input_list <- Setup_Mod_SrvIdx_and_Comps(input_list = input_list,
                                           # data inputs
                                           ObsSrvIdx = three_rg_sable_data$ObsSrvIdx,
                                           ObsSrvIdx_SE = three_rg_sable_data$ObsSrvIdx_SE,
                                           UseSrvIdx =  three_rg_sable_data$UseSrvIdx,
                                           ObsSrvAgeComps = three_rg_sable_data$ObsSrvAgeComps,
                                           ISS_SrvAgeComps = three_rg_sable_data$ISS_SrvAgeComps,
                                           UseSrvAgeComps = three_rg_sable_data$UseSrvAgeComps,
                                           ObsSrvLenComps = three_rg_sable_data$ObsSrvLenComps,
                                           UseSrvLenComps = three_rg_sable_data$UseSrvLenComps,
                                           ISS_SrvLenComps = three_rg_sable_data$ISS_SrvLenComps,

                                           # Model options
                                           srv_idx_type = c("abd", "abd"),
                                           # abundance and biomass for survey fleet 1 and 2
                                           SrvAgeComps_LikeType =
                                             c("Multinomial", "Multinomial"),
                                           # survey age composition likelihood for survey fleet
                                           # 1, and 2
                                           SrvLenComps_LikeType =
                                             c("none", "none"),
                                           #  no length compositions used for survey
                                           SrvAgeComps_Type = c("spltRjntS_Year_1-terminal_Fleet_1",
                                                                "spltRjntS_Year_1-terminal_Fleet_2"),
                                           # survey age comp type
                                           SrvLenComps_Type = c("none_Year_1-terminal_Fleet_1",
                                                                "none_Year_1-terminal_Fleet_2")
                                           # survey length comp type
  )

  # Fishery Selectivity and Catchability
  input_list <- Setup_Mod_Fishsel_and_Q(input_list = input_list,

                                        # Model options
                                        cont_tv_fish_sel = c("none_Fleet_1", "none_Fleet_2"),
                                        # fishery selectivity, whether continuous time-varying

                                        # fishery selectivity blocks
                                        fish_sel_blocks =
                                          c("Block_1_Year_1-56_Fleet_1",
                                            # block 1, fishery ll selex
                                            "Block_2_Year_57-terminal_Fleet_1",
                                            # block 3 fishery ll selex
                                            "none_Fleet_2"),
                                        # no blocks for trawl fishery

                                        # fishery selectivity form
                                        fish_sel_model =
                                          c("logist1_Fleet_1",
                                            "gamma_Fleet_2"),

                                        # fishery catchability blocks
                                        fish_q_blocks =
                                          c("none_Fleet_1",
                                            "none_Fleet_2"),
                                        # no blocks since q is not estimated

                                        # whether to estimate all fixed effects
                                        # for fishery selectivity and later modify
                                        # to fix and share parameters
                                        fish_fixed_sel_pars_spec =
                                          c("est_all", "est_all"),

                                        # whether to estimate all fixed effects
                                        # for fishery catchability
                                        fish_q_spec =
                                          c("fix", "fix")
                                        # fix fishery q since not used
  )


  # Custom parameter sharing for fishery selectivity
  map_ln_fish_fixed_sel_pars <- input_list$par$ln_fish_fixed_sel_pars # mapping fishery selectivity

  # Fixed gear fleet, unique parameters for each sex (time block 1)
  map_ln_fish_fixed_sel_pars[,1,1,1,1] <- 1 # a50, female, time block 1, fixed gear
  map_ln_fish_fixed_sel_pars[,2,1,1,1] <- 2 # delta, female, time block 1, fixed gear
  # (shared with time block 2 and sex)
  map_ln_fish_fixed_sel_pars[,1,1,2,1] <- 3 # a50, male, time block 1, fixed gear
  map_ln_fish_fixed_sel_pars[,2,1,2,1] <- 2 # delta, male, time block 1, fixed gear
  # (shared with time block 2 and sex)

  # time block 2, fixed gear fishery
  map_ln_fish_fixed_sel_pars[,1,2,1,1] <- 4 # a50, female, time block 2, fixed gear
  map_ln_fish_fixed_sel_pars[,2,2,1,1] <- 2 # delta, female, time block 2, fixed gear
  # (shared with time block 1 and sex)
  map_ln_fish_fixed_sel_pars[,1,2,2,1] <- 5 # a50, male, time block 2, fixed gear
  map_ln_fish_fixed_sel_pars[,2,2,2,1] <- 2 # delta, male, time block 2, fixed gear
  # (shared with time block 1 and sex)

  # time block 1 and 2, trawl gear fishery
  map_ln_fish_fixed_sel_pars[,1,1,1,2] <- 6 # amax, female, time block 1, trawl gear
  map_ln_fish_fixed_sel_pars[,2,1,1,2] <- 7 # delta, female, time block 1, trawl gear
  # (shared by sex)
  map_ln_fish_fixed_sel_pars[,1,1,2,2] <- 8 # amax, male, time block 1, trawl gear
  map_ln_fish_fixed_sel_pars[,2,1,2,2] <- 7 # delta, male, time block 1, trawl gear
  # (shared by sex)
  map_ln_fish_fixed_sel_pars[,,2,,2] <- NA # no parameters estimated for time block 2 trawl gear

  input_list$map$ln_fish_fixed_sel_pars <- factor(map_ln_fish_fixed_sel_pars) # input into map list
  input_list$par$ln_fish_fixed_sel_pars[,,,,1] <- log(3) # some more inforamtive starting values
  input_list$par$ln_fish_fixed_sel_pars[,,,,2] <- log(6) # some more inforamtive starting values

  input_list <- Setup_Mod_Srvsel_and_Q(input_list = input_list,

                                       # Model options
                                       # survey selectivity, whether continuous time-varying
                                       cont_tv_srv_sel =
                                         c("none_Fleet_1",
                                           "none_Fleet_2"),

                                       # survey selectivity blocks
                                       srv_sel_blocks =
                                         c("none_Fleet_1",
                                           "none_Fleet_2"
                                         ), # no blocks for jp and domestic survey

                                       # survey selectivity form
                                       srv_sel_model =
                                         c("logist1_Fleet_1",
                                           "logist1_Fleet_2"),

                                       # survey catchability blocks
                                       srv_q_blocks =
                                         c("none_Fleet_1",
                                           "none_Fleet_2"),

                                       # whether to estiamte all fixed effects
                                       # for survey selectivity and later
                                       # modify to fix/share parameters
                                       srv_fixed_sel_pars_spec =
                                         c("est_all",
                                           "est_all"),

                                       # whether to estiamte all
                                       # fixed effects for survey catchability
                                       # spatially-invariant q
                                       srv_q_spec =
                                         c("est_shared_r",
                                           "est_shared_r"),

                                       # Starting values for survey catchability
                                       ln_srv_q = array(8.75,
                                                        dim = c(input_list$data$n_regions, 1,
                                                                input_list$data$n_srv_fleets))
  )

  # Custom mapping survey selectivity stuff
  map_ln_srv_fixed_sel_pars <- input_list$par$ln_srv_fixed_sel_pars # set up mapping factor stuff

  # Coop survey (japanese)
  map_ln_srv_fixed_sel_pars[,1,1,1,1] <- 1 # a50, coop survey, time block 1, female
  map_ln_srv_fixed_sel_pars[,2,1,1,1] <- 2 # delta, coop survey, time block 1, female
  # (sharing with domestic survey)
  map_ln_srv_fixed_sel_pars[,1,1,2,1] <- 3 # a50, coop survey, time block 1, male
  map_ln_srv_fixed_sel_pars[,2,1,2,1] <- 4 # delta, coop survey, time block 1, male
  # (sharing with domestic survey)

  # domestic survey
  map_ln_srv_fixed_sel_pars[,1,1,1,2] <- 5 # a50, domestic survey, time block 1, female
  map_ln_srv_fixed_sel_pars[,2,1,1,2] <- 2 # delta, domestic survey, time block 1, female
  # (sharing with coop survey)
  map_ln_srv_fixed_sel_pars[,1,1,2,2] <- 6 # a50, domestic survey, time block 1, male
  map_ln_srv_fixed_sel_pars[,2,1,2,2] <- 4 # delta, domestic survey, time block 1, male
  # (sharing with coop survey)

  input_list$map$ln_srv_fixed_sel_pars <- factor(map_ln_srv_fixed_sel_pars)  # input into map list
  input_list$par$ln_srv_fixed_sel_pars[] <- log(3) # some more informative starting values

  # set up model weighting stuff
  input_list <- Setup_Mod_Weighting(input_list = input_list,
                                    Wt_Catch = 1,
                                    Wt_FishIdx = 1,
                                    Wt_SrvIdx = 1,
                                    Wt_Rec = 1,
                                    Wt_F = 1,
                                    Wt_Tagging = 1,
                                    # Composition model weighting
                                    Wt_FishAgeComps =
                                      array(1, dim = c(input_list$data$n_regions,
                                                       length(input_list$data$years),
                                                       input_list$data$n_seas,
                                                       input_list$data$n_sexes,
                                                       input_list$data$n_fish_fleets)),
                                    Wt_FishLenComps =
                                      array(1, dim = c(input_list$data$n_regions,
                                                       length(input_list$data$years),
                                                       input_list$data$n_seas,
                                                       input_list$data$n_sexes,
                                                       input_list$data$n_fish_fleets)),
                                    Wt_SrvAgeComps =
                                      array(1, dim = c(input_list$data$n_regions,
                                                       length(input_list$data$years),
                                                       input_list$data$n_seas,
                                                       input_list$data$n_sexes,
                                                       input_list$data$n_srv_fleets)),
                                    Wt_SrvLenComps =
                                      array(1, dim = c(input_list$data$n_regions,
                                                       length(input_list$data$years),
                                                       input_list$data$n_seas,
                                                       input_list$data$n_sexes,
                                                       input_list$data$n_srv_fleets))
  )

  # extract out lists updated with helper functions
  data <- input_list$data
  parameters <- input_list$par
  mapping <- input_list$map

  # set survey timing to midway through year
  data$t_srv[] <- 0.5


  sabie_rtmb_model <- fit_model(data,
                                parameters,
                                mapping,
                                random = NULL,
                                newton_loops = 3,
                                silent = TRUE)

  # Save model results
  sabie_rtmb_model$sd_rep <- RTMB::sdreport(sabie_rtmb_model)

  ssb_expected_vec <- c(39.67343, 107.48289,  92.51175,  39.36165, 106.71604,  91.67044,
                        37.88366, 103.36530,  87.99805,  35.57364,  97.88692,  82.22800,
                        34.04447,  94.17644,  79.61764,  33.74517,  92.57153,  78.73323,
                        33.41889,  90.98589,  77.99490,  32.64240,  88.42775,  76.22456,
                        31.78232,  85.57802,  74.15394,  30.03818,  80.89521,  70.19132,
                        27.93790,  75.43623,  65.32746,  26.02643,  69.99544,  59.61421,
                        23.28865,  63.50973,  53.76333,  19.50736,  55.15801,  46.48303,
                        17.46065,  49.59485,  41.13676,  15.64010,  44.32800,  35.81326,
                        14.16076,  39.77291,  31.34789,  12.36818,  34.71830,  26.62317,
                        10.90454,  32.27983,  23.00080,  10.91353,  30.47202,  22.40091,
                        11.54200,  28.68866,  21.06338,  12.78997,  27.76512,  20.73099,
                        15.80127,  27.80515,  20.51699,  21.35350,  30.42202,  21.49741,
                        28.84921,  35.58794,  23.17528,  34.73732,  42.89803,  27.42985,
                        37.66594,  51.75827,  34.30106,  40.48605,  58.91989,  38.85329,
                        42.07947,  64.73477,  41.99223,  39.42727,  69.39803,  45.47265,
                        34.64654,  72.01091,  48.48852,  30.60451,  71.25642,  48.79042,
                        26.15117,  69.15636,  48.64726,  22.51175,  65.95486,  47.16957,
                        20.07105,  61.19376,  44.34020,  18.50823,  57.28018,  41.11739,
                        17.10786,  53.98262,  38.45956,  16.05590,  51.14422,  36.64825,
                        15.19210,  48.33164,  35.70321,  14.64602,  45.77870,  34.39098,
                        14.22069,  43.65331,  33.11230,  14.15888,  41.74608,  31.15173,
                        14.61611,  41.47845,  29.93888,  15.11204,  42.25352,  29.76285,
                        16.10109,  43.47449,  29.65799,  16.85709,  44.94767,  30.17051,
                        17.38703,  46.23232,  31.30070,  17.87425,  47.22247,  32.09435,
                        16.94653,  47.38386,  32.96218,  15.86883,  46.74063,  33.19290,
                        14.96649,  45.53138,  32.95381,  14.31236,  44.27862,  32.66366,
                        14.01417,  42.63895,  31.70690,  13.77220,  40.83854,  30.49004,
                        13.71247,  38.61034,  28.70865,  13.52290,  36.57973,  27.49647,
                        13.29165,  34.87613,  26.32353,  13.10033,  33.47268,  25.47910,
                        14.17759,  33.13018,  24.49356,  16.70203,  34.41596,  23.92235,
                        21.43838,  38.56011,  24.77585,  28.23219,  46.42713,  27.54328)

  rec_expected_vec <- c( 9.4882457,  5.1430092,  2.2192742,  9.6630398,  5.1961050,
                         2.2291811, 10.0958374,  5.2967932,  2.2484061, 10.9922214,
                         5.4635005,  2.2799843, 12.2623279,  5.6435165,  2.3130034,
                         13.5242196,  5.7626878,  2.3360275, 14.5187579,  5.8033033,
                         2.3465295, 14.2465191,  5.6857193,  2.3322683, 13.2368709,
                         5.4666756,  2.3007542, 11.9777172,  5.1733162,  2.2539899,
                         10.7790350,  4.8350114,  2.1950969,  9.6618306,  4.4499041,
                         2.1204772,  8.6273058,  4.0201232,  2.0285796,  7.8720899,
                         3.5973555,  1.9338144,  7.7642927,  3.2218101,  1.8504594,
                         2.9546454,  0.7340996,  0.6051922,  4.2861911,  0.6259843,
                         0.5790964,  4.2950662,  0.5423907,  0.5643884,  4.4130149,
                         0.6129132,  0.6581187, 37.3789206,  0.9342916, 16.2418030,
                         58.7085464,  1.1014723,  1.6354865,  3.3162862,  0.8741255,
                         0.6377696, 23.3252522,  1.0438339,  0.6227757, 74.8103258,
                         1.0199767,  0.6415209, 37.3534299,  1.1156789,  0.7737032,
                         3.3670105,  1.5305303,  0.7793424, 22.5914344,  8.4847774,
                         0.7974638,  2.9121990,  7.2274555,  0.7033818,  1.2765954,
                         3.2446704,  0.8076420,  1.0703113, 12.0432368,  1.5424296,
                         0.9441979,  2.5853643,  7.8793986,  1.7288712,  5.3667827,
                         12.0724809,  3.7518065,  6.7730932,  1.3388498,  0.8745027,
                         6.8937280,  2.0412569,  0.7034009,  8.7952875,  3.2740892,
                         0.9464316,  2.9811932,  1.3939773,  3.7067847,  3.1176189,
                         3.3451741,  7.4717484,  6.1761534,  1.8894439,  3.4956614,
                         4.4658538,  0.8986409,  8.2310130, 20.9702710,  0.9209671,
                         1.7131843, 19.7925575,  0.9541103,  2.8619230,  5.6611189,
                         1.8745978, 16.9200019, 16.3403548,  1.0912541,  4.2868211,
                         6.8301796,  0.8707011,  2.6256240,  3.6489298,  1.5503650,
                         3.5679318,  6.8753391,  1.7678221,  2.4910274,  3.0371900,
                         1.0510983,  2.7897915,  6.2499170,  2.9613137,  1.6673446,
                         2.7244196,  2.0082918,  5.6112075,  3.2792355,  4.5251297,
                         8.7084097,  4.5978793,  1.1628381,  3.1492820,  1.4751992,
                         0.7748086,  4.3467724,  1.4351274,  1.1514011,  1.0543084,
                         1.1591482,  1.5131083,  3.1544606,  2.0511119,  1.6548438,
                         5.1380380,  4.1908328,  0.8485258, 32.6601177, 20.9888894,
                         1.1995712,  5.5172590,  7.9104128,  1.9855286, 74.7348313,
                         51.2124053,  5.2986286, 16.9730548, 23.2022554,  0.6042292,
                         22.2523710, 11.7528461,  0.5012960, 10.4172817,  5.6120456,
                         2.4140437)

  expect_equal(as.vector(sabie_rtmb_model$rep$SSB), ssb_expected_vec, tolerance = 1e-3)
  expect_equal(as.vector(sabie_rtmb_model$rep$Rec), rec_expected_vec, tolerance = 1e-3)
  expect_true(sabie_rtmb_model$sd_rep$pdHess)

})

