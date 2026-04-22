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
                              use_rec_region_prop_prior = 1,
                              rec_region_prop_prior = data.frame(pop = 1, alpha = I(list(rep(3, input_list$data$n_regions)))),
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
                                  use_conv_fish_tagging = c(1, 0), # using tagging data for fixed gear
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
                                  conv_fish_tag_like = "Multinomial_Release", # Negative Binomial
                                  conv_tag_mixing_period = 2, # Don't fit tagging until release year + 1
                                  conv_tag_t_tagging = 0.5, # tagging happens midway through the year,
                                  # movement does not occur within that year
                                  use_conv_tag_fishrep_prior = 1, # tag reporting rate priors are used
                                  conv_tag_fishrep_prior = tag_prior,
                                  conv_tag_age_pool = as.list(1:30), # whether or
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

  # setting up catch data
  suppressWarnings(
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

  # defining priors
  sex_par <- expand.grid(sex = 1:2, par = 1:2)
  fleet_blocks <- data.frame(
    fleet = c(1, 2),
    block = 1
  )

  # merge together (note that unlike the operational assessment, selectivity
  # blocks are reduced from 3 to 2)
  fish_selex_structure <- merge(fleet_blocks, sex_par)

  # Merge to get all valid combinations
  fish_selex_structure <- merge(fleet_blocks, sex_par) %>%
    dplyr::filter(!(fleet == 1 & block == 1 & sex == 2 & par == 2)) %>%              # remove priors for any unestimated pars -- par1=a50, par2=delta; NEEDS TO MATCH PARAMETER MAPPING
    dplyr::filter(!(fleet == 2 & block == 1 & sex == 2 & par == 1))                  # remove priors for any unestimated pars -- par1=a50, par2=delta; NEEDS TO MATCH PARAMETER MAPPING

  # Add the lognormal prior values - creates a dataframe, each row is a unique parameter combination to apply the prior to
  fish_selex_prior <- cbind(
    region = 1,
    fish_selex_structure,
    mu = 2,                                                                      # All selex means = 1 (means should be defined in normal space)
    sd = 3                                                                       # All selex sd = 5
  )

  fish_selex_prior_tf <- fish_selex_prior %>%                                    # set tighter selex prior for TF
    dplyr::filter((fleet == 2 & par == 1)) %>%
    dplyr::mutate(mu = 2, sd = 1) %>%
    dplyr::full_join(fish_selex_prior %>%  dplyr::filter(!(fleet == 2 & par == 1 )))

  fish_selex_prior_tf <- fish_selex_prior_tf %>%                                    # set tighter selex prior for TF
    dplyr::filter((fleet == 2 & par == 2)) %>%
    dplyr::mutate(mu = 5, sd = 2) %>%
    dplyr::full_join(fish_selex_prior_tf %>%  dplyr::filter(!(fleet == 2 & par == 2)))

  input_list <- Setup_Mod_Fishsel_and_Q(input_list = input_list,

                                        # Model options
                                        cont_tv_fish_sel = c("none_Fleet_1", "none_Fleet_2"),
                                        # fishery selectivity, whether continuous time-varying

                                        # fishery selectivity blocks
                                        fish_sel_blocks =
                                          c("none_Fleet_1",
                                            "none_Fleet_2"),
                                        # no blocks for trawl fishery

                                        # fishery selectivity form
                                        fish_sel_model =
                                          c("logist1_Fleet_1", "gamma_Fleet_2"),

                                        # fishery catchability blocks
                                        fish_q_blocks =
                                          c("none_Fleet_1", "none_Fleet_2"),
                                        # no blocks since q is not estimated

                                        # sharing fishery selex parameters
                                        fish_fixed_sel_pars =
                                          c("est_shared_r", "est_shared_r"),

                                        # whether to estimate all fixed effects
                                        # for fishery catchability
                                        fish_q_spec =
                                          c("fix", "fix"),
                                        Use_fish_selex_prior = 1,
                                        fish_selex_prior = fish_selex_prior
  )


  # setup survey selectivity
  # Define sex and parameter combinations
  sex_par <- expand.grid(sex = 1:2, par = 1:2)

  # Define valid fleet-block combinations (only estimating domestic and jp LLS)
  fleet_blocks <- data.frame(
    fleet = c(1, 2),
    block = c(1, 1)
  )

  # Merge to get all valid combinations
  srv_selex_structure <- merge(fleet_blocks, sex_par)

  # Add the lognormal prior values - creates a dataframe, each row is a unique parameter combination to apply the prior to
  srv_selex_prior <- cbind(
    region = 1,
    srv_selex_structure,
    mu = 1,
    sd = 5
  ) %>%
    filter(!(fleet == 2 & par == 2 & sex == 2)) %>%
    mutate(mu = ifelse(fleet == 2, 2, mu),
           sd = ifelse(fleet == 2, 3, sd))

  input_list <- Setup_Mod_Srvsel_and_Q(input_list = input_list,

                                       # Model options
                                       # survey selectivity, whether continuous time-varying
                                       cont_tv_srv_sel =
                                         c("none_Fleet_1",
                                           "none_Fleet_2"
                                         ),

                                       # survey selectivity blocks
                                       srv_sel_blocks =                          # survey selectivity time blocks if not TV specified above for a given fleet
                                         c("none_Fleet_1",
                                           "none_Fleet_2"                        # No blocks for JPN LLS
                                         ),

                                       # survey selectivity form
                                       srv_sel_model =
                                         c("logist1_Fleet_1",
                                           "logist1_Fleet_2"
                                         ),

                                       # survey catchability blocks
                                       srv_q_blocks =
                                         c("none_Fleet_1",
                                           "none_Fleet_2"
                                         ),

                                       # whether to estiamte all fixed effects
                                       # for survey selectivity and later
                                       # modify to fix/share parameters
                                       srv_fixed_sel_pars_spec =
                                         c("est_shared_r",
                                           "est_shared_r"
                                         ),

                                       # whether to estiamte all
                                       # fixed effects for survey catchability
                                       # spatially-invariant q
                                       srv_q_spec =
                                         c("est_shared_r",
                                           "est_shared_r"
                                         ),
                                       Use_srv_selex_prior = 1,
                                       srv_selex_prior = srv_selex_prior
  )

  # set up model weighting stuff
  input_list <- Setup_Mod_Weighting(input_list = input_list,
                                    Wt_Catch = 1,
                                    Wt_FishIdx = 1,
                                    Wt_SrvIdx = 1,
                                    Wt_Rec = 1,
                                    Wt_F = 1,
                                    Wt_Tagging = 0.5,
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


  # Additional Model Specifications -----------------------------------------

  # Survey Ages (~100 total across all regions)
  data$ISS_SrvAgeComps[] <- 33

  # Fishery Ages
  data$ISS_FishAgeComps[1,,,,] <- 25  # BS
  data$ISS_FishAgeComps[2,,,,] <- 25  # AI
  data$ISS_FishAgeComps[3,,,,] <- 50  # GOA

  # Fishery Lengths - Fixed Gear
  data$ISS_FishLenComps[1,,,,1] <- 13  # BS
  data$ISS_FishLenComps[2,,,,1] <- 13  # AI
  data$ISS_FishLenComps[3,,,,1] <- 18  # GOA

  # Fishery Lengths - Trawl Gear
  data$ISS_FishLenComps[1,,,,1] <- 25  # BS
  data$ISS_FishLenComps[2,,,,1] <- 10  # AI
  data$ISS_FishLenComps[3,,,,1] <- 10  # GOA

  # Map off early delta for fishery
  map_fish_fixed <- array(mapping$ln_fish_fixed_sel_pars, dim = dim(parameters$ln_fish_fixed_sel_pars))
  map_fish_fixed[,2,1,2,1]  <- map_fish_fixed[,2,1,1,1] # share deltas

  # Map off bmax for trawl females
  map_fish_fixed[,1,1,2,2]  <- map_fish_fixed[,1,1,1,2] # share deltas
  mapping$ln_fish_fixed_sel_pars <- factor(map_fish_fixed)

  # Map off delta for JP LLS
  map_srv_fixed <- array(mapping$ln_srv_fixed_sel_pars, dim = dim(parameters$ln_srv_fixed_sel_pars))
  map_srv_fixed[,2,1,2,2]  <- map_srv_fixed[,2,1,1,2] # share deltas
  mapping$ln_srv_fixed_sel_pars <- factor(map_srv_fixed)

  # Some starting values to help out the model
  parameters$ln_srv_fixed_sel_pars[] <- log(2)
  parameters$ln_fish_fixed_sel_pars[,,,,1] <- log(2) # fixed gear
  parameters$ln_fish_fixed_sel_pars[,,,,2] <- log(5) # trawl gear

  # make AD model function
  # Fit model
  st <- Sys.time()
  sabie_rtmb_model <- fit_model(data,
                                parameters,
                                mapping,
                                random = NULL,
                                newton_loops = 3,
                                silent = FALSE
  )
  en <- Sys.time()
  print(en - st)

  sd_rep <- sdreport(sabie_rtmb_model)
  rep <- sabie_rtmb_model$rep

  ssb_expected_vec <- c(42.43286, 113.10831, 93.93930, 42.11580, 112.37242, 93.16166,
                        40.61397, 109.02433, 89.61442, 38.24779, 103.51031, 84.02818,
                        36.57228, 99.69266, 81.48621, 36.14607, 97.99478, 80.60904,
                        35.62111, 96.25595, 79.80418, 34.57893, 93.46656, 77.91364,
                        33.39338, 90.29843, 75.66397, 31.27087, 85.18115, 71.48430,
                        28.82506, 79.21385, 66.37307, 26.71808, 73.28888, 60.45717,
                        23.83731, 66.27699, 54.37803, 20.01292, 57.47868, 46.92137,
                        18.06954, 51.70478, 41.51765, 16.37339, 46.34164, 36.25682,
                        14.96237, 41.74213, 31.91090, 13.13598, 36.61392, 27.31414,
                        11.81574, 34.08713, 23.83993, 11.71187, 32.24100, 23.31492,
                        12.50053, 30.41056, 22.09590, 13.99514, 29.55602, 21.86798,
                        17.10213, 29.57455, 21.57054, 22.63182, 31.91486, 22.33300,
                        28.86516, 36.76074, 24.55732, 35.51773, 43.11028, 27.96888,
                        39.00795, 51.44632, 34.46835, 41.99960, 58.42999, 38.73999,
                        43.45352, 64.43311, 41.88359, 40.34529, 69.27940, 45.46635,
                        35.54768, 71.87370, 48.47891, 31.55745, 71.16268, 48.58312,
                        27.21643, 69.03709, 48.64404, 23.64803, 65.86255, 47.40086,
                        21.21485, 61.57600, 44.61141, 19.69551, 58.02729, 41.41692,
                        18.02174, 55.05570, 39.20299, 17.04518, 52.47211, 37.37824,
                        16.04768, 49.92073, 36.64499, 15.59934, 47.62134, 35.34787,
                        15.21276, 45.79928, 34.13359, 15.19172, 44.10440, 32.33322,
                        15.67008, 43.92173, 31.26167, 16.23841, 44.74078, 31.07301,
                        17.23850, 46.00055, 31.04016, 18.15915, 47.46932, 31.50781,
                        18.55422, 48.64313, 32.84055, 19.04990, 49.51172, 33.63599,
                        18.23541, 49.53954, 34.37894, 17.18566, 48.82853, 34.61992,
                        16.20203, 47.63599, 34.32545, 15.46085, 46.37063, 33.90898,
                        14.98325, 44.72809, 32.88829, 14.62240, 42.85141, 31.61026,
                        14.52832, 40.63486, 29.82018, 14.38454, 38.47667, 28.65332,
                        14.21687, 36.64579, 27.44720, 13.89212, 35.02843, 26.50242,
                        14.67504, 34.32273, 25.36179, 16.87500, 35.18185, 24.53042,
                        21.27853, 38.68389, 25.00855, 27.49620, 45.67836, 27.32340)

  rec_expected_vec <- c(9.3401233, 5.2365995, 2.8938153, 9.4875858, 5.2782822,
                        2.9062696, 9.8085016, 5.3536376, 2.9297097, 10.4339567,
                        5.4726055, 2.9680079, 11.2766678, 5.5972561, 3.0086127,
                        12.0638388, 5.6626820, 3.0341238, 12.6332921, 5.6698740,
                        3.0435884, 12.5435441, 5.5626964, 3.0210939, 12.0228537,
                        5.3831786, 2.9774062, 11.2535941, 5.1411630, 2.9139784,
                        10.4208724, 4.8513431, 2.8343188, 9.5863108, 4.5107587,
                        2.7345369, 8.7293291, 4.1112480, 2.6108274, 7.9631980,
                        3.6769710, 2.4745198, 7.5118961, 3.2328321, 2.3324289,
                        2.5434335, 0.7204771, 0.7340939, 2.9171286, 0.5967620,
                        0.7172988, 3.5995822, 0.4876527, 0.7822116, 23.7685613,
                        0.6472714, 1.0402280, 9.0148663, 0.9133252, 13.9537818,
                        68.8274704, 1.4276178, 2.1858209, 2.8386438, 0.9870487,
                        0.6634661, 23.7476745, 1.1813142, 0.6749995, 80.4819369,
                        1.3348347, 0.9365452, 28.1107967, 2.0704505, 1.6396774,
                        2.5421489, 2.1512058, 0.9544748, 17.1089006, 18.5478090,
                        0.9782643, 5.3282683, 3.4394905, 0.8534229, 0.9832811,
                        3.5062431, 0.9560827, 0.7518336, 8.8211653, 1.9594734,
                        0.6232200, 1.9706499, 19.4931837, 1.1622252, 3.2036475,
                        5.6335716, 2.3238857, 6.6977833, 11.5034086, 0.5368965,
                        2.5023654, 1.6635521, 0.5425474, 9.4573376, 5.9093540,
                        0.7981424, 3.4461510, 1.3785374, 2.9310888, 3.2330424,
                        5.5195879, 8.0702599, 3.4884842, 2.5838121, 2.5051120,
                        8.7923321, 1.2517643, 9.7376861, 13.4100946, 0.8593748,
                        2.1452790, 28.4918022, 1.1357561, 2.2142446, 2.0956307,
                        1.6846894, 16.5174363, 15.8895060, 1.4800144, 4.5487519,
                        7.5476056, 1.0219531, 3.4365326, 3.1217332, 2.0820107,
                        3.3305282, 5.8340355, 2.1513339, 3.0520242, 3.8348572,
                        1.5280776, 2.5149269, 7.1834221, 1.6858094, 1.5839029,
                        2.7144396, 1.7588386, 3.3819137, 2.4957400, 8.9835182,
                        8.0440909, 4.9361704, 1.2238024, 4.4465728, 1.0800822,
                        0.7564184, 4.4697538, 1.0315829, 0.9663743, 0.9249992,
                        0.8817285, 1.2734311, 2.1862950, 1.4514425, 3.4053964,
                        6.1020950, 4.1068450, 0.9936529, 30.6883700, 23.6221809,
                        1.5841907, 3.0493127, 2.6659073, 1.3877151, 70.3381940,
                        48.1133686, 7.2479720, 8.8301509, 19.9198965, 0.7279112,
                        32.2958856, 14.3170497, 0.5901027, 10.0831337, 5.6585413,
                        3.1342959)

  expect_equal(as.vector(sabie_rtmb_model$rep$SSB), ssb_expected_vec, tolerance = 1e-3)
  expect_equal(as.vector(sabie_rtmb_model$rep$Rec), rec_expected_vec, tolerance = 1e-3)
  expect_true(sd_rep$pdHess)

})

