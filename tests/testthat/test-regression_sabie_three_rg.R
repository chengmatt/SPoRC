# Regression test. Expected SSB and recruitment come from a previously validated SPoRC fit, not from
# hand-derived values, so a mismatch means something moved a fitted result. See tests/README.md.
# Regenerated 2026-10-12 after the fishery input sample size fix (see NEWS.md).
#
# Expectations here are flattened with as.vector over a multi-region array, so a change to
# region or season array ordering shows up as a whole-vector mismatch rather than a shift.

library(SPoRC)
library(testthat)
data("three_rg_sable_data")

test_that("Three-region Sablefish RTMB model produces expected results", {


  # Initialize model dimensions and data list
  input_list <- Setup_Mod_Dim(years = seq_along(three_rg_sable_data$years),
                              # years 1 to 62
                              ages = seq_along(three_rg_sable_data$ages),
                              # ages 1 to 30
                              lens = three_rg_sable_data$lens,
                              # lengths 41 to 99
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
                                  # tag releases: region in the first column, year in the second,
                                  # one row per tagged cohort
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
                                  conv_tag_fish_reporting_pars = array(log(0.2 / (1 - 0.2)), dim = c(input_list$data$n_regions, 2, input_list$data$n_fish_fleets))
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
                                                                   input_list$data$n_fish_fleets))
                                        # fixing catch sd at small value
    )
  )

  # ln_F_mean cannot be passed through ... because R partially matches the name
  # to the ln_F_mean_spec formal, so the starting value is assigned post-hoc.
  input_list$par$ln_F_mean[] <- -2

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
    dplyr::filter(!(fleet == 1 & block == 1 & sex == 2 & par == 2)) %>%              # drop priors on parameters that are not estimated, par 1 a50 and par 2 delta, which has to agree with the mapping
    dplyr::filter(!(fleet == 2 & block == 1 & sex == 2 & par == 1))                  # same for the second fleet

  # lognormal prior values, one row per parameter combination the prior applies to
  fish_selex_prior <- cbind(
    region = 1,
    fish_selex_structure,
    mu = 2,                                                                      # All selex means = 1 (means should be defined in normal space)
    sd = 3                                                                       # All selex sd = 5
  )

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

  # fleet and block combinations, estimating only the domestic and Japanese longline surveys
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

                                       # estimate every survey selectivity parameter here,
                                       # then fix or share some below
                                       srv_fixed_sel_pars_spec =
                                         c("est_shared_r",
                                           "est_shared_r"
                                         ),

                                       # estimate every survey catchability,
                                       # one q for every region
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

  # fishery comps: five region values summed over BS, AI, and WGOA for region 1

  # Fishery Ages
  data$ISS_FishAgeComps[1,,,,] <- 59  # BS_AI_WGOA
  data$ISS_FishAgeComps[2,,,,] <- 18  # CGOA
  data$ISS_FishAgeComps[3,,,,] <- 18  # EGOA

  # fishery lengths, fixed gear
  data$ISS_FishLenComps[1,,,,1] <- 31  # BS_AI_WGOA
  data$ISS_FishLenComps[2,,,,1] <-  7  # CGOA
  data$ISS_FishLenComps[3,,,,1] <-  7  # EGOA

  # fishery lengths, trawl gear
  data$ISS_FishLenComps[1,,,,2] <- 36  # BS_AI_WGOA
  data$ISS_FishLenComps[2,,,,2] <-  4  # CGOA
  data$ISS_FishLenComps[3,,,,2] <-  4  # EGOA

  # Map off early delta for fishery
  map_fish_fixed <- array(mapping$fish_fixed_sel_pars, dim = dim(parameters$fish_fixed_sel_pars))
  map_fish_fixed[,2,1,2,1]  <- map_fish_fixed[,2,1,1,1] # share deltas

  # Map off bmax for trawl females
  map_fish_fixed[,1,1,2,2]  <- map_fish_fixed[,1,1,1,2] # share deltas
  mapping$fish_fixed_sel_pars <- factor(map_fish_fixed)

  # map off delta for the Japanese longline survey
  map_srv_fixed <- array(mapping$srv_fixed_sel_pars, dim = dim(parameters$srv_fixed_sel_pars))
  map_srv_fixed[,2,1,2,2]  <- map_srv_fixed[,2,1,1,2] # share deltas
  mapping$srv_fixed_sel_pars <- factor(map_srv_fixed)

  # Some starting values to help out the model
  parameters$srv_fixed_sel_pars[] <- log(2)
  parameters$fish_fixed_sel_pars[,,,,1] <- log(2) # fixed gear
  parameters$fish_fixed_sel_pars[,,,,2] <- log(5) # trawl gear

  # fit model
  suppressWarnings(
    sabie_rtmb_model <- fit_model(
      data,
      parameters,
      mapping,
      random = NULL,
      newton_loops = 3,
      silent = TRUE,
      do_optim = TRUE
    )
  )

  sd_rep <- sdreport(sabie_rtmb_model)
  rep <- sabie_rtmb_model$rep

  ssb_expected_vec <- c(40.95724, 98.03610, 84.51836, 40.64604, 97.11775, 83.47729, 39.20403, 93.78996, 79.76416, 36.99184, 88.57656, 74.15983,
                        35.87366, 85.59352, 71.88787, 35.93806, 84.73012, 71.38898, 36.06477, 84.20178, 71.26249, 35.74253, 82.94798, 70.34777,
                        35.25866, 81.47695, 69.27438, 33.80308, 78.21438, 66.43576, 31.83088, 73.98812, 62.67529, 29.77157, 69.27197, 57.86859,
                        27.04024, 63.58134, 52.86539, 23.28716, 56.05886, 46.41358, 21.00796, 50.82329, 41.63096, 18.89556, 45.68779, 36.71524,
                        17.10166, 41.15956, 32.52981, 14.97115, 36.06680, 27.97108, 13.35635, 33.41974, 24.42802, 13.07590, 31.47857, 23.84725,
                        13.63527, 29.58739, 22.58383, 14.95010, 28.67243, 22.29134, 17.93043, 28.59964, 21.87604, 23.31651, 30.72721, 22.43576,
                        29.53051, 35.15532, 24.33411, 36.09346, 41.01446, 27.58933, 39.50109, 48.71695, 33.90535, 42.51310, 55.12812, 37.98593,
                        43.70079, 60.64232, 41.06156, 40.49966, 65.04619, 44.57443, 35.70849, 67.40359, 47.51256, 31.73026, 66.71550, 47.62125,
                        27.37144, 64.66253, 47.82301, 23.87825, 61.72065, 46.64814, 21.56566, 57.70190, 43.84982, 20.19453, 54.51269, 40.61553,
                        18.60846, 51.84058, 38.40407, 17.68176, 49.46760, 36.63938, 16.77653, 47.13814, 35.87433, 16.39668, 45.00259, 34.55770,
                        16.05722, 43.41568, 33.37656, 16.04472, 41.94771, 31.64416, 16.49790, 41.81436, 30.60595, 17.03990, 42.64145, 30.38453,
                        17.83203, 43.73007, 30.45829, 18.77652, 45.08131, 30.68590, 19.02990, 46.09379, 32.01459, 19.39253, 46.84680, 32.78045,
                        18.45892, 46.80173, 33.52244, 17.31131, 46.08276, 33.78846, 16.30028, 44.98838, 33.51481, 15.57049, 43.88419, 33.17301,
                        15.17109, 42.50237, 32.24408, 14.92709, 40.86461, 31.11098, 14.92905, 38.90784, 29.39680, 14.86433, 36.89429, 28.31336,
                        14.59887, 35.14822, 27.25831, 14.33526, 33.63719, 26.31314, 15.23232, 33.16236, 25.15670, 17.63095, 34.21971, 24.34188,
                        22.21375, 37.60169, 24.79589, 28.58722, 44.18710, 27.09266)

  rec_expected_vec <- c(8.9766835, 5.8090897, 3.3101997, 9.1932561, 5.8983850, 3.3404680, 9.5914622, 6.0365222, 3.3886532, 10.3230657, 6.2422012,
                        3.4630855, 11.2877339, 6.4745067, 3.5385225, 12.3633346, 6.6548656, 3.6089044, 13.1549111, 6.7242936, 3.6450243, 13.1638708,
                        6.6214685, 3.6304436, 12.6056054, 6.4011048, 3.5799490, 11.7328379, 6.0836015, 3.4951052, 10.8119536, 5.7017993, 3.3858769,
                        9.9234018, 5.2669669, 3.2526689, 9.0069478, 4.7520857, 3.0776788, 8.1273993, 4.1706412, 2.8747302, 7.6070791, 3.6059318,
                        2.6709716, 2.7054410, 0.7890830, 0.8200071, 3.0978408, 0.6426694, 0.7807384, 3.7810492, 0.5153928, 0.8223733, 21.1069687,
                        0.7119609, 1.1422600, 12.1297707, 0.9810727, 13.6741287, 62.4996735, 1.8486971, 2.9295814, 2.6587082, 1.0911461, 0.6980897,
                        26.4603994, 1.3749574, 0.7163761, 73.8137098, 1.5865895, 1.0808758, 25.2218613, 3.0504072, 2.2356726, 2.4475013, 2.4503309,
                        1.0992945, 14.9318327, 22.8342030, 1.1360396, 4.0881558, 3.4901320, 0.9337981, 0.8497069, 3.2465683, 1.0457576, 0.6135906,
                        9.5892949, 3.2482393, 0.4960318, 2.3504455, 17.6497554, 0.7336475, 4.2312610, 6.2003816, 1.4382581, 8.7125682, 11.1928584,
                        0.4641959, 2.4017756, 1.6284901, 0.4717690, 8.2223085, 6.4013664, 0.6517815, 5.5040151, 1.2721413, 2.0304442, 3.6049768,
                        5.5735027, 7.4171297, 2.5365286, 3.1488372, 2.1110192, 13.0228943, 1.8238975, 9.2232889, 5.4762323, 1.2067659, 1.7226012,
                        32.8931786, 1.0758114, 1.4603168, 2.8421043, 1.4313116, 14.8229402, 16.0203500, 1.7520944, 3.7048941, 8.3219695, 1.4053290,
                        2.1131116, 3.1904750, 2.4990864, 2.1046756, 7.8366677, 1.8910518, 2.3089779, 4.4429812, 1.4177287, 1.5447270, 9.1789443,
                        1.1644131, 0.9785644, 3.8123131, 1.5241543, 2.4057953, 3.1978284, 10.2835307, 7.7062618, 6.8652061, 0.9874040, 3.0593319,
                        1.0567570, 0.7060449, 3.5889044, 1.1077787, 0.9853811, 0.6934551, 0.9268406, 1.6413394, 1.9136084, 1.6149411, 3.9335149,
                        7.0399933, 5.2150249, 1.1766855, 26.4502331, 24.5888028, 1.2195140, 4.5807652, 3.0198258, 1.2322133, 68.2231486, 53.1271302,
                        7.4360789, 11.2692665, 9.7986471, 0.8310624, 23.8018539, 3.2872105, 0.7403774, 9.6209959, 6.2381043, 3.5722104)

  expect_equal(as.vector(sabie_rtmb_model$rep$SSB), ssb_expected_vec, tolerance = 1e-3)
  expect_equal(as.vector(sabie_rtmb_model$rep$Rec), rec_expected_vec, tolerance = 1e-3)
  expect_true(sd_rep$pdHess)
  expect_jnLL_decomposes(sabie_rtmb_model)

})
