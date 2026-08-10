# Pinned regression test. The expected SSB and recruitment vectors are output from a
# previously validated SPoRC fit of this assessment, not hand-derived values. A mismatch
# means a change moved a fitted result, which is a bug unless the numerical change was
# intended. If it was intended, re-baseline deliberately and say why in NEWS.md. Do not
# paste in fresh output to make the test pass. See tests/README.md.
#
# Expectations here are flattened with as.vector over a multi-region array, so a change to
# region or season array ordering shows up as a whole-vector mismatch rather than a shift.

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
  map_fish_fixed <- array(mapping$fish_fixed_sel_pars, dim = dim(parameters$fish_fixed_sel_pars))
  map_fish_fixed[,2,1,2,1]  <- map_fish_fixed[,2,1,1,1] # share deltas

  # Map off bmax for trawl females
  map_fish_fixed[,1,1,2,2]  <- map_fish_fixed[,1,1,1,2] # share deltas
  mapping$fish_fixed_sel_pars <- factor(map_fish_fixed)

  # Map off delta for JP LLS
  map_srv_fixed <- array(mapping$srv_fixed_sel_pars, dim = dim(parameters$srv_fixed_sel_pars))
  map_srv_fixed[,2,1,2,2]  <- map_srv_fixed[,2,1,1,2] # share deltas
  mapping$srv_fixed_sel_pars <- factor(map_srv_fixed)

  # Some starting values to help out the model
  parameters$srv_fixed_sel_pars[] <- log(2)
  parameters$fish_fixed_sel_pars[,,,,1] <- log(2) # fixed gear
  parameters$fish_fixed_sel_pars[,,,,2] <- log(5) # trawl gear

  # fit model
  suppressWarnings(
    sabie_rtmb_model <- fit_model(data,
                                  parameters,
                                  mapping,
                                  random = NULL,
                                  newton_loops = 3,
                                  silent = TRUE, do_optim = TRUE
    )
  )

  sd_rep <- sdreport(sabie_rtmb_model)
  rep <- sabie_rtmb_model$rep

  ssb_expected_vec <- c(40.42253, 107.86754, 89.73232, 40.11332, 107.08797, 88.89238, 38.65002, 103.72859, 85.30239, 36.38329, 98.27178, 79.71293, 34.88954, 94.61205, 77.21587, 34.71510, 93.17660, 76.42893, 34.46119, 91.82745, 75.83077, 33.68868, 89.53401, 74.25652, 32.74652, 86.91366, 72.39022, 30.84356, 82.37057, 68.64034, 28.57215, 76.95310, 63.96549, 26.58639, 71.51852, 58.46735, 23.79675, 64.94169, 52.77060, 20.03567, 56.51601, 45.66102, 18.11770, 51.01851, 40.53576, 16.43583, 45.86687, 35.51319, 15.03203, 41.42770, 31.36146, 13.20838, 36.42353, 26.92537, 11.88073, 33.97964, 23.57372, 11.76253, 32.19280, 23.13367, 12.53046, 30.40141, 21.98366, 14.00740, 29.56763, 21.80031, 17.10111, 29.59350, 21.53155, 22.62237, 31.93288, 22.31019, 28.85741, 36.77169, 24.54100, 35.50839, 43.11434, 27.96153, 39.00339, 51.43856, 34.46747, 41.99752, 58.41498, 38.74571, 43.45338, 64.41328, 41.89463, 40.34726, 69.25305, 45.48342, 35.55027, 71.84392, 48.49936, 31.55700, 71.13360, 48.60737, 27.21791, 69.00918, 48.66629, 23.64950, 65.83895, 47.42081, 21.21451, 61.55665, 44.62991, 19.69168, 58.01484, 41.43131, 18.01814, 55.04590, 39.21329, 17.03933, 52.46450, 37.38748, 16.04073, 49.91778, 36.64958, 15.58989, 47.62224, 35.35043, 15.20392, 45.80022, 34.13485, 15.18241, 44.10735, 32.33304, 15.66064, 43.92489, 31.26112, 16.22913, 44.74249, 31.07323, 17.22919, 46.00120, 31.04031, 18.14893, 47.46761, 31.50894, 18.54534, 48.63720, 32.84249, 19.04054, 49.50389, 33.63856, 18.22771, 49.52842, 34.38266, 17.17885, 48.81669, 34.62378, 16.19558, 47.62358, 34.33035, 15.45518, 46.35875, 33.91363, 14.97738, 44.71704, 32.89297, 14.61717, 42.84104, 31.61398, 14.52238, 40.62582, 29.82326, 14.37794, 38.47023, 28.65467, 14.21062, 36.63980, 27.44818, 13.88583, 35.02551, 26.50234, 14.66872,
                        34.32218, 25.36115, 16.86998, 35.18300, 24.52997, 21.27319, 38.68657, 25.00824, 27.49360, 45.68241, 27.32455)

  rec_expected_vec <- c(9.6668703, 5.3522270, 2.9355184, 9.8144428, 5.3939010, 2.9477237, 10.1483043, 5.4711338, 2.9713681, 10.8097942, 5.5942765, 3.0104775, 11.6991519, 5.7225873, 3.0517351, 12.5177047, 5.7872946, 3.0771719, 13.0818220, 5.7892758, 3.0855988, 12.9284615, 5.6720905, 3.0611576, 12.3291398, 5.4805084, 3.0150453, 11.4903654, 5.2256246, 2.9488971, 10.6057857, 4.9230777, 2.8663951, 9.7345937, 4.5702954, 2.7636419, 8.8501908, 4.1588812, 2.6367923, 8.0642488, 3.7133642, 2.4973536, 7.5969545, 3.2588516, 2.3520017, 2.5286426, 0.7216194, 0.7375122, 2.8985117, 0.5977473, 0.7223800, 3.5711508, 0.4884500, 0.7926346, 23.6344085, 0.6496813, 1.0582838, 9.1278768, 0.9178031, 13.9497242, 68.6973734, 1.4358499, 2.2026994, 2.8839975, 0.9937650, 0.6668003, 23.7332931, 1.1903495, 0.6783286, 80.4412417, 1.3466377, 0.9418358, 28.0951404, 2.0919651, 1.6475833, 2.5759180, 2.1745673, 0.9587026, 17.0630485, 18.4785953, 0.9834705, 5.3378095, 3.4751614, 0.8569754, 0.9904635, 3.5290780, 0.9593079, 0.7574997, 8.8124826, 1.9656362, 0.6281423, 1.9941372, 19.3983713, 1.1687994, 3.2235826, 5.6809318, 2.3362379, 6.7160409, 11.4076279, 0.5409367, 2.5329440, 1.6827128, 0.5466153, 9.4315802, 5.8828427, 0.8039430, 3.4594225, 1.3889366, 2.9429608, 3.2347831, 5.4912992, 8.0539027, 3.5184332, 2.5832772, 2.5279605, 8.7468393, 1.2532405, 9.7201500, 13.4607263, 0.8624061, 2.1614845, 28.4188375, 1.1420890, 2.2253652, 2.1174356, 1.6904062, 16.4989024, 15.8674898, 1.4856956, 4.5673874, 7.5366375, 1.0260496, 3.4330929, 3.1405466, 2.0794687, 3.3325648, 5.8161917, 2.1491843, 3.0538589, 3.8545847, 1.5282685, 2.5160385, 7.1539882, 1.6896009, 1.5892608, 2.7252600, 1.7688935, 3.3914354,
                        2.5158610, 8.9379266, 8.0431937, 4.9224860, 1.2302717, 4.4516790, 1.0847038, 0.7582844, 4.4592769, 1.0345916, 0.9672864, 0.9304722, 0.8852402, 1.2753402, 2.1936735, 1.4577490, 3.3864268, 6.1019322, 4.1199261, 0.9968865, 30.6710938, 23.5993809, 1.5875769, 3.0831886, 2.6911448, 1.3940958, 70.2675880, 48.1015370, 7.2440726, 8.9549164, 19.9197053,
                        0.7322962, 32.1551112, 14.3071467, 0.5934652, 10.3531844, 5.7595411, 3.1723423)

  expect_equal(as.vector(sabie_rtmb_model$rep$SSB), ssb_expected_vec, tolerance = 1e-3)
  expect_equal(as.vector(sabie_rtmb_model$rep$Rec), rec_expected_vec, tolerance = 1e-3)
  expect_true(sd_rep$pdHess)
  expect_jnLL_decomposes(sabie_rtmb_model)

})

