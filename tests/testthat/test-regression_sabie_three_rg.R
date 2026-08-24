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

  ssb_expected_vec <- c(39.45039, 105.58367, 88.06524, 39.15294, 104.70083, 87.11024, 37.74715, 101.30070, 83.43663, 35.58635, 95.88802, 77.81483, 34.22179, 92.34815, 
                        75.34072, 34.16248, 91.08539, 74.63502, 34.01916, 89.95182, 74.16311, 33.35109, 87.90105, 72.74657, 32.50180, 85.52968, 
                        71.05176, 30.68292, 81.23637, 67.48333, 28.47847, 76.05178, 62.98735, 26.53916, 70.82087, 57.65840, 23.78427, 64.42127, 
                        52.11630, 20.04704, 56.14345, 45.14875, 18.13821, 50.75165, 40.14011, 16.46142, 45.68154, 35.21545, 15.06011, 41.30453, 
                        31.14248, 13.23730, 36.34870, 26.77096, 11.90658, 33.93717, 23.46828, 11.78265, 32.17355, 23.06201, 12.54232, 30.39744, 
                        21.93937, 14.01218, 29.57180, 21.77361, 17.10050, 29.60055, 21.51610, 22.61826, 31.93953, 22.30104, 28.85383, 36.77549, 
                        24.53432, 35.50408, 43.11533, 27.95830, 39.00096, 51.43466, 34.46675, 41.99605, 58.40810, 38.74761, 43.45277, 64.40441, 
                        41.89859, 40.34762, 69.24141, 45.48978, 35.55099, 71.83075, 48.50709, 31.55656, 71.12061, 48.61663, 27.21829, 68.99657, 
                        48.67474, 23.64988, 65.82806, 47.42828, 21.21416, 61.54750, 44.63672, 19.68991, 58.00855, 41.43636, 18.01645, 55.04072, 
                        39.21664, 17.03673, 52.46024, 37.39037, 16.03768, 49.91549, 36.65051, 15.58583, 47.62157, 35.35051, 15.20012, 45.79962, 
                        34.13441, 15.17840, 44.10760, 32.33201, 15.65657, 43.92526, 31.25996, 16.22511, 44.74226, 31.07240, 17.22515, 46.00053, 
                        31.03948, 18.14449, 47.46593, 31.50852, 18.54143, 48.63375, 32.84240, 19.03644, 49.49964, 33.63871, 18.22429, 49.52279, 
                        34.38328, 17.17580, 48.81072, 34.62448, 16.19270, 47.61738, 34.33148, 15.45263, 46.35277, 33.91467, 14.97477, 44.71143, 
                        32.89403, 14.61482, 42.83573, 31.61466, 14.51973, 40.62110, 29.82368, 14.37502, 38.46662, 28.65439, 14.20784, 36.63641, 
                        27.44775, 13.88304, 35.02343, 26.50147, 14.66592, 34.32110, 25.36007, 16.86767, 35.18264, 24.52897, 21.27067, 38.68684, 
                        25.00732, 27.49207, 45.68322, 27.32420)

  rec_expected_vec <- c(9.8009459, 5.3995180, 2.9523811, 9.9477756, 5.4410028, 2.9644338, 10.2864297, 5.5189026, 2.9881350, 10.9622897, 5.6437432, 3.0275782, 11.8702186, 
                        5.7735510, 3.0691051, 12.7007979, 5.8379464, 3.0945241, 13.2612730, 5.8377000, 3.1025230, 13.0807041, 5.7163132, 3.0772798, 
                        12.4491362, 5.5197108, 3.0301706, 11.5826468, 5.2595233, 2.9629082, 10.6777433, 4.9517796, 2.8792469, 9.7923352, 4.5940616, 
                        2.7752871, 8.8973451, 4.1778603, 2.6471622, 8.1037293, 3.7278405, 2.5064505, 7.6301115, 3.2691829, 2.3597779, 2.5227798, 
                        0.7220580, 0.7388365, 2.8910756, 0.5981280, 0.7243600, 3.5598444, 0.4887611, 0.7967551, 23.5809696, 0.6506278, 1.0655011, 
                        9.1719921, 0.9195914, 13.9472601, 68.6460889, 1.4390611, 2.2092983, 2.9020993, 0.9964304, 0.6681210, 23.7259944, 1.1939211, 
                        0.6796377, 80.4256527, 1.3513060, 0.9438988, 28.0896303, 2.1004243, 1.6505632, 2.5893199, 2.1838684, 0.9603608, 17.0452329, 
                        18.4498883, 0.9855226, 5.3415790, 3.4894005, 0.8583768, 0.9933192, 3.5381173, 0.9605651, 0.7597569, 8.8093005, 1.9678562, 
                        0.6301075, 2.0035752, 19.3596659, 1.1713750, 3.2314690, 5.7002413, 2.3412697, 6.7233612, 11.3686270, 0.5425456, 2.5452053, 
                        1.6904434, 0.5482301, 9.4214940, 5.8719426, 0.8062463, 3.4644969, 1.3931283, 2.9476934, 3.2352652, 5.4800950, 8.0473468, 
                        3.5304974, 2.5830664, 2.5370347, 8.7282973, 1.2538072, 9.7133164, 13.4810118, 0.8635928, 2.1678619, 28.3890199, 1.1446465, 
                        2.2296968, 2.1261406, 1.6927462, 16.4915465, 15.8582026, 1.4879725, 4.5750371, 7.5320081, 1.0276415, 3.4315777, 3.1479795, 
                        2.0785337, 3.3333319, 5.8089862, 2.1483859, 3.0546182, 3.8622509, 1.5283494, 2.5164738, 7.1421483, 1.6911606, 1.5913368, 
                        2.7294897, 1.7729419, 3.3952400, 2.5238439, 8.9197334, 8.0428887, 4.9165812, 1.2328484, 4.4536185, 1.0865382, 0.7590346, 
                        4.4551934, 1.0357629, 0.9676603, 0.9326038, 0.8866270, 1.2760992, 2.1966715, 1.4602110, 3.3789738, 6.1016728, 4.1250332, 
                        0.9981187, 30.6640743, 23.5899483, 1.5890103, 3.0965276, 2.7011504, 1.3967008, 70.2402102, 48.0955688, 7.2422607, 9.0033030, 
                        19.9189351, 0.7340583, 32.0991798, 14.3032566, 0.5948105, 10.4622083, 5.8004354, 3.1876287)

  expect_equal(as.vector(sabie_rtmb_model$rep$SSB), ssb_expected_vec, tolerance = 1e-3)
  expect_equal(as.vector(sabie_rtmb_model$rep$Rec), rec_expected_vec, tolerance = 1e-3)
  expect_true(sd_rep$pdHess)
  expect_jnLL_decomposes(sabie_rtmb_model)

})

