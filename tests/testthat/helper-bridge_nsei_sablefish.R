# Purpose: To set up the 2025 NSEI sablefish TMB model (scaa_mod_v23_Aaron2) in SPoRC, with the fishery cpue and mark-recapture as fishery indices
# Creator: Matthew LH. Cheng
# Date: 11/14/26

# Setup Model -------------------------------------------------------------

build_nsei_sablefish_input <- function(tmb_data, tmb_par) {

  # tmb sex index 1 is male and 2 is female, the other way round from SPoRC
  years <- 1975:2025
  ages <- 2:31
  n_yrs <- length(years)
  n_ages <- length(ages)
  n_lens <- length(tmb_data$lenbin)
  n_sexes <- 2

  # Initialize model dimensions and data list
  input_list <- Setup_Mod_Dim(years = years, # vector of years (1975 - 2025)
                              ages = ages, # vector of ages
                              lens = tmb_data$lenbin, # vector of lengths
                              n_regions = 1, # number of regions
                              n_sexes = n_sexes, # number of sexes == 1, female, == 2 male
                              n_fish_fleets = 2, # number of fishery fleets == 1, fixed gear fishery, == 2 mark-recapture (no catch)
                              n_srv_fleets = 1, # number of survey fleets == 1, longline survey
                              bias_correct_pe = "none", # no bias correction in the tmb model
                              verbose = FALSE
  )

  # Setup recruitment stuff
  # the tmb rec penalty is sum((dev - sigma^2 / 2)^2) / sigma^2, a normal at sigmaR / sqrt(2)
  input_list <- Setup_Mod_Rec(input_list = input_list, # input data list from above
                              rec_model = "mean_rec", # recruitment model
                              sigmaR_spec = "fix", # fix sigmaR
                              ln_sigmaR = array(log(exp(tmb_par$log_sigma_r) / sqrt(2)), dim = c(2, 1, 1)), # sigmaR of 1.2 in the tmb model
                              sigmaR_switch = 1, # one sigmaR for all years
                              do_rec_bias_ramp = 0, # no bias ramp
                              init_age_strc = 1, # geometric series to derive initial age structure
                              equil_init_age_strc = "stoch_no_plus", # initial deviations on every age but the plus group
                              init_F_form = "abs", # initial F as a value
                              init_F_spec = "fix", # fix initial F
                              init_F_par = array(log(1e-10), dim = c(1, 1, 2)), # no fishing in the initial age structure
                              use_rinit = 1, # separate recruitment for the initial age structure (log_rinit)
                              t_spawn = tmb_data$spawn_month # spawn timing
  )

  # Setup biological stuff
  # survey weight at age is used for SSB and the survey index, fishery weight at age for the landed catch and fishery cpue
  WAA <- array(0, dim = c(1, 1, n_yrs, 1, n_ages, n_sexes)) # weight at age
  WAA_fish <- array(0, dim = c(1, 1, n_yrs, 1, n_ages, n_sexes, 2)) # fishery weight at age
  WAA_srv <- array(0, dim = c(1, 1, n_yrs, 1, n_ages, n_sexes, 1)) # survey weight at age
  MatAA <- array(0, dim = c(1, 1, n_yrs, 1, n_ages, n_sexes)) # maturity at age
  Fixed_natmort <- array(exp(tmb_par$log_M), dim = c(1, 1, n_yrs, n_ages, n_sexes)) # natural mortality of 0.1
  SizeAgeTrans <- array(0, dim = c(1, 1, n_yrs, 1, n_lens, n_ages, n_sexes)) # size age transition

  for(y in 1:n_yrs) {
    for(s in 1:n_sexes) {
      tmb_sex <- 3 - s # tmb sex index
      WAA[1,1,y,1,,s] <- tmb_data$data_srv_waa[1,,tmb_sex]
      WAA_fish[1,1,y,1,,s,] <- tmb_data$data_fsh_waa[1,,tmb_sex] # same for both fishery fleets, mark-recapture is in numbers
      WAA_srv[1,1,y,1,,s,1] <- tmb_data$data_srv_waa[1,,tmb_sex]
      MatAA[1,1,y,1,,s] <- tmb_data$prop_mature[1,]
      SizeAgeTrans[1,1,y,1,,,s] <- t(tmb_data$agelen_key_fsh[,,tmb_sex]) # tmb key is age by length
    } # end s loop
  } # end y loop

  input_list <- Setup_Mod_Biologicals(input_list = input_list,
                                      # Data inputs
                                      WAA = WAA, # weight at age
                                      WAA_fish = WAA_fish, # fishery weight at age
                                      WAA_srv = WAA_srv, # survey weight at age
                                      MatAA = MatAA, # maturity at age
                                      AgeingError = tmb_data$ageing_error, # ageing error
                                      SizeAgeTrans = SizeAgeTrans, # size age transition matrix
                                      # Model options
                                      fit_lengths = 1, # fitting length compositions
                                      M_spec = "fix", # fix natural mortality
                                      Fixed_natmort = Fixed_natmort
  )

  # Setup movement and tagging stuff (none)
  input_list <- Setup_Mod_Movement(input_list = input_list,
                                   use_fixed_movement = 1,
                                   Fixed_Movement = NA,
                                   do_recruits_move = 0
  )

  input_list <- Setup_Mod_Tagging(input_list = input_list, use_conv_fish_tagging = 0)

  # Setup catch and fishing mortality stuff
  # the mark-recapture fleet has a recorded catch of zero, so its F is zero in every year
  ObsCatch <- array(0, dim = c(1, n_yrs, 1, 2))
  UseCatch <- array(0, dim = c(1, n_yrs, 1, 2))
  ObsCatch[1,,1,1] <- tmb_data$data_catch * 1e3 # landed catch in kg, weight at age is in kg
  UseCatch[1,,1,1] <- 1

  # the tmb F penalty is 0.1 * sum(dev^2), a normal at sd sqrt(5)
  input_list <- Setup_Mod_Catch_and_F(input_list = input_list,
                                      # Data inputs
                                      ObsCatch = ObsCatch,
                                      UseCatch = UseCatch,
                                      # Model options
                                      catch_units = c("biom", "biom"), # catch in biomass
                                      sigmaC_spec = "fix", # fix catch sd
                                      ln_sigmaC = array(log(tmb_data$sigma_catch), dim = c(1, n_yrs, 1, 2)), # catch sd of 0.05
                                      Use_F_pen = 1, # use F penalty
                                      ln_sigmaF = array(log(sqrt(1 / (2 * tmb_data$wt_fpen))), dim = c(1, 1, 2)), # F penalty sd
                                      Fdev_model = "iid", # F deviations are iid
                                      dmr_mean_spec = "fix" # fix discard mortality
  )

  input_list$par$logit_dmr_mean[] <- qlogis(0.16) # discard mortality of 0.16

  # Setup fishery indices and compositions
  # fishery cpue in biomass on fleet 1, mark-recapture abundance (millions of fish) on fleet 2
  ObsFishIdx <- array(NA, dim = c(1, n_yrs, 1, 2))
  ObsFishIdx_SE <- array(NA, dim = c(1, n_yrs, 1, 2))
  UseFishIdx <- array(0, dim = c(1, n_yrs, 1, 2))

  fsh_cpue_yrs <- tmb_data$yrs_fsh_cpue + 1 # tmb years start at 0
  ObsFishIdx[1,fsh_cpue_yrs,1,1] <- tmb_data$data_fsh_cpue
  ObsFishIdx_SE[1,fsh_cpue_yrs,1,1] <- tmb_data$sigma_fsh_cpue
  UseFishIdx[1,fsh_cpue_yrs,1,1] <- 1

  mr_yrs <- tmb_data$yrs_mr + 1
  ObsFishIdx[1,mr_yrs,1,2] <- tmb_data$data_mr
  ObsFishIdx_SE[1,mr_yrs,1,2] <- tmb_data$sigma_mr
  UseFishIdx[1,mr_yrs,1,2] <- 1

  # ages are sexes combined, lengths are by sex, both from the landed catch
  ObsFishAgeComps <- array(0, dim = c(1, n_yrs, 1, n_ages, n_sexes, 2))
  ISS_FishAgeComps <- array(0, dim = c(1, n_yrs, 1, n_sexes, 2))
  UseFishAgeComps <- array(0, dim = c(1, n_yrs, 1, 2))
  ObsFishLenComps <- array(0, dim = c(1, n_yrs, 1, n_lens, n_sexes, 2))
  ISS_FishLenComps <- array(0, dim = c(1, n_yrs, 1, n_sexes, 2))
  UseFishLenComps <- array(0, dim = c(1, n_yrs, 1, 2))

  fsh_age_yrs <- tmb_data$yrs_fsh_age + 1
  ObsFishAgeComps[1,fsh_age_yrs,1,,1,1] <- tmb_data$data_fsh_age
  ISS_FishAgeComps[1,fsh_age_yrs,1,1,1] <- tmb_data$effn_fsh_age
  UseFishAgeComps[1,fsh_age_yrs,1,1] <- 1

  fsh_len_yrs <- tmb_data$yrs_fsh_len + 1
  for(s in 1:n_sexes) {
    ObsFishLenComps[1,fsh_len_yrs,1,,s,1] <- tmb_data$data_fsh_len[,,3 - s]
    ISS_FishLenComps[1,fsh_len_yrs,1,s,1] <- tmb_data$effn_fsh_len[,1,3 - s]
  } # end s loop
  UseFishLenComps[1,fsh_len_yrs,1,1] <- 1

  input_list <- Setup_Mod_FishIdx_and_Comps(input_list = input_list,
                                            # data inputs
                                            ObsFishIdx = ObsFishIdx, # observed fishery index
                                            ObsFishIdx_SE = ObsFishIdx_SE, # lognormal sd
                                            UseFishIdx = UseFishIdx, # whether fishery indices are used
                                            ObsFishAgeComps = ObsFishAgeComps, # observed fishery ages
                                            UseFishAgeComps = UseFishAgeComps, # whether fishery ages are used
                                            ISS_FishAgeComps = ISS_FishAgeComps, # input sample size for fishery ages
                                            ObsFishLenComps = ObsFishLenComps, # observed fishery lengths
                                            UseFishLenComps = UseFishLenComps, # whether fishery lengths are used
                                            ISS_FishLenComps = ISS_FishLenComps, # input sample size for fishery lengths
                                            # Model options
                                            fish_idx_type = c("biom", "abd"), # biomass for fishery cpue, abundance for mark-recapture
                                            FishAgeComps_LikeType = c("Multinomial", "none"), # age comp likelihoods for fishery fleets
                                            FishLenComps_LikeType = c("Multinomial", "none"), # length comp likelihoods for fishery fleets
                                            FishAgeComps_Type = c("agg_Year_1-terminal_Fleet_1",
                                                                  "none_Year_1-terminal_Fleet_2"), # sexes combined
                                            FishLenComps_Type = c("spltRspltS_Year_1-terminal_Fleet_1",
                                                                  "none_Year_1-terminal_Fleet_2") # by sex
  )

  # Setup survey indices and compositions (longline survey only)
  ObsSrvIdx <- array(NA, dim = c(1, n_yrs, 1, 1))
  ObsSrvIdx_SE <- array(NA, dim = c(1, n_yrs, 1, 1))
  UseSrvIdx <- array(0, dim = c(1, n_yrs, 1, 1))

  srv_cpue_yrs <- tmb_data$yrs_srv_cpue + 1
  ObsSrvIdx[1,srv_cpue_yrs,1,1] <- tmb_data$data_srv_cpue
  ObsSrvIdx_SE[1,srv_cpue_yrs,1,1] <- tmb_data$sigma_srv_cpue
  UseSrvIdx[1,srv_cpue_yrs,1,1] <- 1

  # survey ages are sexes combined, lengths are by sex
  ObsSrvAgeComps <- array(0, dim = c(1, n_yrs, 1, n_ages, n_sexes, 1))
  ISS_SrvAgeComps <- array(0, dim = c(1, n_yrs, 1, n_sexes, 1))
  UseSrvAgeComps <- array(0, dim = c(1, n_yrs, 1, 1))
  ObsSrvLenComps <- array(0, dim = c(1, n_yrs, 1, n_lens, n_sexes, 1))
  ISS_SrvLenComps <- array(0, dim = c(1, n_yrs, 1, n_sexes, 1))
  UseSrvLenComps <- array(0, dim = c(1, n_yrs, 1, 1))

  srv_age_yrs <- tmb_data$yrs_srv_age + 1
  ObsSrvAgeComps[1,srv_age_yrs,1,,1,1] <- tmb_data$data_srv_age
  ISS_SrvAgeComps[1,srv_age_yrs,1,1,1] <- tmb_data$effn_srv_age
  UseSrvAgeComps[1,srv_age_yrs,1,1] <- 1

  srv_len_yrs <- tmb_data$yrs_srv_len + 1
  for(s in 1:n_sexes) {
    ObsSrvLenComps[1,srv_len_yrs,1,,s,1] <- tmb_data$data_srv_len[,,3 - s]
    ISS_SrvLenComps[1,srv_len_yrs,1,s,1] <- tmb_data$effn_srv_len[,1,3 - s]
  } # end s loop
  UseSrvLenComps[1,srv_len_yrs,1,1] <- 1

  input_list <- Setup_Mod_SrvIdx_and_Comps(input_list = input_list,
                                           # data inputs
                                           ObsSrvIdx = ObsSrvIdx, # observed survey index
                                           ObsSrvIdx_SE = ObsSrvIdx_SE, # lognormal sd
                                           UseSrvIdx = UseSrvIdx, # whether survey indices are used
                                           ObsSrvAgeComps = ObsSrvAgeComps, # observed survey ages
                                           UseSrvAgeComps = UseSrvAgeComps, # whether survey ages are used
                                           ISS_SrvAgeComps = ISS_SrvAgeComps, # input sample size for survey ages
                                           ObsSrvLenComps = ObsSrvLenComps, # observed survey lengths
                                           UseSrvLenComps = UseSrvLenComps, # whether survey lengths are used
                                           ISS_SrvLenComps = ISS_SrvLenComps, # input sample size for survey lengths
                                           # Model options
                                           srv_idx_type = "abd", # abundance for the longline survey
                                           SrvAgeComps_LikeType = "Multinomial", # survey age composition likelihood
                                           SrvLenComps_LikeType = "Multinomial", # survey length composition likelihood
                                           SrvAgeComps_Type = "agg_Year_1-terminal_Fleet_1", # survey age comp type
                                           SrvLenComps_Type = "spltRspltS_Year_1-terminal_Fleet_1" # survey length comp type
  )

  # Setup fishery selectivity, retention and catchability
  # the mark-recapture fleet takes the same blocks, retention and selectivity parameters as the fishery
  ret_sel_input <- array(0, dim = c(1, 1, n_yrs, 1, n_ages, n_sexes, 2))
  for(y in 1:n_yrs) for(s in 1:n_sexes) ret_sel_input[1,1,y,1,,s,] <- tmb_data$retention[1,,3 - s]

  # mark-recapture index is in millions, so the q prior mean is 1e-6
  fish_q_prior <- data.frame(region = 1,
                             fleet = 2,
                             block = 1,
                             mu = 1e-6,
                             sd = 0.01
  )

  input_list <- Setup_Mod_Fishsel_and_Q(input_list = input_list,
                                        # fishery selectivity blocks (1975 - 1994, 1995 - 2021, 2022 - 2025)
                                        fish_sel_blocks = c("Block_1_Year_1-20_Fleet_1",
                                                            "Block_2_Year_21-47_Fleet_1",
                                                            "Block_3_Year_48-51_Fleet_1",
                                                            "Block_1_Year_1-20_Fleet_2",
                                                            "Block_2_Year_21-47_Fleet_2",
                                                            "Block_3_Year_48-51_Fleet_2"),
                                        # fishery selectivity form
                                        fish_sel_model = c("logist1_Fleet_1_Block_1",
                                                           "logist1_Fleet_1_Block_2",
                                                           "logist1_Fleet_1_Block_3",
                                                           "logist1_Fleet_2_Block_1",
                                                           "logist1_Fleet_2_Block_2",
                                                           "logist1_Fleet_2_Block_3"),
                                        fish_fixed_sel_pars_spec = c("est_all", "est_shared_f_1"), # mark-recapture shares fleet 1 selectivity, mapped off below
                                        # fishery catchability blocks, 3 for the fishery cpue and 1 for mark-recapture
                                        fish_q_blocks = c("Block_1_Year_1-20_Fleet_1",
                                                          "Block_2_Year_21-47_Fleet_1",
                                                          "Block_3_Year_48-51_Fleet_1",
                                                          "none_Fleet_2"),
                                        fish_q_spec = c("est_all", "est_all"), # estimate all catchability
                                        Use_fish_q_prior = 1, # use catchability prior
                                        fish_q_prior = fish_q_prior,
                                        use_fixed_ret_sel = c(1, 1), # fixed retention
                                        ret_sel_input = ret_sel_input # retention at age
  )

  # Setup survey selectivity and catchability
  input_list <- Setup_Mod_Srvsel_and_Q(input_list = input_list,
                                       # survey selectivity blocks (1975 - 1999, 2000 - 2016, 2017 - 2025)
                                       srv_sel_blocks = c("Block_1_Year_1-25_Fleet_1",
                                                          "Block_2_Year_26-42_Fleet_1",
                                                          "Block_3_Year_43-51_Fleet_1"),
                                       # survey selectivity form
                                       srv_sel_model = c("logist1_Fleet_1_Block_1",
                                                         "logist1_Fleet_1_Block_2",
                                                         "logist1_Fleet_1_Block_3"),
                                       srv_fixed_sel_pars_spec = "est_all", # estimate longline survey selectivity
                                       # survey catchability blocks, following the selectivity blocks
                                       srv_q_blocks = c("Block_1_Year_1-25_Fleet_1",
                                                        "Block_2_Year_26-42_Fleet_1",
                                                        "Block_3_Year_43-51_Fleet_1"),
                                       srv_q_spec = "est_all", # estimate catchability
                                       t_srv = array(tmb_data$srv_month, dim = c(1, 1, 1)) # survey timing (july)
  )

  # Setup weighting stuff
  # the tmb model adds 1e-5 inside every log
  input_list <- Setup_Mod_Weighting(input_list = input_list,
                                    addtocomp = 1e-5, # constant added to compositions
                                    comp_const_obs = 1, # added to observed compositions too
                                    addtofishidx = 1e-5, # constant added to fishery indices
                                    addtosrvidx = 1e-5, # constant added to survey indices
                                    Wt_Catch = 1,
                                    Wt_FishIdx = 1,
                                    Wt_SrvIdx = 1,
                                    Wt_Rec = 1,
                                    Wt_F = 1, # F penalty weight is in its sd above
                                    Wt_Tagging = 0
  )

  # fix fishery selectivity, and male longline survey selectivity in the first block
  input_list$map$fish_fixed_sel_pars <- factor(rep(NA, length(input_list$par$fish_fixed_sel_pars)))
  map_srv_sel <- array(as.numeric(input_list$map$srv_fixed_sel_pars), dim = dim(input_list$par$srv_fixed_sel_pars)) # region, par, block, sex, fleet
  map_srv_sel[1,,1,2,1] <- NA
  input_list$map$srv_fixed_sel_pars <- factor(map_srv_sel)

  return(input_list)
} # end function

# Fix Values --------------------------------------------------------------

seed_nsei_sablefish_tmb <- function(input_list, tmb_par) {

  # the tmb rec penalty is centered at sigma^2 / 2, so that goes into mean recruitment
  # and comes off every deviation. the plus group has no tmb deviation, so it is fixed at -sigma^2 / 2
  bias_shift <- exp(tmb_par$log_sigma_r)^2 / 2
  n_init_devs <- length(tmb_par$log_rinit_devs)

  # recruitment and initial age structure
  input_list$par$ln_global_R0[] <- tmb_par$log_rbar + bias_shift
  input_list$par$ln_rinit[] <- tmb_par$log_rinit + bias_shift
  input_list$par$ln_RecDevs[1,1,] <- tmb_par$log_rec_devs - bias_shift
  for(s in 1:2) {
    input_list$par$ln_InitDevs[1,1,1:n_init_devs,s] <- tmb_par$log_rinit_devs - bias_shift # same for both sexes
    input_list$par$ln_InitDevs[1,1,n_init_devs + 1,s] <- -bias_shift # plus group
  } # end s loop

  # fishing mortality (fleet 1 only)
  input_list$par$ln_F_mean[1,1,1] <- tmb_par$log_Fbar
  input_list$par$ln_F_devs[1,,1,1] <- tmb_par$log_F_devs

  # catchability
  input_list$par$ln_srv_q[1,1:3,1] <- tmb_par$srv_logq # longline survey
  input_list$par$ln_fish_q[1,1:3,1] <- tmb_par$fsh_logq # fishery cpue
  input_list$par$ln_fish_q[1,1,2] <- tmb_par$mr_logq + log(1e-6) # mark-recapture, index in millions

  # selectivity, tmb logistic is in age - 2 so a50 moves up by 2
  for(b in 1:3) {
    for(s in 1:2) {
      input_list$par$fish_fixed_sel_pars[1,1,b,s,] <- log(exp(tmb_par$log_fsh_slx_pars[b,1,3 - s]) + 2) # both fishery fleets
      input_list$par$fish_fixed_sel_pars[1,2,b,s,] <- tmb_par$log_fsh_slx_pars[b,2,3 - s]
      input_list$par$srv_fixed_sel_pars[1,1,b,s,1] <- log(exp(tmb_par$log_srv_slx_pars[b,1,3 - s]) + 2)
      input_list$par$srv_fixed_sel_pars[1,2,b,s,1] <- tmb_par$log_srv_slx_pars[b,2,3 - s]
    } # end s loop
  } # end b loop

  return(input_list)
} # end function
