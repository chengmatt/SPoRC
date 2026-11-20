# Purpose: To set up the September 2026 EBS pollock Model 26.0 (Rceattle) in SPoRC and start it at the Rceattle estimate
# Creator: Matthew LH. Cheng
# Date: 11/19/26

# Setup Model -------------------------------------------------------------

build_ebs_pollock_m26_input <- function(inputs, rce) {

  years <- inputs$years
  n_yrs <- length(years)
  n_ages <- length(inputs$ages)
  n_srv <- inputs$n_srv_fleets

  # Initialize model dimensions and data list
  input_list <- Setup_Mod_Dim(years = years, # vector of years (1964 - 2024)
                              ages = inputs$ages, # vector of ages
                              lens = NA, # no lengths
                              n_regions = 1, # number of regions
                              n_sexes = 1, # number of sexes
                              n_fish_fleets = 1, # number of fishery fleets
                              n_srv_fleets = n_srv, # number of survey fleets == 1 bottom trawl, == 2 acoustic trawl, == 3 vessel of opportunity, == 4 acoustic age 1
                              bias_correct_pe = "none", # no bias correction in rceattle
                              verbose = FALSE
  )

  # Setup recruitment stuff
  input_list <- Setup_Mod_Rec(input_list = input_list, # input data list from above
                              rec_model = "mean_rec", # recruitment model
                              sigmaR_spec = "fix", # fix sigmaR
                              ln_sigmaR = array(log(1), dim = c(2, 1, 1)), # sigmaR of 1 for initial ages and recruitment
                              sigmaR_switch = 1, # one sigmaR for all years
                              do_rec_bias_ramp = 0, # no bias ramp
                              init_age_strc = 1, # geometric series to derive initial age structure
                              equil_init_age_strc = 2, # initial deviations on every age including the plus group
                              ln_global_R0 = rce$ln_R0, # mean recruitment
                              t_spawn = 3 / 12 # spawn in march
  )

  # Setup biological stuff
  input_list <- Setup_Mod_Biologicals(input_list = input_list,
                                      WAA = inputs$WAA, # spawning weight at age
                                      WAA_fish = inputs$WAA_fish, # fishery weight at age
                                      WAA_srv = inputs$WAA_srv, # survey weight at age
                                      MatAA = inputs$MatAA, # maturity at age
                                      fit_lengths = 0, # no lengths
                                      M_spec = "fix", # fix natural mortality
                                      Fixed_natmort = inputs$Fixed_natmort # natural mortality at age
  )

  # Setup movement and tagging
  input_list <- Setup_Mod_Movement(input_list = input_list,
                                   use_fixed_movement = 1,
                                   Fixed_Movement = NA,
                                   do_recruits_move = 0
  )

  input_list <- Setup_Mod_Tagging(input_list = input_list, use_conv_fish_tagging = 0)

  # Setup catch and fishing mortality
  input_list <- Setup_Mod_Catch_and_F(input_list = input_list,
                                      ObsCatch = inputs$ObsCatch, # catch in thousand tons
                                      UseCatch = inputs$UseCatch,
                                      Use_F_pen = 0, # no F penalty
                                      ln_F_mean_spec = "fix", # no mean F, so the deviations are log F
                                      sigmaC_spec = "fix", # fix catch sd
                                      ln_sigmaC = array(log(0.05), dim = c(1, n_yrs, 1, 1)) # catch sd of 0.05
  )

  # Setup fishery indices and compositions
  input_list <- Setup_Mod_FishIdx_and_Comps(input_list = input_list,
                                            ObsFishIdx = inputs$ObsFishIdx, # fishery cpue, 1965 - 1976
                                            ObsFishIdx_SE = inputs$ObsFishIdx_SE, # standard error
                                            UseFishIdx = inputs$UseFishIdx,
                                            ObsFishAgeComps = inputs$ObsFishAgeComps,
                                            UseFishAgeComps = inputs$UseFishAgeComps,
                                            ISS_FishAgeComps = inputs$ISS_FishAgeComps,
                                            ObsFishLenComps = array(NA, dim = c(1, n_yrs, 1, 1, 1, 1)),
                                            UseFishLenComps = array(0, dim = c(1, n_yrs, 1, 1)),
                                            ISS_FishLenComps = array(0, dim = c(1, n_yrs, 1, 1, 1)),
                                            fish_idx_type = "biom", # biomass
                                            FishIdx_LikeType = "normal", # normal likelihood
                                            FishAgeComps_LikeType = "Multinomial",
                                            FishLenComps_LikeType = "none",
                                            FishAgeComps_Type = "agg_Year_1-terminal_Fleet_1",
                                            FishLenComps_Type = "none_Year_1-terminal_Fleet_1"
  )

  # Setup survey indices and compositions
  input_list <- Setup_Mod_SrvIdx_and_Comps(input_list = input_list,
                                           ObsSrvIdx = inputs$ObsSrvIdx,
                                           ObsSrvIdx_SE = inputs$ObsSrvIdx_SE,
                                           UseSrvIdx = inputs$UseSrvIdx,
                                           ObsSrvAgeComps = inputs$ObsSrvAgeComps,
                                           UseSrvAgeComps = inputs$UseSrvAgeComps,
                                           ISS_SrvAgeComps = inputs$ISS_SrvAgeComps,
                                           ObsSrvLenComps = array(NA, dim = c(1, n_yrs, 1, 1, 1, n_srv)),
                                           UseSrvLenComps = array(0, dim = c(1, n_yrs, 1, n_srv)),
                                           ISS_SrvLenComps = array(0, dim = c(1, n_yrs, 1, 1, n_srv)),
                                           srv_idx_type = c("biom", "biom", "biom", "abd"), # biomass, and abundance for the age 1 index
                                           srv_idx_ages = list(NULL, 2:15, 2:15, 1), # acoustic indices at ages 2 - 15, age 1 index at age 1
                                           SrvIdx_LikeType = c("mvn", "lognormal", "normal", "lognormal"), # bottom trawl is multivariate normal
                                           SrvIdx_Cov = list(inputs$SrvIdx_Cov, NULL, NULL, NULL), # bottom trawl covariance
                                           SrvAgeComps_LikeType = c("Multinomial", "Multinomial", "none", "none"),
                                           SrvLenComps_LikeType = c("none", "none", "none", "none"),
                                           SrvAgeComps_Type = c("agg_Year_1-terminal_Fleet_1", "agg_Year_1-terminal_Fleet_2",
                                                                "none_Year_1-terminal_Fleet_3", "none_Year_1-terminal_Fleet_4"),
                                           SrvLenComps_Type = c("none_Year_1-terminal_Fleet_1", "none_Year_1-terminal_Fleet_2",
                                                                "none_Year_1-terminal_Fleet_3", "none_Year_1-terminal_Fleet_4"),
                                           SrvAgeComps_bins = list(NULL, 2:15, NULL, NULL) # acoustic ages 2 - 15
  )

  # Setup fishery selectivity and catchability
  # the first year's increment from the fixed base is penalized at weight 1, sd 1 / sqrt(2)
  fish_selex_prior <- data.frame(region = 1,
                                 fleet = 1,
                                 block = 1,
                                 sex = 1,
                                 par = 1:12, # one coefficient per age to 11, and one for 12 - 15
                                 mu = exp(rce$fish_sel_coff), # rceattle base, on the natural scale
                                 sd = 1 / sqrt(2)
  )

  input_list <- Setup_Mod_Fishsel_and_Q(input_list = input_list,
                                        fish_sel_model = "nonparlog_Fleet_1", # non-parametric on the log scale
                                        cont_tv_fish_sel = "rw_Fleet_1", # random walk
                                        fish_sel_blocks = "none_Fleet_1",
                                        fish_q_blocks = "none_Fleet_1",
                                        fish_fixed_sel_pars_spec = "est_all",
                                        fish_q_spec = "est_all",
                                        fishsel_pe_pars_spec = "fix", # fix random walk sd
                                        fish_sel_devs_spec = "est_all",
                                        fishsel_rw_init_sigma = NA, # walk starts at zero
                                        fish_sel_nonpar_est_bins = list(list(list(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12:15))), # ages 12 - 15 share one coefficient
                                        fishsel_devs_shared_bins = list(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12:15),
                                        Use_fish_selex_prior = 1,
                                        fish_selex_prior = fish_selex_prior
  )

  # Setup survey selectivity and catchability
  input_list <- Setup_Mod_Srvsel_and_Q(input_list = input_list,
                                       srv_sel_model = c("logist1_Fleet_1", "nonparlog_Fleet_2", "nonparlog_Fleet_3", "logist1_Fleet_4"),
                                       cont_tv_srv_sel = c("iid_Fleet_1", "rw_Fleet_2", "rw_Fleet_3", "none_Fleet_4"),
                                       srv_sel_blocks = c("none_Fleet_1", "none_Fleet_2", "none_Fleet_3", "none_Fleet_4"),
                                       srv_q_blocks = c("none_Fleet_1", "none_Fleet_2", "none_Fleet_3", "none_Fleet_4"),
                                       srv_fixed_sel_pars_spec = c("est_all", "est_all", "est_shared_f_2", "fix"), # vessel of opportunity shares acoustic selectivity
                                       srv_sel_devs_spec = c("est_all", "est_all", "est_shared_f_2", "fix"),
                                       srvsel_pe_pars_spec = c("fix", "fix", "fix", "fix"),
                                       srvsel_pe_wt = c(0, 0, 0, 0), # survey deviations only enter through the penalties
                                       srv_q_spec = c("est_all", "est_all", "est_all", "est_all"),
                                       srv_q_type = c("arith", "est", "est", "geo"), # solved for the bottom trawl and age 1 index
                                       srv_sel_nonpar_est_bins = list(NULL,
                                                                      list(list(1, 2, 3, 4, 5, 6, 7, 8:15)),
                                                                      list(list(1:2, 3, 4, 5, 6, 7, 8:15)),
                                                                      NULL), # acoustic ages 2 - 8, and 9 - 15 follow age 8
                                       srvsel_devs_shared_bins = list(1, 2, 3, 4, 5, 6, 7, 8:15),
                                       srv_sel_bin_dev_bins = list(1, NULL, NULL, NULL), # bottom trawl age 1 is its own annual value
                                       cont_tv_srvsel_bin_devs = c("rw", "none", "none", "none"),
                                       t_srv = array(c(0.5, 0.5, 0, 0.5), dim = c(1, 1, n_srv)) # july, july, start of year, july
  )

  # Setup model weighting and selectivity penalties
  fish_pen_wts <- list(smooth_dome = 12.5, # declines in log selectivity at every age
                       smooth_bin_curve = rep(1 / 60, n_yrs), # curvature
                       smooth_yr_diff = c(0, rep(1 / (2 * 0.5^2), n_yrs - 1)), # random walk sd of 0.5
                       normalize = FALSE
  )

  bts_rw <- as.numeric(years >= 1983) * 2 # bottom trawl year to year change, ages 3 - 14
  ats_shape <- as.numeric(years >= 1994) # acoustic increases over ages 2 - 15
  ats_curve <- as.numeric(years >= 1995) # acoustic curvature
  ats_rw <- as.numeric(years >= 1995) / (2 * 0.138^2) # acoustic random walk sd of 0.138

  srv_pen_wts <- list(list(smooth_yr_diff = bts_rw, normalize = FALSE,
                           bin_range = list(smooth_yr_diff = c(3, 14))), # bottom trawl
                      list(smooth_bin_diff = ats_shape, smooth_dome = -ats_shape, smooth_bin_curve = ats_curve,
                           smooth_yr_diff = ats_rw, normalize = FALSE,
                           bin_range = list(smooth_bin_diff = c(2, 15), smooth_dome = c(2, 15))), # acoustic trawl
                      list(), # vessel of opportunity shares the acoustic curve
                      list() # acoustic age 1
  )

  input_list <- Setup_Mod_Weighting(input_list = input_list,
                                    addtocomp = 1e-3, # constant added to compositions
                                    addtosrvidx = 0, # nothing added to the indices
                                    addtofishidx = 0,
                                    fish_sel_pen_wts = fish_pen_wts,
                                    srv_sel_pen_wts = srv_pen_wts,
                                    Wt_Catch = 1,
                                    Wt_FishIdx = 1,
                                    Wt_SrvIdx = 1,
                                    Wt_Rec = 1,
                                    Wt_Init_Rec = 1,
                                    Wt_F = 0,
                                    Wt_Tagging = 0
  )

  # Setup mapping
  # bottom trawl deviations from 1982, when the survey starts
  map_srvdev <- array(as.numeric(input_list$map$ln_srvsel_devs), dim = dim(input_list$par$ln_srvsel_devs)) # region, year, par, sex, fleet
  map_srvdev[1,years < 1982,,1,1] <- NA

  # acoustic deviations at ages 2 - 8 every year from 1995, and the vessel of opportunity shares them
  map_srvdev[1,,,1,2:3] <- NA
  k <- max(map_srvdev, na.rm = TRUE) + 1
  for(y in which(years >= 1995)) {
    for(g in list(2, 3, 4, 5, 6, 7, 8:15)) {
      map_srvdev[1,y,g,1,2] <- k
      k <- k + 1
    } # end g loop
  } # end y loop
  map_srvdev[1,,,1,3] <- map_srvdev[1,,,1,2]
  input_list$map$ln_srvsel_devs <- factor(map_srvdev)
  input_list$data$map_ln_srvsel_devs <- map_srvdev

  # fishery deviations every year from 1965, ages 12 - 15 move together
  map_fishdev <- array(NA, dim = dim(input_list$par$ln_fishsel_devs)) # region, year, bin, sex, fleet
  k <- 1
  for(y in 2:n_yrs) {
    for(g in list(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12:15)) {
      map_fishdev[1,y,g,1,1] <- k
      k <- k + 1
    } # end g loop
  } # end y loop
  input_list$map$ln_fishsel_devs <- factor(map_fishdev)
  input_list$data$map_ln_fishsel_devs <- map_fishdev

  # acoustic age 1 kept at zero, and no base logistic pair for the bottom trawl
  map_srvpar <- array(as.numeric(input_list$map$srv_fixed_sel_pars), dim = dim(input_list$par$srv_fixed_sel_pars)) # region, par, block, sex, fleet
  map_srvpar[1,1,1,1,2:3] <- NA
  map_srvpar[1,1:2,1,1,1] <- NA
  input_list$map$srv_fixed_sel_pars <- factor(map_srvpar)

  # bottom trawl age 1 from 1982, with a fixed random walk sd
  map_bindev <- array(as.numeric(input_list$map$ln_srvsel_bin_devs), dim = dim(input_list$par$ln_srvsel_bin_devs))
  map_bindev[1,years < 1982,,1,1] <- NA
  input_list$map$ln_srvsel_bin_devs <- factor(map_bindev)
  input_list$data$map_ln_srvsel_bin_devs <- map_bindev
  input_list$map$srvsel_bin_devs_pe_pars <- factor(rep(NA, length(input_list$par$srvsel_bin_devs_pe_pars)))

  return(input_list)

} # end build_ebs_pollock_m26_input

# Fix Values --------------------------------------------------------------

seed_ebs_pollock_m26_rceattle <- function(input_list, rce) {

  years <- input_list$data$years
  n_yrs <- length(years)
  n_ages <- length(input_list$data$ages)

  # recruitment, initial age structure and fishing mortality
  input_list$par$ln_RecDevs[1,1,] <- rce$rec_dev
  input_list$par$ln_InitDevs[1,1,,1] <- rce$init_dev
  input_list$par$ln_F_devs[1,,1,1] <- rce$log_F

  # fishery: the base plus the first year's increment, then the running sum of the later increments
  bin_f <- c(1:12, 12, 12, 12) # bin of each age
  input_list$par$fish_fixed_sel_pars[1,,1,1,1] <- (rce$fish_sel_coff + rce$fish_sel_coff_dev[1,])[bin_f]
  input_list$par$ln_fishsel_devs[1,2:n_yrs,,1,1] <- apply(rce$fish_sel_coff_dev[-1,], 2, cumsum)[,bin_f]
  input_list$par$fishsel_pe_pars[1,,1,1] <- log(1 / sqrt(2)) # sum of squares of the increments

  # acoustic: log selectivity relative to age 1, where age 1 before it is zeroed is 15 less the sum of ages 2 - 15
  log_sel_ats <- log(cbind(n_ages - rowSums(rce$ats_sel[,2:n_ages]), rce$ats_sel[,2:n_ages]))
  log_sel_ats <- log_sel_ats - log_sel_ats[,1]
  y_94 <- which(years == 1994)
  for(f in 2:3) {
    input_list$par$srv_fixed_sel_pars[1,,1,1,f] <- log_sel_ats[y_94,]
    for(y in which(years > 1994)) input_list$par$ln_srvsel_devs[1,y,,1,f] <- log_sel_ats[y,] - log_sel_ats[y_94,]
  } # end f loop

  # bottom trawl: logistic at age + 0.5 with a multiplicative a50 deviation, so a50 - 0.5 on the ages here
  input_list$par$srv_fixed_sel_pars[1,1:2,1,1,1] <- 0
  input_list$par$ln_srvsel_devs[1,,1,1,1] <- log(rce$bts_a50 * exp(rce$bts_a50_dev) - 0.5)
  input_list$par$ln_srvsel_devs[1,,2,1,1] <- rce$bts_log_slope + rce$bts_log_slope_dev

  # bottom trawl age 1: log selectivity is b1 * exp(dev) in rceattle, with weight 8 on the change in dev
  input_list$par$ln_srvsel_bin_devs[1,,1,1,1] <- rce$bts_age1 * exp(rce$bts_age1_dev)
  input_list$par$srvsel_bin_devs_pe_pars[1,,1,1] <- log(abs(rce$bts_age1) / 4) # weight 8 / b1^2 on the change in log selectivity

  # catchability
  input_list$par$ln_srv_q[1,,2] <- rce$log_q["ats"]
  input_list$par$ln_srv_q[1,,3] <- rce$log_q["avo"]
  input_list$par$ln_fish_q[1,,1] <- rce$log_q["cpue"]

  return(input_list)

} # end seed_ebs_pollock_m26_rceattle
