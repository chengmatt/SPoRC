# Purpose: To bridge to the spatial assessment in Cheng and Marsh et al. 2025  using SPoRC
# Creator: Matthew LH. Cheng
# Date Created: 2/5/25

# Set up ------------------------------------------------------------------
unloadNamespace("SPoRC")
library(here)
library(tidyverse)
library(RTMB)
library(SPoRC)
devtools::load_all(here("R"))

# Read in data inputs
ageing_dat <- dget(here("dev", "dev_data",'2023 Base (23.5)_final model', 'test.rdat')) # for getting ageing error
dat <- readRDS(here("dev", "dev_data", "Spatial Sablefish Model", "data_1_area.RDS"))

# Munging data inputs -----------------------------------------------------
n_seas <- 1

# Setup biological data
waa <- array(0, dim = c(1, dat$n_regions, length(dat$years), n_seas, length(dat$ages), 2))
mataa <- array(0, dim = c(1, dat$n_regions, length(dat$years), n_seas, length(dat$ages), 2))
agelen <- array(0, dim = c(1, dat$n_regions, length(dat$years), n_seas, length(dat$length_bins), length(dat$ages), 2))

for(y in 1:length(dat$years)) {

  # weight at age
  waa[1,1,y,1,,1] <- dat$female_mean_weight_by_age[,y]
  waa[1,1,y,1,,2] <- dat$male_mean_weight_by_age[,y]

  # maturity at age
  mataa[1,1,y,1,,1] <- dat$maturity[,y]
  mataa[1,1,y,1,,2] <- dat$maturity[,y]

  # age length transition
  agelen[1,1,y,1,,,1] <- t(dat$female_age_length_transition[,,y])
  agelen[1,1,y,1,,,2] <- t(dat$male_age_length_transition[,,y])

} # end y loop


# Catch
obscatch <- array(0, dim = c(dat$n_regions, length(dat$years), n_seas, 2))
usecatch <- array(0, dim = c(dat$n_regions, length(dat$years), n_seas, 2))

obscatch[1,,1,1] <- as.vector(dat$fixed_fishery_catch) # fixed gear catches
obscatch[1,,1,2] <- as.vector(dat$trwl_fishery_catch) # trawl gear catches
obscatch[1,1:3,1,2] <- 0 # set first three years catch as 0
usecatch[obscatch != 0] <- 1 # use catch indicator

# fishery index
obsfishidx <- array(0, dim = c(dat$n_regions, length(dat$years), n_seas, 2))
obsfishidx_se <- array(0, dim = c(dat$n_regions, length(dat$years),n_seas, 2))
usefishidx <- array(0, dim = c(dat$n_regions, length(dat$years), n_seas, 2))

# fishery ages
obsfishage <- array(0, dim = c(dat$n_regions, length(dat$years), n_seas, length(dat$ages), 2, 2))
usefishage <- array(0, dim = c(dat$n_regions, length(dat$years),n_seas,  2))
issfishage <- array(0, dim = c(dat$n_regions, length(dat$years), n_seas, 2, 2))

for(y in 1:length(dat$years)) {
  obsfishage[1,y,1,,1,1] <- dat$obs_fixed_catchatage[-c(1:30),,y] # females
  obsfishage[1,y,1,,2,1] <- dat$obs_fixed_catchatage[c(1:30),,y] # males
}

usefishage[1,,1,1] <- dat$fixed_catchatage_indicator # use catch at age indicator
issfishage[1,which(dat$fixed_catchatage_indicator == 1),1,1,1] <- 250 # iss weighting for catch at age

# fishery lengths
obsfishlen <- array(0, dim = c(dat$n_regions, length(dat$years), n_seas, length(dat$length_bins), 2, 2))
usefishlen <- array(0, dim = c(dat$n_regions, length(dat$years), n_seas, 2))
issfishlen <- array(0, dim = c(dat$n_regions, length(dat$years), n_seas, 2, 2))

for(y in 1:length(dat$years)) {

  # fixed gear
  obsfishlen[1,y,1,,1,1] <- dat$obs_fixed_catchatlgth[-c(1:30),,y] # females
  obsfishlen[1,y,1,,2,1] <- dat$obs_fixed_catchatlgth[c(1:30),,y] # males

  # trawl gear
  obsfishlen[1,y,1,,1,2] <- dat$obs_trwl_catchatlgth[-c(1:30),,y] # females
  obsfishlen[1,y,1,,2,2] <- dat$obs_trwl_catchatlgth[c(1:30),,y] # males
}

# fixed gear
usefishlen[1,,1,1] <- dat$fixed_catchatlgth_indicator # use catch at length indicator
issfishlen[1,which(dat$fixed_catchatlgth_indicator == 1),1,1,1] <- 200 # iss weighting for catch at length

# trawl gear
usefishlen[1,,1,2] <- dat$trwl_catchatlgth_indicator # use catch at length indicator
issfishlen[1,which(dat$trwl_catchatlgth_indicator == 1),1,1,2] <- 200 # iss weighting for catch at length

# survey index
obssrvidx <- array(0, dim = c(dat$n_regions, length(dat$years), n_seas, 2))
obssrvidx_se <- array(0, dim = c(dat$n_regions, length(dat$years), n_seas, 2))
usesrvidx <- array(0, dim = c(dat$n_regions, length(dat$years), n_seas, 2))

obssrvidx[] <- dat$obs_srv_bio # survey index
obssrvidx_se[] <- dat$obs_srv_se # survey index se
usesrvidx[] <- dat$srv_bio_indicator # use indicator

# survey ages
obssrvage <- array(0, dim = c(dat$n_regions, length(dat$years), n_seas, length(dat$ages), 2, 2))
usesrvage <- array(0, dim = c(dat$n_regions, length(dat$years), n_seas, 2))
isssrvage <- array(0, dim = c(dat$n_regions, length(dat$years), n_seas, 2, 2))

for(y in 1:length(dat$years)) {

  # domestic
  obssrvage[1,y,1,,1,1] <- dat$obs_srv_catchatage[-c(1:30),,y,1] # females
  obssrvage[1,y,1,,2,1] <- dat$obs_srv_catchatage[c(1:30),,y,1] # males

  # japanese
  obssrvage[1,y,1,,1,2] <- dat$obs_srv_catchatage[-c(1:30),,y,2] # females
  obssrvage[1,y,1,,2,2] <- dat$obs_srv_catchatage[c(1:30),,y,2] # males
}

usesrvage[1,,1,] <- dat$srv_catchatage_indicator # use catch at age indicator
isssrvage[] <- 300 # iss weighting for catch at age

# survey lengths
obssrvlen <- array(0, dim = c(dat$n_regions, length(dat$years), n_seas, length(dat$length_bins), 2, 2))
usesrvlen <- array(0, dim = c(dat$n_regions, length(dat$years), n_seas,2))
isssrvlen <- array(0, dim = c(dat$n_regions, length(dat$years),n_seas,  2, 2))

# ISS weights are kind of arbirtrarly determined -- just trying to match relative weights in spatial assessment given to indices

# Setup model -------------------------------------------------------------

# Initialize model dimensions and data list
input_list <- Setup_Mod_Dim(years = 1:length(dat$years),
                            # vector of years (1 - 62)
                            ages = 1:length(dat$ages),
                            # vector of ages (1 - 30)
                            lens = dat$length_bins,
                            # number of lengths (41 - 99)
                            n_regions = dat$n_regions,
                            # number of regions (5)
                            n_sexes = 2,
                            # number of sexes (2)
                            n_fish_fleets = 2,
                            # number of fishery fleet (2)
                            n_srv_fleets = 2,
                            # number of survey fleets (2)
                            n_seas = 1,
                            n_pop = 1,
                            verbose = TRUE
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
                            ln_sigmaR = array(log(c(0.4, 1.2)), dim = c(2, 1, 1)),
                            # values to fix sigmaR at, or starting values
                            ln_global_R0 = log(30)
)

# Setup biological stuff (using defaults for other stuff)
input_list <- Setup_Mod_Biologicals(input_list = input_list,
                                    WAA = waa, # weight at age
                                    MatAA = mataa, # maturity at age
                                    AgeingError = ageing_dat$age_error,
                                    # ageing error matrix
                                    fit_lengths = 1, # fitting lengths
                                    SizeAgeTrans = agelen,
                                    # size age transition matrix
                                    M_spec = "fix", # fix natural mortality
                                    Fixed_natmort = array(0.104884, dim = c(1, input_list$data$n_regions,
                                                                            length(input_list$data$years),
                                                                            length(input_list$data$ages),
                                                                            input_list$data$n_sexes))
                                    # values to fix natural mortality at
)

# setting up movement parameterization
input_list <- Setup_Mod_Movement(input_list = input_list)

# setting up tagging parameterization
input_list <- Setup_Mod_Tagging(input_list = input_list, use_conv_fish_tagging = c(0,0))

# setting up catch data
input_list <- Setup_Mod_Catch_and_F(input_list = input_list,
                                    # Data inputs
                                    ObsCatch = obscatch,
                                    UseCatch = usecatch,
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
                                          ObsFishIdx = obsfishidx,
                                          ObsFishIdx_SE = obsfishidx_se,
                                          UseFishIdx =  usefishidx,
                                          ObsFishAgeComps = obsfishage,
                                          UseFishAgeComps = usefishage,
                                          ISS_FishAgeComps = issfishage,
                                          ObsFishLenComps = obsfishlen,
                                          UseFishLenComps = usefishlen,
                                          ISS_FishLenComps = issfishlen,

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
)

# Survey Indices and Compositions
input_list <- Setup_Mod_SrvIdx_and_Comps(input_list = input_list,
                                         # data inputs
                                         ObsSrvIdx = obssrvidx,
                                         ObsSrvIdx_SE = obssrvidx_se,
                                         UseSrvIdx =  usesrvidx,
                                         ObsSrvAgeComps = obssrvage,
                                         UseSrvAgeComps = usesrvage,
                                         ISS_SrvAgeComps = isssrvage,
                                         ObsSrvLenComps = obssrvlen,
                                         UseSrvLenComps = usesrvlen,
                                         ISS_SrvLenComps = isssrvlen,

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
                                       c("est_all",
                                         "est_all"),

                                     # Starting values for survey catchability
                                     ln_srv_q = array(9,
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

parameters$ln_fish_fixed_sel_pars[,,,,1] <- log(1) # some more informative starting values
parameters$ln_fish_fixed_sel_pars[,,,,2] <- log(5)

# Fit model
st <- Sys.time()
sabie_rtmb_model <- fit_model(data,
                              parameters,
                              mapping,
                              random = NULL,
                              newton_loops = 5,
                              silent = FALSE,
                              do_optim = T
)
en <- Sys.time()
print(en - st)

sabie_rtmb_model$sd_rep <- sdreport(sabie_rtmb_model)

rep <- sabie_rtmb_model$rep
sd_rep <- sabie_rtmb_model$sd_rep

saveRDS(data, here("dev", "dev_output", "1_Region_Model_Sablefish_SptComparison", "data.RDS"))
saveRDS(sd_rep, here("dev", "dev_output", "1_Region_Model_Sablefish_SptComparison", "sd_rep.RDS"))
saveRDS(rep, here("dev", "dev_output", "1_Region_Model_Sablefish_SptComparison", "rep.RDS"))

