# Purpose: To bridge sabieRTMB v3 to a single region model - 2024 assessment
# Creator: Matthew LH. Cheng
# Date Created: 1/27/25

library(here)
library(R2admb)
library(tidyverse)
library(RTMB)
library(SPoRC)
devtools::load_all(here("R"))

# Version 1 (Original Model) ----------------------------------------------
# Read in data
tem_dat <- dget(here("dev", "dev_data", '2024 Base (23.5)_final model', 'tem.rdat'))
ageing_dat <- dget(here("dev", "dev_data",'2023 Base (23.5)_final model', 'test.rdat')) # for getting ageing error
tem_admb_dat <- readLines(here("dev", "dev_data","2024 Base (23.5)_final model", "tem_2024_na_wh.dat"))
tem_rep <- readLines(here("dev", "dev_data",'2024 Base (23.5)_final model', 'sable.rep'))
tem_par <- read_pars(here("dev", "dev_data",'2024 Base (23.5)_final model', 'tem'))
rtmb_data <- readRDS(here("dev", "dev_data", 'sgl_rg_sable_data.RDS'))

# Prepare Data and Inputs -------------------------------------------------

# Initialize model dimensions and data list
input_list <- Setup_Mod_Dim(years = as.numeric((rownames(tem_dat$t.series))), # vector of years
                            ages = 1:30, # vector of ages
                            lens = seq(41,99,2), # number of lengths
                            n_regions = 1, # number of regions
                            n_sexes = 2, # number of sexes == 1, female, == 2 male
                            n_fish_fleets = 2, # number of fishery fleet == 1, fixed gear, == 2 trawl gear
                            n_srv_fleets = 3, # number of survey fleets
                            verbose = TRUE
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
                            ln_sigmaR = log(c(0.4, 1.2)),
                            rec_model = "mean_rec", # recruitment model
                            sigmaR_spec = "fix_early_est_late", # fix early sigmaR, estiamte late sigmaR
                            InitDevs_spec = NULL, # estimate all initial deviations
                            RecDevs_spec = NULL, # stiamte all recruitment deivations
                            # sexratio = as.vector(c(0.5, 0.5)), # recruitment sex ratio
                            init_age_strc = 2,
                            init_F_prop = 0.1
)

# Specificying natural mortality fixed array
fixed_natmort <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), length(input_list$data$ages), input_list$data$n_sexes))
fixed_natmort[,,,1] <- 0.1134156 # fix female M
fixed_natmort[,,,2] <- 0.1052175 # fix male M

input_list <- Setup_Mod_Biologicals(input_list = input_list,
                                    # Data inputs
                                    WAA = rtmb_data$WAA, # weight-at-age
                                    MatAA = rtmb_data$MatAA, # maturity at age
                                    AgeingError = as.matrix(ageing_dat$age_error), # ageing error
                                    SizeAgeTrans = rtmb_data$SizeAgeTrans, # size age transition matrix
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
input_list <- Setup_Mod_Catch_and_F(input_list = input_list,
                                    # Data inputs
                                    ObsCatch = rtmb_data$ObsCatch,
                                    UseCatch = rtmb_data$UseCatch,
                                    # Model options
                                    Use_F_pen = 1, # whether to use f penalty, == 0 don't use, == 1 use
                                    sigmaC_spec = 'fix'
)

# Setup fishery indices and compositions
input_list <- Setup_Mod_FishIdx_and_Comps(input_list = input_list,
                                          # data inputs
                                          ObsFishIdx = rtmb_data$ObsFishIdx,
                                          ObsFishIdx_SE = rtmb_data$ObsFishIdx_SE,
                                          UseFishIdx =  rtmb_data$UseFishIdx,
                                          ObsFishAgeComps = rtmb_data$ObsFishAgeComps,
                                          UseFishAgeComps = rtmb_data$UseFishAgeComps,
                                          ISS_FishAgeComps = rtmb_data$ISS_FishAgeComps,
                                          ObsFishLenComps = rtmb_data$ObsFishLenComps,
                                          UseFishLenComps = rtmb_data$UseFishLenComps,
                                          ISS_FishLenComps = rtmb_data$ISS_FishLenComps,

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
                                         ObsSrvIdx = rtmb_data$ObsSrvIdx,
                                         ObsSrvIdx_SE = rtmb_data$ObsSrvIdx_SE,
                                         UseSrvIdx =  rtmb_data$UseSrvIdx,
                                         ObsSrvAgeComps = rtmb_data$ObsSrvAgeComps,
                                         ISS_SrvAgeComps = rtmb_data$ISS_SrvAgeComps,
                                         UseSrvAgeComps = rtmb_data$UseSrvAgeComps,
                                         ObsSrvLenComps = rtmb_data$ObsSrvLenComps,
                                         UseSrvLenComps = rtmb_data$UseSrvLenComps,
                                         ISS_SrvLenComps = rtmb_data$ISS_SrvLenComps,

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
input_list$map$ln_fish_fixed_sel_pars <- factor(c(1:7, 2, 8:11, rep(12:13,3), rep(c(14,13),3)))

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
                                     srv_q_spec = c("est_all", "est_all", "est_all")
)

# ll survey, share delta female (index 2) across time blocks and to the coop jp ll survey delta
# ll survey, share delta male (index 5) across time blocks and to the coop jp ll survey delta
# coop jp survey does not estimate parameters and shares deltas with longline survey
# single time block with trawl survey and only one parameter hence, only one parameter estimated across blocks (indices 7 and 8)
input_list$map$ln_srv_fixed_sel_pars <- factor(c(1:3, 2, 4:6, 5,rep(7,4), rep(8, 4), rep(c(NA,2), 2), rep(c(NA, 5), 2)))

# Coop JP Survey (Logistic) Single time block (these estimates are fixed!)
input_list$par$ln_srv_fixed_sel_pars[1,,,1,3] <- c(0.980660760456, tem_par$coefficients[names(tem_par$coefficients) == "log_delta_srv1_f"])
input_list$par$ln_srv_fixed_sel_pars[1,,,2,3] <- c(1.22224502478, tem_par$coefficients[names(tem_par$coefficients) == "log_delta_srv1_m"])

# Setup tagging stuff
input_list <- Setup_Mod_Tagging(input_list = input_list,
                                UseTagging = 0
)

# set up data weighting stuff
Wt_FishAgeComps <- array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_sexes, input_list$data$n_fish_fleets)) # weights for fishery age comps
Wt_FishAgeComps[1,,1,1] <- 0.826107286513784 # Weight for fixed gear age comps

# Fishery length comps
Wt_FishLenComps <- array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_sexes, input_list$data$n_fish_fleets)) # weights for fishery age comps
Wt_FishLenComps[1,,1,1] <- 4.1837057381917 # Weight for fixed gear len comps females
Wt_FishLenComps[1,,2,1] <- 4.26969350917589 # Weight for fixed gear len comps males
Wt_FishLenComps[1,,1,2] <- 0.316485920691651 # Weight for trawl gear len comps females
Wt_FishLenComps[1,,2,2] <- 0.229396580680981 # Weight for trawl gear len comps males

# survey age comps
Wt_SrvAgeComps <- array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_sexes, input_list$data$n_srv_fleets)) # weights for survey age comps
Wt_SrvAgeComps[1,,1,1] <- 3.79224544725927 # Weight for domestic survey ll gear age comps
Wt_SrvAgeComps[1,,1,3] <- 1.31681114024037 # Weight for coop jp survey ll gear age comps

# Survey length comps
Wt_SrvLenComps <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_sexes, input_list$data$n_srv_fleets)) # weights for survey age comps
Wt_SrvLenComps[1,,1,1] <- 1.43792019016567 # Weight for domestic ll survey len comps females
Wt_SrvLenComps[1,,2,1] <- 1.07053763450712 # Weight for domestic ll survey len comps males
Wt_SrvLenComps[1,,1,2] <- 0.670883273592302 # Weight for domestic trawl survey len comps females
Wt_SrvLenComps[1,,2,2] <- 0.465207132450763 # Weight for domestic trawl survey len comps males
Wt_SrvLenComps[1,,1,3] <- 1.27772810174693 # Weight for coop jp ll survey len comps females
Wt_SrvLenComps[1,,2,3] <- 0.857519546948587 # Weight for coop jp ll survey len comps males

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

# Fit Model ---------------------------------------------------------------
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
                              silent = F)

# Save model results
sabie_rtmb_model$sd_rep <- RTMB::sdreport(sabie_rtmb_model)
rep <- sabie_rtmb_model$rep
sd_rep <- sabie_rtmb_model$sd_rep

saveRDS(sd_rep, here("dev", "dev_output", "1_Region_Model_Sablefish", "sd_rep.RDS"))
saveRDS(rep, here("dev", "dev_output", "1_Region_Model_Sablefish", "rep.RDS"))

# Compare to ADMB (23.5b) -------------------------------------------------
# Read in data
tem_dat <- dget(here("dev", "dev_data", '2024 Base (23.5)_final model_v3', 'tem.rdat'))
ageing_dat <- dget(here("dev", "dev_data",'2023 Base (23.5)_final model', 'test.rdat')) # for getting ageing error
tem_admb_dat <- readLines(here("dev", "dev_data","2024 Base (23.5)_final model_v3", "tem_2024_na_wh.dat"))
tem_rep <- readLines(here("dev", "dev_data",'2024 Base (23.5)_final model_v3', 'sable.rep'))
tem_par <- read_pars(here("dev", "dev_data",'2024 Base (23.5)_final model_v3', 'tem'))

# Get population dynamics -------------------------------------------------
ages <- 1:30 # ages
years <- rownames(tem_dat$t.series) # years
ssb <- tem_dat$t.series$spbiom # spawning biomass
total_biom <- tem_dat$t.series$totbiom # total biomass
waa_f <- tem_dat$growthmat$wt.f.block1 # female WAA
waa_m <- tem_dat$growthmat$wt.m.block1 # male WAA
mat <- tem_dat$growthmat$mage.block2 # maturity

### Recruitment -------------------------------------------------------------
rec <- tem_dat$t.series$Recr # total recruitment
init_age <- matrix(cbind(tem_dat$natage.female[1,-1] , tem_dat$natage.male[1,-1] ), ncol = 2) # initial ages
mean_rec <- tem_par$coefficients[names(tem_par$coefficients) == "log_mean_rec"] # mean recruitment
init_rec_devs <- rev(tem_par$coefficients[str_detect(names(tem_par$coefficients), "rec_dev") ][1:28]) # initial devs
rec_devs <- tem_par$coefficients[str_detect(names(tem_par$coefficients), "rec_dev")][-c(1:28)] # rec devs

### Fishing Mortality -------------------------------------------------------
total_f <- tem_dat$t.series$fmort # total F
mean_ll_fish <- tem_par$coefficients[names(tem_par$coefficients) == "log_avg_F_fish1"] # Average F LL
devs_ll_fish <- tem_par$coefficients[str_detect(names(tem_par$coefficients), "log_F_devs_fish1")] # Devs LL
mean_ll_trwl <- tem_par$coefficients[names(tem_par$coefficients) == "log_avg_F_fish3"] # Average F Trawl
devs_ll_trwl <- tem_par$coefficients[str_detect(names(tem_par$coefficients), "log_F_devs_fish3")] # Devs Trawl

# Get time series
rec_series <- data.frame(Par = "Recruitment",
                         Year = 1960:2024,
                         TMB = t(sabie_rtmb_model$rep$Rec),
                         ADMB = rec[1:length(data$years)])

f_series <- data.frame(Par = "Total F",
                       Year = 1960:2024,
                       TMB = apply(sabie_rtmb_model$rep$Fmort,2,sum),
                       ADMB = tem_dat$t.series$fmort[1:length(data$years)])

females_series <- data.frame(Par = "Total Females",
                             Year = 1960:2024,
                             TMB = rowSums(sabie_rtmb_model$rep$NAA[1,-66,,1]),
                             ADMB = tem_dat$t.series$numbers.f[1:length(data$years)])

males_series <- data.frame(Par = "Total Males",
                           Year = 1960:2024,
                           TMB = rowSums(sabie_rtmb_model$rep$NAA[1,-66,,2]),
                           ADMB = tem_dat$t.series$numbers.m[1:length(data$years)])

ssb_series <- data.frame(Par = "SSB",
                         Year = 1960:2024,
                         TMB = t(sabie_rtmb_model$rep$SSB),
                         ADMB = ssb)

ssb_se_series <- data.frame(Par = "SSB (SE)",
                            Year = 1960:2024,
                            TMB = sabie_rtmb_model$sd_rep$sd[names(sabie_rtmb_model$sd_rep$value) == "log(SSB)"] * t(sabie_rtmb_model$rep$SSB),
                            ADMB = tem_par$se[str_detect(names(tem_par$se), "ssb")])

rec_se_series <- data.frame(Par = "Recruitment (SE)",
                            Year = 1960:2024,
                            TMB = sabie_rtmb_model$sd_rep$sd[names(sabie_rtmb_model$sd_rep$value) == "log(Rec)"] * t(sabie_rtmb_model$rep$Rec),
                            ADMB = tem_par$se[str_detect(names(tem_par$se), "pred_rec")])

opt_ts_df_v3 <- rbind(ssb_se_series,rec_series, f_series, females_series, males_series, ssb_series, rec_se_series)

# Get selectivities
dom_ll_fish_f1 <- data.frame(Age = 1:30,
                             TMB = sabie_rtmb_model$rep$fish_sel[1,1,,1,1],
                             ADMB = tem_dat$agesel$fish1sel.f,
                             Type = "Domestic LL Fishery Female Block 1")

dom_ll_fish_m1 <- data.frame(Age = 1:30,
                             TMB = sabie_rtmb_model$rep$fish_sel[1,1,,2,1],
                             ADMB = tem_dat$agesel$fish1sel.m,
                             Type = "Domestic LL Fishery Male Block 1")

dom_ll_fish_f2 <- data.frame(Age = 1:30,
                             TMB = sabie_rtmb_model$rep$fish_sel[1,40,,1,1],
                             ADMB = tem_dat$agesel$fish4sel.f,
                             Type = "Domestic LL Fishery Female Block 2")

dom_ll_fish_m2 <- data.frame(Age = 1:30,
                             TMB = sabie_rtmb_model$rep$fish_sel[1,40,,2,1],
                             ADMB = tem_dat$agesel$fish4sel.m,
                             Type = "Domestic LL Fishery Male Block 2")

dom_ll_fish_f3 <- data.frame(Age = 1:30,
                             TMB = sabie_rtmb_model$rep$fish_sel[1,60,,1,1],
                             ADMB = tem_dat$agesel$fish5sel.f,
                             Type = "Domestic LL Fishery Female Block 3")

dom_ll_fish_m3 <- data.frame(Age = 1:30,
                             TMB = sabie_rtmb_model$rep$fish_sel[1,60,,2,1],
                             ADMB = tem_dat$agesel$fish5sel.m,
                             Type = "Domestic LL Fishery Male Block 3")

dom_trwl_fish_f <- data.frame(Age = 1:30,
                              TMB = sabie_rtmb_model$rep$fish_sel[1,1,,1,2],
                              ADMB = tem_dat$agesel$fish3sel.f,
                              Type = "Domestic Trawl Female")

dom_trwl_fish_m <- data.frame(Age = 1:30,
                              TMB = sabie_rtmb_model$rep$fish_sel[1,1,,2,2],
                              ADMB = tem_dat$agesel$fish3sel.m,
                              Type = "Domestic Trawl Male")

dom_ll_srv_f1 <- data.frame(Age = 1:30,
                            TMB = sabie_rtmb_model$rep$srv_sel[1,1,,1,1],
                            ADMB = tem_dat$agesel$srv1sel.f,
                            Type = "Domestic LL Survey Female Block 1")

dom_ll_srv_m1 <- data.frame(Age = 1:30,
                            TMB = sabie_rtmb_model$rep$srv_sel[1,1,,2,1],
                            ADMB = tem_dat$agesel$srv1sel.m,
                            Type = "Domestic LL Survey Male Block 1")

dom_ll_srv_f2 <- data.frame(Age = 1:30,
                            TMB = sabie_rtmb_model$rep$srv_sel[1,60,,1,1],
                            ADMB = tem_dat$agesel$srv10sel.f,
                            Type = "Domestic LL Survey Female Block 2")

dom_ll_srv_m2 <- data.frame(Age = 1:30,
                            TMB = sabie_rtmb_model$rep$srv_sel[1,60,,2,1],
                            ADMB = tem_dat$agesel$srv10sel.m,
                            Type = "Domestic LL Survey Male Block 2")

dom_trwl_srv_f2 <- data.frame(Age = 1:30,
                              TMB = sabie_rtmb_model$rep$srv_sel[1,60,,1,2],
                              ADMB = tem_dat$agesel$srv7sel.f,
                              Type = "Domestic Trawl Survey Female")

dom_trwl_srv_m2 <- data.frame(Age = 1:30,
                              TMB = sabie_rtmb_model$rep$srv_sel[1,60,,2,2],
                              ADMB = tem_dat$agesel$srv7sel.m,
                              Type = "Domestic Trawl Survey Male")

coop_ll_srv_f2 <- data.frame(Age = 1:30,
                             TMB = sabie_rtmb_model$rep$srv_sel[1,60,,1,3],
                             ADMB = tem_dat$agesel$srv2sel.f,
                             Type = "Coop LL Survey Female")

coop_ll_srv_m2 <- data.frame(Age = 1:30,
                             TMB = sabie_rtmb_model$rep$srv_sel[1,60,,2,3],
                             ADMB = tem_dat$agesel$srv2sel.m,
                             Type = "Coop LL Survey Male")

opt_combined_sel_v3 <- rbind(
  dom_ll_fish_m1,
  dom_ll_fish_f2,
  dom_ll_fish_m2,
  dom_ll_fish_f3,
  dom_ll_fish_m3,
  dom_trwl_fish_f,
  dom_trwl_fish_m,
  dom_ll_srv_f1,
  dom_ll_srv_m1,
  dom_ll_srv_f2,
  dom_ll_srv_m2,
  dom_trwl_srv_f2,
  dom_trwl_srv_m2,
  coop_ll_srv_f2,
  coop_ll_srv_m2
)

# Get parameters
R0_df <- data.frame(Par = "Mean Recruitment",
                    TMB = sabie_rtmb_model$rep$R0,
                    ADMB = exp(parameters$ln_global_R0))

srv_q_df <- data.frame(Par = c('srv_domLL_q', 'srv_trwl_q', 'srv_coopLL_q'),
                       TMB = exp(sabie_rtmb_model$sd_rep$par.fixed[names(sabie_rtmb_model$sd_rep$par.fixed) == "ln_srv_q"]),
                       ADMB = exp(parameters$ln_srv_q[1,,]))

fish_q_df <- data.frame(Par = c('fish_domLL_q1', 'fish_domLL_q2', 'fish_domLL_q3'),
                        TMB = exp(sabie_rtmb_model$sd_rep$par.fixed[names(sabie_rtmb_model$sd_rep$par.fixed) == "ln_fish_q"]),
                        ADMB = exp(parameters$ln_fish_q[,,1]))

par_df_v3 <- rbind(R0_df, srv_q_df, fish_q_df)

# plot!
ggplot(opt_ts_df_v3, aes(x = Year, y = ((TMB - ADMB) / ADMB))) +
  geom_line(size = 2) +
  geom_hline(yintercept = 0, lty = 2, size = 1.3) +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~Par, scales = "free_y") +
  labs(y = "Relative difference (%)")  +
  theme_sablefish()

ggplot() +
  geom_line(opt_ts_df_v3, mapping  = aes(x = Year, y = TMB), size = 1.3, lty = 1) +
  geom_line(opt_ts_df_v3, mapping  = aes(x = Year, y = ADMB), size = 1.3, lty = 2, color = 'red') +
  facet_wrap(~Par, scales = "free_y") +
  labs(x = "Year", color = 'Model', y = "Value") +
  theme_sablefish()

# Vignette Figures
png(here("vignettes", "figures","e_ts_comparison_re.png"), width = 1000, height = 500)
ggplot(opt_ts_df_v3 %>% dplyr::filter(Par %in% c("SSB", "Recruitment")), aes(x = Year, y = ((TMB - ADMB) / ADMB))) +
  geom_line(size = 2) +
  geom_hline(yintercept = 0, lty = 2, size = 1.3) +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~Par, scales = "free_y") +
  labs(y = "Relative difference (%)")  +
  theme_sablefish()
dev.off()

png(here("vignettes", "figures", "e_ts_comparison.png"), width = 1000, height = 500)
ggplot() +
  geom_line(opt_ts_df_v3 %>% dplyr::filter(Par %in% c("SSB", "Recruitment")), mapping  = aes(x = Year, y = TMB), size = 1.3, lty = 1) +
  geom_line(opt_ts_df_v3 %>% dplyr::filter(Par %in% c("SSB", "Recruitment")), mapping  = aes(x = Year, y = ADMB), size = 1.3, lty = 2, color = 'red') +
  facet_wrap(~Par, scales = "free_y") +
  labs(x = "Year", color = 'Model', y = "Value") +
  theme_sablefish()
dev.off()

