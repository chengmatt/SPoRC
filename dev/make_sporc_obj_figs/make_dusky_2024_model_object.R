# Purpose: To bridge to the GOA Dusky ADMB model
# Creator: Matthew LH. Cheng
# Date: 6/30/25

# Setup -------------------------------------------------------------------

library(here)
library(SPoRC)
library(tidyverse)
library(PBSmodelling)
devtools::load_all(here("R"))
data("sgl_rg_dusky_data")

# for comparisons
rep_dat_rtem <- readList(here::here("dev", "dev_data", "dusky_rtem.rep"))
rep_dat <- readLines(here::here("dev", "dev_data", "dusky.rep"))

# Setup Model -------------------------------------------------------------
input_list <- Setup_Mod_Dim(
  # Number of populations
  n_pop = sgl_rg_dusky_data$n_pop,
  natal_region = sgl_rg_dusky_data$natal_region,
  years = sgl_rg_dusky_data$years,
  # vector of years
  ages = sgl_rg_dusky_data$mod_ages,
  # vector of ages
  lens = sgl_rg_dusky_data$lens,
  # number of lengths
  n_regions = sgl_rg_dusky_data$n_regions,
  # number of regions
  n_sexes = sgl_rg_dusky_data$n_sexes,
  # number of sexes
  n_fish_fleets = sgl_rg_dusky_data$n_fish_fleets,
  # number of fishery fleets
  n_srv_fleets = sgl_rg_dusky_data$n_srv_fleets, # number of survey fleets
  n_seas = sgl_rg_dusky_data$n_seas,
  verbose = TRUE # whether to output messages
)

# Setup recruitment stuff (using defaults for other stuff)
input_list <- Setup_Mod_Rec(
  input_list = input_list,

  # Model options
  # Doing bias ramp, but basically setting it so that no lognormal bias correction happens (as in the dusky model)
  do_rec_bias_ramp = 1,
  bias_year = rep(length(sgl_rg_dusky_data$years), 4),
  # do bias ramp (0 == don't do bias ramp, 1 == do bias ramp)
  sigmaR_switch = 1,
  # when to switch from early to late sigmaR (switch in first year)
  ln_sigmaR = array(-0.1068576, dim = c(2, input_list$data$n_pop, input_list$data$n_regions)),
  # Starting values for early and late sigmaR
  rec_model = "mean_rec",
  sigmaR_spec = "fix",
  # fix early sigmaR and late sigmaR
  # recruitment sex ratio
  init_age_strc = 1, # geometric series to derive initial age structure
  ln_global_R0 = log(2.7), # starting value for mean_rec
  t_spawn = sgl_rg_dusky_data$spwn_month
)

input_list <- Setup_Mod_Biologicals(
  input_list = input_list,

  # Data inputs
  WAA = sgl_rg_dusky_data$waa_arr,
  MatAA = sgl_rg_dusky_data$mataa_arr,

  # Model options
  # fit lengths
  fit_lengths = 1,
  SizeAgeTrans = sgl_rg_dusky_data$sizeage,
  AgeingError = sgl_rg_dusky_data$age_error_matrix,
  M_spec = "fix",
  # fixing natural mortality
  Fixed_natmort = sgl_rg_dusky_data$fix_natmort,
  addtocomp = 0.00001
)

# Setup movement stuff (using defaults for other stuff)
input_list <- Setup_Mod_Movement(
  input_list = input_list,
  use_fixed_movement = 1,
  Fixed_Movement = NA,
  do_recruits_move = 0
)

input_list <- Setup_Mod_Tagging(input_list = input_list, UseTagging = 0)


input_list <- Setup_Mod_Catch_and_F(
  input_list = input_list,

  # Data inputs
  ObsCatch = sgl_rg_dusky_data$ObsCatch,
  UseCatch = sgl_rg_dusky_data$UseCatch,

  # Model options
  Use_F_pen = 1,
  # whether to use f penalty, == 0 don't use, == 1 use
  sigmaC_spec = "fix",

  # Fixing sigma C and F
  ln_sigmaC = array(log(sqrt(1/2)), dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets)),
  ln_sigmaF = array(log(sqrt(1/2)), dim = c(input_list$data$n_regions, input_list$data$n_seas, input_list$data$n_fish_fleets))
)

input_list <- Setup_Mod_FishIdx_and_Comps(
  input_list = input_list,

  # data inputs
  ObsFishIdx = sgl_rg_dusky_data$ObsFishIdx, # fishery index
  ObsFishIdx_SE = sgl_rg_dusky_data$ObsFishIdx_SE, # standard errors
  UseFishIdx = sgl_rg_dusky_data$UseFishIdx, # whether fishery indices are used
  ObsFishAgeComps = sgl_rg_dusky_data$ObsFishAgeComps, # observed fishery ages
  UseFishAgeComps = sgl_rg_dusky_data$UseFishAgeComps, # whether fishery ages are used
  ISS_FishAgeComps = sgl_rg_dusky_data$ISS_FishAgeComps, # input sample size for fishery ages
  ObsFishLenComps = sgl_rg_dusky_data$ObsFishLenComps, # observed fishery lengths
  UseFishLenComps = sgl_rg_dusky_data$UseFishLenComps, # whether fishery lengths are used
  ISS_FishLenComps = sgl_rg_dusky_data$ISS_FishLenComps, # input sample size for fishery lengths

  # Model options
  fish_idx_type = c("none"),
  # indices for fishery
  FishAgeComps_LikeType = c("Multinomial"),
  # age comp likelihoods for fishery fleet
  FishLenComps_LikeType = c("Multinomial"),
  # length comp likelihoods for fishery
  FishAgeComps_Type = c("agg_Year_1-terminal_Fleet_1"),
  # age comp structure for fishery
  FishLenComps_Type = c("agg_Year_1-terminal_Fleet_1")
  # length comp structure for fishery
)

# Setup survey indices and compositions
input_list <- Setup_Mod_SrvIdx_and_Comps(
  input_list = input_list,

  # data inputs
  ObsSrvIdx = sgl_rg_dusky_data$ObsSrvIdx, # observed survey index
  ObsSrvIdx_SE = sgl_rg_dusky_data$ObsSrvIdx_SE / sgl_rg_dusky_data$ObsSrvIdx, # lognormal SD
  UseSrvIdx = sgl_rg_dusky_data$UseSrvIdx, # whether survey indices are used
  ObsSrvAgeComps = sgl_rg_dusky_data$ObsSrvAgeComps, # observed survey ages
  ISS_SrvAgeComps = sgl_rg_dusky_data$ISS_SrvAgeComps, # input sample size for survey ages
  UseSrvAgeComps = sgl_rg_dusky_data$UseSrvAgeComps, # whether survey ages are used
  ObsSrvLenComps = sgl_rg_dusky_data$ObsSrvLenComps, # observed survey lengths
  UseSrvLenComps = sgl_rg_dusky_data$UseSrvLenComps, # whether survey lengths are used
  ISS_SrvLenComps = sgl_rg_dusky_data$ISS_SrvLenComps, # input sample size for survey lengths

  # Model options
  srv_idx_type = c("biom"),
  # abundance and biomass for survey fleet 1
  SrvAgeComps_LikeType = c("Multinomial"),
  # survey age composition likelihood for survey fleet 1
  SrvLenComps_LikeType = c("Multinomial"),
  #  survey length composition likelihood for survey fleet 1
  SrvAgeComps_Type = c(
    "agg_Year_1-terminal_Fleet_1"
  ),
  # survey age comp type

  SrvLenComps_Type = c(
    "agg_Year_1-terminal_Fleet_1"
  )
  # survey length comp type
)


input_list <- Setup_Mod_Fishsel_and_Q(

  input_list = input_list,

  # Model options
  # fishery selectivity, whether continuous time-varying
  cont_tv_fish_sel = c("none_Fleet_1"),
  # fishery selectivity blocks
  fish_sel_blocks = c("none_Fleet_1"),
  # fishery selectivity form
  fish_sel_model = c("logist2_Fleet_1"),
  # fishery catchability blocks
  fish_q_blocks = c("none_Fleet_1"),
  # whether to estiamte all fixed effects for fishery selectivity
  fish_fixed_sel_pars_spec = c("est_all"),
  # whether to estiamte all fixed effects for fishery catchability
  fish_q_spec = c("fix")
)

# Setup survey selectivity and catchability
# Set up prior for survey catchability
srv_q_prior <- data.frame(
  region = 1,
  block = 1,
  fleet = 1,
  mu = 1,
  sd = 0.447213595
)

input_list <- Setup_Mod_Srvsel_and_Q(
  input_list = input_list,

  # Model options
  # survey selectivity, whether continuous time-varying
  cont_tv_srv_sel = c("none_Fleet_1"),
  # survey selectivity blocks
  srv_sel_blocks = c("none_Fleet_1"),
  # survey selectivity form
  srv_sel_model = c("logist2_Fleet_1"),
  # survey catchability blocks
  srv_q_blocks = c("none_Fleet_1"),
  # whether to estiamte all fixed effects for survey selectivity
  srv_fixed_sel_pars_spec = c("est_all"),
  # whether to estiamte all fixed effects for survey catchability
  srv_q_spec = c("est_all"),
  Use_srv_q_prior = 1,
  # Use catchability prior
  srv_q_prior = srv_q_prior,
  # survey timing
  t_srv = array(0, dim = c(input_list$data$n_regions,
                           input_list$data$n_seas,
                           input_list$data$n_srv_fleets))
)

# catch weigthing for duskies
Wt_Catch <- array(0, dim = c(sgl_rg_dusky_data$n_regions, length(sgl_rg_dusky_data$years), sgl_rg_dusky_data$n_seas, sgl_rg_dusky_data$n_fish_fleets))
Wt_Catch[,which(sgl_rg_dusky_data$years %in% 1977:1991),1,] <- 2
Wt_Catch[,-which(sgl_rg_dusky_data$years %in% 1977:1991),1,] <- 50

input_list <- Setup_Mod_Weighting(
  input_list = input_list,
  Wt_Catch = Wt_Catch,
  Wt_FishIdx = 1,
  Wt_SrvIdx = 1.66,
  Wt_Rec = 1,
  Wt_F = 2,
  Wt_Tagging = 0,
  Wt_FishAgeComps = array(1, dim = c(input_list$data$n_regions,
                                     length(input_list$data$years),
                                     input_list$data$n_seas,
                                     input_list$data$n_sexes,
                                     input_list$data$n_fish_fleets)),
  Wt_FishLenComps = array(1, dim = c(input_list$data$n_regions,
                                     length(input_list$data$years),
                                     input_list$data$n_seas,
                                     input_list$data$n_sexes,
                                     input_list$data$n_fish_fleets)),
  Wt_SrvAgeComps = array(1, dim = c(input_list$data$n_regions,
                                    length(input_list$data$years),
                                    input_list$data$n_seas,
                                    input_list$data$n_sexes,
                                    input_list$data$n_srv_fleets)),
  Wt_SrvLenComps = array(0, dim = c(input_list$data$n_regions,
                                    length(input_list$data$years),
                                    input_list$data$n_seas,
                                    input_list$data$n_sexes,
                                    input_list$data$n_srv_fleets))
)

data <- input_list$data
parameters <- input_list$par
mapping <- input_list$map

# Optimize Values ---------------------------------------------------------
# n_ages <- length(input_list$data$ages)
# par_dat <- R2admb::read_pars(here::here("dev", "dev_data", "dusky"))
# parameters$ln_global_R0 <- par_dat$coefficients[names(par_dat$coefficients) == 'log_mean_rec']
# parameters$ln_sigmaR[] <- log(par_dat$coefficients[names(par_dat$coefficients) == 'sigr'])
# parameters$ln_F_mean[] <- par_dat$coefficients[names(par_dat$coefficients) == 'log_avg_F']
# parameters$ln_F_devs[] <- par_dat$coefficients[str_detect(names(par_dat$coefficients), "log_F_devs")]
# parameters$ln_InitDevs[,,1:(n_ages - 2)] <- rev(par_dat$coefficients[str_detect(names(par_dat$coefficients), "log_rec_dev")][1:(n_ages - 2)])
# parameters$ln_RecDevs[] <- par_dat$coefficients[str_detect(names(par_dat$coefficients), "log_rec_dev")][-c(1:(n_ages - 2))]
# parameters$ln_srv_q[] <- par_dat$coefficients[str_detect(names(par_dat$coefficients), "q_srv1")]
# parameters$ln_fish_fixed_sel_pars[] <- log(c(par_dat$coefficients[str_detect(names(par_dat$coefficients), "a50")][1], par_dat$coefficients[str_detect(names(par_dat$coefficients), "delta")][1]))
# parameters$ln_srv_fixed_sel_pars[] <- log(c(par_dat$coefficients[str_detect(names(par_dat$coefficients), "a50")][3], par_dat$coefficients[str_detect(names(par_dat$coefficients), "delta")][3]))
#
# # Fit model
# dusky_rtmb_model <- fit_model(data,
#                               parameters,
#                               mapping,
#                               random = NULL,
#                               newton_loops = 3,
#                               silent = FALSE, do_optim = F
# )
#
# plot(dusky_rtmb_model$rep$NAA[1,1,1,1,,1], type = 'l')
# lines(rep_dat_rtem$natage[1,-1])
#
# plot(as.vector(dusky_rtmb_model$rep$SSB), col = 'red') # RTMB
# lines((rep_dat_rtem$tseries[6,]))
#
# SPoRC_rtmb(parameters, data)

# Fit model
dusky_rtmb_model <- fit_model(data,
                              parameters,
                              mapping,
                              random = NULL,
                              newton_loops = 3,
                              silent = FALSE
)

dusky_rtmb_model$sdrep <- RTMB::sdreport(dusky_rtmb_model) # get standard error report

post_optim_sanity_checks(dusky_rtmb_model$sdrep,
                         dusky_rtmb_model$rep,
                         gradient_tol = 1e-10,
                         se_tol = 5,
                         corr_tol = 0.99)

### Comparisons -------------------------------------------------------------

# extract rtmb and admb values

# catch
catch_comp <- data.frame(
  Year = data$years,
  RTMB = as.vector(dusky_rtmb_model$rep$PredCatch),
  ADMB = rep_dat_rtem$tseries[3,],
  Type = 'Predicted Catch'
)

# ssb
ssb_comp <- data.frame(
  Year = data$years,
  RTMB = as.vector(dusky_rtmb_model$rep$SSB),
  ADMB = rep_dat_rtem$tseries[6,],
  Type = 'SSB'
)

# total biomass
totbiom_comp <- data.frame(
  Year = data$years,
  RTMB = as.vector(dusky_rtmb_model$rep$Total_Biom),
  ADMB = rep_dat_rtem$tseries[5,],
  Type = 'Total Biom'
)

# fishing mortality
f_comp <- data.frame(
  Year = data$years,
  RTMB = as.vector(dusky_rtmb_model$rep$Fmort),
  ADMB = rep_dat_rtem$tseries[4,],
  Type = 'F'
)

# recruitment
rec_comp <- data.frame(
  Year = data$years,
  RTMB = as.vector(dusky_rtmb_model$rep$Rec),
  ADMB = rep_dat_rtem$natage[,2],
  Type = 'Recruitment'
)

# compare selectivity
pred_fishsel <- as.numeric(unlist(strsplit(rep_dat[grep("Fishery_Selectivity", rep_dat)], " ")))
pred_fishsel <- pred_fishsel[!is.na(pred_fishsel)]

fishsel_comp <- data.frame(
  Age = data$ages,
  RTMB = dusky_rtmb_model$rep$fish_sel[1,1,,1,1],
  ADMB = pred_fishsel,
  Type = 'Fishery Selectivity'
)

pred_srvsel <- as.numeric(unlist(strsplit(rep_dat[grep("TWL Survey_Selectivity", rep_dat)], " ")))
pred_srvsel <- pred_srvsel[!is.na(pred_srvsel)]

srvsel_comp <- data.frame(
  Age = data$ages,
  RTMB = dusky_rtmb_model$rep$srv_sel[1,1,,1,1],
  ADMB = pred_srvsel,
  Type = 'Survey Selectivity'
)

# combine
sel_comp <- rbind(fishsel_comp, srvsel_comp)
ts_comp <- rbind(totbiom_comp, ssb_comp, catch_comp, f_comp, rec_comp)

ggplot(ts_comp) +
  geom_line(aes(x = Year, y = RTMB, color = 'RTMB'), lty = 2, lwd = 1.3, alpha = 0.95) +
  geom_line(aes(x = Year, y = ADMB, color = 'ADMB'), lty = 1, lwd = 1.3, alpha = 0.95) +
  facet_wrap(~Type, scales = 'free') +
  theme_sablefish() +
  labs(x = ' Year', y = 'Value', color = 'Model')

ggplot(ts_comp) +
  geom_line(aes(x = Year, y = (RTMB - ADMB) / ADMB), size = 1) +
  geom_hline(yintercept = 0, lty = 2, size = 1.3) +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~Type, scales = "free_y") +
  labs(y = "Relative difference (%)") +
  theme_sablefish()

ggplot(sel_comp) +
  geom_line(aes(x = Age, y = RTMB, color = 'RTMB'), lty = 2, lwd = 1.3, alpha = 0.95) +
  geom_line(aes(x = Age, y = ADMB, color = 'ADMB'), lty = 1, lwd = 1.3, alpha = 0.95) +
  facet_wrap(~Type, scales = 'free') +
  theme_sablefish() +
  labs(x = ' Year', y = 'Value', color = 'Model')

# Output this into a data file
dusky_rtmb_model$data <- data
dusky_rtmb_model$parameters <- parameters
dusky_rtmb_model$mapping <- mapping
usethis::use_data(dusky_rtmb_model, internal = FALSE, overwrite = TRUE)
