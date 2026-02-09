# Purpose: To bridge to the GOA Dusky ADMB model
# Creator: Matthew LH. Cheng
# Date: 6/30/25

# Setup -------------------------------------------------------------------

library(here)
library(SPoRC)
library(tidyverse)
library(PBSmodelling)
devtools::load_all(here("R"))

# Read in data file
in_dat <- readLines(here::here("dev", "dev_data", "goa_dusk_2024.dat"))
# Read in report file
rep_dat <- readLines(here::here("dev", "dev_data", "dusky.rep"))
# Read in parameters and stuff for comparison
par_dat <- R2admb::read_pars(here::here("dev", "dev_data", "dusky"))
rep_dat_rtem <- readList(here::here("dev", "dev_data", "dusky_rtem.rep"))

# Extract Data ------------------------------------------------------------
# Sorry for hard coding line numbers ...

# Number of regions
n_regions <- 1

# Number of sexes
n_sexes <- 1

# Number of fishery fleets
n_fish_fleets <- 1

# Number of survey fleets
n_srv_fleets <- 1

# Get model start and end year
years <- in_dat[[18]]:in_dat[[20]]
n_years <- length(years)

# Get model ages (not the same as data ages)
obs_ages <- in_dat[[22]]:(in_dat[[26]])
mod_ages <- c(obs_ages, 31:33) # add more ages because of model ages
n_ages <- length(mod_ages)
n_obs_ages <- length(obs_ages)

# Get length bins
lens <- as.numeric(unlist(strsplit(in_dat[[34]], " ")))
n_lens <- length(lens)

# Get spawning month
spwn_month <- 0

# Get weight at age
waa <- as.numeric(unlist(strsplit(in_dat[[41]], " ")))

# get size age transition
start_line <- grep("Size-age transition matrix", in_dat, ignore.case = TRUE)
data_lines <- in_dat[(start_line + 1):277]
data_lines <- grep("^[ 0-9eE.+\\-]", data_lines, value = TRUE)
flat_values <- scan(text = paste(data_lines, collapse = "\n"), quiet = TRUE)
ncol <- n_lens
nrow <- length(flat_values) / ncol
size_age_matrix <- matrix(flat_values, ncol = ncol, byrow = TRUE)
size_age_matrix[which(size_age_matrix == -2e-4)] <- 2e-4

# get ageing error
# Find where the age error matrix starts
start_line <- grep("age error transition matrix", in_dat, ignore.case = TRUE)
data_lines <- in_dat[(start_line + 1):313]
data_lines <- grep("^[ 0-9eE.+\\-]", data_lines, value = TRUE)
flat_values <- scan(text = paste(data_lines, collapse = "\n"), quiet = TRUE)
ncol <- n_obs_ages
nrow <- length(flat_values) / ncol
age_error_matrix <- matrix(flat_values, ncol = ncol, byrow = TRUE)

# Get maturity
mat_lines <- grep("Maturity", rep_dat)
mat <- as.numeric(unlist(strsplit(rep_dat[mat_lines[1]], " ")))
mat <- mat[!is.na(mat)]

# Fishery Stuff -----------------------------------------------------------
# Get fishery catch
fish_catch_yrs <- as.numeric(unlist(strsplit(in_dat[[47]], " "))[-1])
fish_catch <- as.numeric(unlist(strsplit(in_dat[[48]], " ")))

# Fishery age comps
fish_catch_age_yrs <- as.numeric(unlist(strsplit(in_dat[[101]], " ")))
fish_catch_age_n <- as.numeric(unlist(strsplit(in_dat[[103]], " ")))
fish_catch_age_nhauls <- as.numeric(unlist(strsplit(in_dat[[105]], " ")))

# fishery ages here
# Locate start of oac_fish matrix
start_line <- grep("Observed fishery age compositions.*oac_fish", in_dat, ignore.case = TRUE)
data_lines <- in_dat[(start_line + 1):123]
data_lines <- grep("^[ 0-9eE.+\\-]", data_lines, value = TRUE)
flat_values <- scan(text = paste(data_lines, collapse = "\n"), quiet = TRUE)
ncol <- n_obs_ages
nrow <- length(flat_values) / ncol
oac_fish_age <- matrix(flat_values, ncol = ncol, byrow = TRUE)

# fishery size comps
fish_catch_size_yrs <- as.numeric(unlist(strsplit(in_dat[[166]], " ")))
fish_catch_size_n <- as.numeric(unlist(strsplit(in_dat[[168]], " ")))
fish_catch_size_nhauls <- as.numeric(unlist(strsplit(in_dat[[170]], " ")))

# fishery size comps here
start_line <- grep("Observed fishery size compositions.*proportions at age", in_dat)
data_lines <- in_dat[(start_line + 1):191]
data_lines <- grep("^[ 0-9eE.+\\-]", data_lines, value = TRUE)
flat_values <- scan(text = paste(data_lines, collapse = "\n"), quiet = TRUE)
ncol <- n_lens
nrow <- length(flat_values) / ncol
oac_fish_size <- matrix(flat_values, ncol = ncol, byrow = TRUE)

# Survey Stuff ------------------------------------------------------------
# survey index
bts_srv_yrs <- as.numeric(unlist(strsplit(in_dat[[66]], " ")))
bts_srv_idx <- as.numeric(unlist(strsplit(in_dat[[68]], " ")))
bts_srv_se <- as.numeric(unlist(strsplit(in_dat[[70]], " ")))

# survey age comps
bts_catch_age_yrs <- as.numeric(unlist(strsplit(in_dat[[133]], " ")))
bts_catch_age_n <- as.numeric(unlist(strsplit(in_dat[[135]], " ")))
bts_catch_age_nhauls <- as.numeric(unlist(strsplit(in_dat[[137]], " ")))

# survey ages here
# Locate the start of oac_srv1
start_line <- grep("Observed trawl survey age compositions.*oac_srv1", in_dat)
data_lines <- in_dat[(start_line + 1):156]
data_lines <- grep("^[ 0-9eE.+\\-]", data_lines, value = TRUE)
flat_values <- scan(text = paste(data_lines, collapse = "\n"), quiet = TRUE)
ncol <- n_obs_ages
nrow <- length(flat_values) / ncol
oac_srv1_age <- matrix(flat_values, ncol = ncol, byrow = TRUE)

# survey size comps
bts_catch_size_yrs <- as.numeric(unlist(strsplit(in_dat[[201]], " ")))
bts_catch_size_n <- as.numeric(unlist(strsplit(in_dat[[203]], " ")))
bts_catch_size_nhauls <- as.numeric(unlist(strsplit(in_dat[[205]], " ")))

# survey size comps here
# Locate the start of oac_srv1
start_line <- grep("Observed survey size compositions.*oac_fish", in_dat)
data_lines <- in_dat[(start_line + 1):224]
data_lines <- grep("^[ 0-9eE.+\\-]", data_lines, value = TRUE)
flat_values <- scan(text = paste(data_lines, collapse = "\n"), quiet = TRUE)
ncol <- n_lens
nrow <- length(flat_values) / ncol
oac_srv1_size <- matrix(flat_values, ncol = ncol, byrow = TRUE)


# Setup RTMB inputs -------------------------------------------------------

# biologicals
waa_arr <- array(0, dim = c(n_regions, n_years, n_ages, n_sexes))
fix_natmort <- array(0, dim = c(n_regions, n_years, n_ages, n_sexes))
mataa_arr <- array(0, dim = c(n_regions, n_years, n_ages, n_sexes))
sizeage <- array(0, dim = c(n_regions, n_years, n_lens, n_ages, n_sexes))

for(r in 1:n_regions) {
  for(y in 1:n_years) {
    for(s in 1:n_sexes) {
      waa_arr[r,y,,s] <- waa # weight at age
      fix_natmort[r,y,,s] <- 0.07 # natural mortality
      mataa_arr[r,y,,s] <- mat # maturity at age
      sizeage[r,y,,,s] <- t(size_age_matrix) # size age transition
    } # end s loop
  } # end y loop
} # end r loop

# Catch inputs
ObsCatch <- array(0, dim = c(n_regions, n_years, n_fish_fleets))
UseCatch <- array(1, dim = c(n_regions, n_years, n_fish_fleets))
ObsCatch[] <- fish_catch

# Fishery index inputs -- none used
ObsFishIdx <- array(NA, dim = c(n_regions, n_years, n_fish_fleets))
ObsFishIdx_SE <- array(NA, dim = c(n_regions, n_years, n_fish_fleets))
UseFishIdx <- array(0, dim = c(n_regions, n_years, n_fish_fleets))

# Fishery age comps inputs
ObsFishAgeComps <- array(0, dim = c(n_regions, n_years, n_obs_ages, n_sexes, n_fish_fleets))
UseFishAgeComps <- array(0, dim = c(n_regions, n_years, n_fish_fleets))
ISS_FishAgeComps <- array(0, dim = c(n_regions, n_years, n_sexes, n_fish_fleets))

ObsFishAgeComps[1,which(years %in% fish_catch_age_yrs),,1,1] <- oac_fish_age
UseFishAgeComps[1,which(years %in% fish_catch_age_yrs),1] <- 1
ISS_FishAgeComps[1,which(years %in% fish_catch_age_yrs),1,1] <- sqrt(fish_catch_age_nhauls * fish_catch_age_n)
ISS_FishAgeComps[1,which(years %in% fish_catch_age_yrs),1,1] <- rep_dat_rtem$oac_fish_sample[,1]

# fishery length comps inputs
ObsFishLenComps <- array(0, dim = c(n_regions, n_years, n_lens, n_sexes, n_fish_fleets))
UseFishLenComps <- array(0, dim = c(n_regions, n_years, n_fish_fleets))
ISS_FishLenComps <- array(0, dim = c(n_regions, n_years, n_sexes, n_fish_fleets))

ObsFishLenComps[1,which(years %in% fish_catch_size_yrs),,1,1] <- oac_fish_size
UseFishLenComps[1,which(years %in% fish_catch_size_yrs),1] <- 1
ISS_FishLenComps[1,which(years %in% fish_catch_size_yrs),1,1] <- sqrt(fish_catch_size_nhauls * fish_catch_size_n)
ISS_FishLenComps[1,which(years %in% fish_catch_size_yrs),1,1] <- rep_dat_rtem$osc_fish_sample[,1]

# survey index
ObsSrvIdx <- array(0, dim = c(n_regions, n_years, n_srv_fleets))
ObsSrvIdx_SE <- array(0, dim = c(n_regions, n_years, n_srv_fleets))
UseSrvIdx <- array(0, dim = c(n_regions, n_years, n_srv_fleets))

ObsSrvIdx[1,which(years %in% bts_srv_yrs),1] <- bts_srv_idx

# transform to lognormal sd to be the same as dusky assessment
cv <- bts_srv_se / bts_srv_idx
se <- bts_srv_idx * sqrt(log(1 + cv^2))
ObsSrvIdx_SE[1,which(years %in% bts_srv_yrs),1] <- se
UseSrvIdx[1,which(years %in% bts_srv_yrs),1] <- 1

# survey age comps inputs
ObsSrvAgeComps <- array(0, dim = c(n_regions, n_years, n_obs_ages, n_sexes, n_srv_fleets))
UseSrvAgeComps <- array(0, dim = c(n_regions, n_years, n_srv_fleets))
ISS_SrvAgeComps <- array(0, dim = c(n_regions, n_years, n_sexes, n_srv_fleets))

ObsSrvAgeComps[1,which(years %in% bts_catch_age_yrs),,1,1] <- oac_srv1_age
UseSrvAgeComps[1,which(years %in% bts_catch_age_yrs),1] <- 1
ISS_SrvAgeComps[1,which(years %in% bts_catch_age_yrs),1,1] <- sqrt(bts_catch_age_nhauls * bts_catch_age_n)
ISS_SrvAgeComps[1,which(years %in% bts_catch_age_yrs),1,1] <- rep_dat_rtem$oac_srv1_sample[,1]

# survey length comps inputs
ObsSrvLenComps <- array(0, dim = c(n_regions, n_years, n_lens, n_sexes, n_srv_fleets))
UseSrvLenComps <- array(0, dim = c(n_regions, n_years, n_srv_fleets))
ISS_SrvLenComps <- array(0, dim = c(n_regions, n_years, n_sexes, n_srv_fleets))

ObsSrvLenComps[1,which(years %in% bts_catch_size_yrs),,1,1] <- oac_srv1_size
UseSrvLenComps[1,which(years %in% bts_catch_size_yrs),1] <- 0
ISS_SrvLenComps[1,which(years %in% bts_catch_size_yrs),1,1] <- 0

# Setup Model -------------------------------------------------------------

input_list <- Setup_Mod_Dim(
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
  ln_sigmaR = rep(-0.1068576 , 2), # 2 values for early and late sigma
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
  ln_sigmaC = array(log(sqrt(1/2)), dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets)),
  ln_sigmaF = array(log(sqrt(1/2)), dim = c(input_list$data$n_regions, input_list$data$n_fish_fleets))
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
                           input_list$data$n_srv_fleets))
)

# catch weigthing for duskies
Wt_Catch <- array(0, dim = c(sgl_rg_dusky_data$n_regions, length(sgl_rg_dusky_data$years), sgl_rg_dusky_data$n_fish_fleets))
Wt_Catch[,which(sgl_rg_dusky_data$years %in% 1977:1991),] <- 2
Wt_Catch[,-which(sgl_rg_dusky_data$years %in% 1977:1991),] <- 50

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
                                     input_list$data$n_sexes,
                                     input_list$data$n_fish_fleets)),
  Wt_FishLenComps = array(1, dim = c(input_list$data$n_regions,
                                     length(input_list$data$years),
                                     input_list$data$n_sexes,
                                     input_list$data$n_fish_fleets)),
  Wt_SrvAgeComps = array(1, dim = c(input_list$data$n_regions,
                                    length(input_list$data$years),
                                    input_list$data$n_sexes,
                                    input_list$data$n_srv_fleets)),
  Wt_SrvLenComps = array(0, dim = c(input_list$data$n_regions,
                                    length(input_list$data$years),
                                    input_list$data$n_sexes,
                                    input_list$data$n_srv_fleets))
)

data <- input_list$data
parameters <- input_list$par
mapping <- input_list$map

# Fix Values --------------------------------------------------------------

# Extract parameter values and set starting values
parameters$ln_global_R0 <- par_dat$coefficients[names(par_dat$coefficients) == 'log_mean_rec']
parameters$ln_sigmaR[] <- log(par_dat$coefficients[names(par_dat$coefficients) == 'sigr'])
parameters$ln_F_mean[] <- par_dat$coefficients[names(par_dat$coefficients) == 'log_avg_F']
parameters$ln_F_devs[] <- par_dat$coefficients[str_detect(names(par_dat$coefficients), "log_F_devs")]
parameters$ln_InitDevs[,1:(n_ages - 2)] <- rev(par_dat$coefficients[str_detect(names(par_dat$coefficients), "log_rec_dev")][1:(n_ages - 2)])
parameters$ln_RecDevs[] <- par_dat$coefficients[str_detect(names(par_dat$coefficients), "log_rec_dev")][-c(1:(n_ages - 2))]
parameters$ln_srv_q[] <- par_dat$coefficients[str_detect(names(par_dat$coefficients), "q_srv1")]
parameters$ln_fish_fixed_sel_pars[] <- log(c(par_dat$coefficients[str_detect(names(par_dat$coefficients), "a50")][1], par_dat$coefficients[str_detect(names(par_dat$coefficients), "delta")][1]))
parameters$ln_srv_fixed_sel_pars[] <- log(c(par_dat$coefficients[str_detect(names(par_dat$coefficients), "a50")][3], par_dat$coefficients[str_detect(names(par_dat$coefficients), "delta")][3]))

# Note: We are estimating 131 parameters, the dusky assessment estimates 137
# Disrepancies come from maturity parameters (2), the reference points (3), and the sigmaR (1) which is fixed in this case.

# Check unoptimized values to see if we get back the same values
obj <- RTMB::MakeADFun(SPoRC:::cmb(SPoRC:::SPoRC_rtmb, data), parameters = parameters, map = mapping, random = NULL, silent = F)
obj$rep <- obj$report(obj$env$last.par.best)

plot(obj$rep$NAA[1,1,,1], type = 'l')
lines(rep_dat_rtem$natage[1,-1])

# Compare unoptimized likelihoods
# Catch
sum(data$Wt_Catch * obj$rep$Catch_nLL) # RTMB
rep_dat_rtem$like[1,1]

# survey
sum(data$Wt_SrvIdx * obj$rep$SrvIdx_nLL) # RTMB
rep_dat_rtem$like[3,1]

# Note: Survey index is a bit off because of change in likelihood formulation, as well as calculating survey index midyear

# Fishery ages
sum(obj$rep$FishAgeComps_nLL) # RTMB
rep_dat_rtem$like[5,1]

# Survey ages
sum(obj$rep$SrvAgeComps_nLL) # RTMB
rep_dat_rtem$like[6,1]

# Fishery Sizes
sum(obj$rep$FishLenComps_nLL) # RTMB
rep_dat_rtem$like[7,1]

# Recruitment penalty
sum(obj$rep$Init_Rec_nLL, obj$rep$Rec_nLL) # RTMB
rep_dat_rtem$pen[1,1]

# Fishing mortality penalty
data$Wt_F * sum(obj$rep$Fmort_nLL) # RTMB
rep_dat_rtem$pen[9,1]

# Priors
obj$rep$srv_q_nLL # RTMB
rep_dat_rtem$prior[2,]

# Compare deterministic catch
plot(as.vector(obj$rep$PredCatch), col = 'red') # RTMB
lines((rep_dat_rtem$tseries[3,]))

# Compare deterministic ssb
plot(as.vector(obj$rep$SSB), col = 'red') # RTMB
lines((rep_dat_rtem$tseries[6,]))

# Look at relative difference (minimal)
(as.vector(obj$rep$SSB) - (rep_dat_rtem$tseries[6,])) / (rep_dat_rtem$tseries[6,]) * 100

# Compare deterministic total biomass
plot(as.vector(obj$rep$Total_Biom), col = 'red') # RTMB
lines((rep_dat_rtem$tseries[5,]))

# Look at relative difference (minimal)
(as.vector(obj$rep$Total_Biom) - (rep_dat_rtem$tseries[5,])) / (rep_dat_rtem$tseries[5,]) * 100

# Any differences that remain aren't generally due to population dynamics being different,
# rather, they are due to differences in likelihood formulations

# Optimize Values ---------------------------------------------------------

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

# Compare optimized likelihoods
# Catch
sum(data$Wt_Catch * dusky_rtmb_model$rep$Catch_nLL) # RTMB
rep_dat_rtem$like[1,1]

# survey
sum(data$Wt_SrvIdx * dusky_rtmb_model$rep$SrvIdx_nLL) # RTMB
rep_dat_rtem$like[3,1]

# Note: Survey index is a bit off because of change in likelihood formulation, as well as calculating survey index midyear

# Fishery ages
sum(dusky_rtmb_model$rep$FishAgeComps_nLL) # RTMB
rep_dat_rtem$like[5,1]

# Survey ages
sum(dusky_rtmb_model$rep$SrvAgeComps_nLL) # RTMB
rep_dat_rtem$like[6,1]

# Fishery Sizes
sum(dusky_rtmb_model$rep$FishLenComps_nLL) # RTMB
rep_dat_rtem$like[7,1]

# Recruitment penalty
sum(dusky_rtmb_model$rep$Init_Rec_nLL, dusky_rtmb_model$rep$Rec_nLL) # RTMB
rep_dat_rtem$pen[1,1]

# Fishing mortality penalty
data$Wt_F * sum(dusky_rtmb_model$rep$Fmort_nLL) # RTMB
rep_dat_rtem$pen[9,1]

# Priors
dusky_rtmb_model$rep$srv_q_nLL # RTMB
rep_dat_rtem$prior[2,]

# Plots -------------------------------------------------------------------
get_data_fitted_plot(list(data), "2024_Dusky")
get_biological_plot(list(data), list(dusky_rtmb_model$rep), "2024_Dusky")
get_selex_plot(list(dusky_rtmb_model$rep), "2024_Dusky")
get_ts_plot(list(dusky_rtmb_model$rep), list(dusky_rtmb_model$sdrep), "2024_Dusky")
get_nLL_plot(list(data), list(dusky_rtmb_model$rep), "2024_Dusky")


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
  # coord_cartesian(ylim = c(-0.05, 0.05)) +
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
