# Purpose: Build sgl_rg_goa_nork_data, the inputs, ADMB seed, and ADMB output for the GOA northern rockfish case study
# Creator: Matthew LH. Cheng
# Date Created: 8/7/26
#
# Sources, both under dev/dev_data:
#   goa_northern_dat.RDS      assembled inputs from the 2024 assessment
#   goa_northern_model.rds    ADMB maximum likelihood estimates and outputs,
#                             the comparison target

library(here)

dat <- readRDS(here("dev", "dev_data", "goa_northern_dat.RDS"))
m24 <- readRDS(here("dev", "dev_data", "goa_northern_model.rds"))

# Dimensions -------------------------------------------------------------------
# 50 model ages (2 to 51) against 44 observed (2 to 45): the ageing error matrix maps between them
n_pop <- 1
n_regions <- 1
n_seas <- 1
n_sexes <- 1
n_fish_fleets <- 1
n_srv_fleets <- 1

years <- dat$years
ages <- 2:51
obs_ages <- dat$ages
lens <- dat$length_bins

n_yrs <- length(years)
n_ages <- length(ages)
n_obs_ages <- length(obs_ages)
n_lens <- length(lens)

stopifnot(length(dat$waa) == n_ages, n_obs_ages == 44, n_lens == 31)

# Biologicals ------------------------------------------------------------------
WAA <- array(NA_real_, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes))
WAA[1, 1, , 1, , 1] <- matrix(rep(dat$waa, each = n_yrs), nrow = n_yrs)

# maturity is estimated inside the ADMB template from the binomial maturity data rather than read
# from the inputs, so the fitted ogive is reused and fixed in SPoRC
MatAA <- array(NA_real_, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes))
MatAA[1, 1, , 1, , 1] <- matrix(rep(m24$maa, each = n_yrs), nrow = n_yrs)

# size_age is stored as [age x length]; SPoRC expects [length x age].
SizeAgeTrans <- array(NA_real_, dim = c(n_pop, n_regions, n_yrs, n_seas, n_lens, n_ages, n_sexes))
for(y in seq_len(n_yrs)) SizeAgeTrans[1, 1, y, 1, , , 1] <- t(dat$size_age)

# age_error is [true age x observed age], which is the orientation SPoRC uses.
AgeingError <- array(NA_real_, dim = c(n_yrs, n_ages, n_obs_ages))
for(y in seq_len(n_yrs)) AgeingError[y, , ] <- dat$age_error

# Catch ------------------------------------------------------------------------
# the reconstructed early catches (1961 to 1977) carry a sum of squares weight of 5 against the
# modern series' 50, kept to SPoRC as fixed lognormal catch standard deviations
ObsCatch <- array(dat$catch_obs, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))
UseCatch <- array(dat$catch_ind, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))
ln_sigmaC <- array(NA_real_, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))
ln_sigmaC[1, , 1, 1] <- log(sqrt(1 / (2 * dat$catch_wt)))

# Survey index -----------------------------------------------------------------
# the GOA bottom trawl survey is the only index. its reported standard deviation is arithmetic
# scale, converted to the exact lognormal sigma
srv_ind <- dat$srv_ind
ObsSrvIdx <- array(NA_real_, dim = c(n_regions, n_yrs, n_seas, n_srv_fleets))
ObsSrvIdx_SE <- array(NA_real_, dim = c(n_regions, n_yrs, n_seas, n_srv_fleets))
UseSrvIdx <- array(0L, dim = c(n_regions, n_yrs, n_seas, n_srv_fleets))
ObsSrvIdx[1, srv_ind == 1, 1, 1] <- dat$srv_obs
ObsSrvIdx_SE[1, srv_ind == 1, 1, 1] <- sqrt(log(1 + dat$srv_sd^2 / dat$srv_obs^2))
UseSrvIdx[1, srv_ind == 1, 1, 1] <- 1L

# Composition helpers ----------------------------------------------------------
make_obs_comp <- function(obs_mat, ind_vec, n_bins, n_fl) {
  arr <- array(NA_real_, dim = c(n_regions, n_yrs, n_seas, n_bins, n_sexes, n_fl))
  iy <- 1L
  for(t in seq_len(n_yrs)) {
    if(ind_vec[t] == 1L) {
      arr[1, t, 1, , 1, 1] <- obs_mat[, iy]
      iy <- iy + 1L
    }
  } # end t loop
  arr
}
make_use_comp <- function(ind_vec, n_fl) {
  arr <- array(0L, dim = c(n_regions, n_yrs, n_seas, n_fl))
  arr[1, , 1, 1] <- as.integer(ind_vec)
  arr
}
make_iss_comp <- function(iss_vec, ind_vec, n_fl) {
  arr <- array(0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_fl))
  iy <- 1L
  for(t in seq_len(n_yrs)) {
    if(ind_vec[t] == 1L) {
      arr[1, t, 1, 1, 1] <- iss_vec[iy]
      iy <- iy + 1L
    }
  } # end t loop
  arr
}

# Compositions -----------------------------------------------------------------
# Every composition source holds the assessment's fixed weight of 0.5, which
# is applied over the whole weight array rather than only the observed years.
ObsFishAgeComps <- make_obs_comp(dat$fish_age_obs, dat$fish_age_ind, n_obs_ages, n_fish_fleets)
UseFishAgeComps <- make_use_comp(dat$fish_age_ind, n_fish_fleets)
ISS_FishAgeComps <- make_iss_comp(dat$fish_age_iss, dat$fish_age_ind, n_fish_fleets)
Wt_FishAgeComps <- array(dat$fish_age_wt, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_fish_fleets))

ObsSrvAgeComps <- make_obs_comp(dat$srv_age_obs, dat$srv_age_ind, n_obs_ages, n_srv_fleets)
UseSrvAgeComps <- make_use_comp(dat$srv_age_ind, n_srv_fleets)
ISS_SrvAgeComps <- make_iss_comp(dat$srv_age_iss, dat$srv_age_ind, n_srv_fleets)
Wt_SrvAgeComps <- array(dat$srv_age_wt, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_srv_fleets))

ObsFishLenComps <- make_obs_comp(dat$fish_size_obs, dat$fish_size_ind, n_lens, n_fish_fleets)
UseFishLenComps <- make_use_comp(dat$fish_size_ind, n_fish_fleets)
ISS_FishLenComps <- make_iss_comp(dat$fish_size_iss, dat$fish_size_ind, n_fish_fleets)
Wt_FishLenComps <- array(dat$fish_size_wt, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_fish_fleets))

# The survey contributes no length compositions.
ObsSrvLenComps <- array(NA_real_, dim = c(n_regions, n_yrs, n_seas, n_lens, n_sexes, n_srv_fleets))
UseSrvLenComps <- array(0L, dim = c(n_regions, n_yrs, n_seas, n_srv_fleets))
ISS_SrvLenComps <- array(0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_srv_fleets))
Wt_SrvLenComps <- array(0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_srv_fleets))

# ADMB maximum likelihood estimate ---------------------------------------------
# initial age deviations are back derived from the ADMB first year numbers at age against the
# unfished geometric series, so seeding ln_InitDevs with them reproduces ADMB's start exactly
M <- m24$M
R0 <- exp(m24$log_mean_R)
NAA_equil <- R0 * exp(-(1:(n_ages - 1)) * M)
NAA_equil[n_ages - 1] <- NAA_equil[n_ages - 2] * exp(-M) / (1 - exp(-M))
init_devs <- log(m24$Nat[2:n_ages, 1]) - log(NAA_equil)

mle <- list(
  M = m24$M,
  log_mean_R = m24$log_mean_R,
  log_Rt = m24$log_Rt,
  log_mean_F = m24$log_mean_F,
  log_Ft = m24$log_Ft,
  a50C = m24$a50C,
  deltaC = m24$deltaC,
  a50S = m24$a50S,
  deltaS = m24$deltaS,
  q = m24$q,
  sigmaR = m24$sigmaR,
  init_devs = init_devs
)

# ADMB output ------------------------------------------------------------------
# The comparison target. Nat is [n_ages x n_yrs] over the model ages.
admb <- list(
  NAA = m24$Nat,
  SSB = m24$spawn_bio,
  Rec = m24$recruits,
  Fmort = m24$Ft,
  TotBiom = m24$tot_bio,
  pred_catch = m24$catch_pred,
  pred_srv = m24$srv_pred,
  sel = m24$slx,
  ssqcatch = m24$ssqcatch,
  like_srv = m24$like_srv,
  like_fish_age = m24$like_fish_age,
  like_srv_age = m24$like_srv_age,
  like_fish_size = m24$like_fish_size,
  like_rec = m24$like_rec,
  f_regularity = m24$f_regularity,
  nll_M = m24$nll_M,
  nll_q = m24$nll_q,
  nll = m24$nll,
  B0 = m24$B0,
  B40 = m24$B40,
  B35 = m24$B35,
  F40 = m24$F40,
  F35 = m24$F35
)

sgl_rg_goa_nork_data <- list(
  years = years,
  ages = ages,
  obs_ages = obs_ages,
  lens = lens,
  n_regions = n_regions,
  n_sexes = n_sexes,
  n_fish_fleets = n_fish_fleets,
  n_srv_fleets = n_srv_fleets,
  n_seas = n_seas,
  n_pop = n_pop,
  natal_region = 1,

  WAA = WAA,
  MatAA = MatAA,
  SizeAgeTrans = SizeAgeTrans,
  AgeingError = AgeingError,

  ObsCatch = ObsCatch,
  UseCatch = UseCatch,
  ln_sigmaC = ln_sigmaC,

  ObsSrvIdx = ObsSrvIdx,
  ObsSrvIdx_SE = ObsSrvIdx_SE,
  UseSrvIdx = UseSrvIdx,

  ObsFishAgeComps = ObsFishAgeComps,
  UseFishAgeComps = UseFishAgeComps,
  ISS_FishAgeComps = ISS_FishAgeComps,
  Wt_FishAgeComps = Wt_FishAgeComps,
  ObsFishLenComps = ObsFishLenComps,
  UseFishLenComps = UseFishLenComps,
  ISS_FishLenComps = ISS_FishLenComps,
  Wt_FishLenComps = Wt_FishLenComps,

  ObsSrvAgeComps = ObsSrvAgeComps,
  UseSrvAgeComps = UseSrvAgeComps,
  ISS_SrvAgeComps = ISS_SrvAgeComps,
  Wt_SrvAgeComps = Wt_SrvAgeComps,
  ObsSrvLenComps = ObsSrvLenComps,
  UseSrvLenComps = UseSrvLenComps,
  ISS_SrvLenComps = ISS_SrvLenComps,
  Wt_SrvLenComps = Wt_SrvLenComps,

  yrs_srv = dat$srv_yrs,
  yrs_fish_age = dat$fish_age_yrs,
  yrs_fish_len = dat$fish_size_yrs,
  yrs_srv_age = dat$srv_age_yrs,

  sigmaR = m24$sigmaR,
  mean_M = dat$mean_M,
  cv_M = dat$cv_M,
  mean_q = dat$mean_q,
  cv_q = dat$cv_q,
  catch_wt = dat$catch_wt,
  srv_wt = dat$srv_wt,
  fmort_wt = dat$wt_fmort_reg,

  mle = mle,
  admb = admb
)

usethis::use_data(sgl_rg_goa_nork_data, internal = FALSE, overwrite = TRUE)
