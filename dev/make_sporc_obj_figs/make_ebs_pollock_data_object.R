# Purpose: Build sgl_rg_ebswp_data, the inputs, ADMB seed, and ADMB output for the EBS pollock case study
# Creator: Matthew LH. Cheng
# Date Created: 8/7/26
#
# Sources, all under dev/dev_data:
#   pm_24.dat     2024 assessment input file, read with ebswp::read_dat
#   pm.par        ADMB parameter file
#   pm_base.rds   RTMB reproduction of pm.rep and pm.par, the comparison target
#   cov_2024.dat  bottom trawl survey index covariance matrix

library(here)

in_dat <- ebswp::read_dat(here("dev", "dev_data", "pm_24.dat"))
pm <- readRDS(here("dev", "dev_data", "pm_base.rds"))$report
par_txt <- readLines(here("dev", "dev_data", "pm.par"))
bts_cov <- as.matrix(utils::read.table(here("dev", "dev_data", "cov_2024.dat")))

# Dimensions -----------------------------------------------------------------
# survey fleets: 1 bottom trawl, 2 acoustic trawl, 3 acoustic vessel of opportunity, 4 the
# acoustic survey's age 1 abundance, which the assessment fits as an index in its own right
years <- seq(in_dat$styr, in_dat$endyr)
ages <- 1:in_dat$nages
n_yrs <- length(years)
n_ages <- length(ages)
n_pop <- 1
n_regions <- 1
n_seas <- 1
n_sexes <- 1
n_fish_fleets <- 1
n_srv_fleets <- 4

i_fsh <- which(years %in% in_dat$yrs_fsh_data)
i_bts <- which(years %in% in_dat$yrs_bts_data)
i_ats <- which(years %in% in_dat$yrs_ats_data)
i_avo <- which(years %in% in_dat$yrs_avo)

# Pull a named block out of the ADMB parameter file.
grab <- function(nm) {
  i <- grep(paste0("^# ", nm, ":"), par_txt)
  stopifnot(length(i) == 1)
  j <- i + 1
  v <- numeric(0)
  while(j <= length(par_txt) && !grepl("^#", par_txt[j])) {
    v <- c(v, scan(text = par_txt[j], quiet = TRUE))
    j <- j + 1
  }
  v
}

# Biologicals ----------------------------------------------------------------
# The assessment has a different weight at age matrix for spawning biomass,
# for catch, and for each survey index.
WAA <- array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes))
MatAA <- array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes))
WAA[, 1, , 1, , 1] <- in_dat$wt_ssb[, 1:n_ages]
MatAA[, 1, , 1, , 1] <- rep(in_dat$p_mature[1:n_ages], each = n_yrs)

WAA_fish <- array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_fish_fleets))
WAA_fish[, 1, , 1, , 1, 1] <- in_dat$wt_fsh[, 1:n_ages]

WAA_srv <- array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_srv_fleets))
WAA_srv[, 1, i_bts, 1, , 1, 1] <- in_dat$wt_bts[, 1:n_ages]
WAA_srv[, 1, i_ats, 1, , 1, 2] <- in_dat$wt_ats[, 1:n_ages]
WAA_srv[, 1, i_avo, 1, , 1, 3] <- in_dat$wt_avo[, 1:n_ages]
WAA_srv[, 1, i_ats, 1, , 1, 4] <- in_dat$wt_ats[, 1:n_ages]

AgeingError <- as.matrix(in_dat$age_err)

# Catch ----------------------------------------------------------------------
ObsCatch <- array(0, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))
ObsCatch[1, , 1, 1] <- in_dat$obs_catch
UseCatch <- array(1, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))

# Fishery index and compositions ---------------------------------------------
# The CPUE series covers the second through thirteenth model years and is fit on
# the arithmetic scale, so its standard error is kept untransformed.
ObsFishIdx <- array(NA_real_, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))
ObsFishIdx_SE <- array(NA_real_, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))
UseFishIdx <- array(0, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))
i_cpue <- which(years %in% in_dat$yrs_cpue)
ObsFishIdx[1, i_cpue, 1, 1] <- in_dat$obs_cpue
ObsFishIdx_SE[1, i_cpue, 1, 1] <- in_dat$obs_cpue_std
UseFishIdx[1, i_cpue, 1, 1] <- 1

ObsFishAgeComps <- array(NA_real_, dim = c(n_regions, n_yrs, n_seas, n_ages, n_sexes, n_fish_fleets))
UseFishAgeComps <- array(0, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))
ISS_FishAgeComps <- array(0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_fish_fleets))
ObsFishAgeComps[1, i_fsh, 1, , 1, 1] <- in_dat$oac_fsh_data[, 1:n_ages]
UseFishAgeComps[1, i_fsh, 1, 1] <- 1
ISS_FishAgeComps[1, i_fsh, 1, 1, 1] <- in_dat$sam_fsh

# Survey indices -------------------------------------------------------------
# trawl fits with a full covariance matrix on the arithmetic scale, acoustic lognormally, vessel
# of opportunity normally, age 1 lognormally with a variance folded into its weight
ObsSrvIdx <- array(NA_real_, dim = c(n_regions, n_yrs, n_seas, n_srv_fleets))
ObsSrvIdx_SE <- array(NA_real_, dim = c(n_regions, n_yrs, n_seas, n_srv_fleets))
UseSrvIdx <- array(0, dim = c(n_regions, n_yrs, n_seas, n_srv_fleets))

ObsSrvIdx[1, i_bts, 1, 1] <- in_dat$ob_bts
ObsSrvIdx_SE[1, i_bts, 1, 1] <- in_dat$ob_bts_std
UseSrvIdx[1, i_bts, 1, 1] <- 1

ObsSrvIdx[1, i_ats, 1, 2] <- in_dat$ob_ats
ObsSrvIdx_SE[1, i_ats, 1, 2] <- sqrt(log((in_dat$ob_ats_std / in_dat$ob_ats)^2 + 1))
UseSrvIdx[1, i_ats, 1, 2] <- 1

ObsSrvIdx[1, i_avo, 1, 3] <- in_dat$ob_avo
ObsSrvIdx_SE[1, i_avo, 1, 3] <- in_dat$ob_avo_std
UseSrvIdx[1, i_avo, 1, 3] <- 1

ObsSrvIdx[1, i_ats, 1, 4] <- in_dat$oac_ats[, 1]
ObsSrvIdx_SE[1, i_ats, 1, 4] <- 1
UseSrvIdx[1, i_ats, 1, 4] <- 1

# Survey compositions --------------------------------------------------------
ObsSrvAgeComps <- array(NA_real_, dim = c(n_regions, n_yrs, n_seas, n_ages, n_sexes, n_srv_fleets))
UseSrvAgeComps <- array(0, dim = c(n_regions, n_yrs, n_seas, n_srv_fleets))
ISS_SrvAgeComps <- array(0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_srv_fleets))

ObsSrvAgeComps[1, i_bts, 1, , 1, 1] <- in_dat$oac_bts[, 1:n_ages]
UseSrvAgeComps[1, i_bts, 1, 1] <- 1
ObsSrvAgeComps[1, i_ats, 1, , 1, 2] <- in_dat$oac_ats[, 1:n_ages]
UseSrvAgeComps[1, i_ats, 1, 2] <- 1

# pm.tpl divides the observed trawl composition by its ages 2-15 subtotal while keeping all 15
# entries, which is a per-year inflation of the sample size applied to a normalized composition
bts_c <- rowSums(in_dat$oac_bts[, 1:n_ages]) / rowSums(in_dat$oac_bts[, 2:n_ages])
ISS_SrvAgeComps[1, i_bts, 1, 1, 1] <- floor(in_dat$sam_bts) * bts_c
ISS_SrvAgeComps[1, i_ats, 1, 1, 2] <- floor(in_dat$sam_ats)

# ADMB maximum likelihood estimate -------------------------------------------
# F is reported at age, so it is divided back out by selectivity at the age
# where selectivity peaks to recover the annual scalar.
a_ref <- apply(pm$sel_fsh, 1, which.max)
Fmort <- sapply(1:n_yrs, function(y) pm$F[y, a_ref[y]] / pm$sel_fsh[y, a_ref[y]])

# Fishery selectivity: log scale coefficients over ages 1-12 with a flat tail,
# and deviations as the running sum of pm.tpl's increments at its change years.
yrs_ch_fsh <- 1965:2023
coffs_f <- grab("sel_coffs_fsh")
nsel_f <- length(coffs_f)
inc_f <- matrix(grab("sel_devs_fsh"), nrow = length(yrs_ch_fsh), ncol = nsel_f, byrow = TRUE)
cum_f <- matrix(0, n_yrs, nsel_f)
for(i in 1:(n_yrs - 1)) {
  k <- match(years[i], yrs_ch_fsh)
  cum_f[i + 1, ] <- cum_f[i, ] + if(!is.na(k)) inc_f[k, ] else 0
}
pars_fsh <- c(coffs_f, rep(coffs_f[nsel_f], n_ages - nsel_f))
devs_fsh <- cbind(cum_f, matrix(cum_f[, nsel_f], n_yrs, n_ages - nsel_f))

# acoustic survey: age 1 pinned at log selectivity zero, coefficients over ages 2-8, ages 9-15 flat.
# only within-year differences are identified, so they are read off the reported surface instead
yrs_ch_ats <- 1995:2024
coffs_a <- grab("sel_coffs_ats")
nsel_a <- 8
pars_ats <- c(0, coffs_a, rep(coffs_a[length(coffs_a)], n_ages - nsel_a))
i_ats_all <- which(years >= in_dat$styr_ats)
devs_ats <- matrix(0, n_yrs, n_ages)
for(k in seq_along(i_ats_all)) devs_ats[i_ats_all[k], ] <- log(pm$sel_ats[k, ]) - pars_ats
for(i in seq_len(min(i_ats_all) - 1)) devs_ats[i, ] <- devs_ats[min(i_ats_all), ]

# bottom trawl survey: logistic with a year-specific midpoint and slope plus a free age 1 value.
# SPoRC evaluates at the raw ages while pm.tpl evaluates at 0.5 + age, so the midpoint shifts by 0.5
slp <- grab("sel_slp_bts_dev")
a50 <- grab("sel_a50_bts_dev")
age1 <- grab("sel_age_one_bts_dev")
i_bts_all <- which(years >= in_dat$styr_bts)
bts_b50_dev <- rep(log(exp(a50[1]) - 0.5), n_yrs)
bts_b50_dev[i_bts_all] <- log(exp(a50) - 0.5)
bts_k_dev <- rep(slp[1], n_yrs)
bts_k_dev[i_bts_all] <- slp
bts_age1_dev <- rep(age1[1], n_yrs)
bts_age1_dev[i_bts_all] <- age1

mle <- list(
  Fmort = Fmort,
  Rec = pm$N[, 1],
  log_initdevs = grab("log_initdevs"),
  steepness = pm$steepness,
  ln_global_R0 = 10.0952990352,
  pars_fsh = pars_fsh,
  devs_fsh = devs_fsh,
  pars_ats = pars_ats,
  devs_ats = devs_ats,
  bts_b50_dev = bts_b50_dev,
  bts_k_dev = bts_k_dev,
  bts_age1_dev = bts_age1_dev
)

# ADMB output, the target the bridge is compared against ----------------------
admb <- list(
  SSB = pm$SSB,
  Rec = pm$N[, 1],
  NAA = pm$N,
  Fmort = Fmort,
  sel_fsh = pm$sel_fsh,
  sel_bts = pm$sel_bts,
  sel_ats = pm$sel_ats,
  eb_ats = pm$eb_ats,
  pred_avo = pm$pred_avo,
  pred_cpue = pm$pred_cpue,
  tot_like = pm$tot_like
)

# Write out data -------------------------------------------------------------
sgl_rg_ebswp_data <- list(
  years = years,
  ages = ages,
  n_seas = n_seas,
  n_regions = n_regions,
  n_sexes = n_sexes,
  n_fish_fleets = n_fish_fleets,
  n_srv_fleets = n_srv_fleets,
  n_pop = n_pop,
  natal_region = 1,

  WAA = WAA,
  WAA_fish = WAA_fish,
  WAA_srv = WAA_srv,
  MatAA = MatAA,
  AgeingError = AgeingError,

  ObsCatch = ObsCatch,
  UseCatch = UseCatch,

  ObsFishIdx = ObsFishIdx,
  ObsFishIdx_SE = ObsFishIdx_SE,
  UseFishIdx = UseFishIdx,
  ObsFishAgeComps = ObsFishAgeComps,
  UseFishAgeComps = UseFishAgeComps,
  ISS_FishAgeComps = ISS_FishAgeComps,

  ObsSrvIdx = ObsSrvIdx,
  ObsSrvIdx_SE = ObsSrvIdx_SE,
  UseSrvIdx = UseSrvIdx,
  ObsSrvAgeComps = ObsSrvAgeComps,
  UseSrvAgeComps = UseSrvAgeComps,
  ISS_SrvAgeComps = ISS_SrvAgeComps,
  SrvIdx_Cov = bts_cov,

  yrs_bts = in_dat$yrs_bts_data,
  yrs_ats = in_dat$yrs_ats_data,
  yrs_avo = in_dat$yrs_avo,
  # The stock recruit residuals run over these years only; 1979 is excluded.
  yrs_srr = setdiff(1978:(max(years) - 2), 1979),
  yrs_sel_ch_fsh = yrs_ch_fsh,
  yrs_sel_ch_ats = yrs_ch_ats,

  mle = mle,
  admb = admb
)

usethis::use_data(sgl_rg_ebswp_data, internal = FALSE, overwrite = TRUE)
