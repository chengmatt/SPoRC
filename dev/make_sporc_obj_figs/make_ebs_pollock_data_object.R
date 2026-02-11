library(here)

# Single Region EBS Walleye Pollock ----------------------------------------
# Read in data
in_dat <- ebswp::read_dat(here::here("dev", "dev_data", "pm_24_db.dat"))
rep <- ebswp::read_rep(here::here("dev", "dev_data", "pm.rep"))

# Reformat data list to fit SPoRC
ebswp_SPoRC_data <- list()

# Dimensions ---------------------------------------------------------------
ebswp_SPoRC_data$years <- seq(in_dat$styr, in_dat$endyr)
ebswp_SPoRC_data$ages <- 1:in_dat$nages
ebswp_SPoRC_data$n_seas <- 1
ebswp_SPoRC_data$n_regions <- 1
ebswp_SPoRC_data$n_sexes <- 1
ebswp_SPoRC_data$n_fish_fleets <- 1
ebswp_SPoRC_data$n_srv_fleets <- 3

n_regions <- ebswp_SPoRC_data$n_regions
n_yrs <- length(ebswp_SPoRC_data$years)
n_ages <- length(ebswp_SPoRC_data$ages)
n_seas <- ebswp_SPoRC_data$n_seas
n_sexes <- ebswp_SPoRC_data$n_sexes
n_fish_fleets <- ebswp_SPoRC_data$n_fish_fleets
n_srv_fleets <- ebswp_SPoRC_data$n_srv_fleets
yrs <- in_dat$styr:in_dat$endyr

# Biologicals --------------------------------------------------------------
# WAA: region, year, seas, age, sex
ebswp_SPoRC_data$WAA <- array(0, dim = c(n_regions, n_yrs, n_seas, n_ages, n_sexes))
# MatAA: region, year, seas, age, sex
ebswp_SPoRC_data$MatAA <- array(0, dim = c(n_regions, n_yrs, n_seas, n_ages, n_sexes))
for(y in 1:n_yrs) {
  for(seas in 1:n_seas) {
    ebswp_SPoRC_data$WAA[1, y, seas, , 1] <- in_dat$wt_ssb[y, 1:n_ages]
    ebswp_SPoRC_data$MatAA[1, y, seas, , 1] <- in_dat$p_mature[1:n_ages]
  }
}

# Ageing error
ebswp_SPoRC_data$AgeingError <- as.matrix(in_dat$age_err)

# Catch --------------------------------------------------------------------
# ObsCatch: region, year, seas, fleet
ebswp_SPoRC_data$ObsCatch <- array(0, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))
ebswp_SPoRC_data$ObsCatch[1, , 1, 1] <- in_dat$obs_catch

# UseCatch: region, year, seas, fleet
ebswp_SPoRC_data$UseCatch <- array(1, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))

# Fishery Indices ----------------------------------------------------------
# ObsFishIdx: region, year, seas, fleet
ebswp_SPoRC_data$ObsFishIdx <- array(NA, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))
ebswp_SPoRC_data$ObsFishIdx[1, 2:13, 1, 1] <- in_dat$obs_cpue

# ObsFishIdx_SE: region, year, seas, fleet
ebswp_SPoRC_data$ObsFishIdx_SE <- array(NA, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))
ebswp_SPoRC_data$ObsFishIdx_SE[1, 2:13, 1, 1] <- in_dat$obs_cpue_std / in_dat$obs_cpue

# UseFishIdx: region, year, seas, fleet
ebswp_SPoRC_data$UseFishIdx <- array(0, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))
ebswp_SPoRC_data$UseFishIdx[1, 2:13, 1, 1] <- 1

# Fishery Age Compositions -------------------------------------------------
# ObsFishAgeComps: region, year, seas, age, sex, fleet
ebswp_SPoRC_data$ObsFishAgeComps <- array(NA_real_, dim = c(n_regions, n_yrs, n_seas, n_ages, n_sexes, n_fish_fleets))

# UseFishAgeComps: region, year, seas, fleet
ebswp_SPoRC_data$UseFishAgeComps <- array(0, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))

# ISS_FishAgeComps: region, year, seas, sex, fleet
ebswp_SPoRC_data$ISS_FishAgeComps <- array(0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_fish_fleets))

idx <- which(yrs %in% in_dat$yrs_fsh_data)
for(i in 1:length(idx)) {
  ebswp_SPoRC_data$ObsFishAgeComps[1, idx[i], 1, , 1, 1] <- in_dat$oac_fsh_data[i, 1:n_ages]
  ebswp_SPoRC_data$UseFishAgeComps[1, idx[i], 1, 1] <- 1
  ebswp_SPoRC_data$ISS_FishAgeComps[1, idx[i], 1, 1, 1] <- in_dat$sam_fsh[i]
}

# Fishery Length Compositions (empty) --------------------------------------
# ObsFishLenComps: region, year, seas, len, sex, fleet
# UseFishLenComps: region, year, seas, fleet
# ISS_FishLenComps: region, year, seas, sex, fleet
# (not used for this model - leave unset or add empty arrays if needed)

# Survey Indices -----------------------------------------------------------
# ObsSrvIdx: region, year, seas, fleet
ebswp_SPoRC_data$ObsSrvIdx <- array(NA_real_, dim = c(n_regions, n_yrs, n_seas, n_srv_fleets))

# ObsSrvIdx_SE: region, year, seas, fleet
ebswp_SPoRC_data$ObsSrvIdx_SE <- array(NA_real_, dim = c(n_regions, n_yrs, n_seas, n_srv_fleets))

# UseSrvIdx: region, year, seas, fleet
ebswp_SPoRC_data$UseSrvIdx <- array(0, dim = c(n_regions, n_yrs, n_seas, n_srv_fleets))

# Survey Age Compositions --------------------------------------------------
# ObsSrvAgeComps: region, year, seas, age, sex, fleet
ebswp_SPoRC_data$ObsSrvAgeComps <- array(NA_real_, dim = c(n_regions, n_yrs, n_seas, n_ages, n_sexes, n_srv_fleets))

# UseSrvAgeComps: region, year, seas, fleet
ebswp_SPoRC_data$UseSrvAgeComps <- array(0, dim = c(n_regions, n_yrs, n_seas, n_srv_fleets))

# ISS_SrvAgeComps: region, year, seas, sex, fleet
ebswp_SPoRC_data$ISS_SrvAgeComps <- array(0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_srv_fleets))

# Bottom trawl survey (fleet 1)
idx <- which(yrs %in% in_dat$yrs_bts)
for(i in 1:length(idx)) {
  ebswp_SPoRC_data$ObsSrvIdx[1, idx[i], 1, 1] <- in_dat$ob_bts[i]
  ebswp_SPoRC_data$ObsSrvIdx_SE[1, idx[i], 1, 1] <- in_dat$ob_bts_std[i]
  ebswp_SPoRC_data$UseSrvIdx[1, idx[i], 1, 1] <- 1
  ebswp_SPoRC_data$ObsSrvAgeComps[1, idx[i], 1, , 1, 1] <- in_dat$oac_bts[i, 1:n_ages]
  ebswp_SPoRC_data$UseSrvAgeComps[1, idx[i], 1, 1] <- 1
  ebswp_SPoRC_data$ISS_SrvAgeComps[1, idx[i], 1, 1, 1] <- in_dat$sam_bts[i]
}

# ATS survey (fleet 2)
idx <- which(yrs %in% in_dat$yrs_ats)
for(i in 1:length(idx)) {
  ebswp_SPoRC_data$ObsSrvIdx[1, idx[i], 1, 2] <- in_dat$ob_ats[i]
  ebswp_SPoRC_data$ObsSrvIdx_SE[1, idx[i], 1, 2] <- in_dat$ob_ats_std[i]
  ebswp_SPoRC_data$UseSrvIdx[1, idx[i], 1, 2] <- 1
  ebswp_SPoRC_data$ObsSrvAgeComps[1, idx[i], 1, , 1, 2] <- in_dat$oac_ats[i, 1:n_ages]
  ebswp_SPoRC_data$UseSrvAgeComps[1, idx[i], 1, 2] <- 1
  ebswp_SPoRC_data$ISS_SrvAgeComps[1, idx[i], 1, 1, 2] <- in_dat$sam_ats[i]
}

# AVO survey (fleet 3) - index only, no age comps
idx <- which(yrs %in% in_dat$yrs_avo)
for(i in 1:length(idx)) {
  ebswp_SPoRC_data$ObsSrvIdx[1, idx[i], 1, 3] <- in_dat$ob_avo[i]
  ebswp_SPoRC_data$ObsSrvIdx_SE[1, idx[i], 1, 3] <- in_dat$ob_avo_std[i]
  ebswp_SPoRC_data$UseSrvIdx[1, idx[i], 1, 3] <- 1
}

# Convert SE to CV
ebswp_SPoRC_data$ObsSrvIdx_SE <- ebswp_SPoRC_data$ObsSrvIdx_SE / ebswp_SPoRC_data$ObsSrvIdx

# Save report values from actual assessment
ebswp_SPoRC_data$SSB <- rep$SSB
ebswp_SPoRC_data$R <- rep$R

# Write out data -----------------------------------------------------------
sgl_rg_ebswp_data <- ebswp_SPoRC_data
usethis::use_data(sgl_rg_ebswp_data, internal = FALSE, overwrite = TRUE)
