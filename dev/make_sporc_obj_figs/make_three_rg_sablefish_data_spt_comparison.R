library(here)
library(tidyverse)

# Load in data from previous spatial model and build this
spatial_data <- readRDS(here("dev", 'dev_data', "Spatial Sablefish Model", "data_3_area.RDS"))
ageing_dat <- dget(here("dev", 'dev_data','2023 Base (23.5)_final model', 'test.rdat')) # for getting ageing error

# load in tagging data
tag_rel <- readRDS(here("dev", 'dev_data', "Spatial Sablefish Model", "Tag_release_summarised_3_area.RDS"))
tag_rec <- readRDS(here("dev", 'dev_data', "Spatial Sablefish Model", "Tag_recovery_summarised_3_area.RDS"))

# create empty list to populate
three_rg_sable_data <- list()

# dimensions
years <- three_rg_sable_data$years <- spatial_data$years
ages <- three_rg_sable_data$ages <- spatial_data$ages
lens <- three_rg_sable_data$lens <- spatial_data$length_bins
n_yrs <- length(years)
n_pop <- three_rg_sable_data$n_pop <- 1
three_rg_sable_data$natal_region <- 1
n_ages <- length(ages)
n_lens <- length(three_rg_sable_data$lens)
n_regions <- three_rg_sable_data$n_regions <- 3
n_sexes <- three_rg_sable_data$n_sexes <- 2
n_fish_fleets <- three_rg_sable_data$n_fish_fleets <- 2
n_srv_fleets <- three_rg_sable_data$n_srv_fleets <- 2
n_seas <- three_rg_sable_data$n_seas <- 1
three_rg_sable_data$spawn_seas <- 1
three_rg_sable_data$seasdur <- 1


# biologicals
# Weight at age
three_rg_sable_data$WAA <- array(0, dim = c(n_pop, n_regions,n_yrs, n_seas, n_ages, n_sexes))

# Maturity at age
three_rg_sable_data$MatAA <- array(0, dim = c(n_pop, n_regions,n_yrs, n_seas,n_ages, n_sexes))

for(r in 1:n_regions) {
  # weight at age
  three_rg_sable_data$WAA[,r,,1,,1] <- t(spatial_data$female_mean_weight_by_age) # female weight at age
  three_rg_sable_data$WAA[,r,,1,,2] <- t(spatial_data$male_mean_weight_by_age) # male weight at age

  # maturity
  three_rg_sable_data$MatAA[,r,,1,,1] <- t(spatial_data$maturity) # female maturity at age
  three_rg_sable_data$MatAA[,r,,1,,2] <- t(spatial_data$maturity) # male maturity at age (not used)
}

# ageing error
three_rg_sable_data$AgeingError <- as.matrix(ageing_dat$age_error) # ageing error matrix

# Size Age transition matrix
three_rg_sable_data$SizeAgeTrans <- array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_lens, n_ages, n_sexes)) # size age transition matrix
for(p in 1:n_pop) {
  for(y in 1:n_yrs) {
    for(r in 1:n_regions) {
      three_rg_sable_data$SizeAgeTrans[p,r,y,1,,,1] <- t(spatial_data$female_age_length_transition[,,y])
      three_rg_sable_data$SizeAgeTrans[p,r,y,1,,,1] <- apply(three_rg_sable_data$SizeAgeTrans[p,r,y,1,,,1], 2, function(x) x / sum(x))
      three_rg_sable_data$SizeAgeTrans[p,r,y,1,,,2] <- t(spatial_data$male_age_length_transition[,,y]) / colSums(t(spatial_data$male_age_length_transition[,,y]))
      three_rg_sable_data$SizeAgeTrans[p,r,y,1,,,2] <- apply(three_rg_sable_data$SizeAgeTrans[p,r,y,1,,,2], 2, function(x) x / sum(x))
    } # end r loop
  } # end y loop
}  # end p loop

# Tagging data
tag_rel <- tag_rel %>% mutate(region_num = case_when(
  region_release == "BS_AI_WGOA" ~ 1,
  region_release == "CGOA" ~ 2,
  region_release == "EGOA" ~ 3
)) %>% filter(release_year %in% c(three_rg_sable_data$years))

# name regions correctly
tag_rec <- tag_rec %>% mutate(region_num = case_when(
  recovery_region == "BS_AI_WGOA" ~ 1,
  recovery_region == "CGOA" ~ 2,
  recovery_region == "EGOA" ~ 3
)) %>% filter(recovery_year %in% c(three_rg_sable_data$years)) %>%
  mutate(tag_lib = (recovery_year - release_year) + 1)

# loop through to get tag release indicator
rel_cohorts_df <- tag_rel %>% distinct(release_event_id, release_year, region_release, region_num)
tag_release_ind <- data.frame()
for(i in 1:nrow(rel_cohorts_df)) {
  tmp <- data.frame(regions = rel_cohorts_df[i,4], tag_yrs = rel_cohorts_df[i,2]-1959, tag_seas = 1)
  tag_release_ind <- rbind(tag_release_ind, tmp)
} # end i

three_rg_sable_data$conv_tag_release_indicator <- as.matrix(tag_release_ind) # input releases of cohorts in
three_rg_sable_data$n_conv_tag_cohorts <- nrow(three_rg_sable_data$conv_tag_release_indicator) # number of tag cohorts
three_rg_sable_data$conv_tag_max_tag_liberty <- 15 # maximum liberty to track cohorts

# Get Tagged Fish
three_rg_sable_data$conv_tagged_fish <- array(0, dim = c(three_rg_sable_data$n_conv_tag_cohorts, n_pop, n_ages, n_sexes)) # tagged fish
for(i in 1:nrow(three_rg_sable_data$conv_tag_release_indicator )) {
  tmp_rel_f <- tag_rel %>% filter(release_event_id == i, sex == "F") # filter to a given cohort females
  tmp_rel_m <- tag_rel %>% filter(release_event_id == i, sex == "M") # filter to a given cohort males
  three_rg_sable_data$conv_tagged_fish[i,1,,1] <- tmp_rel_f$Nage_at_release # females
  three_rg_sable_data$conv_tagged_fish[i,1,,2] <- tmp_rel_m$Nage_at_release # males
} # end i for tag cohorts


# Get Observed Recaptures
three_rg_sable_data$obs_conv_tag_fish_recap <- array(0, dim = c(three_rg_sable_data$conv_tag_max_tag_liberty,
                                                              n_seas,
                                                              three_rg_sable_data$n_conv_tag_cohorts,
                                                              n_pop,
                                                              n_regions,
                                                              n_ages,
                                                              n_sexes,
                                                              n_fish_fleets))

# pre filter data frame before looping (only recoveries from fixed gear included)
tag_rec_f <- tag_rec %>% filter(sex == "F")
tag_rec_m <- tag_rec %>% filter(sex == "M")

for(ry in 1:three_rg_sable_data$conv_tag_max_tag_liberty) {
  for(i in 1:three_rg_sable_data$n_conv_tag_cohorts) {
    for(r in 1:n_regions) {
      tmp_rec_f <- tag_rec_f %>% filter(release_event_id == i, region_num == r, tag_lib == ry)
      tmp_rec_m <- tag_rec_m %>% filter(release_event_id == i, region_num == r, tag_lib == ry)
      if(nrow(tmp_rec_f) > 0) three_rg_sable_data$obs_conv_tag_fish_recap[ry,1,i,1,r,,1,1] <- tmp_rec_f$Nage_at_recovery
      if(nrow(tmp_rec_m) > 0) three_rg_sable_data$obs_conv_tag_fish_recap[ry,1,i,1,r,,2,1] <- tmp_rec_m$Nage_at_recovery
    }
  }
}


# Catch data
three_rg_sable_data$ObsCatch <- array(0, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))
three_rg_sable_data$ObsCatch[,,1,1] <- spatial_data$fixed_fishery_catch # fixed gear catch
three_rg_sable_data$ObsCatch[,,1,2] <- spatial_data$trwl_fishery_catch # trawl gear catch
three_rg_sable_data$ObsCatch[,1:3,1,2] <- 0 # input 0 for first three years
three_rg_sable_data$UseCatch <- array(1, c(n_regions, n_yrs, n_seas, n_fish_fleets)) # fit catch data everywhere
three_rg_sable_data$UseCatch[,1:3,1,2] <- 0 # don't use catch for the first three years

# Fishery Indices (not fit to)
three_rg_sable_data$ObsFishIdx <- array(NA, c(n_regions, n_yrs, n_seas, n_fish_fleets))
three_rg_sable_data$ObsFishIdx_SE <- array(NA, c(n_regions, n_yrs, n_seas, n_fish_fleets))
three_rg_sable_data$UseFishIdx <- array(0, c(n_regions, n_yrs, n_seas, n_fish_fleets))

# Fishery Age Compositions (Joint by sex, split by region)
three_rg_sable_data$ObsFishAgeComps <- array(0, dim = c(n_regions, n_yrs, n_seas,
                                                      n_ages, n_sexes, n_fish_fleets))
# loop through to populate
for(y in 1:n_yrs) {
  for(r in 1:n_regions) {
    if(spatial_data$fixed_catchatage_indicator[r,y] == 1) {
      three_rg_sable_data$ObsFishAgeComps[r,y,1,,1,1] <- spatial_data$obs_fixed_catchatage[-(1:n_ages),r,y] / sum(spatial_data$obs_fixed_catchatage[,r,y]) # females
      three_rg_sable_data$ObsFishAgeComps[r,y,1,,2,1] <- spatial_data$obs_fixed_catchatage[1:n_ages,r,y] / sum(spatial_data$obs_fixed_catchatage[,r,y]) # males
    } # end if
  } # end r loop
} # end y loop

# use indicators
three_rg_sable_data$UseFishAgeComps <- array(0, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))
three_rg_sable_data$UseFishAgeComps[,,1,1] <- spatial_data$fixed_catchatage_indicator # fixed gear use indicator

# input sample size
three_rg_sable_data$ISS_FishAgeComps <- array(0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_fish_fleets))
three_rg_sable_data$ISS_FishAgeComps[] <- 50 # ISS for fixed gear fishery

# Fishery Length Compositions (Joint by sex, split by region)
three_rg_sable_data$ObsFishLenComps <- array(0, dim = c(n_regions, n_yrs, n_seas, n_lens, n_sexes, n_fish_fleets))
# loop through to populate
for(y in 1:n_yrs) {
  for(r in 1:n_regions) {

    # fixed gear
    if(spatial_data$fixed_catchatlgth_indicator[r,y] == 1) {
      three_rg_sable_data$ObsFishLenComps[r,y,1,,1,1] <- spatial_data$obs_fixed_catchatlgth[-(1:n_ages),r,y] / sum(spatial_data$obs_fixed_catchatlgth[,r,y]) # females
      three_rg_sable_data$ObsFishLenComps[r,y,1,,2,1] <- spatial_data$obs_fixed_catchatlgth[1:n_ages,r,y] / sum(spatial_data$obs_fixed_catchatlgth[,r,y]) # males
    } # end ifM

    # trawl gear
    if(spatial_data$trwl_catchatlgth_indicator[r,y] == 1) {
      three_rg_sable_data$ObsFishLenComps[r,y,1,,1,2] <- spatial_data$obs_trwl_catchatlgth[-(1:n_ages),r,y] / sum(spatial_data$obs_trwl_catchatlgth[,r,y]) # females
      three_rg_sable_data$ObsFishLenComps[r,y,1,,2,2] <- spatial_data$obs_trwl_catchatlgth[1:n_ages,r,y] / sum(spatial_data$obs_trwl_catchatlgth[,r,y]) # males
    } # end if

  } # end r loop
} # end y loop

# use indicators
three_rg_sable_data$UseFishLenComps <- array(0, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))
three_rg_sable_data$UseFishLenComps[,,1,1] <- spatial_data$fixed_catchatlgth_indicator # fixed gear use indicator
three_rg_sable_data$UseFishLenComps[,,1,2] <- spatial_data$trwl_catchatlgth_indicator # trawl gear use indicator

# input sample size
three_rg_sable_data$ISS_FishLenComps <- array(0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_fish_fleets))
three_rg_sable_data$ISS_FishLenComps[] <- 40 # fixed gear lens ISS

# Survey Indices (fleet 1 = japanese, fleet 2 = domestic)
three_rg_sable_data$ObsSrvIdx_SE <- three_rg_sable_data$ObsSrvIdx <- array(0, dim = c(n_regions, n_yrs, n_seas, n_srv_fleets))
three_rg_sable_data$ObsSrvIdx[,,1,] <- spatial_data$obs_srv_bio
three_rg_sable_data$ObsSrvIdx_SE[,,1,] <- spatial_data$obs_srv_se

three_rg_sable_data$UseSrvIdx <- array(0, dim = c(n_regions, n_yrs, n_seas, n_srv_fleets))
three_rg_sable_data$UseSrvIdx[,,1,] <- spatial_data$srv_bio_indicator # survey indicator

# Survey Age Comps
three_rg_sable_data$ObsSrvAgeComps <- array(0, dim = c(n_regions, n_yrs, n_seas,
                                                     n_ages, n_sexes, n_srv_fleets))

# loop through to populate
for(y in 1:n_yrs) {
  for(r in 1:n_regions) {

    if(spatial_data$srv_catchatage_indicator[r,y,1] == 1) { # coop survey
      three_rg_sable_data$ObsSrvAgeComps[r,y,1,,1,1] <- spatial_data$obs_srv_catchatage[-(1:n_ages),r,y,1] / sum(spatial_data$obs_srv_catchatage[,r,y,1]) # females
      three_rg_sable_data$ObsSrvAgeComps[r,y,1,,2,1] <- spatial_data$obs_srv_catchatage[1:n_ages,r,y,1] / sum(spatial_data$obs_srv_catchatage[,r,y,1]) # males
    } # end if

    if(spatial_data$srv_catchatage_indicator[r,y,2] == 1) { # domestic survey
      three_rg_sable_data$ObsSrvAgeComps[r,y,1,,1,2] <- spatial_data$obs_srv_catchatage[-(1:n_ages),r,y,2] / sum(spatial_data$obs_srv_catchatage[,r,y,2]) # females
      three_rg_sable_data$ObsSrvAgeComps[r,y,1,,2,2] <- spatial_data$obs_srv_catchatage[1:n_ages,r,y,2] / sum(spatial_data$obs_srv_catchatage[,r,y,2]) # males
    } # end if

  } # end r loop
} # end y loop

# use indicators
three_rg_sable_data$UseSrvAgeComps <- array(0, dim = c(n_regions, n_yrs, n_seas, n_srv_fleets))
three_rg_sable_data$UseSrvAgeComps[,,1,] <- spatial_data$srv_catchatage_indicator

# input sample sizes
three_rg_sable_data$ISS_SrvAgeComps <- array(0, dim = c(n_regions, n_yrs, n_seas,
                                                      n_sexes, n_srv_fleets))
three_rg_sable_data$ISS_SrvAgeComps[] <- 60 # input sample size of 60 for all survey age comps

# Survey length compositions (not used)
three_rg_sable_data$ObsSrvLenComps <- array(NA, dim = c(n_regions,n_yrs, n_seas, n_lens, n_sexes, n_srv_fleets))
three_rg_sable_data$UseSrvLenComps <- array(0, dim = c(n_regions, n_yrs, n_seas, n_srv_fleets))
three_rg_sable_data$ISS_SrvLenComps <- array(0, dim = c(n_regions,n_yrs, n_seas, n_sexes, n_srv_fleets))

three_rg_sable_data$do_recruits_move <- 0 # recruit's don't move

# Write data
usethis::use_data(three_rg_sable_data, internal = FALSE, overwrite = TRUE)

