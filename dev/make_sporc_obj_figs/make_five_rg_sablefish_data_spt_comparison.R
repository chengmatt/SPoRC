library(here)

# Load in data from previous spatial model and build this
spatial_data <- readRDS(here("dev", 'dev_data', "Spatial Sablefish Model", "data.RDS"))
ageing_dat <- dget(here("dev", 'dev_data','2023 Base (23.5)_final model', 'test.rdat')) # for getting ageing error

# load in tagging data
tag_rel <- readRDS(here("dev", 'dev_data', "Spatial Sablefish Model", "Tag_release_summarised.RDS"))
tag_rec <- readRDS(here("dev", 'dev_data', "Spatial Sablefish Model", "Tag_recovery_summarised.RDS"))

# create empty list to populate
mlt_rg_sable_data <- list()

# dimensions
years <- mlt_rg_sable_data$years <- spatial_data$years
ages <- mlt_rg_sable_data$ages <- spatial_data$ages
lens <- mlt_rg_sable_data$lens <- spatial_data$length_bins
n_yrs <- length(years)
n_ages <- length(ages)
n_lens <- length(mlt_rg_sable_data$lens)
n_regions <- mlt_rg_sable_data$n_regions <- 5
n_sexes <- mlt_rg_sable_data$n_sexes <- 2
n_fish_fleets <- mlt_rg_sable_data$n_fish_fleets <- 2
n_srv_fleets <- mlt_rg_sable_data$n_srv_fleets <- 2
n_seas <- mlt_rg_sable_data$n_seas <- 1
mlt_rg_sable_data$spawn_seas <- 1
mlt_rg_sable_data$seasdur <- 1

# biologicals
# Weight at age
mlt_rg_sable_data$WAA <- array(0, dim = c(n_regions,n_yrs, n_seas, n_ages, n_sexes))

# Maturity at age
mlt_rg_sable_data$MatAA <- array(0, dim = c(n_regions,n_yrs, n_seas,n_ages, n_sexes))

for(r in 1:n_regions) {
  # weight at age
  mlt_rg_sable_data$WAA[r,,1,,1] <- t(spatial_data$female_mean_weight_by_age) # female weight at age
  mlt_rg_sable_data$WAA[r,,1,,2] <- t(spatial_data$male_mean_weight_by_age) # male weight at age

  # maturity
  mlt_rg_sable_data$MatAA[r,,1,,1] <- t(spatial_data$maturity) # female maturity at age
  mlt_rg_sable_data$MatAA[r,,1,,2] <- t(spatial_data$maturity) # male maturity at age (not used)
}

# ageing error
mlt_rg_sable_data$AgeingError <- as.matrix(ageing_dat$age_error) # ageing error matrix

# Size Age transition matrix
mlt_rg_sable_data$SizeAgeTrans <- array(0, dim = c(n_regions, n_yrs, n_seas, n_lens, n_ages, n_sexes)) # size age transition matrix
for(y in 1:n_yrs) {
  for(r in 1:n_regions) {
    mlt_rg_sable_data$SizeAgeTrans[r,y,1,,,1] <- t(spatial_data$female_age_length_transition[,,y])
    mlt_rg_sable_data$SizeAgeTrans[r,y,1,,,1] <- apply(mlt_rg_sable_data$SizeAgeTrans[r,y,1,,,1], 2, function(x) x / sum(x))
    mlt_rg_sable_data$SizeAgeTrans[r,y,1,,,2] <- t(spatial_data$male_age_length_transition[,,y]) / colSums(t(spatial_data$male_age_length_transition[,,y]))
    mlt_rg_sable_data$SizeAgeTrans[r,y,1,,,2] <- apply(mlt_rg_sable_data$SizeAgeTrans[r,y,1,,,2], 2, function(x) x / sum(x))
  } # end r loop
} # end y loop

# Tagging data
tag_rel <- tag_rel %>% mutate(region_num = case_when(
  region_release == "BS" ~ 1,
  region_release == "AI" ~ 2,
  region_release == "WGOA" ~ 3,
  region_release == "CGOA" ~ 4,
  region_release == "EGOA" ~ 5
)) %>% filter(release_year %in% c(mlt_rg_sable_data$years))

# name regions correctly
tag_rec <- tag_rec %>% mutate(region_num = case_when(
  recovery_region == "BS" ~ 1,
  recovery_region == "AI" ~ 2,
  recovery_region == "WGOA" ~ 3,
  recovery_region == "CGOA" ~ 4,
  recovery_region == "EGOA" ~ 5
)) %>% filter(recovery_year %in% c(mlt_rg_sable_data$years)) %>%
  mutate(tag_lib = (recovery_year - release_year) + 1)

# loop through to get tag release indicator
rel_cohorts_df <- tag_rel %>% distinct(release_event_id, release_year, region_release, region_num)
tag_release_ind <- data.frame()
for(i in 1:nrow(rel_cohorts_df)) {
  tmp <- data.frame(regions = rel_cohorts_df[i,4], tag_yrs = rel_cohorts_df[i,2]-1959, tag_seas = 1)
  tag_release_ind <- rbind(tag_release_ind, tmp)
} # end i

mlt_rg_sable_data$tag_release_indicator <- as.matrix(tag_release_ind) # input releases of cohorts in
mlt_rg_sable_data$n_tag_cohorts <- nrow(mlt_rg_sable_data$tag_release_indicator) # number of tag cohorts
mlt_rg_sable_data$max_tag_liberty <- 15 # maximum liberty to track cohorts

# Get Tagged Fish
mlt_rg_sable_data$Tagged_Fish <- array(0, dim = c(mlt_rg_sable_data$n_tag_cohorts, n_ages, n_sexes)) # tagged fish
for(i in 1:nrow(mlt_rg_sable_data$tag_release_indicator )) {
  tmp_rel_f <- tag_rel %>% filter(release_event_id == i, sex == "F") # filter to a given cohort females
  tmp_rel_m <- tag_rel %>% filter(release_event_id == i, sex == "M") # filter to a given cohort males
  mlt_rg_sable_data$Tagged_Fish[i,,1] <- tmp_rel_f$Nage_at_release # females
  mlt_rg_sable_data$Tagged_Fish[i,,2] <- tmp_rel_m$Nage_at_release # males
} # end i for tag cohorts


# Get Observed Recaptures
mlt_rg_sable_data$Obs_Tag_Recap <- array(0, dim = c(mlt_rg_sable_data$max_tag_liberty,
                                                      n_seas,
                                                      mlt_rg_sable_data$n_tag_cohorts,
                                                      n_regions,
                                                      n_ages,
                                                      n_sexes))

# pre filter data frame before looping
tag_rec_f <- tag_rec %>% filter(sex == "F")
tag_rec_m <- tag_rec %>% filter(sex == "M")

for(ry in 1:mlt_rg_sable_data$max_tag_liberty) {
  for(i in 1:mlt_rg_sable_data$n_tag_cohorts) {
    for(r in 1:n_regions) {
      tmp_rec_f <- tag_rec_f %>% filter(release_event_id == i, region_num == r, tag_lib == ry)
      tmp_rec_m <- tag_rec_m %>% filter(release_event_id == i, region_num == r, tag_lib == ry)
      if(nrow(tmp_rec_f) > 0) mlt_rg_sable_data$Obs_Tag_Recap[ry,1,i,r,,1] <- tmp_rec_f$Nage_at_recovery
      if(nrow(tmp_rec_m) > 0) mlt_rg_sable_data$Obs_Tag_Recap[ry,1,i,r,,2] <- tmp_rec_m$Nage_at_recovery
    }
  }
}


# Catch data
mlt_rg_sable_data$ObsCatch <- array(0, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))
mlt_rg_sable_data$ObsCatch[,,1,1] <- spatial_data$fixed_fishery_catch # fixed gear catch
mlt_rg_sable_data$ObsCatch[,,1,2] <- spatial_data$trwl_fishery_catch # trawl gear catch
mlt_rg_sable_data$ObsCatch[,1:3,1,2] <- 0 # input 0 for first three years
mlt_rg_sable_data$UseCatch <- array(1, c(n_regions, n_yrs, n_seas, n_fish_fleets)) # fit catch data everywhere
mlt_rg_sable_data$UseCatch[,1:3,1,2] <- 0 # don't use catch for the first three years

# Fishery Indices (not fit to)
mlt_rg_sable_data$ObsFishIdx <- array(NA, c(n_regions, n_yrs, n_seas, n_fish_fleets))
mlt_rg_sable_data$ObsFishIdx_SE <- array(NA, c(n_regions, n_yrs, n_seas, n_fish_fleets))
mlt_rg_sable_data$UseFishIdx <- array(0, c(n_regions, n_yrs, n_seas, n_fish_fleets))

# Fishery Age Compositions (Joint by sex, split by region)
mlt_rg_sable_data$ObsFishAgeComps <- array(0, dim = c(n_regions, n_yrs, n_seas,
                                                        n_ages, n_sexes, n_fish_fleets))
# loop through to populate
for(y in 1:n_yrs) {
  for(r in 1:n_regions) {
    if(spatial_data$fixed_catchatage_indicator[r,y] == 1) {
      mlt_rg_sable_data$ObsFishAgeComps[r,y,1,,1,1] <- spatial_data$obs_fixed_catchatage[-(1:n_ages),r,y] / sum(spatial_data$obs_fixed_catchatage[,r,y]) # females
      mlt_rg_sable_data$ObsFishAgeComps[r,y,1,,2,1] <- spatial_data$obs_fixed_catchatage[1:n_ages,r,y] / sum(spatial_data$obs_fixed_catchatage[,r,y]) # males
    } # end if
  } # end r loop
} # end y loop

# use indicators
mlt_rg_sable_data$UseFishAgeComps <- array(0, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))
mlt_rg_sable_data$UseFishAgeComps[,,1,1] <- spatial_data$fixed_catchatage_indicator # fixed gear use indicator

# input sample size
mlt_rg_sable_data$ISS_FishAgeComps <- array(0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_fish_fleets))
mlt_rg_sable_data$ISS_FishAgeComps[] <- 50 # ISS for fixed gear fishery

# Fishery Length Compositions (Joint by sex, split by region)
mlt_rg_sable_data$ObsFishLenComps <- array(0, dim = c(n_regions, n_yrs, n_seas, n_lens, n_sexes, n_fish_fleets))
# loop through to populate
for(y in 1:n_yrs) {
  for(r in 1:n_regions) {

    # fixed gear
    if(spatial_data$fixed_catchatlgth_indicator[r,y] == 1) {
      mlt_rg_sable_data$ObsFishLenComps[r,y,1,,1,1] <- spatial_data$obs_fixed_catchatlgth[-(1:n_ages),r,y] / sum(spatial_data$obs_fixed_catchatlgth[,r,y]) # females
      mlt_rg_sable_data$ObsFishLenComps[r,y,1,,2,1] <- spatial_data$obs_fixed_catchatlgth[1:n_ages,r,y] / sum(spatial_data$obs_fixed_catchatlgth[,r,y]) # males
    } # end ifM

    # trawl gear
    if(spatial_data$trwl_catchatlgth_indicator[r,y] == 1) {
      mlt_rg_sable_data$ObsFishLenComps[r,y,1,,1,2] <- spatial_data$obs_trwl_catchatlgth[-(1:n_ages),r,y] / sum(spatial_data$obs_trwl_catchatlgth[,r,y]) # females
      mlt_rg_sable_data$ObsFishLenComps[r,y,1,,2,2] <- spatial_data$obs_trwl_catchatlgth[1:n_ages,r,y] / sum(spatial_data$obs_trwl_catchatlgth[,r,y]) # males
    } # end if

  } # end r loop
} # end y loop

# use indicators
mlt_rg_sable_data$UseFishLenComps <- array(0, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))
mlt_rg_sable_data$UseFishLenComps[,,1,1] <- spatial_data$fixed_catchatlgth_indicator # fixed gear use indicator
mlt_rg_sable_data$UseFishLenComps[,,1,2] <- spatial_data$trwl_catchatlgth_indicator # trawl gear use indicator

# input sample size
mlt_rg_sable_data$ISS_FishLenComps <- array(0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_fish_fleets))
mlt_rg_sable_data$ISS_FishLenComps[] <- 40 # fixed gear lens ISS

# Survey Indices (fleet 1 = japanese, fleet 2 = domestic)
mlt_rg_sable_data$ObsSrvIdx_SE <- mlt_rg_sable_data$ObsSrvIdx <- array(0, dim = c(n_regions, n_yrs, n_seas, n_srv_fleets))
mlt_rg_sable_data$ObsSrvIdx[,,1,] <- spatial_data$obs_srv_bio
mlt_rg_sable_data$ObsSrvIdx_SE[,,1,] <- spatial_data$obs_srv_se

mlt_rg_sable_data$UseSrvIdx <- array(0, dim = c(n_regions, n_yrs, n_seas, n_srv_fleets))
mlt_rg_sable_data$UseSrvIdx[,,1,] <- spatial_data$srv_bio_indicator # survey indicator

# Survey Age Comps
mlt_rg_sable_data$ObsSrvAgeComps <- array(0, dim = c(n_regions, n_yrs, n_seas,
                                                       n_ages, n_sexes, n_srv_fleets))

# loop through to populate
for(y in 1:n_yrs) {
  for(r in 1:n_regions) {

    if(spatial_data$srv_catchatage_indicator[r,y,1] == 1) { # coop survey
      mlt_rg_sable_data$ObsSrvAgeComps[r,y,1,,1,1] <- spatial_data$obs_srv_catchatage[-(1:n_ages),r,y,1] / sum(spatial_data$obs_srv_catchatage[,r,y,1]) # females
      mlt_rg_sable_data$ObsSrvAgeComps[r,y,1,,2,1] <- spatial_data$obs_srv_catchatage[1:n_ages,r,y,1] / sum(spatial_data$obs_srv_catchatage[,r,y,1]) # males
    } # end if

    if(spatial_data$srv_catchatage_indicator[r,y,2] == 1) { # domestic survey
      mlt_rg_sable_data$ObsSrvAgeComps[r,y,1,,1,2] <- spatial_data$obs_srv_catchatage[-(1:n_ages),r,y,2] / sum(spatial_data$obs_srv_catchatage[,r,y,2]) # females
      mlt_rg_sable_data$ObsSrvAgeComps[r,y,1,,2,2] <- spatial_data$obs_srv_catchatage[1:n_ages,r,y,2] / sum(spatial_data$obs_srv_catchatage[,r,y,2]) # males
    } # end if

  } # end r loop
} # end y loop

# use indicators
mlt_rg_sable_data$UseSrvAgeComps <- array(0, dim = c(n_regions, n_yrs, n_seas, n_srv_fleets))
mlt_rg_sable_data$UseSrvAgeComps[,,1,] <- spatial_data$srv_catchatage_indicator

# input sample sizes
mlt_rg_sable_data$ISS_SrvAgeComps <- array(0, dim = c(n_regions, n_yrs, n_seas,
                                                        n_sexes, n_srv_fleets))
mlt_rg_sable_data$ISS_SrvAgeComps[] <- 60 # input sample size of 60 for all survey age comps

# Survey length compositions (not used)
mlt_rg_sable_data$ObsSrvLenComps <- array(NA, dim = c(n_regions,n_yrs, n_seas, n_lens, n_sexes, n_srv_fleets))
mlt_rg_sable_data$UseSrvLenComps <- array(0, dim = c(n_regions, n_yrs, n_seas, n_srv_fleets))
mlt_rg_sable_data$ISS_SrvLenComps <- array(0, dim = c(n_regions,n_yrs, n_seas, n_sexes, n_srv_fleets))

mlt_rg_sable_data$do_recruits_move <- 0 # recruit's don't move

# Write data
usethis::use_data(mlt_rg_sable_data, internal = FALSE, overwrite = TRUE)
