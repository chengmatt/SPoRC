library(devtools)
library(pkgdown)
library(roxygen2)
library(attachment)
library(usethis)
library(here)

# Build Package -----------------------------------------------------------
document() # document functions
roxygenise() # make sure functions have roxygen documentation
att_amend_desc(update = TRUE)
desc::desc_del_dep("compResidual", "Imports") # move to suggests
desc::desc_set_dep("compResidual", "Suggests") # move to suggests
desc::desc_del_dep("remotes", "Imports") # remove remotes
# usethis::use_news_md() # add news md
Sys.unsetenv("GITHUB_PAT") # may need to unset to build vignettes
# build_vignettes() # build vignettes
build_site(examples = FALSE)
use_build_ignore("dev") # ignore dev folder
usethis::use_build_ignore("_pkgdown.yml") # ignore pkgdown.yml
check() # check package stuff
# devtools::test()
# usethis::edit_git_ignore()
# unignore dev folder
# rbuildignore <- readLines(".Rbuildignore")
# rbuildignore <- rbuildignore[!grepl("^\\^dev\\$", rbuildignore)]
# writeLines(rbuildignore, ".Rbuildignore")
Sys.unsetenv("GITHUB_PAT") # may need to unset to build vignettes
build() # build package
install() # install locally
unloadNamespace('SPoRC')
# devtools::load_all(here("R"))


# Unit Tests --------------------------------------------------------------
# usethis::use_testthat()
# usethis::use_test("dusky_rtmb")
# usethis::use_test("sabie_sgl_rtmb")
# usethis::use_test("ebs_pol_sgl_rtmb")
# usethis::use_test("sabie_three_rg_rtmb")
# usethis::use_test("sgl_rg_simple_sim_test")

# Build Vignettes ---------------------------------------------------------

# usethis::use_vignette("a_model_dimensions")
# usethis::use_vignette("b_model_parameters")
# usethis::use_vignette("c_model_equations")
# usethis::use_vignette("d_model_report")
# usethis::use_vignette("e_single_region_sablefish_case_study")
# usethis::use_vignette("f_single_region_ebs_pollock_case_study")
# usethis::use_vignette("g_spatial_sablefish_case_study")
# usethis::use_vignette("h_closed_loop_simulations")
# usethis::use_vignette("i_reference_points")
# usethis::use_vignette("j_starting_mapping")
# usethis::use_vignette("k_defining_priors")
# usethis::use_vignette("l_simulation_testing")
# usethis::use_vignette("m_simulation_dimensions")
# usethis::use_vignette("n_single_region_ebs_pollock_randomeff_case_study")
# usethis::use_vignette("o_get_started")
# usethis::use_vignette("p_single_region_dusky_alt_mp_testing")
# usethis::use_vignette("q_movement_param")

# # Build Data Objects ------------------------------------------------------

### Single Region Sablefish -------------------------------------------------
# Load in data and build
sgl_rg_sable_data <- readRDS(here("dev", 'dev_data', "sgl_rg_sable_data.RDS")) # read in RTMB data
ageing_dat <- dget(here("dev", 'dev_data','2023 Base (23.5)_final model', 'test.rdat')) # for getting ageing error
sgl_rg_sable_data$age_error <- ageing_dat$age_error # input ageing error
sgl_rg_sable_data$years <- 1:65
sgl_rg_sable_data$ages <- 1:30
sgl_rg_sable_data$n_fish_fleets <- 2
sgl_rg_sable_data$n_srv_fleets <- 3
sgl_rg_sable_data$n_sexes <- 2

# get admb report and plug into rtmb data
tem_dat <- dget(here("dev", 'dev_data', '2024 Base (23.5)_final model_v3', 'tem.rdat'))
sgl_rg_sable_data$admb_recr <- tem_dat$t.series$Recr # recruitment
sgl_rg_sable_data$admb_spbiom <- tem_dat$t.series$spbiom # ssb
usethis::use_data(sgl_rg_sable_data, internal = FALSE, overwrite = TRUE)

# Single Region Sablefish Report
sgl_rg_sable_rep <- readRDS(here("dev", "dev_output", "1_Region_Model_Sablefish", "rep.RDS"))
usethis::use_data(sgl_rg_sable_rep, internal = FALSE, overwrite = TRUE)

### Five Region Sablefish --------------------------------------------------
# Load in data from previous spatial model and build this
spatial_data <- readRDS(here("dev", 'dev_data', "Spatial Sablefish Model", "data.RDS"))
ageing_dat <- dget(here("dev", 'dev_data','2023 Base (23.5)_final model', 'test.rdat')) # for getting ageing error

# load in tagging data
tag_rel <- readRDS(here("dev", 'dev_data', "Spatial Sablefish Model", "Tag_release_summarised.RDS"))
tag_rec <- readRDS(here("dev", 'dev_data', "Spatial Sablefish Model", "Tag_recovery_summarised.RDS"))

# create empty list to populate
mlt_rg_sable_data <- list()

# dimensions
mlt_rg_sable_data$years <- spatial_data$years
mlt_rg_sable_data$ages <- spatial_data$ages
mlt_rg_sable_data$lens <- spatial_data$length_bins
mlt_rg_sable_data$n_regions <- 5
mlt_rg_sable_data$n_sexes <- 2
mlt_rg_sable_data$n_fish_fleets <- 2
mlt_rg_sable_data$n_srv_fleets <- 2

# biologicals
# Weight at age
mlt_rg_sable_data$WAA <- array(0, dim = c(mlt_rg_sable_data$n_regions,length(mlt_rg_sable_data$years),length(mlt_rg_sable_data$ages), mlt_rg_sable_data$n_sexes))

# Maturity at age
mlt_rg_sable_data$MatAA <- array(0, dim = c(mlt_rg_sable_data$n_regions,length(mlt_rg_sable_data$years),length(mlt_rg_sable_data$ages), mlt_rg_sable_data$n_sexes))

for(r in 1:mlt_rg_sable_data$n_regions) {
  # weight at age
  mlt_rg_sable_data$WAA[r,,,1] <- t(spatial_data$female_mean_weight_by_age) # female weight at age
  mlt_rg_sable_data$WAA[r,,,2] <- t(spatial_data$male_mean_weight_by_age) # male weight at age

  # maturity
  mlt_rg_sable_data$MatAA[r,,,1] <- t(spatial_data$maturity) # female maturity at age
  mlt_rg_sable_data$MatAA[r,,,2] <- t(spatial_data$maturity) # male maturity at age (not used)
}

# ageing error
mlt_rg_sable_data$AgeingError <- as.matrix(ageing_dat$age_error) # ageing error matrix

# Size Age transition matrix
mlt_rg_sable_data$SizeAgeTrans <- array(0, dim = c(mlt_rg_sable_data$n_regions, length(mlt_rg_sable_data$years), length(mlt_rg_sable_data$lens),
                                                   length(mlt_rg_sable_data$ages), mlt_rg_sable_data$n_sexes)) # size age transition matrix
for(y in 1:length(mlt_rg_sable_data$years)) {
  for(r in 1:mlt_rg_sable_data$n_regions) {
    mlt_rg_sable_data$SizeAgeTrans[r,y,,,1] <- t(spatial_data$female_age_length_transition[,,y])
    mlt_rg_sable_data$SizeAgeTrans[r,y,,,1] <- apply(mlt_rg_sable_data$SizeAgeTrans[r,y,,,1], 2, function(x) x / sum(x))
    mlt_rg_sable_data$SizeAgeTrans[r,y,,,2] <- t(spatial_data$male_age_length_transition[,,y]) / colSums(t(spatial_data$male_age_length_transition[,,y]))
    mlt_rg_sable_data$SizeAgeTrans[r,y,,,2] <- apply(mlt_rg_sable_data$SizeAgeTrans[r,y,,,2], 2, function(x) x / sum(x))
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
  tmp <- data.frame(regions = rel_cohorts_df[i,4], tag_yrs = rel_cohorts_df[i,2]-1959)
  tag_release_ind <- rbind(tag_release_ind, tmp)
} # end i

mlt_rg_sable_data$tag_release_indicator <- as.matrix(tag_release_ind) # input releases of cohorts in
mlt_rg_sable_data$n_tag_cohorts <- nrow(mlt_rg_sable_data$tag_release_indicator) # number of tag cohorts
mlt_rg_sable_data$max_tag_liberty <- 15 # maximum liberty to track cohorts

# Get Tagged Fish
mlt_rg_sable_data$Tagged_Fish <- array(0, dim = c(mlt_rg_sable_data$n_tag_cohorts, length(mlt_rg_sable_data$ages), mlt_rg_sable_data$n_sexes)) # tagged fish
for(i in 1:nrow(mlt_rg_sable_data$tag_release_indicator )) {
  tmp_rel_f <- tag_rel %>% filter(release_event_id == i, sex == "F") # filter to a given cohort females
  tmp_rel_m <- tag_rel %>% filter(release_event_id == i, sex == "M") # filter to a given cohort males
  mlt_rg_sable_data$Tagged_Fish[i,,1] <- tmp_rel_f$Nage_at_release # females
  mlt_rg_sable_data$Tagged_Fish[i,,2] <- tmp_rel_m$Nage_at_release # males
} # end i for tag cohorts


# Get Observed Recaptures
mlt_rg_sable_data$Obs_Tag_Recap <- array(0, dim = c(mlt_rg_sable_data$max_tag_liberty, mlt_rg_sable_data$n_tag_cohorts, mlt_rg_sable_data$n_regions,
                                                    length(mlt_rg_sable_data$ages), mlt_rg_sable_data$n_sexes))

# pre filter data frame before looping
tag_rec_f <- tag_rec %>% filter(sex == "F")
tag_rec_m <- tag_rec %>% filter(sex == "M")

for(ry in 1:mlt_rg_sable_data$max_tag_liberty) {
  for(i in 1:mlt_rg_sable_data$n_tag_cohorts) {
    for(r in 1:mlt_rg_sable_data$n_regions) {
      tmp_rec_f <- tag_rec_f %>% filter(release_event_id == i, region_num == r, tag_lib == ry)
      tmp_rec_m <- tag_rec_m %>% filter(release_event_id == i, region_num == r, tag_lib == ry)
      if(nrow(tmp_rec_f) > 0) mlt_rg_sable_data$Obs_Tag_Recap[ry,i,r,,1] <- tmp_rec_f$Nage_at_recovery
      if(nrow(tmp_rec_m) > 0) mlt_rg_sable_data$Obs_Tag_Recap[ry,i,r,,2] <- tmp_rec_m$Nage_at_recovery
    }
  }
}

# Catch data
mlt_rg_sable_data$ObsCatch <- array(0, dim = c(mlt_rg_sable_data$n_regions, length(mlt_rg_sable_data$years), mlt_rg_sable_data$n_fish_fleets))
mlt_rg_sable_data$ObsCatch[,,1] <- spatial_data$fixed_fishery_catch # fixed gear catch
mlt_rg_sable_data$ObsCatch[,,2] <- spatial_data$trwl_fishery_catch # trawl gear catch
mlt_rg_sable_data$ObsCatch[,1:3,2] <- 0 # input 0 for first three years
mlt_rg_sable_data$UseCatch <- array(1, c(mlt_rg_sable_data$n_regions, length(mlt_rg_sable_data$years), mlt_rg_sable_data$n_fish_fleets)) # fit catch data everywhere
mlt_rg_sable_data$UseCatch[,1:3,2] <- 0 # don't use catch for the first three years
mlt_rg_sable_data$Catch_Type <- array(1, dim = c(length(mlt_rg_sable_data$years), mlt_rg_sable_data$n_fish_fleets)) # regional catch is availiable

# Fishery Indices (not fit to)
mlt_rg_sable_data$ObsFishIdx <- array(NA, c(mlt_rg_sable_data$n_regions, length(mlt_rg_sable_data$years), mlt_rg_sable_data$n_fish_fleets))
mlt_rg_sable_data$ObsFishIdx_SE <- array(NA, c(mlt_rg_sable_data$n_regions, length(mlt_rg_sable_data$years), mlt_rg_sable_data$n_fish_fleets))
mlt_rg_sable_data$UseFishIdx <- array(0, c(mlt_rg_sable_data$n_regions, length(mlt_rg_sable_data$years), mlt_rg_sable_data$n_fish_fleets))

# Fishery Age Compositions (Joint by sex, split by region)
mlt_rg_sable_data$ObsFishAgeComps <- array(0, dim = c(mlt_rg_sable_data$n_regions, length(mlt_rg_sable_data$years),
                                                      length(mlt_rg_sable_data$ages), mlt_rg_sable_data$n_sexes, mlt_rg_sable_data$n_fish_fleets))
# loop through to populate
for(y in 1:length(mlt_rg_sable_data$years)) {
  for(r in 1:mlt_rg_sable_data$n_regions) {
    if(spatial_data$fixed_catchatage_indicator[r,y] == 1) {
      mlt_rg_sable_data$ObsFishAgeComps[r,y,,1,1] <- spatial_data$obs_fixed_catchatage[-(1:length(mlt_rg_sable_data$ages)),r,y] / sum(spatial_data$obs_fixed_catchatage[,r,y]) # females
      mlt_rg_sable_data$ObsFishAgeComps[r,y,,2,1] <- spatial_data$obs_fixed_catchatage[1:length(mlt_rg_sable_data$ages),r,y] / sum(spatial_data$obs_fixed_catchatage[,r,y]) # males
    } # end if
  } # end r loop
} # end y loop

# use indicators
mlt_rg_sable_data$UseFishAgeComps <- array(0, dim = c(mlt_rg_sable_data$n_regions, length(mlt_rg_sable_data$years), mlt_rg_sable_data$n_fish_fleets))
mlt_rg_sable_data$UseFishAgeComps[,,1] <- spatial_data$fixed_catchatage_indicator # fixed gear use indicator

# input sample size
mlt_rg_sable_data$ISS_FishAgeComps <- array(0, dim = c(mlt_rg_sable_data$n_regions, length(mlt_rg_sable_data$years), mlt_rg_sable_data$n_sexes, mlt_rg_sable_data$n_fish_fleets))
mlt_rg_sable_data$ISS_FishAgeComps[,,1,1] <- 50 # ISS for fixed gear fishery

# Fishery Length Compositions (Joint by sex, split by region)
mlt_rg_sable_data$ObsFishLenComps <- array(0, dim = c(mlt_rg_sable_data$n_regions, length(mlt_rg_sable_data$years),
                                                      length(mlt_rg_sable_data$lens), mlt_rg_sable_data$n_sexes, mlt_rg_sable_data$n_fish_fleets))
# loop through to populate
for(y in 1:length(mlt_rg_sable_data$years)) {
  for(r in 1:mlt_rg_sable_data$n_regions) {

    # fixed gear
    if(spatial_data$fixed_catchatlgth_indicator[r,y] == 1) {
      mlt_rg_sable_data$ObsFishLenComps[r,y,,1,1] <- spatial_data$obs_fixed_catchatlgth[-(1:length(mlt_rg_sable_data$ages)),r,y] / sum(spatial_data$obs_fixed_catchatlgth[,r,y]) # females
      mlt_rg_sable_data$ObsFishLenComps[r,y,,2,1] <- spatial_data$obs_fixed_catchatlgth[1:length(mlt_rg_sable_data$ages),r,y] / sum(spatial_data$obs_fixed_catchatlgth[,r,y]) # males
    } # end if

    # trawl gear
    if(spatial_data$trwl_catchatlgth_indicator[r,y] == 1) {
      mlt_rg_sable_data$ObsFishLenComps[r,y,,1,2] <- spatial_data$obs_trwl_catchatlgth[-(1:length(mlt_rg_sable_data$ages)),r,y] / sum(spatial_data$obs_trwl_catchatlgth[,r,y]) # females
      mlt_rg_sable_data$ObsFishLenComps[r,y,,2,2] <- spatial_data$obs_trwl_catchatlgth[1:length(mlt_rg_sable_data$ages),r,y] / sum(spatial_data$obs_trwl_catchatlgth[,r,y]) # males
    } # end if

  } # end r loop
} # end y loop

# use indicators
mlt_rg_sable_data$UseFishLenComps <- array(0, dim = c(mlt_rg_sable_data$n_regions, length(mlt_rg_sable_data$years), mlt_rg_sable_data$n_fish_fleets))
mlt_rg_sable_data$UseFishLenComps[,,1] <- spatial_data$fixed_catchatlgth_indicator # fixed gear use indicator
mlt_rg_sable_data$UseFishLenComps[,,2] <- spatial_data$trwl_catchatlgth_indicator # trawl gear use indicator

# input sample size
mlt_rg_sable_data$ISS_FishLenComps <- array(0, dim = c(mlt_rg_sable_data$n_regions, length(mlt_rg_sable_data$years), mlt_rg_sable_data$n_sexes, mlt_rg_sable_data$n_fish_fleets))
mlt_rg_sable_data$ISS_FishLenComps[,,1,1] <- 40 # fixed gear lens ISS
mlt_rg_sable_data$ISS_FishLenComps[,,1,2] <- 40 # trawl gear lens ISS

# Survey Indices (fleet 1 = japanese, fleet 2 = domestic)
mlt_rg_sable_data$ObsSrvIdx <- array(0, dim = c(mlt_rg_sable_data$n_regions, length(mlt_rg_sable_data$years), mlt_rg_sable_data$n_srv_fleets))
mlt_rg_sable_data$ObsSrvIdx <- spatial_data$obs_srv_bio
mlt_rg_sable_data$ObsSrvIdx_SE <- spatial_data$obs_srv_se

mlt_rg_sable_data$UseSrvIdx <- array(0, dim = c(mlt_rg_sable_data$n_regions, length(mlt_rg_sable_data$years), mlt_rg_sable_data$n_srv_fleets))
mlt_rg_sable_data$UseSrvIdx <- spatial_data$srv_bio_indicator # survey indicator

# Survey Age Comps
mlt_rg_sable_data$ObsSrvAgeComps <- array(0, dim = c(mlt_rg_sable_data$n_regions, length(mlt_rg_sable_data$years),
                                                     length(mlt_rg_sable_data$ages), mlt_rg_sable_data$n_sexes, mlt_rg_sable_data$n_srv_fleets))

# loop through to populate
for(y in 1:length(mlt_rg_sable_data$years)) {
  for(r in 1:mlt_rg_sable_data$n_regions) {

    if(spatial_data$srv_catchatage_indicator[r,y,1] == 1) { # coop survey
      mlt_rg_sable_data$ObsSrvAgeComps[r,y,,1,1] <- spatial_data$obs_srv_catchatage[-(1:length(mlt_rg_sable_data$ages)),r,y,1] / sum(spatial_data$obs_srv_catchatage[,r,y,1]) # females
      mlt_rg_sable_data$ObsSrvAgeComps[r,y,,2,1] <- spatial_data$obs_srv_catchatage[1:length(mlt_rg_sable_data$ages),r,y,1] / sum(spatial_data$obs_srv_catchatage[,r,y,1]) # males
    } # end if

    if(spatial_data$srv_catchatage_indicator[r,y,2] == 1) { # domestic survey
      mlt_rg_sable_data$ObsSrvAgeComps[r,y,,1,2] <- spatial_data$obs_srv_catchatage[-(1:length(mlt_rg_sable_data$ages)),r,y,2] / sum(spatial_data$obs_srv_catchatage[,r,y,2]) # females
      mlt_rg_sable_data$ObsSrvAgeComps[r,y,,2,2] <- spatial_data$obs_srv_catchatage[1:length(mlt_rg_sable_data$ages),r,y,2] / sum(spatial_data$obs_srv_catchatage[,r,y,2]) # males
    } # end if

  } # end r loop
} # end y loop

# use indicators
mlt_rg_sable_data$UseSrvAgeComps <- array(0, dim = c(mlt_rg_sable_data$n_regions, length(mlt_rg_sable_data$years), mlt_rg_sable_data$n_srv_fleets))
mlt_rg_sable_data$UseSrvAgeComps <- spatial_data$srv_catchatage_indicator

# input sample sizes
mlt_rg_sable_data$ISS_SrvAgeComps <- array(0, dim = c(mlt_rg_sable_data$n_regions, length(mlt_rg_sable_data$years),
                                                      mlt_rg_sable_data$n_sexes, mlt_rg_sable_data$n_srv_fleets))
mlt_rg_sable_data$ISS_SrvAgeComps[,,1,] <- 60 # input sample size of 60 for all survey age comps

# Survey length compositions (not used)
mlt_rg_sable_data$ObsSrvLenComps <- array(NA, dim = c(mlt_rg_sable_data$n_regions,length(mlt_rg_sable_data$years),
                                                      length(mlt_rg_sable_data$lens), mlt_rg_sable_data$n_sexes, mlt_rg_sable_data$n_srv_fleets))
mlt_rg_sable_data$UseSrvLenComps <- array(0, dim = c(mlt_rg_sable_data$n_regions, length(mlt_rg_sable_data$years), mlt_rg_sable_data$n_srv_fleets))
mlt_rg_sable_data$ISS_SrvLenComps <- array(0, dim = c(mlt_rg_sable_data$n_regions,length(mlt_rg_sable_data$years), mlt_rg_sable_data$n_sexes, mlt_rg_sable_data$n_srv_fleets))

mlt_rg_sable_data$do_recruits_move <- 0 # recruit's don't move

# Write data
usethis::use_data(mlt_rg_sable_data, internal = FALSE, overwrite = TRUE)

# Five region sablefish report
mlt_rg_sable_rep <- readRDS(here("dev", "dev_output", "5_Region_Model_Sablefish", "rep.RDS"))
usethis::use_data(mlt_rg_sable_rep, internal = FALSE, overwrite = TRUE, compress = 'xz')

### Three Region Spatial Sablefish ------------------------------------------
# Load in data from previous spatial model and build this
spatial_data <- readRDS(here("dev", 'dev_data', "Spatial Sablefish Model", "data_3_area.RDS"))
ageing_dat <- dget(here("dev", 'dev_data','2023 Base (23.5)_final model', 'test.rdat')) # for getting ageing error

# load in tagging data
tag_rel <- readRDS(here("dev", 'dev_data', "Spatial Sablefish Model", "Tag_release_summarised_3_area.RDS"))
tag_rec <- readRDS(here("dev", 'dev_data', "Spatial Sablefish Model", "Tag_recovery_summarised_3_area.RDS"))

# create empty list to populate
three_rg_sable_data <- list()

# dimensions
three_rg_sable_data$years <- spatial_data$years
three_rg_sable_data$ages <- spatial_data$ages
three_rg_sable_data$lens <- spatial_data$length_bins
three_rg_sable_data$n_regions <- 3
three_rg_sable_data$n_sexes <- 2
three_rg_sable_data$n_fish_fleets <- 2
three_rg_sable_data$n_srv_fleets <- 2

# biologicals
# Weight at age
three_rg_sable_data$WAA <- array(0, dim = c(three_rg_sable_data$n_regions,length(three_rg_sable_data$years),length(three_rg_sable_data$ages), three_rg_sable_data$n_sexes))

# Maturity at age
three_rg_sable_data$MatAA <- array(0, dim = c(three_rg_sable_data$n_regions,length(three_rg_sable_data$years),length(three_rg_sable_data$ages), three_rg_sable_data$n_sexes))

for(r in 1:three_rg_sable_data$n_regions) {
  # weight at age
  three_rg_sable_data$WAA[r,,,1] <- t(spatial_data$female_mean_weight_by_age) # female weight at age
  three_rg_sable_data$WAA[r,,,2] <- t(spatial_data$male_mean_weight_by_age) # male weight at age

  # maturity
  three_rg_sable_data$MatAA[r,,,1] <- t(spatial_data$maturity) # female maturity at age
  three_rg_sable_data$MatAA[r,,,2] <- t(spatial_data$maturity) # male maturity at age (not used)
}

# ageing error
three_rg_sable_data$AgeingError <- as.matrix(ageing_dat$age_error) # ageing error matrix

# Size Age transition matrix
three_rg_sable_data$SizeAgeTrans <- array(0, dim = c(three_rg_sable_data$n_regions, length(three_rg_sable_data$years), length(three_rg_sable_data$lens),
                                                     length(three_rg_sable_data$ages), three_rg_sable_data$n_sexes)) # size age transition matrix
for(y in 1:length(three_rg_sable_data$years)) {
  for(r in 1:three_rg_sable_data$n_regions) {
    three_rg_sable_data$SizeAgeTrans[r,y,,,1] <- t(spatial_data$female_age_length_transition[,,y])
    three_rg_sable_data$SizeAgeTrans[r,y,,,1] <- apply(three_rg_sable_data$SizeAgeTrans[r,y,,,1], 2, function(x) x / sum(x))
    three_rg_sable_data$SizeAgeTrans[r,y,,,2] <- t(spatial_data$male_age_length_transition[,,y]) / colSums(t(spatial_data$male_age_length_transition[,,y]))
    three_rg_sable_data$SizeAgeTrans[r,y,,,2] <- apply(three_rg_sable_data$SizeAgeTrans[r,y,,,2], 2, function(x) x / sum(x))
  } # end r loop
} # end y loop

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
  tmp <- data.frame(regions = rel_cohorts_df[i,4], tag_yrs = rel_cohorts_df[i,2]-1959)
  tag_release_ind <- rbind(tag_release_ind, tmp)
} # end i

three_rg_sable_data$tag_release_indicator <- as.matrix(tag_release_ind) # input releases of cohorts in
three_rg_sable_data$n_tag_cohorts <- nrow(three_rg_sable_data$tag_release_indicator) # number of tag cohorts
three_rg_sable_data$max_tag_liberty <- 15 # maximum liberty to track cohorts

# Get Tagged Fish
three_rg_sable_data$Tagged_Fish <- array(0, dim = c(three_rg_sable_data$n_tag_cohorts, length(three_rg_sable_data$ages), three_rg_sable_data$n_sexes)) # tagged fish
for(i in 1:nrow(three_rg_sable_data$tag_release_indicator )) {
  tmp_rel_f <- tag_rel %>% filter(release_event_id == i, sex == "F") # filter to a given cohort females
  tmp_rel_m <- tag_rel %>% filter(release_event_id == i, sex == "M") # filter to a given cohort males
  three_rg_sable_data$Tagged_Fish[i,,1] <- tmp_rel_f$Nage_at_release # females
  three_rg_sable_data$Tagged_Fish[i,,2] <- tmp_rel_m$Nage_at_release # males
} # end i for tag cohorts


# Get Observed Recaptures
three_rg_sable_data$Obs_Tag_Recap <- array(0, dim = c(three_rg_sable_data$max_tag_liberty, three_rg_sable_data$n_tag_cohorts, three_rg_sable_data$n_regions,
                                                      length(three_rg_sable_data$ages), three_rg_sable_data$n_sexes))

# pre filter data frame before looping
tag_rec_f <- tag_rec %>% filter(sex == "F")
tag_rec_m <- tag_rec %>% filter(sex == "M")

for(ry in 1:three_rg_sable_data$max_tag_liberty) {
  for(i in 1:three_rg_sable_data$n_tag_cohorts) {
    for(r in 1:three_rg_sable_data$n_regions) {
      tmp_rec_f <- tag_rec_f %>% filter(release_event_id == i, region_num == r, tag_lib == ry)
      tmp_rec_m <- tag_rec_m %>% filter(release_event_id == i, region_num == r, tag_lib == ry)
      if(nrow(tmp_rec_f) > 0) three_rg_sable_data$Obs_Tag_Recap[ry,i,r,,1] <- tmp_rec_f$Nage_at_recovery
      if(nrow(tmp_rec_m) > 0) three_rg_sable_data$Obs_Tag_Recap[ry,i,r,,2] <- tmp_rec_m$Nage_at_recovery
    }
  }
}

# Catch data
three_rg_sable_data$ObsCatch <- array(0, dim = c(three_rg_sable_data$n_regions, length(three_rg_sable_data$years), three_rg_sable_data$n_fish_fleets))
three_rg_sable_data$ObsCatch[,,1] <- spatial_data$fixed_fishery_catch # fixed gear catch
three_rg_sable_data$ObsCatch[,,2] <- spatial_data$trwl_fishery_catch # trawl gear catch
three_rg_sable_data$ObsCatch[,1:3,2] <- 0 # input 0 for first three years
three_rg_sable_data$UseCatch <- array(1, c(three_rg_sable_data$n_regions, length(three_rg_sable_data$years), three_rg_sable_data$n_fish_fleets)) # fit catch data everywhere
three_rg_sable_data$UseCatch[,1:3,2] <- 0 # don't use catch for the first three years
three_rg_sable_data$Catch_Type <- array(1, dim = c(length(three_rg_sable_data$years), three_rg_sable_data$n_fish_fleets)) # regional catch is availiable

# Fishery Indices (not fit to)
three_rg_sable_data$ObsFishIdx <- array(NA, c(three_rg_sable_data$n_regions, length(three_rg_sable_data$years), three_rg_sable_data$n_fish_fleets))
three_rg_sable_data$ObsFishIdx_SE <- array(NA, c(three_rg_sable_data$n_regions, length(three_rg_sable_data$years), three_rg_sable_data$n_fish_fleets))
three_rg_sable_data$UseFishIdx <- array(0, c(three_rg_sable_data$n_regions, length(three_rg_sable_data$years), three_rg_sable_data$n_fish_fleets))

# Fishery Age Compositions (Joint by sex, split by region)
three_rg_sable_data$ObsFishAgeComps <- array(0, dim = c(three_rg_sable_data$n_regions, length(three_rg_sable_data$years),
                                                        length(three_rg_sable_data$ages), three_rg_sable_data$n_sexes, three_rg_sable_data$n_fish_fleets))
# loop through to populate
for(y in 1:length(three_rg_sable_data$years)) {
  for(r in 1:three_rg_sable_data$n_regions) {
    if(spatial_data$fixed_catchatage_indicator[r,y] == 1) {
      three_rg_sable_data$ObsFishAgeComps[r,y,,1,1] <- spatial_data$obs_fixed_catchatage[-(1:length(three_rg_sable_data$ages)),r,y] / sum(spatial_data$obs_fixed_catchatage[,r,y]) # females
      three_rg_sable_data$ObsFishAgeComps[r,y,,2,1] <- spatial_data$obs_fixed_catchatage[1:length(three_rg_sable_data$ages),r,y] / sum(spatial_data$obs_fixed_catchatage[,r,y]) # males
    } # end if
  } # end r loop
} # end y loop

# use indicators
three_rg_sable_data$UseFishAgeComps <- array(0, dim = c(three_rg_sable_data$n_regions, length(three_rg_sable_data$years), three_rg_sable_data$n_fish_fleets))
three_rg_sable_data$UseFishAgeComps[,,1] <- spatial_data$fixed_catchatage_indicator # fixed gear use indicator

# input sample size
three_rg_sable_data$ISS_FishAgeComps <- array(0, dim = c(three_rg_sable_data$n_regions, length(three_rg_sable_data$years), three_rg_sable_data$n_sexes, three_rg_sable_data$n_fish_fleets))
three_rg_sable_data$ISS_FishAgeComps[,,1,1] <- 50 # ISS for fixed gear fishery

# Fishery Length Compositions (Joint by sex, split by region)
three_rg_sable_data$ObsFishLenComps <- array(0, dim = c(three_rg_sable_data$n_regions, length(three_rg_sable_data$years),
                                                        length(three_rg_sable_data$lens), three_rg_sable_data$n_sexes, three_rg_sable_data$n_fish_fleets))
# loop through to populate
for(y in 1:length(three_rg_sable_data$years)) {
  for(r in 1:three_rg_sable_data$n_regions) {

    # fixed gear
    if(spatial_data$fixed_catchatlgth_indicator[r,y] == 1) {
      three_rg_sable_data$ObsFishLenComps[r,y,,1,1] <- spatial_data$obs_fixed_catchatlgth[-(1:length(three_rg_sable_data$ages)),r,y] / sum(spatial_data$obs_fixed_catchatlgth[,r,y]) # females
      three_rg_sable_data$ObsFishLenComps[r,y,,2,1] <- spatial_data$obs_fixed_catchatlgth[1:length(three_rg_sable_data$ages),r,y] / sum(spatial_data$obs_fixed_catchatlgth[,r,y]) # males
    } # end ifM

    # trawl gear
    if(spatial_data$trwl_catchatlgth_indicator[r,y] == 1) {
      three_rg_sable_data$ObsFishLenComps[r,y,,1,2] <- spatial_data$obs_trwl_catchatlgth[-(1:length(three_rg_sable_data$ages)),r,y] / sum(spatial_data$obs_trwl_catchatlgth[,r,y]) # females
      three_rg_sable_data$ObsFishLenComps[r,y,,2,2] <- spatial_data$obs_trwl_catchatlgth[1:length(three_rg_sable_data$ages),r,y] / sum(spatial_data$obs_trwl_catchatlgth[,r,y]) # males
    } # end if

  } # end r loop
} # end y loop

# use indicators
three_rg_sable_data$UseFishLenComps <- array(0, dim = c(three_rg_sable_data$n_regions, length(three_rg_sable_data$years), three_rg_sable_data$n_fish_fleets))
three_rg_sable_data$UseFishLenComps[,,1] <- spatial_data$fixed_catchatlgth_indicator # fixed gear use indicator
three_rg_sable_data$UseFishLenComps[,,2] <- spatial_data$trwl_catchatlgth_indicator # trawl gear use indicator

# input sample size
three_rg_sable_data$ISS_FishLenComps <- array(0, dim = c(three_rg_sable_data$n_regions, length(three_rg_sable_data$years), three_rg_sable_data$n_sexes, three_rg_sable_data$n_fish_fleets))
three_rg_sable_data$ISS_FishLenComps[,,1,1] <- 40 # fixed gear lens ISS
three_rg_sable_data$ISS_FishLenComps[,,1,2] <- 40 # trawl gear lens ISS

# Survey Indices (fleet 1 = japanese, fleet 2 = domestic)
three_rg_sable_data$ObsSrvIdx <- array(0, dim = c(three_rg_sable_data$n_regions, length(three_rg_sable_data$years), three_rg_sable_data$n_srv_fleets))
three_rg_sable_data$ObsSrvIdx <- spatial_data$obs_srv_bio
three_rg_sable_data$ObsSrvIdx_SE <- spatial_data$obs_srv_se

three_rg_sable_data$UseSrvIdx <- array(0, dim = c(three_rg_sable_data$n_regions, length(three_rg_sable_data$years), three_rg_sable_data$n_srv_fleets))
three_rg_sable_data$UseSrvIdx <- spatial_data$srv_bio_indicator # survey indicator

# Survey Age Comps
three_rg_sable_data$ObsSrvAgeComps <- array(0, dim = c(three_rg_sable_data$n_regions, length(three_rg_sable_data$years),
                                                       length(three_rg_sable_data$ages), three_rg_sable_data$n_sexes, three_rg_sable_data$n_srv_fleets))

# loop through to populate
for(y in 1:length(three_rg_sable_data$years)) {
  for(r in 1:three_rg_sable_data$n_regions) {

    if(spatial_data$srv_catchatage_indicator[r,y,1] == 1) { # coop survey
      three_rg_sable_data$ObsSrvAgeComps[r,y,,1,1] <- spatial_data$obs_srv_catchatage[-(1:length(three_rg_sable_data$ages)),r,y,1] / sum(spatial_data$obs_srv_catchatage[,r,y,1]) # females
      three_rg_sable_data$ObsSrvAgeComps[r,y,,2,1] <- spatial_data$obs_srv_catchatage[1:length(three_rg_sable_data$ages),r,y,1] / sum(spatial_data$obs_srv_catchatage[,r,y,1]) # males
    } # end if

    if(spatial_data$srv_catchatage_indicator[r,y,2] == 1) { # domestic survey
      three_rg_sable_data$ObsSrvAgeComps[r,y,,1,2] <- spatial_data$obs_srv_catchatage[-(1:length(three_rg_sable_data$ages)),r,y,2] / sum(spatial_data$obs_srv_catchatage[,r,y,2]) # females
      three_rg_sable_data$ObsSrvAgeComps[r,y,,2,2] <- spatial_data$obs_srv_catchatage[1:length(three_rg_sable_data$ages),r,y,2] / sum(spatial_data$obs_srv_catchatage[,r,y,2]) # males
    } # end if

  } # end r loop
} # end y loop

# use indicators
three_rg_sable_data$UseSrvAgeComps <- array(0, dim = c(three_rg_sable_data$n_regions, length(three_rg_sable_data$years), three_rg_sable_data$n_srv_fleets))
three_rg_sable_data$UseSrvAgeComps <- spatial_data$srv_catchatage_indicator

# input sample sizes
three_rg_sable_data$ISS_SrvAgeComps <- array(0, dim = c(three_rg_sable_data$n_regions, length(three_rg_sable_data$years),
                                                        three_rg_sable_data$n_sexes, three_rg_sable_data$n_srv_fleets))
three_rg_sable_data$ISS_SrvAgeComps[,,1,] <- 60 # input sample size of 60 for all survey age comps

# Survey length compositions (not used)
three_rg_sable_data$ObsSrvLenComps <- array(NA, dim = c(three_rg_sable_data$n_regions,length(three_rg_sable_data$years),
                                                        length(three_rg_sable_data$lens), three_rg_sable_data$n_sexes, three_rg_sable_data$n_srv_fleets))
three_rg_sable_data$UseSrvLenComps <- array(0, dim = c(three_rg_sable_data$n_regions, length(three_rg_sable_data$years), three_rg_sable_data$n_srv_fleets))
three_rg_sable_data$ISS_SrvLenComps <- array(0, dim = c(three_rg_sable_data$n_regions,length(three_rg_sable_data$years), three_rg_sable_data$n_sexes, three_rg_sable_data$n_srv_fleets))

three_rg_sable_data$do_recruits_move <- 0 # recruit's don't move

# Write data
usethis::use_data(three_rg_sable_data, internal = FALSE, overwrite = TRUE)
# Three region sablefish report
three_rg_sable_rep <- readRDS(here("dev", "dev_output", "3_Region_Model_Sablefish", "rep.RDS"))
usethis::use_data(three_rg_sable_rep, internal = FALSE, overwrite = TRUE, compress = 'xz')

# Single Region Dusky -----------------------------------------------------

# Read in data file
in_dat <- readLines(here::here("dev", "dev_data", "goa_dusk_2024.dat"))
# Read in report file
rep_dat <- readLines(here::here("dev", "dev_data", "dusky.rep"))
# Read in parameters and stuff for comparison
par_dat <- R2admb::read_pars(here::here("dev", "dev_data", "dusky"))
rep_dat_rtem <- PBSmodelling::readList(here::here("dev", "dev_data", "dusky_rtem.rep"))

sgl_rg_dusky_data <- list()

# Number of regions
sgl_rg_dusky_data$n_regions <- 1

# Number of sexes
sgl_rg_dusky_data$n_sexes <- 1

# Number of fishery fleets
sgl_rg_dusky_data$n_fish_fleets <- 1

# Number of survey fleets
sgl_rg_dusky_data$n_srv_fleets <- 1

# Get model start and end year
sgl_rg_dusky_data$years <- in_dat[[18]]:in_dat[[20]]
n_years <- length(sgl_rg_dusky_data$years)

# Get model ages (not the same as data ages)
obs_ages <- in_dat[[22]]:(in_dat[[26]])
sgl_rg_dusky_data$mod_ages <- c(obs_ages, 31:33) # add more ages because of model ages
n_ages <- length(sgl_rg_dusky_data$mod_ages)
n_obs_ages <- length(obs_ages)

# Get length bins
sgl_rg_dusky_data$lens <- as.numeric(unlist(strsplit(in_dat[[34]], " ")))
n_lens <- length(sgl_rg_dusky_data$lens)

# Get spawning month
sgl_rg_dusky_data$spwn_month <- 0

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

# get ageing error
# Find where the age error matrix starts
start_line <- grep("age error transition matrix", in_dat, ignore.case = TRUE)
data_lines <- in_dat[(start_line + 1):313]
data_lines <- grep("^[ 0-9eE.+\\-]", data_lines, value = TRUE)
flat_values <- scan(text = paste(data_lines, collapse = "\n"), quiet = TRUE)
ncol <- n_obs_ages
nrow <- length(flat_values) / ncol
sgl_rg_dusky_data$age_error_matrix <- matrix(flat_values, ncol = ncol, byrow = TRUE)

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
sgl_rg_dusky_data$waa_arr <- array(0, dim = c(n_regions, n_years, n_ages, n_sexes))
sgl_rg_dusky_data$fix_natmort <- array(0, dim = c(n_regions, n_years, n_ages, n_sexes))
sgl_rg_dusky_data$mataa_arr <- array(0, dim = c(n_regions, n_years, n_ages, n_sexes))
sgl_rg_dusky_data$sizeage <- array(0, dim = c(n_regions, n_years, n_lens, n_ages, n_sexes))

for(r in 1:n_regions) {
  for(y in 1:n_years) {
    for(s in 1:n_sexes) {
      sgl_rg_dusky_data$waa_arr[r,y,,s] <- waa # weight at age
      sgl_rg_dusky_data$fix_natmort[r,y,,s] <- 0.07 # natural mortality
      sgl_rg_dusky_data$mataa_arr[r,y,,s] <- mat # maturity at age
      sgl_rg_dusky_data$sizeage[r,y,,,s] <- t(size_age_matrix) # size age transition
    } # end s loop
  } # end y loop
} # end r loop

# Catch inputs
sgl_rg_dusky_data$ObsCatch <- array(0, dim = c(n_regions, n_years, n_fish_fleets))
sgl_rg_dusky_data$Catch_Type <- array(1, dim = c(n_years, n_fish_fleets))
sgl_rg_dusky_data$UseCatch <- array(1, dim = c(n_regions, n_years, n_fish_fleets))
sgl_rg_dusky_data$ObsCatch[] <- fish_catch

# Fishery index inputs -- none used
sgl_rg_dusky_data$ObsFishIdx <- array(NA, dim = c(n_regions, n_years, n_fish_fleets))
sgl_rg_dusky_data$ObsFishIdx_SE <- array(NA, dim = c(n_regions, n_years, n_fish_fleets))
sgl_rg_dusky_data$UseFishIdx <- array(0, dim = c(n_regions, n_years, n_fish_fleets))

# Fishery age comps inputs
sgl_rg_dusky_data$ObsFishAgeComps <- array(0, dim = c(n_regions, n_years, n_obs_ages, n_sexes, n_fish_fleets))
sgl_rg_dusky_data$UseFishAgeComps <- array(0, dim = c(n_regions, n_years, n_fish_fleets))
sgl_rg_dusky_data$ISS_FishAgeComps <- array(0, dim = c(n_regions, n_years, n_sexes, n_fish_fleets))

sgl_rg_dusky_data$ObsFishAgeComps[1,which(sgl_rg_dusky_data$years %in% fish_catch_age_yrs),,1,1] <- oac_fish_age
sgl_rg_dusky_data$UseFishAgeComps[1,which(sgl_rg_dusky_data$years %in% fish_catch_age_yrs),1] <- 1
sgl_rg_dusky_data$ISS_FishAgeComps[1,which(sgl_rg_dusky_data$years %in% fish_catch_age_yrs),1,1] <- sqrt(fish_catch_age_nhauls * fish_catch_age_n)
sgl_rg_dusky_data$ISS_FishAgeComps[1,which(sgl_rg_dusky_data$years %in% fish_catch_age_yrs),1,1] <- rep_dat_rtem$oac_fish_sample[,1]

# fishery length comps inputs
sgl_rg_dusky_data$ObsFishLenComps <- array(0, dim = c(n_regions, n_years, n_lens, n_sexes, n_fish_fleets))
sgl_rg_dusky_data$UseFishLenComps <- array(0, dim = c(n_regions, n_years, n_fish_fleets))
sgl_rg_dusky_data$ISS_FishLenComps <- array(0, dim = c(n_regions, n_years, n_sexes, n_fish_fleets))

sgl_rg_dusky_data$ObsFishLenComps[1,which(sgl_rg_dusky_data$years %in% fish_catch_size_yrs),,1,1] <- oac_fish_size
sgl_rg_dusky_data$UseFishLenComps[1,which(sgl_rg_dusky_data$years %in% fish_catch_size_yrs),1] <- 1
sgl_rg_dusky_data$ISS_FishLenComps[1,which(sgl_rg_dusky_data$years %in% fish_catch_size_yrs),1,1] <- sqrt(fish_catch_size_nhauls * fish_catch_size_n)
sgl_rg_dusky_data$ISS_FishLenComps[1,which(sgl_rg_dusky_data$years %in% fish_catch_size_yrs),1,1] <- rep_dat_rtem$osc_fish_sample[,1]

# survey index
sgl_rg_dusky_data$ObsSrvIdx <- array(0, dim = c(n_regions, n_years, sgl_rg_dusky_data$n_srv_fleets))
sgl_rg_dusky_data$ObsSrvIdx_SE <- array(0, dim = c(n_regions, n_years, sgl_rg_dusky_data$n_srv_fleets))
sgl_rg_dusky_data$UseSrvIdx <- array(0, dim = c(n_regions, n_years, sgl_rg_dusky_data$n_srv_fleets))

sgl_rg_dusky_data$ObsSrvIdx[1,which(sgl_rg_dusky_data$years %in% bts_srv_yrs),1] <- bts_srv_idx

# transform to lognormal sd to be the same as dusky assessment
cv <- bts_srv_se / bts_srv_idx
se <- bts_srv_idx * sqrt(log(1 + cv^2))
sgl_rg_dusky_data$ObsSrvIdx_SE[1,which(sgl_rg_dusky_data$years %in% bts_srv_yrs),1] <- se
sgl_rg_dusky_data$UseSrvIdx[1,which(sgl_rg_dusky_data$years %in% bts_srv_yrs),1] <- 1

# survey age comps inputs
sgl_rg_dusky_data$ObsSrvAgeComps <- array(0, dim = c(n_regions, n_years, n_obs_ages, n_sexes, sgl_rg_dusky_data$n_srv_fleets))
sgl_rg_dusky_data$UseSrvAgeComps <- array(0, dim = c(n_regions, n_years, sgl_rg_dusky_data$n_srv_fleets))
sgl_rg_dusky_data$ISS_SrvAgeComps <- array(0, dim = c(n_regions, n_years, n_sexes, sgl_rg_dusky_data$n_srv_fleets))

sgl_rg_dusky_data$ObsSrvAgeComps[1,which(sgl_rg_dusky_data$years %in% bts_catch_age_yrs),,1,1] <- oac_srv1_age
sgl_rg_dusky_data$UseSrvAgeComps[1,which(sgl_rg_dusky_data$years %in% bts_catch_age_yrs),1] <- 1
sgl_rg_dusky_data$ISS_SrvAgeComps[1,which(sgl_rg_dusky_data$years %in% bts_catch_age_yrs),1,1] <- sqrt(bts_catch_age_nhauls * bts_catch_age_n)
sgl_rg_dusky_data$ISS_SrvAgeComps[1,which(sgl_rg_dusky_data$years %in% bts_catch_age_yrs),1,1] <- rep_dat_rtem$oac_srv1_sample[,1]

# survey length comps inputs
sgl_rg_dusky_data$ObsSrvLenComps <- array(0, dim = c(n_regions, n_years, n_lens, n_sexes, sgl_rg_dusky_data$n_srv_fleets))
sgl_rg_dusky_data$UseSrvLenComps <- array(0, dim = c(n_regions, n_years, sgl_rg_dusky_data$n_srv_fleets))
sgl_rg_dusky_data$ISS_SrvLenComps <- array(0, dim = c(n_regions, n_years, n_sexes, sgl_rg_dusky_data$n_srv_fleets))

sgl_rg_dusky_data$ObsSrvLenComps[1,which(sgl_rg_dusky_data$years %in% bts_catch_size_yrs),,1,1] <- oac_srv1_size
sgl_rg_dusky_data$UseSrvLenComps[1,which(sgl_rg_dusky_data$years %in% bts_catch_size_yrs),1] <- 0
sgl_rg_dusky_data$ISS_SrvLenComps[1,which(sgl_rg_dusky_data$years %in% bts_catch_size_yrs),1,1] <- 0

# Write data
usethis::use_data(sgl_rg_dusky_data, internal = FALSE, overwrite = TRUE)

