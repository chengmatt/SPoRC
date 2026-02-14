library(here)

# Single Region Sablefish -------------------------------------------------
# Load in data and build
sgl_rg_sable_data <- readRDS(here("dev", 'dev_data', "sgl_rg_sable_data.RDS")) # read in RTMB data
ageing_dat <- dget(here("dev", 'dev_data','2023 Base (23.5)_final model', 'test.rdat')) # for getting ageing error

# Set up dimensions
sgl_rg_sable_data$years <- 1:65
sgl_rg_sable_data$ages <- 1:30
sgl_rg_sable_data$n_fish_fleets <- 2
sgl_rg_sable_data$n_srv_fleets <- 3
sgl_rg_sable_data$n_sexes <- 2
sgl_rg_sable_data$n_seas <- 1
sgl_rg_sable_data$n_regions <- 1
sgl_rg_sable_data$spawn_seas <- 1
sgl_rg_sable_data$seasdur <- 1
sgl_rg_sable_data$lens <- seq(41,99,2)

# Dimensions to use in current R file
n_regions <- sgl_rg_sable_data$n_regions
n_yrs <- length(sgl_rg_sable_data$years)
n_ages <- length(sgl_rg_sable_data$ages)
n_seas <- sgl_rg_sable_data$n_seas
n_sexes <- sgl_rg_sable_data$n_sexes
n_fish_fleets <- sgl_rg_sable_data$n_fish_fleets
n_srv_fleets <- sgl_rg_sable_data$n_srv_fleets
n_lens <- length(sgl_rg_sable_data$lens)


# get ageing error
sgl_rg_sable_data$age_error <- ageing_dat$age_error # input ageing error

# Expand arrays
sgl_rg_sable_data$WAA <- array(sgl_rg_sable_data$WAA, dim = c(n_regions, n_yrs, n_seas, n_ages, n_sexes))
sgl_rg_sable_data$MatAA <- array(sgl_rg_sable_data$MatAA, dim = c(n_regions, n_yrs, n_seas, n_ages, n_sexes))
sgl_rg_sable_data$SizeAgeTrans <- array(sgl_rg_sable_data$SizeAgeTrans, dim = c(n_regions, n_yrs, n_seas, n_lens, n_ages, n_sexes))
sgl_rg_sable_data$ObsCatch <- array(sgl_rg_sable_data$ObsCatch, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))
sgl_rg_sable_data$UseCatch <- array(sgl_rg_sable_data$UseCatch, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))
sgl_rg_sable_data$ObsFishIdx <- array(sgl_rg_sable_data$ObsFishIdx, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))
sgl_rg_sable_data$ObsFishIdx_SE <- array(sgl_rg_sable_data$ObsFishIdx_SE, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))
sgl_rg_sable_data$UseFishIdx <- array(sgl_rg_sable_data$UseFishIdx, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))
sgl_rg_sable_data$ObsSrvIdx <- array(sgl_rg_sable_data$ObsSrvIdx, dim = c(n_regions, n_yrs, n_seas, n_srv_fleets))
sgl_rg_sable_data$ObsSrvIdx_SE <- array(sgl_rg_sable_data$ObsSrvIdx_SE, dim = c(n_regions, n_yrs, n_seas, n_srv_fleets))
sgl_rg_sable_data$UseSrvIdx <- array(sgl_rg_sable_data$UseSrvIdx, dim = c(n_regions, n_yrs, n_seas, n_srv_fleets))
sgl_rg_sable_data$ObsFishAgeComps <- array(sgl_rg_sable_data$ObsFishAgeComps, dim = c(n_regions, n_yrs, n_seas, n_ages, n_sexes, n_fish_fleets))
sgl_rg_sable_data$UseFishAgeComps <- array(sgl_rg_sable_data$UseFishAgeComps, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))
sgl_rg_sable_data$ISS_FishAgeComps <- array(sgl_rg_sable_data$ISS_FishAgeComps, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_fish_fleets))
sgl_rg_sable_data$ObsFishLenComps <- array(sgl_rg_sable_data$ObsFishLenComps, dim = c(n_regions, n_yrs, n_seas, n_lens, n_sexes, n_fish_fleets))
sgl_rg_sable_data$UseFishLenComps <- array(sgl_rg_sable_data$UseFishLenComps, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))
sgl_rg_sable_data$ISS_FishLenComps <- array(sgl_rg_sable_data$ISS_FishLenComps, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_fish_fleets))
sgl_rg_sable_data$ObsSrvAgeComps <- array(sgl_rg_sable_data$ObsSrvAgeComps, dim = c(n_regions, n_yrs, n_seas, n_ages, n_sexes, n_srv_fleets))
sgl_rg_sable_data$UseSrvAgeComps <- array(sgl_rg_sable_data$UseSrvAgeComps, dim = c(n_regions, n_yrs, n_seas, n_srv_fleets))
sgl_rg_sable_data$ISS_SrvAgeComps <- array(sgl_rg_sable_data$ISS_SrvAgeComps, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_srv_fleets))
sgl_rg_sable_data$ObsSrvLenComps <- array(sgl_rg_sable_data$ObsSrvLenComps, dim = c(n_regions, n_yrs, n_seas, n_lens, n_sexes, n_srv_fleets))
sgl_rg_sable_data$UseSrvLenComps <- array(sgl_rg_sable_data$UseSrvLenComps, dim = c(n_regions, n_yrs, n_seas, n_srv_fleets))
sgl_rg_sable_data$ISS_SrvLenComps <- array(sgl_rg_sable_data$ISS_SrvLenComps, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_srv_fleets))

# get admb report and plug into rtmb data
tem_dat <- dget(here("dev", 'dev_data', '2024 Base (23.5)_final model_v3', 'tem.rdat'))
sgl_rg_sable_data$admb_recr <- tem_dat$t.series$Recr # recruitment
sgl_rg_sable_data$admb_spbiom <- tem_dat$t.series$spbiom # ssb
usethis::use_data(sgl_rg_sable_data, internal = FALSE, overwrite = TRUE)
