library(here)

# Single Region Dusky -----------------------------------------------------

# Read in data file
in_dat <- readLines(here::here("dev", "dev_data", "goa_dusk_2024.dat"))
# Read in report file
rep_dat <- readLines(here::here("dev", "dev_data", "dusky.rep"))
# Read in parameters and stuff for comparison
par_dat <- R2admb::read_pars(here::here("dev", "dev_data", "dusky"))
rep_dat_rtem <- PBSmodelling::readList(here::here("dev", "dev_data", "dusky_rtem.rep"))

sgl_rg_dusky_data <- list()

# Number of populations
sgl_rg_dusky_data$n_pop <- 1
n_pop <- 1

# Number of regions
sgl_rg_dusky_data$n_regions <- 1
n_regions <- 1

# Number of sexes
sgl_rg_dusky_data$n_sexes <- 1

# Number of fishery fleets
sgl_rg_dusky_data$n_fish_fleets <- 1

# Number of survey fleets
sgl_rg_dusky_data$n_srv_fleets <- 1

# Number of seasons
sgl_rg_dusky_data$n_seas <- 1
sgl_rg_dusky_data$seasdur <- 1

n_seas <- 1

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
sgl_rg_dusky_data$waa_arr <- array(0, dim = c(n_pop, n_regions, n_years, n_seas, n_ages, n_sexes))
sgl_rg_dusky_data$fix_natmort <- array(0, dim = c(n_pop,n_regions, n_years, n_ages, n_sexes))
sgl_rg_dusky_data$mataa_arr <- array(0, dim = c(n_pop,n_regions, n_years, n_seas, n_ages, n_sexes))
sgl_rg_dusky_data$sizeage <- array(0, dim = c(n_pop,n_regions, n_years, n_seas, n_lens, n_ages, n_sexes))

for(r in 1:n_regions) {
  for(y in 1:n_years) {
    for(s in 1:sgl_rg_dusky_data$n_sexes) {
      sgl_rg_dusky_data$fix_natmort[,r,y,,s] <- 0.07 # natural mortality
      for(seas in 1:n_seas) {
        sgl_rg_dusky_data$waa_arr[,r,y,seas,,s] <- waa # weight at age
        sgl_rg_dusky_data$mataa_arr[,r,y,seas,,s] <- mat # maturity at age
        sgl_rg_dusky_data$sizeage[,r,y,seas,,,s] <- t(size_age_matrix) # size age transition
      }
    } # end s loop
  } # end y loop
} # end r loop

# Catch inputs
sgl_rg_dusky_data$ObsCatch <- base::array(0, dim = c(n_regions, n_years, n_seas, n_fish_fleets))
sgl_rg_dusky_data$Catch_Type <- array(1, dim = c(n_years, n_fish_fleets))
sgl_rg_dusky_data$UseCatch <- array(1, dim = c(n_regions, n_years, n_seas, n_fish_fleets))
sgl_rg_dusky_data$ObsCatch[] <- fish_catch

# Fishery index inputs -- none used
sgl_rg_dusky_data$ObsFishIdx <- array(NA, dim = c(n_regions, n_years, n_seas, n_fish_fleets))
sgl_rg_dusky_data$ObsFishIdx_SE <- array(NA, dim = c(n_regions, n_years, n_seas, n_fish_fleets))
sgl_rg_dusky_data$UseFishIdx <- array(0, dim = c(n_regions, n_years, n_seas, n_fish_fleets))

# Fishery age comps inputs
sgl_rg_dusky_data$ObsFishAgeComps <- array(0, dim = c(n_regions, n_years, n_seas, n_obs_ages, n_sexes, n_fish_fleets))
sgl_rg_dusky_data$UseFishAgeComps <- array(0, dim = c(n_regions, n_years, n_seas, n_fish_fleets))
sgl_rg_dusky_data$ISS_FishAgeComps <- array(0, dim = c(n_regions, n_years, n_seas, n_sexes, n_fish_fleets))

sgl_rg_dusky_data$ObsFishAgeComps[1,which(sgl_rg_dusky_data$years %in% fish_catch_age_yrs),1,,1,1] <- oac_fish_age
sgl_rg_dusky_data$UseFishAgeComps[1,which(sgl_rg_dusky_data$years %in% fish_catch_age_yrs),1,1] <- 1
sgl_rg_dusky_data$ISS_FishAgeComps[1,which(sgl_rg_dusky_data$years %in% fish_catch_age_yrs),1,1,1] <- sqrt(fish_catch_age_nhauls * fish_catch_age_n)
sgl_rg_dusky_data$ISS_FishAgeComps[1,which(sgl_rg_dusky_data$years %in% fish_catch_age_yrs),1,1,1] <- rep_dat_rtem$oac_fish_sample[,1]

# fishery length comps inputs
sgl_rg_dusky_data$ObsFishLenComps <- array(0, dim = c(n_regions, n_years, n_seas, n_lens, n_sexes, n_fish_fleets))
sgl_rg_dusky_data$UseFishLenComps <- array(0, dim = c(n_regions, n_years, n_seas, n_fish_fleets))
sgl_rg_dusky_data$ISS_FishLenComps <- array(0, dim = c(n_regions, n_years, n_seas, n_sexes, n_fish_fleets))

sgl_rg_dusky_data$ObsFishLenComps[1,which(sgl_rg_dusky_data$years %in% fish_catch_size_yrs),1,,1,1] <- oac_fish_size
sgl_rg_dusky_data$UseFishLenComps[1,which(sgl_rg_dusky_data$years %in% fish_catch_size_yrs),1,1] <- 1
sgl_rg_dusky_data$ISS_FishLenComps[1,which(sgl_rg_dusky_data$years %in% fish_catch_size_yrs),1,1,1] <- sqrt(fish_catch_size_nhauls * fish_catch_size_n)
sgl_rg_dusky_data$ISS_FishLenComps[1,which(sgl_rg_dusky_data$years %in% fish_catch_size_yrs),1,1,1] <- rep_dat_rtem$osc_fish_sample[,1]

# survey index
sgl_rg_dusky_data$ObsSrvIdx <- array(0, dim = c(n_regions, n_years, n_seas, sgl_rg_dusky_data$n_srv_fleets))
sgl_rg_dusky_data$ObsSrvIdx_SE <- array(0, dim = c(n_regions, n_years, n_seas, sgl_rg_dusky_data$n_srv_fleets))
sgl_rg_dusky_data$UseSrvIdx <- array(0, dim = c(n_regions, n_years, n_seas, sgl_rg_dusky_data$n_srv_fleets))

sgl_rg_dusky_data$ObsSrvIdx[1,which(sgl_rg_dusky_data$years %in% bts_srv_yrs),1,1] <- bts_srv_idx

# transform to lognormal sd to be the same as dusky assessment
cv <- bts_srv_se / bts_srv_idx
se <- bts_srv_idx * sqrt(log(1 + cv^2))
sgl_rg_dusky_data$ObsSrvIdx_SE[1,which(sgl_rg_dusky_data$years %in% bts_srv_yrs),1,1] <- se
sgl_rg_dusky_data$UseSrvIdx[1,which(sgl_rg_dusky_data$years %in% bts_srv_yrs),1,1] <- 1

# survey age comps inputs
sgl_rg_dusky_data$ObsSrvAgeComps <- array(0, dim = c(n_regions, n_years, n_seas, n_obs_ages, n_sexes, sgl_rg_dusky_data$n_srv_fleets))
sgl_rg_dusky_data$UseSrvAgeComps <- array(0, dim = c(n_regions, n_years, n_seas, sgl_rg_dusky_data$n_srv_fleets))
sgl_rg_dusky_data$ISS_SrvAgeComps <- array(0, dim = c(n_regions, n_years, n_seas, n_sexes, sgl_rg_dusky_data$n_srv_fleets))

sgl_rg_dusky_data$ObsSrvAgeComps[1,which(sgl_rg_dusky_data$years %in% bts_catch_age_yrs),1,,1,1] <- oac_srv1_age
sgl_rg_dusky_data$UseSrvAgeComps[1,which(sgl_rg_dusky_data$years %in% bts_catch_age_yrs),1,1] <- 1
sgl_rg_dusky_data$ISS_SrvAgeComps[1,which(sgl_rg_dusky_data$years %in% bts_catch_age_yrs),1,1,1] <- sqrt(bts_catch_age_nhauls * bts_catch_age_n)
sgl_rg_dusky_data$ISS_SrvAgeComps[1,which(sgl_rg_dusky_data$years %in% bts_catch_age_yrs),1,1,1] <- rep_dat_rtem$oac_srv1_sample[,1]

# survey length comps inputs
sgl_rg_dusky_data$ObsSrvLenComps <- array(0, dim = c(n_regions, n_years, n_seas, n_lens, n_sexes, sgl_rg_dusky_data$n_srv_fleets))
sgl_rg_dusky_data$UseSrvLenComps <- array(0, dim = c(n_regions, n_years, n_seas, sgl_rg_dusky_data$n_srv_fleets))
sgl_rg_dusky_data$ISS_SrvLenComps <- array(0, dim = c(n_regions, n_years, n_seas, n_sexes, sgl_rg_dusky_data$n_srv_fleets))

sgl_rg_dusky_data$ObsSrvLenComps[1,which(sgl_rg_dusky_data$years %in% bts_catch_size_yrs),1,,1,1] <- oac_srv1_size
sgl_rg_dusky_data$UseSrvLenComps[1,which(sgl_rg_dusky_data$years %in% bts_catch_size_yrs),1,1] <- 0
sgl_rg_dusky_data$ISS_SrvLenComps[1,which(sgl_rg_dusky_data$years %in% bts_catch_size_yrs),1,1,1] <- 0

# Write data
usethis::use_data(sgl_rg_dusky_data, internal = FALSE, overwrite = TRUE)

