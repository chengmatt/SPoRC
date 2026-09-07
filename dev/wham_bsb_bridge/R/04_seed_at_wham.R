# Purpose: hold every SPoRC parameter at the wham fit_0 values and report the population dynamics
# Date Created: 9/7/26

library(here)
suppressMessages(pkgload::load_all(here::here(), quiet = TRUE))

out_dir <- here("dev", "wham_bsb_bridge", "output")
dat <- readRDS(file.path(out_dir, "02_sporc_data.rds"))
input_list <- readRDS(file.path(out_dir, "03_input_list.rds"))
wham <- dat$wham

n_pop <- dat$dims$n_pop
n_regions <- dat$dims$n_regions
n_yrs <- dat$dims$n_yrs
n_seas <- dat$dims$n_seas
n_ages <- dat$dims$n_ages
n_fish <- dat$dims$n_fish
n_srv <- dat$dims$n_srv
seasdur <- dat$dims$seasdur
natal <- dat$dims$natal_region

floor_n <- 1e-10 # stands in for the cells wham holds at exactly zero

# fishing mortality. wham F is an annual rate, SPoRC Fmort is already integrated over the season
input_list$par$ln_F_mean[] <- 0
for(f in 1:length(dat$dims$fleet_region)) {
  r <- dat$dims$fleet_region[f]
  gear <- dat$dims$fleet_gear[f]
  for(seas in 1:n_seas) input_list$par$ln_F_devs[r, , seas, gear] <- log(dat$F_ann[, f] * seasdur[seas])
} # end f loop

# aggregate catch observation error, taken from wham's own reported standard deviations
for(f in 1:length(dat$dims$fleet_region)) {
  r <- dat$dims$fleet_region[f]
  gear <- dat$dims$fleet_gear[f]
  for(seas in 1:n_seas) input_list$par$ln_sigmaC[r, , seas, gear] <- log(dat$ObsCatch_SE[r, , dat$obs_seas, gear])
} # end f loop

# survey catchability, time invariant in this fit
for(i in 1:length(dat$dims$index_region)) {
  input_list$par$ln_srv_q[dat$dims$index_region[i], 1, dat$dims$index_gear[i]] <- log(wham$q[1, i])
}

# year one numbers at age, ages 2 and above, read straight off the wham equilibrium
for(p in 1:n_pop) for(r in 1:n_regions) {
  input_list$par$ln_InitDevs[p, r, , 1] <- log(pmax(wham$N1[p, r, 2:n_ages], floor_n))
}

# the state-space numbers at age for years two onward, ages two and above
for(p in 1:n_pop) for(r in 1:n_regions) for(y in 2:n_yrs) {
  input_list$par$ln_NAA[p, r, y, 1, 2:n_ages, 1] <- log(pmax(wham$NAA[p, r, y, 2:n_ages], floor_n))
}

# dirichlet multinomial dispersion. wham uses the saturating form, alpha = p * exp(theta), and
# SPoRC the linear form, alpha = p * exp(theta) * ISS, so the two agree once ISS is divided out
dm_fleets <- which(dat$wham$catch_paa_model == 2)
for(f in dm_fleets) {
  r <- dat$dims$fleet_region[f]
  input_list$par$ln_FishAge_theta[r, 1, dat$dims$fleet_gear[f]] <- dat$wham$catch_paa_pars[f, 1] - log(dat$ISS_FishAgeComps[r, 1, dat$obs_seas, 1, dat$dims$fleet_gear[f]])
} # end f loop

dm_srv <- which(dat$wham$index_paa_model == 2)
for(i in dm_srv) {
  r <- dat$dims$index_region[i]
  seas <- dat$dims$index_seas[i]
  input_list$par$ln_SrvAge_theta[r, 1, dat$dims$index_gear[i]] <- dat$wham$index_paa_pars[i, 1] - log(dat$ISS_SrvAgeComps[r, 1, seas, 1, dat$dims$index_gear[i]])
} # end i loop

# logistic normal dispersion and correlation, taken straight off wham. the sample size
# scaling lives inside the likelihood, so the parameter is the one wham estimates
ln_fleets <- which(dat$wham$catch_paa_model %in% c(5, 6))
for(f in ln_fleets) {
  r <- dat$dims$fleet_region[f]
  gear <- dat$dims$fleet_gear[f]
  input_list$par$ln_FishAge_theta[r, 1, gear] <- dat$wham$catch_paa_pars[f, 1]
  input_list$par$FishAge_corr_pars[r, 1, gear, 1] <- dat$wham$catch_paa_pars[f, 2]
} # end f loop

ln_srv <- which(dat$wham$index_paa_model %in% c(5, 6))
for(i in ln_srv) {
  r <- dat$dims$index_region[i]
  gear <- dat$dims$index_gear[i]
  input_list$par$ln_SrvAge_theta[r, 1, gear] <- dat$wham$index_paa_pars[i, 1]
  input_list$par$SrvAge_corr_pars[r, 1, gear, 1] <- dat$wham$index_paa_pars[i, 2]
} # end i loop

# every parameter is held, nothing is estimated here
map_all <- lapply(input_list$par, function(x) factor(rep(NA, length(x))))
input_list$map[names(map_all)] <- map_all

# first pass with flat recruitment deviations, so the mean recruitment scale can be read off
seed0 <- fit_model(data = input_list$data,
                   parameters = input_list$par,
                   mapping = input_list$map,
                   random = NULL,
                   do_optim = FALSE,
                   silent = TRUE)

rep0 <- seed0$rep

# recruitment deviations that land age one on the wham numbers, natal region only
for(p in 1:n_pop) {
  wham_rec <- wham$NAA[p, natal[p], , 1]
  input_list$par$ln_RecDevs[p, natal[p], ] <- log(wham_rec / rep0$Rec[p, natal[p], ])
}

seed <- fit_model(data = input_list$data,
                  parameters = input_list$par,
                  mapping = input_list$map,
                  random = NULL,
                  do_optim = FALSE,
                  silent = TRUE)

saveRDS(list(input_list = input_list, seed = seed, rep = seed$rep),
        file.path(out_dir, "04_seeded.rds"))
