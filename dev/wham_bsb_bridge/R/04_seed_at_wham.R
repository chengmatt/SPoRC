# Purpose: hold every SPoRC parameter at the wham fit_0 values and report the population dynamics
# Date Created: 9/7/26

library(here)
library(SPoRC)

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
