# Purpose: reshape the wham black sea bass fit_0 inputs and report into SPoRC arrays
# Date Created: 9/7/26

library(here)

out_dir <- here("dev", "wham_bsb_bridge", "output")
wham <- readRDS(file.path(out_dir, "01_wham_fit0.rds"))
dat <- wham$input$data
rep <- wham$rep
par <- wham$parList

# dimensions. wham stocks become SPoRC populations, each homing to its own region
n_pop <- 2
n_regions <- 2
n_yrs <- dat$n_years_model
n_seas <- dat$n_seasons
n_ages <- dat$n_ages
n_sexes <- 1
n_fish <- 4 # north commercial, north recreational, south commercial, south recreational
n_srv <- 4  # north rec cpa, north vast, south rec cpa, south vast

years <- 1989:2021
ages <- 1:n_ages
seasdur <- dat$fracyr_seasons
natal_region <- c(1, 2)

# each wham fleet and index keeps its own SPoRC fleet, so the composition likelihoods,
# which are set per fleet, can differ between the north and the south the way wham has them
fleet_region <- dat$fleet_regions          # 1 1 2 2
fleet_gear <- 1:4
index_region <- dat$index_regions          # 1 1 2 2
index_gear <- 1:4

# selectivity blocks: fleets take 1-4, indices take 5-8
fleet_block <- 1:4
index_block <- 5:8

# spawning happens partway into season 6, t_spawn is a fraction of that season
spawn_seas <- dat$spawn_seasons[1]
t_spawn <- dat$fracyr_SSB[1, 1] / seasdur[spawn_seas]

# survey timing, also a fraction of the season the index sits in
t_srv <- array(0, dim = c(n_regions, n_seas, n_srv))
for(i in 1:dat$n_indices) {
  t_srv[index_region[i], dat$index_seasons[i], index_gear[i]] <- dat$fracyr_indices[1, i] / seasdur[dat$index_seasons[i]]
}

# weight and maturity at age, spawning
WAA <- array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes))
MatAA <- array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes))
for(p in 1:n_pop) {
  for(r in 1:n_regions) {
    for(seas in 1:n_seas) {
      WAA[p, r, , seas, , 1] <- dat$waa[dat$waa_pointer_ssb[p], , ]
      MatAA[p, r, , seas, , 1] <- dat$mature[p, , ]
    } # end seas loop
  } # end r loop
} # end p loop

# fishery and survey weight at age, each wham fleet lands on its own region and gear
WAA_fish <- array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_fish))
for(f in 1:dat$n_fleets) {
  for(p in 1:n_pop) for(seas in 1:n_seas) {
    WAA_fish[p, fleet_region[f], , seas, , 1, fleet_gear[f]] <- dat$waa[dat$waa_pointer_fleets[f], , ]
  }
} # end f loop

WAA_srv <- array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_srv))
for(i in 1:dat$n_indices) {
  for(p in 1:n_pop) for(seas in 1:n_seas) {
    WAA_srv[p, index_region[i], , seas, , 1, index_gear[i]] <- dat$waa[dat$waa_pointer_indices[i], , ]
  }
} # end i loop

# natural mortality, an annual rate that SPoRC scales by season duration
natmort <- array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes))
for(p in 1:n_pop) for(r in 1:n_regions) for(seas in 1:n_seas) natmort[p, r, , seas, , 1] <- rep$MAA[p, r, , ]

# movement, taken straight from the wham transition matrices so must_move comes along
Fixed_Movement <- array(0, dim = c(n_pop, n_regions, n_regions, n_yrs, n_seas, n_ages, n_sexes))
for(p in 1:n_pop) for(y in 1:n_yrs) for(seas in 1:n_seas) for(a in 1:n_ages) {
  Fixed_Movement[p, , , y, seas, a, 1] <- rep$mu[p, a, seas, y, , ]
}

# selectivity at age by block, already holding the wham random effects
selAA <- array(0, dim = c(dat$n_selblocks, n_yrs, n_ages))
for(b in 1:dat$n_selblocks) selAA[b, , ] <- rep$selAA[[b]]

# fishing mortality. wham FAA is F times selectivity at an annual rate
F_ann <- matrix(0, n_yrs, dat$n_fleets)
for(f in 1:dat$n_fleets) {
  for(y in 1:n_yrs) {
    ref_age <- which.max(selAA[fleet_block[f], y, ])
    F_ann[y, f] <- rep$FAA[f, y, ref_age] / selAA[fleet_block[f], y, ref_age]
  }
} # end f loop

# aggregate catch. wham fits one annual total per fleet, so it sits in season 1 and is
# fit against the season total
obs_seas <- 1
ObsCatch <- array(0, dim = c(n_regions, n_yrs, n_seas, n_fish))
UseCatch <- array(0, dim = c(n_regions, n_yrs, n_seas, n_fish))
ObsCatch_SE <- array(0, dim = c(n_regions, n_yrs, n_seas, n_fish))
for(f in 1:dat$n_fleets) {
  ObsCatch[fleet_region[f], , obs_seas, fleet_gear[f]] <- dat$agg_catch[, f]
  UseCatch[fleet_region[f], , obs_seas, fleet_gear[f]] <- dat$use_agg_catch[, f]
  ObsCatch_SE[fleet_region[f], , obs_seas, fleet_gear[f]] <- dat$agg_catch_sigma[, f] * exp(par$log_catch_sig_scale[f])
}

# survey indices sit in one season each, so they map across without reshaping
ObsSrvIdx <- array(0, dim = c(n_regions, n_yrs, n_seas, n_srv))
ObsSrvIdx_SE <- array(0, dim = c(n_regions, n_yrs, n_seas, n_srv))
UseSrvIdx <- array(0, dim = c(n_regions, n_yrs, n_seas, n_srv))
for(i in 1:dat$n_indices) {
  seas <- dat$index_seasons[i]
  ObsSrvIdx[index_region[i], , seas, index_gear[i]] <- dat$agg_indices[, i]
  ObsSrvIdx_SE[index_region[i], , seas, index_gear[i]] <- dat$agg_index_sigma[, i] * exp(par$log_index_sig_scale[i])
  UseSrvIdx[index_region[i], , seas, index_gear[i]] <- dat$use_indices[, i]
}

# survey age compositions, in the same season as their index
ObsSrvAgeComps <- array(0, dim = c(n_regions, n_yrs, n_seas, n_ages, n_sexes, n_srv))
UseSrvAgeComps <- array(0, dim = c(n_regions, n_yrs, n_seas, n_srv))
ISS_SrvAgeComps <- array(0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_srv))
for(i in 1:dat$n_indices) {
  seas <- dat$index_seasons[i]
  ObsSrvAgeComps[index_region[i], , seas, , 1, index_gear[i]] <- dat$index_paa[i, , ]
  UseSrvAgeComps[index_region[i], , seas, index_gear[i]] <- dat$use_index_paa[, i]
  ISS_SrvAgeComps[index_region[i], , seas, 1, index_gear[i]] <- dat$index_Neff[, i]
}

# fishery age compositions, also annual, so they sit alongside the catch in season 1
ObsFishAgeComps <- array(0, dim = c(n_regions, n_yrs, n_seas, n_ages, n_sexes, n_fish))
UseFishAgeComps <- array(0, dim = c(n_regions, n_yrs, n_seas, n_fish))
ISS_FishAgeComps <- array(0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_fish))
for(f in 1:dat$n_fleets) {
  ObsFishAgeComps[fleet_region[f], , obs_seas, , 1, fleet_gear[f]] <- dat$catch_paa[f, , ]
  UseFishAgeComps[fleet_region[f], , obs_seas, fleet_gear[f]] <- dat$use_catch_paa[, f]
  ISS_FishAgeComps[fleet_region[f], , obs_seas, 1, fleet_gear[f]] <- dat$catch_Neff[, f]
}

# wham quantities the comparison reads back
wham_targets <- list(
  NAA = rep$NAA,                   # stock, region, year, age at jan 1
  pred_NAA = rep$pred_NAA,         # one step ahead prediction of the same
  N1 = rep$N1,
  SSB = rep$SSB,
  pred_catch = rep$pred_catch,
  pred_indices = rep$pred_indices,
  pred_CAA = rep$pred_CAA,
  pred_IAA = rep$pred_IAA,
  FAA = rep$FAA,
  q = rep$q,
  MAA = rep$MAA,
  NAA_spawn = rep$NAA_spawn,
  all_NAA = rep$all_NAA,
  nll_agg_catch = rep$nll_agg_catch,
  nll_agg_indices = rep$nll_agg_indices,
  nll_catch_acomp = rep$nll_catch_acomp,
  nll_index_acomp = rep$nll_index_acomp,
  nll_NAA = rep$nll_NAA,
  catch_paa_pars = par$catch_paa_pars,
  index_paa_pars = par$index_paa_pars,
  catch_paa_model = dat$age_comp_model_fleets,
  index_paa_model = dat$age_comp_model_indices,
  mean_rec_pars = par$mean_rec_pars,
  NAA_sigma = par$log_NAA_sigma,
  NAA_rho = par$trans_NAA_rho,
  N1_pars = par$log_N1,
  sel_repars = par$sel_repars
)

saveRDS(list(dims = list(n_pop = n_pop, n_regions = n_regions, n_yrs = n_yrs, n_seas = n_seas,
                         n_ages = n_ages, n_sexes = n_sexes, n_fish = n_fish, n_srv = n_srv,
                         years = years, ages = ages, seasdur = seasdur, natal_region = natal_region,
                         spawn_seas = spawn_seas, t_spawn = t_spawn,
                         fleet_region = fleet_region, fleet_gear = fleet_gear, fleet_block = fleet_block,
                         index_region = index_region, index_gear = index_gear, index_block = index_block,
                         index_seas = dat$index_seasons),
             WAA = WAA,
             MatAA = MatAA,
             WAA_fish = WAA_fish,
             WAA_srv = WAA_srv,
             natmort = natmort,
             Fixed_Movement = Fixed_Movement,
             selAA = selAA,
             F_ann = F_ann,
             t_srv = t_srv,
             obs_seas = obs_seas,
             ObsCatch = ObsCatch,
             ObsCatch_SE = ObsCatch_SE,
             UseCatch = UseCatch,
             ObsSrvIdx = ObsSrvIdx,
             ObsSrvIdx_SE = ObsSrvIdx_SE,
             UseSrvIdx = UseSrvIdx,
             ObsSrvAgeComps = ObsSrvAgeComps,
             UseSrvAgeComps = UseSrvAgeComps,
             ISS_SrvAgeComps = ISS_SrvAgeComps,
             ObsFishAgeComps = ObsFishAgeComps,
             UseFishAgeComps = UseFishAgeComps,
             ISS_FishAgeComps = ISS_FishAgeComps,
             wham = wham_targets),
        file.path(out_dir, "02_sporc_data.rds"))
