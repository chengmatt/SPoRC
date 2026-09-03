# Purpose: certified RE sandeel model, expect jnLL 1749.951 (softplus eps 0.1) with only log_move_diffusion_pars NaN (boundary at zero)
# Date Created: 9/3/26

# nomenclature stuff
# regions are EU1, EU2, UK
# fishery fleets are EU1_s1 EU1_s2 EU2_s1 EU2_s2 UK_s1 UK_s2
# survey fleets are EU1 EU2 UK

library(here)
devtools::load_all(here())

# read in data
dat <- readRDS(here("dev", "spatial_sandeel", "inputs", "north_sea_sandeel_3r.rds"))

# setup model -------------------------------------------------------------

# setup dims
n_pop <- dat$dims$n_pop
n_regions <- dat$dims$n_regions
n_yrs <- dat$dims$n_yrs
n_ages <- dat$dims$n_ages
n_seas <- dat$dims$n_seas
n_sexes <- dat$dims$n_sexes
n_fish <- dat$dims$n_fish_fleets # 6, region x season
n_srv <- dat$dims$n_srv_fleets # 3, EU1 EU2 UK

# setup model dimensions
input_list <- Setup_Mod_Dim(
  years = dat$dims$years,
  ages = dat$dims$ages,
  lens = NA,
  n_pop = n_pop,
  n_regions = n_regions,
  n_seas = n_seas,
  seasdur = dat$dims$seasdur,
  n_sexes = n_sexes,
  n_fish_fleets = n_fish,
  n_srv_fleets = n_srv,
  verbose = TRUE
)

# setup recruitment
input_list <- Setup_Mod_Rec(
  input_list = input_list,
  rec_model = "mean_rec",
  rec_lag = 1,
  spawn_seas = 1,
  t_spawn = 0,
  use_fixed_rec_seas_prop = 1,
  fixed_rec_seas_prop = matrix(c(0, 1), nrow = n_pop, ncol = n_seas), # recruits enter in season 2
  sigmaR_spec = "fix_early_est_late",

  # do initial age strc as equilibrium w/ some devs
  init_age_strc = "matrix",
  equil_init_age_strc = "stoch_all"
)

# starting / fixed values
input_list$par$ln_global_R0 <- log(2e8)
input_list$par$ln_sigmaR[] <- log(4.0128) # same sigR in all regions

# setup biologicals
input_list <- Setup_Mod_Biologicals(
  input_list,
  WAA = dat$WAA,
  WAA_fish = dat$WAA_fish,
  WAA_srv = dat$WAA_srv,
  MatAA = dat$MatAA * 2, # cancels SPoRC's single-sex SSB halving
  fit_lengths = 0,
  M_spec = "fix",
  Fixed_natmort = dat$natmort
)

# setup movement
# EU1 and EU2 do not touch-- both connect through UK
adjacency <- matrix(c(0, 0, 1,
                      0, 0, 1,
                      1, 1, 0), nrow = 3, byrow = TRUE)

# make ctmc dataframe --- must include pop, regions, years, seas, ages, sexes dimensions even if not used
ctmc_data <- expand.grid(
  pop     = 1,
  regions = 1:n_regions,
  years   = 1:n_yrs,
  seas    = 1:n_seas,
  ages    = 1:n_ages,
  sexes   = 1:n_sexes
)

# setup movement
input_list <- Setup_Mod_Movement(
  input_list = input_list,
  use_fixed_movement = 0,
  Fixed_Movement = NA,
  do_recruits_move = 0, # recruits stay put
  move_type = 1,
  move_timing = 2,
  adjacency_mat = as.matrix(adjacency), # adjcacency matrix
  ctmc_move_dat = ctmc_data,
  diffusion_formula  = ~1,  # intercept
  preference_formula = ~0 + factor(regions),
  ctmc_diffusion_bounds = 1,
  area_r = c(38670, 89611, 136781)  # area size
)

# setup catch
input_list <- Setup_Mod_Catch_and_F(
  input_list,

  # not using aggregated catch streams
  ObsCatch = array(0, dim = c(n_regions, n_yrs, n_seas, n_fish)),
  UseCatch = array(0, dim = c(n_regions, n_yrs, n_seas, n_fish)),

  # using catch-at-age data streams
  ObsCatchAA = dat$ObsCatchAA,
  UseCatchAA = dat$UseCatchAA,
  CatchAA_Type = "spltRspltS", # region split
  sigmaCAA_key = dat$sigmaCAA_key, # sharing sigmas for age 1-3, 4-5
  sigmaCAA_spec = "est", # estimate sigma for CAA
  ln_sigmaCAA = array(log(0.5), dim = c(n_ages, n_sexes, n_fish)), # starting value for sigma CAA
  catch_units = array("abd", dim = c(n_fish)),
  Use_F_pen = 0, # free Fs
  sigmaC_spec = "fix",
  sigmaF_spec = "fix"
)

# setup fishery index (not used since using catch-at-age ... )
input_list <- Setup_Mod_FishIdx_and_Comps(
  input_list,

  # no fishery index
  ObsFishIdx = array(NA, dim = c(n_regions, n_yrs, n_seas, n_fish)),
  ObsFishIdx_SE = array(NA, dim = c(n_regions, n_yrs, n_seas, n_fish)),
  UseFishIdx = array(0, dim = c(n_regions, n_yrs, n_seas, n_fish)),

  # no fishery comps
  ObsFishAgeComps = array(0, dim = c(n_regions, n_yrs, n_seas, n_ages, n_sexes, n_fish)),
  UseFishAgeComps = array(0, dim = c(n_regions, n_yrs, n_seas, n_fish)),
  ISS_FishAgeComps = array(0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_fish)),
  ObsFishLenComps = array(0, dim = c(n_regions, n_yrs, n_seas, 1, n_sexes, n_fish)),
  UseFishLenComps = array(0, dim = c(n_regions, n_yrs, n_seas, n_fish)),
  ISS_FishLenComps = array(0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_fish)),
  fish_idx_type = c("none", "none", "none", "none", "none", "none"),
  FishAgeComps_LikeType = c("none", "none", "none", "none", "none", "none"),
  FishLenComps_LikeType = c("none", "none", "none", "none", "none", "none"),
  FishAgeComps_Type = c("none_Year_1-terminal_Fleet_1",
                        "none_Year_1-terminal_Fleet_2",
                        "none_Year_1-terminal_Fleet_3",
                        "none_Year_1-terminal_Fleet_4",
                        "none_Year_1-terminal_Fleet_5",
                        "none_Year_1-terminal_Fleet_6"),
  FishLenComps_Type = c("none_Year_1-terminal_Fleet_1",
                        "none_Year_1-terminal_Fleet_2",
                        "none_Year_1-terminal_Fleet_3",
                        "none_Year_1-terminal_Fleet_4",
                        "none_Year_1-terminal_Fleet_5",
                        "none_Year_1-terminal_Fleet_6")
)

# setup survey indices
input_list <- Setup_Mod_SrvIdx_and_Comps(
  input_list,

  # no aggregate index used
  ObsSrvIdx = dat$ObsSrvIdx,
  ObsSrvIdx_SE = array(0, dim = c(n_regions, n_yrs, n_seas, n_srv)),
  UseSrvIdx = dat$UseSrvIdx,

  # using survey index at age
  ObsSrvIdxAA = dat$ObsSrvIdxAA,
  UseSrvIdxAA = dat$UseSrvIdxAA,

  # survey index time - spatially resolved
  SrvIdxAA_Type = c("spltRspltS", "spltRspltS", "spltRspltS"),

  sigmaSrvIdxAA_key = dat$sigmaSrvIdxAA_key, # sigmas share age 0-1 and then separate for age 2
  sigmaSrvIdxAA_spec = "est",
  ln_sigmaSrvIdxAA = dat$ln_sigmaSrvIdxAA,
  srv_idx_type = c("abd", "abd", "abd"),

  # no marginal age comps
  ObsSrvAgeComps = array(0, dim = c(n_regions, n_yrs, n_seas, n_ages, n_sexes, n_srv)),
  UseSrvAgeComps = array(0, dim = c(n_regions, n_yrs, n_seas, n_srv)),
  ISS_SrvAgeComps = array(0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_srv)),
  ObsSrvLenComps = array(0, dim = c(n_regions, n_yrs, n_seas, 1, n_sexes, n_srv)),
  UseSrvLenComps = array(0, dim = c(n_regions, n_yrs, n_seas, n_srv)),
  ISS_SrvLenComps = array(0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_srv)),

  # marginal age comps - not used
  SrvAgeComps_LikeType = c("none", "none", "none"),
  SrvLenComps_LikeType = c("none", "none", "none"),
  SrvAgeComps_Type = c("none_Year_1-terminal_Fleet_1",
                       "none_Year_1-terminal_Fleet_2",
                       "none_Year_1-terminal_Fleet_3"),
  SrvLenComps_Type = c("none_Year_1-terminal_Fleet_1",
                       "none_Year_1-terminal_Fleet_2",
                       "none_Year_1-terminal_Fleet_3")
)

# setup fishery selex (non-parametric)
input_list <- Setup_Mod_Fishsel_and_Q(
  input_list,
  cont_tv_fish_sel = c("none_Fleet_1",
                       "none_Fleet_2",
                       "none_Fleet_3",
                       "none_Fleet_4",
                       "none_Fleet_5",
                       "none_Fleet_6"),
  fish_sel_model = c("nonparfree_Fleet_1",
                     "nonparfree_Fleet_2",
                     "nonparfree_Fleet_3",
                     "nonparfree_Fleet_4",
                     "nonparfree_Fleet_5",
                     "nonparfree_Fleet_6"),
  fish_q_blocks = c("none_Fleet_1",
                    "none_Fleet_2",
                    "none_Fleet_3",
                    "none_Fleet_4",
                    "none_Fleet_5",
                    "none_Fleet_6"),

  # [[fleet]][[block]][[bin group]] .. ages 0-4 as bins 1-5, ages 3 and 4 grouped
  fish_sel_nonpar_est_bins = list(
    # fleet 1, EU1_s1
    list(
      list(1, 2, 3, c(4, 5))
    ),
    # fleet 2, EU1_s2
    list(
      list(1, 2, 3, c(4, 5))
    ),
    # fleet 3, EU2_s1
    list(
      list(1, 2, 3, c(4, 5))
    ),
    # fleet 4, EU2_s2
    list(
      list(1, 2, 3, c(4, 5))
    ),
    # fleet 5, UK_s1
    list(
      list(1, 2, 3, c(4, 5))
    ),
    # fleet 6, UK_s2
    list(
      list(1, 2, 3, c(4, 5))
    )

  ),
  fish_fixed_sel_pars_spec = c("est_all", "est_shared_f_1", "est_shared_f_1",
                               "est_shared_f_1", "est_shared_f_1", "est_shared_f_1"), # share selex with first fleet - i.e., spatially and seasonally invariant
  fish_q_spec = c("fix", "fix", "fix", "fix", "fix", "fix")
)

# setup survey selex (non-parametric, mapped to the ages each survey sees)
# figure out srv timing
t_srv <- array(0, dim = c(n_regions, n_seas, n_srv))
for(k in 1:n_srv) t_srv[, dat$dims$srv_seas[k], k] <- dat$dims$srv_t[k]

# setup survey selex stuff
input_list <- Setup_Mod_Srvsel_and_Q(
  input_list,
  cont_tv_srv_sel = c("none_Fleet_1",
                      "none_Fleet_2",
                      "none_Fleet_3"),
  srv_sel_blocks = c("none_Fleet_1",
                     "none_Fleet_2",
                     "none_Fleet_3"),
  srv_sel_model = c("nonparfree_Fleet_1",
                    "nonparfree_Fleet_2",
                    "nonparfree_Fleet_3"),
  srv_q_blocks = c("none_Fleet_1",
                   "none_Fleet_2",
                   "none_Fleet_3"),

  # [[fleet]][[block]][[bin group]], one bin per age the survey observes
  srv_sel_nonpar_est_bins = list(
    # fleet 1, EU1, ages 0-2
    list(
      list(1, 2, 3)
    ),
    # fleet 2, EU2, ages 0-2
    list(
      list(1, 2, 3)
    ),
    # fleet 3, UK, ages 0-2
    list(
      list(1, 2, 3)
    )

  ),
  srv_fixed_sel_pars_spec = c("est_all", "est_shared_f_1", "est_shared_f_1"), # share selex here to get consistent scaling
  # nonparfree absorbs the catchability at age, so q fixed at one
  srv_q_spec = c("fix", "fix", "fix"),
  t_srv = t_srv
)

# no tagging
input_list <- Setup_Mod_Tagging(
  input_list = input_list,
  use_conv_fish_tagging = 0
)

# setup weighting
input_list <- Setup_Mod_Weighting(
  input_list,
  Wt_Catch = 1,
  Wt_FishIdx = 0,
  Wt_SrvIdx = 1,
  Wt_Rec = 1,
  Wt_F = 1,
  Wt_Tagging = 0,
  Wt_FishAgeComps = array(0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_fish)),
  Wt_FishLenComps = array(0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_fish)),
  Wt_SrvAgeComps = array(0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_srv)),
  Wt_SrvLenComps = array(0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_srv))
)

# fishery selectivity ------------------------------------------------------------
# setup selx mapping
sel_start <- input_list$par$fish_fixed_sel_pars
sel_map <- array(0, dim = dim(sel_start))
sel_start[, 1, , 1, ] <- -10  # age 0, never fished so map off at small value
for(f in 1:n_fish) for(r in 1:n_regions) sel_map[r, , 1, 1, f] <- c(NA, 1, 2, NA, NA) # what pars to fix / map off ... fix age-0, then map off 4-5 shared at a fixed term so identifiable since free ln_f_mean
input_list$par$fish_fixed_sel_pars <- sel_start
input_list$map$fish_fixed_sel_pars <- factor(sel_map)

# effort driven fishing mortality -------------------------------------------------
ln_F_devs <- input_list$par$ln_F_devs
dev_map <- array(NA, dim = dim(ln_F_devs))
mean_effort <- apply(dat$effort, c(1, 3), mean, na.rm = TRUE)
free <- 0

for(r in 1:n_regions) {
  for(s in 1:n_seas) {
    f <- dat$dims$fish_fleet[r, s]
    for(y in 1:n_yrs) {

      st <- dat$effort_status[r, y, s] # figure out if have effort

      # effort observed, hold the dev at it
      if(st == 1) ln_F_devs[r, y, s, f] <- log(dat$effort[r, y, s])
      # no effort, estimate it, but only where catch is actually fit
      if(st == 2) {
        if(sum(dat$UseCatchAA[r, y, s, , , f]) > 0) {
          free <- free + 1
          ln_F_devs[r, y, s, f] <- log(mean_effort[r, s])
          dev_map[r, y, s, f] <- free
        } else ln_F_devs[r, y, s, f] <- log(1e-12) # zero out Fdev if no catch
      }
      # no landings, no fishing
      if(st == 3) ln_F_devs[r, y, s, f] <- log(1e-12) # zero out fdev
    } # end y loop

  } # end s loop
} # end r loop

# map stuff out
input_list$par$ln_F_devs <- ln_F_devs
input_list$map$ln_F_devs <- factor(dev_map)

# save stuff
saveRDS(input_list, here("dev", "spatial_sandeel", "inputs", "spatial_input_list_re.rds"))



# fit ---------------------------------------------------------------------

# laplace over recruitment deviations
mod_est <- fit_model(
  input_list$data,
  input_list$par,
  input_list$map,
  c('ln_RecDevs'),
  do_optim = TRUE,
  silent = FALSE,
  newton_loops = 3
)

# polish with restarts, keep the best point with finite gradients
best <- mod_est$optim
for(i in 1:3) {
  cand <- tryCatch(stats::nlminb(best$par, mod_est$fn, mod_est$gr,
            control = list(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-15)), error = function(e) NULL)
  if(!is.null(cand)) {
    g <- tryCatch(mod_est$gr(cand$par), error = function(e) NA)
    if(all(is.finite(g)) && cand$objective <= best$objective) best <- cand
  }
} # end i loop
mod_est$optim <- best

# plain sdreport under random effects, he() has no laplace tape
sd_rep <- sdreport(mod_est)

saveRDS(list(mod_est = mod_est$optim, sd_rep = sd_rep, rep = mod_est$report(mod_est$env$last.par.best)),
        here("dev", "spatial_sandeel", "inputs", "re_certified_fit.rds"))
