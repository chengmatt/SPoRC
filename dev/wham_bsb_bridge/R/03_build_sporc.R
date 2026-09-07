# Purpose: build the SPoRC black sea bass model and hold every process at the wham fit_0 values
# Date Created: 9/7/26

library(here)
library(SPoRC)

out_dir <- here("dev", "wham_bsb_bridge", "output")
dat <- readRDS(file.path(out_dir, "02_sporc_data.rds"))

n_pop <- dat$dims$n_pop
n_regions <- dat$dims$n_regions
n_yrs <- dat$dims$n_yrs
n_seas <- dat$dims$n_seas
n_ages <- dat$dims$n_ages
n_sexes <- dat$dims$n_sexes
n_fish <- dat$dims$n_fish
n_srv <- dat$dims$n_srv

# two stocks homing to their own region, eleven seasons, one sex
input_list <- Setup_Mod_Dim(
  years = dat$dims$years,
  ages = dat$dims$ages,
  lens = NA,
  n_pop = n_pop,
  n_regions = n_regions,
  natal_region = dat$dims$natal_region,
  n_seas = n_seas,
  seasdur = dat$dims$seasdur,
  n_sexes = n_sexes,
  n_fish_fleets = n_fish,
  n_srv_fleets = n_srv,
  verbose = FALSE
)

# recruitment about a mean, entering in season 1 of the natal region
input_list <- Setup_Mod_Rec(
  input_list = input_list,
  rec_model = "mean_rec",
  rec_dd = "local",
  rec_lag = 1,
  spawn_seas = dat$dims$spawn_seas,
  t_spawn = dat$dims$t_spawn,
  use_fixed_rec_seas_prop = 1,
  fixed_rec_seas_prop = matrix(rep(c(1, rep(0, n_seas - 1)), each = n_pop), nrow = n_pop, ncol = n_seas),
  rec_region_prop_spec = "no_dispersal",
  sigmaR_spec = "fix",
  do_rec_bias_ramp = 0,
  init_age_strc = "free",
  equil_init_age_strc = "stoch_all",
  InitDevs_spec = "est_shared_r",
  RecDevs_spec = "est_shared_r"
)

# biologicals, with the numbers at age state active over every age and year so wham's
# realized state can be imposed cell by cell
input_list <- Setup_Mod_Biologicals(
  input_list,
  WAA = dat$WAA,
  WAA_fish = dat$WAA_fish,
  WAA_srv = dat$WAA_srv,
  MatAA = dat$MatAA * 2, # cancels the single sex halving so ssb matches wham
  fit_lengths = 0,
  M_spec = "fix",
  Fixed_natmort = dat$natmort,
  NAA_re = "iid",
  NAA_sigma_spec = "fix",
  NAA_re_ages = 2:n_ages,
  NAA_re_years = dat$dims$years[-1]
)

# movement taken from the wham transition matrices, applied after mortality each season
input_list <- Setup_Mod_Movement(
  input_list = input_list,
  use_fixed_movement = 1,
  Fixed_Movement = dat$Fixed_Movement,
  do_recruits_move = 1,
  move_type = 0,
  move_timing = 1
)

# catch. the seasonal split of the observations is a placeholder, F is held at wham values
input_list <- Setup_Mod_Catch_and_F(
  input_list,
  ObsCatch = dat$ObsCatch,
  UseCatch = dat$UseCatch,
  catch_units = c("biom", "biom"),
  Use_F_pen = 0,
  sigmaC_spec = "fix",
  sigmaF_spec = "fix"
)

# no fishery index or fishery comps in this model
input_list <- Setup_Mod_FishIdx_and_Comps(
  input_list,
  ObsFishIdx = array(NA, dim = c(n_regions, n_yrs, n_seas, n_fish)),
  ObsFishIdx_SE = array(NA, dim = c(n_regions, n_yrs, n_seas, n_fish)),
  UseFishIdx = array(0, dim = c(n_regions, n_yrs, n_seas, n_fish)),
  ObsFishAgeComps = dat$ObsFishAgeComps,
  UseFishAgeComps = dat$UseFishAgeComps,
  ISS_FishAgeComps = dat$ISS_FishAgeComps,
  ObsFishLenComps = array(0, dim = c(n_regions, n_yrs, n_seas, 1, n_sexes, n_fish)),
  UseFishLenComps = array(0, dim = c(n_regions, n_yrs, n_seas, n_fish)),
  ISS_FishLenComps = array(0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_fish)),
  fish_idx_type = c("none", "none"),
  FishAgeComps_LikeType = c("none", "none"),
  FishLenComps_LikeType = c("none", "none"),
  FishAgeComps_Type = c("none_Year_1-terminal_Fleet_1",
                        "none_Year_1-terminal_Fleet_2"),
  FishLenComps_Type = c("none_Year_1-terminal_Fleet_1",
                        "none_Year_1-terminal_Fleet_2")
)

# survey indices in numbers, each sitting in its own season
input_list <- Setup_Mod_SrvIdx_and_Comps(
  input_list,
  ObsSrvIdx = dat$ObsSrvIdx,
  ObsSrvIdx_SE = dat$ObsSrvIdx_SE,
  UseSrvIdx = dat$UseSrvIdx,
  srv_idx_type = c("abd", "abd"),
  ObsSrvAgeComps = dat$ObsSrvAgeComps,
  UseSrvAgeComps = dat$UseSrvAgeComps,
  ISS_SrvAgeComps = dat$ISS_SrvAgeComps,
  ObsSrvLenComps = array(0, dim = c(n_regions, n_yrs, n_seas, 1, n_sexes, n_srv)),
  UseSrvLenComps = array(0, dim = c(n_regions, n_yrs, n_seas, n_srv)),
  ISS_SrvLenComps = array(0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_srv)),
  SrvAgeComps_LikeType = c("Multinomial", "Multinomial"),
  SrvLenComps_LikeType = c("none", "none"),
  SrvAgeComps_Type = c("spltRjntS_Year_1-terminal_Fleet_1",
                       "spltRjntS_Year_1-terminal_Fleet_2"),
  SrvLenComps_Type = c("none_Year_1-terminal_Fleet_1",
                       "none_Year_1-terminal_Fleet_2")
)

# fishery selectivity at age read straight off the wham blocks
fish_sel_input <- array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_fish))
for(f in 1:length(dat$dims$fleet_region)) {
  r <- dat$dims$fleet_region[f]
  gear <- dat$dims$fleet_gear[f]
  for(p in 1:n_pop) for(seas in 1:n_seas) fish_sel_input[p, r, , seas, , 1, gear] <- dat$selAA[dat$dims$fleet_block[f], , ]
} # end f loop

input_list <- Setup_Mod_Fishsel_and_Q(
  input_list,
  cont_tv_fish_sel = c("none_Fleet_1", "none_Fleet_2"),
  fish_sel_model = c("logist1_Fleet_1", "logist1_Fleet_2"),
  fish_sel_blocks = c("none_Fleet_1", "none_Fleet_2"),
  fish_q_blocks = c("none_Fleet_1", "none_Fleet_2"),
  fish_q_spec = c("fix", "fix"),
  fish_fixed_sel_pars_spec = c("fix_fish_sel_input", "fix_fish_sel_input"),
  use_fixed_fish_sel = c(1, 1),
  fish_sel_input = fish_sel_input
)

# survey selectivity at age, same route
srv_sel_input <- array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_srv))
for(i in 1:length(dat$dims$index_region)) {
  r <- dat$dims$index_region[i]
  gear <- dat$dims$index_gear[i]
  for(p in 1:n_pop) for(seas in 1:n_seas) srv_sel_input[p, r, , seas, , 1, gear] <- dat$selAA[dat$dims$index_block[i], , ]
} # end i loop

input_list <- Setup_Mod_Srvsel_and_Q(
  input_list,
  cont_tv_srv_sel = c("none_Fleet_1", "none_Fleet_2"),
  srv_sel_model = c("logist1_Fleet_1", "logist1_Fleet_2"),
  srv_sel_blocks = c("none_Fleet_1", "none_Fleet_2"),
  srv_q_blocks = c("none_Fleet_1", "none_Fleet_2"),
  srv_q_spec = c("fix", "fix"),
  srv_fixed_sel_pars_spec = c("fix_srv_sel_input", "fix_srv_sel_input"),
  use_fixed_srv_sel = c(1, 1),
  srv_sel_input = srv_sel_input,
  t_srv = dat$t_srv
)

input_list <- Setup_Mod_Tagging(
  input_list = input_list,
  use_conv_fish_tagging = 0
)

input_list <- Setup_Mod_Weighting(
  input_list,
  Wt_Catch = 1,
  Wt_FishIdx = 0,
  Wt_SrvIdx = 1,
  Wt_Rec = 1,
  Wt_F = 0,
  Wt_Tagging = 0,
  Wt_FishAgeComps = array(0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_fish)),
  Wt_FishLenComps = array(0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_fish)),
  Wt_SrvAgeComps = array(1, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_srv)),
  Wt_SrvLenComps = array(0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_srv))
)

saveRDS(input_list, file.path(out_dir, "03_input_list.rds"))
