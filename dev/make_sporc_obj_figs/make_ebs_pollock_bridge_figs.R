# Purpose: Bridge the 2024 EBS walleye pollock assessment to SPoRC and render
#          the case study figures. Everything is specified the way the pollock
#          assessment specifies it: Ricker recruitment with the stock recruit
#          penalty over its year window, steepness estimated under its beta
#          prior, free initial numbers at age, selectivity through SPoRC's own
#          forms, and deviation penalties centered on the deviations' own mean.
#
#          Three stages: set every parameter to the pollock maximum likelihood
#          estimate and check the objective there, optimize, then compare.
# Creator: Matthew LH. Cheng
# Date Created: 8/7/26

library(here)
library(dplyr)
devtools::load_all(here())
source(here("dev", "make_sporc_obj_figs", "helper-bridge_figs.R"))

dat <- SPoRC::sgl_rg_ebswp_data
yrs <- dat$years
n_yrs <- length(yrs)
n_ages <- length(dat$ages)
n_srv <- dat$n_srv_fleets

i_bts <- which(yrs %in% dat$yrs_bts)
i_ats <- which(yrs %in% dat$yrs_ats)
i_avo <- which(yrs %in% dat$yrs_avo)

# Setup ----------------------------------------------------------------------
input_list <- Setup_Mod_Dim(
  years = yrs,
  ages = dat$ages,
  lens = NA,
  n_regions = dat$n_regions,
  n_sexes = dat$n_sexes,
  n_fish_fleets = dat$n_fish_fleets,
  n_srv_fleets = n_srv,
  n_seas = dat$n_seas,
  n_pop = dat$n_pop,
  natal_region = dat$natal_region,
  verbose = FALSE
)

inv_steepness <- function(s) qlogis((s - 0.2) / 0.8)

input_list <- Setup_Mod_Rec(
  input_list = input_list,
  # SR_ref_yr points spawning biomass per recruit at terminal year biologicals,
  # which is what the assessment's unfished spawning biomass uses. Built from the
  # first year instead, Bzero is 6021.33 rather than pm.tpl's 6109.5576.
  rec_model = "ricker_rec",
  rec_lag = 1,
  SR_ref_yr = n_yrs,
  # do_rec_bias_ramp = 0 sets the ramp to 1 throughout, so the recruitment
  # penalty is centered on -sigmaR^2/2, which is pm.tpl's +sigmaRsq/2 on the
  # residual. Setup_Mod_Rec takes bias_year, and an NA bias_year leaves the ramp
  # at zero and the bias correction absent altogether.
  do_rec_bias_ramp = 0,
  sigmaR_switch = 1,
  sigmaR_spec = "fix",
  # Three separate sigmas in pm.tpl, one per penalty. The stock recruit
  # residuals use sigr, fixed at 1. The initial ages carry weight 0.1 and the
  # recruitment level weight 1, which map onto sigma = 1/sqrt(2w).
  ln_sigmaR = array(c(log(1 / sqrt(0.2)), log(1)), dim = c(2, 1, 1)),
  init_age_strc = 4,
  equil_init_age_strc = 2,
  RecDevs_pen_center = "fixed",
  InitDevs_pen_center = "own_mean",
  ln_global_R0 = dat$mle$ln_global_R0,
  t_spawn = (4 - 1) / 12,
  use_rinit = 0,
  dont_est_recdev_last = 0,
  steepness_h = array(inv_steepness(dat$mle$steepness), dim = c(1, 1)),
  h_spec = "est_shared_pop_r",
  # pm.tpl's beta sits on the unrescaled (0,1) support and is symmetric, so its
  # center is 0.5 rather than the 0.6 its steepnessprior constant declares.
  Use_h_prior = 1,
  h_prior = data.frame(pop = 1, region = 1, mu = 0.5, sd = 0.09, lb = 0, ub = 1),
  # pm.tpl carries a second recruitment statement alongside the stock recruit
  # residuals: a sum of squares on all log recruitments about their own mean,
  # with weight 1, which is sigma = 1/sqrt(2).
  Use_rec_level_pen = 1,
  rec_level_pen_sigma = 1 / sqrt(2),
  rec_level_pen_center = "own_mean"
)

fix_natmort <- array(0, dim = c(1, 1, n_yrs, n_ages, 1))
fix_natmort[, , , 1, ] <- 0.9
fix_natmort[, , , 2, ] <- 0.45
fix_natmort[, , , -c(1, 2), ] <- 0.3

input_list <- Setup_Mod_Biologicals(
  input_list = input_list,
  WAA = dat$WAA,
  WAA_fish = dat$WAA_fish,
  WAA_srv = dat$WAA_srv,
  MatAA = dat$MatAA,
  fit_lengths = 0,
  M_spec = "fix",
  Fixed_natmort = fix_natmort,
  # pm.tpl weights the multinomial by the raw observed composition and puts
  # MN_const inside the logarithm only, which is comp_const_obs = 0. SPoRC's
  # default weights by obs + const, the unbiased choice since its stationary
  # point is p = obs, but that is a different model: on this bridge it moves
  # estimated spawning biomass by a median of 8.8 percent. addtocomp must stay
  # strictly positive, since at zero the offset term evaluates 0*log(0).
  addtocomp = 1e-3,
  comp_const_obs = 0,
  addtosrvidx = 0.01,
  addtofishidx = 0
)

input_list <- Setup_Mod_Movement(
  input_list = input_list,
  use_fixed_movement = 1,
  Fixed_Movement = NA,
  do_recruits_move = 0
)

# pm.tpl has no log_avg_F: fishing mortality is free annual log F, penalized
# about its own mean, which is ln_F_mean_spec = "fix" (mean pinned at zero, so
# the deviations carry all of log F) plus Fdev_pen_center = "own_mean".
input_list <- Setup_Mod_Catch_and_F(
  input_list = input_list,
  ObsCatch = dat$ObsCatch,
  UseCatch = dat$UseCatch,
  Use_F_pen = 1,
  Fdev_model = "iid",
  Fdev_pen_center = "own_mean",
  ln_F_mean_spec = "fix",
  sigmaF_spec = "fix",
  ln_sigmaF = array(log(1 / sqrt(2)), dim = c(1, 1, 1)),
  sigmaC_spec = "fix",
  ln_sigmaC = array(log(0.05), dim = c(1, n_yrs, 1, 1))
)

input_list <- Setup_Mod_FishIdx_and_Comps(
  input_list = input_list,
  ObsFishIdx = dat$ObsFishIdx,
  ObsFishIdx_SE = dat$ObsFishIdx_SE,
  UseFishIdx = dat$UseFishIdx,
  ObsFishAgeComps = dat$ObsFishAgeComps,
  UseFishAgeComps = dat$UseFishAgeComps,
  ISS_FishAgeComps = dat$ISS_FishAgeComps,
  ObsFishLenComps = array(NA_real_, dim = c(1, n_yrs, 1, length(input_list$data$lens), 1, 1)),
  UseFishLenComps = array(0, dim = c(1, n_yrs, 1, 1)),
  ISS_FishLenComps = array(0, dim = c(1, n_yrs, 1, 1, 1)),
  fish_idx_type = "biom",
  FishIdx_LikeType = "normal",
  FishAgeComps_LikeType = "Multinomial",
  FishLenComps_LikeType = "none",
  FishAgeComps_Type = "agg_Year_1-terminal_Fleet_1",
  FishLenComps_Type = "none_Year_1-terminal_Fleet_1"
)

input_list <- Setup_Mod_SrvIdx_and_Comps(
  input_list = input_list,
  ObsSrvIdx = dat$ObsSrvIdx,
  ObsSrvIdx_SE = dat$ObsSrvIdx_SE,
  UseSrvIdx = dat$UseSrvIdx,
  ObsSrvAgeComps = dat$ObsSrvAgeComps,
  ISS_SrvAgeComps = dat$ISS_SrvAgeComps,
  UseSrvAgeComps = dat$UseSrvAgeComps,
  ObsSrvLenComps = array(NA_real_, dim = c(1, n_yrs, 1, length(input_list$data$lens), 1, n_srv)),
  UseSrvLenComps = array(0, dim = c(1, n_yrs, 1, n_srv)),
  ISS_SrvLenComps = array(0, dim = c(1, n_yrs, 1, 1, n_srv)),
  srv_idx_type = c("biom", "biom", "biom", "abd"),
  # Fleet 4 is the acoustic survey's age 1 abundance, which the assessment fits
  # as its own index with its own catchability.
  srv_idx_ages = list(NULL, NULL, NULL, 1),
  # pm.tpl normalizes the acoustic compositions over ages 2-15 only
  SrvAgeComps_bins = list(NULL, 2:15, NULL, NULL),
  SrvIdx_LikeType = c("mvn", "lognormal", "normal", "lognormal"),
  SrvIdx_Cov = list(dat$SrvIdx_Cov, NULL, NULL, NULL),
  SrvAgeComps_LikeType = c("Multinomial", "Multinomial", "none", "none"),
  SrvLenComps_LikeType = rep("none", n_srv),
  SrvAgeComps_Type = c("agg_Year_1-terminal_Fleet_1", "agg_Year_1-terminal_Fleet_2",
                       "none_Year_1-terminal_Fleet_3", "none_Year_1-terminal_Fleet_4"),
  SrvLenComps_Type = paste0("none_Year_1-terminal_Fleet_", 1:n_srv),
  t_srv = array(c(0.5, 0.5, 0, 0.5), dim = c(1, 1, n_srv))
)

input_list <- Setup_Mod_Fishsel_and_Q(
  input_list = input_list,
  fish_sel_model = "nonparlog_Fleet_1",
  cont_tv_fish_sel = "rw_Fleet_1",
  fish_sel_blocks = "none_Fleet_1",
  fish_q_blocks = "none_Fleet_1",
  fish_fixed_sel_pars_spec = "est_all",
  fish_q_spec = "est_all",
  fishsel_pe_pars_spec = "fix",
  fish_sel_devs_spec = "est_all",
  # norm2(sel_devs_fsh) in pm.tpl weights every increment equally, including the
  # first, so the walk starts at zero under its own sigma rather than free.
  fishsel_rw_init_sigma = NA,
  # Nested [[fleet]][[block]] and then the bin groups. Ages 12-15 share one
  # coefficient, which is pm.tpl's flat tail; the deviations share the same
  # grouping so they are not estimated on bins with no free coefficient.
  fish_sel_nonpar_est_bins = list(list(list(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12:15))),
  fishsel_devs_shared_bins = list(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12:15)
)

input_list <- Setup_Mod_Srvsel_and_Q(
  input_list = input_list,
  srv_sel_model = c("logist1_Fleet_1", "nonparlog_Fleet_2", "nonparlog_Fleet_3", "logist1_Fleet_4"),
  cont_tv_srv_sel = c("iid_Fleet_1", "rw_Fleet_2", "rw_Fleet_3", "none_Fleet_4"),
  srv_sel_blocks = paste0("none_Fleet_", 1:n_srv),
  srv_q_blocks = paste0("none_Fleet_", 1:n_srv),
  # the vessel of opportunity index shares the acoustic survey's selectivity
  srv_fixed_sel_pars_spec = c("est_all", "est_all", "est_shared_f_2", "fix"),
  srv_sel_devs_spec = c("est_all", "est_all", "est_shared_f_2", "fix"),
  srvsel_pe_pars_spec = rep("fix", n_srv),
  srv_q_spec = rep("est_all", n_srv),
  # q enters as a multiplicative scalar, so it can be solved for instead of
  # estimated. "geo" is the geometric mean of the observed to predicted ratio,
  # which is the lognormal maximum likelihood value of q and is exact here
  # because the age 1 index has a constant standard error. "arith" is the ratio
  # of means, which is what pm.tpl applies to the trawl survey; it is not the
  # maximum likelihood value under that fleet's multivariate normal likelihood,
  # so it is a bridged convention rather than a concentrated parameter.
  srv_q_type = c("arith", "est", "est", "geo"),
  # A weight of zero makes the objective skip the process error likelihood for
  # ln_srvsel_devs, so the deviations stay estimated but are shaped only by the
  # explicit shape penalties set further down. pm.tpl constrains its survey
  # deviations that way rather than with a distribution, so keeping the default
  # of one would penalize them twice. It is also why srvsel_pe_pars_spec is
  # "fix": with the weight at zero the sigmas never reach the objective. The
  # bin override deviations are unaffected and keep their own process error.
  srvsel_pe_wt = c(0, 0, 0, 0),
  srv_sel_nonpar_est_bins = list(NULL,
                                 list(list(1, 2, 3, 4, 5, 6, 7, 8:15)),
                                 list(list(1:2, 3, 4, 5, 6, 7, 8:15)),
                                 NULL),
  srvsel_devs_shared_bins = list(1, 2, 3, 4, 5, 6, 7, 8:15),
  # the trawl survey's age 1 is a free annual value, not the logistic's
  srv_sel_bin_dev_bins = list(1, NULL, NULL, NULL),
  cont_tv_srvsel_bin_devs = c("rw", "none", "none", "none"),
  t_srv = array(c(0.5, 0.5, 0, 0.5), dim = c(1, 1, n_srv))
)

input_list <- Setup_Mod_Tagging(input_list = input_list, use_conv_fish_tagging = 0)

Wt_Rec <- array(0, dim = dim(input_list$par$ln_RecDevs))
Wt_Rec[1, 1, which(yrs %in% dat$yrs_srr)] <- 1
Wt_SrvIdx <- array(1, dim = c(1, n_yrs, 1, n_srv))
Wt_SrvIdx[1, max(i_ats), 1, 4] <- 0

# pm.tpl's selectivity penalties. The fishery carries a dome penalty over ages
# 6-12, curvature at the first year and every change year, and a random walk
# whose standard deviation is 0.5 except in two years where it opens up to 1.9.
# The trawl survey carries a year to year difference over ages 3-14, and the
# acoustic survey a dome penalty over ages 5-8.
yrs_ch_f <- dat$yrs_sel_ch_fsh
sig_ch <- rep(0.5, length(yrs_ch_f))
sig_ch[55:56] <- 1.9
curve_wt <- rep(0, n_yrs)
curve_wt[1] <- 1 / length(yrs_ch_f)
curve_wt[match(yrs_ch_f, yrs)] <- 1 / length(yrs_ch_f)
rw_wt <- rep(0, n_yrs)
rw_wt[match(yrs_ch_f, yrs)] <- 1 / (2 * sig_ch^2)

fish_pen_wts <- list(smooth_bin_diff = 3, smooth_bin_curve = curve_wt, smooth_yr_diff = rw_wt,
                     normalize = FALSE, bin_range = list(smooth_bin_diff = c(6, 12)))

# The trawl difference runs from 1982, whose predecessor pm.tpl fixes at one.
bts_rw <- rep(0, n_yrs)
bts_rw[match(1982:2024, yrs)] <- 2
yrs_ch_a <- dat$yrs_sel_ch_ats
ats_curve <- rep(0, n_yrs)
ats_curve[match(yrs_ch_a, yrs)] <- 1
ats_rw <- rep(0, n_yrs)
ats_rw[match(yrs_ch_a, yrs)] <- 1 / (2 * 0.138^2)
# pm.tpl applies the acoustic shape penalty only from the survey's first year,
# not across the whole model period.
ats_shape <- rep(0, n_yrs)
ats_shape[match(1994:2024, yrs)] <- 1
ats_spec <- list(smooth_bin_diff = ats_shape, smooth_bin_curve = ats_curve, smooth_yr_diff = ats_rw,
                 normalize = FALSE, bin_range = list(smooth_bin_diff = c(5, 8)))

srv_pen_wts <- list(
  list(smooth_yr_diff = bts_rw, normalize = FALSE, yr_diff_ref = 0,
       bin_range = list(smooth_yr_diff = c(3, 14))),
  ats_spec,
  # the vessel of opportunity index shares the acoustic curve, so it is not
  # penalized twice
  list(),
  list()
)

input_list <- Setup_Mod_Weighting(
  input_list = input_list,
  fish_sel_pen_wts = fish_pen_wts,
  srv_sel_pen_wts = srv_pen_wts,
  Wt_Catch = 1,
  Wt_FishIdx = 1,
  Wt_SrvIdx = Wt_SrvIdx,
  Wt_Rec = Wt_Rec,
  Wt_Init_Rec = 1,
  Wt_F = 1,
  Wt_Tagging = 0,
  Wt_FishAgeComps = array(1, dim = c(1, n_yrs, 1, 1, 1)),
  Wt_FishLenComps = array(1, dim = c(1, n_yrs, 1, 1, 1)),
  Wt_SrvAgeComps = array(1, dim = c(1, n_yrs, 1, 1, n_srv)),
  Wt_SrvLenComps = array(1, dim = c(1, n_yrs, 1, 1, n_srv))
)

data <- input_list$data
parameters <- input_list$par
mapping <- input_list$map

# Mapping --------------------------------------------------------------------
# Setup already leaves a parametric fleet's unused deviation slots and a
# no-time-variation fleet unmapped. What is left is the data window, since
# deviations before a survey exists have nothing to inform them, and the
# deviation groupings, which are levels here rather than pm.tpl's increments.
i_bts_all <- which(yrs >= 1982)
map_srvdev <- array(as.numeric(mapping$ln_srvsel_devs), dim = dim(parameters$ln_srvsel_devs))
map_srvdev[1, -i_bts_all, , 1, 1] <- NA
mapping$ln_srvsel_devs <- factor(map_srvdev)
data$map_ln_srvsel_devs <- map_srvdev

# The first year is redundant with the coefficients and is held at zero, and the
# bins within a group move together. Built explicitly because est_shared_b would
# share a group across years as well, making selectivity time invariant.
build_dev_map <- function(dim_arr, fleet, dev_years, groups) {
  m <- array(NA_real_, dim = dim_arr)
  k <- 1
  for(y in dev_years) {
    for(g in groups) {
      m[1, y, g, 1, fleet] <- k
      k <- k + 1
    }
  } # end y loop
  m
}

grp_f <- list(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12:15)
# pm.tpl applies sel_devs_fsh(ii) going from yrs_ch_fsh(ii) to the NEXT year, so
# the 59 changes land on 1966:2024 and log_sel_fsh is identical in 1964 and 1965.
# Mapping them a year early leaves 1965 estimated with a zero increment and 2024
# fixed, which drops the last increment from the walk.
map_fishdev <- build_dev_map(dim(parameters$ln_fishsel_devs), 1, match(1966:2024, yrs), grp_f)
mapping$ln_fishsel_devs <- factor(map_fishdev)
data$map_ln_fishsel_devs <- map_fishdev

# The acoustic survey's deviations sit on ages 2-8 with 9-15 following age 8;
# age 1 never receives one.
grp_a <- list(2, 3, 4, 5, 6, 7, 8:15)
map_ats <- array(as.numeric(mapping$ln_srvsel_devs), dim = dim(parameters$ln_srvsel_devs))
map_ats[1, , , 1, 2:3] <- NA
k <- max(map_ats, na.rm = TRUE) + 1
for(y in match(yrs_ch_a, yrs)) {
  for(g in grp_a) {
    map_ats[1, y, g, 1, 2] <- k
    k <- k + 1
  }
} # end y loop
# pm.tpl predicts the vessel of opportunity index from sel_ats itself, so that
# fleet shares the acoustic deviations, not just its coefficients. Leaving fleet
# 3 without any freezes its selectivity while fleet 2 moves, which is invisible
# at the pollock MLE and only appears after optimization.
map_ats[1, , , 1, 3] <- map_ats[1, , , 1, 2]
mapping$ln_srvsel_devs <- factor(map_ats)
data$map_ln_srvsel_devs <- map_ats

# pm.tpl holds acoustic age 1 log selectivity at zero rather than estimating it,
# so the objective is not stationary there unless it is fixed here too. It also
# has no base logistic parameters for the trawl survey: its slope and midpoint
# ARE the deviations, so the base is held at zero.
map_srvpar <- array(as.numeric(mapping$srv_fixed_sel_pars), dim = dim(parameters$srv_fixed_sel_pars))
map_srvpar[1, 1, 1, 1, 2:3] <- NA
map_srvpar[1, 1:2, 1, 1, 1] <- NA
mapping$srv_fixed_sel_pars <- factor(map_srvpar)

# pm.tpl fixes the age 1 series' penalty weight at 8 rather than estimating it.
mapping$srvsel_bin_devs_pe_pars <- factor(rep(NA, length(parameters$srvsel_bin_devs_pe_pars)))

map_bindev <- array(as.numeric(mapping$ln_srvsel_bin_devs), dim = dim(parameters$ln_srvsel_bin_devs))
map_bindev[1, -i_bts_all, , 1, 1] <- NA
mapping$ln_srvsel_bin_devs <- factor(map_bindev)
data$map_ln_srvsel_bin_devs <- map_bindev

# exp(par + dev) rescaled by its own mean is invariant to shifting par and dev
# together, so the level has to be pinned or the likelihood is flat along it.
# This is what pm.tpl's avgsel penalty does.
data$Use_fish_selex_penalty <- 1
fish_pen <- data.frame(region = 1, fleet = 1, block = 1, sex = 1, wt = 10)
fish_pen$par <- list(1:12)
data$fish_selex_penalty <- fish_pen

data$Use_srv_selex_penalty <- 1
srv_pen <- data.frame(region = 1, fleet = 2, block = 1, sex = 1, wt = 10)
srv_pen$par <- list(2:8)
data$srv_selex_penalty <- srv_pen

# Stage 1: set every parameter to the pollock MLE -----------------------------
# ln_F_mean_spec = "fix" already pinned the mean at zero and mapped it off, so
# the deviations are log F outright.
parameters$ln_F_devs[1, , 1, 1] <- log(dat$mle$Fmort)
parameters$ln_InitDevs[1, 1, ] <- dat$mle$log_initdevs

parameters$fish_fixed_sel_pars[1, 1:n_ages, 1, 1, 1] <- dat$mle$pars_fsh
parameters$ln_fishsel_devs[1, 1:n_yrs, , 1, 1] <- dat$mle$devs_fsh
for(f in 2:3) {
  parameters$srv_fixed_sel_pars[1, 1:n_ages, 1, 1, f] <- dat$mle$pars_ats
  parameters$ln_srvsel_devs[1, 1:n_yrs, , 1, f] <- dat$mle$devs_ats
} # end f loop
parameters$srv_fixed_sel_pars[1, 1:2, 1, 1, 1] <- 0
parameters$ln_srvsel_devs[1, 1:n_yrs, 1, 1, 1] <- dat$mle$bts_b50_dev
parameters$ln_srvsel_devs[1, 1:n_yrs, 2, 1, 1] <- dat$mle$bts_k_dev
parameters$ln_srvsel_bin_devs[1, 1:n_yrs, 1, 1, 1] <- dat$mle$bts_age1_dev

# pm.tpl penalizes the fishery selectivity increments with weight 1, which is
# sigma = 1/sqrt(2) on the random walk over the deviation levels, and the trawl
# age 1 series' first differences with weight 8.
parameters$fishsel_pe_pars[1, , 1, 1] <- log(1 / sqrt(2))
parameters$srvsel_bin_devs_pe_pars[1, , 1, 1] <- log(1 / sqrt(16))

obj <- fit_model(data, parameters, mapping, do_optim = FALSE, silent = TRUE)

cat("=== Stage 1: selectivity surfaces at the pollock MLE ===\n")
cat("fishery max pct diff:",
    100 * max(abs(obj$rep$fish_sel[1, 1, 1:n_yrs, 1, , 1, 1] - dat$admb$sel_fsh) / dat$admb$sel_fsh), "\n")
# sel_bts carries one row per year from 1982 to 2024 while there are only 42
# survey years, because 2020 has no trawl survey.
bts_rows <- match(yrs[i_bts], 1982:2024)
cat("trawl   max pct diff:",
    100 * max(abs(obj$rep$srv_sel[1, 1, i_bts, 1, , 1, 1] - dat$admb$sel_bts[bts_rows, ]) / dat$admb$sel_bts[bts_rows, ]), "\n")

# Recruitment deviations are stock recruit residuals here, and each year's
# residual depends on the previous year's spawning biomass, so they are solved
# by forward substitution rather than assigned.
free <- obj$par
idx_rec <- which(names(free) == "ln_RecDevs")
for(it in 1:40) {
  r <- obj$report(free)
  gap <- log(dat$mle$Rec) - log(as.vector(r$Rec[1, 1, 1:n_yrs]))
  if(max(abs(gap)) < 1e-12) break
  free[idx_rec] <- free[idx_rec] + gap
} # end it loop
cat("recruitment matched after", it, "passes; max abs log gap:", max(abs(gap)), "\n")

# The two estimated catchabilities are recovered from the ratio of the pollock
# predicted index to SPoRC's, which is exact because q enters multiplicatively.
# Fleets 1 and 4 solve theirs analytically and carry no free parameter.
r <- obj$report(free)
idx_q <- which(names(free) == "ln_srv_q")
free[idx_q] <- free[idx_q] + log(c(
  mean(dat$admb$eb_ats / as.vector(r$PredSrvIdx[1, 1, i_ats, 1, 2])),
  mean(dat$admb$pred_avo / as.vector(r$PredSrvIdx[1, 1, i_avo, 1, 3]))))
idx_fq <- which(names(free) == "ln_fish_q")
free[idx_fq] <- free[idx_fq] + log(mean(dat$admb$pred_cpue / as.vector(r$PredFishIdx[1, 1, 2:13, 1, 1])))

rep1 <- obj$report(free)
cat("\n=== Stage 1: population at the pollock MLE ===\n")
cat("N   max pct diff:", 100 * max(abs(rep1$NAA[1, 1, 1:n_yrs, 1, , 1] - dat$admb$NAA) / dat$admb$NAA), "\n")
cat("SSB max pct diff:", 100 * max(abs(rep1$SSB[1, 1, 1:n_yrs] - dat$admb$SSB) / dat$admb$SSB), "\n")
cat("jnLL at the pollock MLE:", obj$fn(free), "\n")
cat("max |gradient| there   :", max(abs(obj$gr(free))), "\n")

# Stage 2: optimize -----------------------------------------------------------
parameters2 <- obj$env$parList(free)
est <- fit_model(data, parameters2, mapping, do_optim = TRUE, newton_loops = 2, silent = TRUE)
cat("\n=== Stage 2: optimized ===\n")
cat("free parameters:", length(est$par), "\n")
cat("final jnLL:", est$optim$objective, "  max |gradient|:", max(abs(est$gr(est$optim$par))), "\n")
est$sdrep <- RTMB::sdreport(est)

# Stage 3: compare ------------------------------------------------------------
sdr <- est$sdrep
rep <- est$rep

ssb <- as.vector(rep$SSB[1, 1, 1:n_yrs])
rec <- as.vector(rep$Rec[1, 1, 1:n_yrs])

# The two series overplot exactly, so the difference gets its own panel, which is
# what bridge_ts_figure lays out.
ggplot2::ggsave(here("vignettes", "figures", "f_ebs_pol_ts_comparison.png"),
                bridge_ts_figure(yrs = yrs, ssb = ssb, rec = rec,
                                 admb_ssb = dat$admb$SSB, admb_rec = dat$admb$Rec,
                                 label = "2024 Pollock Assessment",
                                 ssb_se = bridge_se(sdr, "log_SSB", ssb, exact = TRUE),
                                 rec_se = bridge_se(sdr, "log_Rec", rec, exact = TRUE),
                                 legend_nrow = NULL),
                width = 17, height = 9, dpi = 150)

# Fishery selectivity, normalized to its maximum so the shape is readable.
sel_df <- reshape2::melt(rep$fish_sel[1, 1, 1:n_yrs, 1, , 1, 1]) %>%
  dplyr::rename(Year = Var1, Age = Var2) %>%
  dplyr::mutate(Year = yrs[Year], value = value / max(value), type = "SPoRC") %>%
  bind_rows(reshape2::melt(dat$admb$sel_fsh) %>%
              dplyr::rename(Year = Var1, Age = Var2) %>%
              dplyr::mutate(Year = yrs[Year], value = value / max(value),
                            type = "2024 Pollock Assessment")) %>%
  dplyr::filter(Age %in% 3:11)

p_sel <- ggplot2::ggplot(sel_df, ggplot2::aes(x = Year, y = value, color = type)) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::facet_wrap(~Age, labeller = ggplot2::label_both) +
  ggthemes::scale_color_colorblind() +
  ggplot2::theme_bw(base_size = 15) +
  ggplot2::theme(legend.position = "top") +
  ggplot2::labs(x = "Year", y = "Relative Selectivity", color = "Type")

ggplot2::ggsave(here("vignettes", "figures", "f_ebs_pol_fishsel_comparison.png"),
                p_sel, width = 12, height = 7, dpi = 150)

cat("\n=== Stage 3: optimized SPoRC against the pollock assessment ===\n")
print(rbind(bridge_cmp("SSB", ssb, dat$admb$SSB, signed = TRUE),
            bridge_cmp("Recruitment", rec, dat$admb$Rec, signed = TRUE)),
      row.names = FALSE, digits = 4)

saveRDS(list(rep = rep, sdrep = est$sdrep, opt = est$optim,
             data = data, parameters = parameters2, mapping = mapping),
        here("dev", "dev_output", "ebs_pollock_bridge.rds"))
