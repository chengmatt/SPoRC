# Purpose: Bridge SAM's Northeast Arctic cod assessment to SPoRC, self test it, and render its figures
# Creator: Matthew LH. Cheng
# Date Created: 10/5/26

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
devtools::load_all(here())
source(here("dev", "make_sporc_obj_figs", "helper-bridge_figs.R"))

d <- sgl_rg_neacod_data$inputs
sam <- sgl_rg_neacod_data$sam
years <- d$years
ages <- d$ages
n_yrs <- d$n_yrs
n_ages <- d$n_ages
n_srv <- d$n_srv_fleets
n_sims <- 50 # replicates in the conditioned self test
n_sims_joint <- 150 # replicates in the joint self test, enough to read a bias of about one percent

# Build -------------------------------------------------------------------------

input_list <- Setup_Mod_Dim(
  n_pop = 1,
  years = years,
  ages = ages,
  lens = NA,
  n_regions = 1,
  n_sexes = 1,
  n_seas = 1,
  n_fish_fleets = 1,
  n_srv_fleets = n_srv,
  verbose = FALSE
)

input_list <- Setup_Mod_Rec(
  input_list,
  rec_model = "mean_rec",
  RecDevs_model = "rw",
  dont_pen_recdev_first = 1,
  sigmaR_spec = "fix_early_est_late",
  sigmaR_switch = 1,
  do_rec_bias_ramp = 0,
  init_age_strc = "free",
  equil_init_age_strc = "stoch_all_no_pen",
  t_spawn = 0,
  ln_global_R0_spec = "fix",
  ln_global_R0 = sam$logN[1, 1]
)

input_list <- Setup_Mod_Biologicals(
  input_list,
  WAA = d$WAA,
  WAA_fish = d$WAA_fish,
  WAA_srv = d$WAA_srv,
  MatAA = d$MatAA,
  fit_lengths = 0,
  M_spec = "fix",
  Fixed_natmort = d$natmort,
  NAA_re = "iid",
  NAA_re_ages = ages[-1],
  NAA_re_years = years[-1],
  NAA_sigma_spec = "est",
  NAA_sigma_ageblk_spec_vals = list(seq_len(n_ages))
)

input_list <- Setup_Mod_Movement(input_list, use_fixed_movement = 1, do_recruits_move = 0, Fixed_Movement = NA)
input_list <- Setup_Mod_Tagging(input_list, use_conv_fish_tagging = 0)

# missing aggregate catch keeps f estimated every year, the catch at age is the fishery likelihood
input_list <- Setup_Mod_Catch_and_F(
  input_list,
  ObsCatch = array(NA, dim = c(1, n_yrs, 1, 1)),
  UseCatch = array(0, dim = c(1, n_yrs, 1, 1)),
  catch_units = "abd",
  ObsCatchAA = d$ObsCatchAA,
  UseCatchAA = d$UseCatchAA,
  CatchAA_LikeType = "lognormal",
  CatchAA_Type = "spltRaggS",
  sigmaCAA_key = d$sigmaCAA_key,
  sigmaCAA_spec = "est",
  sigmaC_spec = "fix",
  ln_F_mean_spec = "fix",
  Fdev_model = "rw",
  sigmaF_spec = "est_all",
  Use_F_pen = 1
)

input_list <- Setup_Mod_FishIdx_and_Comps(
  input_list,
  ObsFishIdx = array(1, dim = c(1, n_yrs, 1, 1)),
  ObsFishIdx_SE = array(0.1, dim = c(1, n_yrs, 1, 1)),
  UseFishIdx = array(0, dim = c(1, n_yrs, 1, 1)),
  ObsFishAgeComps = array(0, dim = c(1, n_yrs, 1, n_ages, 1, 1)),
  UseFishAgeComps = array(0, dim = c(1, n_yrs, 1, 1)),
  ISS_FishAgeComps = array(0, dim = c(1, n_yrs, 1, 1, 1)),
  ObsFishLenComps = array(0, dim = c(1, n_yrs, 1, 1, 1, 1)),
  UseFishLenComps = array(0, dim = c(1, n_yrs, 1, 1)),
  ISS_FishLenComps = array(0, dim = c(1, n_yrs, 1, 1, 1)),
  fish_idx_type = "none",
  FishAgeComps_LikeType = "none",
  FishLenComps_LikeType = "none",
  FishAgeComps_Type = "none_Year_1-terminal_Fleet_1",
  FishLenComps_Type = "none_Year_1-terminal_Fleet_1"
)

input_list <- Setup_Mod_SrvIdx_and_Comps(
  input_list,
  ObsSrvIdx = array(1, dim = c(1, n_yrs, 1, n_srv)),
  ObsSrvIdx_SE = array(0.1, dim = c(1, n_yrs, 1, n_srv)),
  UseSrvIdx = array(0, dim = c(1, n_yrs, 1, n_srv)),
  ObsSrvIdxAA = d$ObsSrvIdxAA,
  UseSrvIdxAA = d$UseSrvIdxAA,
  SrvIdxAA_LikeType = "lognormal",
  SrvIdxAA_Type = "spltRaggS",
  sigmaSrvIdxAA_key = d$sigmaSrvIdxAA_key,
  sigmaSrvIdxAA_spec = "est",
  ObsSrvAgeComps = array(0, dim = c(1, n_yrs, 1, n_ages, 1, n_srv)),
  UseSrvAgeComps = array(0, dim = c(1, n_yrs, 1, n_srv)),
  ISS_SrvAgeComps = array(0, dim = c(1, n_yrs, 1, 1, n_srv)),
  ObsSrvLenComps = array(0, dim = c(1, n_yrs, 1, 1, 1, n_srv)),
  UseSrvLenComps = array(0, dim = c(1, n_yrs, 1, n_srv)),
  ISS_SrvLenComps = array(0, dim = c(1, n_yrs, 1, 1, n_srv)),
  srv_idx_type = c("none", "none", "none", "none", "none"),
  SrvAgeComps_LikeType = c("none", "none", "none", "none", "none"),
  SrvLenComps_LikeType = c("none", "none", "none", "none", "none"),
  SrvAgeComps_Type = c("none_Year_1-terminal_Fleet_1", "none_Year_1-terminal_Fleet_2", "none_Year_1-terminal_Fleet_3",
                       "none_Year_1-terminal_Fleet_4", "none_Year_1-terminal_Fleet_5"),
  SrvLenComps_Type = c("none_Year_1-terminal_Fleet_1", "none_Year_1-terminal_Fleet_2", "none_Year_1-terminal_Fleet_3",
                       "none_Year_1-terminal_Fleet_4", "none_Year_1-terminal_Fleet_5")
)

# log f at age, ages 14 and 15+ sharing one value and one walk
fish_bin_groups <- list(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12:13)
input_list <- Setup_Mod_Fishsel_and_Q(
  input_list,
  fish_sel_model = "nonparfree_Fleet_1",
  fish_sel_blocks = "none_Fleet_1",
  fish_q_blocks = "none_Fleet_1",
  fish_fixed_sel_pars_spec = "est_all",
  fish_q_spec = "fix",
  fish_sel_nonpar_est_bins = list(list(fish_bin_groups)),
  cont_tv_fish_sel = "rw_Fleet_1",
  fishsel_pe_pars_spec = "est_shared_b",
  fish_sel_devs_spec = "est_shared_b",
  fishsel_devs_shared_bins = fish_bin_groups,
  fishsel_dont_est_dev_first = 1
)

# survey catchability at age, ages 11 and up sharing one value
srv_bin_groups <- list(1, 2, 3, 4, 5, 6, 7, 8, 9:13)
input_list <- Setup_Mod_Srvsel_and_Q(
  input_list,
  srv_sel_model = c("nonparfree_Fleet_1", "nonparfree_Fleet_2", "nonparfree_Fleet_3", "nonparfree_Fleet_4", "nonparfree_Fleet_5"),
  srv_sel_blocks = c("none_Fleet_1", "none_Fleet_2", "none_Fleet_3", "none_Fleet_4", "none_Fleet_5"),
  srv_q_blocks = c("none_Fleet_1", "none_Fleet_2", "none_Fleet_3", "none_Fleet_4", "none_Fleet_5"),
  srv_fixed_sel_pars_spec = c("est_all", "est_all", "est_all", "est_all", "est_all"),
  srv_q_spec = c("fix", "fix", "fix", "fix", "fix"),
  cont_tv_srv_sel = c("none_Fleet_1", "none_Fleet_2", "none_Fleet_3", "none_Fleet_4", "none_Fleet_5"),
  srv_sel_nonpar_est_bins = list(list(srv_bin_groups), list(srv_bin_groups), list(srv_bin_groups), list(srv_bin_groups), list(srv_bin_groups)),
  t_srv = array(d$srv_time, dim = c(1, 1, n_srv))
)

input_list <- Setup_Mod_Weighting(
  input_list,
  Wt_Catch = 1,
  Wt_FishIdx = 1,
  Wt_SrvIdx = 1,
  Wt_Rec = 1,
  Wt_F = 1,
  Wt_Tagging = 0,
  addtosrvidx = 0,
  addtofishidx = 0,
  Wt_FishAgeComps = array(1, dim = c(1, n_yrs, 1, 1, 1)),
  Wt_FishLenComps = array(1, dim = c(1, n_yrs, 1, 1, 1)),
  Wt_SrvAgeComps = array(1, dim = c(1, n_yrs, 1, 1, n_srv)),
  Wt_SrvLenComps = array(1, dim = c(1, n_yrs, 1, 1, n_srv))
)

# Start at SAM's estimate -----------------------------------------------------------

# year one's log f at age in the base values, later years as their difference from it
input_list$par$ln_F_mean[] <- 0
input_list$par$ln_F_devs[] <- 0
input_list$par$fish_fixed_sel_pars[1, , 1, 1, 1] <- sam$logFF[1, ]
for(y in 2:n_yrs) input_list$par$ln_fishsel_devs[1, y, , 1, 1] <- sam$logFF[y, ] - sam$logFF[1, ]

# the f increment sd split into the part shared across ages and the part each age has alone
sd_F <- exp(sam$logsdF[1])
input_list$par$ln_sigmaF[] <- log(sd_F * sqrt(sam$rho))
input_list$par$fishsel_pe_pars[1, , 1, 1] <- log(sd_F * sqrt(1 - sam$rho))

# survey catchability at age, unobserved ages reading the oldest observed one
input_list$par$ln_srv_q[] <- 0
for(sf in 1:n_srv) {
  q_idx <- d$srv_q_key[, sf]
  q_idx[is.na(q_idx)] <- max(q_idx, na.rm = TRUE)
  input_list$par$srv_fixed_sel_pars[1, , 1, 1, sf] <- sam$logQ[q_idx]
} # end sf loop

# observation and process sds, then the states
input_list$par$ln_sigmaCAA[] <- sam$logSdLogObs[1]
for(sf in 1:n_srv) input_list$par$ln_sigmaSrvIdxAA[, 1, sf] <- sam$logSdLogObs[1 + sf]
input_list$par$ln_sigmaNAA[] <- sam$logSdLogN[2]
input_list$par$ln_sigmaR[2, 1, 1] <- sam$logSdLogN[1]
input_list$par$ln_InitDevs[1, 1, , 1] <- sam$logN[1, 2:n_ages]
input_list$par$ln_NAA[1, 1, , 1, , 1] <- sam$logN
input_list$par$ln_RecDevs[] <- 0
for(y in 2:n_yrs) input_list$par$ln_RecDevs[1, 1, y] <- sam$logN[y, 1] - sam$logN[1, 1]

# Fit -------------------------------------------------------------------------------

# year one's log f at age is integrated out too, as SAM integrates its first f state
re_names <- c("ln_NAA", "ln_RecDevs", "ln_InitDevs", "ln_fishsel_devs", "ln_F_devs", "fish_fixed_sel_pars")
fit <- fit_model(input_list$data, input_list$par, input_list$map, random = re_names,
                 do_optim = TRUE, newton_loops = 3, silent = TRUE)
sd_rep <- RTMB::sdreport(fit, getJointPrecision = TRUE)
rep <- fit$rep
best <- fit$env$last.par.best[-fit$env$random]
bridge_nll <- c(sporc = as.numeric(fit$fn(best)), sam = sam$nll, max_grad = max(abs(fit$gr(best))))

# Bridge figures ----------------------------------------------------------------------

ggplot2::ggsave(
  here("vignettes", "figures", "ai_neacod_ts.png"),
  bridge_ts_figure(years, rep$SSB[1, 1, ], rep$Rec[1, 1, ], sam$ssb, exp(sam$logN[, 1]), "SAM", ref_name = "SAM"),
  width = 17,
  height = 9,
  dpi = 150
)

# f at age, six ages
faa_ages <- c(3, 5, 7, 9, 11, 13)
faa_df <- bind_rows(
  data.frame(Year = rep(years, n_ages), Age = rep(ages, each = n_yrs), value = as.vector(rep$tot_FAA[1, 1, , 1, , 1, 1]), type = "SPoRC"),
  data.frame(Year = rep(years, n_ages), Age = rep(ages, each = n_yrs), value = as.vector(exp(sam$logFF)), type = "SAM")
) %>% filter(Age %in% faa_ages) %>% mutate(Age = paste("Age", Age))

p_faa <- ggplot(faa_df, aes(Year, value, color = type, linetype = type)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~Age, ncol = 3) +
  ggthemes::scale_color_colorblind() +
  theme_bw(base_size = 18) +
  theme(legend.position = "top") +
  labs(x = "Year", y = "Fishing mortality", color = "Type", linetype = "Type")
ggplot2::ggsave(here("vignettes", "figures", "ai_neacod_faa.png"), p_faa, width = 14, height = 8, dpi = 150)

# Self tests --------------------------------------------------------------------------

what <- c("SSB", "Rec", "tot_FAA")
what_par <- c("ln_sigmaCAA", "ln_sigmaSrvIdxAA", "ln_sigmaNAA", "ln_sigmaR", "ln_sigmaF", "fishsel_pe_pars")

# conditioned on the fit's processes, so only the observations are drawn
set.seed(123)
st_cond <- simulation_self_test(data = fit$data, parameters = input_list$par, mapping = input_list$map, random = re_names,
                                rep = rep, sd_rep = sd_rep, obj = fit, n_sims = n_sims, newton_loops = 3,
                                what = what, what_par = what_par)

# parameters drawn from the fit's joint precision, then every process fresh at them
set.seed(456)
st_joint <- simulation_self_test(data = fit$data, parameters = input_list$par, mapping = input_list$map, random = re_names,
                                 rep = rep, sd_rep = sd_rep, obj = fit, n_sims = n_sims_joint, newton_loops = 3, sim_type = "joint",
                                 what = what, what_par = what_par)

# perfect data, conditioned on the fit's processes
set.seed(789)
st_perfect <- simulation_self_test(data = fit$data, parameters = input_list$par, mapping = input_list$map, random = re_names,
                                   rep = rep, sd_rep = sd_rep, obj = fit, n_sims = 3, newton_loops = 3, perfect_data = TRUE,
                                   what = what, what_par = what_par)

# perfect data under the joint design, so each population is one the fit has not seen
set.seed(1011)
st_perfect_joint <- simulation_self_test(data = fit$data, parameters = input_list$par, mapping = input_list$map, random = re_names,
                                         rep = rep, sd_rep = sd_rep, obj = fit, n_sims = 5, newton_loops = 3, perfect_data = TRUE,
                                         sim_type = "joint", what = what, what_par = what_par)

saveRDS(list(bridge_nll = bridge_nll, st_cond = st_cond, st_joint = st_joint, st_perfect = st_perfect, st_perfect_joint = st_perfect_joint),
        here("dev", "dev_output", "nea_cod_selftest.rds"))

# Self test figures -------------------------------------------------------------------

# relative error by year: spawning biomass, recruitment and mean f over ages 5-10
fbar_ages <- which(ages %in% 5:10)
designs <- list("Conditional" = st_cond, "Joint" = st_joint)
err_df <- NULL
for(design in names(designs)) {
  st <- designs[[design]]
  n_rep <- dim(st$SSB)[4]
  fbar_est <- apply(st$tot_FAA[1, 1, , 1, fbar_ages, 1, 1, , drop = FALSE], c(3, 8), mean)
  fbar_true <- apply(st$truth$tot_FAA[1, 1, , 1, fbar_ages, 1, 1, , drop = FALSE], c(3, 8), mean)
  err_df <- bind_rows(
    err_df,
    data.frame(Design = design, Par = "Spawning biomass", Year = rep(years, n_rep), value = as.vector(st$SSB[1, 1, , ] / st$truth$SSB[1, 1, , ] - 1)),
    data.frame(Design = design, Par = "Recruitment", Year = rep(years, n_rep), value = as.vector(st$Rec[1, 1, , ] / st$truth$Rec[1, 1, , ] - 1)),
    data.frame(Design = design, Par = "Mean F, ages 5-10", Year = rep(years, n_rep), value = as.vector(fbar_est / fbar_true - 1))
  )
} # end design loop

err_sum <- err_df %>%
  group_by(Design, Par, Year) %>%
  summarize(lwr = quantile(value, 0.05, na.rm = TRUE), med = median(value, na.rm = TRUE), upr = quantile(value, 0.95, na.rm = TRUE), .groups = "drop")

p_err <- ggplot(err_sum, aes(Year, 100 * med, ymin = 100 * lwr, ymax = 100 * upr, color = Design, fill = Design)) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_ribbon(alpha = 0.25, color = NA) +
  geom_line(linewidth = 0.9) +
  facet_grid(Par ~ Design) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  theme_bw(base_size = 16) +
  theme(legend.position = "none") +
  labs(x = "Year", y = "Relative error (%)")
ggplot2::ggsave(here("vignettes", "figures", "ai_neacod_selftest.png"), p_err, width = 14, height = 10, dpi = 150)

# observation and process sds, each replicate's estimate over the value it ran on: one catch sd, one per survey
# (read at its youngest age), the numbers at age sd, the late recruitment sd, and the two f walk sds
sd_df <- NULL
for(design in names(designs)) {
  st <- designs[[design]]
  n_rep <- dim(st$SSB)[4]
  caa_est <- array(st$ln_sigmaCAA, dim = c(n_ages, n_rep)) # age by replicate
  caa_true <- array(st$truth$ln_sigmaCAA, dim = c(n_ages, n_rep))
  srv_est <- array(st$ln_sigmaSrvIdxAA, dim = c(n_ages, n_srv, n_rep)) # age by survey by replicate
  srv_true <- array(st$truth$ln_sigmaSrvIdxAA, dim = c(n_ages, n_srv, n_rep))
  naa_est <- array(st$ln_sigmaNAA, dim = c(length(st$ln_sigmaNAA) / n_rep, n_rep))
  naa_true <- array(st$truth$ln_sigmaNAA, dim = dim(naa_est))
  rec_est <- array(st$ln_sigmaR, dim = c(2, n_rep)) # early and late
  rec_true <- array(st$truth$ln_sigmaR, dim = c(2, n_rep))
  f_est <- array(st$ln_sigmaF, dim = c(length(st$ln_sigmaF) / n_rep, n_rep))
  f_true <- array(st$truth$ln_sigmaF, dim = dim(f_est))
  fsel_est <- array(st$fishsel_pe_pars, dim = c(length(st$fishsel_pe_pars) / n_rep, n_rep))
  fsel_true <- array(st$truth$fishsel_pe_pars, dim = dim(fsel_est))

  sd_df <- bind_rows(sd_df, data.frame(Design = design, Par = "Catch at age", value = exp(caa_est[1, ] - caa_true[1, ])))
  for(sf in 1:n_srv) sd_df <- bind_rows(sd_df, data.frame(Design = design, Par = paste("Survey", sf), value = exp(srv_est[1, sf, ] - srv_true[1, sf, ])))

  # a conditional replicate runs on the fit's smoothed states, not on draws at the fitted process sds, so those are left out
  if(design == "Conditional") next
  sd_df <- bind_rows(
    sd_df,
    data.frame(Design = design, Par = "Numbers at age process", value = exp(naa_est[1, ] - naa_true[1, ])),
    data.frame(Design = design, Par = "Recruitment process", value = exp(rec_est[2, ] - rec_true[2, ])),
    data.frame(Design = design, Par = "F level process", value = exp(f_est[1, ] - f_true[1, ])),
    data.frame(Design = design, Par = "F at age process", value = exp(fsel_est[1, ] - fsel_true[1, ]))
  )
} # end design loop
sd_df$Par <- factor(sd_df$Par, levels = c("Catch at age", "Survey 1", "Survey 2", "Survey 3", "Survey 4", "Survey 5",
                                          "Numbers at age process", "Recruitment process", "F level process", "F at age process"))

p_sd <- ggplot(sd_df, aes(Par, 100 * (value - 1), fill = Design)) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_boxplot(outlier.size = 0.8) +
  ggthemes::scale_fill_colorblind() +
  theme_bw(base_size = 16) +
  theme(legend.position = "top", axis.text.x = element_text(angle = 35, hjust = 1)) +
  labs(x = NULL, y = "Relative error in the sd (%)", fill = NULL)
ggplot2::ggsave(here("vignettes", "figures", "ai_neacod_selftest_sd.png"), p_sd, width = 14, height = 8, dpi = 150)
