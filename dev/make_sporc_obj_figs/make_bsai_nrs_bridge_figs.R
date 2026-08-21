# Purpose: Bridge the 2024 BSAI northern rock sole assessment (flatfish model
#          fm.tpl, Model 24.2) to SPoRC from the packaged data object, and
#          render the case study figures. The model is specified the way the
#          assessment specifies it: two sexes with sex specific natural
#          mortality, free initial numbers at age estimated separately for each
#          sex about one shared level and tied to each other, mean recruitment
#          with a Ricker curve fitted as a penalty, time varying logistic
#          fishery selectivity with a male curve offset, and a survey read as a
#          July biomass index against January 1 compositions.
#
#          Three stages: set every parameter to the assessment's maximum
#          likelihood estimate and check the objective there, optimize, then
#          compare.
# Creator: Matthew LH. Cheng
# Date Created: 8/21/26

library(here)
library(dplyr)
devtools::load_all(here())
source(here("dev", "make_sporc_obj_figs", "helper-bridge_figs.R"))

dat <- SPoRC::sgl_rg_bsai_nrs_data
yrs <- dat$years
n_yrs <- length(yrs)
n_ages <- length(dat$ages)
n_sexes <- dat$n_sexes
n_srv <- dat$n_srv_fleets
mle <- dat$mle

i_srv <- match(dat$yrs_srv, yrs)
inv_steepness <- function(s) qlogis((s - 0.2) / 0.8)

# Setup ----------------------------------------------------------------------
input_list <- Setup_Mod_Dim(
  years = yrs, ages = dat$ages, lens = NA,
  n_regions = dat$n_regions, n_sexes = n_sexes,
  n_fish_fleets = dat$n_fish_fleets, n_srv_fleets = n_srv,
  n_seas = dat$n_seas, n_pop = dat$n_pop, natal_region = dat$natal_region,
  verbose = FALSE
)

# The Ricker in the assessment is R = A S exp(-B S). SPoRC writes it in
# depletion form, R = R0 (S/S0) exp(a (1 - S/S0)) with a = log(4h / (1 - h)),
# so the two are the same curve at a = log(A phi0) and R0 = a / (B phi0), with
# phi0 the unfished female spawning biomass per recruit at the reference year.
Nspr <- numeric(n_ages)
Nspr[1] <- 0.5
for(a in 2:n_ages) Nspr[a] <- Nspr[a-1] * exp(-mle$M_f)
Nspr[n_ages] <- Nspr[n_ages] / (1 - exp(-mle$M_f))
phi0 <- sum(Nspr * exp(-dat$t_spawn * mle$M_f) * dat$WAA[1,1,n_yrs,1,,1] * dat$MatAA[1,1,n_yrs,1,,1])
a_sr <- log(exp(mle$R_logalpha) * phi0)
h_sr <- exp(a_sr) / (4 + exp(a_sr))
sr_R0 <- a_sr / (exp(mle$R_logbeta) * phi0)

input_list <- Setup_Mod_Rec(
  input_list = input_list,
  rec_model = "mean_rec", rec_lag = 1, SR_ref_yr = n_yrs,
  sr_penalty = "ricker", sr_pen_sigma = dat$sr_pen_sigma,
  sr_pen_yrs = dat$sr_pen_yrs, sr_R0_spec = "est",
  steepness_h = array(inv_steepness(h_sr), dim = c(1, 1)), h_spec = "est_shared_pop_r",
  ln_sr_R0 = array(log(sr_R0), dim = 1),
  # the ramp is turned on with every break past the last year, which centres the
  # penalty on zero; do_rec_bias_ramp = 0 would centre it on -sigma^2/2
  do_rec_bias_ramp = 1, bias_year = rep(n_yrs + 1, 4),
  sigmaR_switch = 1, sigmaR_spec = "fix",
  ln_sigmaR = array(log(dat$sigmaR), dim = c(2, 1, 1)),
  RecDevs_pen_center = "fixed", dont_est_recdev_last = 0,
  init_age_strc = 4, equil_init_age_strc = 2,
  InitDevs_spec = NULL, InitDevs_sex_spec = "est_all", InitDevs_pen_center = "own_mean",
  Use_init_sex_pen = 1, init_sex_pen_sigma = dat$sigmaR,
  ln_global_R0 = mle$mean_log_rec, t_spawn = dat$t_spawn, use_rinit = 0
)

input_list <- Setup_Mod_Biologicals(
  input_list = input_list,
  WAA = dat$WAA, WAA_fish = dat$WAA_fish, WAA_srv = dat$WAA_srv, MatAA = dat$MatAA,
  AgeingError = dat$AgeingError, fit_lengths = 0,
  M_spec = "est_ln_M",
  M_popblk_spec = list(1), M_regionblk_spec = list(1), M_yearblk_spec = list(1:n_yrs),
  M_ageblk_spec = list(1:n_ages), M_sexblk_spec = list(1, 2),
  Use_M_prior = 1,
  M_prior = data.frame(popblk = 1, regionblk = 1, yearblk = 1, ageblk = 1, sexblk = 1,
                       mu = dat$m_prior$mu, sd = dat$m_prior$sd),
  addtocomp = 1e-3, comp_const_obs = 1, addtosrvidx = 0, addtofishidx = 0
)

input_list <- Setup_Mod_Movement(input_list = input_list, use_fixed_movement = 1,
                                 Fixed_Movement = NA, do_recruits_move = 0)

input_list <- Setup_Mod_Tagging(input_list = input_list, use_conv_fish_tagging = 0)

input_list <- Setup_Mod_Catch_and_F(
  input_list = input_list, ObsCatch = dat$ObsCatch, UseCatch = dat$UseCatch,
  Use_F_pen = 0, ln_F_mean_spec = "fix",
  sigmaC_spec = "fix", ln_sigmaC = array(log(dat$sigmaC), dim = c(1, n_yrs, 1, 1))
)

input_list <- Setup_Mod_FishIdx_and_Comps(
  input_list = input_list,
  ObsFishIdx = array(NA_real_, dim = c(1, n_yrs, 1, 1)),
  ObsFishIdx_SE = array(NA_real_, dim = c(1, n_yrs, 1, 1)),
  UseFishIdx = array(0, dim = c(1, n_yrs, 1, 1)),
  ObsFishAgeComps = dat$ObsFishAgeComps, UseFishAgeComps = dat$UseFishAgeComps,
  ISS_FishAgeComps = dat$ISS_FishAgeComps,
  ObsFishLenComps = array(NA_real_, dim = c(1, n_yrs, 1, length(input_list$data$lens), n_sexes, 1)),
  UseFishLenComps = array(0, dim = c(1, n_yrs, 1, 1)),
  ISS_FishLenComps = array(0, dim = c(1, n_yrs, 1, n_sexes, 1)),
  fish_idx_type = "none", FishIdx_LikeType = "lognormal",
  FishAgeComps_LikeType = "Multinomial", FishLenComps_LikeType = "none",
  FishAgeComps_Type = "spltRjntS_Year_1-terminal_Fleet_1",
  FishLenComps_Type = "none_Year_1-terminal_Fleet_1"
)

input_list <- Setup_Mod_SrvIdx_and_Comps(
  input_list = input_list,
  ObsSrvIdx = dat$ObsSrvIdx, ObsSrvIdx_SE = dat$ObsSrvIdx_SE, UseSrvIdx = dat$UseSrvIdx,
  ObsSrvAgeComps = dat$ObsSrvAgeComps, UseSrvAgeComps = dat$UseSrvAgeComps,
  ISS_SrvAgeComps = dat$ISS_SrvAgeComps,
  ObsSrvLenComps = array(NA_real_, dim = c(1, n_yrs, 1, length(input_list$data$lens), n_sexes, n_srv)),
  UseSrvLenComps = array(0, dim = c(1, n_yrs, 1, n_srv)),
  ISS_SrvLenComps = array(0, dim = c(1, n_yrs, 1, n_sexes, n_srv)),
  srv_idx_type = c("biom", "none"), SrvIdx_LikeType = rep("lognormal", n_srv),
  SrvAgeComps_LikeType = c("none", "Multinomial"), SrvLenComps_LikeType = rep("none", n_srv),
  SrvAgeComps_Type = c("none_Year_1-terminal_Fleet_1", "spltRjntS_Year_1-terminal_Fleet_2"),
  SrvLenComps_Type = paste0("none_Year_1-terminal_Fleet_", 1:n_srv),
  t_srv = array(dat$t_srv, dim = c(1, 1, n_srv))
)

input_list <- Setup_Mod_Fishsel_and_Q(
  input_list = input_list,
  fish_sel_model = paste0("logist1_Fleet_1_NSelBins_", dat$nselages),
  cont_tv_fish_sel = "iid_Fleet_1",
  fish_sel_blocks = "none_Fleet_1", fish_q_blocks = "none_Fleet_1",
  fish_fixed_sel_pars_spec = "est_all", fish_sel_devs_spec = "est_all",
  fish_sel_sex_offset = "scale",
  fishsel_pe_pars_spec = "fix", fish_q_spec = "fix"
)

input_list <- Setup_Mod_Srvsel_and_Q(
  input_list = input_list,
  srv_sel_model = paste0("logist1_Fleet_", 1:n_srv, "_NSelBins_", dat$nselages),
  cont_tv_srv_sel = paste0("none_Fleet_", 1:n_srv),
  srv_sel_blocks = paste0("none_Fleet_", 1:n_srv), srv_q_blocks = paste0("none_Fleet_", 1:n_srv),
  srv_fixed_sel_pars_spec = c("est_all", "est_shared_f_1"),
  srv_sel_sex_offset = rep("par", n_srv),
  srv_q_spec = c("est_all", "fix"),
  Use_srv_q_prior = 1,
  srv_q_prior = data.frame(region = 1, fleet = 1, block = 1, mu = dat$q_prior$mu, sd = dat$q_prior$sd),
  t_srv = array(dat$t_srv, dim = c(1, 1, n_srv))
)

input_list <- Setup_Mod_Weighting(
  input_list = input_list,
  Wt_Catch = 1, Wt_FishIdx = 0, Wt_SrvIdx = 1, Wt_Rec = 1, Wt_Init_Rec = 1, Wt_F = 1, Wt_Tagging = 0,
  Wt_FishAgeComps = array(1, dim = c(1, n_yrs, 1, n_sexes, 1)),
  Wt_FishLenComps = array(1, dim = c(1, n_yrs, 1, n_sexes, 1)),
  Wt_SrvAgeComps = array(1, dim = c(1, n_yrs, 1, n_sexes, n_srv)),
  Wt_SrvLenComps = array(1, dim = c(1, n_yrs, 1, n_sexes, n_srv))
)

data <- input_list$data
parameters <- input_list$par
mapping <- input_list$map

# Stage 1: start at the assessment's estimate --------------------------------
parameters$ln_F_mean[1,1,1] <- mle$log_avg_fmort
parameters$ln_F_devs[1,,1,1] <- mle$fmort_dev
parameters$ln_RecDevs[1,1,] <- mle$rec_dev
# free numbers at age absorb the assessment's shared mean_log_init
parameters$ln_InitDevs[1,1,,1] <- mle$mean_log_init + mle$init_dev_f
parameters$ln_InitDevs[1,1,,2] <- mle$mean_log_init + mle$init_dev_m
parameters$ln_M[1] <- log(mle$M_f)
parameters$ln_M[2] <- log(mle$M_m)
parameters$ln_srv_q[1,1,1] <- mle$ln_q

parameters$fish_fixed_sel_pars[1,1,1,1,1] <- log(mle$sel50_fsh_f)
parameters$fish_fixed_sel_pars[1,2,1,1,1] <- log(mle$sel_slope_fsh_f)
parameters$fish_fixed_sel_pars[1,1,1,2,1] <- log(mle$sel50_fsh_m)
parameters$fish_fixed_sel_pars[1,2,1,2,1] <- log(mle$sel_slope_fsh_m)
parameters$ln_fishsel_devs[1,,1,1,1] <- mle$sel50_devs_f
parameters$ln_fishsel_devs[1,,2,1,1] <- mle$slope_devs_f
parameters$ln_fishsel_devs[1,,1,2,1] <- mle$sel50_devs_m
parameters$ln_fishsel_devs[1,,2,2,1] <- mle$slope_devs_m
parameters$ln_fishsel_sex_scale[1,1,2,1] <- mle$male_sel_offset
parameters$fishsel_pe_pars[1,1,,1] <- log(dat$a50_sigma)
parameters$fishsel_pe_pars[1,2,,1] <- log(dat$slp_sigma)

# The two survey fleets share these through the map, and a collapsed parameter
# starts at the mean of its cells, so both fleets carry the value.
for(sf in 1:n_srv) {
  parameters$srv_fixed_sel_pars[1,1,1,1,sf] <- log(mle$sel50_srv)
  parameters$srv_fixed_sel_pars[1,2,1,1,sf] <- log(mle$sel_slope_srv)
  parameters$srv_fixed_sel_pars[1,1,1,2,sf] <- mle$sel50_srv_m
  parameters$srv_fixed_sel_pars[1,2,1,2,sf] <- mle$sel_slope_srv_m
} # end sf loop

seed <- fit_model(data, parameters, mapping, do_optim = FALSE, silent = TRUE)
rep1 <- seed$rep

cat("=== Stage 1: at the assessment's estimate ===\n")
pd <- function(a, b) 100 * (a - b) / b
totb_jan1 <- sapply(1:n_yrs, function(y) sum(rep1$NAA[1,1,y,1,,1] * dat$WAA[1,1,y,1,,1]) +
                      sum(rep1$NAA[1,1,y,1,,2] * dat$WAA[1,1,y,1,,2]))
for(row in list(c("spawning biomass", "SSB"), c("recruitment", "Rec"))) {
  v <- as.vector(rep1[[row[2]]][1,1,1:n_yrs])
  target <- if(row[2] == "SSB") dat$fm$SSB else dat$fm$Rec
  cat(sprintf("%-24s max %.3g %%\n", row[1], max(abs(pd(v, target)))))
} # end row loop
cat(sprintf("%-24s max %.3g %%\n", "total biomass (Jan 1)", max(abs(pd(totb_jan1, dat$fm$TotBiom)))))
cat(sprintf("%-24s max %.3g %%\n", "predicted catch", max(abs(pd(as.vector(rep1$PredCatch[1,1,1:n_yrs,1,1]), dat$fm$pred_catch)))))

# Likelihood crosswalk. SPoRC writes each component as a proper density while
# the assessment drops normalising constants, so each Gaussian block is compared
# net of exactly the constants the assessment omits.
lc <- function(sigma, n) n * (log(sigma) + 0.5 * log(2 * pi))
lik <- data.frame(
  component = c("survey index", "catch", "fishery ages", "survey ages", "rec devs",
                "init devs", "init sex tie", "SR (Ricker) penalty", "selex devs", "q prior", "M prior"),
  SPoRC_net = c(
    sum(rep1$SrvIdx_nLL) - sum(lc(dat$ObsSrvIdx_SE[1, i_srv, 1, 1], 1)),
    sum(rep1$Catch_nLL) - lc(dat$sigmaC, n_yrs),
    sum(rep1$FishAgeComps_nLL),
    sum(rep1$SrvAgeComps_nLL),
    sum(rep1$Rec_nLL) - lc(dat$sigmaR, n_yrs),
    sum(rep1$Init_Rec_nLL) - lc(dat$sigmaR, n_sexes * (n_ages - 1)),
    sum(rep1$Init_Sex_nLL) - lc(dat$sigmaR, n_ages - 1),
    sum(rep1$SR_pen_nLL) - lc(dat$sr_pen_sigma, sum(data$sr_pen_yrs)),
    rep1$sel_nLL - lc(dat$a50_sigma, n_sexes * n_yrs) - lc(dat$slp_sigma, n_sexes * n_yrs),
    rep1$srv_q_nLL - lc(dat$q_prior$sd, 1),
    rep1$M_nLL - lc(dat$m_prior$sd, 1)),
  assessment = with(dat$fm$Like_Comp,
                    c(srv, catch * 300, fsh_age, srv_age, rec, init, init_like, sr,
                      sel_a50 + sel_slope, q_prior, m_prior))
)
lik$difference <- lik$SPoRC_net - lik$assessment
cat("\n=== Stage 1: likelihood components ===\n")
print(lik, row.names = FALSE, digits = 10)
cat("\njnLL here:", seed$rep$jnLL, "  max |gradient|:", max(abs(seed$gr(seed$par))), "\n")

# Stage 2: optimize ----------------------------------------------------------
est <- fit_model(data, parameters, mapping, do_optim = TRUE, newton_loops = 3, silent = TRUE)
est$sdrep <- RTMB::sdreport(est, hessian.fixed = est$he(est$optim$par))
cat("\n=== Stage 2: optimized ===\n")
cat("free parameters:", length(est$par), "  final jnLL:", est$optim$objective,
    "  max |gradient|:", max(abs(est$gr(est$optim$par))), "\n")

# Stage 3: compare and render ------------------------------------------------
rep <- est$rep
ssb <- as.vector(rep$SSB[1,1,1:n_yrs])
rec <- as.vector(rep$Rec[1,1,1:n_yrs])

p_ts <- bridge_ts_figure(yrs = yrs, ssb = ssb, rec = rec,
                         admb_ssb = dat$fm$SSB, admb_rec = dat$fm$Rec,
                         label = "ADMB fm.tpl (Model 24.2)",
                         ssb_se = bridge_se(est$sdrep, "log_SSB", ssb),
                         rec_se = bridge_se(est$sdrep, "log_Rec", rec))
ggplot2::ggsave(here("vignettes", "figures", "ac_bsai_nrs_ts_comparison.png"), p_ts,
                width = 17, height = 9, dpi = 200)

# Selectivity in the terminal year, by sex and gear
adm_sel <- function(b50, k, scale = 1) {
  s <- 1 / (1 + exp(-k * ((1:n_ages) - b50)))
  s[(dat$nselages + 1):n_ages] <- s[dat$nselages]
  s * scale
}
sel_df <- dplyr::bind_rows(
  bridge_sel_rows(dat$ages, rep$fish_sel[1,1,n_yrs,1,,1,1],
                  adm_sel(mle$sel50_fsh_f * exp(mle$sel50_devs_f[n_yrs]), mle$sel_slope_fsh_f * exp(mle$slope_devs_f[n_yrs])),
                  "Fishery, female", "ADMB fm.tpl (Model 24.2)"),
  bridge_sel_rows(dat$ages, rep$fish_sel[1,1,n_yrs,1,,2,1],
                  adm_sel(mle$sel50_fsh_m * exp(mle$sel50_devs_m[n_yrs]), mle$sel_slope_fsh_m * exp(mle$slope_devs_m[n_yrs]), exp(mle$male_sel_offset)),
                  "Fishery, male", "ADMB fm.tpl (Model 24.2)"),
  bridge_sel_rows(dat$ages, rep$srv_sel[1,1,n_yrs,1,,1,1],
                  adm_sel(mle$sel50_srv, mle$sel_slope_srv),
                  "Survey, female", "ADMB fm.tpl (Model 24.2)"),
  bridge_sel_rows(dat$ages, rep$srv_sel[1,1,n_yrs,1,,2,1],
                  adm_sel(mle$sel50_srv * exp(mle$sel50_srv_m), mle$sel_slope_srv * exp(mle$sel_slope_srv_m)),
                  "Survey, male", "ADMB fm.tpl (Model 24.2)")
)
ggplot2::ggsave(here("vignettes", "figures", "ac_bsai_nrs_sel_comparison.png"),
                bridge_sel_figure(sel_df, nrow = 2), width = 15, height = 9, dpi = 200)

cat("\n=== Stage 3: optimized against the assessment ===\n")
print(rbind(bridge_cmp("Spawning biomass", ssb, dat$fm$SSB),
            bridge_cmp("Recruitment", rec, dat$fm$Rec)), row.names = FALSE, digits = 4)
pl <- est$env$parList(est$optim$par)
print(data.frame(
  quantity = c("M female", "M male", "survey q", "male fishery selectivity scale"),
  SPoRC = c(exp(pl$ln_M[1]), exp(pl$ln_M[2]), exp(pl$ln_srv_q[1,1,1]), exp(pl$ln_fishsel_sex_scale[1,1,2,1])),
  assessment = c(mle$M_f, mle$M_m, exp(mle$ln_q), exp(mle$male_sel_offset))
), row.names = FALSE, digits = 5)
cat("\nfigures written to vignettes/figures\n")
