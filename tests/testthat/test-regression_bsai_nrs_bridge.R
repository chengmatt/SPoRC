# Self-validating bridge test: the expectations are the 2024 BSAI northern rock sole assessment's own
# reported quantities in sgl_rg_bsai_nrs_data$fm, evaluated at its estimate without optimizing.

library(SPoRC)
library(testthat)
data("sgl_rg_bsai_nrs_data")

test_that("BSAI northern rock sole bridges to the 2024 ADMB assessment at its own MLE", {

  dat <- sgl_rg_bsai_nrs_data
  yrs <- dat$years
  n_yrs <- length(yrs)
  n_ages <- length(dat$ages)
  n_sexes <- dat$n_sexes
  n_srv <- dat$n_srv_fleets
  mle <- dat$mle
  i_srv <- match(dat$yrs_srv, yrs)

  inv_steepness <- function(s) qlogis((s - 0.2) / 0.8)

  input_list <- Setup_Mod_Dim(
    years = yrs,
    ages = dat$ages,
    lens = NA,
    n_regions = dat$n_regions,
    n_sexes = n_sexes,
    n_fish_fleets = dat$n_fish_fleets,
    n_srv_fleets = n_srv,
    n_seas = dat$n_seas,
    n_pop = dat$n_pop,
    natal_region = dat$natal_region,
    verbose = FALSE
  )

  # the assessment's Ricker, R = A S exp(-B S), in SPoRC's depletion form
  Nspr <- numeric(n_ages)
  Nspr[1] <- 0.5
  for(a in 2:n_ages) Nspr[a] <- Nspr[a-1] * exp(-mle$M_f)
  Nspr[n_ages] <- Nspr[n_ages] / (1 - exp(-mle$M_f))
  phi0 <- sum(Nspr * exp(-dat$t_spawn * mle$M_f) * dat$WAA[1,1,n_yrs,1,,1] * dat$MatAA[1,1,n_yrs,1,,1])
  a_sr <- log(exp(mle$R_logalpha) * phi0)
  h_sr <- exp(a_sr) / (4 + exp(a_sr))

  input_list <- Setup_Mod_Rec(
    input_list = input_list,
    rec_model = "mean_rec",
    rec_lag = 1,
    SR_ref_yr = n_yrs,
    sr_penalty = "ricker",
    sr_pen_sigma = dat$sr_pen_sigma,
    sr_pen_yrs = dat$sr_pen_yrs,
    sr_R0_spec = "est",
    steepness_h = array(inv_steepness(h_sr), dim = c(1, 1)),
    h_spec = "est_shared_pop_r",
    ln_sr_R0 = array(log(a_sr / (exp(mle$R_logbeta) * phi0)), dim = 1),
    do_rec_bias_ramp = 1,
    bias_year = rep(n_yrs + 1, 4),
    sigmaR_switch = 1,
    sigmaR_spec = "fix",
    ln_sigmaR = array(log(dat$sigmaR), dim = c(2, 1, 1)),
    RecDevs_pen_center = "fixed",
    dont_est_recdev_last = 0,
    init_age_strc = 4,
    equil_init_age_strc = 2,
    InitDevs_spec = NULL,
    InitDevs_sex_spec = "est_all",
    InitDevs_pen_center = "own_mean",
    Use_init_sex_pen = 1,
    init_sex_pen_sigma = dat$sigmaR,
    ln_global_R0 = mle$mean_log_rec,
    t_spawn = dat$t_spawn,
    use_rinit = 0
  )

  input_list <- Setup_Mod_Biologicals(
    input_list = input_list,
    WAA = dat$WAA,
    WAA_fish = dat$WAA_fish,
    WAA_srv = dat$WAA_srv,
    MatAA = dat$MatAA,
    AgeingError = dat$AgeingError,
    fit_lengths = 0,
    M_spec = "est_ln_M",
    M_popblk_spec = list(1),
    M_regionblk_spec = list(1),
    M_yearblk_spec = list(1:n_yrs),
    M_ageblk_spec = list(1:n_ages),
    M_sexblk_spec = list(1, 2),
    Use_M_prior = 1,
    M_prior = data.frame(
      popblk = 1,
      regionblk = 1,
      yearblk = 1,
      ageblk = 1,
      sexblk = 1,
      mu = dat$m_prior$mu,
      sd = dat$m_prior$sd
    ),
    addtocomp = 1e-3,
    comp_const_obs = 1,
    addtosrvidx = 0,
    addtofishidx = 0
  )

  input_list <- Setup_Mod_Movement(
    input_list = input_list,
    use_fixed_movement = 1,
    Fixed_Movement = NA,
    do_recruits_move = 0
  )
  input_list <- Setup_Mod_Tagging(input_list = input_list, use_conv_fish_tagging = 0)

  suppressWarnings(input_list <- Setup_Mod_Catch_and_F(
    input_list = input_list,
    ObsCatch = dat$ObsCatch,
    UseCatch = dat$UseCatch,
    Use_F_pen = 0,
    ln_F_mean_spec = "fix",
    sigmaC_spec = "fix",
    ln_sigmaC = array(log(dat$sigmaC), dim = c(1, n_yrs, 1, 1))
  ))

  input_list <- Setup_Mod_FishIdx_and_Comps(
    input_list = input_list,
    ObsFishIdx = array(NA_real_, dim = c(1, n_yrs, 1, 1)),
    ObsFishIdx_SE = array(NA_real_, dim = c(1, n_yrs, 1, 1)),
    UseFishIdx = array(0, dim = c(1, n_yrs, 1, 1)),
    ObsFishAgeComps = dat$ObsFishAgeComps,
    UseFishAgeComps = dat$UseFishAgeComps,
    ISS_FishAgeComps = dat$ISS_FishAgeComps,
    ObsFishLenComps = array(NA_real_, dim = c(1, n_yrs, 1, length(input_list$data$lens), n_sexes, 1)),
    UseFishLenComps = array(0, dim = c(1, n_yrs, 1, 1)),
    ISS_FishLenComps = array(0, dim = c(1, n_yrs, 1, n_sexes, 1)),
    fish_idx_type = "none",
    FishIdx_LikeType = "lognormal",
    FishAgeComps_LikeType = "Multinomial",
    FishLenComps_LikeType = "none",
    FishAgeComps_Type = "spltRjntS_Year_1-terminal_Fleet_1",
    FishLenComps_Type = "none_Year_1-terminal_Fleet_1"
  )

  input_list <- Setup_Mod_SrvIdx_and_Comps(
    input_list = input_list,
    ObsSrvIdx = dat$ObsSrvIdx,
    ObsSrvIdx_SE = dat$ObsSrvIdx_SE,
    UseSrvIdx = dat$UseSrvIdx,
    ObsSrvAgeComps = dat$ObsSrvAgeComps,
    UseSrvAgeComps = dat$UseSrvAgeComps,
    ISS_SrvAgeComps = dat$ISS_SrvAgeComps,
    ObsSrvLenComps = array(NA_real_, dim = c(1, n_yrs, 1, length(input_list$data$lens), n_sexes, n_srv)),
    UseSrvLenComps = array(0, dim = c(1, n_yrs, 1, n_srv)),
    ISS_SrvLenComps = array(0, dim = c(1, n_yrs, 1, n_sexes, n_srv)),
    srv_idx_type = c("biom", "none"),
    SrvIdx_LikeType = rep("lognormal", n_srv),
    SrvAgeComps_LikeType = c("none", "Multinomial"),
    SrvLenComps_LikeType = rep("none", n_srv),
    SrvAgeComps_Type = c("none_Year_1-terminal_Fleet_1", "spltRjntS_Year_1-terminal_Fleet_2"),
    SrvLenComps_Type = paste0("none_Year_1-terminal_Fleet_", 1:n_srv),
    t_srv = array(dat$t_srv, dim = c(1, 1, n_srv))
  )

  input_list <- Setup_Mod_Fishsel_and_Q(
    input_list = input_list,
    fish_sel_model = paste0("logist1_Fleet_1_NSelBins_", dat$nselages),
    cont_tv_fish_sel = "iid_Fleet_1",
    fish_sel_blocks = "none_Fleet_1",
    fish_q_blocks = "none_Fleet_1",
    fish_fixed_sel_pars_spec = "est_all",
    fish_sel_devs_spec = "est_all",
    fish_sel_sex_offset = "scale",
    fishsel_pe_pars_spec = "fix",
    fish_q_spec = "fix"
  )

  input_list <- Setup_Mod_Srvsel_and_Q(
    input_list = input_list,
    srv_sel_model = paste0("logist1_Fleet_", 1:n_srv, "_NSelBins_", dat$nselages),
    cont_tv_srv_sel = paste0("none_Fleet_", 1:n_srv),
    srv_sel_blocks = paste0("none_Fleet_", 1:n_srv),
    srv_q_blocks = paste0("none_Fleet_", 1:n_srv),
    srv_fixed_sel_pars_spec = c("est_all", "est_shared_f_1"),
    srv_sel_sex_offset = rep("par", n_srv),
    srv_q_spec = c("est_all", "fix"),
    Use_srv_q_prior = 1,
    srv_q_prior = data.frame(region = 1, fleet = 1, block = 1, mu = dat$q_prior$mu, sd = dat$q_prior$sd),
    t_srv = array(dat$t_srv, dim = c(1, 1, n_srv))
  )

  input_list <- Setup_Mod_Weighting(
    input_list = input_list,
    Wt_Catch = 1,
    Wt_FishIdx = 0,
    Wt_SrvIdx = 1,
    Wt_Rec = 1,
    Wt_Init_Rec = 1,
    Wt_F = 1,
    Wt_Tagging = 0,
    Wt_FishAgeComps = array(1, dim = c(1, n_yrs, 1, n_sexes, 1)),
    Wt_FishLenComps = array(1, dim = c(1, n_yrs, 1, n_sexes, 1)),
    Wt_SrvAgeComps = array(1, dim = c(1, n_yrs, 1, n_sexes, n_srv)),
    Wt_SrvLenComps = array(1, dim = c(1, n_yrs, 1, n_sexes, n_srv))
  )

  data <- input_list$data
  parameters <- input_list$par
  mapping <- input_list$map

  # Set every parameter to the assessment's MLE ----
  parameters$ln_F_mean[1,1,1] <- mle$log_avg_fmort
  parameters$ln_F_devs[1,,1,1] <- mle$fmort_dev
  parameters$ln_RecDevs[1,1,] <- mle$rec_dev
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
  for(sf in 1:n_srv) {
    parameters$srv_fixed_sel_pars[1,1,1,1,sf] <- log(mle$sel50_srv)
    parameters$srv_fixed_sel_pars[1,2,1,1,sf] <- log(mle$sel_slope_srv)
    parameters$srv_fixed_sel_pars[1,1,1,2,sf] <- mle$sel50_srv_m
    parameters$srv_fixed_sel_pars[1,2,1,2,sf] <- mle$sel_slope_srv_m
  } # end sf loop

  obj <- fit_model(data, parameters, mapping, do_optim = FALSE, silent = TRUE)
  rep <- obj$rep

  # Population, against the assessment's reported series ----
  pd <- function(a, b) max(abs(100 * (a - b) / b))
  expect_lt(pd(as.vector(rep$SSB[1,1,1:n_yrs]), dat$fm$SSB), 0.01)
  expect_lt(pd(as.vector(rep$Rec[1,1,1:n_yrs]), dat$fm$Rec), 0.01)
  expect_lt(pd(as.vector(rep$PredCatch[1,1,1:n_yrs,1,1]), dat$fm$pred_catch), 0.01)
  totb_jan1 <- sapply(1:n_yrs, function(y) sum(rep$NAA[1,1,y,1,,1] * dat$WAA[1,1,y,1,,1]) +
                        sum(rep$NAA[1,1,y,1,,2] * dat$WAA[1,1,y,1,,2]))
  expect_lt(pd(totb_jan1, dat$fm$TotBiom), 0.01)

  # Selectivity, against the assessment's own form including its age-17 plateau ----
  adm_sel <- function(b50, k, scale = 1) {
    s <- 1 / (1 + exp(-k * ((1:n_ages) - b50)))
    s[(dat$nselages + 1):n_ages] <- s[dat$nselages]
    s * scale
  }
  sel_f <- t(sapply(1:n_yrs, function(y) adm_sel(mle$sel50_fsh_f * exp(mle$sel50_devs_f[y]),
                                                 mle$sel_slope_fsh_f * exp(mle$slope_devs_f[y]))))
  sel_m <- t(sapply(1:n_yrs, function(y) adm_sel(mle$sel50_fsh_m * exp(mle$sel50_devs_m[y]),
                                                 mle$sel_slope_fsh_m * exp(mle$slope_devs_m[y]),
                                                 exp(mle$male_sel_offset))))
  expect_lt(max(abs(rep$fish_sel[1,1,1:n_yrs,1,,1,1] - sel_f)), 1e-12)
  expect_lt(max(abs(rep$fish_sel[1,1,1:n_yrs,1,,2,1] - sel_m)), 1e-12)
  expect_lt(max(abs(rep$srv_sel[1,1,1,1,,1,1] - adm_sel(mle$sel50_srv, mle$sel_slope_srv))), 1e-12)
  expect_lt(max(abs(rep$srv_sel[1,1,1,1,,2,1] - adm_sel(mle$sel50_srv * exp(mle$sel50_srv_m),
                                                        mle$sel_slope_srv * exp(mle$sel_slope_srv_m)))), 1e-12)
  # both survey fleets read the one curve
  expect_equal(as.numeric(rep$srv_sel[1,1,1,1,,1,2]), as.numeric(rep$srv_sel[1,1,1,1,,1,1]), tolerance = 1e-12)

  # Likelihood components, net of the constants the assessment omits ----
  lc <- function(sigma, n) n * (log(sigma) + 0.5 * log(2 * pi))
  cmp <- dat$fm$Like_Comp
  near_admb <- function(sporc, admb) expect_lt(abs(sporc - admb) / max(abs(admb), 1e-8), 1e-5)
  near_admb(sum(rep$SrvIdx_nLL) - sum(lc(dat$ObsSrvIdx_SE[1, i_srv, 1, 1], 1)), cmp$srv)
  near_admb(sum(rep$Catch_nLL) - lc(dat$sigmaC, n_yrs), cmp$catch * 300)
  near_admb(sum(rep$FishAgeComps_nLL), cmp$fsh_age)
  near_admb(sum(rep$SrvAgeComps_nLL), cmp$srv_age)
  near_admb(sum(rep$Rec_nLL) - lc(dat$sigmaR, n_yrs), cmp$rec)
  near_admb(sum(rep$Init_Rec_nLL) - lc(dat$sigmaR, n_sexes * (n_ages - 1)), cmp$init)
  near_admb(sum(rep$Init_Sex_nLL) - lc(dat$sigmaR, n_ages - 1), cmp$init_like)
  near_admb(sum(rep$SR_pen_nLL) - lc(dat$sr_pen_sigma, sum(data$sr_pen_yrs)), cmp$sr)
  near_admb(rep$sel_nLL - lc(dat$a50_sigma, n_sexes * n_yrs) - lc(dat$slp_sigma, n_sexes * n_yrs),
            cmp$sel_a50 + cmp$sel_slope)
  near_admb(rep$srv_q_nLL - lc(dat$q_prior$sd, 1), cmp$q_prior)
  near_admb(rep$M_nLL - lc(dat$m_prior$sd, 1), cmp$m_prior)
})
