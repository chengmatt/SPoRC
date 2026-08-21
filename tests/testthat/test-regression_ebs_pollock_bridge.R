# Self-validating bridge test. The expectations are not stored SPoRC output: they are
# the 2024 EBS walleye pollock assessment's own reported quantities, shipped in
# sgl_rg_ebswp_data$admb. Every parameter is set to the assessment's maximum likelihood
# estimate and the model is evaluated there without optimizing, so a failure means the
# population dynamics, a likelihood, or a selectivity form no longer reproduces the
# assessment at a point where it is known to. Do not loosen the tolerances to make this
# pass. See tests/README.md and the EBS pollock case study vignette.

library(SPoRC)
library(testthat)
data("sgl_rg_ebswp_data")

test_that("EBS pollock bridges exactly to the 2024 ADMB assessment at its own MLE", {

  dat <- sgl_rg_ebswp_data
  yrs <- dat$years
  n_yrs <- length(yrs)
  n_ages <- length(dat$ages)
  n_srv <- dat$n_srv_fleets

  i_bts <- which(yrs %in% dat$yrs_bts)
  i_ats <- which(yrs %in% dat$yrs_ats)
  i_avo <- which(yrs %in% dat$yrs_avo)

  # Dimensions ----
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

  # Recruitment: Ricker with a one year lag, steepness estimated under a beta
  # prior, free initial numbers at age, and three separate penalty variances.
  input_list <- Setup_Mod_Rec(
    input_list = input_list,
    rec_model = "ricker_rec",
    rec_lag = 1,
    SR_ref_yr = n_yrs,
    do_rec_bias_ramp = 0,
    sigmaR_switch = 1,
    sigmaR_spec = "fix",
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
    Use_h_prior = 1,
    h_prior = data.frame(pop = 1, region = 1, mu = 0.5, sd = 0.09, lb = 0, ub = 1),
    Use_rec_level_pen = 1,
    rec_level_pen_sigma = 1 / sqrt(2),
    rec_level_pen_center = "own_mean"
  )

  # Biologicals: age specific fixed M, separate weight at age for spawning
  # biomass, catch and each survey, and the assessment's multinomial constant convention.
  fix_natmort <- array(0, dim = c(1, 1, n_yrs, n_ages, 1))
  fix_natmort[, , , 1, ] <- 0.9
  fix_natmort[, , , 2, ] <- 0.45
  fix_natmort[, , , -c(1, 2), ] <- 0.3

  suppressWarnings(
    input_list <- Setup_Mod_Biologicals(
      input_list = input_list,
      WAA = dat$WAA,
      WAA_fish = dat$WAA_fish,
      WAA_srv = dat$WAA_srv,
      MatAA = dat$MatAA,
      fit_lengths = 0,
      M_spec = "fix",
      Fixed_natmort = fix_natmort,
      addtocomp = 1e-3,
      comp_const_obs = 0,
      addtosrvidx = 0.01,
      addtofishidx = 0
    )
  )

  input_list <- Setup_Mod_Movement(
    input_list = input_list,
    use_fixed_movement = 1,
    Fixed_Movement = NA,
    do_recruits_move = 0
  )

  # Catch and F: no mean F parameter, so the deviations carry all of log F and
  # are penalised about their own mean.
  suppressWarnings(
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

  # Four survey indices with four different error structures. Fleet 4 is the
  # acoustic survey's age 1 abundance, which the assessment fits as its own
  # index with its own catchability.
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
    srv_idx_ages = list(NULL, NULL, NULL, 1),
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

  # Fishery selectivity: non-parametric on the log scale, ages 12-15 sharing one
  # coefficient, evolving as a random walk that starts at zero.
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
    fishsel_rw_init_sigma = NA,
    fish_sel_nonpar_est_bins = list(list(list(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12:15))),
    fishsel_devs_shared_bins = list(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12:15)
  )

  input_list <- Setup_Mod_Srvsel_and_Q(
    input_list = input_list,
    srv_sel_model = c("logist1_Fleet_1", "nonparlog_Fleet_2", "nonparlog_Fleet_3", "logist1_Fleet_4"),
    cont_tv_srv_sel = c("iid_Fleet_1", "rw_Fleet_2", "rw_Fleet_3", "none_Fleet_4"),
    srv_sel_blocks = paste0("none_Fleet_", 1:n_srv),
    srv_q_blocks = paste0("none_Fleet_", 1:n_srv),
    srv_fixed_sel_pars_spec = c("est_all", "est_all", "est_shared_f_2", "fix"),
    srv_sel_devs_spec = c("est_all", "est_all", "est_shared_f_2", "fix"),
    srvsel_pe_pars_spec = rep("fix", n_srv),
    srv_q_spec = rep("est_all", n_srv),
    srv_q_type = c("arith", "est", "est", "geo"),
    # A weight of zero makes the objective skip the process error likelihood for
    # ln_srvsel_devs, so the deviations stay estimated but are shaped only by the
    # explicit shape penalties set further down. pm.tpl constrains its survey
    # deviations that way rather than with a distribution, so the default of one
    # would penalize them twice and the bridge would not match. It is also why
    # srvsel_pe_pars_spec is "fix": with the weight at zero the sigmas never
    # reach the objective. The bin override deviations keep their own process
    # error, which this weight does not touch.
    srvsel_pe_wt = c(0, 0, 0, 0),
    srv_sel_nonpar_est_bins = list(NULL,
                                   list(list(1, 2, 3, 4, 5, 6, 7, 8:15)),
                                   list(list(1:2, 3, 4, 5, 6, 7, 8:15)),
                                   NULL),
    srvsel_devs_shared_bins = list(1, 2, 3, 4, 5, 6, 7, 8:15),
    srv_sel_bin_dev_bins = list(1, NULL, NULL, NULL),
    cont_tv_srvsel_bin_devs = c("rw", "none", "none", "none"),
    t_srv = array(c(0.5, 0.5, 0, 0.5), dim = c(1, 1, n_srv))
  )

  input_list <- Setup_Mod_Tagging(input_list = input_list, use_conv_fish_tagging = 0)

  # Weighting and the assessment's selectivity shape penalties ----
  Wt_Rec <- array(0, dim = dim(input_list$par$ln_RecDevs))
  Wt_Rec[1, 1, which(yrs %in% dat$yrs_srr)] <- 1
  Wt_SrvIdx <- array(1, dim = c(1, n_yrs, 1, n_srv))
  Wt_SrvIdx[1, max(i_ats), 1, 4] <- 0

  yrs_ch_f <- dat$yrs_sel_ch_fsh
  sig_ch <- rep(0.5, length(yrs_ch_f))
  sig_ch[55:56] <- 1.9
  curve_wt <- rep(0, n_yrs)
  curve_wt[1] <- 1 / length(yrs_ch_f)
  curve_wt[match(yrs_ch_f, yrs)] <- 1 / length(yrs_ch_f)
  rw_wt <- rep(0, n_yrs)
  rw_wt[match(yrs_ch_f, yrs)] <- 1 / (2 * sig_ch^2)

  fish_pen_wts <- list(smooth_bin_diff = 3, smooth_bin_curve = curve_wt,
                       smooth_yr_diff = rw_wt, normalize = FALSE,
                       bin_range = list(smooth_bin_diff = c(6, 12)))

  bts_rw <- rep(0, n_yrs)
  bts_rw[match(1982:2024, yrs)] <- 2
  yrs_ch_a <- dat$yrs_sel_ch_ats
  ats_curve <- rep(0, n_yrs)
  ats_curve[match(yrs_ch_a, yrs)] <- 1
  ats_rw <- rep(0, n_yrs)
  ats_rw[match(yrs_ch_a, yrs)] <- 1 / (2 * 0.138^2)
  ats_shape <- rep(0, n_yrs)
  ats_shape[match(1994:2024, yrs)] <- 1

  srv_pen_wts <- list(
    list(smooth_yr_diff = bts_rw, normalize = FALSE, yr_diff_ref = 0,
         bin_range = list(smooth_yr_diff = c(3, 14))),
    list(smooth_bin_diff = ats_shape, smooth_bin_curve = ats_curve,
         smooth_yr_diff = ats_rw, normalize = FALSE,
         bin_range = list(smooth_bin_diff = c(5, 8))),
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

  # Mapping. SPoRC parameterises selectivity deviations as levels while the
  # assessment uses increments, so the first year's level is redundant with the
  # coefficients and is held at zero, and bins within a group move together.
  i_bts_all <- which(yrs >= 1982)
  map_srvdev <- array(as.numeric(mapping$ln_srvsel_devs), dim = dim(parameters$ln_srvsel_devs))
  map_srvdev[1, -i_bts_all, , 1, 1] <- NA
  mapping$ln_srvsel_devs <- factor(map_srvdev)
  data$map_ln_srvsel_devs <- map_srvdev

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
  map_fishdev <- build_dev_map(dim(parameters$ln_fishsel_devs), 1, match(1966:2024, yrs), grp_f)
  mapping$ln_fishsel_devs <- factor(map_fishdev)
  data$map_ln_fishsel_devs <- map_fishdev

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
  map_ats[1, , , 1, 3] <- map_ats[1, , , 1, 2]
  mapping$ln_srvsel_devs <- factor(map_ats)
  data$map_ln_srvsel_devs <- map_ats

  map_srvpar <- array(as.numeric(mapping$srv_fixed_sel_pars), dim = dim(parameters$srv_fixed_sel_pars))
  map_srvpar[1, 1, 1, 1, 2:3] <- NA
  map_srvpar[1, 1:2, 1, 1, 1] <- NA
  mapping$srv_fixed_sel_pars <- factor(map_srvpar)

  mapping$srvsel_bin_devs_pe_pars <- factor(rep(NA, length(parameters$srvsel_bin_devs_pe_pars)))

  map_bindev <- array(as.numeric(mapping$ln_srvsel_bin_devs), dim = dim(parameters$ln_srvsel_bin_devs))
  map_bindev[1, -i_bts_all, , 1, 1] <- NA
  mapping$ln_srvsel_bin_devs <- factor(map_bindev)
  data$map_ln_srvsel_bin_devs <- map_bindev

  # exp(par + dev) rescaled by its own mean is invariant to shifting par and dev
  # together, so the level is pinned. This is the assessment's avgsel penalty.
  data$Use_fish_selex_penalty <- 1
  fish_pen <- data.frame(region = 1, fleet = 1, block = 1, sex = 1, wt = 10)
  fish_pen$par <- list(1:12)
  data$fish_selex_penalty <- fish_pen

  data$Use_srv_selex_penalty <- 1
  srv_pen <- data.frame(region = 1, fleet = 2, block = 1, sex = 1, wt = 10)
  srv_pen$par <- list(2:8)
  data$srv_selex_penalty <- srv_pen

  # Set every parameter to the assessment's MLE ----
  parameters$ln_F_devs[1, , 1, 1] <- log(dat$mle$Fmort)
  parameters$ln_InitDevs[1, 1, , ] <- dat$mle$log_initdevs

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

  parameters$fishsel_pe_pars[1, , 1, 1] <- log(1 / sqrt(2))
  parameters$srvsel_bin_devs_pe_pars[1, , 1, 1] <- log(1 / sqrt(16))

  obj <- fit_model(data, parameters, mapping, do_optim = FALSE, silent = TRUE)

  # Selectivity is built by SPoRC's own forms, so matching the assessment's
  # surfaces is a check on the forms rather than on the data.
  expect_equal(obj$rep$fish_sel[1, 1, 1:n_yrs, 1, , 1, 1], dat$admb$sel_fsh,
               tolerance = 1e-10, ignore_attr = TRUE)
  # sel_bts carries one row per year from 1982 to 2024, while there are only 42
  # survey years, because 2020 has no trawl survey.
  bts_rows <- match(yrs[i_bts], 1982:2024)
  expect_equal(obj$rep$srv_sel[1, 1, i_bts, 1, , 1, 1], dat$admb$sel_bts[bts_rows, ],
               tolerance = 1e-10, ignore_attr = TRUE)

  # Recruitment deviations are stock recruit residuals, and each year's residual
  # depends on the previous year's spawning biomass, so they are solved by
  # forward substitution rather than assigned.
  free <- obj$par
  idx_rec <- which(names(free) == "ln_RecDevs")
  for(it in 1:40) {
    r <- obj$report(free)
    gap <- log(dat$mle$Rec) - log(as.vector(r$Rec[1, 1, 1:n_yrs]))
    if(max(abs(gap)) < 1e-12) break
    free[idx_rec] <- free[idx_rec] + gap
  } # end it loop
  expect_lt(max(abs(gap)), 1e-12)

  # Catchability enters multiplicatively, so it is recovered exactly from the
  # ratio of the assessment's predicted index to SPoRC's.
  r <- obj$report(free)
  idx_q <- which(names(free) == "ln_srv_q")
  free[idx_q] <- free[idx_q] + log(c(
    mean(dat$admb$eb_ats / as.vector(r$PredSrvIdx[1, 1, i_ats, 1, 2])),
    mean(dat$admb$pred_avo / as.vector(r$PredSrvIdx[1, 1, i_avo, 1, 3]))))
  idx_fq <- which(names(free) == "ln_fish_q")
  free[idx_fq] <- free[idx_fq] + log(mean(dat$admb$pred_cpue / as.vector(r$PredFishIdx[1, 1, 2:13, 1, 1])))

  rep1 <- obj$report(free)

  # The population at the assessment's MLE ----
  expect_equal(rep1$NAA[1, 1, 1:n_yrs, 1, , 1], dat$admb$NAA,
               tolerance = 1e-9, ignore_attr = TRUE)
  expect_equal(as.vector(rep1$SSB[1, 1, 1:n_yrs]), as.vector(dat$admb$SSB),
               tolerance = 1e-9, ignore_attr = TRUE)

  # The assessment is at its own optimum, so SPoRC's gradient must vanish there
  # too. A specification difference shows up here before it shows up in SSB.
  expect_lt(max(abs(obj$gr(free))), 1e-3)

  # Pinned only because the assessment and SPoRC drop different normalising constants, so
  # the totals differ by a constant rather than matching outright.
  expect_equal(obj$fn(free), 1547.851, tolerance = 1e-5)

  expect_jnLL_decomposes(obj)
})
