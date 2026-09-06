# Self-test of the sex-linked selectivity options and sex-specific initial age deviations, end to end.
#
# A two-sex OM whose truth IS the offset structure (male fishery curve scaled by exp(-0.163), survey male
# parameters as log offsets, a plateau from bin 4, per-sex initial age curves) generates joint-sex data.
#
# The estimation model holds that structure through fish_sel_sex_offset / srv_sel_sex_offset / NSelBins /
# InitDevs_sex_spec.
#
# simulation_self_test then conditions the OM on the fit, so the fitted surfaces and deviations travel back
# through Setup_Sim_Fishing / Setup_Sim_Survey / Setup_Sim_Rec.

library(testthat)
library(RTMB)

test_that("sex offsets, the selectivity plateau, and per-sex initial deviations recover through the simulation", {

  n_yrs <- 25; n_ages <- 6; n_sexes <- 2; n_sims <- 10
  waa <- 5 / (1 + exp(-3 * ((1:n_ages) - 3)))
  mat <- 1 / (1 + exp(-3 * ((1:n_ages) - 3)))
  f_ramp <- c(seq(0.05, 0.45, length.out = 15), seq(0.45, 0.15, length.out = 10))

  # truth: the male fishery curve is the male logistic times exp(gamma) with a
  # plateau from bin 4; survey male parameters are log offsets on the female's
  fish_scale_true <- -0.163
  srv_off_true <- c(-0.2, 0.25) # on ln_b50, ln_k
  plat <- function(s) { s[5:n_ages] <- s[4]; s }
  logi <- function(b50, k) 1 / (1 + exp(-k * ((1:n_ages) - b50)))
  fsel_f <- plat(logi(2.5, 2.0)); fsel_m <- plat(logi(2.9, 1.6)) * exp(fish_scale_true)
  ssel_f <- logi(3.0, 1.5); ssel_m <- logi(3.0 * exp(srv_off_true[1]), 1.5 * exp(srv_off_true[2]))
  init_devs_true <- cbind(c(0.25, -0.2, 0.15, -0.1, 0.05), c(-0.15, 0.25, -0.1, 0.2, -0.05))

  sel_arr <- function(f_curve, m_curve) {
    a <- array(0, dim = c(1, 1, n_yrs, 1, n_ages, n_sexes, 1))
    for(y in 1:n_yrs) { a[1,1,y,1,,1,1] <- f_curve; a[1,1,y,1,,2,1] <- m_curve }
    a
  }

  sim_list <- Setup_Sim_Dim(
    n_sims = 1,
    n_yrs = n_yrs,
    n_regions = 1,
    n_ages = n_ages,
    n_lens = NULL,
    n_sexes = n_sexes,
    n_fish_fleets = 1,
    n_srv_fleets = 1,
    n_pop = 1
  )
  sim_list <- Setup_Sim_Containers(sim_list)
  sim_list <- Setup_Sim_Fishing(
    sim_list = sim_list,
    fish_sel_input = replicate(1, sel_arr(fsel_f, fsel_m)),
    ret_sel_input = replicate(1, array(1, dim = c(1, 1, n_yrs, 1, n_ages, n_sexes, 1))),
    dmr_input = array(0, dim = c(1, n_yrs, 1, 1, 1)),
    Fmort_input = array(rep(f_ramp, each = 1), dim = c(1, n_yrs, 1, 1, 1)),
    ISS_FishAgeComps = array(1000, dim = c(1, n_yrs, 1, n_sexes, 1, 1)),
    FishAgeComps_Type = matrix("spltRjntS", n_yrs, 1)
  )
  sim_list <- Setup_Sim_Survey(
    sim_list = sim_list,
    srv_sel_input = replicate(1, sel_arr(ssel_f, ssel_m)),
    ObsSrvIdx_SE = array(0.05, dim = c(1, n_yrs, 1, 1)),
    ISS_SrvAgeComps = array(1000, dim = c(1, n_yrs, 1, n_sexes, 1, 1)),
    SrvAgeComps_Type = matrix("spltRjntS", n_yrs, 1)
  )
  biol6 <- function(val) {
    a <- array(0, dim = c(1, 1, n_yrs, 1, n_ages, n_sexes))
    for(s in 1:n_sexes) a[1,1,,1,,s] <- matrix(rep(val, each = n_yrs), n_yrs, n_ages)
    a
  }
  waa_fleet <- array(0, dim = c(1, 1, n_yrs, 1, n_ages, n_sexes, 1)); waa_fleet[1,1,,1,,,1] <- biol6(waa)[1,1,,1,,]
  suppressWarnings(sim_list <- Setup_Sim_Biologicals(
    sim_list = sim_list,
    natmort_input = replicate(1, array(0.3, dim = c(1, 1, n_yrs, n_ages, n_sexes))),
    WAA_input = replicate(1, biol6(waa)),
    WAA_fish_input = replicate(1, waa_fleet),
    WAA_srv_input = replicate(1, waa_fleet),
    MatAA_input = replicate(1, biol6(mat))
  ))
  sim_list <- Setup_Sim_Tagging(sim_list = sim_list, use_conv_fish_tagging = 0)
  sim_list$Movement <- array(1, dim = c(1, 1, 1, n_yrs, 1, n_ages, n_sexes, 1))
  sim_list <- Setup_Sim_Rec(
    sim_list = sim_list,
    R0_input = replicate(1, array(5, dim = c(1, 1, n_yrs))),
    ln_sigmaR = array(log(0.1), dim = c(2, 1, 1)),
    recruitment_opt = "mean_rec",
    init_age_strc = 1,
    # one curve per sex, the new 5-D layout
    ln_InitDevs_input = array(init_devs_true, dim = c(1, 1, n_ages - 1, n_sexes, 1))
  )
  set.seed(77)
  sim_out <- Simulate_Pop_Static(sim_list = sim_list, output_path = NULL)

  # the operating model's initial numbers really are sex-specific
  expect_gt(max(abs(sim_out$ln_InitDevs[1,1,,1,1] - sim_out$ln_InitDevs[1,1,,2,1])), 0.2)

  # slice the simulated environment into estimation-shaped data
  sim_data <- simulation_data_to_SPoRC(sim_env = sim_out, y = n_yrs, sim = 1)

  # Estimation model with the same structure expressed through the options
  input_list <- Setup_Mod_Dim(
    years = 1:n_yrs,
    ages = 1:n_ages,
    lens = NULL,
    n_regions = 1,
    n_sexes = n_sexes,
    n_fish_fleets = 1,
    n_srv_fleets = 1,
    n_pop = 1,
    natal_region = 1,
    verbose = FALSE
  )
  input_list <- Setup_Mod_Rec(
    input_list = input_list,
    do_rec_bias_ramp = 0,
    sigmaR_switch = 1,
    ln_sigmaR = array(log(0.1), c(2, 1, 1)),
    rec_model = "mean_rec",
    sigmaR_spec = "fix",
    init_age_strc = 1,
    equil_init_age_strc = 2,
    InitDevs_sex_spec = "est_all",
    ln_global_R0 = log(5)
  )
  input_list <- Setup_Mod_Biologicals(
    input_list = input_list,
    WAA = sim_data$WAA,
    MatAA = sim_data$MatAA,
    WAA_fish = sim_data$WAA_fish,
    WAA_srv = sim_data$WAA_srv,
    fit_lengths = 0,
    AgeingError = sim_data$AgeingError,
    M_spec = "fix",
    Fixed_natmort = array(0.3, dim = c(1, 1, n_yrs, n_ages, n_sexes))
  )
  input_list <- Setup_Mod_Tagging(input_list = input_list, use_conv_fish_tagging = 0)
  input_list <- Setup_Mod_Movement(
    input_list = input_list,
    use_fixed_movement = 1,
    Fixed_Movement = NA,
    do_recruits_move = 0
  )
  suppressWarnings(input_list <- Setup_Mod_Catch_and_F(
    input_list = input_list,
    ObsCatch = sim_data$ObsCatch,
    UseCatch = sim_data$UseCatch,
    Use_F_pen = 1,
    sigmaC_spec = "fix",
    ln_sigmaC = sim_data$ln_sigmaC,
    ln_sigmaF = array(log(1), dim = c(1, 1, 1)),
    ObsDiscard = sim_data$ObsDiscard,
    UseDiscard = sim_data$UseDiscard,
    sigma_dmr_spec = "fix",
    dmr_mean_spec = "fix",
    ln_sigmaD = sim_data$ln_sigmaD
  ))
  input_list <- Setup_Mod_FishIdx_and_Comps(
    input_list = input_list,
    ObsFishIdx = sim_data$ObsFishIdx,
    ObsFishIdx_SE = sim_data$ObsFishIdx_SE,
    UseFishIdx = array(0, dim = dim(sim_data$UseFishIdx)),
    ObsFishAgeComps = sim_data$ObsFishAgeComps,
    ObsFishLenComps = NULL,
    UseFishAgeComps = sim_data$UseFishAgeComps,
    UseFishLenComps = array(0, dim = dim(sim_data$UseFishAgeComps)),
    ISS_FishAgeComps = sim_data$ISS_FishAgeComps,
    ISS_FishLenComps = NULL,
    fish_idx_type = "biom",
    FishAgeComps_LikeType = "Multinomial",
    FishLenComps_LikeType = "none",
    FishAgeComps_Type = "spltRjntS_Year_1-terminal_Fleet_1",
    FishLenComps_Type = "none_Year_1-terminal_Fleet_1"
  )
  input_list <- Setup_Mod_SrvIdx_and_Comps(
    input_list = input_list,
    ObsSrvIdx = sim_data$ObsSrvIdx,
    ObsSrvIdx_SE = sim_data$ObsSrvIdx_SE,
    UseSrvIdx = sim_data$UseSrvIdx,
    ObsSrvAgeComps = sim_data$ObsSrvAgeComps,
    ObsSrvLenComps = NULL,
    UseSrvAgeComps = sim_data$UseSrvAgeComps,
    UseSrvLenComps = array(0, dim = dim(sim_data$UseSrvAgeComps)),
    ISS_SrvAgeComps = sim_data$ISS_SrvAgeComps,
    ISS_SrvLenComps = NULL,
    srv_idx_type = "biom",
    SrvAgeComps_LikeType = "Multinomial",
    SrvLenComps_LikeType = "none",
    SrvAgeComps_Type = "spltRjntS_Year_1-terminal_Fleet_1",
    SrvLenComps_Type = "none_Year_1-terminal_Fleet_1"
  )
  input_list <- Setup_Mod_Fishsel_and_Q(
    input_list = input_list,
    fish_sel_model = "logist1_Fleet_1_NSelBins_4",
    fish_fixed_sel_pars_spec = "est_all",
    fish_q_spec = "fix",
    fish_sel_sex_offset = "scale",
    use_fixed_ret_sel = 1
  )
  input_list <- Setup_Mod_Srvsel_and_Q(
    input_list = input_list,
    srv_sel_model = "logist1_Fleet_1",
    srv_fixed_sel_pars_spec = "est_all",
    srv_q_spec = "est_all",
    srv_sel_sex_offset = "par"
  )
  input_list <- Setup_Mod_Weighting(
    input_list = input_list,
    Wt_Catch = 1,
    Wt_FishIdx = 0,
    Wt_SrvIdx = 1,
    Wt_Rec = 1,
    Wt_Init_Rec = 1,
    Wt_F = 1,
    Wt_FishAgeComps = array(1, dim = c(1, n_yrs, 1, n_sexes, 1)),
    Wt_SrvAgeComps = array(1, dim = c(1, n_yrs, 1, n_sexes, 1))
  )

  fit <- fit_model(input_list$data, input_list$par, input_list$map, random = NULL, silent = TRUE)
  expect_lt(max(abs(fit$gr(fit$env$last.par.best))), 1e-3)
  sd_rep <- RTMB::sdreport(fit)

  # single-fit recovery of the offset parameters themselves
  pl <- fit$env$parList(fit$env$last.par.best)
  expect_lt(abs(pl$ln_fishsel_sex_scale[1,1,2,1] - fish_scale_true), 0.1)
  expect_lt(max(abs(pl$srv_fixed_sel_pars[1,,1,2,1] - c(srv_off_true[1], srv_off_true[2]))), 0.15)
  expect_lt(max(abs(pl$ln_InitDevs[1,1,,1] - init_devs_true[,1])), 0.25)
  expect_lt(max(abs(pl$ln_InitDevs[1,1,,2] - init_devs_true[,2])), 0.25)
  # both sexes' fitted curves hold the plateau (bin 5 rides bin 4)
  for(s in 1:n_sexes) expect_equal(as.numeric(fit$rep$fish_sel[1,1,1,1,5,s,1]), as.numeric(fit$rep$fish_sel[1,1,1,1,4,s,1]), tolerance = 1e-10)

  # the self-test conditions the operating model on the fit, so the per-sex
  # surfaces and initial curves travel back through the simulation wiring
  set.seed(107)
  res <- simulation_self_test(
    data = fit$data,
    parameters = input_list$par,
    mapping = input_list$map,
    random = NULL,
    rep = fit$rep,
    sd_rep = sd_rep,
    n_sims = n_sims,
    what = c("SSB", "fish_sel", "srv_sel"),
    sim_recruitment = "input"
  )
  for(w in c("SSB", "fish_sel", "srv_sel")) { # res also has an sd_rep slot, so iterate the names
    truth <- as.numeric(fit$rep[[w]])
    est_mat <- matrix(res[[w]], nrow = length(truth), ncol = n_sims)
    keep <- truth > 1e-8 # selectivity arrays have structural zeros outside 1:n_yrs
    re <- sweep(est_mat[keep, , drop = FALSE], 1, truth[keep], "-") / truth[keep]
    expect_lt(abs(stats::median(re, na.rm = TRUE)), 0.05)
  } # end w loop
})
