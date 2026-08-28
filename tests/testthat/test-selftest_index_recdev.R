# Self-test of the recruitment deviation index, end to end through the
# simulation module. The operating model carries two survey fleets: an ordinary
# biomass survey, and one whose observation IS the year class strength, drawn as
# q times the deviation with normal error. The estimation model fits the second
# through srv_idx_type = "recdev", so the option is exercised on both sides
# rather than only in the estimation model.
#
# The two sides line up because the operating model draws its deviation about
# zero and carries the bias correction inside recruitment, while the estimation
# model carries the deviation and centers its penalty; the anomaly the index
# reads is the same quantity in both.

library(SPoRC)
library(testthat)
library(RTMB)

test_that("a recruitment deviation index is simulated and recovered", {

  n_yrs <- 40; n_ages <- 8; n_srv <- 2
  sigmaR <- 0.5
  q_rec_true <- 0.8
  waa <- 5 / (1 + exp(-1.2 * ((1:n_ages) - 3)))
  mat <- 1 / (1 + exp(-1.5 * ((1:n_ages) - 3)))
  f_ramp <- c(seq(0.03, 0.35, length.out = 25), seq(0.35, 0.12, length.out = 15))
  curve <- function(slope, infl, n_fleets) {
    array(rep(1 / (1 + exp(-slope * ((1:n_ages) - infl))), each = n_yrs),
          dim = c(1, 1, n_yrs, 1, n_ages, 1, n_fleets))
  }

  sim_list <- Setup_Sim_Dim(n_sims = 1, n_yrs = n_yrs, n_regions = 1, n_ages = n_ages,
                            n_lens = NULL, n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = n_srv, n_pop = 1)
  sim_list <- Setup_Sim_Containers(sim_list)
  sim_list <- Setup_Sim_Fishing(
    sim_list = sim_list,
    fish_sel_input = replicate(1, curve(3, 2, 1)),
    # everything caught is retained, so the estimation model's default full
    # retention is the same statement
    ret_sel_input = replicate(1, array(1, dim = c(1, 1, n_yrs, 1, n_ages, 1, 1))),
    dmr_input = array(0, dim = c(1, n_yrs, 1, 1, 1)),
    Fmort_input = array(f_ramp, dim = c(1, n_yrs, 1, 1, 1)),
    ISS_FishAgeComps = array(300, dim = c(1, n_yrs, 1, 1, 1, 1))
  )

  # fleet 1 is an ordinary biomass survey; fleet 2 observes the deviations
  srv_q <- array(1, dim = c(1, n_yrs, n_srv, 1))
  srv_q[1, , 2, 1] <- q_rec_true
  sim_list <- Setup_Sim_Survey(
    sim_list = sim_list,
    srv_sel_input = replicate(1, curve(1, 3, n_srv)),
    ObsSrvIdx_SE = array(rep(c(0.15, 0.30), each = n_yrs), dim = c(1, n_yrs, 1, n_srv)),
    srv_q_input = srv_q,
    srv_idx_type = c("biom", "recdev"),
    SrvIdx_LikeType = c("lognormal", "normal"),
    ISS_SrvAgeComps = array(300, dim = c(1, n_yrs, 1, 1, n_srv, 1))
  )
  biol <- function(val) array(rep(val, each = n_yrs), dim = c(1, 1, n_yrs, 1, n_ages, 1))
  suppressWarnings(sim_list <- Setup_Sim_Biologicals(
    sim_list = sim_list,
    natmort_input = replicate(1, array(0.25, dim = c(1, 1, n_yrs, n_ages, 1))),
    WAA_input = replicate(1, biol(waa)),
    WAA_fish_input = replicate(1, array(rep(waa, each = n_yrs), dim = c(1, 1, n_yrs, 1, n_ages, 1, 1))),
    WAA_srv_input = replicate(1, array(rep(waa, each = n_yrs), dim = c(1, 1, n_yrs, 1, n_ages, 1, n_srv))),
    MatAA_input = replicate(1, biol(mat))
  ))
  sim_list <- Setup_Sim_Tagging(sim_list = sim_list, use_conv_fish_tagging = 0)
  sim_list$Movement <- array(1, dim = c(1, 1, 1, n_yrs, 1, n_ages, 1, 1))
  sim_list <- Setup_Sim_Rec(
    sim_list = sim_list, recruitment_opt = "mean_rec", do_recruits_move = 0,
    R0_input = replicate(1, array(30, dim = c(1, 1, n_yrs))),
    ln_sigmaR = array(log(sigmaR), dim = c(2, 1, 1)),
    init_age_strc = 1
  )

  set.seed(1234)
  om <- Simulate_Pop_Static(sim_list = sim_list, output_path = NULL)

  # the operating model's index IS the deviation times catchability
  expect_equal(as.vector(om$TrueSrvIdx[1, , 1, 2, 1]),
               q_rec_true * as.vector(om$ln_RecDevs[1, 1, , 1]), tolerance = 1e-10)
  # and it is not the biomass survey, which reads the population
  expect_gt(stats::sd(as.vector(om$TrueSrvIdx[1, , 1, 1, 1])), 0)
  expect_true(any(om$TrueSrvIdx[1, , 1, 2, 1] < 0)) # deviations are signed

  sim_data <- simulation_data_to_SPoRC(sim_env = om, y = n_yrs, sim = 1)

  build <- function(dat) {
    il <- Setup_Mod_Dim(years = 1:n_yrs, ages = 1:n_ages, lens = NA, n_regions = 1, n_sexes = 1,
                        n_fish_fleets = 1, n_srv_fleets = n_srv, n_seas = 1, n_pop = 1,
                        natal_region = 1, verbose = FALSE)
    il <- Setup_Mod_Rec(input_list = il, rec_model = "mean_rec", rec_lag = 1, SR_ref_yr = 1,
                        do_rec_bias_ramp = 0, sigmaR_switch = 1, sigmaR_spec = "fix",
                        ln_sigmaR = array(log(sigmaR), dim = c(2, 1, 1)),
                        # the operating model draws initial age deviations, so the
                        # estimation model estimates them rather than starting from
                        # a bare equilibrium
                        init_age_strc = 1, equil_init_age_strc = 2,
                        RecDevs_pen_center = "fixed",
                        ln_global_R0 = log(30), t_spawn = 0)
    il <- Setup_Mod_Biologicals(input_list = il, WAA = dat$WAA, WAA_fish = dat$WAA_fish,
                                WAA_srv = dat$WAA_srv, MatAA = dat$MatAA, fit_lengths = 0,
                                M_spec = "fix", Fixed_natmort = array(0.25, dim = c(1, 1, n_yrs, n_ages, 1)),
                                AgeingError = dat$AgeingError,
                                addtocomp = 1e-3, comp_const_obs = 1, addtosrvidx = 0, addtofishidx = 0)
    il <- Setup_Mod_Movement(input_list = il, use_fixed_movement = 1, Fixed_Movement = NA, do_recruits_move = 0)
    il <- Setup_Mod_Tagging(input_list = il, use_conv_fish_tagging = 0)
    il <- Setup_Mod_Catch_and_F(input_list = il, ObsCatch = dat$ObsCatch, UseCatch = dat$UseCatch,
                                Use_F_pen = 0, sigmaC_spec = "fix",
                                ln_sigmaC = array(log(0.05), dim = c(1, n_yrs, 1, 1)))
    il <- Setup_Mod_FishIdx_and_Comps(
      input_list = il,
      ObsFishIdx = array(NA_real_, dim = c(1, n_yrs, 1, 1)), ObsFishIdx_SE = array(NA_real_, dim = c(1, n_yrs, 1, 1)),
      UseFishIdx = array(0, dim = c(1, n_yrs, 1, 1)), fish_idx_type = "none", FishIdx_LikeType = "lognormal",
      ObsFishAgeComps = dat$ObsFishAgeComps, UseFishAgeComps = dat$UseFishAgeComps,
      ISS_FishAgeComps = dat$ISS_FishAgeComps,
      ObsFishLenComps = array(NA_real_, dim = c(1, n_yrs, 1, 1, 1, 1)),
      UseFishLenComps = array(0, dim = c(1, n_yrs, 1, 1)), ISS_FishLenComps = array(0, dim = c(1, n_yrs, 1, 1, 1)),
      FishAgeComps_LikeType = "Multinomial", FishLenComps_LikeType = "none",
      FishAgeComps_Type = "agg_Year_1-terminal_Fleet_1", FishLenComps_Type = "none_Year_1-terminal_Fleet_1")
    il <- Setup_Mod_SrvIdx_and_Comps(
      input_list = il,
      ObsSrvIdx = dat$ObsSrvIdx, ObsSrvIdx_SE = dat$ObsSrvIdx_SE, UseSrvIdx = dat$UseSrvIdx,
      srv_idx_type = c("biom", "recdev"), SrvIdx_LikeType = c("lognormal", "normal"),
      ObsSrvAgeComps = dat$ObsSrvAgeComps, UseSrvAgeComps = dat$UseSrvAgeComps,
      ISS_SrvAgeComps = dat$ISS_SrvAgeComps,
      ObsSrvLenComps = array(NA_real_, dim = c(1, n_yrs, 1, 1, 1, n_srv)),
      UseSrvLenComps = array(0, dim = c(1, n_yrs, 1, n_srv)), ISS_SrvLenComps = array(0, dim = c(1, n_yrs, 1, 1, n_srv)),
      SrvAgeComps_LikeType = rep("Multinomial", n_srv), SrvLenComps_LikeType = rep("none", n_srv),
      SrvAgeComps_Type = paste0("agg_Year_1-terminal_Fleet_", 1:n_srv),
      SrvLenComps_Type = paste0("none_Year_1-terminal_Fleet_", 1:n_srv),
      t_srv = array(1, dim = c(1, 1, n_srv)))
    il <- Setup_Mod_Fishsel_and_Q(input_list = il, fish_sel_model = "logist1_Fleet_1",
                                  cont_tv_fish_sel = "none_Fleet_1", fish_sel_blocks = "none_Fleet_1",
                                  fish_q_blocks = "none_Fleet_1", fish_fixed_sel_pars_spec = "est_all",
                                  fish_q_spec = "fix")
    il <- Setup_Mod_Srvsel_and_Q(input_list = il, srv_sel_model = paste0("logist1_Fleet_", 1:n_srv),
                                 cont_tv_srv_sel = paste0("none_Fleet_", 1:n_srv),
                                 srv_sel_blocks = paste0("none_Fleet_", 1:n_srv),
                                 srv_q_blocks = paste0("none_Fleet_", 1:n_srv),
                                 srv_fixed_sel_pars_spec = c("est_all", "fix"),
                                 srv_q_spec = rep("est_all", n_srv),
                                 t_srv = array(1, dim = c(1, 1, n_srv)))
    il <- Setup_Mod_Weighting(input_list = il, Wt_Catch = 1, Wt_FishIdx = 0, Wt_SrvIdx = 1,
                              Wt_Rec = 1, Wt_F = 1, Wt_Tagging = 0,
                              Wt_FishAgeComps = array(1, dim = c(1, n_yrs, 1, 1, 1)),
                              Wt_FishLenComps = array(1, dim = c(1, n_yrs, 1, 1, 1)),
                              Wt_SrvAgeComps = array(1, dim = c(1, n_yrs, 1, 1, n_srv)),
                              Wt_SrvLenComps = array(1, dim = c(1, n_yrs, 1, 1, n_srv)))
    il
  }

  il <- build(sim_data)
  il$par$ln_srv_q[1, 1, 2] <- log(q_rec_true)
  fit <- fit_model(il$data, il$par, il$map, do_optim = TRUE, newton_loops = 3, silent = TRUE)
  pl <- fit$env$parList(fit$optim$par)

  # the estimation model's anomaly is the operating model's own deviation
  expect_equal(as.vector(fit$rep$RecDev_anom[1, 1, ]),
               as.vector(pl$ln_RecDevs[1, 1, ]) + 0.5 * sigmaR^2, tolerance = 1e-10)

  # the deviation index's fleet predicts that anomaly times its catchability, and
  # reads nothing from the population
  expect_equal(as.vector(fit$rep$PredSrvIdx[1, 1, , 1, 2]),
               exp(pl$ln_srv_q[1, 1, 2]) * as.vector(fit$rep$RecDev_anom[1, 1, ]), tolerance = 1e-10)

  # the deviations the index carries come back
  dev_true <- as.vector(om$ln_RecDevs[1, 1, , 1])
  dev_hat <- as.vector(fit$rep$RecDev_anom[1, 1, ])
  expect_gt(stats::cor(dev_true, dev_hat), 0.9)
  q_hat <- exp(pl$ln_srv_q[1, 1, 2])
  expect_lt(abs(q_hat - q_rec_true) / q_rec_true, 0.25)

  # and the fit recovers itself through the simulation, which is what puts the
  # operating model's own recdev path under test: simulation_self_test hands
  # srv_idx_type back to the simulator, so every replicate redraws this fleet
  # as a deviation index rather than as an abundance
  expect_lt(max(abs(fit$gr(fit$env$last.par.best))), 1e-3)
  sd_rep <- RTMB::sdreport(fit)
  set.seed(202)
  n_sims <- 20
  res <- simulation_self_test(data = fit$data, parameters = il$par, mapping = il$map,
                              random = NULL, rep = fit$rep, sd_rep = sd_rep,
                              n_sims = n_sims, what = c("SSB", "srv_q"),
                              sim_recruitment = "input")
  med_re <- function(w) {
    truth <- as.numeric(fit$rep[[w]])
    est <- matrix(res[[which(c("SSB", "srv_q") == w)]], nrow = length(truth), ncol = n_sims)
    stats::median(sweep(est, 1, truth, "-") / matrix(truth, length(truth), n_sims), na.rm = TRUE)
  }
  expect_lt(abs(med_re("SSB")), 0.05)
  expect_lt(abs(med_re("srv_q")), 0.05)
})
