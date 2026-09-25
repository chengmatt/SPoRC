# A self test of the dsem end to end: the operating model states its own arrows and values, recruitment follows an
# AR1 covariate, and the estimation model declares RecDevs_model = "dsem". With near-exact data and almost no
# process error the effect, R0 and every deviation come back to the third decimal; with real process error the
# effect is recovered within its standard error over replicates.

dsem_selftest_om <- function(arrows, values, cov_obs_sd = c(env = NA), n_sims = 6, n_yrs = 40, seed = 321, mod_var_logscale = FALSE, exact = FALSE) {

  sim_list <- Setup_Sim_Dim(n_sims = n_sims, n_yrs = n_yrs, n_regions = 1, n_ages = 10, n_lens = NULL, n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = 1, n_pop = 1)
  sim_list <- Setup_Sim_Containers(sim_list)
  logistic <- function(slope, infl) 1 / (1 + exp(-slope * ((1:10) - infl)))
  yearly <- function(v, dims) replicate(sim_list$n_sims, array(rep(v, each = n_yrs), dim = dims))
  # exact: catch and indices measured to a tenth of a percent, compositions from a hundred thousand fish
  obs_sd <- if(exact) 0.001 else 0.2
  iss <- if(exact) 1e5 else 100
  sim_list <- Setup_Sim_Fishing(sim_list, fish_sel_input = yearly(logistic(3, 5), c(1, 1, n_yrs, 1, 10, 1, 1)),
                                ln_sigmaC = array(log(if(exact) 0.001 else 0.02), dim = c(1, n_yrs, 1, 1)), ObsFishIdx_SE = array(obs_sd, dim = c(1, n_yrs, 1, 1)),
                                ISS_FishAgeComps = array(iss, dim = c(1, n_yrs, 1, 1, 1, n_sims)))
  sim_list <- Setup_Sim_Survey(sim_list, srv_sel_input = yearly(logistic(1, 3), c(1, 1, n_yrs, 1, 10, 1, 1)),
                               ObsSrvIdx_SE = array(obs_sd, dim = c(1, n_yrs, 1, 1)), ISS_SrvAgeComps = array(iss, dim = c(1, n_yrs, 1, 1, 1, n_sims)))
  sim_list <- suppressWarnings(Setup_Sim_Biologicals(sim_list,
                                                     natmort_input = replicate(sim_list$n_sims, array(0.3, dim = c(1, 1, n_yrs, 10, 1))),
                                                     WAA_input = yearly(5 * logistic(3, 3), c(1, 1, n_yrs, 1, 10, 1)),
                                                     WAA_fish_input = yearly(5 * logistic(3, 3), c(1, 1, n_yrs, 1, 10, 1, 1)),
                                                     WAA_srv_input = yearly(5 * logistic(3, 3), c(1, 1, n_yrs, 1, 10, 1, 1)),
                                                     MatAA_input = yearly(logistic(3, 3), c(1, 1, n_yrs, 1, 10, 1))))
  sim_list <- Setup_Sim_Tagging(sim_list, use_conv_fish_tagging = 0)
  sim_list$Movement <- array(1, dim = c(1, 1, 1, n_yrs, 1, 10, 1, sim_list$n_sims))
  sim_list <- Setup_Sim_Rec(sim_list,
                            R0_input = replicate(sim_list$n_sims, array(5, dim = c(1, 1, n_yrs))),
                            rinit_input = array(2, dim = c(1, 1, sim_list$n_sims)),
                            use_rinit = 1,
                            ln_sigmaR = array(log(1), dim = c(2, 1, 1)),
                            recruitment_opt = "mean_rec",
                            init_age_strc = 1,
                            ln_InitDevs_input = array(0, dim = c(1, 1, 9, 1, sim_list$n_sims))) # the population starts in equilibrium

  # the dsem written from scratch: these arrows and values are the truth the refit is judged against
  sim_list <- Setup_Sim_DSEM(sim_list, dsem_arrows = arrows, dsem_values = values, dsem_cov_obs_sd = cov_obs_sd, mod_var_logscale = mod_var_logscale)
  set.seed(seed)
  Simulate_Pop_Static(sim_list = sim_list, output_path = NULL)

}

# the simple self test's estimation model, with the recruitment deviations handed to the dsem. the covariate mean is
# fixed at the truth (zero), the convention the operating model drew under, so R0 means the same thing on both sides
dsem_selftest_em <- function(sim_obj, sim, arrows, cov_family = c(env = "fixed"), hold = NULL, mod_var_logscale = FALSE, init_devs = TRUE, dsem = TRUE) {

  sim_data <- simulation_data_to_SPoRC(sim_env = sim_obj, y = sim_obj$n_years, sim = sim)
  il <- Setup_Mod_Dim(years = 1:sim_obj$n_years, ages = 1:sim_obj$n_ages, lens = sim_obj$n_lens, n_regions = 1, n_sexes = 1,
                      n_fish_fleets = 1, n_srv_fleets = 1, n_pop = 1, natal_region = 1, verbose = FALSE)
  # with the dsem, the declaration: the deviations' density is the arrows and sigmaR is fixed. without it, plain
  # independent deviations with the one sigmaR the years read estimated, which is the marginal sd the covariate explains part of
  il <- if(dsem) Setup_Mod_Rec(il, do_rec_bias_ramp = 0, sigmaR_switch = 1, ln_sigmaR = array(log(1), c(2, 1, 1)), rec_model = "mean_rec", use_rinit = 1,
                               init_age_strc = 1, equil_init_age_strc = 2, ln_global_R0 = log(5), ln_rinit = log(2), RecDevs_model = "dsem")
        else Setup_Mod_Rec(il, do_rec_bias_ramp = 0, sigmaR_switch = 1, ln_sigmaR = array(log(1), c(2, 1, 1)), rec_model = "mean_rec", use_rinit = 1,
                           init_age_strc = 1, equil_init_age_strc = 2, ln_global_R0 = log(5), ln_rinit = log(2), RecDevs_model = "iid", sigmaR_spec = "fix_early_est_late")
  il <- Setup_Mod_Biologicals(il, WAA = sim_data$WAA, MatAA = sim_data$MatAA, WAA_fish = sim_data$WAA_fish, WAA_srv = sim_data$WAA_srv, fit_lengths = 0,
                              AgeingError = sim_data$AgeingError, M_spec = "fix", Fixed_natmort = array(0.3, dim = c(1, 1, sim_obj$n_years, 10, 1)))
  il <- Setup_Mod_Tagging(il, use_conv_fish_tagging = 0)
  il <- Setup_Mod_Movement(il, use_fixed_movement = 1, Fixed_Movement = NA, do_recruits_move = 0)
  il <- suppressWarnings(Setup_Mod_Catch_and_F(il, ObsCatch = sim_data$ObsCatch, UseCatch = sim_data$UseCatch, Use_F_pen = 1, sigmaC_spec = "fix",
                                               ln_sigmaC = sim_data$ln_sigmaC, ln_sigmaF = array(log(1), dim = c(1, 1, 1))))
  il <- Setup_Mod_FishIdx_and_Comps(il, ObsFishIdx = sim_data$ObsFishIdx, ObsFishIdx_SE = sim_data$ObsFishIdx_SE, UseFishIdx = sim_data$UseFishIdx,
                                    ObsFishAgeComps = sim_data$ObsFishAgeComps, ObsFishLenComps = sim_data$ObsFishLenComps, UseFishAgeComps = sim_data$UseFishAgeComps,
                                    UseFishLenComps = sim_data$UseFishLenComps, ISS_FishAgeComps = sim_data$ISS_FishAgeComps, ISS_FishLenComps = sim_data$ISS_FishLenComps,
                                    fish_idx_type = "biom", FishAgeComps_LikeType = "Multinomial", FishLenComps_LikeType = "none",
                                    FishAgeComps_Type = "agg_Year_1-terminal_Fleet_1", FishLenComps_Type = "none_Year_1-terminal_Fleet_1")
  il <- Setup_Mod_SrvIdx_and_Comps(il, ObsSrvIdx = sim_data$ObsSrvIdx, ObsSrvIdx_SE = sim_data$ObsSrvIdx_SE, UseSrvIdx = sim_data$UseSrvIdx,
                                   ObsSrvAgeComps = sim_data$ObsSrvAgeComps, ObsSrvLenComps = sim_data$ObsSrvLenComps, UseSrvAgeComps = sim_data$UseSrvAgeComps,
                                   UseSrvLenComps = sim_data$UseSrvLenComps, ISS_SrvAgeComps = sim_data$ISS_SrvAgeComps, ISS_SrvLenComps = sim_data$ISS_SrvLenComps,
                                   srv_idx_type = "biom", SrvAgeComps_LikeType = "Multinomial", SrvLenComps_LikeType = "none",
                                   SrvAgeComps_Type = "agg_Year_1-terminal_Fleet_1", SrvLenComps_Type = "none_Year_1-terminal_Fleet_1")
  il <- Setup_Mod_Fishsel_and_Q(il, fish_sel_model = "logist2_Fleet_1", fish_fixed_sel_pars_spec = "est_all", fish_q_spec = "est_all")
  il <- Setup_Mod_Srvsel_and_Q(il, srv_sel_model = "logist2_Fleet_1", srv_fixed_sel_pars_spec = "est_all", srv_q_spec = "est_all")
  il <- Setup_Mod_Weighting(il, Wt_Catch = 1, Wt_FishIdx = 1, Wt_SrvIdx = 1, Wt_Rec = 1, Wt_F = 1, Wt_Tagging = 0,
                            Wt_FishAgeComps = array(1, dim = c(1, sim_obj$n_years, 1, 1, 1)), Wt_FishLenComps = array(1, dim = c(1, sim_obj$n_years, 1, 1, 1)),
                            Wt_SrvAgeComps = array(1, dim = c(1, sim_obj$n_years, 1, 1, 1)), Wt_SrvLenComps = array(1, dim = c(1, sim_obj$n_years, 1, 1, 1)))
  # the operating model started in equilibrium, so the initial deviations may be taken as known
  if(!init_devs) {
    il$map$ln_InitDevs <- factor(rep(NA, length(il$par$ln_InitDevs)))
    il$par$ln_InitDevs[] <- 0
  }

  if(!dsem) return(il)

  # the covariate this replicate observed, with its measurement error, and the same arrows; any other series the
  # arrows name that is not a deviation series is a latent column, all NA
  env <- data.frame(year = 1:sim_obj$n_years, env = sim_obj$dsem_cov_obs_sim[,1,sim])
  for(nm in setdiff(names(cov_family), "env")) env[[nm]] <- NA_real_
  Setup_Mod_DSEM(il, hold_arrows(arrows, hold), env, dsem_family = cov_family, dsem_mu_spec = c(env = 0), mod_var_logscale = mod_var_logscale)

}

test_that("with near-exact data and almost no process error the refit returns the truth", {

  arrows <- c("env -> rec, 0, b_env", "env -> env, 1, rho_env", "env <-> env, 0, sd_env", "rec <-> rec, 0, sd_rec")
  truth <- c(b_env = 0.6, rho_env = 0.5, sd_env = 1, sd_rec = 0.02) # recruitment is the covariate's effect and almost nothing else
  om <- dsem_selftest_om(arrows, truth, n_sims = 3, exact = TRUE)

  # the truth went in: recruitment deviations are the drawn grid column, the covariate is the state itself, no initial deviations
  expect_equal(om$ln_RecDevs[1,1,,1], om$dsem_x_sim[,2,1], tolerance = 1e-12)
  expect_equal(om$dsem_cov_obs_sim[,1,1], om$dsem_x_sim[,1,1], tolerance = 1e-12)
  expect_true(all(om$dsem_drawn$ln_RecDevs))
  expect_true(all(om$ln_InitDevs == 0))

  for(i in 1:om$n_sims) {
    il <- dsem_selftest_em(om, i, arrows, hold = c(sd_rec = 0.02, sd_env = 1, rho_env = 0.5), init_devs = FALSE)
    expect_equal(il$data$dsem_declared, "rec")
    expect_true(all(is.na(il$map$ln_sigmaR))) # the declaration fixed sigmaR
    fit <- fit_model(il$data, il$par, il$map, random = c("ln_RecDevs", "dsem_x"), newton_loops = 3, silent = TRUE)
    pl <- fit$env$parList()
    expect_lt(max(abs(fit$gr(fit$env$last.par.best[-fit$env$random]))), 1e-4)
    expect_equal(pl$dsem_beta[1], 0.6, tolerance = 0.02) # the effect
    expect_equal(as.numeric(exp(pl$ln_global_R0)), 5, tolerance = 0.01) # R0 means the same thing on both sides
    expect_equal(as.numeric(exp(pl$ln_rinit)), 2, tolerance = 0.01)
    expect_lt(max(abs(pl$ln_RecDevs[1,1,] - om$ln_RecDevs[1,1,,i])), 0.03) # every deviation
    expect_lt(max(abs(fit$rep$SSB[1,1,] / om$SSB[1,1,,i] - 1)), 0.005) # and the biomass
  } # end i loop

})

test_that("with real process error the effect is recovered within its standard error over replicates", {

  arrows <- c("env -> rec, 0, b_env", "env -> env, 1, rho_env", "env <-> env, 0, sd_env", "rec <-> rec, 0, sd_rec")
  truth <- c(b_env = 0.6, rho_env = 0.5, sd_env = 1, sd_rec = 0.5)
  om <- dsem_selftest_om(arrows, truth, n_sims = 6)

  # the covariate's own process and the innovation sd kept at the truth, so the effect is what the refit has to find
  est <- matrix(NA_real_, om$n_sims, 2, dimnames = list(NULL, c("b_env", "se_b")))
  for(i in 1:om$n_sims) {
    il <- dsem_selftest_em(om, i, arrows, hold = c(sd_rec = 0.5, sd_env = 1, rho_env = 0.5))
    fit <- fit_model(il$data, il$par, il$map, random = c("ln_RecDevs", "dsem_x"), newton_loops = 1, silent = TRUE)
    s <- summary(RTMB::sdreport(fit), "fixed")
    est[i, ] <- c(fit$env$parList()$dsem_beta[match("b_env", il$data$dsem_model$beta_names)], s[rownames(s) == "dsem_beta", 2])
  } # end i loop

  # forty years of deviations with an sd of 0.5 identify a slope to about 0.08, and that is what comes back
  expect_lt(abs(mean(est[, "b_env"]) - 0.6), 2 * sd(est[, "b_env"]) / sqrt(om$n_sims) + 0.05)
  expect_true(all(is.finite(est[, "se_b"])))
  expect_lt(abs(median(est[, "se_b"]) / sd(est[, "b_env"]) - 1), 0.75)

})

test_that("the covariate takes variance out of recruitment: the dsem's innovation sd sits below the plain sigmaR", {

  # the same replicates fit twice: with the arrows (innovation sd estimated) and as plain independent deviations
  # (sigmaR estimated). the marginal sd the plain fit sees is sqrt(b^2 var(env) + sd_rec^2), about 0.85 here
  arrows <- c("env -> rec, 0, b_env", "env -> env, 1, rho_env", "env <-> env, 0, sd_env", "rec <-> rec, 0, sd_rec")
  truth <- c(b_env = 0.6, rho_env = 0.5, sd_env = 1, sd_rec = 0.5)
  om <- dsem_selftest_om(arrows, truth, n_sims = 3, exact = TRUE, seed = 11)

  est <- matrix(NA_real_, om$n_sims, 3, dimnames = list(NULL, c("sd_rec", "sigmaR", "R0_iid")))
  for(i in 1:om$n_sims) {
    il <- dsem_selftest_em(om, i, arrows, hold = c(sd_env = 1, rho_env = 0.5), init_devs = FALSE)
    fit <- fit_model(il$data, il$par, il$map, random = c("ln_RecDevs", "dsem_x"), newton_loops = 1, silent = TRUE)
    il_iid <- dsem_selftest_em(om, i, arrows, init_devs = FALSE, dsem = FALSE)
    fit_iid <- fit_model(il_iid$data, il_iid$par, il_iid$map, random = "ln_RecDevs", newton_loops = 1, silent = TRUE)
    est[i, ] <- c(exp(fit$env$parList()$ln_dsem_sd[match("sd_rec", il$data$dsem_model$ln_sd_names)]), exp(fit_iid$env$parList()$ln_sigmaR[2,1,1]), exp(fit_iid$env$parList()$ln_global_R0))
  } # end i loop

  expect_true(all(est[, "sd_rec"] < est[, "sigmaR"]))
  expect_equal(mean(est[, "sd_rec"]), 0.5, tolerance = 0.25)
  expect_equal(mean(est[, "sigmaR"]), sqrt(0.36 * 4 / 3 + 0.25), tolerance = 0.25)
  # both fits run the full lognormal correction. the dsem's is taken given the covariate, so its R0 is mean
  # recruitment at the covariate's mean, 5; the plain fit's is over the covariate too, 5 exp(0.36 var(env) / 2)
  expect_equal(mean(est[, "R0_iid"]), 5 * exp(0.36 * 4 / 3 / 2), tolerance = 0.3)

})

test_that("an effect that changes by year is written and refitted as a moderated arrow", {

  # the effect is its own AR1 series, b_env, and the arrow env -> rec reads it each year
  arrows <- c("b_env -> b_env, 1, rho_b", "b_env <-> b_env, 0, sd_b", "env -> rec, 0, b_env", "env -> env, 1, rho_env", "env <-> env, 0, sd_env", "rec <-> rec, 0, sd_rec")
  truth <- c(rho_b = 0.8, sd_b = 0.15, rho_env = 0.5, sd_env = 1, sd_rec = 0.5)
  om <- dsem_selftest_om(arrows, truth, cov_obs_sd = c(env = 0.2, b_env = NA), n_sims = 2)

  # the drawn coefficient series is what the recruitment deviations were built with
  x <- om$dsem_x_sim[,,1] # columns env, b_env, rec
  eps <- x[-1, 3] - x[-1, 2] * x[-1, 1] # the innovation of rec is its value less the moderated effect
  expect_equal(sd(eps), 0.5, tolerance = 0.25)
  expect_gt(cor(x[-1, 2], x[-nrow(x), 2]), 0) # the effect series is autocorrelated, loosely on forty years

  # the refit takes the same arrows, with the effect series latent (all NA), and evaluates
  il <- dsem_selftest_em(om, 1, arrows, cov_family = c(env = "normal", b_env = "fixed"), hold = c(sd_rec = 0.5, sd_env = 1, sd_b = 0.15, rho_b = 0.8))
  expect_true(any(il$data$dsem_model$arrows$mod_idx > 0))
  obj <- fit_model(il$data, il$par, il$map, random = c("ln_RecDevs", "dsem_x"), do_optim = FALSE, silent = TRUE)
  expect_true(is.finite(obj$fn(obj$par)))

})
