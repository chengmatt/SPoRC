# Catchability deviations in the self test: past the conditioning years the operating model draws the
# process the fit penalizes, at each replicate's own sigma under joint, and only where the fit estimates.

library(SPoRC)
library(testthat)

test_that("past the conditioning years the self test draws the fit's catchability walk", {

  rw <- q_rw_fit()
  sigma_fit <- exp(as.numeric(rw$fit$env$parList(par = rw$fit$env$last.par.best)$ln_sigma_srv_q)) # the estimate the self test reads

  # conditional: every replicate walks at the fitted sigma
  sl <- q_selftest_capture(rw$fit, rw$sd_rep, n_sims = 2, n_cond_yrs = 30, sim_type = "conditional")
  expect_equal(sl$srv_q_model, 3) # a random walk
  expect_equal(as.numeric(sl$sigma_srv_q), rep(sigma_fit, 2))
  expect_true(all(sl$srv_q_devs_est))

  # joint: each replicate has its own sigma and its own deviations over the conditioning years
  sl <- q_selftest_capture(rw$fit, rw$sd_rep, n_sims = 2, n_cond_yrs = 30, sim_type = "joint")
  expect_false(isTRUE(all.equal(sl$sigma_srv_q[1,1,1], sl$sigma_srv_q[1,1,2])))
  expect_false(isTRUE(all.equal(sl$ln_srv_q_devs[1,1:30,1,1], sl$ln_srv_q_devs[1,1:30,1,2])))

  env <- Setup_sim_env(sl)
  set.seed(1)
  for(sim in 1:2) draw_sim_q_devs(sim, env)
  devs <- env$ln_srv_q_devs[1,,1,]
  expect_equal(devs[1:30,], sl$ln_srv_q_devs[1,1:30,1,]) # the conditioning years are each replicate's own
  for(sim in 1:2) expect_equal(stats::sd(diff(devs[30:60,sim])), sl$sigma_srv_q[1,1,sim], tolerance = 0.3) # steps at its own sigma

})

test_that("a deviation the fit keeps fixed keeps its value, and a walk starts at the first estimated year", {

  # the first 15 years and years 26 to 30 are left out of the fit's deviations, the first 15 fixed at 2
  n_yrs <- 40
  set.seed(5)
  fit_devs <- c(rep(2, 15), stats::rnorm(10, 0, 0.2), rep(0, 5), stats::rnorm(10, 0, 0.2))
  sl <- Setup_Sim_q_devs(q_cond_sl(n_yrs = n_yrs, devs = fit_devs, n_cond = 0), srv_q_model = "rw", sigma_srv_q = 0.01)
  sl$srv_q_devs_est <- array(c(rep(FALSE, 15), rep(TRUE, 10), rep(FALSE, 5), rep(TRUE, 10)), dim = c(1, n_yrs, 1))
  env <- Setup_sim_env(sl)
  set.seed(2)
  draw_sim_q_devs(1, env)
  devs <- env$ln_srv_q_devs[1,,1,1]

  expect_equal(devs[c(1:15, 26:30)], fit_devs[c(1:15, 26:30)]) # the fixed years keep the fit's values
  expect_lt(abs(devs[16]), 0.05) # the walk starts at zero, as the penalty's first estimated year does, not from the fixed 2
  expect_lt(abs(devs[31]), 0.05) # and steps on from the fixed zero in year 30

})

test_that("a self test drawn without conditioning recovers the catchability sigma", {

  skip_on_cran()

  # with no conditioning years every deviation is drawn from the fitted walk. before the walk was routed
  # the refits saw flat catchability and put sigma at zero
  rw <- q_rw_fit()
  sigma_fit <- exp(as.numeric(rw$fit$env$parList(par = rw$fit$env$last.par.best)$ln_sigma_srv_q)) # the estimate the self test reads
  set.seed(11)
  st <- suppressWarnings(suppressMessages(
    simulation_self_test(data = rw$fit$data, parameters = rw$fit$env$parList(), mapping = rw$fit$mapping,
                         random = "ln_srv_q_devs", rep = rw$fit$rep, sd_rep = rw$sd_rep, obj = rw$fit,
                         n_sims = 4, n_cond_yrs = 0, newton_loops = 1, what = "SSB", what_par = "ln_sigma_srv_q")))
  sigma_hat <- exp(as.numeric(st$ln_sigma_srv_q))

  expect_true(all(sigma_hat > 0.5 * sigma_fit))
  expect_equal(mean(sigma_hat), sigma_fit, tolerance = 0.2)

})
