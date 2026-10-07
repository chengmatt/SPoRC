# Catchability deviations in the self test: conditional keeps the fit's, and joint draws the process the fit
# penalizes fresh, at each replicate's own sigma and only where the fit estimates.

library(SPoRC)
library(testthat)

test_that("conditional keeps the fit's catchability walk and joint draws a fresh one at each replicate's own sigma", {

  rw <- q_rw_fit()

  # conditional: every replicate runs on the fit's deviations
  sl <- q_selftest_capture(rw$fit, rw$sd_rep, n_sims = 2, sim_type = "conditional")
  expect_equal(sl$ln_srv_q_devs[1,,1,1], sl$ln_srv_q_devs[1,,1,2])

  # joint: a random walk at each replicate's own sigma
  sl <- q_selftest_capture(rw$fit, rw$sd_rep, n_sims = 2, sim_type = "joint")
  expect_equal(sl$srv_q_model, 3) # a random walk
  expect_true(all(sl$srv_q_devs_est))
  expect_false(isTRUE(all.equal(sl$sigma_srv_q[1,1,1], sl$sigma_srv_q[1,1,2])))

  env <- Setup_sim_env(sl)
  set.seed(1)
  for(sim in 1:2) draw_sim_q_devs(sim, env)
  devs <- env$ln_srv_q_devs[1,,1,]
  for(sim in 1:2) {
    expect_false(isTRUE(all.equal(devs[,sim], sl$ln_srv_q_devs[1,,1,sim]))) # drawn, not the deviations it was handed
    expect_equal(stats::sd(diff(devs[,sim])), sl$sigma_srv_q[1,1,sim], tolerance = 0.3) # steps at its own sigma
  } # end sim loop

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

test_that("a walk the fit starts diffusely keeps its first estimated year and steps on from it", {

  # the first estimated year, 16, sits at 0.7, a level the replicate's own fit or draw set
  n_yrs <- 30
  fit_devs <- c(rep(0, 15), 0.7, rep(0, 14))
  draw_devs <- function(init_sigma) {
    sl <- Setup_Sim_q_devs(q_cond_sl(n_yrs = n_yrs, devs = fit_devs, n_cond = 0), srv_q_model = "rw", sigma_srv_q = 0.2,
                           srv_q_rw_init_sigma = init_sigma)
    sl$srv_q_devs_est <- array(c(rep(FALSE, 15), rep(TRUE, 15)), dim = c(1, n_yrs, 1))
    env <- Setup_sim_env(sl)
    replicate(3000, { draw_sim_q_devs(1, env); env$ln_srv_q_devs[1,,1,1] })
  }

  set.seed(3)
  diffuse <- draw_devs(5)
  expect_true(all(diffuse[16,] == 0.7)) # the diffuse start keeps its value
  expect_equal(mean(diffuse[17,]), 0.7, tolerance = 0.02) # and the walk steps on from it
  expect_equal(sd(diffuse[17,]), 0.2, tolerance = 0.05)

  # under NA the penalty starts the walk at zero under its own sigma, and so does the draw
  own <- draw_devs(NA)
  expect_equal(mean(own[16,]), 0, tolerance = 0.02)
  expect_equal(sd(own[16,]), 0.2, tolerance = 0.05)

})

test_that("a joint self test recovers the catchability sigma", {

  skip_on_cran()

  # every deviation is drawn fresh from each replicate's walk. before the walk was routed the refits saw flat
  # catchability and put sigma at zero
  rw <- q_rw_fit()
  set.seed(11)
  st <- suppressWarnings(suppressMessages(
    simulation_self_test(data = rw$fit$data, parameters = rw$fit$env$parList(), mapping = rw$fit$mapping,
                         random = "ln_srv_q_devs", rep = rw$fit$rep, sd_rep = rw$sd_rep, obj = rw$fit,
                         n_sims = 4, sim_type = "joint", newton_loops = 1, what = "SSB", what_par = "ln_sigma_srv_q")))
  sigma_hat <- exp(as.numeric(st$ln_sigma_srv_q))
  sigma_true <- exp(as.numeric(st$truth$ln_sigma_srv_q)) # each replicate's own draw

  expect_true(all(sigma_hat > 0.5 * sigma_true))
  expect_equal(mean(sigma_hat / sigma_true), 1, tolerance = 0.2)

})
