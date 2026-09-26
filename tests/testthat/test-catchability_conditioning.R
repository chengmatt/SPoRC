# The operating model takes the fit's block mean catchability and its deviations separately. These cover
# the split, the conditioning years reproducing the fit, and what a projection walks from.

test_that("the split inverts the fit's reported catchability", {

  n_yrs <- 12
  set.seed(4)
  ln_q <- array(log(c(0.05, 0.12)), dim = c(1, 2, 1))
  blocks <- array(rep(c(1, 2), each = n_yrs / 2), dim = c(1, n_yrs, 1))
  devs <- array(stats::rnorm(n_yrs, 0, 0.3), dim = c(1, n_yrs, 1))
  rep_q <- array(exp(ln_q[1, blocks[1,,1], 1]), dim = c(1, n_yrs, 1)) * exp(devs)

  out <- split_reported_q(rep_q, ln_q, blocks, q_type = 0)

  expect_equal(as.numeric(out$q_mean), exp(ln_q[1, blocks[1,,1], 1]), tolerance = 1e-12)
  expect_equal(as.numeric(out$devs), as.numeric(devs), tolerance = 1e-12)
  expect_equal(as.numeric(out$q_mean * exp(out$devs)), as.numeric(rep_q), tolerance = 1e-12)
})

test_that("an analytically solved fleet is its own mean and takes no deviations", {

  n_yrs <- 6
  rep_q <- array(seq(0.01, 0.06, length.out = n_yrs), dim = c(1, n_yrs, 1))
  ln_q <- array(log(999), dim = c(1, 1, 1)) # never read, so a wrong value must not reach the mean
  blocks <- array(1, dim = c(1, n_yrs, 1))

  out <- split_reported_q(rep_q, ln_q, blocks, q_type = 1)

  expect_equal(as.numeric(out$q_mean), as.numeric(rep_q))
  expect_equal(as.numeric(out$devs), rep(0, n_yrs))
})

test_that("conditioning years reproduce the fit and a projection walk starts from the last deviation", {

  n_yrs <- 30
  n_cond <- 20
  set.seed(3)
  fit_devs <- c(stats::rnorm(n_cond, 0, 0.2), rep(0, n_yrs - n_cond))
  sl <- q_cond_sl(n_yrs = n_yrs, devs = fit_devs, n_cond = n_cond)
  env <- Setup_sim_env(Setup_Sim_q_devs(sl, srv_q_model = "rw", sigma_srv_q = 0.05))

  set.seed(8)
  draw_sim_q_devs(1, env)

  # the mean times the fit's own deviation is the catchability the fit reported
  expect_equal(as.numeric(env$srv_q[1, seq_len(n_cond), 1, 1]),
               0.05 * exp(fit_devs[seq_len(n_cond)]), tolerance = 1e-12)
  expect_equal(as.numeric(env$ln_srv_q_devs[1, seq_len(n_cond), 1, 1]), fit_devs[seq_len(n_cond)])

  # and the walk continues from the fit's last deviation rather than restarting at zero
  expect_lt(abs(env$ln_srv_q_devs[1, n_cond + 1, 1, 1] - fit_devs[n_cond]), 0.25)
  expect_gt(max(abs(env$ln_srv_q_devs[1, (n_cond + 1):n_yrs, 1, 1])), 0.01)
})

test_that("a process error is refused when every year is a conditioning year", {

  sl <- q_cond_sl(n_yrs = 20, n_cond = 20)
  for(form in c("iid", "rw", "ar1", "dsem")) {
    expect_error(Setup_Sim_q_devs(sl, srv_q_model = form, sigma_srv_q = 0.2), "would never be drawn")
  } # end form loop

  expect_error(Setup_Sim_q_devs(sl, fish_q_model = "rw", sigma_fish_q = 0.2), "would never be drawn")
  expect_no_error(Setup_Sim_q_devs(q_cond_sl(n_yrs = 30, n_cond = 20), srv_q_model = "rw", sigma_srv_q = 0.2))
})

test_that("deviations on a different grid than the catchability are refused", {

  sl <- q_cond_sl(n_yrs = 20)
  sl$ln_srv_q_devs <- array(0, dim = c(1, 15, 1, 1))
  expect_error(Setup_sim_env(sl), "same region, year, fleet and replicate grid")
})

test_that("a dsem cell reaches catchability on a fleet with no process error of its own", {

  n_yrs <- 20
  env <- Setup_sim_env(q_cond_sl(n_yrs = n_yrs))
  expect_equal(as.numeric(env$srv_q_model), 1) # "none", so nothing of the fleet's own is drawn

  set.seed(6)
  drawn <- stats::rnorm(n_yrs, 0, 0.3)
  env$ln_srv_q_devs[1,,1,1] <- drawn
  env$dsem_drawn$ln_srv_q_devs[1,,1] <- TRUE

  draw_sim_q_devs(1, env)
  expect_equal(as.numeric(env$srv_q[1,,1,1]), 0.05 * exp(drawn), tolerance = 1e-12)

  # a year the dsem did not write keeps the catchability it came with
  env2 <- Setup_sim_env(q_cond_sl(n_yrs = n_yrs))
  env2$ln_srv_q_devs[1,,1,1] <- drawn
  env2$dsem_drawn$ln_srv_q_devs[1, 1:5, 1] <- TRUE

  draw_sim_q_devs(1, env2)
  expect_equal(as.numeric(env2$srv_q[1, 1:5, 1, 1]), 0.05 * exp(drawn[1:5]), tolerance = 1e-12)
  expect_equal(as.numeric(env2$srv_q[1, 6:n_yrs, 1, 1]), rep(0.05, n_yrs - 5))
})

test_that("a self test reproduces the fit's catchability, dsem series and all", {

  fit <- q_dsem_fit()
  sl <- q_selftest_simlist(fit)

  expect_equal(sl$n_cond_yrs, length(fit$data$years))

  # the mean is the block value and the deviations are what the fit put on top of it
  q_blk <- exp(as.numeric(fit$env$parList()$ln_srv_q[1,1,1]))
  expect_equal(as.numeric(sl$srv_q[1,,1,1]), rep(q_blk, length(fit$data$years)), tolerance = 1e-12)
  expect_gt(stats::sd(as.numeric(sl$ln_srv_q_devs[1,,1,1])), 0.05)

  env <- Setup_sim_env(sl)
  expect_true(all(env$dsem_drawn$ln_srv_q_devs)) # the dsem wrote every fitted year

  draw_sim_q_devs(1, env)
  expect_equal(as.numeric(env$srv_q[1,,1,1]), as.numeric(fit$rep$srv_q[1,,1]), tolerance = 1e-12)
})

test_that("a closed loop reads the fitted years back and lets a dsem drive the projection", {

  fit <- q_dsem_fit()
  pars <- fit$env$parList()
  n_fit <- length(fit$data$years)
  n_cl <- 8

  sl <- suppressWarnings(suppressMessages(
    condition_closed_loop_simulations(closed_loop_yrs = n_cl, n_sims = 2, data = fit$data, parameters = pars,
                                      mapping = fit$mapping, sd_rep = q_dsem_sdrep(fit), rep = fit$rep,
                                      random = "ln_srv_q_devs")))
  sl <- suppressWarnings(suppressMessages(Setup_Sim_DSEM(sl, fit$data, pars, rep = fit$rep, condition_on_fit = TRUE)))

  q_blk <- exp(as.numeric(pars$ln_srv_q[1,1,1]))
  expect_equal(as.numeric(sl$srv_q[1,,1,1]), rep(q_blk, n_fit + n_cl), tolerance = 1e-12)

  env <- Setup_sim_env(sl)
  draw_sim_q_devs(1, env)
  draw_sim_q_devs(2, env)
  proj <- (n_fit + 1):(n_fit + n_cl)

  # the fitted years come back to the fit's own catchability, drawn cells and all
  expect_equal(as.numeric(env$srv_q[1, seq_len(n_fit), 1, 1]), as.numeric(fit$rep$srv_q[1,,1]), tolerance = 1e-3)

  # and the projection is the block value times the dsem's own draw, a different path in each replicate
  expect_equal(as.numeric(env$srv_q[1, proj, 1, 1]),
               q_blk * exp(as.numeric(env$ln_srv_q_devs[1, proj, 1, 1])), tolerance = 1e-12)
  expect_gt(max(abs(env$srv_q[1, proj, 1, 1] - env$srv_q[1, proj, 1, 2])), 0.1)
})

test_that("a projected catchability series is read over the conditioning years and refused past them", {

  fit <- q_dsem_projected_fit()
  pars <- fit$env$parList()
  n_fit <- length(fit$data$years)
  expect_true(any(fit$data$dsem_model$project_k)) # the sd is fixed at zero, so nothing is integrated

  # a self test runs the fitted years alone, where the series is read from the fit rather than drawn
  sl <- q_selftest_simlist(fit)
  env <- expect_no_error(Setup_sim_env(sl))
  draw_sim_q_devs(1, env)
  expect_equal(as.numeric(env$srv_q[1,,1,1]), as.numeric(fit$rep$srv_q[1,,1]), tolerance = 1e-12)

  # past them there is no innovation to draw, so the cells would come back NaN rather than a catchability
  sl_cl <- suppressWarnings(suppressMessages(
    condition_closed_loop_simulations(closed_loop_yrs = 8, n_sims = 2, data = fit$data, parameters = pars,
                                      mapping = fit$mapping, sd_rep = RTMB::sdreport(fit), rep = fit$rep,
                                      random = NULL)))
  add_dsem <- function(sl) suppressWarnings(suppressMessages(Setup_Sim_DSEM(sl, fit$data, pars, rep = fit$rep, condition_on_fit = TRUE)))
  expect_error(add_dsem(sl_cl), "sd fixed at zero") # said at the setup call, before any cell is drawn

  # and the refusal is about the years left to draw, not the catchability the operating model was handed
  sl_cl$n_cond_yrs <- sl_cl$n_yrs
  sl_env <- add_dsem(sl_cl)
  expect_silent(check_q_dsem_drawable(sl_env))
  expect_no_error(Setup_sim_env(sl_env))
})
