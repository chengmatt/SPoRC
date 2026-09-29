# Under sim_type = "joint" every replicate's operating model has to run at that replicate's own
# draw. Checks the draws differ, that they reach the simulation inputs, and that a parameter the
# data pin down returns at its own draw rather than at the fit.

library(SPoRC)
library(testthat)

# near-exact data, so a refit recovers whatever generated it and tracking is visible
routing_fit <- function(seed = 41) {
  om <- selftest_make_om(idx_se_om = 0.005, iss_om = 2e5, seed = seed)
  sim_data <- simulation_data_to_SPoRC(sim_env = om, y = selftest_cfg$n_yrs, sim = 1)
  inp <- selftest_build_input(sim_data)
  fit <- fit_model(inp$data, inp$par, inp$map, random = NULL, silent = TRUE)
  list(fit = fit, inp = inp, sd_rep = RTMB::sdreport(fit, getJointPrecision = TRUE))
}

# the sim_list Simulate_Pop_Static is handed, which is what the operating model runs on
capture_sim_list <- function(f) {
  on.exit(suppressMessages(untrace("Simulate_Pop_Static", where = asNamespace("SPoRC"))), add = TRUE)
  suppressMessages(trace("Simulate_Pop_Static", where = asNamespace("SPoRC"),
                         tracer = quote(assign("..sim_list_seen", sim_list, envir = globalenv())),
                         print = FALSE))
  f()
  get("..sim_list_seen", envir = globalenv())
}

# largest spread across replicates of any single cell; zero means every replicate got the same value
rep_spread <- function(a) {
  d <- dim(a)
  if(is.null(d)) return(0)
  m <- matrix(a, ncol = d[length(d)])
  max(apply(m, 1, function(v) diff(range(v))))
}

run_self_test <- function(fx, type, n_sims = 4, what_par = NULL, seed = 2024) {
  set.seed(seed)
  simulation_self_test(fx$fit$data, fx$inp$par, fx$inp$map, random = NULL,
                       rep = fx$fit$rep, sd_rep = fx$sd_rep, obj = fx$fit,
                       n_sims = n_sims, newton_loops = 1, what = "SSB",
                       what_par = what_par, sim_type = type)
}


test_that("joint gives every replicate its own parameters and conditional gives them all the fit", {

  fx <- routing_fit()
  set.seed(1)
  vj <- SPoRC:::sim_draw_views("joint", 4, fx$fit$rep, fx$inp$par, fx$inp$map, fx$sd_rep, NULL, fx$fit)
  vc <- SPoRC:::sim_draw_views("conditional", 4, fx$fit$rep, fx$inp$par, fx$inp$map, fx$sd_rep, NULL, fx$fit)

  # conditional hands the same fitted list to every replicate
  for(i in 2:4) expect_identical(vc$pars[[i]], vc$pars[[1]])
  expect_identical(vc$pars[[1]], vc$fit_pars)

  # joint hands each a different one, and the fitted list is still available separately
  for(i in 2:4) expect_false(isTRUE(all.equal(vj$pars[[i]]$ln_srv_q, vj$pars[[1]]$ln_srv_q)))
  expect_equal(vj$fit_pars, vc$fit_pars)

  # a fit with no random effects has no joint precision, so the fixed-effect covariance stands in
  expect_null(fx$sd_rep$jointPrecision)

})


test_that("the drawn parameters reach the operating model's own inputs", {

  fx <- routing_fit()
  sl_j <- capture_sim_list(function() run_self_test(fx, "joint"))
  sl_c <- capture_sim_list(function() run_self_test(fx, "conditional"))

  # R0 and the seasonal split were built from the fitted report until 2026-09-28, so every replicate
  # ran at the same recruitment scale whatever it drew
  for(nm in c("R0", "Rec_input", "Fmort", "fish_sel", "srv_sel")) {
    expect_gt(rep_spread(sl_j[[nm]]), 0)
    expect_equal(rep_spread(sl_c[[nm]]), 0)
  }

  # rinit is mapped off here under use_rinit = 0, so it has nothing to draw and stays put
  expect_equal(rep_spread(sl_j$rinit), 0)

})


test_that("a parameter the data pin down returns at its own draw", {

  fx <- routing_fit()
  res <- run_self_test(fx, "joint", n_sims = 6, what_par = "ln_srv_q")

  est <- as.numeric(res$ln_srv_q)
  tru <- as.numeric(res$truth$ln_srv_q)
  fitted_q <- as.numeric(fx$fit$env$parList(par = fx$fit$env$last.par.best)$ln_srv_q)

  expect_gt(stats::sd(tru), 0) # the draws moved at all
  # each estimate sits closer to the draw that generated it than to the fitted value
  expect_lt(mean(abs(est - tru)), mean(abs(est - fitted_q)))

})


test_that("conditional scores against the fit and joint against each replicate's draw", {

  fx <- routing_fit()
  rc <- run_self_test(fx, "conditional", what_par = "ln_srv_q")
  rj <- run_self_test(fx, "joint", what_par = "ln_srv_q")

  expect_equal(stats::sd(as.numeric(rc$truth$ln_srv_q)), 0) # one truth throughout
  expect_gt(stats::sd(as.numeric(rj$truth$ln_srv_q)), 0)    # one per replicate

  expect_equal(dim(rc$SSB)[length(dim(rc$SSB))], 4)
  expect_equal(length(as.numeric(rc$ln_srv_q)), 4)
  expect_true(all(c("SSB", "ln_srv_q") %in% names(rc$truth)))

})


test_that("what_par refuses a name that is not a parameter", {

  fx <- routing_fit()
  expect_error(run_self_test(fx, "conditional", what_par = "not_a_parameter"), "what_par")

})
