# Checks fit_model(model = ) on a model of the user's own. A data list Setup_Mod_Dim did not build used to stop in
# fill_optional_modules, which read n_regions it does not have; such a list now skips the SPoRC setup steps, while a
# SPoRC data list with an edited copy of the objective still takes them.

library(testthat)
library(SPoRC)
data("dusky_rtmb_model")

# random walk state observed with error, nothing from SPoRC
set.seed(1)
n_obs <- 50
toy_data <- list(y = cumsum(rnorm(n_obs, 0, 0.3)) + rnorm(n_obs, 0, 0.5), # observations
                 map_u = rep(0, n_obs)) # named like a SPoRC map copy, which must be left alone
toy_pars <- list(u = rep(0, n_obs), ln_sd_u = 0, ln_sd_y = 0)
toy_model <- function(pars, data) {
  RTMB::getAll(pars, data)
  nll = -sum(dnorm(diff(u), 0, exp(ln_sd_u), log = TRUE)) # random walk
  nll = nll - sum(dnorm(y, u, exp(ln_sd_y), log = TRUE)) # observations
  RTMB::REPORT(u)
  nll
}

test_that("a model of the user's own with random effects fits as RTMB fits it by hand", {

  ref <- RTMB::MakeADFun(cmb(toy_model, toy_data), toy_pars, random = "u", silent = TRUE)
  opt <- stats::nlminb(ref$par, ref$fn, ref$gr)

  for(toy_map in list(list(), NULL)) {
    fit <- suppressMessages(fit_model(toy_data, toy_pars, mapping = toy_map, random = "u", model = toy_model,
                                      silent = TRUE, do_sdrep = TRUE))
    expect_equal(fit$optim$objective, opt$objective, tolerance = 1e-8)
    expect_equal(as.numeric(fit$optim$par), as.numeric(opt$par), tolerance = 1e-5) # Newton steps refine past nlminb
    expect_length(fit$rep$u, n_obs)
    expect_true(fit$post_optim_sanity_checks)
    expect_identical(fit$data, toy_data) # no setup fields added, and map_u not rewritten from the map
  } # end toy_map loop

})

test_that("a model of the user's own with fixed effects takes the Newton steps", {

  # both sds fixed, so the objective is quadratic in u and one Newton step lands on the optimum
  toy_map <- list(ln_sd_u = factor(NA), ln_sd_y = factor(NA))
  fit <- fit_model(toy_data, toy_pars, mapping = toy_map, model = toy_model, silent = TRUE, newton_loops = 3)

  expect_length(fit$par, n_obs)
  expect_lt(max(abs(fit$gr(fit$optim$par))), 1e-8)
  expect_false(is.null(fit$hessian)) # the tape Hessian was used

})

test_that("a SPoRC data list with an edited copy of the objective still takes the setup steps", {

  data <- dusky_rtmb_model$data
  data$use_fixed_movement <- NULL # as if Setup_Mod_Movement had not been called

  copied_objective <- function(pars, data) SPoRC_rtmb(pars, data) # stands in for an edited copy, as in profiling
  obj <- fit_model(data, dusky_rtmb_model$parameters, dusky_rtmb_model$mapping, model = copied_objective,
                   do_optim = FALSE, silent = TRUE)

  expect_identical(obj$data$use_fixed_movement, 1) # one region, so movement filled with its off setting
  expect_equal(obj$fn(obj$par), dusky_rtmb_model$fn(dusky_rtmb_model$par), tolerance = 1e-10)

})
