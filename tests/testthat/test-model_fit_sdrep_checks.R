# Checks fit_model(do_sdrep = ): FALSE (the default) attaches no sdreport and runs no checks, TRUE attaches the sdreport and
# the result of post_optim_sanity_checks on it, and post_optim_sanity_checks_control reaches the checks, with defaults for any
# tolerance left out and an error for any other name.

library(testthat)
library(SPoRC)
data("dusky_rtmb_model")

data <- dusky_rtmb_model$data
parameters <- dusky_rtmb_model$parameters
mapping <- dusky_rtmb_model$mapping

test_that("do_sdrep = FALSE attaches no sdreport and prints nothing", {

  expect_no_message(fit <- fit_model(data, parameters, mapping, newton_loops = 1, silent = TRUE))
  expect_null(fit$sd_rep)
  expect_null(fit$post_optim_sanity_checks)

})

test_that("do_sdrep = TRUE attaches the sdreport and the checks on it", {

  expect_message(fit <- fit_model(data, parameters, mapping, newton_loops = 1, silent = TRUE, do_sdrep = TRUE),
                 "Successfully passed")
  expect_s3_class(fit$sd_rep, "sdreport")
  expect_equal(length(fit$sd_rep$par.fixed), length(fit$par))
  expect_equal(unname(fit$sd_rep$par.fixed), unname(fit$env$last.par.best)) # same parameters $rep is evaluated at
  expect_true(fit$post_optim_sanity_checks)

})

test_that("post_optim_sanity_checks_control reaches the checks", {

  # every fit has a gradient above zero, so a zero tolerance must fail and name the parameter
  tight <- list(gradient_tol = 0, se_tol = 100, corr_tol = 0.99)
  msgs <- capture_messages(fit <- fit_model(data, parameters, mapping, newton_loops = 1, silent = TRUE, do_sdrep = TRUE,
                                            post_optim_sanity_checks_control = tight))
  expect_false(fit$post_optim_sanity_checks)
  expect_true(any(grepl("absolute gradient", msgs)))
  expect_false(any(grepl("Successfully passed", msgs)))

})

test_that("a control list naming only some tolerances takes the defaults for the rest", {

  # se_tol and corr_tol are left out, so they reach the checks as their defaults rather than NULL
  msgs <- capture_messages(fit <- fit_model(data, parameters, mapping, newton_loops = 1, silent = TRUE, do_sdrep = TRUE,
                                            post_optim_sanity_checks_control = list(gradient_tol = 0)))
  expect_false(fit$post_optim_sanity_checks)
  expect_true(any(grepl("absolute gradient", msgs)))
  expect_false(any(grepl("standard error|Parameter pairs", msgs))) # the default se and correlation tolerances pass

})

test_that("a control list element that is not one of the three tolerances stops the fit", {

  expect_error(fit_model(data, parameters, mapping, silent = TRUE, do_sdrep = TRUE,
                         post_optim_sanity_checks_control = list(gradient_tol = 1e-4, grad_tol = 1e-4)), "'grad_tol'")
  expect_error(fit_model(data, parameters, mapping, silent = TRUE, do_sdrep = TRUE,
                         post_optim_sanity_checks_control = list(1e-4)), "takes only gradient_tol, se_tol and corr_tol")
  expect_error(fit_model(data, parameters, mapping, silent = TRUE, do_optim = FALSE, # stops before the fit, sdreport or not
                         post_optim_sanity_checks_control = list(se_tolerance = 10)), "'se_tolerance'")

})
