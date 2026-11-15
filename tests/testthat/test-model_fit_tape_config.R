# Checks fit_model's tape_config and tmb_config: each setting gives the default model's objective and gradient, the
# session's settings are restored after the fit and after a build that stops, NULL leaves them alone, and a misspelled
# name or value stops before anything is built.

library(testthat)
library(SPoRC)
data("dusky_rtmb_model")

data <- dusky_rtmb_model$data
parameters <- dusky_rtmb_model$parameters
mapping <- dusky_rtmb_model$mapping

test_that("each tape_config setting gives the default objective and gradient and restores the session settings", {

  tape_before <- RTMB::TapeConfig() # no arguments reads the current settings
  obj_default <- fit_model(data, parameters, mapping, do_optim = FALSE, silent = TRUE)
  fn_default <- obj_default$fn(obj_default$par)
  gr_default <- as.numeric(obj_default$gr(obj_default$par))

  tape_settings <- list(list(vectorize = "enable"),
                        list(atomic = "enable"),
                        list(atomic = "disable"),
                        list(comparison = "tape"),
                        list(atomic = "enable", vectorize = "enable"))

  for(i in seq_along(tape_settings)) {
    obj <- fit_model(data, parameters, mapping, do_optim = FALSE, silent = TRUE, tape_config = tape_settings[[i]])
    expect_equal(obj$fn(obj$par), fn_default, tolerance = 1e-12) # recorded differently, same objective
    expect_equal(as.numeric(obj$gr(obj$par)), gr_default, tolerance = 1e-8) # and gradient, to rounding
    expect_identical(RTMB::TapeConfig(), tape_before) # put back on exit
  } # end i loop

})

test_that("tmb_config sparse Hessian compression gives the same Laplace objective and restores TMB's settings", {

  tmb_before <- TMB::config(DLL = "RTMB")
  obj_default <- fit_model(data, parameters, mapping, random = "ln_RecDevs", do_optim = FALSE, silent = TRUE)
  obj_compress <- fit_model(data, parameters, mapping, random = "ln_RecDevs", do_optim = FALSE, silent = TRUE,
                            tmb_config = list(tmbad.sparse_hessian_compress = 1))

  expect_equal(obj_compress$fn(obj_compress$par), obj_default$fn(obj_default$par), tolerance = 1e-12)
  expect_equal(as.numeric(obj_compress$gr(obj_compress$par)), as.numeric(obj_default$gr(obj_default$par)), tolerance = 1e-8)
  expect_identical(TMB::config(DLL = "RTMB"), tmb_before)

})

test_that("settings are restored when the build stops", {

  tape_before <- RTMB::TapeConfig()
  tmb_before <- TMB::config(DLL = "RTMB")

  # a model that stops while it is being recorded
  stop_model <- function(pars, data) stop("stop inside the model")
  expect_error(fit_model(data, parameters, mapping, model = stop_model, silent = TRUE,
                         tape_config = list(vectorize = "enable"),
                         tmb_config = list(tmbad.sparse_hessian_compress = 1)),
               "stop inside the model")

  expect_identical(RTMB::TapeConfig(), tape_before)
  expect_identical(TMB::config(DLL = "RTMB"), tmb_before)

})

test_that("NULL leaves a setting made for the session alone", {

  tape_before <- RTMB::TapeConfig()
  RTMB::TapeConfig(vectorize = "enable") # set by the user before fitting
  tape_user <- RTMB::TapeConfig()

  fit_model(data, parameters, mapping, do_optim = FALSE, silent = TRUE)
  expect_identical(RTMB::TapeConfig(), tape_user)

  RTMB::TapeConfig(tape_before) # put the session back for the other tests
  expect_identical(RTMB::TapeConfig(), tape_before)

})

test_that("a misspelled name or value stops before the model is built, naming the argument", {

  expect_error(fit_model(data, parameters, mapping, do_optim = FALSE, silent = TRUE, tape_config = list(vectorise = "enable")),
               "tape_config takes only comparison, atomic and vectorize, not: 'vectorise'")
  expect_error(fit_model(data, parameters, mapping, do_optim = FALSE, silent = TRUE, tape_config = list(vectorize = "on")),
               "tape_config\\$vectorize must be one of")
  expect_error(fit_model(data, parameters, mapping, do_optim = FALSE, silent = TRUE, tape_config = list("enable")),
               "tape_config takes only")
  expect_error(fit_model(data, parameters, mapping, do_optim = FALSE, silent = TRUE, tape_config = "enable"),
               "tape_config must be NULL or a named list")
  expect_error(fit_model(data, parameters, mapping, do_optim = FALSE, silent = TRUE, tmb_config = list(sparse_hessian_compress = 1)),
               "tmb_config takes only .* not: 'sparse_hessian_compress'")
  expect_error(fit_model(data, parameters, mapping, do_optim = FALSE, silent = TRUE, tmb_config = list(tmbad.sparse_hessian_compress = "yes")),
               "tmb_config\\$tmbad.sparse_hessian_compress must be a single number")

})
