# Checks that asking fit_model and do_jitter for zero Newton loops takes no Newton step. Both loops were written
# 1:n, which at zero runs over c(1, 0) and took two steps, so a fit asked for none came back refined.

library(testthat)
library(SPoRC)
data("dusky_rtmb_model")

data <- dusky_rtmb_model$data
parameters <- dusky_rtmb_model$parameters
mapping <- dusky_rtmb_model$mapping

test_that("fit_model with newton_loops = 0 returns the nlminb optimum untouched", {

  fit <- fit_model(data, parameters, mapping, random = NULL, newton_loops = 0, silent = TRUE)

  # the same optimization by hand, on an object built the way fit_model builds it and with its default control
  seed <- fit_model(data, parameters, mapping, random = NULL, do_optim = FALSE, silent = TRUE)
  opt <- stats::nlminb(seed$par, seed$fn, seed$gr,
                       control = list(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-12),
                       lower = rep(-Inf, length(seed$par)),
                       upper = rep(Inf, length(seed$par)))

  expect_identical(fit$optim$par, opt$par)
  expect_identical(fit$optim$objective, opt$objective)

})

test_that("do_jitter with n_newton_loops = 0 reports the nlminb optimum untouched", {

  # sd = 0 starts the one run from the model's own start
  jit <- do_jitter(data, parameters, mapping, NULL, sd = 0, n_jitter = 1, n_newton_loops = 0, do_par = FALSE)

  # the same optimization by hand, with do_jitter's control
  obj <- RTMB::MakeADFun(cmb(SPoRC_rtmb, data), parameters = parameters, map = mapping, silent = TRUE)
  opt <- stats::nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-15))
  sd_rep <- RTMB::sdreport(obj)

  expect_identical(unique(jit$jnLL), obj$report(obj$env$last.par.best)$jnLL)
  expect_identical(unique(jit$Max_Gradient), max(abs(sd_rep$gradient.fixed))) # a Newton step moves this most

})
