# Checks the standard errors get_par_est_info reports against RTMB's own summary of the sdreport, for fixed effects and for
# random effects, whose variances the sdreport stores.

library(testthat)
library(SPoRC)
data("sgl_rg_dusky_data")

test_that("random effect standard errors match RTMB's summary of the sdreport", {

  input_list <- build_goa_dusky_input(sgl_rg_dusky_data)
  fit <- suppressMessages(fit_model(input_list$data, input_list$par, input_list$map, random = c("ln_RecDevs", "ln_InitDevs"),
                                    newton_loops = 0, silent = TRUE, do_sdrep = TRUE))

  est_pars <- get_par_est_info(fit$parameters, fit$mapping, fit$sd_rep)$est_pars
  sd_fixed <- summary(fit$sd_rep, "fixed")
  sd_random <- summary(fit$sd_rep, "random")

  # rows are fixed effects first, then random effects, in the sdreport's order
  n_fixed <- nrow(sd_fixed)
  n_random <- nrow(sd_random)
  expect_equal(nrow(est_pars), n_fixed + n_random)
  expect_equal(est_pars$Par, c(rownames(sd_fixed), rownames(sd_random)))

  expect_equal(est_pars$SE_Val[1:n_fixed], unname(sd_fixed[, "Std. Error"]))
  expect_equal(est_pars$SE_Val[n_fixed + 1:n_random], unname(sd_random[, "Std. Error"]))
  expect_equal(est_pars$Est_Val[n_fixed + 1:n_random], unname(sd_random[, "Estimate"]))

  # the variance is not the standard error here, so the old value would fail the line above
  expect_true(max(abs(sd_random[, "Std. Error"]^2 - sd_random[, "Std. Error"])) > 0.01)

})

test_that("a fit without random effects still returns one standard error per estimate", {

  input_list <- build_goa_dusky_input(sgl_rg_dusky_data)
  fit <- suppressMessages(fit_model(input_list$data, input_list$par, input_list$map, newton_loops = 0, silent = TRUE, do_sdrep = TRUE))

  est_pars <- get_par_est_info(fit$parameters, fit$mapping, fit$sd_rep)$est_pars
  expect_equal(nrow(est_pars), length(fit$par))
  expect_equal(est_pars$SE_Val, unname(summary(fit$sd_rep, "fixed")[, "Std. Error"]))

})
