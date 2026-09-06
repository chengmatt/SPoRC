library(SPoRC)
library(testthat)

# Operating model to estimation model self test for the state-space numbers at age. The population
# is advanced with a known process error, index and composition data are generated from it, and the
# estimating model has to recover the process from those data alone with the states integrated out.
#
# This asks a harder question than test-selftest_naa_state.R, which handed over the true states and
# estimated only the process parameters. Here the process variance has to be separated from
# observation error, which is the part of a state-space fit that fails in practice.

test_that("the estimation model recovers a simulated numbers-at-age process error", {
  om <- naaom_make_om(NAA_re = "iid", sigmaNAA = 0.25, rho_age = 0, rho_year = 0, seed = 808)
  dat <- naaom_om_data(om)

  il_det <- naaom_build_em(dat, NAA_re = "none")
  det <- fit_model(il_det$data, il_det$par, il_det$map, random = NULL, silent = TRUE)

  il <- naaom_build_em(dat, NAA_re = "iid")
  fit <- suppressWarnings(fit_model(il$data, il$par, il$map, random = "ln_NAA", silent = TRUE))

  expect_lt(max(abs(fit$gr(naaom_fixed(fit)))), 1e-3)
  expect_equal(exp(naaom_fixed(fit, "ln_sigmaNAA")), 0.25, tolerance = 0.25)
  expect_equal(length(fit$env$random), il$data$n_est_naa_re)
})

test_that("a correlated process error is recovered and improves the biomass trajectory", {
  om <- naaom_make_om(NAA_re = "2dar1", sigmaNAA = 0.25, rho_age = 0.5, rho_year = 0.4, seed = 808)
  dat <- naaom_om_data(om)

  il_det <- naaom_build_em(dat, NAA_re = "none")
  det <- fit_model(il_det$data, il_det$par, il_det$map, random = NULL, silent = TRUE)

  il <- naaom_build_em(dat, NAA_re = "2dar1")
  fit <- suppressWarnings(fit_model(il$data, il$par, il$map, random = "ln_NAA", silent = TRUE))

  expect_lt(max(abs(fit$gr(naaom_fixed(fit)))), 1e-3)
  expect_equal(exp(naaom_fixed(fit, "ln_sigmaNAA")), 0.25, tolerance = 0.25)

  # the correlations are weakly informed from data alone, so this is a sign-and-scale check rather
  # than a tight one: both are estimated positive and short of the boundary
  rho <- 2 / (1 + exp(-2 * naaom_fixed(fit, "NAA_pe_pars"))) - 1
  expect_true(all(rho[1:2] > 0))
  expect_true(all(rho[1:2] < 0.95))

  # a population with real process error is badly served by a model that denies it
  expect_lt(naaom_ssb_rmse(fit, om), naaom_ssb_rmse(det, om) / 2)
})

test_that("the state does not damage a trajectory that has no process error to find", {
  # turning the state on when the operating model is deterministic should cost little, which is
  # what makes it safe to try: the penalty pulls the innovations to zero and the fit stands
  om <- naaom_make_om(NAA_re = "none", seed = 404)
  dat <- naaom_om_data(om)

  il_det <- naaom_build_em(dat, NAA_re = "none")
  det <- fit_model(il_det$data, il_det$par, il_det$map, random = NULL, silent = TRUE)

  il <- naaom_build_em(dat, NAA_re = "iid")
  fit <- suppressWarnings(fit_model(il$data, il$par, il$map, random = "ln_NAA", silent = TRUE))

  expect_lt(max(abs(fit$gr(naaom_fixed(fit)))), 1e-3)
  expect_lt(naaom_ssb_rmse(fit, om), naaom_ssb_rmse(det, om) * 1.5)
})
