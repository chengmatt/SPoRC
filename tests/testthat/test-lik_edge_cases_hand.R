# Likelihood pieces at their edges, against values worked out by hand (1e-10): the logistic normal with empty bins
# dropped when one, two or three bins are left, and the total index standard deviation under each form.

library(SPoRC)
library(testthat)

test_that("the logistic normal with empty bins dropped is fit from two positive bins up", {

  pred <- c(0.1, 0.2, 0.3, 0.4) # expected proportions
  sigma <- 0.8 # per fish standard deviation
  ISS <- 50 # input sample size
  variance <- sigma^2 / ISS

  # one positive bin leaves no log ratio, so it adds nothing
  expect_equal(SPoRC:::get_logistnormal_miss0_nLL(obs = c(0, 1, 0, 0), pred = pred, ln_sigma = log(sigma), ISS = ISS), 0)

  # two positive bins, 2 and 4: one log ratio against the last positive bin, plus the change of variables
  obs <- c(0, 0.3, 0, 0.7)
  resid <- log(0.3 / 0.7) - log(pred[2] / pred[4])
  by_hand <- -stats::dnorm(resid, mean = 0, sd = sqrt(variance), log = TRUE) + log(0.3) + log(0.7)
  expect_equal(as.numeric(SPoRC:::get_logistnormal_miss0_nLL(obs = obs, pred = pred, ln_sigma = log(sigma), ISS = ISS)), by_hand, tolerance = 1e-10)
  expect_gt(abs(by_hand), 0.1) # far from the zero a skipped cell would give

  # three positive bins, 1, 2 and 4, independent
  obs <- c(0.2, 0.3, 0, 0.5)
  resid <- c(log(0.2 / 0.5) - log(pred[1] / pred[4]), log(0.3 / 0.5) - log(pred[2] / pred[4]))
  by_hand <- -sum(stats::dnorm(resid, mean = 0, sd = sqrt(variance), log = TRUE)) + log(0.2) + log(0.3) + log(0.5)
  expect_equal(as.numeric(SPoRC:::get_logistnormal_miss0_nLL(obs = obs, pred = pred, ln_sigma = log(sigma), ISS = ISS)), by_hand, tolerance = 1e-10)

  # the same three bins under an AR1 across bins: the second residual given the first, one bin apart
  trans_rho <- 0.4
  rho <- stats::plogis(trans_rho)
  by_hand <- -(stats::dnorm(resid[1], mean = 0, sd = sqrt(variance), log = TRUE) +
                 stats::dnorm(resid[2], mean = rho * resid[1], sd = sqrt(variance * (1 - rho^2)), log = TRUE)) + log(0.2) + log(0.3) + log(0.5)
  expect_equal(as.numeric(SPoRC:::get_logistnormal_miss0_nLL(obs = obs, pred = pred, ln_sigma = log(sigma), ISS = ISS, corr_type = 1, trans_rho = trans_rho)),
               by_hand, tolerance = 1e-10)

})

test_that("the total index standard deviation combines the reported and estimated parts under each form", {

  # reported standard errors [region, year, season, fleet] and one estimated value per fleet
  se <- array(c(0.10, 0.20, 0.30, 0.15, 0.25, 0.35), dim = c(1, 3, 1, 2))
  extra <- c(0.05, 0.40)
  ln_sigma <- log(extra)

  # form 0 keeps the reported values
  expect_equal(SPoRC:::build_idx_sd(se = se, ln_sigma = ln_sigma, form = 0), se)

  for(f in 1:2) {
    expect_equal(SPoRC:::build_idx_sd(se = se, ln_sigma = ln_sigma, form = 1)[1,,1,f], se[1,,1,f] + extra[f], tolerance = 1e-12) # added
    expect_equal(SPoRC:::build_idx_sd(se = se, ln_sigma = ln_sigma, form = 2)[1,,1,f], sqrt(se[1,,1,f]^2 + extra[f]^2), tolerance = 1e-12) # added in variance
    expect_equal(SPoRC:::build_idx_sd(se = se, ln_sigma = ln_sigma, form = 3)[1,,1,f], rep(extra[f], 3), tolerance = 1e-12) # replaced
  } # end f loop

  # a data source with no standard errors stays empty under every form
  for(form in 0:3) expect_null(SPoRC:::build_idx_sd(se = NULL, ln_sigma = ln_sigma, form = form))

})
