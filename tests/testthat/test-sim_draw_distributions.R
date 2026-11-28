# Operating model draws, measured over many replicates: AR1 catchability deviations start at the stationary sd and
# step at sigma (5 pct), resampled recruitment keeps the arithmetic and harmonic means of its series (1 pct), AR1
# recruitment deviations follow the bias ramp (5 pct), fishing mortality deviations keep the fitted years and step
# on from the last one (5 pct), and a bias-corrected lognormal index is centered on the truth (0.5 pct).

library(SPoRC)
library(testthat)

test_that("AR1 catchability deviations are drawn at the stationary sd from the first year on", {

  n_yrs <- 6
  n_sims <- 4000
  sigma <- 0.3 # sd of each year's innovation
  rho <- 0.8 # correlation between years

  # what draw_sim_q_devs reads, for one region and one survey fleet with no years kept from a fit
  sim_env <- new.env()
  sim_env$n_regions <- 1
  sim_env$n_cond_yrs <- 0
  sim_env$srv_q_model <- 4 # ar1
  sim_env$ln_srv_q_devs <- array(0, dim = c(1, n_yrs, 1, n_sims))
  sim_env$srv_q <- array(2, dim = c(1, n_yrs, 1, n_sims)) # mean catchability
  sim_env$sigma_srv_q <- matrix(sigma, nrow = 1, ncol = 1)
  sim_env$srv_q_rho <- matrix(rho, nrow = 1, ncol = 1)

  set.seed(11)
  for(sim in 1:n_sims) SPoRC:::draw_sim_q_devs(sim, sim_env)
  devs <- sim_env$ln_srv_q_devs[1,,1,] # [year, replicate]

  # every year has the stationary sd, the first included
  stationary_sd <- sigma / sqrt(1 - rho^2)
  expect_equal(stats::sd(devs[1,]), stationary_sd, tolerance = 0.05)
  expect_equal(stats::sd(devs[n_yrs,]), stationary_sd, tolerance = 0.05)

  # each later year is rho times the year before plus an innovation at sigma
  expect_equal(stats::cor(devs[1,], devs[2,]), rho, tolerance = 0.03)
  expect_equal(stats::sd(devs[2,] - rho * devs[1,]), sigma, tolerance = 0.05)

  # catchability is the mean times the exponentiated deviation
  expect_equal(sim_env$srv_q[1,,1,], 2 * exp(devs), tolerance = 1e-12)

})

test_that("resampled recruitment keeps the arithmetic and harmonic means of the series", {

  recruitment <- c(5, 8, 12, 20, 9, 14, 30, 7) # a right-skewed series

  set.seed(12)
  draws <- SPoRC:::rinvgauss_rec(sims = 200000, recruitment = recruitment)

  expect_true(all(draws > 0))
  expect_equal(mean(draws), mean(recruitment), tolerance = 0.01)
  expect_equal(1 / mean(1 / draws), 1 / mean(1 / recruitment), tolerance = 0.01)

  # the two means differ by a third here, so a draw centered on the wrong one would fail
  expect_gt(mean(recruitment) / (1 / mean(1 / recruitment)), 1.3)

})

test_that("AR1 recruitment deviations are drawn about a mean that follows the bias ramp", {

  n_yrs <- 5
  n_sims <- 6000
  sigma <- 0.8 # sd of each year's innovation
  rho <- 0.6 # correlation between years
  ramp <- c(0, 0.3, 0.6, 1, 1) # bias ramp by year

  # what rec_devs_past_fit reads from a fit, for one population and region with no years kept
  data <- list(RecDevs_model = 3) # ar1
  pars <- list(ln_RecDevs = array(0, dim = c(1, 1, n_yrs)),
               ln_sigmaR = array(log(sigma), dim = c(2, 1, 1)),
               RecDevs_rho = array(atanh(rho), dim = c(1, 1))) # the model stores rho before its tanh transform
  rep <- list(bias_ramp = ramp, rec_region_prop = matrix(1, nrow = 1, ncol = 1))

  set.seed(13)
  devs <- matrix(0, n_yrs, n_sims) # [year, replicate]
  for(sim in 1:n_sims) devs[,sim] <- SPoRC:::rec_devs_past_fit(data = data, pars = pars, rep = rep, n_cond_yrs = 0, n_yrs = n_yrs)[1,1,]

  # each year's mean is minus half the stationary variance times that year's ramp, at the stationary sd
  stationary_var <- sigma^2 / (1 - rho^2)
  center <- -stationary_var / 2 * ramp
  for(y in 1:n_yrs) {
    expect_lt(abs(mean(devs[y,]) - center[y]), 0.05)
    expect_equal(stats::sd(devs[y,]), sqrt(stationary_var), tolerance = 0.05)
  } # end y loop

  # each later year is rho times last year's distance from its own mean, plus an innovation at sigma
  for(y in 2:n_yrs) {
    innovation <- (devs[y,] - center[y]) - rho * (devs[y - 1,] - center[y - 1])
    expect_lt(abs(mean(innovation)), 0.04)
    expect_equal(stats::sd(innovation), sigma, tolerance = 0.05)
  } # end y loop

  # the ramp moves the mean by half a unit across years here, far more than the tolerance
  expect_gt(abs(center[n_yrs] - center[1]), 0.4)

})

test_that("fishing mortality deviations keep the fitted years and are drawn on from the last of them", {

  fit_devs <- c(0.5, -0.2, 0.3, 0, 0, 0) # the fit's deviations, with the last three years to draw
  n_cond <- 3 # years that keep the fit's
  n_sims <- 4000
  sigma <- 0.25
  rho <- 0.6
  all_years <- rep(TRUE, length(fit_devs))

  # draws [year, replicate] under one process
  draw_series <- function(PE_model, is_est, seed) {
    set.seed(seed)
    draws <- matrix(0, length(fit_devs), n_sims)
    for(sim in 1:n_sims) draws[,sim] <- SPoRC:::draw_F_dev_series(devs = fit_devs, is_est = is_est, PE_model = PE_model, sigma = sigma, rho = rho, n_cond = n_cond)
    draws
  }

  # independent: the fitted years are returned as they were, the last of them included, and the rest are drawn about zero
  draws <- draw_series(PE_model = 1, is_est = all_years, seed = 51)
  for(y in 1:n_cond) expect_equal(draws[y,], rep(fit_devs[y], n_sims))
  expect_lt(abs(mean(draws[4,])), 4 * sigma / sqrt(n_sims))
  expect_equal(stats::sd(draws[4,]), sigma, tolerance = 0.05)

  # random walk: year 4 steps from the fit's year 3, and year 5 from year 4
  draws <- draw_series(PE_model = 2, is_est = all_years, seed = 52)
  expect_equal(draws[3,], rep(fit_devs[3], n_sims))
  expect_lt(abs(mean(draws[4,]) - fit_devs[3]), 4 * sigma / sqrt(n_sims))
  expect_equal(stats::sd(draws[4,]), sigma, tolerance = 0.05)
  expect_equal(stats::sd(draws[5,] - draws[4,]), sigma, tolerance = 0.05)

  # AR1: year 4 is rho times the fit's year 3, plus an innovation at sigma
  draws <- draw_series(PE_model = 3, is_est = all_years, seed = 53)
  expect_lt(abs(mean(draws[4,]) - rho * fit_devs[3]), 4 * sigma / sqrt(n_sims))
  expect_equal(stats::sd(draws[4,]), sigma, tolerance = 0.05)

  # a year the fit does not estimate is skipped, and the next one is drawn across the two year gap
  with_gap <- c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE)
  draws <- draw_series(PE_model = 2, is_est = with_gap, seed = 54)
  expect_equal(draws[4,], rep(fit_devs[4], n_sims)) # not estimated, so left as it was
  expect_equal(stats::sd(draws[5,]), sigma * sqrt(2), tolerance = 0.05)
  draws <- draw_series(PE_model = 3, is_est = with_gap, seed = 55)
  expect_lt(abs(mean(draws[5,]) - rho^2 * fit_devs[3]), 4 * sigma * sqrt(1 + rho^2) / sqrt(n_sims))
  expect_equal(stats::sd(draws[5,]), sigma * sqrt(1 + rho^2), tolerance = 0.05)

})

test_that("a lognormal index draw is centered on the truth when bias corrected, and above it when not", {

  truth <- 100
  se <- 0.5 # log scale
  n_draws <- 200000

  set.seed(56)
  corrected <- SPoRC:::draw_index_obs(true = rep(truth, n_draws), se = se, like_type = 0, bias_correct_oe = 1)
  set.seed(56)
  uncorrected <- SPoRC:::draw_index_obs(true = rep(truth, n_draws), se = se, like_type = 0, bias_correct_oe = 0)

  expect_equal(mean(corrected) / truth, 1, tolerance = 0.005)
  expect_equal(mean(uncorrected) / truth, exp(0.5 * se^2), tolerance = 0.005) # 13 pct above the truth
  expect_equal(stats::sd(log(corrected)), se, tolerance = 0.01)

  # a normal index is centered on the truth with that standard deviation
  set.seed(57)
  normal <- SPoRC:::draw_index_obs(true = rep(truth, n_draws), se = 8, like_type = 1)
  expect_equal(mean(normal), truth, tolerance = 0.001)
  expect_equal(stats::sd(normal), 8, tolerance = 0.01)

})
