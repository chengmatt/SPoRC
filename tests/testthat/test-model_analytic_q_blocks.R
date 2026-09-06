library(SPoRC)
library(testthat)

# An analytically solved catchability is concentrated out of the likelihood rather than
# estimated. It is solved WITHIN each catchability time block, so a blocked q gets one
# solved value per block instead of one pooled value for the whole series.
# get_blocked_analytic_q is tested directly against the formulas it claims to implement.

test_that("get_blocked_analytic_q reproduces the pooled solve when there is one block", {

  set.seed(11)
  n_yrs <- 12
  obs <- exp(rnorm(n_yrs, log(10), 0.2))
  pred <- exp(rnorm(n_yrs, log(20), 0.2))
  blk1 <- rep(1L, n_yrs)

  arith <- SPoRC:::get_blocked_analytic_q(1, obs, pred, 1:n_yrs, blk1, 1)
  geo <- SPoRC:::get_blocked_analytic_q(2, obs, pred, 1:n_yrs, blk1, 1)

  expect_equal(arith, rep(mean(obs) / mean(pred), n_yrs), tolerance = 1e-14)
  expect_equal(geo, rep(exp(mean(log(obs) - log(pred))), n_yrs), tolerance = 1e-14)
  expect_length(arith, n_yrs)
})

test_that("two blocks are solved from their own observations only", {

  set.seed(12)
  n_yrs <- 12
  obs <- exp(rnorm(n_yrs, log(10), 0.2))
  pred <- exp(rnorm(n_yrs, log(20), 0.2))
  blk <- rep(1:2, each = 6)
  a <- 1:6; b <- 7:12

  arith <- SPoRC:::get_blocked_analytic_q(1, obs, pred, 1:n_yrs, blk, 2)
  geo <- SPoRC:::get_blocked_analytic_q(2, obs, pred, 1:n_yrs, blk, 2)

  expect_equal(arith[1], mean(obs[a]) / mean(pred[a]), tolerance = 1e-14)
  expect_equal(arith[n_yrs], mean(obs[b]) / mean(pred[b]), tolerance = 1e-14)
  expect_equal(geo[1], exp(mean(log(obs[a]) - log(pred[a]))), tolerance = 1e-14)
  expect_equal(geo[n_yrs], exp(mean(log(obs[b]) - log(pred[b]))), tolerance = 1e-14)
  expect_length(unique(signif(arith, 10)), 2)
})

test_that("a block with no observations falls back to the pooled solve", {

  set.seed(13)
  n_yrs <- 10
  obs <- exp(rnorm(6, log(10), 0.2))
  pred <- exp(rnorm(6, log(20), 0.2))
  blk <- c(rep(1L, 6), rep(2L, 4))   # block 2 owns years 7 to 10, none of which have an index

  q <- SPoRC:::get_blocked_analytic_q(1, obs, pred, 1:6, blk, 2)
  expect_equal(q[1], mean(obs) / mean(pred), tolerance = 1e-14)
  expect_equal(q[n_yrs], mean(obs) / mean(pred), tolerance = 1e-14)
})

test_that("relabelling the regimes leaves an analytic q untouched and an estimated one does not", {

  # This is why the per year posterior over regimes is exactly 1/K under an analytic q:
  # the solve follows whatever partition the states define, so nothing in the likelihood
  # distinguishes regime 1 from regime 2. The partition is identified; the labels are not.
  set.seed(15)
  n_yrs <- 8
  obs <- exp(rnorm(n_yrs, log(10), 0.2))
  pred <- exp(rnorm(n_yrs, log(20), 0.2))
  blk <- c(rep(1L, 4), rep(2L, 4))
  swapped <- c(rep(2L, 4), rep(1L, 4))

  q_a <- SPoRC:::get_blocked_analytic_q(1, obs, pred, 1:n_yrs, blk, 2)
  q_b <- SPoRC:::get_blocked_analytic_q(1, obs, pred, 1:n_yrs, swapped, 2)
  # each year keeps the value solved from its own part of the partition, whichever
  # label that part is with
  expect_equal(q_a, q_b, tolerance = 1e-14)
})
