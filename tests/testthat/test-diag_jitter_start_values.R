# Checks how jitter_start_values() splits a starting vector into fixed and random parts, and that a
# wrong-length par_vec errors instead of recycling. Uses a toy RTMB model, not a SPoRC one.

library(testthat)
library(SPoRC)
library(RTMB)   # dnorm has to resolve to RTMB's overload, not stats'

make_toy_re_obj <- function(seed = 1) {
  set.seed(seed)
  dat <- list(y = stats::rnorm(20))
  f <- function(p) {
    -sum(dnorm(p$u, 0, 1, log = TRUE)) -
      sum(dnorm(dat$y, p$mu + p$u, exp(p$ln_sd), log = TRUE))
  }
  RTMB::MakeADFun(f, list(mu = 0, ln_sd = 0, u = rep(0, 20)), random = "u", silent = TRUE)
}

test_that("jitter_start_values() returns the random start unperturbed unless asked to jitter it", {

  obj <- make_toy_re_obj()
  rand_idx <- obj$env$random

  no_re <- jitter_start_values(obj, par_vec = NULL, sd = 0.1, jitter_random = FALSE)
  expect_length(no_re$fixed, length(obj$par))
  expect_equal(unname(no_re$random), unname(obj$env$par[rand_idx]))

  with_re <- jitter_start_values(obj, par_vec = NULL, sd = 0.1, jitter_random = TRUE)
  expect_length(with_re$fixed, length(obj$par))
  expect_length(with_re$random, length(rand_idx))
  expect_false(isTRUE(all.equal(unname(with_re$random), unname(obj$env$par[rand_idx]))))
})

test_that("jitter_start_values() keeps par_vec's random values as the start when jitter_random is FALSE", {

  obj <- make_toy_re_obj()
  rand_idx <- obj$env$random

  # a joint par_vec sets where the inner solve starts, whether or not it is perturbed
  joint <- obj$env$par
  joint[rand_idx] <- 0.3
  out <- jitter_start_values(obj, par_vec = joint, sd = 0.1, jitter_random = FALSE)
  expect_equal(unname(out$random), rep(0.3, length(rand_idx)))
})

test_that("jitter_start_values() accepts a fixed-effect par_vec and a joint one, splitting the joint one on the random indices", {

  obj <- make_toy_re_obj()
  rand_idx <- obj$env$random

  fixed_only <- jitter_start_values(obj, par_vec = obj$par, sd = 0, jitter_random = TRUE)
  expect_equal(unname(fixed_only$fixed), unname(obj$par))

  joint <- obj$env$par
  joint[rand_idx] <- 0.5
  split <- jitter_start_values(obj, par_vec = joint, sd = 0, jitter_random = TRUE)
  expect_equal(unname(split$fixed), unname(joint[-rand_idx]))
  expect_equal(unname(split$random), rep(0.5, length(rand_idx)))
})

test_that("jitter_start_values() rejects a par_vec that is neither the fixed nor the joint length", {

  obj <- make_toy_re_obj()
  # the old code recycled here, producing a vector nlminb could not hand back to obj$fn
  expect_error(jitter_start_values(obj, par_vec = rep(0, 7), sd = 0.1, jitter_random = TRUE),
               "par_vec has length 7")
})

test_that("jitter_start_values() returns no random values when the model has no random effects", {

  set.seed(2)
  dat <- list(y = stats::rnorm(20))
  obj <- RTMB::MakeADFun(function(p) -sum(dnorm(dat$y, p$mu, exp(p$ln_sd), log = TRUE)),
                         list(mu = 0, ln_sd = 0), silent = TRUE)

  out <- jitter_start_values(obj, par_vec = NULL, sd = 0.1, jitter_random = TRUE)
  expect_length(out$fixed, length(obj$par))
  expect_length(out$random, 0)
})

test_that("the drawn random effects can be seeded into last.par, where the inner solve starts", {

  obj <- make_toy_re_obj()
  rand_idx <- obj$env$random
  draws <- jitter_start_values(obj, par_vec = NULL, sd = 0.5, jitter_random = TRUE)$random

  obj$env$last.par[rand_idx] <- draws
  expect_equal(unname(obj$env$last.par[rand_idx]), unname(draws))
})
