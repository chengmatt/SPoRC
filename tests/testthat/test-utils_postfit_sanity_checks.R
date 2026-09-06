library(SPoRC)
library(testthat)

# post_optim_sanity_checks() is the gate that decides whether a fit is reported
# as converged. Its four criteria only ever fire on a bad fit, so they are
# driven here with hand-built sdreport stand-ins rather than a real model.

make_sd_rep <- function(
  gradient = c(1e-6, 1e-6),
  pdHess = TRUE,
  cov = diag(c(0.04, 0.09)),
  par_names = c("ln_M", "ln_q")
) {
  names(gradient) <- par_names
  dimnames(cov) <- list(par_names, par_names)
  list(
    gradient.fixed = gradient,
    par.fixed = stats::setNames(c(0.1, 0.2), par_names),
    pdHess = pdHess,
    cov.fixed = cov
  )
}

good_rep <- list(jnLL = 123.4)

test_that("a clean fit passes every check", {
  expect_true(suppressMessages(SPoRC::post_optim_sanity_checks(make_sd_rep(), good_rep)))
  expect_message(SPoRC::post_optim_sanity_checks(make_sd_rep(), good_rep),
                 "Successfully passed")
})

test_that("a non-finite joint likelihood fails the check", {
  for(bad in list(Inf, -Inf, NA_real_, NaN)) {
    expect_false(suppressMessages(
      SPoRC::post_optim_sanity_checks(make_sd_rep(), list(jnLL = bad))))
  } # end bad loop
  expect_message(SPoRC::post_optim_sanity_checks(make_sd_rep(), list(jnLL = Inf)),
                 "joint log-likelihood")
})

test_that("a gradient above tolerance fails and names the parameter", {
  sd_rep <- make_sd_rep(gradient = c(1e-6, 0.5))
  expect_false(suppressMessages(SPoRC::post_optim_sanity_checks(sd_rep, good_rep)))
  # The message must identify the worst parameter, not just report a failure.
  expect_message(SPoRC::post_optim_sanity_checks(sd_rep, good_rep), "ln_q")
  expect_message(SPoRC::post_optim_sanity_checks(sd_rep, good_rep), "non-convergence")
})

test_that("the gradient tolerance is configurable", {
  sd_rep <- make_sd_rep(gradient = c(1e-6, 0.5))
  # The same fit passes once the tolerance is loosened past the max gradient.
  expect_true(suppressMessages(
    SPoRC::post_optim_sanity_checks(sd_rep, good_rep, gradient_tol = 1)))
  expect_false(suppressMessages(
    SPoRC::post_optim_sanity_checks(sd_rep, good_rep, gradient_tol = 1e-8)))
})

test_that("a non-positive-definite Hessian fails the check", {
  sd_rep <- make_sd_rep(pdHess = FALSE)
  expect_false(suppressMessages(SPoRC::post_optim_sanity_checks(sd_rep, good_rep)))
  expect_message(SPoRC::post_optim_sanity_checks(sd_rep, good_rep),
                 "Hessian is not positive definite")
})

test_that("non-finite standard errors fail the check", {
  # A negative variance yields NaN under sqrt(), which is the signature of a
  # covariance matrix that could not be inverted properly.
  sd_rep <- make_sd_rep(cov = diag(c(0.04, -1)))
  expect_false(suppressMessages(
    suppressWarnings(SPoRC::post_optim_sanity_checks(sd_rep, good_rep))))
  expect_message(suppressWarnings(SPoRC::post_optim_sanity_checks(sd_rep, good_rep)),
                 "non finite elements in standard errors")
})

test_that("an oversized standard error fails and names the parameter", {
  sd_rep <- make_sd_rep(cov = diag(c(0.04, 1e6)))
  expect_false(suppressMessages(SPoRC::post_optim_sanity_checks(sd_rep, good_rep)))
  expect_message(SPoRC::post_optim_sanity_checks(sd_rep, good_rep), "ln_q")
  expect_message(SPoRC::post_optim_sanity_checks(sd_rep, good_rep), "standard error")
})

test_that("the standard error tolerance is configurable", {
  sd_rep <- make_sd_rep(cov = diag(c(0.04, 1e6)))
  expect_true(suppressMessages(
    SPoRC::post_optim_sanity_checks(sd_rep, good_rep, se_tol = 1e5)))
})

test_that("a near-perfect parameter correlation fails and names both parameters", {
  cov <- matrix(c(1, 0.995, 0.995, 1), nrow = 2)
  sd_rep <- make_sd_rep(cov = cov)
  expect_false(suppressMessages(SPoRC::post_optim_sanity_checks(sd_rep, good_rep)))
  msgs <- capture_messages(SPoRC::post_optim_sanity_checks(sd_rep, good_rep))
  expect_true(any(grepl("Parameter pairs", msgs)))
  expect_true(any(grepl("ln_M", msgs) & grepl("ln_q", msgs)))
})

test_that("the correlation tolerance is configurable", {
  cov <- matrix(c(1, 0.995, 0.995, 1), nrow = 2)
  sd_rep <- make_sd_rep(cov = cov)
  expect_true(suppressMessages(
    SPoRC::post_optim_sanity_checks(sd_rep, good_rep, corr_tol = 0.999)))
  expect_false(suppressMessages(
    SPoRC::post_optim_sanity_checks(sd_rep, good_rep, corr_tol = 0.5)))
})

test_that("a strongly negative correlation is caught as well", {
  # The check is on absolute correlation, so sign must not matter.
  cov <- matrix(c(1, -0.995, -0.995, 1), nrow = 2)
  expect_false(suppressMessages(
    SPoRC::post_optim_sanity_checks(make_sd_rep(cov = cov), good_rep)))
})

test_that("non-finite standard errors skip the correlation branch", {
  # cov2cor() on a matrix with a negative diagonal would itself be meaningless,
  # so the SE check must short-circuit before the correlation check runs.
  sd_rep <- make_sd_rep(cov = diag(c(0.04, -1)))
  msgs <- suppressWarnings(capture_messages(
    SPoRC::post_optim_sanity_checks(sd_rep, good_rep)))
  expect_false(any(grepl("Parameter pairs", msgs)))
})

test_that("several simultaneous failures are all reported", {
  sd_rep <- make_sd_rep(gradient = c(1e-6, 0.5), pdHess = FALSE)
  msgs <- capture_messages(
    SPoRC::post_optim_sanity_checks(sd_rep, list(jnLL = Inf)))
  expect_true(any(grepl("joint log-likelihood", msgs)))
  expect_true(any(grepl("absolute gradient", msgs)))
  expect_true(any(grepl("Hessian", msgs)))
  expect_false(any(grepl("Successfully passed", msgs)))
})

# ---------------------------------------------------------------------------
# marg_AIC
# ---------------------------------------------------------------------------

test_that("marg_AIC reads the objective from nlminb output", {
  opt <- list(par = c(1, 2, 3), objective = 100)
  expect_equal(SPoRC::marg_AIC(opt), 2 * 3 + 2 * 100)
})

test_that("marg_AIC reads the objective from optim output", {
  opt <- list(par = c(1, 2, 3), value = 100)
  expect_equal(SPoRC::marg_AIC(opt), 2 * 3 + 2 * 100)
})

test_that("marg_AIC applies the small-sample correction for finite n", {
  opt <- list(par = c(1, 2, 3), objective = 100)
  k <- 3
  n <- 50
  expect_equal(SPoRC::marg_AIC(opt, n = n),
               2 * k + 2 * 100 + 2 * k * (k + 1) / (n - k - 1))
  # A finite n must penalize more than the uncorrected form.
  expect_gt(SPoRC::marg_AIC(opt, n = n), SPoRC::marg_AIC(opt))
})

test_that("marg_AIC respects a custom penalty multiplier", {
  opt <- list(par = c(1, 2, 3), objective = 100)
  expect_equal(SPoRC::marg_AIC(opt, p = 4), 4 * 3 + 2 * 100)
})

test_that("marg_AIC scales the penalty with the number of parameters", {
  small <- list(par = rep(0, 2), objective = 100)
  large <- list(par = rep(0, 10), objective = 100)
  expect_equal(SPoRC::marg_AIC(large) - SPoRC::marg_AIC(small), 2 * (10 - 2))
})

# ---------------------------------------------------------------------------
# get_optim_param_list
# ---------------------------------------------------------------------------

test_that("get_optim_param_list fills a mapped vector parameter", {
  parameters <- list(ln_M = c(0, 0))
  mapping <- list(ln_M = factor(c(1, 2)))
  sd_rep <- list(par.fixed = c(ln_M = 0.1, ln_M = 0.2), par.random = NULL)

  out <- SPoRC:::get_optim_param_list(parameters, mapping, sd_rep, random = NULL)
  expect_equal(unname(out$ln_M), c(0.1, 0.2))
})

test_that("get_optim_param_list fills an unmapped array parameter", {
  parameters <- list(ln_q = matrix(0, nrow = 2, ncol = 2))
  sd_rep <- list(par.fixed = c(ln_q = 1, ln_q = 2, ln_q = 3, ln_q = 4),
                 par.random = NULL)

  out <- SPoRC:::get_optim_param_list(parameters, list(), sd_rep, random = NULL)
  expect_equal(dim(out$ln_q), c(2, 2))
  expect_equal(unname(as.vector(out$ln_q)), c(1, 2, 3, 4))
})

test_that("get_optim_param_list leaves NA-mapped elements at their starting values", {
  # An NA map entry means the parameter was fixed, so it must keep its start.
  parameters <- list(ln_M = c(-99, 0))
  mapping <- list(ln_M = factor(c(NA, 1), exclude = NULL))
  mapping$ln_M <- factor(c(NA, 1))
  sd_rep <- list(par.fixed = c(ln_M = 0.7), par.random = NULL)

  out <- SPoRC:::get_optim_param_list(parameters, mapping, sd_rep, random = NULL)
  expect_equal(unname(out$ln_M[1]), -99)
  expect_equal(unname(out$ln_M[2]), 0.7)
})

test_that("get_optim_param_list shares one estimate across mapped-together elements", {
  # Two elements mapped to the same level are a single estimated parameter.
  parameters <- list(ln_M = c(0, 0, 0))
  mapping <- list(ln_M = factor(c(1, 1, 2)))
  sd_rep <- list(par.fixed = c(ln_M = 0.3, ln_M = 0.8), par.random = NULL)

  out <- SPoRC:::get_optim_param_list(parameters, mapping, sd_rep, random = NULL)
  expect_equal(unname(out$ln_M), c(0.3, 0.3, 0.8))
})

test_that("get_optim_param_list sources random effects from par.random", {
  parameters <- list(ln_rec_devs = c(0, 0))
  mapping <- list(ln_rec_devs = factor(c(1, 2)))
  sd_rep <- list(par.fixed = c(ln_M = 0.1),
                 par.random = c(ln_rec_devs = -0.4, ln_rec_devs = 0.6))

  out <- SPoRC:::get_optim_param_list(parameters, mapping, sd_rep,
                                      random = "ln_rec_devs")
  expect_equal(unname(out$ln_rec_devs), c(-0.4, 0.6))
})

test_that("get_optim_param_list leaves parameters absent from the fit untouched", {
  parameters <- list(ln_M = c(0, 0), ln_unused = c(5, 5))
  mapping <- list(ln_M = factor(c(1, 2)))
  sd_rep <- list(par.fixed = c(ln_M = 0.1, ln_M = 0.2), par.random = NULL)

  out <- SPoRC:::get_optim_param_list(parameters, mapping, sd_rep, random = NULL)
  expect_equal(out$ln_unused, c(5, 5))
})
