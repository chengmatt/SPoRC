# The operating model runs on the exact matrix exponential and the estimating model on the
# implicit solve, which approximates it rather than reaching the same numbers faster.
#
# The question is what the assessment gets wrong, and how fast that goes away with substeps.
#
# The operating model is the one from test-integration_move_timing_diffusion.R, built by
# helper-spatial_ctmc_om.R: everything the estimating model needs is known exactly, the generator is
# passed to the simulator directly, and the diffusion rate is the only free parameter.

library(SPoRC)
library(testthat)
library(Matrix)

# Test setups. The operating models are built once; only the estimation model varies. ----

oms <- list("2" = build_om(2), "0" = build_om(0))

# Report SSB at the generating parameters, with no estimation involved. This isolates the
# forward distortion: how far the implicit operator moves the population itself.
ssb_rel_err <- function(om, move_timing, nsub) {
  em <- build_em(om, move_timing, em_expm_nsub = nsub)
  pinned <- pin_at_truth(em, TRUE_LOG_THETA)
  obj <- fit_model(
    pinned$data,
    pinned$par,
    pinned$map,
    random = NULL,
    silent = TRUE,
    do_optim = FALSE
  )
  truth <- as.vector(om$SSB[, , , 1])
  max(abs(as.vector(obj$report(obj$par)$SSB) - truth) / truth)
}

# 1. Forward distortion of the population, at the generating parameters ------

test_that("the implicit solve leaves SSB alone at nsub = 0 and converges back to it", {
  # no substeps is the exact exponential, so this is the agreement the move_timing test asserts.
  # one substep is a plain solve and is badly off, checked here rather than left unstated
  err <- vapply(c(0, 1, 8, 512), function(n) ssb_rel_err(oms[["2"]], 2, n), numeric(1))

  expect_lt(err[1], 1e-3)                       # exact
  expect_gt(err[2], 0.05)                       # one backward Euler step is a different model
  expect_true(all(diff(err[-1]) < 0))           # more substeps, less distortion
  expect_lt(err[4], 5e-3)                       # and it converges back
})

test_that("under move_timing 0 only the movement fractions are approximated", {
  # the discrete timings never exponentiate the generator net of mortality, so the error is far
  # smaller than at the same substeps under continuous movement, where it compounds
  err_mt0 <- ssb_rel_err(oms[["0"]], 0, 1)
  err_mt2 <- ssb_rel_err(oms[["2"]], 2, 1)

  expect_lt(err_mt0, 0.05)
  expect_gt(err_mt2, 5 * err_mt0)
})

# 2. Bias in the estimated movement rate -------------------------------------

test_that("an implicit estimation model inflates diffusion, and substeps remove the bias", {
  # the implicit solve spreads fish further per step at the same rate, so the likelihood answers
  # by raising the rate, and that direction matters: a rate biased high is not a safe error
  fit_theta <- function(nsub) {
    em <- build_em(oms[["2"]], 2, em_expm_nsub = nsub)
    pinned <- pin_at_truth(em, START_LOG_THETA)
    fit <- fit_model(pinned$data, pinned$par, pinned$map, random = NULL, silent = TRUE)
    best <- fit$env$last.par.best
    expect_lt(max(abs(fit$gr(best))), 1e-3)
    exp(as.numeric(best[1]))
  }

  truth <- exp(TRUE_LOG_THETA)
  theta_exact <- fit_theta(0)
  theta_be1 <- fit_theta(1)
  theta_be64 <- fit_theta(64)

  # started at log(0.10), so recovery at nsub = 0 is a real recovery
  expect_equal(theta_exact, truth, tolerance = 0.02)

  # one backward Euler step costs about 12% on the movement rate, biased high
  expect_gt(theta_be1 / truth, 1.05)

  # 64 substeps put the estimate back within the exact model's own recovery error
  expect_lt(abs(theta_be64 - truth), 3 * abs(theta_exact - truth) + 0.005)
  expect_lt(abs(theta_be64 - truth), abs(theta_be1 - truth) / 5)
})
