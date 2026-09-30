# A joint simulation on a fit with random effects draws from sd_rep$jointPrecision, which the
# cov.fixed route in test-selftest_joint_routing.R never reaches. Checks the precision is usable,
# that the branch runs, and that every replicate gets its own parameter vector.

library(SPoRC)
library(testthat)

# the same small model the feature self-tests use, with the recruitment deviations integrated
joint_prec_fit <- function(seed = 41) {
  om <- selftest_make_om(idx_se_om = 0.05, iss_om = 1e4, seed = seed)
  sim_data <- simulation_data_to_SPoRC(sim_env = om, y = selftest_cfg$n_yrs, sim = 1)
  inp <- selftest_build_input(sim_data)
  fit <- fit_model(inp$data, inp$par, inp$map, random = "ln_RecDevs", silent = TRUE)
  list(fit = fit, inp = inp, sd_rep = RTMB::sdreport(fit, getJointPrecision = TRUE))
}

test_that("a fit with random effects gives a joint precision a draw can be taken from", {

  skip_on_cran()
  fm <- joint_prec_fit()
  prec <- fm$sd_rep$jointPrecision

  # the branch only exists when sdreport assembled one
  expect_false(is.null(prec))
  expect_true(all(is.finite(as.matrix(prec))))

  # both blocks are there, so the matrix spans the fixed effects and the deviations
  expect_equal(nrow(prec), length(fm$fit$env$last.par.best))
  expect_true(all(c("ln_RecDevs") %in% rownames(prec)))

  # rmvnorm_prec factors it, which is what failed when a fit had not converged
  expect_no_error(Matrix::Cholesky(prec, super = TRUE))
  expect_true(all(is.finite(fm$sd_rep$cov.fixed)))
  expect_true(fm$sd_rep$pdHess)

  # and the draw comes back the right shape and finite
  set.seed(1)
  draws <- rmvnorm_prec(fm$fit$env$last.par.best, prec, n_sims = 4)
  expect_equal(dim(draws), c(nrow(prec), 4))
  expect_true(all(is.finite(draws)))
})

test_that("a joint self test draws each replicate its own parameters through that branch", {

  skip_on_cran()
  fm <- joint_prec_fit()
  n_sims <- 3

  set.seed(2)
  res <- simulation_self_test(
    fm$inp$data,
    fm$inp$par,
    fm$inp$map,
    random = "ln_RecDevs",
    obj = fm$fit,
    rep = fm$fit$rep,
    sd_rep = fm$sd_rep,
    n_sims = n_sims,
    newton_loops = 0,
    do_par = FALSE,
    what = c("SSB"),
    what_par = c("ln_global_R0"),
    sim_type = "joint"
  )

  # the truth is a draw per replicate, not the fit repeated
  truth <- res$truth$SSB
  expect_equal(dim(truth)[length(dim(truth))], n_sims)
  spread <- apply(matrix(truth, ncol = n_sims), 1, function(z) diff(range(z)))
  expect_true(max(spread) > 0)

  # and the refits are finite, so the draws stayed in a region the model can be fit at
  expect_true(all(is.finite(res$SSB)))
})

test_that("a joint self test refuses a fit whose sdreport has no joint precision", {

  skip_on_cran()
  om <- selftest_make_om(idx_se_om = 0.05, iss_om = 1e4, seed = 41)
  sim_data <- simulation_data_to_SPoRC(sim_env = om, y = selftest_cfg$n_yrs, sim = 1)
  inp <- selftest_build_input(sim_data)
  fit <- fit_model(inp$data, inp$par, inp$map, random = "ln_RecDevs", silent = TRUE)

  # sdreport without getJointPrecision, so the branch has nothing to draw from
  expect_error(
    simulation_self_test(inp$data, inp$par, inp$map, random = "ln_RecDevs", obj = fit,
                         rep = fit$rep, sd_rep = RTMB::sdreport(fit), n_sims = 2,
                         newton_loops = 0, do_par = FALSE, what = c("SSB"),
                         what_par = c("ln_global_R0"), sim_type = "joint"),
    "getJointPrecision"
  )
})
