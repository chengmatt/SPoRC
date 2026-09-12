library(SPoRC)
library(testthat)

# The self test has to hand its at-age settings to the operating model and then hand the draws
# back to the refit. With either half missing the operating model draws nothing, the reset blocks
# copy that nothing over the observations, and every replicate refits a block of zeros.

# a model fitting catch at age on one fishery fleet alongside an aggregated survey index,
# over whichever observed ages use_aa flags
at_age_selftest_model <- function(n_yrs, n_ages, sigma_caa, use_aa = NULL) {

  aa_dim <- c(1, n_yrs, 1, n_ages, 1, 1)
  sigma_dim <- c(n_ages, 1, 1)
  if(is.null(use_aa)) use_aa <- array(1, dim = aa_dim)

  build <- function(obs_aa) {
    build_at_age(
      n_yrs = n_yrs,
      n_ages = n_ages,
      ObsCatchAA = obs_aa,
      UseCatchAA = use_aa,
      sigmaCAA_spec = "fix",
      sigmaCAA_key = array(1L, dim = sigma_dim),
      ln_sigmaCAA = array(log(sigma_caa), dim = sigma_dim)
    )
  }

  # the observations are seeded from the model's own prediction, so the fit the self test
  # is conditioned on sits at a sensible point rather than chasing an arbitrary series
  il <- build(array(1e3, dim = aa_dim))
  return(build(array(at_age_rep(il)$PredCatchAA, dim = aa_dim)))
}


test_that("a self test on a catch-at-age model refits the operating model's draws", {

  n_yrs <- 15; n_ages <- 5
  use_aa <- array(1, dim = c(1, n_yrs, 1, n_ages, 1, 1))
  il <- at_age_selftest_model(n_yrs, n_ages, 0.2)

  fit <- fit_model(il$data, il$par, il$map, newton_loops = 3, silent = TRUE)
  sd_rep <- RTMB::sdreport(fit)

  sim_file <- tempfile(fileext = ".RDS")
  set.seed(4021)
  res <- simulation_self_test(
    data = fit$data,
    parameters = il$par,
    mapping = il$map,
    random = NULL,
    rep = fit$rep,
    sd_rep = sd_rep,
    n_sims = 1,
    output_path = sim_file,
    what = c("SSB", "CatchAA_nLL")
  )

  # the operating model drew catch at age wherever the fit observes it, and carries the
  # discards at age it drew rather than dropping them on the way out
  om <- readRDS(sim_file)
  expect_equal(dim(om$ObsCatchAA), c(1, n_yrs, 1, n_ages, 1, 1, 1))
  expect_true(all(as.numeric(om$ObsCatchAA)[as.numeric(use_aa) == 1] > 0))
  expect_equal(dim(om$TrueDiscardAA), dim(om$TrueCatchAA))
  expect_equal(dim(om$ObsDiscardAA), dim(om$ObsCatchAA))

  # and the refit fit those draws. an undrawn observation is zero, which the lognormal
  # takes the log of, so a finite likelihood over fitted cells is only reachable from
  # observations that arrived, and a replicate is a fresh draw rather than the input data
  caa_nLL <- as.numeric(res$CatchAA_nLL)
  expect_true(all(is.finite(caa_nLL)))
  expect_true(any(caa_nLL != 0))
  expect_false(isTRUE(all.equal(caa_nLL, as.numeric(fit$rep$CatchAA_nLL))))
  expect_true(all(is.finite(as.numeric(res$SSB))))
})


test_that("a self test on catch at age recovers spawning biomass without median bias", {

  n_yrs <- 30; n_ages <- 8; n_sims <- 20
  il <- at_age_selftest_model(n_yrs, n_ages, 0.2)

  fit <- fit_model(il$data, il$par, il$map, newton_loops = 3, silent = TRUE)
  expect_lt(max(abs(fit$gr(fit$env$last.par.best))), 1e-3)
  sd_rep <- RTMB::sdreport(fit)

  set.seed(918)
  res <- simulation_self_test(
    data = fit$data,
    parameters = il$par,
    mapping = il$map,
    random = NULL,
    rep = fit$rep,
    sd_rep = sd_rep,
    n_sims = n_sims,
    what = "SSB"
  )

  # every replicate refits, and the series they return is centered on the one they came from
  truth <- as.numeric(fit$rep$SSB)
  est <- matrix(res$SSB, nrow = length(truth), ncol = n_sims)
  expect_equal(sum(is.na(est)), 0)
  rel_err <- sweep(est, 1, truth, "-") / matrix(truth, length(truth), n_sims)
  expect_lt(abs(stats::median(rel_err)), 0.05)
  expect_lt(stats::median(abs(rel_err)), 0.15)
})


test_that("an age the fit leaves out is left out of the draws as well", {

  # the at-age sources have no bin argument because the use flags already run over the
  # observed ages, so this is what a bin restriction is on those streams
  n_yrs <- 15; n_ages <- 5
  use_aa <- array(1, dim = c(1, n_yrs, 1, n_ages, 1, 1))
  use_aa[, , , 1, , ] <- 0  # age one is never aged
  il <- at_age_selftest_model(n_yrs, n_ages, 0.2, use_aa = use_aa)

  fit <- fit_model(il$data, il$par, il$map, newton_loops = 3, silent = TRUE)
  sd_rep <- RTMB::sdreport(fit)

  sim_file <- tempfile(fileext = ".RDS")
  set.seed(77)
  res <- simulation_self_test(
    data = fit$data,
    parameters = il$par,
    mapping = il$map,
    random = NULL,
    rep = fit$rep,
    sd_rep = sd_rep,
    n_sims = 1,
    output_path = sim_file,
    what = c("SSB", "CatchAA_nLL")
  )

  om <- readRDS(sim_file)
  expect_true(all(om$ObsCatchAA[, , , 1, , , ] == 0))    # dropped age never drawn
  expect_true(all(om$ObsCatchAA[, , , -1, , , ] > 0))    # the rest still are
  expect_true(all(as.numeric(res$CatchAA_nLL)[as.numeric(use_aa) == 0] == 0))
  expect_true(all(is.finite(as.numeric(res$CatchAA_nLL))))
  expect_true(all(is.finite(as.numeric(res$SSB))))
})
