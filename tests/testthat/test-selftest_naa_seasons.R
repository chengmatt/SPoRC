library(SPoRC)
library(testthat)

# Operating-model and estimation-model agreement for the seasonal state-space numbers at age.
# The simulator draws innovations at every active season boundary and the penalty is the density of
# that same covariance, so the two have to agree on which dim is the season and on what
# sigmaNAA means. A disagreement here produces a plausible number rather than an error.

test_that("a simulated season correlation appears across seasons and leaves ages and years alone", {
  set.seed(31)
  nk <- 3; ny <- 40; na <- 10
  env <- naaseas_env(
    n_yrs = ny,
    n_ages = na,
    n_seas = nk,
    NAA_re = 1,
    sigmaNAA = 0.4,
    NAA_re_season = 1,
    naa_season_corr = c(0.7, 0.3, 0.5)
  )
  D <- replicate(600, SPoRC:::draw_naa_innovations(env), simplify = FALSE)

  kc <- function(i, j) mean(vapply(D, function(x)
    stats::cor(as.vector(x[1,1,,i,,1]), as.vector(x[1,1,,j,,1])), numeric(1)))
  expect_equal(kc(1, 2), 0.7, tolerance = 0.03)
  expect_equal(kc(1, 3), 0.3, tolerance = 0.03)
  expect_equal(kc(2, 3), 0.5, tolerance = 0.03)

  # the correlation is on the season dim and nowhere else
  cor_age <- function(x) stats::cor(as.vector(x[1,1,,1,1:(na-1),1]), as.vector(x[1,1,,1,2:na,1]))
  cor_yr <- function(x) stats::cor(as.vector(x[1,1,1:(ny-1),1,,1]), as.vector(x[1,1,2:ny,1,,1]))
  expect_equal(mean(vapply(D, cor_age, numeric(1))), 0, tolerance = 0.03)
  expect_equal(mean(vapply(D, cor_yr, numeric(1))), 0, tolerance = 0.03)
  expect_equal(mean(vapply(D, function(x) stats::sd(as.vector(x)), numeric(1))), 0.4, tolerance = 0.02)
})

test_that("the simulator draws only in the seasons the state is active over", {
  set.seed(32)
  env <- naaseas_env(n_yrs = 20, n_ages = 6, n_seas = 3, naa_re_seas = c(1L, 3L))
  eta <- SPoRC:::draw_naa_innovations(env)
  expect_equal(dim(eta), c(1, 1, 20, 3, 6, 1))
  expect_equal(max(abs(eta[,,,2,,])), 0)
  expect_gt(stats::sd(as.vector(eta[,,,c(1, 3),,])), 0)
})

test_that("the penalty recovers the seasonal process the simulator drew from", {
  # Given the true states, maximizing the penalty must return the process parameters. This is what
  # catches a simulator and a penalty that disagree about which dim the season is, since a
  # swapped dim still returns a finite optimum, just the wrong one.
  set.seed(33)
  nk <- 3; ny <- 40; na <- 10
  env <- naaseas_env(
    n_yrs = ny,
    n_ages = na,
    n_seas = nk,
    NAA_re = 1,
    sigmaNAA = 0.4,
    NAA_re_season = 1,
    naa_season_corr = c(0.7, 0.3, 0.5)
  )
  pred <- array(exp(5), dim = c(1, 1, ny, nk, na, 1))

  est <- t(vapply(1:12, function(i) {
    eta <- SPoRC:::draw_naa_innovations(env)
    nll <- function(th) SPoRC:::Get_NAA_state_penalty(
      log(pred) + eta,
      pred,
      array(exp(th[1]), dim = dim(pred)),
      1:na,
      1:ny,
      1:nk,
      NAA_re = 1,
      NAA_pe_pars = array(0, dim = c(1, 1, 3, 1)),
      NAA_re_season = 1,
      NAA_season_corr_pars = array(th[2:4], dim = c(1, 3, 1))
    )
    f <- stats::nlminb(c(log(0.25), rep(0, 3)), nll)
    C <- SPoRC:::build_us_corr(f$par[2:4], nk)
    c(exp(f$par[1]), C[2,1], C[3,1], C[3,2])
  }, numeric(4)))

  m <- colMeans(est)
  expect_equal(m[1], 0.4, tolerance = 0.03)
  expect_equal(m[2], 0.7, tolerance = 0.06)
  expect_equal(m[3], 0.3, tolerance = 0.06)
  expect_equal(m[4], 0.5, tolerance = 0.06)
})

test_that("the operating model applies the state at every active season boundary", {
  om_all <- naaseas_make_om(NAA_re = "iid", sigmaNAA = 0.35, NAA_re_seasons = "all", seed = 707)
  om_ann <- naaseas_make_om(NAA_re = "iid", sigmaNAA = 0.35, NAA_re_seasons = "annual", seed = 707)

  expect_equal(dim(om_all$naa_eta)[4], naaseas_cfg$n_seas)
  expect_gt(stats::sd(as.vector(om_all$naa_eta[,,,2,,,])), 0)
  # the annual operating model leaves the within-year step deterministic
  expect_equal(max(abs(om_ann$naa_eta[,,,2,,,])), 0)
  expect_gt(stats::sd(as.vector(om_ann$naa_eta[,,,1,,,])), 0)
})

test_that("the estimation model recovers a seasonal process error from data alone", {
  # Nothing about the states is known here: the variance has to be separated from observation error
  # with the states integrated out, which is the part of a state-space fit that fails in practice.
  # Both seasons are fished and surveyed, so both have information about their own state.
  sd_true <- 0.35
  out <- t(vapply(1:4, function(i) {
    om <- naaseas_make_om(NAA_re = "iid", sigmaNAA = sd_true, NAA_re_seasons = "all", seed = 700 + i)
    il <- naaseas_build_em(naaseas_om_data(om), NAA_re = "iid", NAA_re_seasons = "all")
    fit <- suppressWarnings(fit_model(il$data, il$par, il$map, random = "ln_NAA", silent = TRUE))
    p <- fit$env$last.par.best
    c(exp(unname(p[-fit$env$random][names(p[-fit$env$random]) == "ln_sigmaNAA"])),
      naaseas_ssb_rmse(fit, om),
      length(fit$env$random) == il$data$n_est_naa_re)
  }, numeric(3)))

  expect_true(all(out[, 3] == 1))                                  # every state is a random effect
  expect_equal(stats::median(out[, 1]), sd_true, tolerance = 0.06) # process error recovered
  expect_lt(stats::median(out[, 2]), 0.15)                         # and the biomass tracks
})

test_that("the annual state on a seasonal model estimates one boundary per year", {
  om <- naaseas_make_om(NAA_re = "iid", sigmaNAA = 0.35, NAA_re_seasons = "annual", seed = 707)
  dat <- naaseas_om_data(om)
  ann <- naaseas_build_em(dat, NAA_re = "iid")
  all_k <- naaseas_build_em(dat, NAA_re = "iid", NAA_re_seasons = "all")

  expect_equal(ann$data$naa_re_seas, 1L)
  expect_equal(all_k$data$n_est_naa_re, naaseas_cfg$n_seas * ann$data$n_est_naa_re)

  fit <- suppressWarnings(fit_model(ann$data, ann$par, ann$map, random = "ln_NAA", silent = TRUE))
  expect_equal(length(fit$env$random), ann$data$n_est_naa_re)
  p <- fit$env$last.par.best
  sd_hat <- exp(unname(p[-fit$env$random][names(p[-fit$env$random]) == "ln_sigmaNAA"]))
  expect_equal(sd_hat, 0.35, tolerance = 0.08)
  expect_lt(naaseas_ssb_rmse(fit, om), 0.15)
})
