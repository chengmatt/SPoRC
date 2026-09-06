library(SPoRC)
library(testthat)

# Two settings added together: a stock-recruit curve fitted as a penalty on the
# recruitment residual rather than generating recruitment, and a standardization
# window for non-parametric log-scale selectivity. Both default to the previous
# behavior, so the defaults are pinned alongside the new paths.

test_that("get_sr_penalty is a normal density on the log residual over the named years", {

  set.seed(1)
  Rec <- array(exp(rnorm(12)), dim = c(1, 1, 12))
  SR_pred <- array(exp(rnorm(12)), dim = c(1, 1, 12))
  sigma <- 0.4

  all_yrs <- SPoRC:::get_sr_penalty(Rec, SR_pred, sigma)
  expect_equal(as.vector(all_yrs),
               -dnorm(log(as.vector(Rec)) - log(as.vector(SR_pred)), 0, sigma, log = TRUE))

  # Years outside the window contribute nothing at all, rather than contributing
  # a residual that happens to be small.
  win <- 3:9
  part <- SPoRC:::get_sr_penalty(Rec, SR_pred, sigma, yrs = win)
  expect_equal(as.vector(part)[-win], rep(0, 12 - length(win)))
  expect_equal(as.vector(part)[win], as.vector(all_yrs)[win])

})

test_that("Setup_Mod_Rec encodes sr_R0_spec and maps ln_sr_R0 to match", {

  base <- Setup_Mod_Dim(
    years = 1:20,
    ages = 1:10,
    lens = NA,
    n_regions = 1,
    n_sexes = 1,
    n_fish_fleets = 1,
    n_srv_fleets = 1,
    n_seas = 1,
    n_pop = 1,
    natal_region = 1,
    verbose = FALSE
  )
  rec <- function(...) suppressWarnings(Setup_Mod_Rec(
    input_list = base,
    rec_model = "mean_rec",
    rec_lag = 1,
    SR_ref_yr = 1,
    init_age_strc = 1,
    equil_init_age_strc = 2,
    t_spawn = 0,
    ...
  ))

  # ln_sr_R0 is only ever in the parameter vector when the curve owns a scale of
  # its own; under the other settings it is mapped off so it cannot move.
  shared <- rec(sr_penalty = "bh", sr_pen_sigma = 0.5, sr_R0_spec = "shared", use_rinit = 0)
  expect_equal(shared$data$sr_R0_spec, 0)
  expect_true(all(is.na(shared$map$ln_sr_R0)))

  est <- rec(sr_penalty = "bh", sr_pen_sigma = 0.5, sr_R0_spec = "est", use_rinit = 0)
  expect_equal(est$data$sr_R0_spec, 1)
  expect_false(any(is.na(est$map$ln_sr_R0)))

  rinit <- rec(sr_penalty = "bh", sr_pen_sigma = 0.5, sr_R0_spec = "rinit", use_rinit = 1)
  expect_equal(rinit$data$sr_R0_spec, 2)
  expect_true(all(is.na(rinit$map$ln_sr_R0)))

  # No curve at all leaves the penalty off and steepness unreachable.
  none <- rec(use_rinit = 0)
  expect_equal(none$data$sr_penalty, 0)
  expect_true(all(is.na(none$map$steepness_h)))

})

test_that("Setup_Mod_Rec rejects a penalty that would be counted twice or has no scale", {

  base <- Setup_Mod_Dim(
    years = 1:20,
    ages = 1:10,
    lens = NA,
    n_regions = 1,
    n_sexes = 1,
    n_fish_fleets = 1,
    n_srv_fleets = 1,
    n_seas = 1,
    n_pop = 1,
    natal_region = 1,
    verbose = FALSE
  )

  # Under bh_rec the deviation already IS the residual, so penalizing it as well
  # would penalize the same quantity twice.
  expect_error(
    suppressWarnings(Setup_Mod_Rec(
      base,
      rec_model = "bh_rec",
      rec_lag = 1,
      SR_ref_yr = 1,
      sr_penalty = "bh",
      sr_pen_sigma = 0.5,
      init_age_strc = 1,
      equil_init_age_strc = 2,
      t_spawn = 0,
      use_rinit = 0
    )),
    "only valid with rec_model")

  # "rinit" takes the curve's scale from a parameter that use_rinit = 0 leaves out.
  expect_error(
    suppressWarnings(Setup_Mod_Rec(
      base,
      rec_model = "mean_rec",
      rec_lag = 1,
      SR_ref_yr = 1,
      sr_penalty = "bh",
      sr_pen_sigma = 0.5,
      sr_R0_spec = "rinit",
      init_age_strc = 1,
      equil_init_age_strc = 2,
      t_spawn = 0,
      use_rinit = 0
    )),
    "needs use_rinit = 1")

  expect_error(
    suppressWarnings(Setup_Mod_Rec(
      base,
      rec_model = "mean_rec",
      rec_lag = 1,
      SR_ref_yr = 1,
      sr_penalty = "bh",
      sr_pen_sigma = 0.5,
      sr_R0_spec = "average",
      init_age_strc = 1,
      equil_init_age_strc = 2,
      t_spawn = 0,
      use_rinit = 0
    )),
    "shared, est, or rinit")

})

test_that("sel_norm_bins sets which bins the nonparlog standardization averages over", {

  pars <- c(-1.2, -0.4, 0.1, 0.5, 0.7, 0.75, 0.7, 0.6, 0.4, 0.2)
  bins <- seq_along(pars)
  sel <- function(nb) SPoRC:::Get_Selex(
    Selex_Model = 9,
    TimeVary_Model = 0,
    pars = pars,
    ln_seldevs = NULL,
    Region = 1,
    Year = 1,
    Bin = bins,
    Sex = 1,
    sel_norm_bins = nb
  )

  # Not supplying a window averages over every bin, which is the behavior the
  # argument replaced.
  expect_equal(as.vector(sel(NULL)), exp(pars) / mean(exp(pars)))
  expect_equal(mean(as.vector(sel(NULL))), 1)
  expect_equal(as.vector(sel(integer(0))), as.vector(sel(NULL)))
  expect_equal(as.vector(sel(bins)), as.vector(sel(NULL)))

  # A window makes the mean one over that window, not over the whole range, and
  # every bin is rescaled by the same constant.
  win <- 4:8
  s <- as.vector(sel(win))
  expect_equal(mean(s[win]), 1)
  expect_false(isTRUE(all.equal(mean(s), 1)))
  expect_equal(s / as.vector(sel(NULL)), rep(s[1] / as.vector(sel(NULL))[1], length(bins)))

})

test_that("sel_norm_bins defaults to every bin through the setup functions", {

  il <- Setup_Mod_Dim(
    years = 1:5,
    ages = 1:8,
    lens = NA,
    n_regions = 1,
    n_sexes = 1,
    n_fish_fleets = 1,
    n_srv_fleets = 1,
    n_seas = 1,
    n_pop = 1,
    natal_region = 1,
    verbose = FALSE
  )
  il <- SPoRC:::setup_sel_norm_bins(il, NULL, prefix = "srv", n_fleets = 1, bins = 8)
  expect_equal(dim(il$data$srv_sel_norm_bins), c(8L, 1L))
  expect_true(all(il$data$srv_sel_norm_bins == 1))

  il <- SPoRC:::setup_sel_norm_bins(il, list(3:6), prefix = "fish", n_fleets = 1, bins = 8)
  expect_equal(which(il$data$fish_sel_norm_bins[, 1] == 1), 3:6)

  expect_error(SPoRC:::setup_sel_norm_bins(il, list(3:20), prefix = "fish", n_fleets = 1, bins = 8),
               "outside 1:8")
  expect_error(SPoRC:::setup_sel_norm_bins(il, list(integer(0)), prefix = "fish", n_fleets = 1, bins = 8),
               "at least one bin")
  expect_error(SPoRC:::setup_sel_norm_bins(il, list(1:3, 1:3), prefix = "fish", n_fleets = 1, bins = 8),
               "one element per fleet")

})
