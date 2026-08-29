# Multi-population models, which nothing else in the suite evaluates.
#
# Every dataset that ships with the package is single-population: mlt_rg_sable
# and three_rg_sable are multi-*region* with n_pop = 1. The tests that build
# n_pop > 1 exercise truncate_yr, which reshapes data and never touches the
# objective. And the one worked multi-population example,
# vignettes/r_natal-homing-pop-lrgr-rg.Rmd, has eval = FALSE on twelve of its
# thirteen chunks including the one that fits, so it has never run in a build.
#
# So the natal-homing machinery had an entry in the API, a vignette, and no
# execution. These are the first checks that evaluate and differentiate it.

multi_pop_input <- function(n_pop, n_regions, natal_region, n_yrs = 10, n_ages = 6) {
  off <- array(0, dim = c(n_regions, n_yrs, 1, 1))
  sweep_input(
    dims = list(n_regions = n_regions, n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = 1,
                n_yrs = n_yrs, n_ages = n_ages, n_pop = n_pop, natal_region = natal_region),
    # more than one population requires local density dependence
    rec = if(n_pop > 1) list(rec_dd = "local") else list(),
    # the observation streams are off: what is under test is the process model
    # once populations are stacked
    catch = list(UseCatch = off),
    fishidx = list(UseFishIdx = off, fish_idx_type = "none",
                   UseFishAgeComps = off, FishAgeComps_LikeType = "none"),
    srvidx = list(UseSrvIdx = off, srv_idx_type = "none",
                  UseSrvAgeComps = off, SrvAgeComps_LikeType = "none"))
}


test_that("stacked identical populations multiply the likelihood", {
  # Populations that do not stray and share every rate are independent copies of
  # one stock, so the joint negative log likelihood of k of them across k regions
  # is k^2 times a single-region single-population stock. It holds only if the
  # population margin is walked correctly everywhere the objective indexes on it,
  # which makes it a collapse relation like the region and sex ones.
  one <- as.numeric(SPoRC_rtmb(multi_pop_input(1, 1, NA)$par, multi_pop_input(1, 1, NA)$data))
  expect_true(is.finite(one))

  for(k in 2:3) {
    il <- multi_pop_input(k, k, seq_len(k))
    expect_equal(as.numeric(SPoRC_rtmb(il$par, il$data)), one * k^2, tolerance = 1e-10,
                 label = sprintf("objective at n_pop = %d", k))
  }
})


test_that("a multi-population model builds an AD tape", {
  # The estimation path, which nothing exercised before. Evaluating the objective
  # on ordinary doubles is not enough: an operation that drops the tape leaves
  # the value right and the model unfittable.
  for(cfg in list(list(k = 2, nr = 1, nat = c(1, 1)),
                  list(k = 2, nr = 2, nat = c(1, 2)),
                  list(k = 3, nr = 3, nat = 1:3))) {
    il <- multi_pop_input(cfg$k, cfg$nr, cfg$nat)
    expect_no_error(fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE),
                    message = sprintf("taping at n_pop = %d, n_regions = %d", cfg$k, cfg$nr))
  }
})


test_that("a multi-population model has a finite gradient", {
  # Taping and differentiating are separate claims: a tape can build and still
  # produce a gradient full of NaN.
  il <- multi_pop_input(2, 2, c(1, 2))
  obj <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)

  expect_true(is.finite(obj$fn(obj$par)))
  expect_true(all(is.finite(obj$gr(obj$par))))
})


test_that("unfished recruitment must be given one value per population", {
  # A single starting value on a model carrying several populations is read
  # position by position further in and indexes past its own end, which reaches
  # RTMB as "not a valid advector" rather than as a problem with the argument.
  #
  # Caught by the same guard every other starting value goes through, rather than
  # by a check special to this one argument.
  expect_error(
    sweep_input(dims = list(n_regions = 1, n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = 1,
                            n_yrs = 10, n_ages = 6, n_pop = 2, natal_region = c(1, 1)),
                rec = list(rec_dd = "local", ln_global_R0 = log(1e6))),
    "starting value for ln_global_R0 is length 1 where the model expects 2")

  expect_no_error(
    sweep_input(dims = list(n_regions = 1, n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = 1,
                            n_yrs = 10, n_ages = 6, n_pop = 2, natal_region = c(1, 1)),
                rec = list(rec_dd = "local", ln_global_R0 = rep(log(1e6), 2))))
})
