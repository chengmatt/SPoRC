library(SPoRC)
library(testthat)

# Minimal 1-pop, 1-region, 1-year, 1-season, 2-age, 1-sex, 1-fleet fixture.
# Dimensions match the real model's convention even when degenerate (extent-1
# dims), so this exercises the same array-indexing/dropping behavior as a
# real (if tiny) model run.
make_fishery_obs_input <- function(catch_units = 0, discard_units = 0, fish_idx_type = 0, fit_lengths = 0) {

  n_pop <- 1; n_regions <- 1; n_yrs <- 1; n_seas <- 1; n_fish_fleets <- 1; n_sexes <- 1; n_ages <- 2; n_lens <- 2

  ZAA <- array(c(0.5, 0.3), dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes))
  NAA <- array(c(100, 50), dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes))
  ret_FAA <- array(c(0.2, 0.1), dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_fish_fleets))
  disc_FAA <- array(c(0.05, 0.02), dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_fish_fleets))
  WAA_fish <- array(c(1.5, 3.0), dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_fish_fleets))
  fish_sel <- array(1, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_fish_fleets))
  ret_sel <- array(1, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_fish_fleets))
  dmr <- array(1, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))

  list(
    n_pop = n_pop, n_regions = n_regions, n_yrs = n_yrs, n_seas = n_seas,
    n_fish_fleets = n_fish_fleets, n_sexes = n_sexes,
    fish_q_blocks = array(1L, dim = c(n_regions, n_yrs, n_fish_fleets)),
    ln_fish_q = array(log(0.001), dim = c(n_regions, 1, n_fish_fleets)),
    fish_q = array(0, dim = c(n_regions, n_yrs, n_fish_fleets)),
    ret_FAA = ret_FAA, disc_FAA = disc_FAA, ZAA = ZAA, NAA = NAA,
    CAA = array(0, dim = dim(ret_FAA)), DAA = array(0, dim = dim(ret_FAA)),
    CAL = array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_lens, n_sexes, n_fish_fleets)),
    DAL = array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_lens, n_sexes, n_fish_fleets)),
    PredCatch = array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_fish_fleets)),
    PredDiscard = array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_fish_fleets)),
    PredFishIdx = array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_fish_fleets)),
    fit_lengths = fit_lengths,
    SizeAgeTrans = array(rep(diag(n_ages), n_pop * n_regions * n_yrs * n_seas * n_sexes),
                         dim = c(n_pop, n_regions, n_yrs, n_seas, n_lens, n_ages, n_sexes)),
    catch_units = catch_units, discard_units = discard_units, WAA_fish = WAA_fish, dmr = dmr,
    fish_idx_type = fish_idx_type, fish_sel = fish_sel, ret_sel = ret_sel
  )
}

test_that("get_fishery_observation_model implements the Baranov catch equation", {

  il <- make_fishery_obs_input()
  out <- do.call(SPoRC:::get_fishery_observation_model, il)

  ZAA_vec <- c(0.5, 0.3); NAA_vec <- c(100, 50)
  ret_FAA_vec <- c(0.2, 0.1); disc_FAA_vec <- c(0.05, 0.02)

  expected_CAA <- ret_FAA_vec / ZAA_vec * NAA_vec * (1 - exp(-ZAA_vec))
  expected_DAA <- disc_FAA_vec / ZAA_vec * NAA_vec * (1 - exp(-ZAA_vec))

  expect_equal(as.numeric(out$CAA[1,1,1,1,,1,1]), expected_CAA, tolerance = 1e-8)
  expect_equal(as.numeric(out$DAA[1,1,1,1,,1,1]), expected_DAA, tolerance = 1e-8)
})

test_that("get_fishery_observation_model: catch_units and discard_units switches", {

  il_abd <- make_fishery_obs_input(catch_units = 0, discard_units = 0)
  out_abd <- do.call(SPoRC:::get_fishery_observation_model, il_abd)
  expect_equal(out_abd$PredCatch[1,1,1,1,1], sum(out_abd$CAA[1,1,1,1,,1,1]), tolerance = 1e-8)
  expect_equal(out_abd$PredDiscard[1,1,1,1,1], sum(out_abd$DAA[1,1,1,1,,1,1]), tolerance = 1e-8) # dmr = 1

  il_biom <- make_fishery_obs_input(catch_units = 1, discard_units = 1)
  out_biom <- do.call(SPoRC:::get_fishery_observation_model, il_biom)
  WAA_vec <- c(1.5, 3.0)
  expect_equal(out_biom$PredCatch[1,1,1,1,1], sum(out_biom$CAA[1,1,1,1,,1,1] * WAA_vec), tolerance = 1e-8)
  expect_equal(out_biom$PredDiscard[1,1,1,1,1], sum(out_biom$DAA[1,1,1,1,,1,1] * WAA_vec), tolerance = 1e-8)

  # abundance/biomass fraction variants must fall strictly within (0, 1)
  il_frac <- make_fishery_obs_input(discard_units = 2)
  out_frac <- do.call(SPoRC:::get_fishery_observation_model, il_frac)
  expect_true(out_frac$PredDiscard[1,1,1,1,1] > 0 && out_frac$PredDiscard[1,1,1,1,1] < 1)

  il_frac_b <- make_fishery_obs_input(discard_units = 3)
  out_frac_b <- do.call(SPoRC:::get_fishery_observation_model, il_frac_b)
  expect_true(out_frac_b$PredDiscard[1,1,1,1,1] > 0 && out_frac_b$PredDiscard[1,1,1,1,1] < 1)
})

test_that("get_fishery_observation_model: fish_q and fishery index (abundance vs biomass)", {

  il <- make_fishery_obs_input(fish_idx_type = 0)
  out <- do.call(SPoRC:::get_fishery_observation_model, il)

  expect_equal(out$fish_q[1,1,1], 0.001, tolerance = 1e-10) # exp(log(0.001))
  # full selectivity (fish_sel = ret_sel = 1) -> index is just q * sum(NAA)
  expect_equal(out$PredFishIdx[1,1,1,1,1], 0.001 * sum(c(100, 50)), tolerance = 1e-8)

  il_b <- make_fishery_obs_input(fish_idx_type = 1)
  out_b <- do.call(SPoRC:::get_fishery_observation_model, il_b)
  expect_equal(out_b$PredFishIdx[1,1,1,1,1], 0.001 * sum(c(100, 50) * c(1.5, 3.0)), tolerance = 1e-8)
})

test_that("get_fishery_observation_model: length compositions via SizeAgeTrans (identity here)", {

  il <- make_fishery_obs_input(fit_lengths = 1)
  out <- do.call(SPoRC:::get_fishery_observation_model, il)

  # SizeAgeTrans is the identity matrix in this fixture, so CAL should equal CAA elementwise
  expect_equal(as.numeric(out$CAL[1,1,1,1,,1,1]), as.numeric(out$CAA[1,1,1,1,,1,1]), tolerance = 1e-8)
  expect_equal(as.numeric(out$DAL[1,1,1,1,,1,1]), as.numeric(out$DAA[1,1,1,1,,1,1]), tolerance = 1e-8)
})

make_survey_obs_input <- function(srv_idx_type = 0, srv_selex_type = 0, do_srv_q_cov = 0, fit_lengths = 0) {

  n_pop <- 1; n_regions <- 1; n_yrs <- 1; n_seas <- 1; n_srv_fleets <- 1; n_sexes <- 1; n_ages <- 2; n_lens <- 2; n_cov <- 1

  ZAA <- array(c(0.5, 0.3), dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes))
  NAA <- array(c(100, 50), dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes))
  WAA_srv <- array(c(1.5, 3.0), dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_srv_fleets))
  srv_sel <- array(c(0.4, 0.9), dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_srv_fleets))
  srv_sel_l <- array(1, dim = c(n_regions, n_yrs, n_lens, n_sexes, n_srv_fleets))

  list(
    n_pop = n_pop, n_regions = n_regions, n_yrs = n_yrs, n_seas = n_seas,
    n_srv_fleets = n_srv_fleets, n_sexes = n_sexes,
    srv_q_blocks = array(1L, dim = c(n_regions, n_yrs, n_srv_fleets)),
    ln_srv_q = array(log(0.002), dim = c(n_regions, 1, n_srv_fleets)),
    srv_q = array(0, dim = c(n_regions, n_yrs, n_srv_fleets)),
    do_srv_q_cov = do_srv_q_cov,
    srv_q_cov = array(0.1, dim = c(n_regions, n_yrs, n_srv_fleets, n_cov)),
    srv_q_coeff = array(0.5, dim = c(n_regions, n_srv_fleets, n_cov)),
    srv_selex_type = srv_selex_type, srv_sel = srv_sel, srv_sel_l = srv_sel_l,
    SizeAgeTrans = array(rep(diag(n_ages), n_pop * n_regions * n_yrs * n_seas * n_sexes),
                         dim = c(n_pop, n_regions, n_yrs, n_seas, n_lens, n_ages, n_sexes)),
    NAA = NAA, ZAA = ZAA,
    t_srv = array(0.5, dim = c(n_regions, n_seas, n_srv_fleets)),
    SrvIAA = array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_srv_fleets)),
    fit_lengths = fit_lengths,
    SrvIAL = array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_lens, n_sexes, n_srv_fleets)),
    srv_idx_type = srv_idx_type, WAA_srv = WAA_srv,
    PredSrvIdx = array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_srv_fleets))
  )
}

test_that("get_survey_observation_model implements the mid-season-survival index-at-age formula", {

  il <- make_survey_obs_input()
  out <- do.call(SPoRC:::get_survey_observation_model, il)

  NAA_vec <- c(100, 50); ZAA_vec <- c(0.5, 0.3); sel_vec <- c(0.4, 0.9); t_srv <- 0.5
  expected_SrvIAA <- NAA_vec * sel_vec * exp(-t_srv * ZAA_vec)

  expect_equal(as.numeric(out$SrvIAA[1,1,1,1,,1,1]), expected_SrvIAA, tolerance = 1e-8)
  expect_equal(out$srv_q[1,1,1], 0.002, tolerance = 1e-10)
  expect_equal(out$PredSrvIdx[1,1,1,1,1], 0.002 * sum(expected_SrvIAA), tolerance = 1e-8) # srv_idx_type = 0 (abundance)
})

test_that("get_survey_observation_model: biomass index and catchability covariate effect", {

  il_b <- make_survey_obs_input(srv_idx_type = 1)
  out_b <- do.call(SPoRC:::get_survey_observation_model, il_b)
  WAA_vec <- c(1.5, 3.0)
  expect_equal(out_b$PredSrvIdx[1,1,1,1,1], 0.002 * sum(out_b$SrvIAA[1,1,1,1,,1,1] * WAA_vec), tolerance = 1e-8)

  il_cov <- make_survey_obs_input(do_srv_q_cov = 1)
  out_cov <- do.call(SPoRC:::get_survey_observation_model, il_cov)
  expect_equal(out_cov$srv_q[1,1,1], 0.002 * exp(0.1 * 0.5), tolerance = 1e-10) # q * exp(cov * coeff)
})

test_that("get_survey_observation_model: length-based selectivity is converted via SizeAgeTrans", {

  il <- make_survey_obs_input(srv_selex_type = 1) # srv_sel_l is all 1s in this fixture
  out <- do.call(SPoRC:::get_survey_observation_model, il)

  # SizeAgeTrans is the identity here, so age-based srv_sel should end up all 1s too
  expect_equal(as.numeric(out$srv_sel[1,1,1,1,,1,1]), c(1, 1), tolerance = 1e-8)
})
