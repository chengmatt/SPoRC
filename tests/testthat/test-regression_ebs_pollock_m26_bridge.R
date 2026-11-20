# checks SPoRC against the September 2026 EBS pollock Model 26.0 (Rceattle) at its estimate, without optimizing
# the fishery cpue (mean numbers over the year here, start of year in rceattle) and the bottom trawl age 1 penalty
# differ, so the gradient is only checked on the survey selectivity and catchability
library(SPoRC)
library(testthat)
data("sgl_rg_ebswp_m26_data")

test_that("EBS pollock reproduces Rceattle Model 26.0 at its estimate", {

  inputs <- sgl_rg_ebswp_m26_data$inputs
  rce <- sgl_rg_ebswp_m26_data$rceattle
  n_yrs <- length(inputs$years)
  y_94 <- which(inputs$years == 1994)

  # setup model and fix values at the rceattle estimate
  input_list <- suppressWarnings(build_ebs_pollock_m26_input(inputs = inputs, rce = rce))
  input_list <- seed_ebs_pollock_m26_rceattle(input_list = input_list, rce = rce)

  obj <- fit_model(data = input_list$data,
                   parameters = input_list$par,
                   mapping = input_list$map,
                   random = NULL,
                   do_optim = FALSE,
                   silent = TRUE)
  sporc_rep <- obj$report(obj$par)

  # 1218 estimated parameters
  expect_equal(length(obj$par), 1218)

  # Compare dynamics
  expect_equal(sporc_rep$NAA[1,1,1:n_yrs,1,,1], rce$naa, tolerance = 1e-8)
  expect_equal(as.vector(sporc_rep$Rec), rce$rec, tolerance = 1e-8)
  expect_equal(as.vector(sporc_rep$SSB), rce$ssb, tolerance = 1e-8)

  # Compare selectivity (acoustic age 1 is not used here and zero in rceattle)
  expect_equal(sporc_rep$fish_sel[1,1,1:n_yrs,1,,1,1], rce$fish_sel, tolerance = 1e-8)
  expect_equal(sporc_rep$srv_sel[1,1,1:n_yrs,1,,1,1], rce$bts_sel, tolerance = 1e-8)
  expect_equal(sporc_rep$srv_sel[1,1,y_94:n_yrs,1,2:15,1,2], rce$ats_sel[y_94:n_yrs,2:15], tolerance = 1e-8)
  expect_equal(sporc_rep$srv_sel[,,,,,,3], sporc_rep$srv_sel[,,,,,,2]) # vessel of opportunity shares the acoustic curve

  # Compare likelihoods
  expect_equal(sum(sporc_rep$Catch_nLL), rce$catch_nll, tolerance = 1e-8)
  expect_equal(sum(sporc_rep$Rec_nLL), rce$rec_nll, tolerance = 1e-8)
  expect_equal(sum(sporc_rep$Init_Rec_nLL), rce$init_nll, tolerance = 1e-8)
  expect_equal(sum(sporc_rep$FishAgeComps_nLL) + sum(sporc_rep$SrvAgeComps_nLL), rce$comp_nll, tolerance = 1e-8)

  # rceattle's estimate is a stationary point of the survey selectivity and catchability here
  grad <- obj$gr(obj$par)
  expect_lt(max(abs(grad[names(obj$par) %in% c("ln_srvsel_devs", "srv_fixed_sel_pars", "ln_srv_q")])), 1e-3)

  # the fishery cpue from mean numbers over the year is below rceattle's start of year prediction, more so as F increases
  cpue_ratio <- sporc_rep$PredFishIdx[1,1,2:13,1,1] / rce$cpue_hat[2:13]
  expect_true(all(cpue_ratio < 1))
  expect_lt(cpue_ratio[12], cpue_ratio[1])

})
