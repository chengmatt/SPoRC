# refits SPoRC from the September 2026 EBS pollock Model 26.0 (Rceattle) estimate and checks it against Rceattle
# the fishery cpue is predicted from mean numbers over the year here, so SSB differs by up to 7.2% before 1980
# (1966) and by 1.4% after
library(SPoRC)
library(testthat)
data("sgl_rg_ebswp_m26_data")

test_that("EBS pollock refit lands near Rceattle Model 26.0", {

  inputs <- sgl_rg_ebswp_m26_data$inputs
  rce <- sgl_rg_ebswp_m26_data$rceattle
  years <- inputs$years

  # setup model and start at the rceattle estimate
  input_list <- suppressWarnings(build_ebs_pollock_m26_input(inputs = inputs, rce = rce))
  input_list <- seed_ebs_pollock_m26_rceattle(input_list = input_list, rce = rce)

  ebs_pollock_model <- fit_model(data = input_list$data,
                                 parameters = input_list$par,
                                 mapping = input_list$map,
                                 random = NULL,
                                 newton_loops = 2,
                                 silent = TRUE)

  expect_lt(max(abs(ebs_pollock_model$gr(ebs_pollock_model$optim$par))), 1e-6)
  expect_equal(ebs_pollock_model$optim$objective, 1481.62735844, tolerance = 1e-7)

  # compare SSB and recruitment
  sporc_rep <- ebs_pollock_model$report(ebs_pollock_model$env$last.par.best)
  ssb_rel <- as.vector(sporc_rep$SSB) / rce$ssb - 1
  rec_rel <- as.vector(sporc_rep$Rec) / rce$rec - 1

  expect_lt(max(abs(ssb_rel)), 0.075)
  expect_lt(max(abs(ssb_rel[years >= 1980])), 0.015)
  expect_lt(max(abs(ssb_rel[years >= 2000])), 0.005)
  expect_lt(max(abs(rec_rel)), 0.045)
  expect_lt(max(abs(rec_rel[years >= 1980])), 0.015)

  expect_jnLL_decomposes(ebs_pollock_model)

})
