# refits SPoRC from the saved 2025 NSEI sablefish TMB estimate and checks it against the converged TMB fit
# the fishery cpue and mark-recapture are fishery indices here (mean numbers over the year) and tmb predicts
# them in august, so the fit sits about 2% from tmb (recruitment in 1981, informed only by the fishery cpue, 51%)
library(SPoRC)
library(testthat)
data("sgl_rg_nsei_sablefish_data")

test_that("NSEI sablefish refit lands near the converged 2025 TMB assessment", {

  tmb_data <- sgl_rg_nsei_sablefish_data$inputs
  tmb <- sgl_rg_nsei_sablefish_data$tmb_converged
  tmb_par <- sgl_rg_nsei_sablefish_data$tmb_saved$par
  n_yrs <- tmb_data$nyr

  # setup model and start at the saved tmb estimate
  input_list <- suppressWarnings(build_nsei_sablefish_input(tmb_data = tmb_data, tmb_par = tmb_par))
  input_list <- seed_nsei_sablefish_tmb(input_list = input_list, tmb_par = tmb_par)

  # fit model, no newton steps since the 2017 - 2025 survey selectivity is knife edge
  nsei_model <- fit_model(data = input_list$data,
                          parameters = input_list$par,
                          mapping = input_list$map,
                          random = NULL,
                          newton_loops = 0,
                          silent = TRUE)

  # same 150 parameters as tmb, less its five SPR rates
  expect_equal(length(nsei_model$par), 150)
  expect_lt(max(abs(nsei_model$gr(nsei_model$optim$par))), 1e-2)
  expect_equal(nsei_model$optim$objective, 1955.4876, tolerance = 1e-7)

  # compare SSB, F and recruitment
  sporc_rep <- nsei_model$report(nsei_model$env$last.par.best)
  ssb_rel <- as.vector(sporc_rep$SSB) / tmb$ssb - 1
  F_rel <- sporc_rep$Fmort[1,,1,1] / tmb$Fmort - 1
  rec_rel <- as.vector(sporc_rep$Rec) / tmb$pred_rec - 1

  expect_lt(max(abs(ssb_rel)), 0.025)
  expect_lt(abs(ssb_rel[n_yrs]), 0.025)
  expect_lt(max(abs(F_rel)), 0.025)
  expect_lt(max(abs(rec_rel[26:n_yrs])), 0.025) # 2000 on
  expect_lt(abs(rec_rel[n_yrs]), 0.02)
  expect_equal(max(sporc_rep$Fmort[1,,1,2]), 0) # mark-recapture fleet doesn't fish

  expect_jnLL_decomposes(nsei_model)
})
