# Regression test. Expected SSB and recruitment come from a previously validated SPoRC fit, not from
# hand-derived values, so a mismatch means something moved a fitted result. See tests/README.md.

library(SPoRC)
library(testthat)
data("sgl_rg_dusky_data")

test_that("Dusky RTMB model produces expected results", {

  # The specification lives in helper-bridge_goa_dusky.R, which the case study
  # figure script sources as well, so a change to the specification moves this
  # test and the figures together.
  input_list <- build_goa_dusky_input(sgl_rg_dusky_data)

  data <- input_list$data
  parameters <- input_list$par
  mapping <- input_list$map

  # Fit model
  dusky_rtmb_model <- fit_model(data,
                                parameters,
                                mapping,
                                random = NULL,
                                newton_loops = 3,
                                silent = TRUE
  )

  dusky_rtmb_model$sdrep <- RTMB::sdreport(dusky_rtmb_model) # get standard error report

  ssb_expected_vec <- c(
    12797.54, 12250.75, 11909.20, 11599.08, 11192.37,
    10790.43, 10507.33, 10338.19, 10565.86, 11306.05,
    12302.58, 13455.75, 14217.22, 14785.66, 15393.80,
    15934.92, 15768.58, 16105.93, 17027.92, 18489.64,
    20460.80, 22348.89, 23729.31, 24165.39, 24860.83,
    25933.28, 26978.12, 28242.44, 29758.03, 31503.33,
    33118.19, 34220.20, 35018.81, 35787.45, 36175.29,
    36467.88, 35801.61, 35360.99, 35004.50, 34987.86,
    35038.66, 35694.16, 36410.03, 37326.83, 38208.01,
    38414.28, 38341.49, 37409.68
  )

  rec_expected_vec <- c(
    1.719130, 1.821452, 2.325670, 5.221616,
    6.055137, 6.385826, 3.798459, 4.855644,
    3.303878, 3.093645, 2.270344, 8.978767,
    5.840166, 17.847508, 12.477896, 10.773610,
    3.022077, 7.685572, 5.743625, 16.754576,
    3.156288, 9.194834, 18.978940, 2.694959,
    12.000002, 14.767476, 6.897940, 9.725435,
    9.246636, 4.339117, 4.744178, 5.182619,
    6.320576, 7.264755, 11.971134, 9.814706,
    7.959109, 17.536878, 4.838579, 7.376432,
    6.538344, 6.301670, 2.189419, 2.404738,
    2.061323, 2.393618, 2.529818, 2.838106
  )

  # Testing to see if we obtian same values back
  expect_equal(dusky_rtmb_model$rep$SSB[1,1,], ssb_expected_vec, tolerance = 1e-5)
  expect_equal(dusky_rtmb_model$rep$Rec[1,1,], rec_expected_vec, tolerance = 1e-5)
  expect_true(dusky_rtmb_model$sdrep$pdHess)
  expect_jnLL_decomposes(dusky_rtmb_model)

})
