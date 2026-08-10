# Pinned regression test. The expected SSB and recruitment vectors are output from a
# previously validated SPoRC fit of this assessment, not hand-derived values. A mismatch
# means a change moved a fitted result, which is a bug unless the numerical change was
# intended. If it was intended, re-baseline deliberately and say why in NEWS.md. Do not
# paste in fresh output to make the test pass. See tests/README.md.
#
# The companion test-regression_rebs_bridge.R checks the same configuration against the
# ADMB assessment's own output without optimizing, so a failure here that does not also
# fail there is an optimizer or setup change rather than a model change.

library(SPoRC)
library(testthat)
data("sgl_rg_rebs_data")

test_that("Single-region BSAI rougheye RTMB model produces expected results", {

  dat <- sgl_rg_rebs_data
  n_yrs <- length(dat$years)

  input_list <- seed_rebs_mle(build_rebs_input(dat), dat)

  rebs_rtmb_model <- fit_model(input_list$data,
                               input_list$par,
                               input_list$map,
                               random = NULL,
                               newton_loops = 3,
                               silent = TRUE
  )

  rebs_rtmb_model$sdrep <- RTMB::sdreport(rebs_rtmb_model)

  ssb_expected_vec <- c(
    4965.5592, 4936.9389, 4304.1979, 3645.9628, 3573.7630,
    3541.1350, 3617.2411, 3732.6991, 3861.9195, 4003.2285,
    4146.5810, 4279.1812, 4398.6338, 4371.2817, 4095.3453,
    4130.3425, 3894.1081, 3718.7256, 3596.2786, 3535.2058,
    3340.8467, 3121.9685, 3019.7667, 2956.2071, 2900.5613,
    2778.4887, 2737.3102, 2721.0101, 2705.3356, 2712.1847,
    2684.7409, 2667.7099, 2646.3977, 2622.5667, 2602.0940,
    2604.5752, 2596.1153, 2577.7748, 2604.0343, 2653.3378,
    2727.1488, 2806.5566, 2896.8879, 2985.1564, 3075.9058,
    3208.4183, 3377.8036, 3554.3408
  )

  rec_expected_vec <- c(
    1.262295, 1.202792, 1.156209, 1.118094, 1.100478,
    1.112351, 1.135931, 1.137713, 1.069551, 0.941126,
    0.805440, 0.696163, 0.612661, 0.547931, 0.503272,
    0.481688, 0.485209, 0.509930, 0.558100, 0.625700,
    0.720575, 0.851462, 1.018484, 1.347845, 2.177848,
    2.425607, 2.268003, 2.486581, 3.270275, 2.198060,
    1.839679, 1.947456, 2.295805, 2.432463, 2.583417,
    2.258798, 1.925996, 27.860708, 1.731920, 2.094088,
    2.655659, 2.978946, 2.266719, 1.626486, 1.400226,
    1.380554, 1.380554, 1.380554
  )

  expect_equal(as.vector(rebs_rtmb_model$rep$SSB)[1:n_yrs], ssb_expected_vec, tolerance = 1e-2)
  expect_equal(as.vector(rebs_rtmb_model$rep$Rec)[1:n_yrs], rec_expected_vec, tolerance = 1e-2)
  expect_true(rebs_rtmb_model$sdrep$pdHess)
  expect_jnLL_decomposes(rebs_rtmb_model)

  # The refit stays on the assessment. This is not a pinned number: it is the
  # assessment's own spawning biomass, so it holds the refit to the bridge.
  expect_lt(max(abs(as.vector(rebs_rtmb_model$rep$SSB)[1:n_yrs] / dat$admb$SSB - 1)), 1e-2)
})
