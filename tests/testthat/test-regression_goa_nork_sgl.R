# Regression test. Expected SSB and recruitment come from a previously validated SPoRC fit, not from
# hand-derived values, so a mismatch means something moved a fitted result. See tests/README.md.
#
# test-regression_goa_nork_bridge.R checks the same configuration without optimizing, so a failure here
# that does not also fail there is an optimizer or setup change, not a model change.

library(SPoRC)
library(testthat)
data("sgl_rg_goa_nork_data")

test_that("Single-region GOA northern rockfish RTMB model produces expected results", {

  dat <- sgl_rg_goa_nork_data
  n_yrs <- length(dat$years)

  input_list <- seed_goa_nork_mle(build_goa_nork_input(dat), dat)

  goa_nork_rtmb_model <- fit_model(input_list$data,
                                   input_list$par,
                                   input_list$map,
                                   random = NULL,
                                   newton_loops = 3,
                                   silent = TRUE
  )

  goa_nork_rtmb_model$sdrep <- RTMB::sdreport(goa_nork_rtmb_model)

  ssb_expected_vec <- c(
    52739.7815, 51381.9664, 48965.0920, 44970.3088, 38448.7819,
    29286.0799, 24559.2916, 22200.3607, 20639.6941, 19985.3438,
    20143.1075, 19720.1564, 19366.0254, 19565.6134, 20053.1045,
    20847.2107, 22197.3108, 24837.5684, 28243.6679, 32258.6517,
    36553.1279, 40550.0565, 43358.3540, 46281.1425, 50630.9437,
    55811.6301, 61241.2189, 66386.8070, 70787.3271, 74603.1505,
    78202.6437, 80632.6637, 81749.6912, 84116.5213, 85628.6311,
    86548.3499, 87635.6767, 88135.2216, 87979.7562, 86330.9976,
    85458.0187, 84730.3243, 84251.2578, 83475.5291, 83404.8227,
    83674.4176, 83687.0817, 83781.7173, 83463.5912, 82533.8701,
    80890.6962, 78756.6685, 75240.7775, 71393.0952, 67613.5375,
    63938.5823, 60553.6834, 58016.4290, 55295.3409, 52446.4223,
    49828.8565, 47282.9249, 45060.6398, 43227.9327
  )

  rec_expected_vec <- c(
    11.292193, 12.647912, 12.443073, 11.138207, 10.488693, 11.050747,
    11.131869, 11.464699, 15.419174, 21.304245, 16.187207, 62.540293,
    18.638468, 18.102690, 29.313303, 19.069898, 20.961201, 70.266253,
    40.965609, 22.718040, 23.258926, 27.313346, 34.268332, 41.091350,
    24.132269, 64.867869, 36.534747, 16.487876, 23.007692, 25.228578,
    10.310866, 24.250671, 16.284170, 14.370434, 9.860093, 51.418988,
    34.365523, 22.941036, 23.409333, 41.364471, 17.548764, 13.308349,
    16.827808, 8.672576, 3.806243, 3.789356, 6.268051, 5.064427,
    6.081825, 6.216343, 7.025763, 4.534942, 3.466788, 5.756939, 3.324092,
    6.283911, 3.941711, 5.236298, 6.357308, 7.687046, 9.225481,
    10.088509, 10.472627, 10.951033
  )

  expect_equal(as.vector(goa_nork_rtmb_model$rep$SSB)[1:n_yrs], ssb_expected_vec, tolerance = 1e-2)
  expect_equal(as.vector(goa_nork_rtmb_model$rep$Rec)[1:n_yrs], rec_expected_vec, tolerance = 1e-2)
  expect_true(goa_nork_rtmb_model$sdrep$pdHess)
  expect_jnLL_decomposes(goa_nork_rtmb_model)

  # The refit stays on the assessment. This is not a pinned number: it is the
  # assessment's own spawning biomass, so it holds the refit to the bridge.
  expect_lt(max(abs(as.vector(goa_nork_rtmb_model$rep$SSB)[1:n_yrs] / dat$admb$SSB - 1)), 1e-2)
})
