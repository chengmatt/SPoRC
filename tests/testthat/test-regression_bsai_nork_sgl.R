# Pinned regression test. The expected SSB and recruitment vectors are output from a
# previously validated SPoRC fit of this assessment, not hand-derived values. A mismatch
# means a change moved a fitted result, which is a bug unless the numerical change was
# intended. If it was intended, re-baseline deliberately and say why in NEWS.md. Do not
# paste in fresh output to make the test pass. See tests/README.md.
#
# The companion test-regression_bsai_nork_bridge.R checks the same configuration against
# the ADMB assessment's own output without optimizing, so a failure here that does not
# also fail there is an optimizer or setup change rather than a model change.
#
# This refit estimates the uncapped logistic survey selectivity, where the bridge
# supplies the assessment's age 30 edge hold as a fixed input. That is the only
# specification difference between the two.

library(SPoRC)
library(testthat)
data("sgl_rg_bsai_nork_data")

test_that("Single-region BSAI northern rockfish RTMB model produces expected results", {

  dat <- sgl_rg_bsai_nork_data
  n_yrs <- length(dat$years)

  input_list <- seed_bsai_nork_mle(build_bsai_nork_input(dat), dat)

  bsai_nork_rtmb_model <- fit_model(input_list$data,
                                    input_list$par,
                                    input_list$map,
                                    random = NULL,
                                    newton_loops = 3,
                                    silent = TRUE
  )

  bsai_nork_rtmb_model$sdrep <- RTMB::sdreport(bsai_nork_rtmb_model)

  ssb_expected_vec <- c(
    57740.4999, 60696.0217, 64269.6049, 68368.9543, 71943.2275,
    75107.7206, 77967.2874, 80746.3118, 83477.3643, 86199.0976,
    88846.4073, 91380.3512, 93744.3794, 96163.0604, 98874.5680,
    102606.2898, 106634.5214, 110154.2458, 113484.5160, 116347.9750,
    118536.7101, 121661.2015, 123702.3996, 124972.7874, 126117.1736,
    127664.1909, 131106.2219, 135204.1799, 139368.5212, 143349.4946,
    144299.5936, 144766.5836, 144870.7333, 144739.1928, 144568.6980,
    145362.5493, 147302.4516, 149622.2562, 150146.5025, 148851.7429,
    147695.6014, 146974.3804, 145001.0458, 141765.1221, 138853.7910,
    136540.3010, 131852.1707
  )

  rec_expected_vec <- c(
    47.561863, 41.901182, 42.413738, 49.350764, 58.946852,
    50.144371, 41.606933, 37.270746, 30.388878, 34.004510,
    103.046052, 196.385359, 35.992298, 29.510389, 39.433211,
    167.734090, 30.904662, 26.174863, 23.225069, 107.219911,
    47.348523, 128.799357, 127.469704, 74.322824, 69.011592,
    25.198595, 30.639545, 25.156454, 58.007707, 49.517558,
    40.141810, 149.959080, 37.598547, 41.061950, 43.420095,
    49.473703, 22.580788, 38.658575, 21.690772, 28.570015,
    36.021755, 36.208181, 39.553889, 23.170632, 44.669352,
    44.669352, 44.669352
  )

  expect_equal(as.vector(bsai_nork_rtmb_model$rep$SSB)[1:n_yrs], ssb_expected_vec,
               tolerance = 1e-2)
  expect_equal(as.vector(bsai_nork_rtmb_model$rep$Rec)[1:n_yrs], rec_expected_vec,
               tolerance = 1e-2)
  expect_true(bsai_nork_rtmb_model$sdrep$pdHess)
  expect_jnLL_decomposes(bsai_nork_rtmb_model)

  # The refit stays on the assessment. This is not a pinned number: it is the
  # assessment's own spawning biomass, so it holds the refit to the bridge.
  ssb_fit <- as.vector(bsai_nork_rtmb_model$rep$SSB)[1:n_yrs]
  expect_lt(max(abs(ssb_fit / dat$admb$SSB - 1)), 2e-2)

  # Recruitment likewise, but only over the years the deviations are estimated. The
  # assessment's three terminal recruits are dev free and it builds them as the MEAN
  # recruitment, exp(mean_log_rec + sigmaR^2 / 2), while SPoRC builds them as the
  # median, exp(R0). With sigmaR fixed at 0.75 that is a factor of exp(sigmaR^2 / 2)
  # = 1.32 between the two by construction, which the refit's own shift in the mean
  # level carries to 1.36. This is a convention difference, not a wiring bug: the
  # bridge test reproduces those same years exactly because it seeds R0 with the bias
  # correction already in it. Asserted here so that the gap stays understood rather
  # than drifting.
  rec_fit <- as.vector(bsai_nork_rtmb_model$rep$Rec)[1:n_yrs]
  n_fixed <- 3
  est_yrs <- seq_len(n_yrs - n_fixed)
  expect_lt(max(abs(rec_fit[est_yrs] / dat$admb$Rec[est_yrs] - 1)), 2e-2)

  term_yrs <- seq(n_yrs - n_fixed + 1, n_yrs)
  term_ratio <- dat$admb$Rec[term_yrs] / rec_fit[term_yrs]
  expect_equal(term_ratio, rep(term_ratio[1], n_fixed), tolerance = 1e-8,
               ignore_attr = TRUE)
  expect_gt(term_ratio[1], exp(dat$sigmaR^2 / 2) * 0.98)
  expect_lt(term_ratio[1], exp(dat$sigmaR^2 / 2) * 1.05)
})
