# Pinned regression test. The expected SSB and recruitment vectors are output from a
# previously validated SPoRC fit of this assessment, not hand-derived values. A mismatch
# means a change moved a fitted result, which is a bug unless the numerical change was
# intended. If it was intended, re-baseline deliberately and say why in NEWS.md. Do not
# paste in fresh output to make the test pass. See tests/README.md.
#
# The companion test-regression_bsai_pop_bridge.R checks the same configuration against
# the ADMB assessment's own output without optimizing, so a failure here that does not
# also fail there is an optimizer or setup change rather than a model change. The two
# share one builder, so there is no specification difference between them.

library(SPoRC)
library(testthat)
data("sgl_rg_bsai_pop_data")

test_that("Single-region BSAI Pacific ocean perch RTMB model produces expected results", {

  dat <- sgl_rg_bsai_pop_data
  n_yrs <- length(dat$years)
  yr_ind <- match(dat$years, dat$admb$yrs)

  input_list <- seed_bsai_pop_mle(build_bsai_pop_input(dat), dat)

  bsai_pop_rtmb_model <- fit_model(input_list$data,
                                   input_list$par,
                                   input_list$map,
                                   random = NULL,
                                   newton_loops = 3,
                                   silent = TRUE
  )

  bsai_pop_rtmb_model$sdrep <- RTMB::sdreport(bsai_pop_rtmb_model)

  ssb_expected_vec <- c(
    379741.4498, 378411.6115, 366122.0445, 360354.5342, 337634.9144,
    291031.1825, 244841.2630, 209111.4235, 183449.1484, 162121.5182,
    146908.3902, 129881.6343, 123055.2531, 114140.5136, 108194.2815,
    95653.0633, 85682.2163, 76674.8701, 74072.0375, 72837.9318,
    72358.8422, 72872.3284, 73990.1369, 76585.9634, 79867.7976,
    83980.4334, 89295.0812, 95702.0298, 103329.1484, 112613.5250,
    122223.8500, 130953.6221, 144371.0773, 157640.7533, 171417.4809,
    186835.9350, 201748.3435, 214891.5203, 227792.3202, 239221.9527,
    247216.9563, 254420.0814, 261156.3667, 268215.3219, 276780.8551,
    288897.6717, 303012.3660, 316622.8044, 329008.7260, 341905.7266,
    354110.9466, 362801.0271, 367599.6780, 370964.7683, 371981.1440,
    372989.8992, 374286.1331, 375189.5801, 374842.7446, 371659.3058,
    364978.2189, 359443.9768, 355441.1404, 351181.7864, 345907.2989
  )

  rec_expected_vec <- c(
    207.166524, 55.124140, 50.365531, 62.646135, 195.519347,
    258.768670, 44.551354, 30.301247, 28.947194, 33.521711,
    30.425322, 22.419861, 28.292917, 44.499482, 35.016903,
    41.448102, 20.462845, 23.426039, 21.337636, 52.902785,
    46.153291, 28.376367, 54.256843, 97.909887, 115.566339,
    69.211869, 53.753334, 289.589066, 55.124829, 137.773875,
    59.437653, 209.677456, 109.809514, 66.042301, 31.953609,
    49.103885, 41.067666, 132.573048, 100.758868, 223.646003,
    49.087592, 194.691429, 62.315521, 298.795738, 40.901887,
    105.328096, 43.388179, 149.069664, 180.437661, 92.419886,
    80.713681, 220.323866, 71.406103, 76.817899, 129.334654,
    116.377270, 84.520032, 155.134548, 71.335479, 151.408426,
    32.287011, 32.349114, 87.204672, 87.204672, 87.204672
  )

  expect_equal(as.vector(bsai_pop_rtmb_model$rep$SSB)[1:n_yrs], ssb_expected_vec,
               tolerance = 1e-2)
  expect_equal(as.vector(bsai_pop_rtmb_model$rep$Rec)[1:n_yrs], rec_expected_vec,
               tolerance = 1e-2)
  expect_true(bsai_pop_rtmb_model$sdrep$pdHess)
  expect_jnLL_decomposes(bsai_pop_rtmb_model)

  # The refit stays on the assessment. These are not pinned numbers: they are the
  # assessment's own estimates, so they hold the refit to the bridge. The initial
  # equilibrium recruitment, the mean F and both catchabilities are recovered to
  # better than 0.05 percent, and spawning biomass to better than 0.2 percent over
  # all 65 years, which is what makes the recruitment level shift below the only
  # place the two models part company.
  fit_par <- bsai_pop_rtmb_model$env$parList(bsai_pop_rtmb_model$env$last.par.best)
  s2 <- dat$sigmaR^2 / 2
  expect_equal(exp(as.vector(fit_par$ln_M)[1]), dat$mle$M, tolerance = 1e-3)
  expect_equal(as.vector(fit_par$ln_rinit)[1] - s2, dat$mle$log_rinit, tolerance = 1e-3)
  expect_equal(as.vector(fit_par$ln_F_mean)[1], dat$mle$log_avg_fmort, tolerance = 1e-3)
  expect_equal(exp(as.vector(fit_par$ln_srv_q)[1:2]), dat$mle$q_srv, tolerance = 1e-3)

  ssb_fit <- as.vector(bsai_pop_rtmb_model$rep$SSB)[1:n_yrs]
  expect_lt(max(abs(ssb_fit / dat$admb$SSB[yr_ind] - 1)), 3e-3)

  # Recruitment likewise, but only over the years the deviations are estimated. The
  # assessment's three terminal recruits are dev free and it builds them as the mean
  # recruitment, exp(mean_log_rec + sigmaR^2 / 2). SPoRC builds them as exp(R0), and
  # R0 is seeded with that same bias correction, so unlike the BSAI northern rockfish
  # case the two do not differ by exp(sigmaR^2 / 2) here: this configuration runs the
  # bias ramp at one, which centres the recruitment penalty on the shift the seeds
  # carry and leaves R0 where it was put.
  #
  # What is left is a level shift. The assessment declares its recruitment deviations
  # as a dev_vector, which is constrained to sum to zero, so it cannot slide the
  # overall level between mean_log_rec and the deviations. SPoRC's deviations are
  # free and it does slide, by about 1.6 percent in log space, which lands as a 6.8
  # percent gap in the three terminal recruits and 0.17 percent in terminal spawning
  # biomass. The gap is asserted rather than tolerated so that it stays understood: it
  # is exactly the shift in R0, and it is the same number in all three terminal years.
  rec_fit <- as.vector(bsai_pop_rtmb_model$rep$Rec)[1:n_yrs]
  n_fixed <- dat$fixedrec
  est_yrs <- seq_len(n_yrs - n_fixed)
  expect_lt(max(abs(rec_fit[est_yrs] / dat$admb$Rec[yr_ind][est_yrs] - 1)), 2e-2)

  term_yrs <- seq(n_yrs - n_fixed + 1, n_yrs)
  term_ratio <- dat$admb$Rec[yr_ind][term_yrs] / rec_fit[term_yrs]
  expect_equal(term_ratio, rep(term_ratio[1], n_fixed), tolerance = 1e-8,
               ignore_attr = TRUE)

  # The gap is the shift in the recruitment level and nothing else.
  r0_shift <- (dat$mle$mean_log_rec + s2) - as.vector(fit_par$ln_global_R0)[1]
  expect_equal(term_ratio[1], exp(r0_shift), tolerance = 1e-6, ignore_attr = TRUE)
  expect_lt(abs(r0_shift), 0.08)
})
