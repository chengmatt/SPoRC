# Ricker recruitment regenerated from the fitted curve each replicate, so the stock-recruit
# relationship itself is what has to come back. Run in the deterministic limit, where the data
# essentially trace the curve and any bias left is structural, so the tolerances are near-exact.
#
# Under realistic process error the same design gives steepness medians a few percent high with R0
# correspondingly low, which is the shape of the stock-recruit ridge rather than a defect.

library(SPoRC)
library(testthat)

test_that("Ricker recruitment recovers SSB, recruitment, R0, and steepness", {

  sigR <- 0.05
  om <- selftest_make_om(
    recruitment_opt = "ricker_rec",
    sigmaR = sigR,
    idx_se_om = 0.02,
    iss_om = 2000
  )
  sim_data <- simulation_data_to_SPoRC(sim_env = om, y = selftest_cfg$n_yrs, sim = 1)

  out <- selftest_run(
    selftest_build_input(sim_data, rec_model = "ricker_rec", sigmaR = sigR),
    what = c("SSB", "Rec", "R0", "h_trans"),
    sim_recruitment = "model",
    seed = 105
  )

  expect_lt(abs(out$summ$SSB[["median_RE"]]), 0.01)
  expect_lt(abs(out$summ$Rec[["median_RE"]]), 0.01)
  expect_lt(abs(out$summ$R0[["median_RE"]]), 0.01)
  expect_lt(abs(out$summ$h_trans[["median_RE"]]), 0.01)

})
