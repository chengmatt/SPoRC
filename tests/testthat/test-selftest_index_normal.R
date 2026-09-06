library(SPoRC)
library(testthat)

# Simulate under an arithmetic-scale normal index error and refit: the operating
# model draws through the same error structure the estimation model fits, so
# median relative error should sit at zero. Routines shared through
# helper-selftest_features.R.

test_that("a normal index likelihood recovers SSB and catchability", {

  om <- selftest_make_om(SrvIdx_LikeType = "normal")
  sim_data <- simulation_data_to_SPoRC(sim_env = om, y = selftest_cfg$n_yrs, sim = 1)

  out <- selftest_run(selftest_build_input(sim_data, SrvIdx_LikeType = "normal"),
                      what = c("SSB", "srv_q"), seed = 101)

  expect_lt(abs(out$summ$SSB[["median_RE"]]), 0.025)
  expect_lt(abs(out$summ$srv_q[["median_RE"]]), 0.025)

})
