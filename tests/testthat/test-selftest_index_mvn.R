library(SPoRC)
library(testthat)

# Simulate under a multivariate normal index error (survey-wide scaling error
# with a strong common factor) and refit with the same covariance. The shared
# factor shifts a replicate's whole series, which catchability absorbs, so the
# catchability median converges slower than under independent error and has
# a looser tolerance. Routines shared through helper-selftest_features.R.

test_that("an mvn index likelihood recovers SSB and catchability", {

  om_ln <- selftest_make_om()
  S <- selftest_mvn_cov(mean(om_ln$TrueSrvIdx[1, , 1, 1, 1]))
  use_all <- array(0, dim = c(1, selftest_cfg$n_yrs, 1, 1)); use_all[1, , 1, 1] <- 1

  om <- selftest_make_om(SrvIdx_LikeType = "mvn", SrvIdx_Cov = list(S), UseSrvIdx = use_all)
  sim_data <- simulation_data_to_SPoRC(sim_env = om, y = selftest_cfg$n_yrs, sim = 1)

  out <- selftest_run(selftest_build_input(sim_data, SrvIdx_LikeType = "mvn", SrvIdx_Cov = list(S)),
                      what = c("SSB", "srv_q"), seed = 102)

  expect_lt(abs(out$summ$SSB[["median_RE"]]), 0.025)
  expect_lt(abs(out$summ$srv_q[["median_RE"]]), 0.05)

})
