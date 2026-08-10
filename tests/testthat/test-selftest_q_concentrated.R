library(SPoRC)
library(testthat)

# Concentrated catchability: the analytic arithmetic and geometric solves fix
# ln_srv_q and derive q from the data inside the likelihood, so the self-test
# checks the solved q and the biomass scale it implies are both recovered.
# Machinery shared through helper-selftest_features.R.

test_that("an arithmetically concentrated catchability recovers SSB and q", {

  om <- selftest_make_om()
  sim_data <- simulation_data_to_SPoRC(sim_env = om, y = selftest_cfg$n_yrs, sim = 1)

  out <- selftest_run(selftest_build_input(sim_data, srv_q_type = "arith"),
                      what = c("SSB", "srv_q"), seed = 103)

  expect_lt(abs(out$summ$SSB[["median_RE"]]), 0.025)
  expect_lt(abs(out$summ$srv_q[["median_RE"]]), 0.025)

})

test_that("a geometrically concentrated catchability recovers SSB and q", {

  om <- selftest_make_om()
  sim_data <- simulation_data_to_SPoRC(sim_env = om, y = selftest_cfg$n_yrs, sim = 1)

  out <- selftest_run(selftest_build_input(sim_data, srv_q_type = "geo"),
                      what = c("SSB", "srv_q"), seed = 104)

  expect_lt(abs(out$summ$SSB[["median_RE"]]), 0.025)
  expect_lt(abs(out$summ$srv_q[["median_RE"]]), 0.025)

})
