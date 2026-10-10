# simulation_data_to_SPoRC with one tagging flag per fishery fleet, as Setup_Mod_Tagging stores them.
# a closed loop conditioned on such a fit reads the whole vector back from the operating model.

library(SPoRC)
library(testthat)

test_that("simulation_data_to_SPoRC reads one tagging flag per fishery fleet", {
  sim_obj <- objective_setup_sim()
  n_yrs <- sim_obj$n_years

  # two fishery fleets, neither fit to tag recaptures
  sim_obj$use_conv_fish_tagging <- c(0, 0)
  sim_data <- simulation_data_to_SPoRC(sim_env = sim_obj, y = n_yrs, sim = 1)

  expect_equal(sim_data$use_conv_fish_tagging, c(0, 0))
  expect_null(sim_data$conv_tag_release_indicator) # no fleet tagged, so no tag data
  expect_null(sim_data$obs_conv_tag_fish_recap)
  expect_null(sim_data$conv_tagged_fish)
  expect_null(sim_data$conv_tagged_fish_attr)
  expect_null(sim_data$n_tag_cohorts)

  # second fleet tagged, so releases up to the last assessment year are kept
  sim_obj$use_conv_fish_tagging <- c(0, 1)
  sim_data <- simulation_data_to_SPoRC(sim_env = sim_obj, y = n_yrs - 5, sim = 1)

  expect_equal(sim_data$n_tag_cohorts, sum(sim_obj$conv_tag_release_indicator[, 2] <= n_yrs - 5))
  expect_true(all(sim_data$conv_tag_release_indicator[, 2] <= n_yrs - 5))
})
