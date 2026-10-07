# A recovery-conditioned tag likelihood conditions on the recaptures it fits, so a fleet that reports no
# tags adds nothing to its total. The Dirichlet-multinomial reads that total directly.

library(testthat)

# one cohort, one year at liberty, two regions and two fleets, the second reporting no tags
two_fleet_tag_nLL <- function(fleet2_obs, like_type = 5) {

  arr_dim <- c(1, 1, 1, 1, 2, 1, 1, 2) # liberty, season, cohort, pop, region, age, sex, fleet
  obs <- array(0, dim = arr_dim)
  pred <- array(0, dim = arr_dim)
  obs[1,1,1,1,,1,1,1] <- c(3, 5)
  obs[1,1,1,1,,1,1,2] <- fleet2_obs
  pred[1,1,1,1,,1,1,1] <- c(2.5, 4.5)
  pred[1,1,1,1,,1,1,2] <- c(6, 6)

  out <- get_conv_tag_likelihoods(n_conv_tag_cohorts = 1,
                                  conv_tag_release_indicator = matrix(c(1, 1, 1), nrow = 1),
                                  conv_tag_max_liberty = 1,
                                  n_yrs = 1,
                                  n_seas = 1,
                                  conv_tag_mixing_period = 1,
                                  n_fish_fleets = 2,
                                  use_conv_fish_tagging = c(1, 0),
                                  n_conv_tag_pop_pool = 1,
                                  n_regions = 2,
                                  n_conv_tag_age_pool = 1,
                                  n_conv_tag_sex_pool = 1,
                                  conv_tag_pop_pool = list(1),
                                  conv_tag_age_pool = list(1),
                                  conv_tag_sex_pool = list(1),
                                  conv_fish_tag_like = like_type,
                                  conv_fish_tag_nLL = array(0, dim = c(1, 1, 1, 2, 2)),
                                  obs_conv_tag_fish_recap = obs,
                                  pred_conv_tag_fish_recap = pred,
                                  addtotag = 0.001,
                                  ln_conv_fish_tag_theta = log(8),
                                  conv_tagged_fish = array(100, dim = c(1, 1, 1, 1)))
  sum(out)
}

test_that("a fleet reporting no tags adds nothing to the recovery-conditioned total", {

  for(like_type in c(3, 5)) {
    expect_equal(two_fleet_tag_nLL(c(0, 0), like_type), two_fleet_tag_nLL(c(40, 25), like_type), info = like_type)
  } # end like_type loop

  # the Dirichlet-multinomial over the reporting fleet's two cells
  obs_fit <- c(3, 5) + 0.001
  pred_fit <- c(2.5, 4.5) + 0.001
  expected <- -ddirmult(obs = obs_fit / sum(obs_fit), pred = pred_fit / sum(pred_fit), Ntotal = sum(obs_fit), ln_theta = log(8), TRUE)
  expect_equal(two_fleet_tag_nLL(c(40, 25), 5), expected, tolerance = 1e-10)
})
