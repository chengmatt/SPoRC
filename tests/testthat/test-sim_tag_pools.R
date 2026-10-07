# Tag recapture draws against the likelihood the estimation model fits them with: a negative binomial over
# a pooled group of ages is one draw, and a cohort with nothing to recapture writes zeros rather than failing.

library(SPoRC)
library(testthat)

n_ages_tag <- 4

# one cohort, liberty year and season; one population, region, sex and fleet; four ages
tag_draw <- function(like, pred_per_age, pools = NULL, attr = "a", theta = 1) {
  dims <- c(1, 1, 1, 1, 1, n_ages_tag, 1, 1, 1) # liberty, season, cohort, pop, region, age, sex, fleet, sim
  pred <- array(pred_per_age, dim = dims)
  simulate_conv_tag_fish_recaptures(conv_fish_tag_like = like,
                                    tag_recaptures_attr = attr,
                                    conv_tagged_fish = array(1000, dim = c(1, 1, n_ages_tag, 1, 1)),
                                    pred_conv_tag_fish_recap = pred,
                                    obs_conv_tag_fish_recap = array(0, dim = dims),
                                    ln_conv_fish_tag_theta = log(theta),
                                    ry = 1, rseas = 1, tc = 1, sim = 1,
                                    n_pop = 1, n_regions = 1, n_ages = n_ages_tag, n_sexes = 1, n_fish_fleets = 1,
                                    tag_pools = pools)
}

test_that("a negative binomial over pooled ages is one draw at the pooled mean", {

  pools <- list(pop = list(1), age = list(1:n_ages_tag), sex = list(1))
  set.seed(14)
  pooled <- replicate(20000, tag_draw(1, 5, pools))
  set.seed(14)
  by_age <- replicate(20000, sum(tag_draw(1, 5)))

  # one draw, written into the group's first age
  expect_true(all(pooled[1,1,1,1,1,-1,1,1,1,] == 0))
  group_count <- pooled[1,1,1,1,1,1,1,1,1,]

  # mean 20 either way, but variance mu + mu^2 / theta = 420 for one draw and 120 for four summed
  expect_equal(mean(group_count), 20, tolerance = 0.03)
  expect_equal(stats::var(group_count), 420, tolerance = 0.06)
  expect_equal(stats::var(by_age), 120, tolerance = 0.06)
})

test_that("pools of one level draw exactly what they drew before", {

  set.seed(3)
  with_pools <- tag_draw(1, 5, list(pop = list(1), age = as.list(1:n_ages_tag), sex = list(1)))
  set.seed(3)
  without <- tag_draw(1, 5)
  expect_identical(with_pools, without)
})

test_that("a cohort with nothing to recapture writes zeros rather than failing", {

  for(like in c(3, 5)) expect_true(all(tag_draw(like, 0.01) == 0)) # under half a recapture expected in total
})

test_that("a recovery-conditioned draw spreads its total over the fleets that report tags", {

  dims <- c(1, 1, 1, 1, 1, 1, 1, 2, 1) # liberty, season, cohort, pop, region, age, sex, fleet, sim
  pred <- array(c(30, 70), dim = dims) # fleet one expects 30, fleet two 70
  for(like in c(3, 5)) {
    set.seed(6)
    drawn <- simulate_conv_tag_fish_recaptures(conv_fish_tag_like = like,
                                               tag_recaptures_attr = "none",
                                               conv_tagged_fish = array(1000, dim = c(1, 1, 1, 1, 1)),
                                               pred_conv_tag_fish_recap = pred,
                                               obs_conv_tag_fish_recap = array(0, dim = dims),
                                               ln_conv_fish_tag_theta = log(5),
                                               ry = 1, rseas = 1, tc = 1, sim = 1,
                                               n_pop = 1, n_regions = 1, n_ages = 1, n_sexes = 1, n_fish_fleets = 2,
                                               tag_fleets = 1)
    expect_equal(drawn[1,1,1,1,1,1,1,1,1], 30, info = like) # the reporting fleet takes the whole of its own total
    expect_equal(drawn[1,1,1,1,1,1,1,2,1], 0, info = like)
  } # end like loop
})
