library(SPoRC)
library(testthat)

# set_data_indicator_unused() is how MSE closed-loop runs and retrospective
# peels withhold data from the estimation model. Silently failing to zero an
# indicator leaks future data into the fit, so every data type is checked.

n_pop <- 2
n_regions <- 3
n_years <- 5
n_seas <- 2
n_ages <- 4
n_sexes <- 2
n_fish_fleets <- 2
n_srv_fleets <- 1

# Names grouped by the array shape the function expects.
fish_agg <- c("UseCatch", "UseDiscard", "UseFishIdx",
              "UseFishAgeComps", "UseFishLenComps",
              "UseFishAgeComps_discard", "UseFishLenComps_discard")
fish_pop <- paste0(fish_agg, "_pop")
srv_agg <- c("UseSrvIdx", "UseSrvAgeComps", "UseSrvLenComps")
srv_pop <- paste0(srv_agg, "_pop")

make_data <- function(n_cohorts = 4, use_tagging = 1) {
  d <- list(years = 1:n_years)

  for(nm in fish_agg) d[[nm]] <- array(1, dim = c(n_regions, n_years, n_seas, n_fish_fleets))
  for(nm in fish_pop) d[[nm]] <- array(1, dim = c(n_pop, n_regions, n_years, n_seas, n_fish_fleets))
  for(nm in srv_agg) d[[nm]] <- array(1, dim = c(n_regions, n_years, n_seas, n_srv_fleets))
  for(nm in srv_pop) d[[nm]] <- array(1, dim = c(n_pop, n_regions, n_years, n_seas, n_srv_fleets))

  d$use_conv_fish_tagging <- use_tagging
  # Column 2 holds the release year; cohorts are released in years 1, 2, 3, 4.
  d$conv_tag_release_indicator <- cbind(rep(1, n_cohorts), seq_len(n_cohorts))
  d$n_conv_tag_cohorts <- n_cohorts
  d$conv_tagged_fish <- array(1, dim = c(n_cohorts, n_pop, n_ages, n_sexes))
  d$obs_conv_tag_fish_recap <- array(1, dim = c(2, n_seas, n_cohorts, n_pop,
                                                n_regions, n_ages, n_sexes, n_fish_fleets))
  return(d)
}

# Pull the year slice of an indicator array regardless of how many leading
# dimensions it carries.
year_slice <- function(arr, yrs, yr_dim) {
  idx <- rep(list(quote(expr = )), length(dim(arr)))
  idx[[yr_dim]] <- yrs
  do.call(`[`, c(list(arr), idx, drop = FALSE))
}

yr_dim_of <- function(nm) if(grepl("_pop$", nm)) 3 else 2

all_names <- c(fish_agg, fish_pop, srv_agg, srv_pop)

test_that("the default call zeroes every indicator in the named years", {
  d <- make_data()
  out <- SPoRC::set_data_indicator_unused(d, unused_years = c(4, 5))

  for(nm in all_names) {
    yd <- yr_dim_of(nm)
    expect_true(all(year_slice(out[[nm]], c(4, 5), yd) == 0), info = nm)
    expect_true(all(year_slice(out[[nm]], 1:3, yd) == 1), info = nm)
  } # end name loop
})

test_that("only the requested data types are modified", {
  d <- make_data()
  out <- SPoRC::set_data_indicator_unused(d, unused_years = 3, what = c("Catch", "SrvIdx"))

  expect_true(all(year_slice(out$UseCatch, 3, 2) == 0))
  expect_true(all(year_slice(out$UseSrvIdx, 3, 2) == 0))
  # Everything else keeps its original values.
  for(nm in setdiff(all_names, c("UseCatch", "UseSrvIdx"))) {
    expect_true(all(out[[nm]] == 1), info = nm)
  } # end name loop
})

test_that("population-level indicators are independent of aggregated ones", {
  d <- make_data()
  out <- SPoRC::set_data_indicator_unused(d, unused_years = 2, what = "Catch_pop")
  expect_true(all(year_slice(out$UseCatch_pop, 2, 3) == 0))
  expect_true(all(out$UseCatch == 1))
})

test_that("years outside the model range are ignored", {
  d <- make_data()
  out <- SPoRC::set_data_indicator_unused(d, unused_years = c(4, 99, -1))

  for(nm in all_names) {
    yd <- yr_dim_of(nm)
    expect_true(all(year_slice(out[[nm]], 4, yd) == 0), info = nm)
    expect_true(all(year_slice(out[[nm]], c(1, 2, 3, 5), yd) == 1), info = nm)
  } # end name loop
})

test_that("an entirely out-of-range year vector leaves indicators untouched", {
  d <- make_data()
  out <- SPoRC::set_data_indicator_unused(d, unused_years = c(50, 60))
  for(nm in all_names) {
    expect_true(all(out[[nm]] == 1), info = nm)
  } # end name loop
  # Tagging cohorts are keyed on the same year filter, so none should drop.
  expect_equal(out$n_conv_tag_cohorts, 4)
})

test_that("tag cohorts released in unused years are dropped", {
  d <- make_data(n_cohorts = 4)
  out <- SPoRC::set_data_indicator_unused(d, unused_years = c(3, 4))

  expect_equal(out$n_conv_tag_cohorts, 2)
  expect_equal(nrow(out$conv_tag_release_indicator), 2)
  expect_equal(out$conv_tag_release_indicator[, 2], c(1, 2))
  expect_equal(dim(out$conv_tagged_fish)[1], 2)
  # Cohorts sit on the third axis of the recapture array.
  expect_equal(dim(out$obs_conv_tag_fish_recap)[3], 2)
})

test_that("dropping a middle cohort keeps the surviving release years", {
  d <- make_data(n_cohorts = 4)
  out <- SPoRC::set_data_indicator_unused(d, unused_years = 2)
  expect_equal(out$n_conv_tag_cohorts, 3)
  expect_equal(out$conv_tag_release_indicator[, 2], c(1, 3, 4))
})

test_that("tagging is untouched when conventional tagging is off", {
  d <- make_data(n_cohorts = 4, use_tagging = 0)
  out <- SPoRC::set_data_indicator_unused(d, unused_years = c(3, 4))
  expect_equal(out$n_conv_tag_cohorts, 4)
  expect_equal(dim(out$conv_tagged_fish)[1], 4)
})

test_that("tagging is untouched when conv_tagging is not requested", {
  d <- make_data(n_cohorts = 4)
  out <- SPoRC::set_data_indicator_unused(d, unused_years = c(3, 4), what = "Catch")
  expect_equal(out$n_conv_tag_cohorts, 4)
  expect_equal(dim(out$obs_conv_tag_fish_recap)[3], 4)
})

test_that("array shapes are preserved apart from the tagging cohort axis", {
  d <- make_data()
  out <- SPoRC::set_data_indicator_unused(d, unused_years = c(2, 3))
  for(nm in all_names) {
    expect_equal(dim(out[[nm]]), dim(d[[nm]]), info = nm)
  } # end name loop
})

test_that("withholding every year zeroes all indicators", {
  # This is the peel-everything corner of a retrospective run.
  d <- make_data()
  out <- SPoRC::set_data_indicator_unused(d, unused_years = 1:n_years)
  for(nm in all_names) {
    expect_true(all(out[[nm]] == 0), info = nm)
  } # end name loop
  expect_equal(out$n_conv_tag_cohorts, 0)
})
