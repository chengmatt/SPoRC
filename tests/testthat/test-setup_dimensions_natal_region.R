library(SPoRC)
library(testthat)

# Setup_Sim_Dim and Setup_Mod_Dim both infer natal_region when it is not
# supplied. Getting that mapping wrong silently misassigns recruitment to the
# wrong region in a multi-population model, so each inference branch and its
# error case is covered here.

sim_dim <- function(n_pop = 1, n_regions = 1, natal_region = NULL, n_sexes = 2, ...) {
  SPoRC::Setup_Sim_Dim(n_sims = 1, n_yrs = 5, n_pop = n_pop, n_regions = n_regions,
                       natal_region = natal_region, n_ages = 10, n_sexes = n_sexes,
                       n_fish_fleets = 1, n_srv_fleets = 1, ...)
}

mod_dim <- function(n_pop = 1, n_regions = 1, natal_region = NULL, lens = NULL, ...) {
  SPoRC::Setup_Mod_Dim(years = 1:5, ages = 1:10, lens = lens, n_pop = n_pop,
                       natal_region = natal_region, n_regions = n_regions,
                       n_sexes = 2, n_fish_fleets = 1, n_srv_fleets = 1, ...)
}

test_that("Setup_Sim_Dim sends every population to region 1 in a single-region model", {
  out <- sim_dim(n_pop = 3, n_regions = 1)
  expect_equal(out$natal_region, rep(1, 3))
})

test_that("Setup_Sim_Dim maps populations one-to-one when counts match", {
  out <- sim_dim(n_pop = 3, n_regions = 3)
  expect_equal(out$natal_region, 1:3)
})

test_that("Setup_Sim_Dim homes a single population to region 1 across many regions", {
  out <- sim_dim(n_pop = 1, n_regions = 4)
  expect_equal(out$natal_region, 1)
})

test_that("Setup_Sim_Dim requires natal_region when the mapping is ambiguous", {
  # Two populations across three regions has no defensible default.
  expect_error(sim_dim(n_pop = 2, n_regions = 3), "natal_region must be specified")
})

test_that("Setup_Sim_Dim keeps an explicit natal_region", {
  out <- sim_dim(n_pop = 2, n_regions = 3, natal_region = c(3, 1))
  expect_equal(out$natal_region, c(3, 1))
})

test_that("Setup_Sim_Dim rejects more than two sexes", {
  expect_error(sim_dim(n_sexes = 3), "number of sexes")
})

test_that("Setup_Sim_Dim accepts one or two sexes", {
  expect_equal(sim_dim(n_sexes = 1)$n_sexes, 1)
  expect_equal(sim_dim(n_sexes = 2)$n_sexes, 2)
})

test_that("Setup_Sim_Dim derives season durations that sum to one", {
  expect_equal(sim_dim()$seasdur, 1)
  out <- sim_dim(n_seas = 4)
  expect_equal(out$seasdur, rep(0.25, 4))
  expect_equal(sum(out$seasdur), 1)
})

test_that("Setup_Sim_Dim sets the initialisation iterations from the age range", {
  expect_equal(sim_dim()$init_iter, 10 * 10)
})

test_that("Setup_Mod_Dim sends every population to region 1 in a single-region model", {
  out <- mod_dim(n_pop = 3, n_regions = 1)
  expect_equal(out$data$natal_region, rep(1, 3))
})

test_that("Setup_Mod_Dim maps populations one-to-one when counts match", {
  out <- mod_dim(n_pop = 3, n_regions = 3)
  expect_equal(out$data$natal_region, 1:3)
})

test_that("Setup_Mod_Dim homes a single population to region 1 across many regions", {
  out <- mod_dim(n_pop = 1, n_regions = 4)
  expect_equal(out$data$natal_region, 1)
})

test_that("Setup_Mod_Dim requires natal_region when the mapping is ambiguous", {
  expect_error(mod_dim(n_pop = 2, n_regions = 3), "natal_region must be specified")
})

test_that("Setup_Mod_Dim keeps an explicit natal_region", {
  out <- mod_dim(n_pop = 2, n_regions = 3, natal_region = c(2, 3))
  expect_equal(out$data$natal_region, c(2, 3))
})

test_that("Setup_Mod_Dim returns the empty par and map sublists downstream setup expects", {
  out <- mod_dim()
  expect_named(out, c("data", "par", "map", "verbose", "store_config", "version"))
  expect_equal(out$par, list())
  expect_equal(out$map, list())
})

test_that("Setup_Mod_Dim stores lens as 1 when none are supplied", {
  # Length bins are optional, and downstream array construction needs a size.
  expect_equal(mod_dim()$data$lens, 1)
  expect_equal(mod_dim(lens = seq(10, 50, 10))$data$lens, seq(10, 50, 10))
})

test_that("Setup_Mod_Dim prints its dimension summary only when verbose", {
  expect_silent(mod_dim(verbose = FALSE))
  msgs <- capture_messages(mod_dim(verbose = TRUE))
  expect_true(any(grepl("Number of Years: 5", msgs)))
  expect_true(any(grepl("Number of Regions: 1", msgs)))
  expect_true(any(grepl("Number of Age Bins: 10", msgs)))
  expect_true(any(grepl("Number of Fishery Fleets: 1", msgs)))
})

test_that("Setup_Mod_Dim reports one duration message per season when verbose", {
  msgs <- capture_messages(mod_dim(n_seas = 4, verbose = TRUE))
  expect_equal(sum(grepl("^Duration of season", msgs)), 4)
})

test_that("Setup_Mod_Dim records its own call only when store_config is set", {
  expect_null(mod_dim()$config)
  cfg <- mod_dim(n_pop = 2, n_regions = 2, store_config = TRUE)$config
  expect_equal(cfg$Setup_Mod_Dim$n_pop, 2)
  expect_equal(cfg$Setup_Mod_Dim$n_regions, 2)
})
