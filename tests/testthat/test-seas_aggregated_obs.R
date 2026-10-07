# Data set to "aggSeas" are fit once a year against the season total rather than once a season.
# Checks the parser, the one observation a year rule, that a fleet reporting annually still fishes
# every season, and that the likelihood reads the year total.

library(SPoRC)
library(testthat)

# Specification parsing ------------------------------------------------------

test_that("parse_seas_agg_spec resolves values, recycles and rejects bad input", {

  expect_equal(SPoRC:::parse_seas_agg_spec(NULL, "x", 3), rep(0L, 3))
  expect_equal(SPoRC:::parse_seas_agg_spec("aggSeas", "x", 3), rep(1L, 3))
  expect_equal(SPoRC:::parse_seas_agg_spec(c("spltSeas", "aggSeas"), "x", 2), c(0L, 1L))

  # a fitted model passes its resolved codes straight back to an operating model
  expect_equal(SPoRC:::parse_seas_agg_spec(c(0, 1), "x", 2), c(0L, 1L))

  expect_error(SPoRC:::parse_seas_agg_spec("annual", "Catch_seas_Type", 2), "Catch_seas_Type has invalid value")
  expect_error(SPoRC:::parse_seas_agg_spec(c("aggSeas", "aggSeas"), "Catch_seas_Type", 3), "3 fleets")
})

test_that("check_seas_agg_use allows one observed season a year and refuses more", {

  use <- array(0, dim = c(2, 3, 4, 2))
  use[1, 1, 2, 1] <- 1
  use[1, 2, 3, 1] <- 1
  expect_null(SPoRC:::check_seas_agg_use(use, c(1, 0), "UseCatch"))

  use[1, 1, 3, 1] <- 1 # a second season in the same year
  expect_error(SPoRC:::check_seas_agg_use(use, c(1, 0), "UseCatch"), "more than one season")

  # a fleet left seasonal is never checked
  expect_null(SPoRC:::check_seas_agg_use(use, c(0, 0), "UseCatch"))

  # population-specific arrays put population first
  use_pop <- array(0, dim = c(2, 2, 3, 4, 2))
  use_pop[1, 1, 1, 2, 1] <- 1
  expect_null(SPoRC:::check_seas_agg_use(use_pop, c(1, 0), "UseCatch_pop"))
  use_pop[1, 1, 1, 4, 1] <- 1
  expect_error(SPoRC:::check_seas_agg_use(use_pop, c(1, 0), "UseCatch_pop"), "more than one season")
})

# Objective ------------------------------------------------------------------

# One two season model, reported both ways, so the difference is only the reporting
seas_agg_pair <- function() {

  sim <- seasonal_M_sim()
  seasonal_input <- seasonal_M_input(sim, list(M_spec = "fix", Fixed_natmort = seasonal_M_fixed(seasonal_M_cfg$M)))

  annual_input <- seasonal_input
  obs <- annual_input$data$ObsCatch
  use <- annual_input$data$UseCatch

  # the year's catch becomes one observation sitting in season one
  year_total <- apply(obs, c(1, 2, 4), sum)
  obs[] <- 0
  obs[, , 1, ] <- year_total
  use[] <- 0
  use[, , 1, ] <- 1

  annual_input$data$ObsCatch <- obs
  annual_input$data$UseCatch <- use
  annual_input$data$Catch_seas_Type <- 1L

  list(seasonal = seasonal_input, annual = annual_input)
}

test_that("a fleet reporting once a year is still fished in every season", {

  pair <- seas_agg_pair()

  fit <- fit_model(data = pair$annual$data,
                   parameters = pair$annual$par,
                   mapping = pair$annual$map,
                   random = NULL,
                   do_optim = FALSE,
                   silent = TRUE)

  # season two has no observation of its own, and would be a closure without the setting
  expect_true(all(fit$rep$Fmort[, , 2, 1] > 0))
  expect_equal(dim(fit$rep$Fmort)[3], seasonal_M_cfg$n_seas)
})

test_that("the aggregated catch likelihood reads the year total", {

  pair <- seas_agg_pair()

  fit <- fit_model(data = pair$annual$data,
                   parameters = pair$annual$par,
                   mapping = pair$annual$map,
                   random = NULL,
                   do_optim = FALSE,
                   silent = TRUE)

  pred_year <- apply(fit$rep$PredCatch, c(2, 3, 5), sum) # summed over populations and seasons
  obs_year <- pair$annual$data$ObsCatch[, , 1, 1]
  sigma <- exp(fit$rep$ln_sigmaC[, , 1, 1])

  expect_equal(as.vector(fit$rep$Catch_nLL[, , 1, 1]),
               as.vector(-stats::dnorm(log(obs_year), log(pred_year[, , 1]), sigma, log = TRUE)),
               tolerance = 1e-10)

  # the other seasons hold no likelihood of their own
  expect_true(all(fit$rep$Catch_nLL[, , -1, ] == 0))
})

test_that("aggregating the catch leaves the population dynamics alone", {

  pair <- seas_agg_pair()

  seasonal_fit <- fit_model(data = pair$seasonal$data, parameters = pair$seasonal$par,
                            mapping = pair$seasonal$map, random = NULL, do_optim = FALSE, silent = TRUE)

  annual_fit <- fit_model(data = pair$annual$data, parameters = pair$annual$par,
                          mapping = pair$annual$map, random = NULL, do_optim = FALSE, silent = TRUE)

  # the same parameters, so the numbers at age and the catch they imply cannot move
  expect_equal(annual_fit$rep$NAA, seasonal_fit$rep$NAA, tolerance = 1e-12)
  expect_equal(annual_fit$rep$SSB, seasonal_fit$rep$SSB, tolerance = 1e-12)
  expect_equal(apply(annual_fit$rep$PredCatch, c(2, 3, 5), sum),
               apply(seasonal_fit$rep$PredCatch, c(2, 3, 5), sum), tolerance = 1e-12)
})

test_that("the setting is inert when the model has one season", {

  seas_agg <- SPoRC:::parse_seas_agg_spec("aggSeas", "Catch_seas_Type", 1)
  use <- array(1, dim = c(1, 5, 1, 1))
  expect_null(SPoRC:::check_seas_agg_use(use, seas_agg, "UseCatch"))

  # nothing to collapse when the year is one season long
  true_arr <- array(stats::runif(5), dim = c(1, 5, 1, 1, 1))
  obs_arr <- true_arr
  out <- SPoRC:::collapse_seas_obs(true_arr, obs_arr, array(0.1, dim = c(1, 5, 1, 1)),
                                   seas_agg, 0, y = 1, sim = 1, n_seas = 1,
                                   n_regions = 1, n_fleets = 1)
  expect_identical(out$true, true_arr)
  expect_identical(out$obs, obs_arr)
})

# Simulation -----------------------------------------------------------------

test_that("collapse_seas_obs puts the year in season one and draws once from that total", {

  set.seed(11)
  true_arr <- array(0, dim = c(1, 2, 3, 1, 1)) # region, year, season, fleet, sim
  true_arr[1, 1, , 1, 1] <- c(2, 3, 5)
  true_arr[1, 2, , 1, 1] <- c(1, 1, 1)
  obs_arr <- true_arr

  out <- SPoRC:::collapse_seas_obs(true_arr, obs_arr, array(0, dim = c(1, 2, 3, 1)),
                                   seas_agg = 1L, like_type = 0L, y = 1, sim = 1,
                                   n_seas = 3, n_regions = 1, n_fleets = 1)

  expect_equal(out$true[1, 1, 1, 1, 1], 10) # the year's total
  expect_true(all(out$true[1, 1, -1, 1, 1] == 0))
  expect_equal(out$obs[1, 1, 1, 1, 1], 10) # zero error, so the draw is the total
  expect_true(all(out$obs[1, 1, -1, 1, 1] == 0))

  # a year that was not asked for is untouched
  expect_equal(out$true[1, 2, , 1, 1], c(1, 1, 1))
})

test_that("collapse_seas_obs writes the year into the season the fit holds it in", {

  true_arr <- array(0, dim = c(1, 1, 3, 1, 1)) # region, year, season, fleet, sim
  true_arr[1, 1, , 1, 1] <- c(2, 3, 5)

  out <- SPoRC:::collapse_seas_obs(true_arr, true_arr, array(0, dim = c(1, 1, 3, 1)),
                                   seas_agg = 1, like_type = 0, y = 1, sim = 1,
                                   n_seas = 3, n_regions = 1, n_fleets = 1, slot = 2)

  expect_equal(as.vector(out$true[1, 1, , 1, 1]), c(0, 10, 0))
  expect_equal(as.vector(out$obs[1, 1, , 1, 1]), c(0, 10, 0))
})

test_that("simulate_comps draws a year total from the numbers summed over its seasons", {

  # pop, region, year, season, age, sex, fleet, sim
  numbers <- array(stats::runif(2 * 3 * 4), dim = c(1, 1, 2, 3, 4, 1, 1, 1))
  year_summed <- numbers
  year_summed[1, 1, 1, 1, , 1, 1, 1] <- apply(numbers[1, 1, 1, , , 1, 1, 1], 2, sum)

  comp_args <- list(r = 1, y = 1, f = 1, seas = 1, sim = 1,
                    ISS = array(100, dim = c(1, 2, 3, 1, 1, 1)),
                    AgeingError = array(diag(4), dim = c(2, 4, 4, 1)),
                    comp_like = 0,
                    comp_type = matrix(1, nrow = 2, ncol = 1),
                    n_sexes = 1,
                    n_regions = 1,
                    n_cat = 4,
                    Obs = array(0, dim = c(1, 2, 3, 4, 1, 1, 1)),
                    age_or_len = 0)

  set.seed(5)
  over_seasons <- do.call(SPoRC:::simulate_comps, c(comp_args, list(Exp = numbers, exp_seas = 1:3)))
  set.seed(5)
  presummed <- do.call(SPoRC:::simulate_comps, c(comp_args, list(Exp = year_summed)))

  expect_equal(over_seasons, presummed)
  expect_equal(sum(over_seasons[1, 1, 1, , 1, 1, 1]), 100)
})
