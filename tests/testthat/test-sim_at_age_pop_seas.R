# Population-specific at-age data and at-age year totals in the operating model. Each is drawn from the
# quantity the estimation model predicts, so at the fit the operating model's truth is the prediction.

library(SPoRC)
library(testthat)

# One self test at the starting values, keeping the operating model it drew from
pop_seas_self_test <- function(n_sims = 2) {

  built <- suppressWarnings(suppressMessages(pop_seas_at_age_input()))
  il <- built$input
  obj <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)

  sim_file <- tempfile(fileext = ".RDS")
  set.seed(11)
  res <- suppressWarnings(suppressMessages(simulation_self_test(
    data = obj$data,
    parameters = obj$parameters,
    mapping = obj$mapping,
    random = NULL,
    rep = obj$rep,
    sd_rep = list(par.fixed = obj$par, par.random = NULL),
    n_sims = n_sims,
    newton_loops = 0,
    output_path = sim_file,
    what = c("CatchAA_nLL", "CatchAA_pop_nLL", "SrvIdxAA_pop_nLL")
  )))

  c(built, list(obj = obj, res = res, om = readRDS(sim_file)))
}

test_that("the operating model's at-age truth is the estimation model's prediction at the fit", {

  st <- pop_seas_self_test()
  om <- st$om
  rep <- st$obj$rep
  rel_diff <- function(a, b) max(abs(a - b) / abs(b))

  # each population drawn from its own catch at age, never summed over populations
  on <- st$use_caa_pop == 1
  expect_lt(rel_diff(om$TrueCatchAA_pop[,,,,,,,1][on], rep$PredCatchAA_pop[on]), 1e-12)

  # a year total, summed over seasons and populations, drawn into the season its flags name
  on <- st$use_caa == 1
  expect_lt(rel_diff(om$TrueCatchAA[,,,,,,1][on], rep$PredCatchAA[on]), 1e-12)
  expect_true(all(om$ObsCatchAA[,,2,,,2,] == 0))

  # the survey's year total by population sits in season two and nowhere else
  on <- st$use_saa_pop == 1
  expect_lt(rel_diff(om$TrueSrvIdxAA_pop[,,,,,,,1][on], rep$PredSrvIdxAA_pop[on]), 1e-12)
  expect_true(all(om$ObsSrvIdxAA_pop[,,,1,,,,] == 0))
  expect_true(all(om$ObsSrvIdxAA_pop[,,,2,,,,] > 0))

  # the refit read the draws rather than the data it was conditioned on
  expect_true(all(is.finite(st$res$CatchAA_pop_nLL)))
  expect_false(isTRUE(all.equal(sum(st$res$SrvIdxAA_pop_nLL[,,,,,,,1]), sum(rep$SrvIdxAA_pop_nLL))))
  expect_false(isTRUE(all.equal(sum(st$res$CatchAA_nLL[,,,,,,1]), sum(rep$CatchAA_nLL))))
})

test_that("simulation_data_to_SPoRC returns the at-age data with the operating model's use flags", {

  st <- pop_seas_self_test()
  d <- pop_seas_dims
  extracted <- simulation_data_to_SPoRC(st$om, y = d$n_yrs, sim = 2)

  expect_equal(extracted$UseCatchAA_pop, st$use_caa_pop)
  expect_equal(extracted$UseCatchAA, st$use_caa)
  expect_equal(extracted$UseSrvIdxAA_pop, st$use_saa_pop)
  expect_equal(as.vector(extracted$ObsCatchAA_pop), as.vector(st$om$ObsCatchAA_pop[,,,,,,,2]))
  expect_equal(dim(extracted$ObsSrvIdxAA_pop_SE), dim(st$use_saa_pop))

  # a data source the operating model never drew comes back NULL
  expect_null(extracted$ObsDiscardAA)
  expect_null(extracted$UseSrvIdxAA)

  # a shorter year range keeps the first years
  early <- simulation_data_to_SPoRC(st$om, y = 3, sim = 1)
  expect_equal(dim(early$ObsCatchAA_pop)[3], 3) # population, region, year
  expect_equal(early$UseCatchAA_pop, st$use_caa_pop[,,1:3,,,,,drop = FALSE])
})

test_that("check_seas_agg_use reads the season dim of an at-age array", {

  # region, year, season, observed age, sex, fleet, with one season observed a year
  use <- array(0, dim = c(1, 3, 2, 4, 1, 2))
  use[1,,1,,1,1] <- 1
  expect_null(check_seas_agg_use(use, c(1, 0), "UseCatchAA"))

  # every age observed in both seasons of a year
  use[1,2,2,,1,1] <- 1
  expect_error(check_seas_agg_use(use, c(1, 0), "UseCatchAA"), "more than one season")

  # the population-specific array puts population first
  use_pop <- array(0, dim = c(2, 1, 3, 2, 4, 1, 2))
  use_pop[,1,,2,,1,1] <- 1
  expect_null(check_seas_agg_use(use_pop, c(1, 0), "UseCatchAA_pop"))
  use_pop[2,1,3,1,,1,1] <- 1
  expect_error(check_seas_agg_use(use_pop, c(1, 0), "UseCatchAA_pop"), "more than one season")
})
