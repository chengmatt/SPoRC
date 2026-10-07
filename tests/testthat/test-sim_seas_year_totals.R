# Data sources reported once a year in a seasonal model. The operating model draws each year total from
# the year's numbers into the season the fit holds it in, and a discard fraction is the year's discards
# over the year's catch on both sides, never a sum of seasonal or population fractions.

library(SPoRC)
library(testthat)

test_that("each year total is drawn from the year's numbers into the season the fit holds it in", {

  il <- year_total_input()
  obj <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)
  rep <- obj$rep
  d <- year_total_dims
  terms <- c("Catch_nLL", "Discard_nLL", "Discard_pop_nLL", "FishAgeComps_pop_nLL", "FishAgeComps_discard_nLL", "SrvAgeComps_pop_nLL")

  sim_file <- tempfile(fileext = ".RDS")
  set.seed(3)
  res <- suppressWarnings(suppressMessages(simulation_self_test(
    data = obj$data,
    parameters = obj$parameters,
    mapping = obj$mapping,
    random = NULL,
    rep = rep,
    sd_rep = list(par.fixed = obj$par, par.random = NULL),
    n_sims = 1,
    newton_loops = 0,
    output_path = sim_file,
    what = terms
  )))
  om <- readRDS(sim_file)

  # the catch, summed over populations and seasons, sits in season two
  pred_catch <- apply(rep$PredCatch[,,,,1, drop = FALSE], c(2, 3), sum)
  expect_true(all(om$ObsCatch[,,1,1,1] == 0))
  expect_lt(max(abs(om$TrueCatch[,,2,1,1] / pred_catch - 1)), 1e-12)

  # the year's discard fraction is the estimation model's, from the same sums
  pred_discard <- array(0, dim = c(d$n_regions, d$n_yrs))
  for(r in seq_len(d$n_regions)) for(y in seq_len(d$n_yrs))
    pred_discard[r,y] <- get_discard_pred(rep$PredDiscard, rep$CAA, rep$DAA, rep$dmr, obj$data$WAA_fish, 2, seq_len(d$n_pop), r, y, seq_len(d$n_seas), 1)
  expect_true(all(om$ObsDiscard[,,1,1,1] == 0))
  expect_lt(max(abs(om$TrueDiscard[,,2,1,1] / pred_discard - 1)), 1e-12)
  expect_lt(max(abs(om$TrueDiscard_pop[,,,,1,1] / rep$PredDiscard[,,,,1] - 1)), 1e-12) # population by season

  # each composition sits in its own season alone, drawn at that season's sample size
  expect_true(all(om$ObsFishAgeComps_pop[,,,1,,,,] == 0))
  expect_true(all(apply(om$ObsFishAgeComps_pop[,,,2,,,1,1], c(1, 2, 3), sum) == 100))
  expect_true(all(om$ObsFishAgeComps_discard[,,2,,,,] == 0))
  expect_true(all(apply(om$ObsFishAgeComps_discard[,,1,,,1,1], c(1, 2), sum) == 80))
  expect_true(all(om$ObsSrvAgeComps_pop[,,,1,,,,] == 0))
  expect_true(all(apply(om$ObsSrvAgeComps_pop[,,,2,,,1,1], c(1, 2, 3), sum) == 60))

  # and the refit reads every one of them
  for(term_name in terms) {
    expect_true(all(is.finite(res[[term_name]])), info = term_name)
    expect_true(sum(res[[term_name]]) != 0, info = term_name)
  } # end term_name loop
})

test_that("a discard fraction over populations or seasons is a ratio of sums", {

  # two populations each discarding 30 percent of their catch discard 30 percent together
  n_ages <- 3
  caa <- array(7, dim = c(2, 1, 1, 2, n_ages, 1, 1)) # pop, region, year, season, age, sex, fleet
  daa <- array(3, dim = dim(caa))
  dmr <- array(1, dim = c(1, 1, 2, 1))
  waa <- array(1, dim = dim(caa))
  pred <- array(0.3, dim = c(2, 1, 1, 2, 1))
  expect_equal(get_discard_pred(pred, caa, daa, dmr, waa, 2, 1:2, 1, 1, 1, 1), 0.3)
  expect_equal(get_discard_pred(pred, caa, daa, dmr, waa, 2, 1:2, 1, 1, 1:2, 1), 0.3)

  # unequal seasons weigh by their catch: 18 of 60 fish discarded in season one and 1 of 40 in season
  # two make 19 of 100 for the year, where averaging the two fractions would give 0.1625
  caa[,,,2,,,] <- 39 / (2 * n_ages)
  daa[,,,2,,,] <- 1 / (2 * n_ages)
  expect_equal(get_discard_pred(pred, caa, daa, dmr, waa, 2, 1:2, 1, 1, 1:2, 1), 19 / 100)

  # one population in one season is the fraction the model always had
  expect_equal(get_discard_pred(pred, caa, daa, dmr, waa, 2, 1, 1, 1, 1, 1), 1 - sum(caa[1,1,1,1,,,]) / sum(caa[1,1,1,1,,,] + daa[1,1,1,1,,,]))

  # numbers and weight still add up
  expect_equal(get_discard_pred(pred, caa, daa, dmr, waa, 0, 1:2, 1, 1, 1:2, 1), sum(pred))
})

test_that("seas_agg_slot_list reads the fit's season and refuses two in one year", {

  use <- array(0, dim = c(2, 4, 3, 1)) # region, year, season, fleet
  use[,1,2,1] <- 1
  use[,2,3,1] <- 1
  data <- list(n_fish_fleets = 1, UseCatch = use, Catch_seas_Type = 1)
  slot <- seas_agg_slot_list(data, n_yrs = 6, platform = "fish")$Catch

  expect_equal(as.vector(slot), c(2, 3, 1, 1, 1, 1)) # no observation in years 3 and 4, and the last fitted year after

  data$UseCatch[2,1,2,1] <- 0
  data$UseCatch[2,1,3,1] <- 1 # region two holds year one's total in another season
  expect_error(seas_agg_slot_list(data, n_yrs = 4, platform = "fish"), "draws a fleet's year total into one season")
})

test_that("the operating model refuses a 2d logistic normal that is not joint by sex", {

  expect_error(check_sim_2d_comp(c(0, 4), matrix(c(1, 1, 1, 1), 2, 2), "comp_fishage_like", "FishAgeComps_Type"), "joint by sex")
  expect_error(check_sim_2d_comp(7, matrix(c(2, 0), 2, 1), "comp_srvage_like", "SrvAgeComps_Type"), "joint by sex")
  expect_null(check_sim_2d_comp(c(4, 1), matrix(c(2, 2, 1, 1), 2, 2), "comp_fishage_like", "FishAgeComps_Type"))
})
