library(SPoRC)
library(testthat)

# Population-specific discard observations. Every other fixture supplies discards only in
# aggregate, so the Discard_pop_nLL block of SPoRC_rtmb, and the Wt_Discard_pop term it
# feeds into jnLL, are otherwise never reached.

build <- function(...) suppressWarnings(suppressMessages(objective_fixture_input(...)))

sim_discards <- function() {
  sim_obj <- objective_fixture_sim()
  simulation_data_to_SPoRC(sim_env = sim_obj, y = sim_obj$n_years, sim = 1)
}

with_pop_discards <- function() {
  sim_data <- sim_discards()
  build(catch_f = list(ObsDiscard_pop = sim_data$ObsDiscard_pop,
                       UseDiscard_pop = sim_data$UseDiscard_pop))
}


test_that("population-specific discards produce a lognormal likelihood per observation", {
  aggregated_only <- evaluate_input(build())
  with_pop <- evaluate_input(with_pop_discards())

  input <- with_pop_discards()
  use_pop <- input$data$UseDiscard_pop
  obs_pop <- input$data$ObsDiscard_pop
  ln_sigma <- input$data$ln_sigmaD_pop
  if(is.null(ln_sigma)) ln_sigma <- input$par$ln_sigmaD_pop

  expect_gt(sum(use_pop), 0)

  # with no population-specific discard data the whole block is skipped
  expect_equal(sum(aggregated_only$rep$Discard_pop_nLL), 0)

  # lognormal on the observation against the predicted discard for that population
  pred <- with_pop$rep$PredDiscard
  expected <- array(0, dim = dim(use_pop))
  idx <- which(use_pop == 1, arr.ind = TRUE)
  for(i in seq_len(nrow(idx))) {
    cell <- idx[i, ]
    expected[cell[1], cell[2], cell[3], cell[4], cell[5]] <-
      -stats::dnorm(log(obs_pop[cell[1], cell[2], cell[3], cell[4], cell[5]]),
                    log(pred[cell[1], cell[2], cell[3], cell[4], cell[5]]),
                    exp(ln_sigma[cell[1], cell[2], cell[3], cell[4], cell[5]]), log = TRUE)
  }

  expect_equal(as.numeric(with_pop$rep$Discard_pop_nLL), as.numeric(expected), tolerance = 1e-8)
})


test_that("population-specific discards enter jnLL through their own weight", {
  input <- with_pop_discards()
  baseline <- evaluate_input(input)

  # doubling the weight adds one more copy of what the term already contributed, which is
  # only true if the term is in the jnLL sum with that weight and no other
  doubled_data <- input$data
  doubled_data$Wt_Discard_pop <- doubled_data$Wt_Discard_pop * 2
  doubled <- fit_model(doubled_data, input$par, input$map, random = NULL,
                       silent = TRUE, do_optim = FALSE)

  contribution <- sum(input$data$Wt_Discard_pop * baseline$rep$Discard_pop_nLL)
  expect_true(contribution != 0)
  expect_equal(doubled$rep$jnLL - baseline$rep$jnLL, contribution, tolerance = 1e-8)

  expect_jnLL_decomposes(baseline, label = "population-specific discards")
})
