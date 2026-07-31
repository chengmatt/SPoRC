library(SPoRC)
library(testthat)

# Each prior and penalty below is a single line of SPoRC_rtmb that no other test switches
# on. The pattern throughout is the same: evaluate the model with the term off, evaluate it
# again with the term on, and check that jnLL moved by exactly the value worked out by hand
# from the prior's own mean and standard deviation. Nothing is optimised, so both
# evaluations sit at the same parameter values and the difference is the term itself.

build <- function(...) suppressWarnings(suppressMessages(objective_fixture_input(...)))

# every likelihood other than the one under test must be untouched by switching a prior on,
# since a prior adds to jnLL without feeding back into the population dynamics
expect_other_likelihoods_unchanged <- function(with_term, without_term, except) {
  components <- setdiff(grep("_nLL$", names(without_term$rep), value = TRUE), except)
  for(nm in components) {
    expect_equal(with_term$rep[[nm]], without_term$rep[[nm]],
                 info = paste0("switching on ", except, " changed ", nm))
  }
}


test_that("the R0 prior adds a normal penalty on ln_global_R0 and nothing else", {
  r0_prior <- data.frame(pop = 1, mu = 4, sd = 0.2)

  off <- evaluate_input(build())
  on  <- evaluate_input(build(rec = list(use_r0_prior = 1, r0_prior = r0_prior)))

  expected <- -stats::dnorm(log(5), log(r0_prior$mu), r0_prior$sd, log = TRUE)

  expect_equal(off$rep$R0_nLL, 0)
  expect_equal(as.numeric(on$rep$R0_nLL), expected, tolerance = 1e-8)
  expect_equal(on$rep$jnLL - off$rep$jnLL, expected, tolerance = 1e-8)
  expect_other_likelihoods_unchanged(on, off, except = "R0_nLL")
})


test_that("the fishery catchability prior adds a normal penalty on ln_fish_q and nothing else", {
  fish_q_prior <- data.frame(region = 1, block = 1, fleet = 1, mu = 0.5, sd = 0.3)

  off <- evaluate_input(build())
  on  <- evaluate_input(build(fishsel = list(Use_fish_q_prior = 1, fish_q_prior = fish_q_prior)))

  ln_fish_q <- build()$par$ln_fish_q[1, 1, 1]
  expected <- -stats::dnorm(ln_fish_q, log(fish_q_prior$mu), fish_q_prior$sd, log = TRUE)

  expect_equal(off$rep$fish_q_nLL, 0)
  expect_equal(as.numeric(on$rep$fish_q_nLL), expected, tolerance = 1e-8)
  expect_equal(on$rep$jnLL - off$rep$jnLL, expected, tolerance = 1e-8)
  expect_other_likelihoods_unchanged(on, off, except = "fish_q_nLL")
})


test_that("the retained selectivity prior adds a normal penalty on ret_fixed_sel_pars and nothing else", {
  ret_selex_prior <- data.frame(region = 1, par = 1, block = 1, sex = 1, fleet = 1, mu = 2, sd = 0.4)

  off <- evaluate_input(build())
  on  <- evaluate_input(build(fishsel = list(Use_ret_selex_prior = 1, ret_selex_prior = ret_selex_prior)))

  ret_par <- build()$par$ret_fixed_sel_pars[1, 1, 1, 1, 1]
  expected <- -stats::dnorm(ret_par, log(ret_selex_prior$mu), ret_selex_prior$sd, log = TRUE)

  # sel_nLL also carries the smoothness penalties, which are zero here, so the prior is the
  # whole of it
  expect_equal(off$rep$sel_nLL, 0)
  expect_equal(as.numeric(on$rep$sel_nLL), expected, tolerance = 1e-8)
  expect_equal(on$rep$jnLL - off$rep$jnLL, expected, tolerance = 1e-8)
  expect_other_likelihoods_unchanged(on, off, except = "sel_nLL")
})


test_that("the discard mortality rate penalty applies to every fished year, with or without discard data", {
  # dmr is identified through total mortality, so a deviation is estimated and
  # penalized in every fished cell. Discard observations are dropped from the back
  # half of the series while catch continues, so a penalty keyed on discard data
  # would cover only the first 15 years and this test would catch it
  use_discard <- build()$data$UseDiscard
  use_discard[, 16:30, , ] <- 0

  # the deviations are estimated in both models so that only the penalty itself differs
  catch_f_off <- list(Use_dmr_pen = 0, dmr_dev_spec = "est_all", UseDiscard = use_discard)
  catch_f_on  <- list(Use_dmr_pen = 1, dmr_dev_spec = "est_all", UseDiscard = use_discard)

  off <- evaluate_input(build(catch_f = catch_f_off))
  on  <- evaluate_input(build(catch_f = catch_f_on))

  input <- build(catch_f = catch_f_on)
  devs <- input$par$logit_dmr_devs
  sigma_dmr <- exp(input$par$ln_sigma_dmr)
  use_catch <- input$data$UseCatch

  # one penalty per fished region-year-season-fleet cell
  expected_cells <- array(0, dim = dim(devs))
  for(f in seq_len(dim(devs)[4])) for(y in seq_len(dim(devs)[2])) {
    for(r in seq_len(dim(devs)[1])) for(seas in seq_len(dim(devs)[3])) {
      if(use_catch[r, y, seas, f] == 1) {
        expected_cells[r, y, seas, f] <- -stats::dnorm(devs[r, y, seas, f], 0, sigma_dmr[r, seas, f], log = TRUE)
      }
    }
  }
  expected <- sum(expected_cells)

  expect_gt(sum(use_catch), sum(use_discard)) # the two sets genuinely differ
  expect_gt(expected, 0)
  expect_equal(sum(off$rep$dmr_nLL), 0)
  expect_equal(as.numeric(on$rep$dmr_nLL), as.numeric(expected_cells), tolerance = 1e-8)

  # every estimated deviation is penalized, and no deviation fixed at zero contributes
  expect_equal(as.vector(on$rep$dmr_nLL != 0), !is.na(input$map$logit_dmr_devs))

  # dmr_nLL enters jnLL weighted by Wt_D
  expect_equal(on$rep$jnLL - off$rep$jnLL, input$data$Wt_D * expected, tolerance = 1e-8)
  expect_other_likelihoods_unchanged(on, off, except = "dmr_nLL")
})


test_that("jnLL still decomposes with these priors and penalties switched on", {
  input <- build(rec = list(use_r0_prior = 1, r0_prior = data.frame(pop = 1, mu = 4, sd = 0.2)),
                 catch_f = list(Use_dmr_pen = 1, dmr_dev_spec = "est_all"),
                 fishsel = list(Use_fish_q_prior = 1,
                                fish_q_prior = data.frame(region = 1, block = 1, fleet = 1, mu = 0.5, sd = 0.3),
                                Use_ret_selex_prior = 1,
                                ret_selex_prior = data.frame(region = 1, par = 1, block = 1, sex = 1,
                                                             fleet = 1, mu = 2, sd = 0.4)))
  model <- evaluate_input(input)

  expect_jnLL_decomposes(model, label = "all priors on")

  expect_true(as.numeric(model$rep$R0_nLL) != 0)
  expect_true(as.numeric(model$rep$fish_q_nLL) != 0)
  expect_true(as.numeric(model$rep$sel_nLL) != 0)
  expect_true(sum(model$rep$dmr_nLL) != 0)
})
