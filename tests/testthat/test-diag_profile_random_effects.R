# A profile with random effects warns that the data components cannot be read on their own, and its
# jnLL_df is the marginal: profiling R0 at the fitted value gives back the fit's own objective.

library(SPoRC)
library(testthat)

data("sgl_rg_dusky_data")

test_that("random effects raise the warning and name the random effects", {
  expect_warning(SPoRC:::check_profile_random(c("ln_RecDevs", "ln_InitDevs")),
                 "ln_RecDevs, ln_InitDevs")
})

test_that("no random effects, no warning", {
  expect_no_warning(SPoRC:::check_profile_random(NULL))
  expect_no_warning(SPoRC:::check_profile_random(character(0)))
})

test_that("with random effects, jnLL_df at the fitted R0 is the fit's marginal, and the laplace term is the rest", {

  input_list <- build_goa_dusky_input(sgl_rg_dusky_data)
  random <- c("ln_RecDevs", "ln_InitDevs")
  fit <- fit_model(input_list$data, input_list$par, input_list$map, random = random,
                   newton_loops = 0, silent = TRUE)
  fit_par <- fit$env$parList()
  R0 <- as.numeric(fit_par$ln_global_R0)

  expect_warning(
    prof <- do_likelihood_profile(input_list$data, fit_par, input_list$map, random, what = "ln_global_R0",
                                  min_val = R0, max_val = R0),
    "Random effects are estimated"
  )

  expect_equal(prof$jnLL_df$value, as.numeric(fit$optim$objective), tolerance = 1e-6)
  expect_equal(prof$laplace_nLL_df$value, as.numeric(fit$optim$objective) - fit$rep$jnLL, tolerance = 1e-6)
  expect_true(abs(prof$laplace_nLL_df$value) > 0.1) # the marginal and jnLL at the mode differ here
  expect_true("Laplace" %in% prof$agg_nLL$type)

})

test_that("without random effects the laplace term is zero and stays out of agg_nLL", {

  input_list <- build_goa_dusky_input(sgl_rg_dusky_data)
  fit <- fit_model(input_list$data, input_list$par, input_list$map, random = NULL,
                   newton_loops = 0, silent = TRUE)
  fit_par <- fit$env$parList()
  R0 <- as.numeric(fit_par$ln_global_R0)

  prof <- do_likelihood_profile(input_list$data, fit_par, input_list$map, NULL, what = "ln_global_R0",
                                min_val = R0, max_val = R0)

  expect_identical(prof$laplace_nLL_df$value, 0)
  expect_false("Laplace" %in% prof$agg_nLL$type)

})
