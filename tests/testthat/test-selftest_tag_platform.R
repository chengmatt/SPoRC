# A tagging model set up without a release platform keeps the default the operating model reads, in the
# estimation model's data and in a self test built from a data list that predates the default.

library(SPoRC)
library(testthat)

test_that("a tagging model without a release platform stores the default, and the self test falls back to it", {

  om <- build_om(move_timing = 0)
  em <- build_em(om, 0)
  default <- default_tag_release_platform(em$data$conv_tag_release_indicator)
  expect_equal(dim(default), c(nrow(em$data$conv_tag_release_indicator), 2))
  expect_true(all(default[, "platform"] == "survey") && all(default[, "fleet"] == "1"))
  expect_false(is.null(em$data$conv_tag_release_platform)) # the setup stored one

  # a data list from before the default, as a saved fit may hold; full attribution never reads the platform, so it fits
  em$data$conv_tag_release_platform <- NULL
  obj <- fit_model(em$data, em$par, em$map, random = NULL, do_optim = FALSE, silent = TRUE)
  seen <- NULL
  testthat::with_mocked_bindings(
    Simulate_Pop_Static = function(sim_list, ...) { seen <<- sim_list; stop("captured") },
    try(suppressWarnings(suppressMessages(simulation_self_test(data = obj$data, parameters = em$par, mapping = em$map, random = NULL,
                                                               rep = obj$rep, sd_rep = list(par.fixed = obj$par, par.random = NULL),
                                                               n_sims = 1, newton_loops = 0, what = "SSB"))), silent = TRUE),
    .package = "SPoRC"
  )
  expect_equal(seen$conv_tag_release_platform, default)

})
