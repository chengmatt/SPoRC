# A "par_normal" selectivity prior row is a normal on the parameter's natural scale, exp(par), where a
# "par" row is a normal on the log parameter.

prior_rows <- function(type, mu, sd) {
  data.frame(region = 1, fleet = 1, block = 1, sex = 1, par = 1, mu = mu, sd = sd, type = type)
}

test_that("a par_normal row is a normal density on exp(par)", {
  pars <- array(log(c(3.7, 2.6)), dim = c(1, 2, 1, 1, 1)) # region, par, block, sex, fleet
  blocks <- array(1, dim = c(1, 5, 1))
  for(mu in c(0, 10)) {
    got <- SPoRC:::get_selex_prior(prior_rows("par_normal", mu, 3), pars, NULL, NULL, 0, blocks)
    expect_equal(as.numeric(got), -dnorm(3.7, mu, 3, log = TRUE), tolerance = 1e-12)
  }
  # the default row type is unchanged
  got <- SPoRC:::get_selex_prior(prior_rows("par", 4, 0.5), pars, NULL, NULL, 0, blocks)
  expect_equal(as.numeric(got), -dnorm(log(3.7), log(4), 0.5, log = TRUE), tolerance = 1e-12)
})

test_that("a model built with a par_normal prior reports it in the selectivity penalty", {
  build <- function(use_prior) {
    il <- suppressMessages(sweep_input(dims = list(n_regions = 1, n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = 1),
                                       fishsel = list(Use_fish_selex_prior = use_prior, fish_selex_prior = prior_rows("par_normal", 10, 3))))
    il$par$fish_fixed_sel_pars[1,1,1,1,1] <- log(4)
    suppressWarnings(fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE))$rep$sel_nLL
  }
  expect_equal(build(1) - build(0), -dnorm(4, 10, 3, log = TRUE), tolerance = 1e-10)
})

test_that("an unknown type and a non-positive sd are refused", {
  build <- function(prior) {
    suppressMessages(sweep_input(dims = list(n_regions = 1, n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = 1),
                                 fishsel = list(Use_fish_selex_prior = 1, fish_selex_prior = prior)))
  }
  expect_error(build(prior_rows("par_natural", 10, 3)), "must be 'par', 'par_normal' or 'value'")
  expect_error(build(prior_rows("par_normal", 10, 0)), "finite and positive")
})
