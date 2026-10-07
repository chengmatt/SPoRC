# Correlated at-age residuals in the operating model: each form draws at the correlations the estimation model builds from
# the same parameters, iid draws are unchanged, and the self test hands the operating model the fit's settings.

n_yrs_c <- 6
n_obs_c <- 5

# standardized residuals of one region, season, sex and fleet, drawn for many replicates
at_age_corr_draws <- function(code, rho_age, rho_year = 0, us_pars = NULL, use = NULL, n_sims = 4000) {
  env <- new.env()
  env$n_sims <- n_sims
  env$UseCatchAA <- if(is.null(use)) array(1, dim = c(1, n_yrs_c, 1, n_obs_c, 1, 1)) else use
  env$AgeObsCorr_catch <- code
  env$trans_rho_catch <- array(atanh(rho_age), dim = c(1, 1, 1)) # rho_trans is tanh
  env$trans_rho_catch_year <- array(atanh(rho_year), dim = c(1, 1, 1))
  env$trans_rho_catch_us <- array(if(is.null(us_pars)) 0 else us_pars, dim = c(n_obs_c * (n_obs_c - 1) / 2, 1, 1, 1))
  draw_sim_at_age_corr(env)
  env$CatchAA_std_resid[1,,1,,1,1,] # year by observed age by replicate
}

test_that("an ar1 across ages draws at rho to the age gap, a gap included, and leaves years independent", {
  set.seed(1)
  z <- at_age_corr_draws(1, 0.6)
  expect_equal(cor(z[3,1,], z[3,2,]), 0.6, tolerance = 0.05)
  expect_equal(cor(z[3,1,], z[3,3,]), 0.36, tolerance = 0.1)
  expect_lt(abs(cor(z[3,2,], z[4,2,])), 0.05)
  expect_equal(sd(z[3,2,]), 1, tolerance = 0.05)

  use_gap <- array(1, dim = c(1, n_yrs_c, 1, n_obs_c, 1, 1))
  use_gap[1,,1,3,1,1] <- 0 # age three never observed
  z <- at_age_corr_draws(1, 0.6, use = use_gap)
  expect_equal(cor(z[3,2,], z[3,4,]), 0.36, tolerance = 0.1) # two ages apart across the gap
  expect_true(all(z[,3,] == 0))
})

test_that("the unstructured form draws at the correlation build_us_corr gives", {
  us_pars <- c(0.5, -0.3, 0.2, 0.1, 0.4, -0.2, 0.3, 0.1, 0.2, -0.1)
  set.seed(2)
  z <- at_age_corr_draws(2, 0, us_pars = us_pars)
  expect_lt(max(abs(cor(t(z[2,,])) - build_us_corr(us_pars, n_obs_c))), 0.06)
})

test_that("the separable ar1 draws the observed years by ages together", {
  set.seed(3)
  z <- at_age_corr_draws(3, 0.5, rho_year = 0.7)
  expect_equal(cor(z[2,2,], z[3,2,]), 0.7, tolerance = 0.05)
  expect_equal(cor(z[2,2,], z[2,3,]), 0.5, tolerance = 0.05)
  expect_equal(cor(z[2,2,], z[3,3,]), 0.35, tolerance = 0.1)
})

test_that("an iid fleet draws exactly what it drew before the correlated forms existed", {
  numbers <- array(c(100, 50, 25, 10), dim = c(1, 1, 4, 1))
  use <- array(1, dim = c(1, 4, 1))
  ln_sigma <- array(log(0.3), dim = c(4, 1))
  for(oe in c(0, 1)) {
    set.seed(4)
    drawn <- sim_at_age_cell(numbers, numbers, use, array(0, dim = c(1, 4, 1)), ln_sigma, 1, 0, 0, FALSE, 1, bias_correct_oe = oe)
    set.seed(4)
    by_hand <- vapply(1:4, function(a) numbers[1,1,a,1] * exp(stats::rnorm(1, if(oe == 1) -0.5 * 0.3^2 else 0, 0.3)), numeric(1))
    expect_identical(as.vector(drawn$obs), by_hand)
  } # end oe loop
})

test_that("the self test hands the operating model the fit's at-age correlation", {
  d <- c(1, 20, 1, 5, 1, 1) # region, year, season, observed age, sex, fleet
  il <- suppressWarnings(suppressMessages(build_at_age(ObsCatchAA = array(1e3, dim = d), UseCatchAA = array(1, dim = d), AgeObsCorr_catch = "1dar1")))
  expect_true(any(!is.na(il$map$trans_rho_catch))) # estimated, so it reaches the self test from the fit
  il$par$trans_rho_catch[] <- 0.4
  obj <- suppressWarnings(fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE))
  captured <- new.env()
  testthat::with_mocked_bindings(
    Simulate_Pop_Static = function(sim_list, ...) { captured$sim_list <- sim_list; stop("captured") },
    try(suppressWarnings(suppressMessages(
      simulation_self_test(data = obj$data, parameters = obj$parameters, mapping = obj$mapping, random = NULL, rep = obj$rep,
                           sd_rep = list(par.fixed = obj$par, par.random = NULL), n_sims = 1, newton_loops = 0, what = "SSB"))), silent = TRUE),
    .package = "SPoRC")
  expect_equal(captured$sim_list$AgeObsCorr_catch, 1)
  expect_equal(as.vector(captured$sim_list$trans_rho_catch), rep(0.4, length(il$par$trans_rho_catch)))
  expect_true(all(captured$sim_list$AgeObsCorr_srv_idx == 0)) # the survey's own setting, iid here
})
