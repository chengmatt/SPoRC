library(SPoRC)
library(testthat)

# Fleet-specific ageing error. A fishery reading otoliths and a survey reading
# scales do not misclassify the same way, so each fleet may carry its own matrix.
# Both default to the shared AgeingError, and a model written before these
# existed has to come out bit for bit unchanged.

# a row-stochastic smearing matrix: some of each true age is read one age either side
smear <- function(n, p = 0.2) {
  m <- diag(1 - p, n)
  for(i in seq_len(n)) {
    nb <- c(i - 1, i + 1)
    nb <- nb[nb >= 1 & nb <= n]
    m[i, nb] <- m[i, nb] + p / length(nb)
    m[i, i] <- m[i, i] + (1 - sum(m[i, ]))   # keep the row summing to one exactly
  } # end i loop
  m
}

test_that("expand_fleet_ageing_error defaults every fleet to the shared matrix", {
  n_years <- 4; n_ages <- 6
  shared <- array(0, dim = c(n_years, n_ages, n_ages))
  for(i in seq_len(n_years)) shared[i,,] <- smear(n_ages)

  out <- expand_fleet_ageing_error(NULL, shared, 3, "AgeingError_fish")
  expect_equal(dim(out), c(n_years, n_ages, n_ages, 3))
  for(f in 1:3) expect_equal(out[,,,f], shared, tolerance = 1e-12)
})

test_that("expand_fleet_ageing_error accepts the time-invariant and time-varying forms", {
  n_years <- 3; n_ages <- 5; n_fleets <- 2
  shared <- array(0, dim = c(n_years, n_ages, n_ages))
  for(i in seq_len(n_years)) shared[i,,] <- diag(1, n_ages)

  # time-invariant, one matrix per fleet, broadcast over years
  inv <- array(0, dim = c(n_ages, n_ages, n_fleets))
  inv[,,1] <- diag(1, n_ages)
  inv[,,2] <- smear(n_ages, 0.3)
  out <- expand_fleet_ageing_error(inv, shared, n_fleets, "AgeingError_fish")
  expect_equal(dim(out), c(n_years, n_ages, n_ages, n_fleets))
  for(i in seq_len(n_years)) expect_equal(out[i,,,2], inv[,,2], tolerance = 1e-12)

  # time-varying, taken as given
  tv <- array(0, dim = c(n_years, n_ages, n_ages, n_fleets))
  for(i in seq_len(n_years)) for(f in seq_len(n_fleets)) tv[i,,,f] <- smear(n_ages, 0.1 * i)
  out2 <- expand_fleet_ageing_error(tv, shared, n_fleets, "AgeingError_fish")
  expect_equal(out2, tv, tolerance = 1e-12)
})

test_that("expand_fleet_ageing_error refuses shapes that would misalign the comps", {
  n_years <- 3; n_ages <- 5
  shared <- array(0, dim = c(n_years, n_ages, n_ages))
  for(i in seq_len(n_years)) shared[i,,] <- diag(1, n_ages)

  wrong_fleets <- array(diag(1, n_ages), dim = c(n_ages, n_ages, 3))
  expect_error(expand_fleet_ageing_error(wrong_fleets, shared, 2, "AgeingError_fish"), "3 fleets")

  # a different observed age count per fleet cannot work: the observed comp arrays
  # carry one age dimension shared across fleets
  wrong_obs <- array(0, dim = c(n_ages, n_ages - 1, 2))
  for(f in 1:2) wrong_obs[,,f] <- diag(1, n_ages)[, 2:n_ages]
  expect_error(expand_fleet_ageing_error(wrong_obs, shared, 2, "AgeingError_fish"), "matching AgeingError")

  # a row that does not sum to one is reported through the setup messages rather
  # than rejected, because the likelihood renormalizes the expectation after the
  # multiply and models supplied such a matrix long before it was ever checked
  bad_rows <- array(0, dim = c(n_ages, n_ages, 2))
  for(f in 1:2) bad_rows[,,f] <- diag(1, n_ages)
  bad_rows[1, 1, 2] <- 0.5
  messages_list <<- character(0)
  out <- expand_fleet_ageing_error(bad_rows, shared, 2, "AgeingError_fish")
  expect_equal(out[1,,,2], bad_rows[,,2], tolerance = 1e-12)   # passed through untouched
  expect_true(any(grepl("sum to neither", messages_list)))

  expect_error(expand_fleet_ageing_error(diag(1, n_ages), shared, 2, "AgeingError_fish"), "must be NULL")
})

test_that("fleet_ageing_error falls back for a model fitted before the arrays existed", {
  n_years <- 2; n_ages <- 4
  shared <- array(0, dim = c(n_years, n_ages, n_ages))
  for(i in seq_len(n_years)) shared[i,,] <- smear(n_ages)

  old_data <- list(n_fish_fleets = 2, n_srv_fleets = 1)   # no AgeingError_fish at all
  out <- fleet_ageing_error(old_data, shared, "fish")
  expect_equal(dim(out), c(n_years, n_ages, n_ages, 2))
  for(f in 1:2) expect_equal(out[,,,f], shared, tolerance = 1e-12)

  own <- array(0, dim = c(n_years, n_ages, n_ages, 2))
  for(i in seq_len(n_years)) { own[i,,,1] <- diag(1, n_ages); own[i,,,2] <- smear(n_ages, 0.4) }
  new_data <- list(n_fish_fleets = 2, n_srv_fleets = 1, AgeingError_fish = own)
  expect_equal(fleet_ageing_error(new_data, shared, "fish"), own, tolerance = 1e-12)
})

test_that("a fleet's own ageing error changes only that fleet's composition likelihood", {
  input <- objective_fixture_input()
  base <- evaluate_input(input)
  n_ages <- length(input$data$ages)
  n_yrs <- length(input$data$years)

  # handing every fleet the shared matrix must reproduce the original exactly,
  # which is the backwards-compatibility guarantee
  same <- input
  same$data$AgeingError_fish <- fleet_ageing_error(input$data, input$data$AgeingError, "fish")
  same$data$AgeingError_srv <- fleet_ageing_error(input$data, input$data$AgeingError, "srv")
  expect_equal(evaluate_input(same)$fn(), base$fn(), tolerance = 1e-12)

  # smearing the survey's ages moves the survey age composition likelihood and
  # leaves the fishery's alone
  srv_smeared <- same
  for(i in seq_len(n_yrs)) srv_smeared$data$AgeingError_srv[i,,,1] <- smear(n_ages, 0.3)
  out <- evaluate_input(srv_smeared)
  expect_false(isTRUE(all.equal(out$fn(), base$fn())))
  expect_false(isTRUE(all.equal(sum(out$rep$SrvAgeComps_nLL), sum(base$rep$SrvAgeComps_nLL))))
  expect_equal(sum(out$rep$FishAgeComps_nLL), sum(base$rep$FishAgeComps_nLL), tolerance = 1e-10)
})

# The operating model simulates compositions through the per-fleet matrices, and
# the estimation model has to be handed the same ones back. Simulate_Pop_Static
# returns a curated list rather than the whole environment, so anything it does
# not name is lost, and simulation_data_to_SPoRC then quietly gives every fleet
# the shared matrix. That fails silently: the self-test still runs, it just
# estimates against an ageing error the data were never generated with.

test_that("Simulate_Pop_Static carries the per-fleet matrices through to the estimation model", {
  skip_if_not(exists("caal_make_om"), "helper-selftest_caal.R not loaded")

  om <- caal_make_om()

  expect_true("AgeingError_fish" %in% names(om),
              info = "Simulate_Pop_Static must name AgeingError_fish in its returned list")
  expect_true("AgeingError_srv" %in% names(om),
              info = "Simulate_Pop_Static must name AgeingError_srv in its returned list")
  # every fleet defaults to the shared matrix, so the slabs agree with it
  expect_equal(dim(om$AgeingError_fish)[1:3], dim(om$AgeingError)[1:3])
  expect_equal(om$AgeingError_fish[,,,1,1], om$AgeingError[,,,1], tolerance = 1e-12)
  expect_equal(om$AgeingError_srv[,,,1,1], om$AgeingError[,,,1], tolerance = 1e-12)

  # and they survive the peel into the estimation model's data list
  sd <- simulation_data_to_SPoRC(sim_env = om, y = dim(om$AgeingError)[1], sim = 1)
  expect_false(is.null(sd$AgeingError_fish))
  expect_false(is.null(sd$AgeingError_srv))
  expect_equal(dim(sd$AgeingError_fish)[1], dim(om$AgeingError)[1])
  expect_equal(sd$AgeingError_fish[,,,1], sd$AgeingError, tolerance = 1e-12)
})
