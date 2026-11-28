# Setup refuses a fixed movement matrix or a spawning movement matrix whose rows do not sum to one, or that is only
# partly given, naming the argument, and builds when every row sums to one.

# Setup ----------------------------------------------------------------------

# movement fractions [pop, from, to, year, season, age, sex] for the sweep model with two regions and one fleet
movement_rows_input <- function(stay_region_1, leave_region_1) {
  Movement <- array(0, dim = c(1, 2, 2, sweep_dims$n_yrs, 1, sweep_dims$n_ages, sweep_dims$n_sexes))
  Movement[1,1,1,,,,] <- stay_region_1
  Movement[1,1,2,,,,] <- leave_region_1
  Movement[1,2,1,,,,] <- 0.2
  Movement[1,2,2,,,,] <- 0.8
  Movement
}

# spawning movement [pop, from, to, year, age, sex] for two populations, one born in each region
spawn_rows_input <- function(returning, staying_away) {
  spawn_move <- array(0, dim = c(2, 2, 2, sweep_dims$n_yrs, sweep_dims$n_ages, sweep_dims$n_sexes))
  for(p in 1:2) {
    home <- p
    away <- 3 - p
    spawn_move[p,home,home,,,] <- 1
    spawn_move[p,away,home,,,] <- returning
    spawn_move[p,away,away,,,] <- staying_away
  } # end p loop
  spawn_move
}

# Tests ----------------------------------------------------------------------

test_that("a fixed movement matrix is refused when a row does not sum to one", {

  dims <- list(n_regions = 2, n_fish_fleets = 1)

  # every row sums to one
  expect_no_error(suppressWarnings(sweep_input(dims = dims, move = list(Fixed_Movement = movement_rows_input(0.7, 0.3)), stop_after = "move")))

  # region 1's fish sum to 1.2
  expect_error(suppressWarnings(sweep_input(dims = dims, move = list(Fixed_Movement = movement_rows_input(0.7, 0.5)), stop_after = "move")),
               "Fixed_Movement has rows that do not sum to 1.*population 1, region 1, which sums to 1.2")

  # the matrix is not read when movement is estimated, and the default is an identity matrix
  expect_no_error(suppressWarnings(sweep_input(dims = dims, move = list(use_fixed_movement = 0, Fixed_Movement = NA), stop_after = "move")))
  expect_no_error(suppressWarnings(sweep_input(dims = dims, move = list(Fixed_Movement = NA), stop_after = "move")))

})

test_that("a spawning movement matrix is refused when a row does not sum to one", {

  dims <- list(n_pop = 2, natal_region = c(1, 2), n_regions = 2, n_fish_fleets = 1)

  # 60 pct return to spawn and 40 pct stay away
  expect_no_error(suppressWarnings(sweep_input(dims = dims, rec = list(rec_dd = "local", sgl_seas_spawning_movement = spawn_rows_input(0.6, 0.4)),
                                               stop_after = "rec")))

  # 60 pct return and 60 pct stay away
  expect_error(suppressWarnings(sweep_input(dims = dims, rec = list(rec_dd = "local", sgl_seas_spawning_movement = spawn_rows_input(0.6, 0.6)),
                                            stop_after = "rec")),
               "sgl_seas_spawning_movement has rows that do not sum to 1.*which sums to 1.2")

})

test_that("rows within rounding of one are accepted", {

  # a third each, written to seven decimals
  Movement <- array(0.3333333, dim = c(1, 3, 3, sweep_dims$n_yrs, 1, sweep_dims$n_ages, sweep_dims$n_sexes))
  expect_no_error(SPoRC:::check_movement_rows(Movement, "Fixed_Movement"))

  # written to three decimals a row sums to 0.999, which loses a fish in a thousand every step
  expect_error(SPoRC:::check_movement_rows(array(0.333, dim = dim(Movement)), "Fixed_Movement"), "sums to 0.999")

})

test_that("a movement matrix with only some values missing is refused, and one not given takes the default", {

  dims <- list(n_regions = 2, n_fish_fleets = 1)

  # one fraction missing from a matrix that is used
  with_missing <- movement_rows_input(0.7, 0.3)
  with_missing[1,2,1,3,1,2,1] <- NA
  expect_error(suppressWarnings(sweep_input(dims = dims, move = list(Fixed_Movement = with_missing), stop_after = "move")),
               "Fixed_Movement has 1 missing values")

  # the same matrix is not read when movement is estimated
  expect_no_error(suppressWarnings(sweep_input(dims = dims, move = list(use_fixed_movement = 0, Fixed_Movement = with_missing), stop_after = "move")))

  # not given at all: no movement
  il <- suppressWarnings(sweep_input(dims = dims, move = list(Fixed_Movement = NA), stop_after = "move"))
  expect_equal(il$data$Fixed_Movement[1,,,1,1,1,1], diag(2))

  # spawning movement, one fraction missing
  dims_pop <- list(n_pop = 2, natal_region = c(1, 2), n_regions = 2, n_fish_fleets = 1)
  spawn_missing <- spawn_rows_input(0.6, 0.4)
  spawn_missing[2,1,2,4,3,1] <- NA
  expect_error(suppressWarnings(sweep_input(dims = dims_pop, rec = list(rec_dd = "local", sgl_seas_spawning_movement = spawn_missing), stop_after = "rec")),
               "sgl_seas_spawning_movement has 1 missing values")

  # not given at all: every fish spawns in its natal region
  il <- suppressWarnings(sweep_input(dims = dims_pop, rec = list(rec_dd = "local"), stop_after = "rec"))
  expect_equal(il$data$sgl_seas_spawning_movement[2,,,1,1,1], matrix(c(0, 0, 1, 1), nrow = 2))

})
