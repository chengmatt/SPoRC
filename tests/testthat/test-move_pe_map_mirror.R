# Checks map_move_devs takes cells out of the movement penalty by exactly their own contribution
# (1e-12), the way a dsem takes a series over, and leaves the rest alone.

move_pen_setup <- function(seed = 9, n_yrs = 8) {
  set.seed(seed)
  d <- c(1, 2, 2, n_yrs, 1, 3, 1) # pop, from, to, year, season, age, sex
  list(d = d, n_yrs = n_yrs,
       move_devs = array(rnorm(prod(d), 0, 0.3), dim = d),
       PE_pars = array(log(0.25), dim = c(1, 2, 1, 3, 1)), # pop, from, season, age, sex
       adjacency = matrix(1, 2, 2))
}

move_pen <- function(s, mirror) Get_move_PE_loglik(cont_vary_movement = "iid_y", PE_pars = s$PE_pars, move_devs = s$move_devs, map_move_devs = mirror,
                                                   do_recruits_move = 1, adjacency_collapsed = s$adjacency, move_type = 0)

test_that("a mirror keeping every cell matches the penalty by hand", {

  s <- move_pen_setup()
  full <- array(1, dim = s$d)
  # iid_y keys on year only, so the penalty reads pop 1, season 1, age 1, sex 1 across from, to and year
  by_hand <- sum(dnorm(as.vector(s$move_devs[1,,,,1,1,1]), 0, 0.25, log = TRUE))
  expect_equal(move_pen(s, full), by_hand, tolerance = 1e-12)

})

test_that("blanked cells drop exactly their own contribution", {

  s <- move_pen_setup()
  full <- array(1, dim = s$d)
  mirror <- full
  mirror[1,1,2,3:5,1,1,1] <- NA # region 1 to region 2, years 3 to 5, the cells the penalty reads

  dropped <- sum(dnorm(s$move_devs[1,1,2,3:5,1,1,1], 0, 0.25, log = TRUE))
  expect_equal(move_pen(s, mirror), move_pen(s, full) - dropped, tolerance = 1e-12)

  # a blank on a cell the penalty never reads (age 2 under a year-only model) changes nothing
  elsewhere <- full; elsewhere[1,1,2,3,1,2,1] <- NA
  expect_equal(move_pen(s, elsewhere), move_pen(s, full), tolerance = 1e-12)

  expect_equal(move_pen(s, array(NA_real_, dim = s$d)), 0, tolerance = 1e-12) # everything out

})
