library(SPoRC)
library(testthat)

# ── adjacency_mat validation for CTMC movement ───────────────────────────────

# The generator lays the diffusion rate wherever adjacency is 1, then overwrites its
# own diagonal with the negative column sums to conserve abundance. Anything other
# than a 0/1 matrix with a zero diagonal leaves the generator columns summing to
# something other than zero, and the movement matrix then loses abundance instead of
# redistributing it. These checks are the setup-side guard against that.

make_dim_input_list <- function(n_regions = 2) {
  Setup_Mod_Dim(years = 1:5, ages = 1:6, lens = 1,
                n_regions = n_regions, n_sexes = 1,
                n_fish_fleets = 1, n_srv_fleets = 1,
                n_seas = 1, n_pop = 1, natal_region = 1, verbose = FALSE)
}

make_ctmc_dat <- function(n_regions = 2) {
  expand.grid(pop = 1, regions = 1:n_regions, years = 1:5, sexes = 1, ages = 1:6, seas = 1)
}

run_ctmc_setup <- function(adjacency_mat, use_fixed_movement = 0, move_timing = 2, n_regions = 2) {
  Setup_Mod_Movement(input_list = make_dim_input_list(n_regions),
                     use_fixed_movement = use_fixed_movement,
                     do_recruits_move = 0,
                     move_type = 1,
                     adjacency_mat = adjacency_mat,
                     diffusion_formula = ~1,
                     preference_formula = ~0,
                     move_timing = move_timing,
                     ctmc_move_dat = make_ctmc_dat(n_regions))
}

test_that("a non-zero adjacency diagonal is rejected", {
  # diag(1, 2) is the identity, which connects nothing and puts the diffusion rate on
  # the generator diagonal, draining exp(-theta) of the population every time step
  expect_error(run_ctmc_setup(diag(1, 2)), "must have a zero diagonal")
  expect_error(run_ctmc_setup(diag(1, 2)), "region\\(s\\): 1, 2")

  # a single contaminated region is caught and named
  A <- matrix(c(1, 1, 1, 0), 2, 2)
  expect_error(run_ctmc_setup(A), "region\\(s\\): 1$|region\\(s\\): 1\\.")
})

test_that("non-binary and missing adjacency entries are rejected", {
  expect_error(run_ctmc_setup(matrix(c(0, 0.5, 0.5, 0), 2, 2)), "must be 0 \\(not connected\\) or 1")
  expect_error(run_ctmc_setup(matrix(c(0, NA, 1, 0), 2, 2)), "cannot contain NA")
})

test_that("a fully disconnected adjacency matrix is rejected when movement is estimated", {
  expect_error(run_ctmc_setup(matrix(0, 2, 2)), "no off-diagonal connections")

  # adjacency is unused when the movement matrix is fixed, so it is not an error there
  expect_no_error(run_ctmc_setup(matrix(0, 2, 2), use_fixed_movement = 1, move_timing = 0))
})

test_that("a valid adjacency matrix passes and is stored unchanged", {
  A <- matrix(c(0, 1, 1, 0), 2, 2)
  il <- run_ctmc_setup(A)
  expect_identical(il$data$adjacency_mat, A)

  # partial connectivity is legal: region 3 reachable only from region 2
  A3 <- matrix(c(0, 1, 0,
                 1, 0, 1,
                 0, 1, 0), 3, 3, byrow = TRUE)
  expect_no_error(run_ctmc_setup(A3, n_regions = 3))
})

test_that("a valid adjacency matrix yields a conservative generator", {
  # the payoff of the guard: rows of Movement sum to 1 and columns of the generator sum to 0
  A <- matrix(c(0, 1, 1, 0), 2, 2)
  mv <- SPoRC:::Get_Movement(
    move_type = 1, do_recruits_move = 1,
    n_pop = 1, n_regions = 2, n_yrs = 1, n_proj_yrs_devs = 0,
    n_ages = 1, n_sexes = 1, n_seas = 1,
    move_pars = NULL,
    move_devs = array(0, dim = c(1, 2, 1, 1, 1, 1, 1)),
    use_fixed_movement = 0,
    ctmc_move_dat = expand.grid(pop = 1, regions = 1:2, years = 1, seas = 1, ages = 1, sexes = 1),
    preference_formula = ~0, diffusion_formula = ~1,
    log_move_diffusion_pars = array(log(0.3), dim = c(1, 1)),
    move_preference_pars = array(0, dim = c(1, 1)),
    area_r = c(1, 1), adjacency_mat = A, ctmc_diffusion_bounds = 0
  )
  M <- matrix(mv$Movement[1, , , 1, 1, 1, 1], 2, 2)
  Q <- matrix(mv$Mrate[1, , , 1, 1, 1, 1], 2, 2)
  expect_equal(rowSums(M), c(1, 1), tolerance = 1e-10)
  expect_equal(rowSums(Q), c(0, 0), tolerance = 1e-10) # Mrate is stored transposed
  expect_true(all(M > 0)) # both regions genuinely exchange
})
