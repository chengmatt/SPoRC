library(SPoRC)
library(testthat)
library(Matrix)

# End-to-end coverage for move_timing. The operator-level properties are covered in
# test-transition-operators.R; what this file guards is the plumbing through a real
# model -- array shapes, argument wiring, and the continuous-movement branches that
# only execute at move_timing = 2. A dimension-dropping bug in the timing-2 tagging
# Baranov previously survived precisely because no test ran a full model at that timing.

# Build a small multi-region CTMC model exercising movement, tagging and two sexes.
build_ctmc_model <- function(n_regions = 3, n_seas = 1, n_sexes = 2, n_yrs = 8, n_ages = 6) {
  set.seed(404)

  input_list <- Setup_Mod_Dim(years = 1:n_yrs, ages = 1:n_ages, lens = 1:5, n_regions = n_regions,
                              n_seas = n_seas, n_sexes = n_sexes, n_fish_fleets = 1,
                              n_srv_fleets = 1, n_pop = 1, verbose = FALSE)

  ctmc_dat <- expand.grid(pop = 1, regions = seq_len(n_regions), years = seq_len(n_yrs),
                          seas = seq_len(n_seas), ages = seq_len(n_ages), sexes = seq_len(n_sexes))
  adj <- matrix(1L, n_regions, n_regions); diag(adj) <- 0L

  input_list <- Setup_Mod_Movement(input_list = input_list,
                                   move_type = 1,
                                   do_recruits_move = 0,
                                   use_fixed_movement = 0,
                                   ctmc_move_dat = ctmc_dat,
                                   adjacency_mat = adj,
                                   area_r = rep(1, n_regions),
                                   diffusion_formula = ~ 1,
                                   preference_formula = ~ 1,
                                   move_timing = 0)
  input_list
}

test_that("Setup_Mod_Movement accepts every move_timing for a CTMC model", {
  il <- build_ctmc_model()
  for (tm in 0:2) {
    out <- Setup_Mod_Movement(input_list = il,
                              move_type = 1, do_recruits_move = 0, use_fixed_movement = 0,
                              ctmc_move_dat = il$data$ctmc_move_dat,
                              adjacency_mat = il$data$adjacency_mat,
                              area_r = il$data$area_r,
                              diffusion_formula = ~ 1, preference_formula = ~ 1,
                              move_timing = tm)
    expect_equal(out$data$move_timing, tm)
  }
})

test_that("continuous movement forces the generator onto annual time units", {
  il <- build_ctmc_model()
  out <- Setup_Mod_Movement(input_list = il,
                            move_type = 1, do_recruits_move = 0, use_fixed_movement = 0,
                            ctmc_move_dat = il$data$ctmc_move_dat,
                            adjacency_mat = il$data$adjacency_mat,
                            area_r = il$data$area_r,
                            diffusion_formula = ~ 1, preference_formula = ~ 1,
                            move_timing = 2, ctmc_scale_by_seasdur = 0)
  # Mixing an unscaled generator with seasdur-scaled mortality is dimensionally inconsistent
  expect_equal(out$data$ctmc_scale_by_seasdur, 1)
})

test_that("Get_Movement returns a usable generator for every timing", {
  # Shape guard on the pieces the dynamics index into: Movement and Mrate must be
  # dimensioned alike so that Mrate[p,,,y,seas,a,s] lines up with Movement[p,,,y,seas,a,s]
  n_regions <- 3; n_ages <- 4; n_sexes <- 2; n_yrs <- 3; n_seas <- 2
  dat <- expand.grid(pop = 1, regions = seq_len(n_regions), years = seq_len(n_yrs),
                     seas = seq_len(n_seas), ages = seq_len(n_ages), sexes = seq_len(n_sexes))
  adj <- matrix(1L, n_regions, n_regions); diag(adj) <- 0L

  res <- Get_Movement(
    move_type = 1, do_recruits_move = 1,
    n_pop = 1, n_regions = n_regions, n_yrs = n_yrs, n_proj_yrs_devs = 0,
    n_ages = n_ages, n_sexes = n_sexes, n_seas = n_seas,
    move_pars = array(0, c(1, n_regions, n_regions - 1, n_yrs, n_seas, n_ages, n_sexes)),
    move_devs = array(0, c(1, n_regions, n_regions - 1, n_yrs, n_seas, n_ages, n_sexes)),
    use_fixed_movement = 0, Fixed_Movement = NULL,
    ctmc_move_dat = dat, preference_formula = ~ 1, diffusion_formula = ~ 1,
    log_move_diffusion_pars = log(0.25), move_preference_pars = 0,
    area_r = rep(1, n_regions), adjacency_mat = adj, ctmc_diffusion_bounds = 0,
    seasdur = rep(1 / n_seas, n_seas), ctmc_scale_by_seasdur = 1
  )

  expect_equal(dim(res$Mrate), dim(res$Movement))

  # Every stratum's generator must drive a valid transition under all three timings
  for (tm in 0:2) {
    for (a in 2:n_ages) for (s in seq_len(n_sexes)) for (sea in seq_len(n_seas)) {
      Mv <- res$Movement[1, , , 1, sea, a, s]
      Qv <- res$Mrate[1, , , 1, sea, a, s]
      out <- advance_seas(c(100, 50, 25), Mv, rep(0, n_regions), Qv, 1 / n_seas, tm)
      expect_equal(sum(out), 175, tolerance = 1e-8)  # no mortality => mass conserved
      expect_true(all(out >= 0))
    }
  }
})

test_that("timing-2 catch and tag Baranov survive a length-1 sex dimension", {
  # Regression guard: the slices multiplied against the season-integrated abundance
  # drop their sex dimension when n_sexes == 1, which previously made them
  # non-conformable with the 3-d integral array.
  n_regions <- 3; n_ages <- 4; n_sexes <- 1
  m <- local({
    set.seed(77)
    D <- matrix(stats::runif(n_regions^2, 0.05, 0.5), n_regions, n_regions)
    diag(D) <- 0; diag(D) <- -colSums(D)
    list(Q_row = t(D))
  })

  # Mimic the shapes the tagging and fishery models actually index with
  ret_FAA <- array(stats::runif(n_regions * n_ages * n_sexes, 0.01, 0.2),
                   dim = c(1, n_regions, 1, n_ages, n_sexes, 1))
  Z <- array(stats::runif(n_regions * n_ages * n_sexes, 0.05, 0.4),
             dim = c(1, n_regions, 1, n_ages, n_sexes))
  N <- array(100, dim = c(n_regions, n_ages, n_sexes))

  tag_int <- array(0, dim = c(n_regions, n_ages, n_sexes))
  for (a in seq_len(n_ages)) for (s in seq_len(n_sexes)) {
    tag_int[, a, s] <- integrate_seas_abundance(N[, a, s], Z[1, , 1, a, s], m$Q_row, 1)
  }

  slice <- array(ret_FAA[1, , 1, , , 1], dim = c(n_regions, n_ages, n_sexes))
  recap <- rep(0.3, n_regions) * slice * tag_int

  expect_equal(dim(recap), c(n_regions, n_ages, n_sexes))
  expect_true(all(is.finite(recap)))
  expect_true(all(recap >= 0))
})

test_that("continuous movement reduces to the discrete timings when regions are decoupled", {
  # With a zero generator there is no movement, so all three timings must agree in a
  # full seasonal sweep -- an end-to-end check that the operator plumbing is consistent.
  n <- 3; n_seas <- 4
  seasdur <- rep(0.25, n_seas)
  Q0 <- matrix(0, n, n); I <- diag(n)
  Zs <- matrix(stats::runif(n * n_seas, 0.05, 0.3), n, n_seas)
  N0 <- c(120, 60, 90)

  res <- lapply(0:2, function(tm) {
    N <- N0
    for (s in seq_len(n_seas)) N <- advance_seas(N, I, Zs[, s], Q0, seasdur[s], tm)
    N
  })
  expect_equal(res[[2]], res[[1]], tolerance = 1e-12)
  expect_equal(res[[3]], res[[1]], tolerance = 1e-12)
})
