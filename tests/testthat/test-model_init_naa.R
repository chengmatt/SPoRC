library(testthat)
library(RTMB)

## Scenario 1: single pop/region/sex/season, UNFISHED
make_unfished_single <- function(n_ages = 4, M = 0.2, R0 = 1000) {
  n_pop <- 1; n_regions <- 1; n_sexes <- 1; n_seas <- 1; n_fish_fleets <- 1

  list(
    n_regions = n_regions,
    n_pop = n_pop,
    n_sexes = n_sexes,
    n_ages = n_ages,
    n_seas = n_seas,
    n_fish_fleets = n_fish_fleets,
    seasdur = 1,
    rec_seas_prop = matrix(1, nrow = n_pop, ncol = n_seas),
    natmort = array(M, dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes)),
    init_F = array(0, dim = c(n_regions, n_seas, n_fish_fleets)),
    dmr = array(0, dim = c(n_regions, n_seas, n_fish_fleets)),
    fish_sel = array(0.5, dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes, n_fish_fleets)),
    ret_sel = array(1, dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes, n_fish_fleets)),
    R0_r = matrix(R0, nrow = n_pop, ncol = n_regions),
    sexratio = array(1, dim = c(n_pop, n_regions, n_sexes)),
    Movement = array(1, dim = c(n_pop, n_regions, n_regions, n_seas, n_ages, n_sexes)),
    do_recruits_move = 0,
    ln_InitDevs = array(0, dim = c(n_pop, n_regions, n_ages - 1))
  )
}

expected_unfished_naa <- function(n_ages, M, R0) {
  naa <- numeric(n_ages)
  naa[1] <- R0
  for (a in 2:(n_ages - 1)) naa[a] <- R0 * exp(-M * (a - 1))
  naa[n_ages] <- naa[n_ages - 1] * exp(-M) / (1 - exp(-M))
  naa
}

test_that("init_age_strc = 1 (scalar geometric series) matches the hand-derived unfished closed form", {

  d <- make_unfished_single(n_ages = 4, M = 0.2, R0 = 1000)
  out <- Get_Init_NAA(
    init_age_strc = 1,
    init_iter = 0,
    n_regions = d$n_regions,
    n_pop = d$n_pop,
    n_sexes = d$n_sexes,
    n_ages = d$n_ages,
    n_seas = d$n_seas,
    n_fish_fleets = d$n_fish_fleets,
    seasdur = d$seasdur,
    rec_seas_prop = d$rec_seas_prop,
    natmort = d$natmort,
    init_F = d$init_F,
    dmr = d$dmr,
    fish_sel = d$fish_sel,
    ret_sel = d$ret_sel,
    R0_r = d$R0_r,
    sexratio = d$sexratio,
    Movement = d$Movement,
    do_recruits_move = d$do_recruits_move,
    ln_InitDevs = d$ln_InitDevs
  )

  expected <- expected_unfished_naa(4, 0.2, 1000)
  expect_equal(as.numeric(out[1, 1, , 1]), expected, tolerance = 1e-8)
})

test_that("init_age_strc = 2 (matrix geometric series) reduces to the scalar solution when n_regions = 1", {

  d <- make_unfished_single(n_ages = 4, M = 0.2, R0 = 1000)
  out <- Get_Init_NAA(
    init_age_strc = 2,
    init_iter = 0,
    n_regions = d$n_regions,
    n_pop = d$n_pop,
    n_sexes = d$n_sexes,
    n_ages = d$n_ages,
    n_seas = d$n_seas,
    n_fish_fleets = d$n_fish_fleets,
    seasdur = d$seasdur,
    rec_seas_prop = d$rec_seas_prop,
    natmort = d$natmort,
    init_F = d$init_F,
    dmr = d$dmr,
    fish_sel = d$fish_sel,
    ret_sel = d$ret_sel,
    R0_r = d$R0_r,
    sexratio = d$sexratio,
    Movement = d$Movement,
    do_recruits_move = d$do_recruits_move,
    ln_InitDevs = d$ln_InitDevs
  )

  expected <- expected_unfished_naa(4, 0.2, 1000)
  expect_equal(as.numeric(out[1, 1, , 1]), expected, tolerance = 1e-8)
})

test_that("init_age_strc = 3 (hybrid) reduces to the scalar solution when n_regions = 1", {

  d <- make_unfished_single(n_ages = 4, M = 0.2, R0 = 1000)
  out <- Get_Init_NAA(
    init_age_strc = 3,
    init_iter = 0,
    n_regions = d$n_regions,
    n_pop = d$n_pop,
    n_sexes = d$n_sexes,
    n_ages = d$n_ages,
    n_seas = d$n_seas,
    n_fish_fleets = d$n_fish_fleets,
    seasdur = d$seasdur,
    rec_seas_prop = d$rec_seas_prop,
    natmort = d$natmort,
    init_F = d$init_F,
    dmr = d$dmr,
    fish_sel = d$fish_sel,
    ret_sel = d$ret_sel,
    R0_r = d$R0_r,
    sexratio = d$sexratio,
    Movement = d$Movement,
    do_recruits_move = d$do_recruits_move,
    ln_InitDevs = d$ln_InitDevs
  )

  expected <- expected_unfished_naa(4, 0.2, 1000)
  expect_equal(as.numeric(out[1, 1, , 1]), expected, tolerance = 1e-8)
})

test_that("init_age_strc = 0 (iterative) converges to the scalar closed form given enough iterations", {

  d <- make_unfished_single(n_ages = 4, M = 0.2, R0 = 1000)
  out <- Get_Init_NAA(
    init_age_strc = 0,
    init_iter = 50,
    n_regions = d$n_regions,
    n_pop = d$n_pop,
    n_sexes = d$n_sexes,
    n_ages = d$n_ages,
    n_seas = d$n_seas,
    n_fish_fleets = d$n_fish_fleets,
    seasdur = d$seasdur,
    rec_seas_prop = d$rec_seas_prop,
    natmort = d$natmort,
    init_F = d$init_F,
    dmr = d$dmr,
    fish_sel = d$fish_sel,
    ret_sel = d$ret_sel,
    R0_r = d$R0_r,
    sexratio = d$sexratio,
    Movement = d$Movement,
    do_recruits_move = d$do_recruits_move,
    ln_InitDevs = d$ln_InitDevs
  )

  expected <- expected_unfished_naa(4, 0.2, 1000)
  expect_equal(as.numeric(out[1, 1, , 1]), expected, tolerance = 1e-4)
})

test_that("init_age_strc = 0 with too few iterations has NOT yet converged (sanity check on the convergence test itself)", {

  d <- make_unfished_single(n_ages = 4, M = 0.2, R0 = 1000)
  out <- Get_Init_NAA(
    init_age_strc = 0,
    init_iter = 1,
    n_regions = d$n_regions,
    n_pop = d$n_pop,
    n_sexes = d$n_sexes,
    n_ages = d$n_ages,
    n_seas = d$n_seas,
    n_fish_fleets = d$n_fish_fleets,
    seasdur = d$seasdur,
    rec_seas_prop = d$rec_seas_prop,
    natmort = d$natmort,
    init_F = d$init_F,
    dmr = d$dmr,
    fish_sel = d$fish_sel,
    ret_sel = d$ret_sel,
    R0_r = d$R0_r,
    sexratio = d$sexratio,
    Movement = d$Movement,
    do_recruits_move = d$do_recruits_move,
    ln_InitDevs = d$ln_InitDevs
  )

  expected <- expected_unfished_naa(4, 0.2, 1000)
  # with only 1 iteration the plus group in particular should NOT have converged
  expect_false(isTRUE(all.equal(as.numeric(out[1, 1, , 1]), expected, tolerance = 1e-4)))
})

## Scenario 2: single pop/region/sex/season, FISHED (retained + discard F)
test_that("init_age_strc = 1 correctly decomposes fishing mortality into retained + discard-mortality components", {

  n_pop <- 1; n_regions <- 1; n_sexes <- 1; n_ages <- 3; n_seas <- 1; n_fish_fleets <- 1
  M <- 0.2; R0 <- 1000; Finit <- 2.0; ret <- 0.8; dmr_val <- 0.3

  natmort <- array(M, dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes))
  init_F <- array(Finit, dim = c(n_regions, n_seas, n_fish_fleets))
  dmr <- array(dmr_val, dim = c(n_regions, n_seas, n_fish_fleets))

  fish_sel <- array(0, dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes, n_fish_fleets))
  fish_sel[, , , 2, , ] <- 1 # age 2
  fish_sel[, , , 3, , ] <- 1 # plus group
  ret_sel <- array(ret, dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes, n_fish_fleets))

  out <- Get_Init_NAA(
    init_age_strc = 1,
    init_iter = 0,
    n_regions = n_regions,
    n_pop = n_pop,
    n_sexes = n_sexes,
    n_ages = n_ages,
    n_seas = n_seas,
    n_fish_fleets = n_fish_fleets,
    seasdur = 1,
    rec_seas_prop = matrix(1, nrow = n_pop, ncol = n_seas),
    natmort = natmort,
    init_F = init_F,
    dmr = dmr,
    fish_sel = fish_sel,
    ret_sel = ret_sel,
    R0_r = matrix(R0, nrow = n_pop, ncol = n_regions),
    sexratio = array(1, dim = c(n_pop, n_regions, n_sexes)),
    Movement = array(1, dim = c(n_pop, n_regions, n_regions, n_seas, n_ages, n_sexes)),
    do_recruits_move = 0,
    ln_InitDevs = array(0, dim = c(n_pop, n_regions, n_ages - 1))
  )

  F_age <- Finit * (ret + (1 - ret) * dmr_val) # F for ages with sel = 1 (both age 2 and the plus group here)
  Z1 <- M                # age 1 unselected, F = 0
  Z2 <- M + F_age
  Z3 <- M + F_age        # plus group, same selectivity as age 2 in this setup

  N1 <- R0
  N2 <- N1 * exp(-Z1)
  N3 <- N2 * exp(-Z2) / (1 - exp(-Z3))

  expect_equal(as.numeric(out[1, 1, , 1]), c(N1, N2, N3), tolerance = 1e-8)
})

## Sex ratio handling
test_that("Get_Init_NAA() allocates age-1 recruits across sexes according to sexratio", {

  n_pop <- 1; n_regions <- 1; n_sexes <- 2; n_ages <- 3; n_seas <- 1; n_fish_fleets <- 1
  M <- 0.2; R0 <- 1000

  sexratio <- array(0, dim = c(n_pop, n_regions, n_sexes))
  sexratio[, , 1] <- 0.6
  sexratio[, , 2] <- 0.4

  out <- Get_Init_NAA(
    init_age_strc = 1,
    init_iter = 0,
    n_regions = n_regions,
    n_pop = n_pop,
    n_sexes = n_sexes,
    n_ages = n_ages,
    n_seas = n_seas,
    n_fish_fleets = n_fish_fleets,
    seasdur = 1,
    rec_seas_prop = matrix(1, nrow = n_pop, ncol = n_seas),
    natmort = array(M, dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes)),
    init_F = array(0, dim = c(n_regions, n_seas, n_fish_fleets)),
    dmr = array(0, dim = c(n_regions, n_seas, n_fish_fleets)),
    fish_sel = array(0.5, dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes, n_fish_fleets)),
    ret_sel = array(1, dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes, n_fish_fleets)),
    R0_r = matrix(R0, nrow = n_pop, ncol = n_regions),
    sexratio = sexratio,
    Movement = array(1, dim = c(n_pop, n_regions, n_regions, n_seas, n_ages, n_sexes)),
    do_recruits_move = 0,
    ln_InitDevs = array(0, dim = c(n_pop, n_regions, n_ages - 1))
  )

  expect_equal(out[1, 1, 1, 1], R0 * 0.6) # age 1, sex 1
  expect_equal(out[1, 1, 1, 2], R0 * 0.4) # age 1, sex 2
  # both sexes should decay at the same rate here since M is sex-invariant in this setup
  expect_equal(out[1, 1, 2, 1] / out[1, 1, 1, 1], out[1, 1, 2, 2] / out[1, 1, 1, 2], tolerance = 1e-8)
})

## ln_InitDevs application
test_that("ln_InitDevs multiplicatively scales ages 2:n_ages but leaves age 1 untouched", {

  d <- make_unfished_single(n_ages = 4, M = 0.2, R0 = 1000)

  out_baseline <- Get_Init_NAA(
    init_age_strc = 1,
    init_iter = 0,
    n_regions = d$n_regions,
    n_pop = d$n_pop,
    n_sexes = d$n_sexes,
    n_ages = d$n_ages,
    n_seas = d$n_seas,
    n_fish_fleets = d$n_fish_fleets,
    seasdur = d$seasdur,
    rec_seas_prop = d$rec_seas_prop,
    natmort = d$natmort,
    init_F = d$init_F,
    dmr = d$dmr,
    fish_sel = d$fish_sel,
    ret_sel = d$ret_sel,
    R0_r = d$R0_r,
    sexratio = d$sexratio,
    Movement = d$Movement,
    do_recruits_move = d$do_recruits_move,
    ln_InitDevs = d$ln_InitDevs
  )

  d$ln_InitDevs[1, 1, ] <- log(2) # double ages 2:4
  out_dev <- Get_Init_NAA(
    init_age_strc = 1,
    init_iter = 0,
    n_regions = d$n_regions,
    n_pop = d$n_pop,
    n_sexes = d$n_sexes,
    n_ages = d$n_ages,
    n_seas = d$n_seas,
    n_fish_fleets = d$n_fish_fleets,
    seasdur = d$seasdur,
    rec_seas_prop = d$rec_seas_prop,
    natmort = d$natmort,
    init_F = d$init_F,
    dmr = d$dmr,
    fish_sel = d$fish_sel,
    ret_sel = d$ret_sel,
    R0_r = d$R0_r,
    sexratio = d$sexratio,
    Movement = d$Movement,
    do_recruits_move = d$do_recruits_move,
    ln_InitDevs = d$ln_InitDevs
  )

  expect_equal(out_dev[1, 1, 1, 1], out_baseline[1, 1, 1, 1]) # age 1 untouched
  expect_equal(as.numeric(out_dev[1, 1, 2:4, 1]), as.numeric(out_baseline[1, 1, 2:4, 1]) * 2, tolerance = 1e-8)
})

## Multi-region movement (init_age_strc = 2): self-consistency of the
## plus-group linear solve, rather than an independently re-derived value.
test_that("init_age_strc = 2 plus-group solution satisfies its own equilibrium equation under 2-region movement", {

  n_pop <- 1; n_regions <- 2; n_sexes <- 1; n_ages <- 4; n_seas <- 1; n_fish_fleets <- 1
  M <- 0.2; R0 <- 500

  Movement <- array(0, dim = c(n_pop, n_regions, n_regions, n_seas, n_ages, n_sexes))
  Movement[, 1, 1, , , ] <- 0.7; Movement[, 1, 2, , , ] <- 0.3
  Movement[, 2, 1, , , ] <- 0.3; Movement[, 2, 2, , , ] <- 0.7

  natmort <- array(M, dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes))
  init_F <- array(0, dim = c(n_regions, n_seas, n_fish_fleets))
  dmr <- array(0, dim = c(n_regions, n_seas, n_fish_fleets))
  fish_sel <- array(0, dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes, n_fish_fleets))
  ret_sel <- array(1, dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes, n_fish_fleets))
  R0_r <- matrix(R0, nrow = n_pop, ncol = n_regions)

  out <- Get_Init_NAA(
    init_age_strc = 2,
    init_iter = 0,
    n_regions = n_regions,
    n_pop = n_pop,
    n_sexes = n_sexes,
    n_ages = n_ages,
    n_seas = n_seas,
    n_fish_fleets = n_fish_fleets,
    seasdur = 1,
    rec_seas_prop = matrix(1, nrow = n_pop, ncol = n_seas),
    natmort = natmort,
    init_F = init_F,
    dmr = dmr,
    fish_sel = fish_sel,
    ret_sel = ret_sel,
    R0_r = R0_r,
    sexratio = array(1, dim = c(n_pop, n_regions, n_sexes)),
    Movement = Movement,
    do_recruits_move = 0,
    ln_InitDevs = array(0, dim = c(n_pop, n_regions, n_ages - 1))
  )

  expect_equal(dim(out), c(n_pop, n_regions, n_ages, n_sexes))
  expect_true(all(is.finite(out)))
  expect_true(all(out > 0))

  # Self-consistency: N_plus should satisfy N_plus = T_plus %*% N_plus + T_penult %*% N_penult,
  # i.e. (I - T_plus) %*% N_plus == T_penult %*% N_penult, using the same annual
  # transition matrices the source builds (survival + movement, single season here).
  S <- diag(exp(-M), n_regions)
  T_mat <- t(Movement[1, , , 1, n_ages, 1]) %*% S # same transpose convention as source
  N_plus <- out[1, , n_ages, 1]
  N_penult <- out[1, , n_ages - 1, 1]
  lhs <- (diag(n_regions) - T_mat) %*% N_plus
  rhs <- T_mat %*% N_penult # T_penult == T_plus here since fish_sel = 0 and natmort is age-invariant
  expect_equal(as.numeric(lhs), as.numeric(rhs), tolerance = 1e-6)
})

## Cross-method consistency smoke test (n_regions = 1: all 4 methods agree)
test_that("all four init_age_strc methods agree with each other when n_regions = 1", {

  d <- make_unfished_single(n_ages = 5, M = 0.15, R0 = 800)

  run <- function(method, init_iter = 0) {
    Get_Init_NAA(
      init_age_strc = method,
      init_iter = init_iter,
      n_regions = d$n_regions,
      n_pop = d$n_pop,
      n_sexes = d$n_sexes,
      n_ages = d$n_ages,
      n_seas = d$n_seas,
      n_fish_fleets = d$n_fish_fleets,
      seasdur = d$seasdur,
      rec_seas_prop = d$rec_seas_prop,
      natmort = d$natmort,
      init_F = d$init_F,
      dmr = d$dmr,
      fish_sel = d$fish_sel,
      ret_sel = d$ret_sel,
      R0_r = d$R0_r,
      sexratio = d$sexratio,
      Movement = d$Movement,
      do_recruits_move = d$do_recruits_move,
      ln_InitDevs = d$ln_InitDevs
    )
  }

  out0 <- run(0, init_iter = 60)
  out1 <- run(1)
  out2 <- run(2)
  out3 <- run(3)

  expect_equal(as.numeric(out1), as.numeric(out2), tolerance = 1e-8)
  expect_equal(as.numeric(out1), as.numeric(out3), tolerance = 1e-8)
  expect_equal(as.numeric(out1), as.numeric(out0), tolerance = 1e-3)
})

## Basic dimension / sanity smoke test with multiple pops, regions, sexes, seasons
test_that("Get_Init_NAA() returns correctly-dimensioned, finite, non-negative output for a larger multi-dimensional setup", {

  n_pop <- 2; n_regions <- 2; n_sexes <- 2; n_ages <- 6; n_seas <- 2; n_fish_fleets <- 1
  M <- 0.25

  natmort <- array(M, dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes))
  init_F <- array(0.1, dim = c(n_regions, n_seas, n_fish_fleets))
  dmr <- array(0.2, dim = c(n_regions, n_seas, n_fish_fleets))
  fish_sel <- array(0.6, dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes, n_fish_fleets))
  ret_sel <- array(0.7, dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes, n_fish_fleets))
  R0_r <- matrix(c(300, 200, 250, 150), nrow = n_pop, ncol = n_regions)

  sexratio <- array(0.5, dim = c(n_pop, n_regions, n_sexes))

  Movement <- array(0, dim = c(n_pop, n_regions, n_regions, n_seas, n_ages, n_sexes))
  Movement[, 1, 1, , , ] <- 0.8; Movement[, 1, 2, , , ] <- 0.2
  Movement[, 2, 1, , , ] <- 0.2; Movement[, 2, 2, , , ] <- 0.8

  rec_seas_prop <- matrix(c(0.6, 0.4), nrow = n_pop, ncol = n_seas, byrow = TRUE)
  ln_InitDevs <- array(0, dim = c(n_pop, n_regions, n_ages - 1))

  for (method in 0:3) {
    out <- Get_Init_NAA(
      init_age_strc = method,
      init_iter = if (method == 0) 30 else 0,
      n_regions = n_regions,
      n_pop = n_pop,
      n_sexes = n_sexes,
      n_ages = n_ages,
      n_seas = n_seas,
      n_fish_fleets = n_fish_fleets,
      seasdur = c(0.5, 0.5),
      rec_seas_prop = rec_seas_prop,
      natmort = natmort,
      init_F = init_F,
      dmr = dmr,
      fish_sel = fish_sel,
      ret_sel = ret_sel,
      R0_r = R0_r,
      sexratio = sexratio,
      Movement = Movement,
      do_recruits_move = 1,
      ln_InitDevs = ln_InitDevs
    )

    expect_equal(dim(out), c(n_pop, n_regions, n_ages, n_sexes), info = paste("method", method))
    expect_true(all(is.finite(out)), info = paste("method", method))
    expect_true(all(out >= 0), info = paste("method", method))
  }
})
