library(SPoRC)
library(testthat)
library(Matrix)

# The implicit option swaps Matrix::expm for (I - A/n)^-n inside mat_exp, which every
# move_timing = 2 operator and Get_Movement route through. These tests pin the three
# things that make that swap safe: expm_nsub = 0 is bit-for-bit the old code path, the
# implicit operator keeps the structural properties the dynamics rely on (non-negative
# fractions, exact abundance accounting), and it converges to the exponential as the
# substep count grows.

# Column-convention generator net of mortality, as build_seas_operator forms it
make_A <- function(n, Z, seed) {
  set.seed(seed)
  D <- matrix(stats::runif(n * n, 0.05, 0.5), n, n)
  diag(D) <- 0
  diag(D) <- -colSums(D)
  list(Q_row = t(D), A = D - diag(Z, n))
}

test_that("expm_nsub = 0 is bit-for-bit the exact matrix exponential", {
  # The default must not perturb any existing model, so this is identity, not tolerance.
  for (n in c(1, 2, 4, 6)) {
    Z <- stats::runif(n, 0.05, 0.6)
    m <- make_A(n, Z, seed = 100 + n)
    ref <- as.matrix(Matrix::expm(methods::as(m$A, "sparseMatrix")))
    expect_identical(SPoRC:::mat_exp(m$A, 0), ref,
                     label = sprintf("mat_exp exact, n = %d", n))
  }
})

test_that("expm_nsub = 1 is exactly solve(I - A)", {
  # This is the plain backward Euler step the option was asked for; the substep
  # machinery must degenerate to it with no extra matrix products.
  n <- 5
  Z <- stats::runif(n, 0.05, 0.6)
  m <- make_A(n, Z, seed = 11)
  expect_equal(SPoRC:::mat_exp(m$A, 1), solve(diag(n) - m$A), tolerance = 0)
})

test_that("repeated squaring matches repeated multiplication", {
  # mat_exp raises the substep operator by squaring, which is why nsub is restricted to
  # powers of two. Check it against the sequential product it stands in for.
  n <- 4
  Z <- stats::runif(n, 0.05, 0.6)
  m <- make_A(n, Z, seed = 12)
  for (nsub in c(1, 2, 4, 8, 64, 512)) {
    B <- solve(diag(n) - m$A / nsub)
    naive <- diag(n)
    for (i in 1:nsub) naive <- naive %*% B
    expect_equal(SPoRC:::mat_exp(m$A, nsub), naive, tolerance = 1e-12,
                 label = sprintf("repeated squaring, nsub = %d", nsub))
  }
})

test_that("the implicit operator stays non-negative and converges to the exponential", {
  n <- 5
  Z <- stats::runif(n, 0.05, 0.6)
  m <- make_A(n, Z, seed = 13)
  ref <- as.matrix(Matrix::expm(methods::as(m$A, "sparseMatrix")))

  errs <- sapply(c(1, 4, 16, 64, 256), function(nsub) {
    E <- SPoRC:::mat_exp(m$A, nsub)
    expect_true(min(E) >= 0, label = sprintf("non-negative at nsub = %d", nsub))
    max(abs(E - ref))
  })
  # first order in 1/n: each 4x increase in substeps must cut the error, and the
  # ratio must approach 4 rather than stalling
  expect_true(all(diff(errs) < 0))
  expect_gt(errs[1] / errs[5], 50)
  expect_lt(errs[5], 1e-2)
})

test_that("implicit movement fractions are exactly column-stochastic without mortality", {
  # Get_Movement exponentiates a pure generator, so the approximation must not leak
  # abundance: (I - Q/n)^-1 is the resolvent, which is stochastic for any n.
  n <- 6
  m <- make_A(n, rep(0, n), seed = 14)
  for (nsub in c(0, 1, 4, 32)) {
    E <- SPoRC:::mat_exp(m$A, nsub)
    expect_equal(unname(colSums(E)), rep(1, n), tolerance = 1e-12,
                 label = sprintf("colSums at nsub = %d", nsub))
  }
})

test_that("the fused operator and integral conserve abundance exactly under the implicit scheme", {
  # get_population_projection takes survivors from T and catch from Integral. Under the
  # exponential these balance by construction; under backward Euler they still do,
  # because 1'(I - A/n) = 1' + z'/n. Without that, catch and numbers at age would drift
  # apart by the discretisation error rather than agreeing to machine precision.
  for (nsub in c(1, 4, 32)) {
    for (dur in c(1, 0.4)) {
      n <- 4
      Z <- stats::runif(n, 0.05, 0.6) * dur
      m <- make_A(n, rep(0, n), seed = 15)
      N <- stats::runif(n, 50, 500)

      fused <- SPoRC:::seas_operator_and_integral(Z, m$Q_row, dur, expm_nsub = nsub)
      survivors <- as.vector(t(N) %*% fused$T)
      removals <- Z * as.vector(fused$Integral %*% N)

      expect_equal(sum(survivors) + sum(removals), sum(N), tolerance = 1e-10,
                   label = sprintf("accounting, nsub = %d, dur = %g", nsub, dur))
    } # end dur loop
  } # end nsub loop
})

test_that("the transition helpers thread expm_nsub through to the same operator", {
  # advance_seas, spawn_state and catch_at_age must all use the operator the flag asks
  # for, not silently fall back to the exact exponential.
  n <- 4
  dur <- 0.7
  Z <- stats::runif(n, 0.05, 0.6) * dur
  m <- make_A(n, rep(0, n), seed = 16)
  N <- stats::runif(n, 50, 500)
  nsub <- 8

  Tm <- build_seas_operator(NULL, Z, m$Q_row, dur, move_timing = 2, expm_nsub = nsub)
  expect_false(isTRUE(all.equal(Tm, build_seas_operator(NULL, Z, m$Q_row, dur, move_timing = 2))))

  expect_equal(advance_seas(N, NULL, Z, m$Q_row, dur, move_timing = 2, expm_nsub = nsub),
               as.vector(t(N) %*% Tm), tolerance = 1e-12)

  fused <- SPoRC:::seas_operator_and_integral(Z, m$Q_row, dur, expm_nsub = nsub)
  expect_equal(Tm, fused$T, tolerance = 1e-12)
  expect_equal(integrate_seas_abundance(N, Z, m$Q_row, dur, expm_nsub = nsub),
               as.vector(fused$Integral %*% N), tolerance = 1e-12)

  F_landed <- stats::runif(n, 0.02, 0.2)
  expect_equal(catch_at_age(N, NULL, Z, m$Q_row, dur, F_landed, move_timing = 2, expm_nsub = nsub),
               F_landed * as.vector(fused$Integral %*% N), tolerance = 1e-12)

  # spawn_state and survey_state propagate part of the way through the season
  expect_equal(spawn_state(N, NULL, Z, m$Q_row, dur, t_spawn = 0.3, move_timing = 2, expm_nsub = nsub),
               as.vector(SPoRC:::mat_exp((t(m$Q_row) * dur - diag(Z, n)) * 0.3, nsub) %*% N),
               tolerance = 1e-12)
  expect_equal(survey_state(N, NULL, Z, m$Q_row, dur, t_srv = 0.5, move_timing = 2, expm_nsub = nsub),
               as.vector(SPoRC:::mat_exp((t(m$Q_row) * dur - diag(Z, n)) * 0.5, nsub) %*% N),
               tolerance = 1e-12)
})

test_that("expm_nsub leaves move_timing 0 and 1 untouched", {
  # Those timings never take an exponential, so the flag must be inert for them.
  n <- 4
  Z <- stats::runif(n, 0.05, 0.6)
  m <- make_A(n, rep(0, n), seed = 17)
  Move <- t(as.matrix(Matrix::expm(methods::as(t(m$Q_row), "sparseMatrix"))))
  N <- stats::runif(n, 50, 500)
  for (mt in c(0, 1)) {
    expect_identical(build_seas_operator(Move, Z, m$Q_row, 1, mt, expm_nsub = 0),
                     build_seas_operator(Move, Z, m$Q_row, 1, mt, expm_nsub = 32))
    expect_identical(advance_seas(N, Move, Z, m$Q_row, 1, mt, expm_nsub = 0),
                     advance_seas(N, Move, Z, m$Q_row, 1, mt, expm_nsub = 32))
  } # end mt loop
})

test_that("mat_exp is differentiable and its gradient converges to the exact one", {
  # The whole point of the option is a cheaper adjoint, so the tape has to build and
  # the derivative has to be the derivative of the thing it approximates.
  n <- 4
  adj <- matrix(0, n, n)
  for (i in 1:(n - 1)) { adj[i, i + 1] <- 1; adj[i + 1, i] <- 1 }
  Z <- c(0.3, 0.25, 0.4, 0.2)

  f <- function(nsub) function(p) {
    Q <- adj * exp(p$lq)
    diag(Q) <- -colSums(Q)
    sum(SPoRC:::mat_exp(Q - diag(Z), nsub) %*% (1:n))
  }
  at <- log(0.2)
  exact <- RTMB::MakeADFun(f(0), list(lq = at), silent = TRUE)
  impl <- RTMB::MakeADFun(f(2048), list(lq = at), silent = TRUE)

  expect_equal(impl$fn(at), exact$fn(at), tolerance = 1e-3)
  expect_equal(as.vector(impl$gr(at)), as.vector(exact$gr(at)), tolerance = 1e-3)
})

test_that("build_plus_group_T threads expm_nsub into the reference point operators", {
  # The per-recruit and MSY routines compose the plus group transition from
  # build_seas_operator, so the flag has to survive that composition. Reference points
  # built on the exact exponential while the fit used the implicit solve would be
  # internally inconsistent, which is the failure this guards.
  n <- 3
  n_seas <- 2
  seasdur <- c(0.4, 0.6)
  set.seed(21)
  D <- matrix(stats::runif(n * n, 0.05, 0.5), n, n)
  diag(D) <- 0
  diag(D) <- -colSums(D)
  Q_row <- t(D)

  Mrate <- array(rep(Q_row, n_seas), dim = c(n, n, n_seas))
  Mov <- array(rep(diag(n), n_seas), dim = c(n, n, n_seas))
  M_penult <- rep(0.2, n)   # natural mortality is a per-region vector
  M_plus <- rep(0.25, n)
  F_penult <- matrix(0.1, n, n_seas)
  F_plus <- matrix(0.12, n, n_seas)

  args <- list(M_penult = M_penult, M_plus = M_plus, F_penult = F_penult, F_plus = F_plus,
               Mov_penult = Mov, Mov_plus = Mov, n_regions = n, n_seas = n_seas,
               seasdur = seasdur, Mrate_penult = Mrate, Mrate_plus = Mrate,
               move_timing = 2)

  exact <- do.call(SPoRC:::build_plus_group_T, c(args, list(expm_nsub = 0)))
  impl1 <- do.call(SPoRC:::build_plus_group_T, c(args, list(expm_nsub = 1)))
  impl512 <- do.call(SPoRC:::build_plus_group_T, c(args, list(expm_nsub = 512)))

  # the flag must actually change the operator, then converge back onto the exact one
  expect_false(isTRUE(all.equal(exact, impl1)))
  for (nm in names(exact)) {
    expect_equal(impl512[[nm]], exact[[nm]], tolerance = 1e-3, label = nm)
    expect_gt(max(abs(impl1[[nm]] - exact[[nm]])), max(abs(impl512[[nm]] - exact[[nm]])))
  } # end nm loop
})

test_that("Setup_Mod_Movement rejects a substep count that is not a power of two", {
  # Substeps are applied by repeated squaring, so anything else would silently be rounded
  # to a different scheme than the user asked for. Reject it at setup instead.
  il <- Setup_Mod_Dim(years = 1:5, ages = 1:6, lens = NULL, n_regions = 3, n_sexes = 1,
                      n_fish_fleets = 1, n_srv_fleets = 1, n_seas = 1, n_pop = 1,
                      verbose = FALSE)
  adj <- matrix(1L, 3, 3); diag(adj) <- 0L
  mv <- function(nsub) Setup_Mod_Movement(
    input_list = il, move_type = 1, do_recruits_move = 0, use_fixed_movement = 0,
    ctmc_move_dat = expand.grid(pop = 1, regions = 1:3, years = 1:5, seas = 1,
                                ages = 1:6, sexes = 1),
    adjacency_mat = adj, area_r = rep(1, 3), diffusion_formula = ~ 1,
    preference_formula = ~ 0, move_timing = 2, move_expm_nsub = nsub, verbose = FALSE)

  for (good in c(0, 1, 2, 4, 8, 512)) {
    expect_equal(mv(good)$data$move_expm_nsub, as.integer(good))
  } # end good loop

  for (bad in c(3, 5, 6, 13, 100)) {
    expect_error(mv(bad), "power of two", info = paste("nsub =", bad))
  } # end bad loop

  # and the earlier guards still hold
  expect_error(mv(-1), "not correctly specified")
  expect_error(mv(2.5), "not correctly specified")
  expect_error(mv("8"), "not correctly specified")
})
