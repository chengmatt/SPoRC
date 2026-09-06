library(SPoRC)
library(testthat)

test_that("the fused operator/integral agrees with computing each separately", {
  # get_population_projection takes both the transition operator and the catch integral
  # from one Van Loan exponential. That is only sound because the operator really is the
  # top-left block of the same matrix, so pin it: any change here would silently alter
  # either the numbers at age or the catch.
  set.seed(808)
  for (n in c(1, 2, 4)) {
    for (dur in c(1, 0.25, 0.6)) {
      Q <- matrix(stats::runif(n * n, 0.05, 0.5), n, n)
      diag(Q) <- 0; diag(Q) <- -rowSums(Q)
      Z <- stats::runif(n, 0.02, 0.6) * dur
      N <- stats::runif(n, 50, 500)

      fused <- SPoRC:::seas_operator_and_integral(Z, Q, dur)

      expect_equal(fused$T, build_seas_operator(NULL, Z, Q, dur, move_timing = 2),
                   tolerance = 1e-12,
                   label = sprintf("fused operator, n = %d, dur = %g", n, dur))
      expect_equal(as.vector(t(N) %*% fused$T),
                   advance_seas(N, NULL, Z, Q, dur, move_timing = 2),
                   tolerance = 1e-12)
      expect_equal(as.vector(fused$Integral %*% N),
                   integrate_seas_abundance(N, Z, Q, dur),
                   tolerance = 1e-12,
                   label = sprintf("fused integral, n = %d, dur = %g", n, dur))
    }
  }
})

test_that("the fused integral reduces to the ordinary Baranov without movement", {
  # With a zero generator the regions decouple, so the integral must collapse to the
  # scalar (1 - exp(-Z)) / Z that the timing 0 and 1 catch equations use.
  n <- 3
  Z <- c(0.35, 0.12, 0.5)
  N <- c(400, 250, 90)
  fused <- SPoRC:::seas_operator_and_integral(Z, matrix(0, n, n), dur = 1)
  expect_equal(as.vector(fused$Integral %*% N), N * (1 - exp(-Z)) / Z, tolerance = 1e-12)
  expect_equal(as.vector(t(N) %*% fused$T), N * exp(-Z), tolerance = 1e-12)
})
library(Matrix)

# Helpers ---------------------------------------------------------------------

# Build a valid CTMC generator in COLUMN convention (colSums 0), matching how
# Get_Movement constructs Q_ss, and return both the stored row-convention
# generator (as Mrate holds it) and the corresponding movement fractions
# (as Movement holds them).
make_move <- function(n, dur = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  D <- matrix(stats::runif(n * n, 0.05, 0.5), n, n)
  diag(D) <- 0
  diag(D) <- -colSums(D)
  list(Q_row = t(D), Move_row = t(as.matrix(Matrix::expm(D * dur))))
}

test_that("movement matrices built for testing are row-stochastic", {
  m <- make_move(4, seed = 1)
  expect_equal(unname(rowSums(m$Move_row)), rep(1, 4), tolerance = 1e-12)
  expect_equal(unname(rowSums(m$Q_row)), rep(0, 4), tolerance = 1e-12)
})

# move_timing = 0 must reproduce the historical expressions exactly ------------

test_that("move_timing = 0 reproduces the legacy move-then-die expression exactly", {
  m <- make_move(4, seed = 2)
  N <- c(100, 250, 30, 80)
  Z <- c(0.2, 0.5, 0.1, 0.9)

  legacy <- as.vector(t(N) %*% m$Move_row) * exp(-Z)
  expect_identical(
    advance_seas(N, m$Move_row, Z, m$Q_row, dur = 1, move_timing = 0),
    legacy
  )
})

test_that("build_seas_operator at move_timing = 0 equals Move %*% diag(surv)", {
  m <- make_move(3, seed = 3)
  Z <- c(0.3, 0.15, 0.6)
  expect_equal(
    build_seas_operator(m$Move_row, Z, m$Q_row, dur = 1, move_timing = 0),
    m$Move_row %*% diag(exp(-Z), 3)
  )
})

test_that("move_timing = 1 equals the die-then-move expression exactly", {
  m <- make_move(4, seed = 4)
  N <- c(100, 250, 30, 80)
  Z <- c(0.2, 0.5, 0.1, 0.9)
  expect_identical(
    advance_seas(N, m$Move_row, Z, m$Q_row, dur = 1, move_timing = 1),
    as.vector(t(N * exp(-Z)) %*% m$Move_row)
  )
})

# Cases where all three timings must coincide ---------------------------------

test_that("all three timings agree when Z is constant across regions", {
  # diag(z * I) commutes with the generator, so ordering cannot matter
  m <- make_move(4, seed = 5)
  N <- c(100, 250, 30, 80)
  Z <- rep(0.37, 4)

  t0 <- advance_seas(N, m$Move_row, Z, m$Q_row, 1, 0)
  t1 <- advance_seas(N, m$Move_row, Z, m$Q_row, 1, 1)
  t2 <- advance_seas(N, m$Move_row, Z, m$Q_row, 1, 2)

  expect_equal(t1, t0, tolerance = 1e-10)
  expect_equal(t2, t0, tolerance = 1e-10)
})

test_that("all three timings agree when there is no mortality", {
  m <- make_move(4, seed = 6)
  N <- c(100, 250, 30, 80)
  Z <- rep(0, 4)

  t0 <- advance_seas(N, m$Move_row, Z, m$Q_row, 1, 0)
  expect_equal(advance_seas(N, m$Move_row, Z, m$Q_row, 1, 1), t0, tolerance = 1e-12)
  expect_equal(advance_seas(N, m$Move_row, Z, m$Q_row, 1, 2), t0, tolerance = 1e-12)
})

test_that("all three timings agree when there is no movement", {
  n <- 4
  N <- c(100, 250, 30, 80)
  Z <- c(0.2, 0.5, 0.1, 0.9)
  I <- diag(n); Q0 <- matrix(0, n, n)

  expected <- N * exp(-Z)
  expect_equal(advance_seas(N, I, Z, Q0, 1, 0), expected, tolerance = 1e-12)
  expect_equal(advance_seas(N, I, Z, Q0, 1, 1), expected, tolerance = 1e-12)
  expect_equal(advance_seas(N, I, Z, Q0, 1, 2), expected, tolerance = 1e-12)
})

test_that("continuous movement conserves abundance when mortality is zero", {
  m <- make_move(4, seed = 7)
  N <- c(100, 250, 30, 80)
  out <- advance_seas(N, m$Move_row, rep(0, 4), m$Q_row, 1, 2)
  expect_equal(sum(out), sum(N), tolerance = 1e-10)
  expect_true(all(out >= 0))
})

# The timings differ when they should -------------------------------

test_that("timings differ under region-varying mortality, with continuous in between", {
  m <- make_move(4, seed = 8)
  N <- c(100, 250, 30, 80)
  Z <- c(0.2, 0.5, 0.1, 0.9)

  d0 <- advance_seas(N, m$Move_row, Z, m$Q_row, 1, 0)
  d1 <- advance_seas(N, m$Move_row, Z, m$Q_row, 1, 1)
  d2 <- advance_seas(N, m$Move_row, Z, m$Q_row, 1, 2)

  expect_false(isTRUE(all.equal(d0, d1)))
  expect_false(isTRUE(all.equal(d0, d2)))
  # continuous movement should sit between the two discrete orderings
  expect_true(all(d2 >= pmin(d0, d1) - 1e-9 & d2 <= pmax(d0, d1) + 1e-9))
})

# Spatial Baranov / Van Loan integral -----------------------------------------

test_that("season-integrated abundance reduces to the Baranov form without movement", {
  n <- 4
  N <- c(100, 250, 30, 80)
  Z <- c(0.2, 0.5, 0.1, 0.9)
  expect_equal(
    integrate_seas_abundance(N, Z, matrix(0, n, n), dur = 1),
    N * (1 - exp(-Z)) / Z,
    tolerance = 1e-10
  )
})

test_that("season-integrated abundance matches numerical integration under movement", {
  m <- make_move(3, seed = 9)
  N <- c(100, 250, 30)
  Z <- c(0.2, 0.5, 0.1)

  A <- t(m$Q_row) - diag(Z, 3)
  grid <- seq(0, 1, length.out = 20001)[-1]
  brute <- rowSums(vapply(grid, function(tau) {
    as.vector(as.matrix(Matrix::expm(A * tau)) %*% N)
  }, numeric(3))) / length(grid)

  vl <- integrate_seas_abundance(N, Z, m$Q_row, dur = 1)
  expect_equal(vl, brute, tolerance = 1e-4)
})

# Spawning-state propagator ---------------------------------------------------

test_that("spawn_state at move_timing = 0 reproduces the legacy expression", {
  m <- make_move(3, seed = 10)
  N <- c(100, 250, 30)
  Z <- c(0.2, 0.5, 0.1)
  t_spawn <- 0.35
  expect_identical(
    spawn_state(N, m$Move_row, Z, m$Q_row, 1, t_spawn, move_timing = 0),
    as.vector(t(N) %*% m$Move_row) * exp(-t_spawn * Z)
  )
})

test_that("spawn_state with t_spawn = 0 leaves continuous state unchanged", {
  m <- make_move(3, seed = 11)
  N <- c(100, 250, 30)
  Z <- c(0.2, 0.5, 0.1)
  expect_equal(spawn_state(N, m$Move_row, Z, m$Q_row, 1, 0, move_timing = 2), N,
               tolerance = 1e-12)
})

test_that("spawn_state with t_spawn = 1 equals a full continuous season step", {
  m <- make_move(3, seed = 12)
  N <- c(100, 250, 30)
  Z <- c(0.2, 0.5, 0.1)
  expect_equal(
    spawn_state(N, m$Move_row, Z, m$Q_row, 1, 1, move_timing = 2),
    advance_seas(N, m$Move_row, Z, m$Q_row, 1, 2),
    tolerance = 1e-10
  )
})

# Season-duration scaling of the generator ------------------------------------

test_that("continuous steps compose across seasons to the full-year operator", {
  # Four quarter-length seasons under a constant generator and constant mortality
  # must equal one full-year step, which is the property that makes the seasdur
  # scaling of Q meaningful.
  n <- 3
  m <- make_move(n, dur = 1, seed = 13)
  N <- c(100, 250, 30)
  Z_annual <- c(0.2, 0.5, 0.1)
  seasdur <- rep(0.25, 4)

  stepwise <- N
  for (s in 1:4) {
    stepwise <- advance_seas(stepwise, m$Move_row, Z_annual * seasdur[s],
                             m$Q_row, dur = seasdur[s], move_timing = 2)
  }
  oneshot <- advance_seas(N, m$Move_row, Z_annual, m$Q_row, dur = 1, move_timing = 2)
  expect_equal(stepwise, oneshot, tolerance = 1e-10)
})

test_that("Get_Movement scales the generator by seasdur when asked", {
  n_regions <- 3; n_seas <- 2
  dat <- expand.grid(
    pop = 1,
    regions = 1:n_regions,
    years = 1,
    seas = 1:n_seas,
    ages = 1:2,
    sexes = 1
  )
  adj <- matrix(1L, n_regions, n_regions); diag(adj) <- 0L
  seasdur <- c(0.25, 0.75)

  call_gm <- function(scale) {
    Get_Movement(
      move_type = 1,
      do_recruits_move = 1,
      n_pop = 1,
      n_regions = n_regions,
      n_yrs = 1,
      n_proj_yrs_devs = 0,
      n_ages = 2,
      n_sexes = 1,
      n_seas = n_seas,
      move_pars = array(0, c(1, n_regions, n_regions - 1, 1, n_seas, 2, 1)),
      move_devs = array(0, c(1, n_regions, n_regions - 1, 1, n_seas, 2, 1)),
      use_fixed_movement = 0,
      Fixed_Movement = NULL,
      ctmc_move_dat = dat,
      preference_formula = ~ 1,
      diffusion_formula = ~ 1,
      log_move_diffusion_pars = log(0.3),
      move_preference_pars = 0,
      area_r = rep(1, n_regions),
      adjacency_mat = adj,
      ctmc_diffusion_bounds = 0,
      seasdur = seasdur,
      ctmc_scale_by_seasdur = scale
    )
  }

  legacy <- call_gm(0)
  scaled <- call_gm(1)

  # Mrate is stored unscaled under both settings
  expect_equal(scaled$Mrate, legacy$Mrate, tolerance = 1e-12)

  # Legacy exponentiates Q once per season regardless of duration, so both seasons match
  expect_equal(legacy$Movement[1, , , 1, 1, 2, 1],
               legacy$Movement[1, , , 1, 2, 2, 1], tolerance = 1e-12)

  # Scaled movement must differ between unequal-duration seasons, and match expm(Q * dur)
  expect_false(isTRUE(all.equal(scaled$Movement[1, , , 1, 1, 2, 1],
                                scaled$Movement[1, , , 1, 2, 2, 1])))
  for (s in 1:n_seas) {
    Q_row <- scaled$Mrate[1, , , 1, s, 2, 1]
    expected <- t(as.matrix(Matrix::expm(t(Q_row) * seasdur[s])))
    expect_equal(scaled$Movement[1, , , 1, s, 2, 1], expected,
                 tolerance = 1e-10, ignore_attr = TRUE)
  }

  # Movement stays row-stochastic under scaling
  for (s in 1:n_seas) {
    expect_equal(rowSums(scaled$Movement[1, , , 1, s, 2, 1]), rep(1, n_regions),
                 tolerance = 1e-10, ignore_attr = TRUE)
  }
})

test_that("Get_Movement defaults to legacy unscaled generator behavior", {
  # Default ctmc_scale_by_seasdur = 0 keeps existing model results reproducible
  expect_equal(formals(Get_Movement)$ctmc_scale_by_seasdur, 0)
  expect_equal(formals(build_seas_operator)$move_timing, 0)
})

# Plus-group transition composition -------------------------------------------

test_that("composed annual transition equals stepping season by season", {
  # The plus-group recursions in Get_Init_NAA (init_age_strc = 2), Get_Det_Recruitment
  # and build_plus_group_T all build an annual operator by composing per-season ones.
  # That operator must reproduce what sequential seasonal stepping gives, for every
  # timing. init_age_strc = 2 previously composed t(Movement) %*% S right-to-left,
  # which both inverted the movement/mortality order and traversed seasons backwards.
  set.seed(11); n <- 3; n_seas <- 3
  seasdur <- c(0.2, 0.5, 0.3)
  M <- array(0, c(n, n, n_seas)); Q <- array(0, c(n, n, n_seas))
  for (s in seq_len(n_seas)) {
    mv <- make_move(n, dur = seasdur[s])
    M[, , s] <- mv$Move_row; Q[, , s] <- mv$Q_row
  }
  Zs <- matrix(stats::runif(n * n_seas, 0.05, 0.3), n, n_seas)
  N0 <- c(100, 40, 70)

  for (tm in 0:2) {
    stepwise <- N0
    for (s in seq_len(n_seas)) stepwise <- advance_seas(stepwise, M[, , s], Zs[, s], Q[, , s], seasdur[s], tm)

    Tann <- diag(n)
    for (s in seq_len(n_seas)) Tann <- t(build_seas_operator(M[, , s], Zs[, s], Q[, , s], seasdur[s], tm)) %*% Tann

    expect_equal(as.vector(Tann %*% N0), stepwise, tolerance = 1e-10,
                 label = sprintf("composed annual operator, move_timing = %d", tm))
  }
})

test_that("right-composing seasons is not equivalent to left-composing", {
  # Guards the season-ordering half of the init_age_strc = 2 fix: if these two agreed,
  # the test above could not distinguish the traversal order.
  set.seed(12); n <- 3
  m1 <- make_move(n); m2 <- make_move(n)
  S1 <- diag(exp(-c(0.2, 0.5, 0.1)), n); S2 <- diag(exp(-c(0.6, 0.1, 0.4)), n)
  left  <- (S2 %*% t(m2$Move_row)) %*% (S1 %*% t(m1$Move_row))
  right <- (S1 %*% t(m1$Move_row)) %*% (S2 %*% t(m2$Move_row))
  expect_false(isTRUE(all.equal(left, right)))
})

# Guard rails -----------------------------------------------------------------

test_that("build_seas_operator rejects an unknown move_timing", {
  m <- make_move(3, seed = 14)
  expect_error(build_seas_operator(m$Move_row, rep(0.1, 3), m$Q_row, 1, 3),
               "move_timing must be")
})

test_that("continuous movement requires a rate matrix", {
  m <- make_move(3, seed = 15)
  expect_error(build_seas_operator(m$Move_row, rep(0.1, 3), NULL, 1, 2),
               "required when move_timing = 2")
})

test_that("spawn_state rejects an unknown move_timing only via its callers", {
  # spawn_state falls through to the continuous branch for any timing > 1, so a
  # missing generator is the failure mode rather than a timing check
  m <- make_move(3, seed = 16)
  expect_error(spawn_state(c(1, 1, 1), m$Move_row, rep(0.1, 3), NULL, 1, 0.5, 2))
})

test_that("Setup_Mod_Movement validates move_timing", {
  # A bad value must be rejected before the not-yet-wired guard fires
  dummy <- list(
    data = list(
      n_pop = 1,
      n_regions = 2,
      n_seas = 1,
      years = 1:3,
      ages = 1:4,
      n_sexes = 1,
      n_proj_yrs_devs = 0
    ),
    par = list(),
    map = list(),
    verbose = FALSE,
    store_config = FALSE
  )
  expect_error(Setup_Mod_Movement(dummy, move_timing = 7), "move_timing is not correctly specified")
})

test_that("continuous movement requires an estimated CTMC generator", {
  # Unstructured multinomial-logit movement has no instantaneous rate matrix, and one
  # cannot in general be recovered from the fractions (the Markov embedding problem)
  dummy <- list(
    data = list(
      n_pop = 1,
      n_regions = 2,
      n_seas = 1,
      years = 1:3,
      ages = 1:4,
      n_sexes = 1,
      n_proj_yrs_devs = 0
    ),
    par = list(),
    map = list(),
    verbose = FALSE,
    store_config = FALSE
  )
  expect_error(Setup_Mod_Movement(dummy, move_timing = 2, move_type = 0),
               "requires move_type == 1")
})

# Survey timing --------------------------------------------------------------

test_that("survey_state leaves the discrete timings unchanged", {
  # Under move_timing 0 and 1 movement is already resolved for the season, so the
  # historical elementwise discount is exact and must not be perturbed
  m <- make_move(3, seed = 20)
  N <- c(100, 60, 30); Z <- c(0.3, 0.15, 0.45)
  for (tm in 0:1) {
    expect_equal(survey_state(N, m$Move_row, Z, m$Q_row, 1, 0.4, tm), N * exp(-0.4 * Z))
  }
})

test_that("survey_state brackets correctly under continuous movement", {
  m <- make_move(3, seed = 21)
  N <- c(100, 60, 30); Z <- c(0.3, 0.15, 0.45)
  # t_srv = 0 is the start of the season, t_srv = 1 a full seasonal step
  expect_equal(survey_state(N, m$Move_row, Z, m$Q_row, 1, 0, 2), N, tolerance = 1e-12)
  expect_equal(survey_state(N, m$Move_row, Z, m$Q_row, 1, 1, 2),
               advance_seas(N, m$Move_row, Z, m$Q_row, 1, 2), tolerance = 1e-10)
  # with no movement it must collapse onto the elementwise discount
  expect_equal(survey_state(N, m$Move_row, Z, matrix(0, 3, 3), 1, 0.4, 2),
               N * exp(-0.4 * Z), tolerance = 1e-10)
  # and with movement it must actually differ from it
  expect_false(isTRUE(all.equal(survey_state(N, m$Move_row, Z, m$Q_row, 1, 0.4, 2),
                                N * exp(-0.4 * Z))))
})

test_that("region-varying survey timing reads each region's own propagation", {
  # t_srv can differ by region, which no single operator represents; the convention is
  # that region r observes the population propagated to that region's survey time
  m <- make_move(3, seed = 22)
  N <- c(100, 60, 30); Z <- c(0.3, 0.15, 0.45)
  tv <- c(0.2, 0.5, 0.8)
  mixed <- survey_state(N, m$Move_row, Z, m$Q_row, 1, tv, 2)
  ref <- vapply(seq_along(tv), function(r) survey_state(N, m$Move_row, Z, m$Q_row, 1, tv[r], 2)[r],
                numeric(1))
  expect_equal(mixed, ref, tolerance = 1e-12)
  # constant t_srv must take the single-propagation branch and agree with it
  expect_equal(survey_state(N, m$Move_row, Z, m$Q_row, 1, rep(0.4, 3), 2),
               survey_state(N, m$Move_row, Z, m$Q_row, 1, 0.4, 2), tolerance = 1e-12)
})
