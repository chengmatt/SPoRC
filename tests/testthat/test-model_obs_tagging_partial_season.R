test_that("partial-season tag exposure keeps F/Z a valid death fraction", {

  # A tag cohort released partway through a season is at liberty for a fraction t of it, so it
  # accrues F*t of fishing mortality and Z*t of total mortality. Baranov's ratio must therefore
  # stay F/Z. Scaling Z alone -- the pre-fix behavior -- gives F/(Z*t), which is not a fraction.
  F_seas <- 0.40
  M_seas <- 0.30
  Z_seas <- F_seas + M_seas

  for (t in c(0.2222, 0.3333, 0.5, 0.7778, 1)) {
    ratio_scaled_both <- (F_seas * t) / (Z_seas * t)
    ratio_scaled_z    <- F_seas / (Z_seas * t)

    expect_lte(ratio_scaled_both, 1)
    expect_equal(ratio_scaled_both, F_seas / Z_seas)

    # the pre-fix ratio is inflated by exactly 1/t, and leaves [0,1] once t < F/Z
    expect_equal(ratio_scaled_z, (F_seas / Z_seas) / t)
    if (t < F_seas / Z_seas) expect_gt(ratio_scaled_z, 1)
  }
})

test_that("catch_at_age reproduces the analytic partial-interval Baranov under both timings", {

  # Single region, so movement drops out and the answer is closed form. Both the discrete branch
  # and the move_timing = 2 (spatial Baranov) branch must agree with it once F and Z are scaled
  # together by the time at liberty.
  N <- 1000
  F_seas <- 0.40
  M_seas <- 0.30
  Z_seas <- F_seas + M_seas
  seasdur <- 0.75

  for (t in c(0.2222, 0.5, 1)) {
    F_eff <- F_seas * t
    Z_eff <- Z_seas * t
    analytic <- (F_eff / Z_eff) * (1 - exp(-Z_eff)) * N

    disc <- catch_at_age(
      N = N,
      Move = diag(1),
      Z = Z_eff,
      Q = NULL,
      dur = seasdur * t,
      F_landed = F_eff,
      move_timing = 0
    )

    cont <- catch_at_age(
      N = N,
      Move = diag(1),
      Z = Z_eff,
      Q = matrix(0, 1, 1),
      dur = seasdur * t,
      F_landed = F_eff,
      move_timing = 2
    )

    expect_equal(as.vector(disc), analytic, tolerance = 1e-10)
    expect_equal(as.vector(cont), analytic, tolerance = 1e-8)
  }
})

# Mid-season releases under continuous movement -------------------------------
# get_tagging_observation_model scales a mid-season cohort's mortality by tag_frac and
# its generator by tag_dur = seasdur * tag_frac. Under move_timing = 2 the cohort must
# therefore still move, for exactly the fraction of the season it was at liberty for --
# unlike timings 0 and 1, which skip the movement step because a full-season transition
# matrix cannot represent a partial interval.

tag_gen <- function(n, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  D <- matrix(stats::runif(n * n, 0.05, 0.5), n, n)
  diag(D) <- 0
  diag(D) <- -colSums(D)
  t(D) # row convention, as stored in Mrate
}

test_that("a mid-season release still redistributes under move_timing = 2", {

  n <- 3
  Q <- tag_gen(n, seed = 4242)
  seasdur <- 0.5
  Z_full <- c(0.20, 0.45, 0.70)
  N <- c(1000, 0, 0) # every tag released into region 1

  for (tag_frac in c(0.25, 0.5, 0.9)) {
    moved <- advance_seas(N, NULL, Z_full * tag_frac, Q, seasdur * tag_frac, move_timing = 2)
    frozen <- advance_seas(N, diag(n), Z_full * tag_frac, matrix(0, n, n),
                           seasdur * tag_frac, move_timing = 2)

    # tags reach the other regions rather than sitting in the release region
    expect_gt(moved[2], 0)
    expect_gt(moved[3], 0)
    expect_equal(as.vector(frozen[2:3]), c(0, 0), tolerance = 1e-12)
    expect_lt(moved[1], frozen[1])
  }
})

test_that("the partial interval is exactly a fraction of the full-season generator", {

  # tag_dur and the tag_frac-scaled Z together must reproduce exp(Lambda * tag_frac),
  # where Lambda is the ordinary full-season generator. If either were scaled without
  # the other, this identity would fail.
  n <- 3
  Q <- tag_gen(n, seed = 99)
  seasdur <- 0.75
  Z_full <- c(0.15, 0.40, 0.62)
  N <- c(600, 250, 90)

  Lambda <- t(Q) * seasdur - diag(Z_full, n) # column convention, full season

  for (tag_frac in c(0.1, 0.3333, 0.75, 1)) {
    got <- advance_seas(N, NULL, Z_full * tag_frac, Q, seasdur * tag_frac, move_timing = 2)
    want <- as.vector(as.matrix(Matrix::expm(methods::as(Lambda * tag_frac, "sparseMatrix"))) %*% N)
    expect_equal(got, want, tolerance = 1e-10,
                 label = sprintf("partial propagation, tag_frac = %g", tag_frac))
  }
})

test_that("two partial tag steps compose into one full-season step", {

  # Semigroup check: being at liberty for two halves of a season must land in the same
  # place as being at liberty for the whole of it. This only holds if the generator is
  # scaled by the at-liberty fraction rather than frozen.
  n <- 4
  Q <- tag_gen(n, seed = 777)
  seasdur <- 0.6
  Z_full <- c(0.10, 0.35, 0.55, 0.22)
  N <- c(400, 300, 200, 100)

  half <- advance_seas(N, NULL, Z_full * 0.5, Q, seasdur * 0.5, move_timing = 2)
  twice <- advance_seas(half, NULL, Z_full * 0.5, Q, seasdur * 0.5, move_timing = 2)
  once <- advance_seas(N, NULL, Z_full, Q, seasdur, move_timing = 2)

  expect_equal(twice, once, tolerance = 1e-10)
})

test_that("partial-interval tag recaptures match direct integration of the true dynamics", {

  # The recapture block computes F_full * tag_frac times the unit-interval integral over
  # the tag_frac-scaled generator. That should equal the exact partial-interval catch,
  # integral over [0, tag_frac] of F_full * exp(Lambda * u) * N du. Verified against
  # quadrature on the untouched full-season generator.
  n <- 3
  Q <- tag_gen(n, seed = 31415)
  seasdur <- 0.8
  Z_full <- c(0.18, 0.44, 0.66)
  F_full <- c(0.09, 0.21, 0.30)
  N <- c(800, 150, 50)

  Lambda <- t(Q) * seasdur - diag(Z_full, n)

  for (tag_frac in c(0.2, 0.5, 1)) {
    got <- (F_full * tag_frac) *
      integrate_seas_abundance(N, Z_full * tag_frac, Q, seasdur * tag_frac)

    grid <- seq(0, tag_frac, length.out = 8001)
    step <- grid[2] - grid[1]
    vals <- lapply(grid, function(u) {
      as.vector(as.matrix(Matrix::expm(methods::as(Lambda * u, "sparseMatrix"))) %*% N)
    })
    quad <- (Reduce(`+`, vals) - (vals[[1]] + vals[[length(vals)]]) / 2) * step
    want <- F_full * quad

    expect_equal(got, want, tolerance = 1e-6,
                 label = sprintf("partial recapture integral, tag_frac = %g", tag_frac))
  }
})

test_that("scaling every mortality component is equivalent to scaling the total", {

  # The patch scales tmp_ZAA, tmp_FAA, tmp_ret_FAA and tmp_disc_DAA by tag_frac. Scaling the
  # components and summing must equal scaling the summed total, or the tag Z would silently
  # change meaning relative to the population Z.
  t <- 0.3333
  natmort <- 0.30; seasdur <- 0.75; shed <- 0.05
  ret_F <- 0.25; disc_F <- 0.10

  Z_then_scale <- ((natmort * seasdur) + (ret_F + disc_F) + (shed * seasdur)) * t
  Z_from_parts <- ((natmort * seasdur) * t) + ((ret_F * t) + (disc_F * t)) + ((shed * seasdur) * t)

  expect_equal(Z_then_scale, Z_from_parts)
})
