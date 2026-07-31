library(SPoRC)
library(testthat)

# Get_Fdev_PE_loglik keys on the deviation map mirrored into the data list rather
# than recomputing the fished set, so that a deviation mapped off by hand is
# neither estimated nor penalized. This builds the map do_Fmort_mapping would.
fdev_map <- function(UseCatch, UseCatch_pop, ObsCatch) {
  has_catch <- UseCatch == 1 | apply(UseCatch_pop == 1, c(2,3,4,5), any) | is.na(ObsCatch)
  map <- array(NA_real_, dim = dim(UseCatch))
  map[has_catch] <- seq_len(sum(has_catch))
  map
}

# ── do_Fdev_rho_mapping ──────────────────────────────────────────────────────

test_that("do_Fdev_rho_mapping only activates Fdev_rho under Fdev_model = 'ar1'", {

  make_il <- function(Fdev_model_code, n_regions = 2, n_seas = 2, n_fish_fleets = 2) {
    messages_list <<- character(0) # collect_message() writes to this global
    list(
      data = list(n_regions = n_regions, n_seas = n_seas, n_fish_fleets = n_fish_fleets,
                 Fdev_model = Fdev_model_code),
      par = list(Fdev_rho = array(0, dim = c(n_regions, n_seas, n_fish_fleets))),
      map = list()
    )
  }

  test_that("Fdev_model = iid (1) maps Fdev_rho entirely to NA regardless of spec", {
    il <- make_il(1)
    il <- SPoRC:::do_Fdev_rho_mapping(il, "est_all")
    expect_true(all(is.na(il$map$Fdev_rho)))
  })

  test_that("Fdev_model = rw (2) maps Fdev_rho entirely to NA regardless of spec", {
    il <- make_il(2)
    il <- SPoRC:::do_Fdev_rho_mapping(il, "est_all")
    expect_true(all(is.na(il$map$Fdev_rho)))
  })

  test_that("Fdev_model = ar1 (3) follows the requested sharing spec", {
    il <- make_il(3)
    il <- SPoRC:::do_Fdev_rho_mapping(il, "est_shared_r_seas_f")
    expect_equal(length(unique(il$map$Fdev_rho[!is.na(il$map$Fdev_rho)])), 1)

    il2 <- make_il(3)
    il2 <- SPoRC:::do_Fdev_rho_mapping(il2, "est_all")
    expect_equal(length(unique(il2$map$Fdev_rho[!is.na(il2$map$Fdev_rho)])), 2 * 2 * 2)
  })
})

# ── Get_Fdev_PE_loglik ───────────────────────────────────────────────────────

test_that("Get_Fdev_PE_loglik matches hand-computed values for iid/rw/ar1", {

  n_regions <- 2; n_yrs <- 5; n_seas <- 1; n_fish_fleets <- 2

  # region 1, fleet 1: catch starts late (years 2:5) -- "first active" is NOT
  # calendar year 1, exercising the general case.
  # region 2, fleet 1: catch active all years.
  # fleet 2: entirely inactive (no catch anywhere) -- should contribute 0.
  UseCatch <- array(0, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))
  UseCatch[1, 2:5, 1, 1] <- 1
  UseCatch[2, 1:5, 1, 1] <- 1
  UseCatch_pop <- array(0, dim = c(1, n_regions, n_yrs, n_seas, n_fish_fleets))
  ObsCatch <- array(0, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets)) # true recorded zeros, no missing data here

  set.seed(1)
  ln_F_devs <- array(rnorm(n_regions * n_yrs * n_seas * n_fish_fleets, sd = 0.4),
                     dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))
  ln_sigmaF <- array(log(seq(0.3, 0.6, length.out = n_regions * n_seas * n_fish_fleets)),
                     dim = c(n_regions, n_seas, n_fish_fleets))
  Fdev_rho <- array(qlogis((seq(0.2, 0.7, length.out = n_regions * n_seas * n_fish_fleets) + 1) / 2) / 2,
                    dim = c(n_regions, n_seas, n_fish_fleets)) # arbitrary unconstrained values

  rho_trans <- function(x) 2 / (1 + exp(-2 * x)) - 1

  hand_ll_array <- function(PE_model) {
    out <- array(0, dim = dim(ln_F_devs))
    for (f in 1:n_fish_fleets) {
      for (r in 1:n_regions) {
        for (seas in 1:n_seas) {
          sigma <- exp(ln_sigmaF[r, seas, f])
          rho <- rho_trans(Fdev_rho[r, seas, f])
          active <- which(UseCatch[r, , seas, f] == 1)
          if (length(active) == 0) next
          first <- TRUE
          for (y in active) {
            if (PE_model == 1) {
              out[r, y, seas, f] <- -dnorm(ln_F_devs[r, y, seas, f], 0, sigma, log = TRUE)
            } else if (PE_model == 2) {
              if (first) out[r, y, seas, f] <- -dnorm(ln_F_devs[r, y, seas, f], 0, 5, log = TRUE)
              else out[r, y, seas, f] <- -dnorm(ln_F_devs[r, y, seas, f], ln_F_devs[r, y - 1, seas, f], sigma, log = TRUE)
            } else if (PE_model == 3) {
              if (first) out[r, y, seas, f] <- -dnorm(ln_F_devs[r, y, seas, f], 0, sigma / sqrt(1 - rho^2), log = TRUE)
              else out[r, y, seas, f] <- -dnorm(ln_F_devs[r, y, seas, f], rho * ln_F_devs[r, y - 1, seas, f], sigma, log = TRUE)
            }
            first <- FALSE
          }
        }
      }
    }
    out
  }

  for (PE_model in 1:3) {
    got <- SPoRC:::Get_Fdev_PE_loglik(PE_model = PE_model, ln_sigmaF = ln_sigmaF, Fdev_rho = Fdev_rho,
                                      ln_F_devs = ln_F_devs, map_ln_F_devs = fdev_map(UseCatch, UseCatch_pop, ObsCatch))
    expect_equal(got, hand_ll_array(PE_model), tolerance = 1e-10, info = paste("PE_model:", PE_model))

    # entirely-inactive fleet 2 contributes nothing
    expect_true(all(got[, , , 2] == 0), info = paste("inactive fleet zero, PE_model:", PE_model))

    # skipped years (region 1, fleet 1, year 1) contribute nothing
    expect_equal(got[1, 1, 1, 1], 0, info = paste("skipped year zero, PE_model:", PE_model))
  }
})

test_that("Get_Fdev_PE_loglik handles multi-year gaps via the closed-form marginal transition", {

  # region 1: catch active in years 1, 2, 3, then closed for 3 years (4, 5, 6),
  # resuming in year 7 -- a gap of d = 4 between the last active year (3) and
  # the next (7). No sharing across regions/seasons/fleets needed here.
  n_regions <- 1; n_yrs <- 7; n_seas <- 1; n_fish_fleets <- 1
  UseCatch <- array(0, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))
  UseCatch[1, c(1,2,3,7), 1, 1] <- 1
  UseCatch_pop <- array(0, dim = c(1, n_regions, n_yrs, n_seas, n_fish_fleets))
  ObsCatch <- array(0, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets)) # true recorded zeros during the closure, not missing

  set.seed(3)
  ln_F_devs <- array(rnorm(n_yrs, sd = 0.4), dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))
  sigma <- 0.45
  ln_sigmaF <- array(log(sigma), dim = c(n_regions, n_seas, n_fish_fleets))
  rho_raw <- 0.3
  rho <- 2 / (1 + exp(-2 * rho_raw)) - 1
  Fdev_rho <- array(rho_raw, dim = c(n_regions, n_seas, n_fish_fleets))

  d <- 7 - 3 # elapsed gap between the last active year (3) and the next (7)

  # random walk: variance inflates linearly with the elapsed gap
  got_rw <- SPoRC:::Get_Fdev_PE_loglik(PE_model = 2, ln_sigmaF = ln_sigmaF, Fdev_rho = Fdev_rho,
                                       ln_F_devs = ln_F_devs, map_ln_F_devs = fdev_map(UseCatch, UseCatch_pop, ObsCatch))
  expected_rw_y7 <- -dnorm(ln_F_devs[1,7,1,1], ln_F_devs[1,3,1,1], sigma * sqrt(d), log = TRUE)
  expect_equal(got_rw[1,7,1,1], expected_rw_y7, tolerance = 1e-10)

  # ar1: mean decays by rho^d, variance is sigma^2 * sum_{i=0}^{d-1} rho^(2i)
  got_ar1 <- SPoRC:::Get_Fdev_PE_loglik(PE_model = 3, ln_sigmaF = ln_sigmaF, Fdev_rho = Fdev_rho,
                                       ln_F_devs = ln_F_devs, map_ln_F_devs = fdev_map(UseCatch, UseCatch_pop, ObsCatch))
  geom_sum <- sum(rho^(2 * (0:(d-1))))
  expected_ar1_y7 <- -dnorm(ln_F_devs[1,7,1,1], rho^d * ln_F_devs[1,3,1,1], sigma * sqrt(geom_sum), log = TRUE)
  expect_equal(got_ar1[1,7,1,1], expected_ar1_y7, tolerance = 1e-10)

  # closed-form geometric sum matches the closed-form (1 - rho^(2d)) / (1 - rho^2)
  expect_equal(geom_sum, (1 - rho^(2*d)) / (1 - rho^2), tolerance = 1e-10)

  # d = 1 (contiguous) reduces exactly to the standard single-step transition
  expected_rw_y2 <- -dnorm(ln_F_devs[1,2,1,1], ln_F_devs[1,1,1,1], sigma, log = TRUE)
  expect_equal(got_rw[1,2,1,1], expected_rw_y2, tolerance = 1e-10)
  expected_ar1_y2 <- -dnorm(ln_F_devs[1,2,1,1], rho * ln_F_devs[1,1,1,1], sigma, log = TRUE)
  expect_equal(got_ar1[1,2,1,1], expected_ar1_y2, tolerance = 1e-10)
})

test_that("Get_Fdev_PE_loglik iid model matches the pre-refactor inline dnorm formula", {
  n_regions <- 2; n_yrs <- 3; n_seas <- 1; n_fish_fleets <- 1
  UseCatch <- array(1, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))
  UseCatch_pop <- array(0, dim = c(1, n_regions, n_yrs, n_seas, n_fish_fleets))
  ObsCatch <- array(0, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))
  set.seed(2)
  ln_F_devs <- array(rnorm(n_regions * n_yrs), dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))
  ln_sigmaF <- array(log(0.5), dim = c(n_regions, n_seas, n_fish_fleets))
  Fdev_rho <- array(0, dim = c(n_regions, n_seas, n_fish_fleets))

  got <- SPoRC:::Get_Fdev_PE_loglik(PE_model = 1, ln_sigmaF = ln_sigmaF, Fdev_rho = Fdev_rho,
                                    ln_F_devs = ln_F_devs, map_ln_F_devs = fdev_map(UseCatch, UseCatch_pop, ObsCatch))

  expected <- array(0, dim = dim(ln_F_devs))
  for (f in 1:n_fish_fleets) for (y in 1:n_yrs) for (r in 1:n_regions) for (seas in 1:n_seas) {
    expected[r, y, seas, f] <- -dnorm(ln_F_devs[r, y, seas, f], 0, exp(ln_sigmaF[r, seas, f]), log = TRUE)
  }
  expect_equal(got, expected, tolerance = 1e-10)
})

# ── Missing (NA) vs. true-zero (0) ObsCatch ─────────────────────────────────

test_that("do_Fmort_mapping estimates a deviation for missing (NA) ObsCatch but not for a true recorded zero", {

  n_regions <- 1; n_yrs <- 4; n_seas <- 1; n_fish_fleets <- 1
  UseCatch <- array(0, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets)) # nothing fit anywhere
  UseCatch_pop <- array(0, dim = c(1, n_regions, n_yrs, n_seas, n_fish_fleets))

  # year 1: true recorded zero (closure) -> no deviation estimated
  # year 2: missing (NA) -> deviation IS estimated
  # years 3-4: also true recorded zeros
  ObsCatch <- array(0, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))
  ObsCatch[1, 2, 1, 1] <- NA

  il <- list(
    data = list(n_regions = n_regions, n_fish_fleets = n_fish_fleets,
               years = 1:n_yrs, n_seas = n_seas,
               UseCatch = UseCatch, UseCatch_pop = UseCatch_pop, ObsCatch = ObsCatch),
    par = list(ln_F_devs = array(0, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))),
    map = list()
  )
  il <- SPoRC:::do_Fmort_mapping(il)
  map_arr <- il$map$ln_F_devs

  expect_true(is.na(map_arr[1])) # year 1: true zero -> NA (not estimated)
  expect_false(is.na(map_arr[2])) # year 2: missing -> estimated
  expect_true(is.na(map_arr[3])) # year 3: true zero -> NA
  expect_true(is.na(map_arr[4])) # year 4: true zero -> NA
})

test_that("Get_Fdev_PE_loglik treats a missing (NA) year as an ordinary active year, not a gap", {

  n_regions <- 1; n_yrs <- 5; n_seas <- 1; n_fish_fleets <- 1

  # years 1, 2 fit normally; year 3 is a true closure (recorded zero, UseCatch = 0);
  # year 4 is missing (NA, UseCatch = 0) -- should behave like an ordinary active
  # year (d = 1 relative to year 2, NOT a gap of d = 2 skipping over it); year 5
  # fit normally again, with d = 1 relative to year 4.
  UseCatch <- array(0, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))
  UseCatch[1, c(1,2,5), 1, 1] <- 1
  UseCatch_pop <- array(0, dim = c(1, n_regions, n_yrs, n_seas, n_fish_fleets))
  ObsCatch <- array(0, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))
  ObsCatch[1, 3, 1, 1] <- 0   # true closure
  ObsCatch[1, 4, 1, 1] <- NA  # missing

  set.seed(4)
  ln_F_devs <- array(rnorm(n_yrs, sd = 0.4), dim = c(n_regions, n_yrs, n_seas, n_fish_fleets))
  sigma <- 0.4
  ln_sigmaF <- array(log(sigma), dim = c(n_regions, n_seas, n_fish_fleets))
  Fdev_rho <- array(0, dim = c(n_regions, n_seas, n_fish_fleets))

  # random walk
  got_rw <- SPoRC:::Get_Fdev_PE_loglik(PE_model = 2, ln_sigmaF = ln_sigmaF, Fdev_rho = Fdev_rho,
                                       ln_F_devs = ln_F_devs, map_ln_F_devs = fdev_map(UseCatch, UseCatch_pop, ObsCatch))

  # year 3 (true closure) contributes nothing and is skipped from the active sequence
  expect_equal(got_rw[1,3,1,1], 0)

  # year 4 (missing) is estimated as an ordinary single-step transition from year 2
  # (its true previous active/missing neighbor), i.e. d = 4 - 2 = 2 since year 3 was
  # skipped as a closure -- NOT treated as itself starting a new gap
  expected_y4 <- -dnorm(ln_F_devs[1,4,1,1], ln_F_devs[1,2,1,1], sigma * sqrt(2), log = TRUE)
  expect_equal(got_rw[1,4,1,1], expected_y4, tolerance = 1e-10)

  # year 5 continues from year 4 (the missing year) at d = 1, NOT from year 2 at d = 3,
  # confirming the missing year itself becomes the new "last active" reference
  expected_y5 <- -dnorm(ln_F_devs[1,5,1,1], ln_F_devs[1,4,1,1], sigma, log = TRUE)
  expect_equal(got_rw[1,5,1,1], expected_y5, tolerance = 1e-10)
})
