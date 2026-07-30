library(SPoRC)
library(testthat)

# Minimal input_list carrying only what do_Fmort_mapping reads
mk_input <- function(UseCatch, ObsCatch, UseCatch_pop = NULL) {
  d <- dim(UseCatch)
  if (is.null(UseCatch_pop)) UseCatch_pop <- array(0, dim = c(1, d))
  list(data = list(n_regions = d[1], years = seq_len(d[2]), n_seas = d[3],
                   n_fish_fleets = d[4], UseCatch = UseCatch,
                   UseCatch_pop = UseCatch_pop, ObsCatch = ObsCatch),
       par = list(), map = list())
}

test_that("ln_F_mean is estimated only where a cell is fished in some year", {

  # 3 regions x 2 seasons -> 6 fleets, each fishing only its own region/season
  n_reg <- 3; n_seas <- 2; n_yrs <- 10; n_ff <- n_reg * n_seas
  UseCatch <- array(0, dim = c(n_reg, n_yrs, n_seas, n_ff))
  ObsCatch <- array(0, dim = c(n_reg, n_yrs, n_seas, n_ff))
  for (r in 1:n_reg) for (s in 1:n_seas) {
    f <- (r - 1) * n_seas + s
    UseCatch[r, , s, f] <- 1
    ObsCatch[r, , s, f] <- 5
  }

  out <- SPoRC:::do_Fmort_mapping(mk_input(UseCatch, ObsCatch))
  map <- array(as.integer(out$map$ln_F_mean), dim = c(n_reg, n_seas, n_ff))

  expect_equal(nlevels(out$map$ln_F_mean), n_ff)
  for (r in 1:n_reg) for (s in 1:n_seas) for (f in 1:n_ff) {
    is_own <- f == (r - 1) * n_seas + s
    if (is_own) expect_false(is.na(map[r, s, f]))
    else expect_true(is.na(map[r, s, f]))
  }

  # every estimated cell gets its own index -- no unintended sharing
  expect_equal(sort(map[!is.na(map)]), 1:n_ff)
})

test_that("a single fleet fishing everywhere leaves all ln_F_mean estimated", {

  n_reg <- 3; n_seas <- 2; n_yrs <- 10
  UseCatch <- array(1, dim = c(n_reg, n_yrs, n_seas, 1))
  ObsCatch <- array(5, dim = c(n_reg, n_yrs, n_seas, 1))

  out <- SPoRC:::do_Fmort_mapping(mk_input(UseCatch, ObsCatch))

  expect_false(any(is.na(out$map$ln_F_mean)))
  expect_equal(nlevels(out$map$ln_F_mean), n_reg * n_seas)
})

test_that("a cell fished in only one year still gets an estimated mean", {

  n_reg <- 2; n_seas <- 1; n_yrs <- 10; n_ff <- 1
  UseCatch <- array(0, dim = c(n_reg, n_yrs, n_seas, n_ff))
  ObsCatch <- array(0, dim = c(n_reg, n_yrs, n_seas, n_ff))
  UseCatch[1, , 1, 1] <- 1; ObsCatch[1, , 1, 1] <- 5
  UseCatch[2, 4, 1, 1] <- 1; ObsCatch[2, 4, 1, 1] <- 5   # region 2 fished once

  out <- SPoRC:::do_Fmort_mapping(mk_input(UseCatch, ObsCatch))

  expect_equal(nlevels(out$map$ln_F_mean), 2)
  expect_false(any(is.na(out$map$ln_F_mean)))
})

test_that("missing (NA) catch counts as fished, matching the ln_F_devs rule", {

  n_reg <- 2; n_seas <- 1; n_yrs <- 5; n_ff <- 1
  UseCatch <- array(0, dim = c(n_reg, n_yrs, n_seas, n_ff))
  ObsCatch <- array(0, dim = c(n_reg, n_yrs, n_seas, n_ff))
  UseCatch[1, , 1, 1] <- 1; ObsCatch[1, , 1, 1] <- 5
  ObsCatch[2, 3, 1, 1] <- NA  # region 2: no catch fit, but one missing obs

  out <- SPoRC:::do_Fmort_mapping(mk_input(UseCatch, ObsCatch))
  expect_false(any(is.na(out$map$ln_F_mean)))

  # a true closure (recorded zero, never used) stays fixed
  ObsCatch[2, 3, 1, 1] <- 0
  out2 <- SPoRC:::do_Fmort_mapping(mk_input(UseCatch, ObsCatch))
  map2 <- array(as.integer(out2$map$ln_F_mean), dim = c(n_reg, n_seas, n_ff))
  expect_true(is.na(map2[2, 1, 1]))
  expect_false(is.na(map2[1, 1, 1]))
})

test_that("pop-specific catch alone activates ln_F_mean", {

  n_pop <- 2; n_reg <- 2; n_seas <- 1; n_yrs <- 5; n_ff <- 1
  UseCatch <- array(0, dim = c(n_reg, n_yrs, n_seas, n_ff))
  ObsCatch <- array(0, dim = c(n_reg, n_yrs, n_seas, n_ff))
  UseCatch_pop <- array(0, dim = c(n_pop, n_reg, n_yrs, n_seas, n_ff))
  UseCatch_pop[1, 2, , 1, 1] <- 1   # only region 2, only via pop channel

  out <- SPoRC:::do_Fmort_mapping(mk_input(UseCatch, ObsCatch, UseCatch_pop))
  map <- array(as.integer(out$map$ln_F_mean), dim = c(n_reg, n_seas, n_ff))

  expect_true(is.na(map[1, 1, 1]))
  expect_false(is.na(map[2, 1, 1]))
})

test_that("ln_F_mean map matches the year-collapsed ln_F_devs map", {

  set.seed(42)
  n_reg <- 3; n_seas <- 2; n_yrs <- 8; n_ff <- 4
  UseCatch <- array(rbinom(n_reg * n_yrs * n_seas * n_ff, 1, 0.3),
                    dim = c(n_reg, n_yrs, n_seas, n_ff))
  ObsCatch <- array(5, dim = c(n_reg, n_yrs, n_seas, n_ff))

  out <- SPoRC:::do_Fmort_mapping(mk_input(UseCatch, ObsCatch))
  dev_map <- array(as.integer(out$map$ln_F_devs), dim = c(n_reg, n_yrs, n_seas, n_ff))
  mean_map <- array(as.integer(out$map$ln_F_mean), dim = c(n_reg, n_seas, n_ff))

  # a mean is estimated iff at least one dev in that cell is estimated
  dev_active <- apply(!is.na(dev_map), c(1, 3, 4), any)
  expect_equal(!is.na(mean_map), dev_active)
})
