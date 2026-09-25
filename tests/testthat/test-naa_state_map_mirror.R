# Checks map_ln_NAA takes cells out of the numbers at age penalty without changing the rest (1e-12),
# and that a partial blank is refused when the state is correlated.

naa_penalty_setup <- function(seed = 7, ny = 6, na = 4) {
  set.seed(seed)
  d <- c(1, 1, ny, 1, na, 1)
  list(d = d, ny = ny, na = na,
       ln_NAA = array(rnorm(prod(d), 5, 0.3), dim = d),
       NAA_pred = array(exp(rnorm(prod(d), 5, 0.1)), dim = d),
       sigmaNAA = array(0.4, dim = d))
}

naa_penalty <- function(s, mirror, NAA_re = 1) Get_NAA_state_penalty(
  ln_NAA = s$ln_NAA, NAA_pred = s$NAA_pred, sigmaNAA = s$sigmaNAA,
  naa_re_ages = 1:s$na, naa_re_yrs = 1:s$ny, naa_re_seas = 1,
  NAA_re = NAA_re, NAA_pe_pars = array(0, dim = c(1, 1, 3, 1)), map_ln_NAA = mirror,
  NAA_re_region = 0, NAA_region_corr_pars = array(0, dim = c(1, 1, 1)),
  NAA_re_pop = 0, NAA_pop_corr_pars = 0, NAA_re_sex = 0, NAA_sex_corr_pars = 0,
  NAA_re_season = 0, NAA_season_corr_pars = array(0, dim = c(1, 1, 1)),
  naa_re_where = matrix(1, 1, 1))

test_that("no mirror and a mirror keeping every cell give the same penalty", {

  s <- naa_penalty_setup()
  by_hand <- -sum(dnorm(as.vector(s$ln_NAA - log(s$NAA_pred)), 0, 0.4, log = TRUE))

  expect_equal(naa_penalty(s, NULL), by_hand, tolerance = 1e-12)
  expect_equal(naa_penalty(s, array(1, dim = s$d)), by_hand, tolerance = 1e-12)

})

test_that("a blanked cell drops exactly its own contribution", {

  s <- naa_penalty_setup()
  mirror <- array(1, dim = s$d)
  mirror[1,1,5:6,1,3,1] <- NA # last two years of age 3

  eta <- as.vector(s$ln_NAA - log(s$NAA_pred))
  dropped <- -sum(dnorm(eta[(3 - 1) * s$ny + 5:6], 0, 0.4, log = TRUE))

  expect_equal(naa_penalty(s, mirror), naa_penalty(s, array(1, dim = s$d)) - dropped, tolerance = 1e-12)
  expect_equal(naa_penalty(s, array(NA_real_, dim = s$d)), 0, tolerance = 1e-12) # every cell out

})

test_that("a partial blank is refused when the state is correlated", {

  s <- naa_penalty_setup()
  mirror <- array(1, dim = s$d)
  mirror[1,1,5:6,1,3,1] <- NA

  # a cell left out sits inside the joint density, the same reason naa_re_where is refused there
  expect_error(naa_penalty(s, mirror, NAA_re = 3), "cannot be dropped one at a time")
  expect_error(naa_penalty(s, mirror, NAA_re = 2), "cannot be dropped one at a time")

  # but taking every cell out is fine under any form, since the whole penalty goes
  expect_equal(naa_penalty(s, array(NA_real_, dim = s$d), NAA_re = 3), 0, tolerance = 1e-12)

})
