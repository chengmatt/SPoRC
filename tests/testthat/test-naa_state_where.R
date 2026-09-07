library(SPoRC)
library(testthat)

# A natal homing population holds no fish in a region it never reaches, so it has no
# numbers at age state there and the penalty would otherwise take the logarithm of zero.
# Checks the argument, the map, the penalty and the guard against correlated regions.

# a small state on two populations and two regions, the second population homing to region two
where_pars <- function(naa_re_where = NULL, n_pop = 2, n_regions = 2, ny = 4, na = 3, seed = 8) {

  set.seed(seed)
  d <- c(n_pop, n_regions, ny, 1, na, 1)
  pred <- array(stats::runif(prod(d), 50, 200), dim = d)
  state <- array(log(pred) + stats::rnorm(prod(d), 0, 0.2), dim = d)

  # the second population never occupies the first region, so its numbers there are zero
  if(!is.null(naa_re_where)) for(p in 1:n_pop) for(r in 1:n_regions) {
    if(naa_re_where[p, r] == 0) {
      pred[p, r, , , , ] <- 0
      state[p, r, , , , ] <- log(1e-10)
    }
  }

  list(ln_NAA = state,
       NAA_pred = pred,
       sigmaNAA = array(0.3, dim = d),
       naa_re_yrs = 2:ny,
       naa_re_ages = 2:na,
       naa_re_seas = 1)
}


test_that("a structurally empty cell would make the penalty infinite and does not", {

  occupancy <- matrix(c(1, 0, 1, 1), nrow = 2, ncol = 2)
  p <- where_pars(occupancy)

  # without saying where the state runs, the empty cell takes log(0)
  unguarded <- SPoRC:::Get_NAA_state_penalty(
    ln_NAA = p$ln_NAA, NAA_pred = p$NAA_pred, sigmaNAA = p$sigmaNAA,
    naa_re_ages = p$naa_re_ages, naa_re_yrs = p$naa_re_yrs, naa_re_seas = p$naa_re_seas,
    NAA_re = 1)
  expect_true(is.infinite(unguarded))

  guarded <- SPoRC:::Get_NAA_state_penalty(
    ln_NAA = p$ln_NAA, NAA_pred = p$NAA_pred, sigmaNAA = p$sigmaNAA,
    naa_re_ages = p$naa_re_ages, naa_re_yrs = p$naa_re_yrs, naa_re_seas = p$naa_re_seas,
    NAA_re = 1, naa_re_where = occupancy)
  expect_true(is.finite(guarded))
})


test_that("the penalty is the sum over the cells the state runs on and nothing else", {

  occupancy <- matrix(c(1, 0, 1, 1), nrow = 2, ncol = 2)
  p <- where_pars(occupancy)

  got <- SPoRC:::Get_NAA_state_penalty(
    ln_NAA = p$ln_NAA, NAA_pred = p$NAA_pred, sigmaNAA = p$sigmaNAA,
    naa_re_ages = p$naa_re_ages, naa_re_yrs = p$naa_re_yrs, naa_re_seas = p$naa_re_seas,
    NAA_re = 1, naa_re_where = occupancy)

  want <- 0
  for(pop in 1:2) for(r in 1:2) {
    if(occupancy[pop, r] == 0) next
    eta <- p$ln_NAA[pop, r, p$naa_re_yrs, 1, p$naa_re_ages, 1] -
           log(p$NAA_pred[pop, r, p$naa_re_yrs, 1, p$naa_re_ages, 1])
    want <- want - sum(stats::dnorm(as.vector(eta), 0, 0.3, log = TRUE))
  }

  expect_equal(got, want, tolerance = 1e-10)

  # every cell occupied is the same answer the argument-free call gives
  full <- where_pars(NULL)
  expect_equal(
    SPoRC:::Get_NAA_state_penalty(ln_NAA = full$ln_NAA, NAA_pred = full$NAA_pred,
      sigmaNAA = full$sigmaNAA, naa_re_ages = full$naa_re_ages, naa_re_yrs = full$naa_re_yrs,
      naa_re_seas = full$naa_re_seas, NAA_re = 1, naa_re_where = matrix(1, 2, 2)),
    SPoRC:::Get_NAA_state_penalty(ln_NAA = full$ln_NAA, NAA_pred = full$NAA_pred,
      sigmaNAA = full$sigmaNAA, naa_re_ages = full$naa_re_ages, naa_re_yrs = full$naa_re_yrs,
      naa_re_seas = full$naa_re_seas, NAA_re = 1),
    tolerance = 1e-12)
})


test_that("a correlated form skips the dropped cell too", {

  occupancy <- matrix(c(1, 0, 1, 1), nrow = 2, ncol = 2)
  p <- where_pars(occupancy)
  pe <- array(0.4, dim = c(2, 2, 3, 1))

  for(form in c(2, 3, 4)) { # ar1 over ages, over years, and separable over both
    got <- SPoRC:::Get_NAA_state_penalty(
      ln_NAA = p$ln_NAA, NAA_pred = p$NAA_pred, sigmaNAA = p$sigmaNAA,
      naa_re_ages = p$naa_re_ages, naa_re_yrs = p$naa_re_yrs, naa_re_seas = p$naa_re_seas,
      NAA_re = form, NAA_pe_pars = pe, naa_re_where = occupancy)
    expect_true(is.finite(got))
  } # end form loop
})


test_that("a dropped cell is refused inside a correlation that spans it", {

  occupancy <- matrix(c(1, 0, 1, 1), nrow = 2, ncol = 2)
  p <- where_pars(occupancy)

  expect_error(
    SPoRC:::Get_NAA_state_penalty(
      ln_NAA = p$ln_NAA, NAA_pred = p$NAA_pred, sigmaNAA = p$sigmaNAA,
      naa_re_ages = p$naa_re_ages, naa_re_yrs = p$naa_re_yrs, naa_re_seas = p$naa_re_seas,
      NAA_re = 1, NAA_re_region = 1, NAA_region_corr_pars = array(0.2, dim = c(2, 1, 1)),
      naa_re_where = occupancy),
    "correlates")

  expect_error(
    SPoRC:::Get_NAA_state_penalty(
      ln_NAA = p$ln_NAA, NAA_pred = p$NAA_pred, sigmaNAA = p$sigmaNAA,
      naa_re_ages = p$naa_re_ages, naa_re_yrs = p$naa_re_yrs, naa_re_seas = p$naa_re_seas,
      NAA_re = 1, NAA_re_pop = 1, NAA_pop_corr_pars = array(0.2, dim = 1),
      naa_re_where = occupancy),
    "correlates")
})


test_that("setup drops the cell from the map as well as from the penalty", {

  sim <- seasonal_M_sim()
  base <- seasonal_M_input(sim, list(M_spec = "fix", Fixed_natmort = seasonal_M_fixed(seasonal_M_cfg$M)))

  # one population here, so the only thing to check is that the field lands and the map follows
  expect_equal(dim(base$data$naa_re_where), c(1, 1))
  expect_equal(as.vector(base$data$naa_re_where), 1L)

  expect_error(SPoRC:::do_NAAstate_mapping(base, NAA_re = "iid", NAA_re_ages = NULL,
                                           NAA_re_years = NULL, NAA_sigma_spec = "fix",
                                           NAA_sigma_popblk_spec_vals = list(1),
                                           NAA_sigma_regionblk_spec_vals = list(1),
                                           NAA_sigma_yearblk_spec_vals = list(1),
                                           NAA_sigma_seasblk_spec_vals = list(1),
                                           NAA_sigma_ageblk_spec_vals = list(1),
                                           NAA_sigma_sexblk_spec_vals = list(1),
                                           NAA_re_where = matrix(2, 1, 1)),
               "other than 0 and 1")
})
