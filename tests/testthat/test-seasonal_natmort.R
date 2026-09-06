library(SPoRC)
library(testthat)

# M has a season dim and holds a rate per year, so mortality in a season is that
# rate times seasdur. Two things to check: one season block is numerically what
# the package did before, and more than one block puts the rate in the season it
# was asked for.

# Helpers -------------------------------------------------------------------

test_that("expand_natmort_seasons holds a rate across seasons and leaves a seasonal array alone", {

  x <- array(stats::runif(2 * 3 * 4 * 5 * 2, 0.1, 0.4), dim = c(2, 3, 4, 5, 2))
  e <- SPoRC:::expand_natmort_seasons(x, 3)

  expect_equal(dim(e), c(2, 3, 4, 3, 5, 2))
  for(seas in 1:3) expect_equal(e[,,,seas,,], x)

  # already has seasons, passes through
  expect_identical(SPoRC:::expand_natmort_seasons(e, 3), e)
  expect_null(SPoRC:::expand_natmort_seasons(NULL, 3))
  expect_error(SPoRC:::expand_natmort_seasons(array(1, dim = c(2, 2)), 3), "must have 6 dimensions")

  # the season dim can go anywhere, which is how the recruitment routines read it
  z <- array(stats::runif(2 * 3 * 5, 0.1, 0.4), dim = c(2, 3, 5))
  g <- SPoRC:::expand_natmort_seasons(z, 4, seas_dim = 3, n_dim = 4)
  expect_equal(dim(g), c(2, 3, 4, 5))
  for(seas in 1:4) expect_equal(array(g[,,seas,,drop = FALSE], dim = dim(z)), z)
})


test_that("collapse_natmort_annual is the duration weighted sum over seasons", {

  seasdur <- c(0.2, 0.5, 0.3)
  y <- array(stats::runif(2 * 3 * 4 * 3 * 5 * 2, 0.1, 0.4), dim = c(2, 3, 4, 3, 5, 2))

  expect_equal(SPoRC:::collapse_natmort_annual(y, seasdur),
               y[,,,1,,] * 0.2 + y[,,,2,,] * 0.5 + y[,,,3,,] * 0.3)

  # constant rate collapses back to itself, since durations sum to one
  x <- array(stats::runif(2 * 3 * 4 * 5 * 2, 0.1, 0.4), dim = c(2, 3, 4, 5, 2))
  expect_equal(SPoRC:::collapse_natmort_annual(SPoRC:::expand_natmort_seasons(x, 3), seasdur), x)

  # seasons can sit anywhere, which is how refpts read it
  z <- array(stats::runif(3 * 3 * 5, 0.1, 0.4), dim = c(3, 3, 5))
  expect_equal(SPoRC:::collapse_natmort_annual(z, seasdur, seas_dim = 2),
               z[,1,] * 0.2 + z[,2,] * 0.5 + z[,3,] * 0.3)

  expect_error(SPoRC:::collapse_natmort_annual(y, c(0.5, 0.5)), "season dimension of natmort")
})


# Block structure -----------------------------------------------------------

test_that("M_seasblk_spec builds the block index and the parameter array over seasons", {

  sim_obj <- seasonal_M_sim()

  const <- suppressWarnings(seasonal_M_input(sim_obj, list(M_spec = "est_ln_M")))
  split <- suppressWarnings(seasonal_M_input(sim_obj, list(M_spec = "est_ln_M",
                                                           M_seasblk_spec = list(1, 2))))
  n_yrs <- seasonal_M_cfg$n_yrs; n_ages <- seasonal_M_cfg$n_ages

  # season dim in slot 4, same as WAA and MatAA
  expect_equal(dim(const$data$M_blocks), c(1, 1, n_yrs, 2, n_ages, 1))
  expect_equal(dim(split$data$M_blocks), c(1, 1, n_yrs, 2, n_ages, 1))

  # one block is one shared par, two blocks are two
  expect_equal(length(const$par$ln_M), 1L)
  expect_equal(length(unique(as.vector(const$data$M_blocks))), 1L)
  expect_equal(length(split$par$ln_M), 2L)
  expect_equal(dim(split$par$ln_M), c(1, 1, 1, 2, 1, 1))

  # each block covers the season it was given
  expect_equal(unique(as.vector(split$data$M_blocks[,,,1,,])), 1L)
  expect_equal(unique(as.vector(split$data$M_blocks[,,,2,,])), 2L)

  # both estimated, not mapped off
  expect_equal(sum(!is.na(split$map$ln_M)), 2L)
})


test_that("season blocks may group seasons, and are validated where they are set", {

  sim_obj <- seasonal_M_sim()

  # both seasons in one block is the default model
  grouped <- suppressWarnings(seasonal_M_input(sim_obj, list(M_spec = "est_ln_M",
                                                             M_seasblk_spec = list(1:2))))
  const <- suppressWarnings(seasonal_M_input(sim_obj, list(M_spec = "est_ln_M")))
  expect_identical(grouped$data$M_blocks, const$data$M_blocks)
  expect_identical(dim(grouped$par$ln_M), dim(const$par$ln_M))

  expect_error(suppressWarnings(seasonal_M_input(sim_obj, list(M_spec = "est_ln_M",
                                                               M_seasblk_spec = "seasonal"))),
               "M_seasblk_spec must be")
})


# Backwards compatibility ---------------------------------------------------

test_that("a five dimensional fixed mortality array is the same model as the six dimensional one", {

  sim_obj <- seasonal_M_sim()

  five <- suppressWarnings(seasonal_M_input(sim_obj, list(M_spec = "fix",
                                                          Fixed_natmort = seasonal_M_fixed(0.3, six_d = FALSE))))
  six <- suppressWarnings(seasonal_M_input(sim_obj, list(M_spec = "fix",
                                                         Fixed_natmort = seasonal_M_fixed(0.3, six_d = TRUE))))

  # 5d is expanded across seasons at setup, so the data lists match before the objective
  expect_equal(dim(five$data$Fixed_natmort), c(1, 1, seasonal_M_cfg$n_yrs, 2, seasonal_M_cfg$n_ages, 1))
  expect_identical(five$data$Fixed_natmort, six$data$Fixed_natmort)

  o5 <- fit_model(five$data, five$par, five$map, random = NULL, silent = TRUE, do_optim = FALSE)
  o6 <- fit_model(six$data, six$par, six$map, random = NULL, silent = TRUE, do_optim = FALSE)
  expect_identical(o5$fn(o5$par), o6$fn(o6$par))
})


test_that("an input list with five dimensional mortality is promoted inside the objective", {

  sim_obj <- seasonal_M_sim()
  il <- suppressWarnings(seasonal_M_input(sim_obj, list(M_spec = "fix",
                                                        Fixed_natmort = seasonal_M_fixed(0.3))))
  ref <- fit_model(il$data, il$par, il$map, random = NULL, silent = TRUE, do_optim = FALSE)

  # older input lists have no season dim on either array
  legacy <- il
  legacy$data$Fixed_natmort <- seasonal_M_fixed(0.3, six_d = FALSE)
  legacy$data$M_blocks <- array(1L, dim = c(1, 1, seasonal_M_cfg$n_yrs, seasonal_M_cfg$n_ages, 1))
  old <- fit_model(legacy$data, legacy$par, legacy$map, random = NULL, silent = TRUE, do_optim = FALSE)

  expect_identical(old$fn(old$par), ref$fn(ref$par))
  expect_equal(dim(old$report(old$par)$natmort), c(1, 1, seasonal_M_cfg$n_yrs, 2, seasonal_M_cfg$n_ages, 1))
})


# What the rate means -------------------------------------------------------

test_that("the mortality applied in a season is the rate times that season's duration", {

  sim_obj <- seasonal_M_sim()
  M_by_seas <- c(0.5, 0.2)

  il <- suppressWarnings(seasonal_M_input(sim_obj, list(M_spec = "fix",
                                                        Fixed_natmort = seasonal_M_fixed(M_by_seas))))
  o <- fit_model(il$data, il$par, il$map, random = NULL, silent = TRUE, do_optim = FALSE)
  r <- o$report(o$par)

  # reported array holds the rate, not the amount
  expect_equal(unique(as.vector(r$natmort[,,,1,,])), M_by_seas[1])
  expect_equal(unique(as.vector(r$natmort[,,,2,,])), M_by_seas[2])

  # and Z holds rate * seasdur
  seasdur <- seasonal_M_cfg$seasdur
  FAA <- apply(r$ret_FAA, 1:6, sum) + apply(r$disc_FAA, 1:6, sum)
  for(seas in 1:2) {
    expect_equal(as.vector(r$ZAA[,,,seas,,] - FAA[,,,seas,,]),
                 rep(M_by_seas[seas] * seasdur[seas], seasonal_M_cfg$n_yrs * seasonal_M_cfg$n_ages))
  }
})


test_that("splitting an annual rate between seasons leaves the start of year numbers unchanged", {

  # this is what lets a stock be bridged on its annual total. Two rates summing
  # to the same annual accumulate the same over the year, so start of year
  # numbers match while within-year numbers don't.
  sim_obj <- seasonal_M_sim()
  seasdur <- seasonal_M_cfg$seasdur
  M_const <- 0.3
  M_split <- c(0.5, (M_const - 0.5 * seasdur[1]) / seasdur[2])
  expect_equal(sum(M_split * seasdur), M_const)

  rep_of <- function(M) {
    il <- suppressWarnings(seasonal_M_input(sim_obj, list(M_spec = "fix", Fixed_natmort = seasonal_M_fixed(M))))
    o <- fit_model(il$data, il$par, il$map, random = NULL, silent = TRUE, do_optim = FALSE)
    o$report(o$par)
  }
  r_const <- rep_of(M_const)
  r_split <- rep_of(M_split)

  # start of year, before any season is stepped
  expect_equal(r_const$NAA[,,,1,,], r_split$NAA[,,,1,,])

  # season 2 differs, which is what a within-year obs sees
  expect_false(isTRUE(all.equal(r_const$NAA[,,,2,,], r_split$NAA[,,,2,,])))

  # so the split isn't a relabelling, it changes the fit
  expect_false(isTRUE(all.equal(r_const$jnLL, r_split$jnLL)))
})


# Guard rails ---------------------------------------------------------------

test_that("season blocks warn on a single season model, where they cannot do anything", {

  input_list <- Setup_Mod_Dim(
    years = 1:10,
    ages = 1:5,
    lens = NULL,
    n_regions = 1,
    n_sexes = 1,
    n_fish_fleets = 1,
    n_srv_fleets = 1,
    n_pop = 1,
    natal_region = 1,
    verbose = FALSE
  )
  biol <- array(1, dim = c(1, 1, 10, 1, 5, 1))

  expect_warning(
    Setup_Mod_Biologicals(
      input_list = input_list,
      WAA = biol,
      MatAA = biol,
      AgeingError = diag(5),
      M_spec = "est_ln_M",
      M_seasblk_spec = list(1)
    ),
    "single season"
  )
})


test_that("a fixed mortality array of the wrong shape is rejected by name", {

  input_list <- Setup_Mod_Dim(
    years = 1:10,
    ages = 1:5,
    lens = NULL,
    n_regions = 1,
    n_sexes = 1,
    n_fish_fleets = 1,
    n_srv_fleets = 1,
    n_seas = 2,
    seasdur = c(0.4, 0.6),
    n_pop = 1,
    natal_region = 1,
    verbose = FALSE
  )
  biol <- array(1, dim = c(1, 1, 10, 2, 5, 1))

  expect_error(
    Setup_Mod_Biologicals(input_list = input_list, WAA = biol, MatAA = biol,
                          AgeingError = diag(5), M_spec = "fix",
                          # three seasons in a two season model
                          Fixed_natmort = array(0.3, dim = c(1, 1, 10, 3, 5, 1))),
    "Fixed_natmort"
  )
})


# Estimation ----------------------------------------------------------------

fit_seasonal_M <- function(M_true, sigmaR, idx_se, iss, seed = 411) {
  sim_obj <- seasonal_M_sim(M_by_seas = M_true, seed = seed, sigmaR = sigmaR, idx_se = idx_se, iss = iss)
  il <- suppressWarnings(seasonal_M_input(sim_obj, list(M_spec = "est_ln_M", M_seasblk_spec = list(1, 2)),
                                          sigmaR = sigmaR))
  il$par$ln_M[] <- log(0.3) # start away from both truths, and at one value for both seasons
  fit <- fit_model(il$data, il$par, il$map, random = NULL, silent = TRUE, do_optim = TRUE, newton_loops = 1)
  list(fit = fit,
       M_hat = c(unique(as.vector(fit$rep$natmort[,,,1,,])), unique(as.vector(fit$rep$natmort[,,,2,,]))))
}


test_that("two season blocks recover the rate in each season when the data are near noiseless", {

  M_true <- c(0.45, 0.2)
  out <- fit_seasonal_M(M_true, sigmaR = 0.02, idx_se = 0.01, iss = 5000)

  expect_lt(max(abs(out$fit$gr(out$fit$optim$par))), 1e-3)

  # both seasons come back, so the block structure puts the rate the user asked
  # for into the season they asked for it in
  expect_equal(out$M_hat[1], M_true[1], tolerance = 0.05)
  expect_equal(out$M_hat[2], M_true[2], tolerance = 0.05)

  # and the two blocks separate rather than collapsing onto the shared start
  expect_gt(out$M_hat[1] - out$M_hat[2], 0.15)
})


test_that("at realistic noise the annual total is identified but the seasonal split is not", {

  # The only thing separating the seasons is the change in numbers within the
  # year, which just the seasonal comps and index see. The annual total is
  # informed by everything, so it comes back much more sharply than the split.
  # This is the confounding Setup_Mod_Biologicals warns about. If a change ever
  # seems to sharpen the split at realistic noise, go look at it.
  M_true <- c(0.45, 0.2)
  seasdur <- seasonal_M_cfg$seasdur
  out <- fit_seasonal_M(M_true, sigmaR = 0.3, idx_se = 0.1, iss = 300)

  expect_lt(max(abs(out$fit$gr(out$fit$optim$par))), 1e-3)

  # the duration weighted annual total lands within a few percent
  expect_equal(sum(out$M_hat * seasdur), sum(M_true * seasdur), tolerance = 0.05)

  # while at least one of the seasonal rates is well outside that
  expect_gt(max(abs(out$M_hat - M_true) / M_true), 0.1)
})


# Priors --------------------------------------------------------------------

test_that("the M prior reads a season block, and defaults to the first season without one", {

  sim_obj <- seasonal_M_sim()

  # a prior written before mortality kept seasons names no season block
  legacy_prior <- data.frame(
    popblk = 1,
    regionblk = 1,
    yearblk = 1,
    ageblk = 1,
    sexblk = 1,
    mu = 0.3,
    sd = 0.1
  )
  seas_prior <- cbind(legacy_prior, seasblk = 1)

  build <- function(prior) {
    suppressWarnings(seasonal_M_input(sim_obj, list(
      M_spec = "est_ln_M",
      Use_M_prior = 1,
      M_prior = prior
    )))
  }
  a <- build(legacy_prior)
  b <- build(seas_prior)

  oa <- fit_model(a$data, a$par, a$map, random = NULL, silent = TRUE, do_optim = FALSE)
  ob <- fit_model(b$data, b$par, b$map, random = NULL, silent = TRUE, do_optim = FALSE)

  # with one block, naming the season or not is the same prior
  expect_identical(oa$fn(oa$par), ob$fn(ob$par))
  expect_gt(oa$report(oa$par)$M_nLL, 0)

  # a season that doesn't exist gets rejected where the prior is set
  expect_error(build(cbind(legacy_prior, seasblk = 3)), "seasblk")
})


test_that("the prior lands on the season block it names when the seasons are split", {

  sim_obj <- seasonal_M_sim()
  base <- data.frame(popblk = 1, regionblk = 1, yearblk = 1, ageblk = 1, sexblk = 1)

  # tight prior far from the start, one season at a time
  nLL_on <- function(seas) {
    il <- suppressWarnings(seasonal_M_input(
      sim_obj, list(
        M_spec = "est_ln_M",
        M_seasblk_spec = list(1, 2),
        Use_M_prior = 1,
        M_prior = cbind(base, seasblk = seas, mu = 0.9, sd = 0.01)
      )))
    il$par$ln_M[] <- log(c(0.9, 0.3)) # season one sits at the prior mean, season two does not
    o <- fit_model(il$data, il$par, il$map, random = NULL, silent = TRUE, do_optim = FALSE)
    o$report(o$par)$M_nLL
  }

  # season 1 sits at the prior mean so the penalty is near its floor, season 2
  # doesn't
  expect_lt(nLL_on(1), nLL_on(2))
})
