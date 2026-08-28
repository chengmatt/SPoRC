# The at-age streams are stored over regions and sexes whatever a fleet reports,
# and each fleet names the margins it sums over, the density it uses, where its
# observation error comes from and how its ages covary. These tests pin those
# four choices against what the model actually computes.

library(SPoRC)
library(testthat)


# Aggregation margins -----------------------------------------------------

test_that("a sex-split stream reads each sex's own catch", {
  ny <- 20; na <- 5
  d <- c(1, ny, 1, na, 2, 1)
  il <- build_at_age(n_sexes = 2, ObsCatchAA = array(1e3, dim = d),
                     UseCatchAA = array(1, dim = d), CatchAA_Type = "spltRspltS")
  r <- at_age_rep(il)

  expect_true(is.finite(r$jnLL))
  expect_equal(dim(r$CatchAA_nLL), c(1, ny, 1, na, 2, 1))
  for(s in 1:2) {
    expect_equal(as.numeric(r$PredCatchAA[1, 10, 1, , s, 1]),
                 as.numeric(r$CAA[1, 1, 10, 1, , s, 1]))
  }
})

test_that("a sex-summed stream reads the sum over sexes, in sex slot one", {
  ny <- 20; na <- 5
  d <- c(1, ny, 1, na, 2, 1)
  use <- array(0, dim = d); use[, , , , 1, ] <- 1
  il <- build_at_age(n_sexes = 2, ObsCatchAA = array(2e3, dim = d), UseCatchAA = use,
                     CatchAA_Type = "spltRaggS")
  r <- at_age_rep(il)

  expect_equal(as.numeric(r$PredCatchAA[1, 10, 1, , 1, 1]),
               as.numeric(apply(r$CAA[1, 1, 10, 1, , , 1], 1, sum)))
  expect_true(all(r$PredCatchAA[, , , , 2, ] == 0))
})

test_that("a region-summed stream reads the sum over regions, in region one", {
  ny <- 20; na <- 5
  d <- c(2, ny, 1, na, 1, 1)
  use <- array(0, dim = d); use[1, , , , , ] <- 1
  il <- build_at_age(n_regions = 2, ObsCatchAA = array(2e3, dim = d), UseCatchAA = use,
                     CatchAA_Type = "agg")
  r <- at_age_rep(il)

  expect_equal(as.numeric(r$PredCatchAA[1, 10, 1, , 1, 1]),
               as.numeric(apply(r$CAA[1, , 10, 1, , , 1], 2, sum)))
  expect_true(all(r$PredCatchAA[2, , , , , ] == 0))
})

test_that("flagging a summed margin outside slot one is refused", {
  ny <- 20; na <- 5
  dr <- c(2, ny, 1, na, 1, 1)
  expect_error(build_at_age(n_regions = 2, ObsCatchAA = array(1, dim = dr),
                            UseCatchAA = array(1, dim = dr), CatchAA_Type = "agg"),
               "sums over regions")
  ds <- c(1, ny, 1, na, 2, 1)
  expect_error(build_at_age(n_sexes = 2, ObsCatchAA = array(1, dim = ds),
                            UseCatchAA = array(1, dim = ds), CatchAA_Type = "spltRaggS"),
               "sums over sexes")
  expect_error(build_at_age(n_sexes = 2, ObsCatchAA = array(1, dim = ds),
                            UseCatchAA = array(1, dim = ds), CatchAA_Type = "nope"),
               "Valid options: agg, spltRaggS, aggRspltS, spltRspltS")
})

test_that("an array missing the sex margin is refused rather than reinterpreted", {
  ny <- 20; na <- 5
  flat <- c(1, ny, 1, na, 1)                       # no sex margin
  full <- c(1, ny, 1, na, 2, 1)

  expect_error(build_at_age(n_sexes = 2, ObsCatchAA = array(2e3, dim = flat),
                            UseCatchAA = array(1, dim = flat)),
               "not the correct dimension")
  # the standard errors are held to the same shape as the observations
  expect_error(build_at_age(n_sexes = 2, ObsCatchAA = array(2e3, dim = full),
                            UseCatchAA = array(1, dim = full),
                            ObsCatchAA_SE = array(0.2, dim = flat)),
               "not the correct dimension")
  # a sex-summed stream says so with a Type, and still carries the full array
  use <- array(0, dim = full); use[, , , , 1, ] <- 1
  il <- build_at_age(n_sexes = 2, ObsCatchAA = array(2e3, dim = full), UseCatchAA = use)
  expect_equal(dim(il$data$UseCatchAA), full)
  expect_equal(sum(il$data$UseCatchAA[, , , , 2, ]), 0)
})

test_that("the observation error key and its starting values need the sex margin", {
  ny <- 20; na <- 5
  d <- c(1, ny, 1, na, 2, 1)
  obs <- array(1e3, dim = d); use <- array(1, dim = d)

  expect_error(build_at_age(n_sexes = 2, ObsCatchAA = obs, UseCatchAA = use,
                            CatchAA_Type = "spltRspltS", sigmaCAA_key = array(1:na, dim = c(na, 1))),
               "key is not the correct dimension")
  expect_error(build_at_age(n_sexes = 2, ObsCatchAA = obs, UseCatchAA = use,
                            CatchAA_Type = "spltRspltS", ln_sigmaCAA = array(log(0.3), dim = c(na, 1))),
               "not the correct dimension")
})


# Observation error keys --------------------------------------------------

test_that("the observation error key carries a sex margin and can share over it", {
  ny <- 20; na <- 5
  d <- c(1, ny, 1, na, 2, 1)
  obs <- array(1e3, dim = d); use <- array(1, dim = d)

  # a key repeating its entries across the sexes couples them together
  il <- build_at_age(n_sexes = 2, ObsCatchAA = obs, UseCatchAA = use, CatchAA_Type = "spltRspltS",
                     sigmaCAA_key = array(rep(1:na, times = 2), dim = c(na, 2, 1)))
  expect_equal(dim(il$par$ln_sigmaCAA), c(na, 2L, 1L))
  expect_equal(length(unique(stats::na.omit(as.integer(il$map$ln_sigmaCAA)))), na)

  # given the margin, each sex can carry its own
  key <- array(0L, dim = c(na, 2, 1))
  key[, 1, 1] <- 1:na; key[, 2, 1] <- na + (1:na)
  il2 <- build_at_age(n_sexes = 2, ObsCatchAA = obs, UseCatchAA = use, CatchAA_Type = "spltRspltS",
                      sigmaCAA_key = key)
  expect_equal(length(unique(stats::na.omit(as.integer(il2$map$ln_sigmaCAA)))), 2 * na)
})

test_that("a sex a fleet never observes carries no observation error parameter", {
  ny <- 20; na <- 5
  d <- c(1, ny, 1, na, 2, 1)
  use <- array(0, dim = d); use[, , , , 1, ] <- 1
  il <- build_at_age(n_sexes = 2, ObsCatchAA = array(2e3, dim = d), UseCatchAA = use,
                     CatchAA_Type = "spltRaggS")
  m <- array(as.integer(il$map$ln_sigmaCAA), dim = c(na, 2, 1))
  expect_true(all(!is.na(m[, 1, 1])))
  expect_true(all(is.na(m[, 2, 1])))
})


# Density and where the error comes from ----------------------------------

test_that("a normal stream is fit on the natural scale", {
  ny <- 20; na <- 5
  d <- c(1, ny, 1, na, 1, 1)
  sig <- 0.3
  il <- build_at_age(ObsCatchAA = array(2e3, dim = d), UseCatchAA = array(1, dim = d),
                     CatchAA_LikeType = "normal", sigmaCAA_spec = "fix",
                     ln_sigmaCAA = array(log(sig), dim = c(na, 1, 1)))
  r <- at_age_rep(il)
  expect_equal(as.numeric(r$CatchAA_nLL[1, 10, 1, , 1, 1]),
               -stats::dnorm(2e3, as.numeric(r$PredCatchAA[1, 10, 1, , 1, 1]), sig, log = TRUE))
})

test_that("reported standard errors enter the way the aggregated index lets them", {
  ny <- 20; na <- 5
  d <- c(1, ny, 1, na, 1, 1)
  obs <- array(2e3, dim = d); use <- array(1, dim = d); se <- array(0.4, dim = d)
  sig <- 0.3

  # the reported errors alone, which leaves the parameter with nothing to read
  il <- build_at_age(ObsCatchAA = obs, UseCatchAA = use, ObsCatchAA_SE = se,
                     CatchAA_sigma_form = "data")
  r <- at_age_rep(il)
  expect_equal(as.numeric(r$CatchAA_nLL[1, 10, 1, , 1, 1]),
               -stats::dnorm(log(2e3), log(as.numeric(r$PredCatchAA[1, 10, 1, , 1, 1])), 0.4, log = TRUE))
  expect_true(all(is.na(as.integer(il$map$ln_sigmaCAA))))

  # and both, added or in quadrature
  for(form in c("est_additive", "est_quadrature")) {
    ilf <- build_at_age(ObsCatchAA = obs, UseCatchAA = use, ObsCatchAA_SE = se,
                        CatchAA_sigma_form = form, sigmaCAA_spec = "fix",
                        ln_sigmaCAA = array(log(sig), dim = c(na, 1, 1)))
    rf <- at_age_rep(ilf)
    want_sd <- if(form == "est_additive") 0.4 + sig else sqrt(0.4^2 + sig^2)
    expect_equal(as.numeric(rf$CatchAA_nLL[1, 10, 1, , 1, 1]),
                 -stats::dnorm(log(2e3), log(as.numeric(rf$PredCatchAA[1, 10, 1, , 1, 1])),
                               want_sd, log = TRUE))
  } # end form loop

  expect_error(build_at_age(ObsCatchAA = obs, UseCatchAA = use, CatchAA_sigma_form = "nope"),
               "Valid options: none, data, est_additive, est_quadrature")
  expect_error(build_at_age(ObsCatchAA = obs, UseCatchAA = use, CatchAA_LikeType = "poisson"),
               "Valid options: lognormal, normal")
})


# Correlation across ages -------------------------------------------------

test_that("an AR(1) across ages is spaced by age, not by position", {
  ny <- 20; na <- 6
  d <- c(1, ny, 1, na, 1, 1)
  ages <- c(1, 2, 4, 5)                      # a gap between the second and third
  sig <- 0.3; rho <- 0.5
  use <- array(0, dim = d); use[1, , 1, ages, 1, 1] <- 1
  obs <- array(0, dim = d); obs[1, , 1, ages, 1, 1] <- 2e3

  il <- build_at_age(n_ages = na, ObsCatchAA = obs, UseCatchAA = use, AgeObsCorr_catch = "1dar1",
                     sigmaCAA_spec = "fix", ln_sigmaCAA = array(log(sig), dim = c(na, 1, 1)),
                     trans_rho_catch = array(atanh(rho), dim = c(1, 1)))
  r <- at_age_rep(il)
  resid <- log(2e3) - log(as.numeric(r$PredCatchAA[1, 10, 1, ages, 1, 1]))

  spaced <- (rho^abs(outer(ages, ages, "-"))) * sig^2
  packed <- (rho^abs(outer(seq_along(ages), seq_along(ages), "-"))) * sig^2
  expect_equal(sum(r$CatchAA_nLL[1, 10, 1, , 1, 1]),
               -dmvn_ref(resid, spaced))
  expect_false(isTRUE(all.equal(sum(r$CatchAA_nLL[1, 10, 1, , 1, 1]),
                                -dmvn_ref(resid, packed))))
})

test_that("an unstructured correlation matches its matrix and is guarded", {
  ny <- 20; na <- 6
  d <- c(1, ny, 1, na, 1, 1)
  sig <- 0.3
  set.seed(7); pars <- stats::rnorm(na * (na - 1) / 2, 0, 0.3)

  il <- build_at_age(n_ages = na, ObsCatchAA = array(2e3, dim = d), UseCatchAA = array(1, dim = d),
                     AgeObsCorr_catch = "us", sigmaCAA_spec = "fix",
                     ln_sigmaCAA = array(log(sig), dim = c(na, 1, 1)),
                     trans_rho_catch_us = array(pars, dim = c(length(pars), 1, 1)))
  r <- at_age_rep(il)
  resid <- log(2e3) - log(as.numeric(r$PredCatchAA[1, 10, 1, , 1, 1]))
  R <- build_us_corr(pars, na)

  expect_equal(sum(r$CatchAA_nLL[1, 10, 1, , 1, 1]),
               -dmvn_ref(resid, R * sig^2))
  expect_equal(length(unique(stats::na.omit(as.integer(il$map$trans_rho_catch_us)))),
               na * (na - 1) / 2)

  # more correlations than cells is refused rather than estimated
  few <- c(1, 5, 1, na, 1, 1)
  expect_error(build_at_age(n_yrs = 5, n_ages = na, ObsCatchAA = array(2e3, dim = few),
                            UseCatchAA = array(1, dim = few), AgeObsCorr_catch = "us"),
               "needs more cells than it has parameters")
})

test_that("a separable correlation over ages and years needs a complete grid", {
  ny <- 8; na <- 4
  d <- c(1, ny, 1, na, 1, 1)
  sig <- 0.3; rho_a <- 0.4; rho_y <- 0.6

  il <- build_at_age(n_yrs = ny, n_ages = na, ObsCatchAA = array(2e3, dim = d),
                     UseCatchAA = array(1, dim = d), AgeObsCorr_catch = "2dar1",
                     sigmaCAA_spec = "fix", ln_sigmaCAA = array(log(sig), dim = c(na, 1, 1)),
                     trans_rho_catch = array(atanh(rho_a), dim = c(1, 1)),
                     trans_rho_catch_year = array(atanh(rho_y), dim = c(1, 1)))
  r <- at_age_rep(il)
  resid <- log(2e3) - log(as.numeric(r$PredCatchAA[1, , 1, , 1, 1]))   # year varies fastest
  Ra <- rho_a^abs(outer(1:na, 1:na, "-")); Ry <- rho_y^abs(outer(1:ny, 1:ny, "-"))

  expect_equal(sum(r$CatchAA_nLL[1, , 1, , 1, 1]),
               -dmvn_ref(resid, kronecker(Ra, Ry) * sig^2))
  expect_false(all(is.na(as.integer(il$map$trans_rho_catch_year))))

  gap <- array(1, dim = d); gap[1, 3, 1, 2, 1, 1] <- 0
  expect_error(build_at_age(n_yrs = ny, n_ages = na, ObsCatchAA = array(2e3, dim = d),
                            UseCatchAA = gap, AgeObsCorr_catch = "2dar1"),
               "complete grid")
})

test_that("the correlation structure is chosen per fleet", {
  ny <- 20; na <- 5
  d <- c(1, ny, 1, na, 1, 2)
  il <- build_at_age(n_fleets = 2, ObsCatchAA = array(2e3, dim = d), UseCatchAA = array(1, dim = d),
                     AgeObsCorr_catch = c("iid", "1dar1"))

  expect_equal(as.numeric(il$data$AgeObsCorr_catch), c(0, 1))
  expect_true(is.na(as.integer(il$map$trans_rho_catch)[1]))
  expect_false(is.na(as.integer(il$map$trans_rho_catch)[2]))

  r <- at_age_rep(il)
  # an independent fleet contributes age by age; a correlated one contributes once
  expect_equal(sum(r$CatchAA_nLL[1, 10, 1, , 1, 1] != 0), na)
  expect_equal(sum(r$CatchAA_nLL[1, 10, 1, , 1, 2] != 0), 1)

  expect_error(build_at_age(n_fleets = 2, ObsCatchAA = array(1, dim = d), UseCatchAA = array(1, dim = d),
                            AgeObsCorr_catch = c("iid", "iid", "iid")),
               "Supply one setting per fleet")
})

test_that("correlations share through the package's spec strings", {
  ny <- 20; na <- 5
  d <- c(2, ny, 1, na, 2, 2)
  obs <- array(1e3, dim = d); use <- array(1, dim = d)
  mk <- function(corr, spec) build_at_age(n_regions = 2, n_sexes = 2, n_fleets = 2,
                                          ObsCatchAA = obs, UseCatchAA = use,
                                          CatchAA_Type = "spltRspltS",
                                          AgeObsCorr_catch = corr, rho_catch_spec = spec)
  n_est <- function(corr, spec) length(unique(stats::na.omit(as.integer(mk(corr, spec)$map$trans_rho_catch))))

  # the correlations sit over region, sex and fleet
  expect_equal(dim(mk("1dar1", NULL)$par$trans_rho_catch), c(2L, 2L, 2L))

  expect_equal(n_est("1dar1", NULL), 2)                  # default: one per fleet
  expect_equal(n_est("1dar1", "est_shared_r_s"), 2)      # which is what that spec says
  expect_equal(n_est("1dar1", "est_shared_s"), 4)        # region by fleet
  expect_equal(n_est("1dar1", "est_shared_r"), 4)        # sex by fleet
  expect_equal(n_est("1dar1", "est_shared_r_s_f"), 1)    # one correlation throughout
  expect_equal(n_est("1dar1", "est_all"), 8)
  expect_true(all(is.na(as.integer(mk("1dar1", "fix")$map$trans_rho_catch))))

  # the same spec shares whole matrices under "us", each pair its own parameter
  n_pairs <- na * (na - 1) / 2
  us_pars <- function(spec) length(unique(stats::na.omit(as.integer(mk("us", spec)$map$trans_rho_catch_us))))
  expect_equal(us_pars("est_shared_r_s"), 2 * n_pairs)   # one matrix per fleet
  expect_equal(us_pars("est_shared_r_s_f"), n_pairs)     # one matrix throughout
  expect_equal(us_pars("est_all"), 8 * n_pairs)
  expect_equal(us_pars("fix"), 0)

  # a dimension the correlation does not have is refused by name
  expect_error(mk("1dar1", "est_shared_y"), "not recognized")
})

test_that("a margin a fleet never observes carries no correlation", {
  ny <- 20; na <- 5
  d <- c(2, ny, 1, na, 2, 1)

  # summed over regions: the observation lives in region one, so region two has
  # nothing to inform a correlation with
  ur <- array(0, dim = d); ur[1, , , , 1, ] <- 1
  il <- build_at_age(n_regions = 2, n_sexes = 2, ObsCatchAA = array(2e3, dim = d), UseCatchAA = ur,
                     CatchAA_Type = "agg", AgeObsCorr_catch = "1dar1", rho_catch_spec = "est_all")
  m <- array(as.integer(il$map$trans_rho_catch), dim = c(2, 2, 1))
  expect_false(is.na(m[1, 1, 1]))
  expect_true(all(is.na(m[2, , 1])))   # region two
  expect_true(all(is.na(m[, 2, 1])))   # sex two

  # and the population-specific stream carries a population margin of its own
  dp <- c(1, d)
  il2 <- build_at_age(n_regions = 2, n_sexes = 2, ObsCatchAA_pop = array(2e3, dim = dp),
                      UseCatchAA_pop = array(1, dim = dp), CatchAA_pop_Type = "spltRspltS",
                      AgeObsCorr_catch_pop = "1dar1")
  expect_equal(dim(il2$par$trans_rho_catch_pop), c(1L, 2L, 2L, 1L))
  expect_true(is.finite(at_age_rep(il2)$jnLL))
})

test_that("the population-specific stream carries its own correlation", {
  ny <- 20; na <- 5
  d <- c(1, ny, 1, na, 1, 1)
  dp <- c(1, d)
  il <- build_at_age(ObsCatchAA_pop = array(2e3, dim = dp), UseCatchAA_pop = array(1, dim = dp),
                     AgeObsCorr_catch = "iid", AgeObsCorr_catch_pop = "1dar1")

  expect_equal(as.numeric(il$data$AgeObsCorr_catch), 0)
  expect_equal(as.numeric(il$data$AgeObsCorr_catch_pop), 1)
  expect_true(all(is.na(as.integer(il$map$trans_rho_catch))))
  expect_false(all(is.na(as.integer(il$map$trans_rho_catch_pop))))

  r <- at_age_rep(il)
  expect_true(is.finite(r$jnLL))
  expect_equal(sum(r$CatchAA_pop_nLL[1, 1, 10, 1, , 1, 1] != 0), 1)  # one correlated cell
})


# Discards ----------------------------------------------------------------

test_that("discards at age are the total discarded, in their own units", {
  ny <- 20; na <- 5
  d <- c(1, ny, 1, na, 1, 1)
  il <- build_at_age(ObsCatchAA = array(2e3, dim = d), UseCatchAA = array(1, dim = d),
                     ObsDiscardAA = array(50, dim = d), UseDiscardAA = array(1, dim = d),
                     discard_units = "abd", dmr_spec = "fix")
  r <- at_age_rep(il)

  # the dead discards the model tracks, raised by the discard mortality rate,
  # which is what the aggregated discard stream is also stated as
  expect_equal(as.numeric(r$PredDiscardAA[1, 10, 1, , 1, 1]),
               as.numeric(apply(r$DAA[1, 1, 10, 1, , , 1, drop = FALSE], 5, sum)) / r$dmr[1, 10, 1, 1])

  # a fraction of the catch is not a statement about one age
  expect_error(build_at_age(ObsDiscardAA = array(50, dim = d), UseDiscardAA = array(1, dim = d),
                            discard_units = "biom_frac"),
               "property of the catch as a whole")
})


# Simulation --------------------------------------------------------------

sim_at_age_om <- function(type, n_sexes = 2, n_sims = 40, seed = 42) {

  sl <- Setup_Sim_Dim(n_sims = n_sims, n_yrs = 20, n_regions = 1, n_ages = 8, n_lens = NULL,
                      n_sexes = n_sexes, n_fish_fleets = 1, n_srv_fleets = 1, n_pop = 1)
  sl <- Setup_Sim_Containers(sl)

  use <- array(0, dim = c(sl$n_regions, sl$n_yrs, sl$n_seas, sl$n_ages, sl$n_sexes, sl$n_fish_fleets))
  if(type == "spltRspltS") use[] <- 1 else use[, , , , 1, ] <- 1

  sl <- Setup_Sim_Fishing(
    sim_list = sl,
    fish_sel_input = replicate(n = sl$n_sims,
      array(rep(1 / (1 + exp(-1.5 * ((1:sl$n_ages) - 3))), each = sl$n_yrs),
            dim = c(sl$n_pop, sl$n_regions, sl$n_yrs, sl$n_seas, sl$n_ages, sl$n_sexes, sl$n_fish_fleets))),
    UseCatchAA = use, use_catch_aa = 1, CatchAA_Type = type,
    ln_sigmaCAA = array(log(0.2), dim = c(sl$n_ages, sl$n_sexes, sl$n_fish_fleets)),
    catch_units = array(0, dim = sl$n_fish_fleets))

  sl <- Setup_Sim_Survey(sim_list = sl,
    srv_sel_input = replicate(n = sl$n_sims,
      array(rep(1 / (1 + exp(-1 * ((1:sl$n_ages) - 3))), each = sl$n_yrs),
            dim = c(sl$n_pop, sl$n_regions, sl$n_yrs, sl$n_seas, sl$n_ages, sl$n_sexes, sl$n_srv_fleets))))

  sl <- suppressWarnings(Setup_Sim_Biologicals(
    sim_list = sl,
    natmort_input = replicate(n = sl$n_sims, array(0.3, dim = c(sl$n_pop, sl$n_regions, sl$n_yrs, sl$n_ages, sl$n_sexes))),
    WAA_input = replicate(n = sl$n_sims, array(1, dim = c(sl$n_pop, sl$n_regions, sl$n_yrs, sl$n_seas, sl$n_ages, sl$n_sexes))),
    WAA_fish_input = replicate(n = sl$n_sims, array(1, dim = c(sl$n_pop, sl$n_regions, sl$n_yrs, sl$n_seas, sl$n_ages, sl$n_sexes, sl$n_fish_fleets))),
    WAA_srv_input = replicate(n = sl$n_sims, array(1, dim = c(sl$n_pop, sl$n_regions, sl$n_yrs, sl$n_seas, sl$n_ages, sl$n_sexes, sl$n_srv_fleets))),
    MatAA_input = replicate(n = sl$n_sims, array(1, dim = c(sl$n_pop, sl$n_regions, sl$n_yrs, sl$n_seas, sl$n_ages, sl$n_sexes)))))

  sl <- Setup_Sim_Tagging(sim_list = sl, use_conv_fish_tagging = 0)
  sl$Movement <- array(1, dim = c(sl$n_pop, sl$n_regions, sl$n_regions, sl$n_yrs, sl$n_seas, sl$n_ages, sl$n_sexes, sl$n_sims))
  sl <- Setup_Sim_Rec(sim_list = sl,
    R0_input = replicate(n = sl$n_sims, array(5, dim = c(sl$n_pop, sl$n_regions, sl$n_yrs))),
    rinit_input = array(2, dim = c(sl$n_pop, sl$n_regions, sl$n_sims)), use_rinit = 1,
    ln_sigmaR = array(log(0.5), dim = c(2, sl$n_pop, sl$n_regions)),
    recruitment_opt = "mean_rec", init_age_strc = 1)

  set.seed(seed)
  Simulate_Pop_Static(sim_list = sl, output_path = NULL)
}

test_that("the operating model draws catch at age sex by sex", {
  om <- sim_at_age_om("spltRspltS")

  # both sexes are drawn, around their own true value
  expect_true(all(om$TrueCatchAA[1, , 1, , 1, 1, ] > 0))
  expect_true(all(om$TrueCatchAA[1, , 1, , 2, 1, ] > 0))
  expect_equal(as.numeric(om$TrueCatchAA[1, 10, 1, , 1, 1, 1]),
               as.numeric(apply(om$CAA[, 1, 10, 1, , 1, 1, 1, drop = FALSE], 5, sum)))

  dev <- log(om$ObsCatchAA[1, , 1, , , 1, ] / om$TrueCatchAA[1, , 1, , , 1, ])
  expect_lt(abs(mean(dev)), 0.02)              # lognormal draws, centered
  expect_equal(stats::sd(as.numeric(dev)), 0.2, tolerance = 0.05)
})

test_that("a sex-summed operating model draws one observation, in sex slot one", {
  om <- sim_at_age_om("spltRaggS")

  expect_true(all(om$TrueCatchAA[1, , 1, , 1, 1, ] > 0))
  expect_true(all(om$TrueCatchAA[1, , 1, , 2, 1, ] == 0))
  expect_equal(as.numeric(om$TrueCatchAA[1, 10, 1, , 1, 1, 1]),
               as.numeric(apply(om$CAA[, 1, 10, 1, , , 1, 1, drop = FALSE], 5, sum)))
})


# One-step-ahead residuals ------------------------------------------------

test_that("at-age OSA residuals carry the age and sex the observation came from", {
  ny <- 12; na <- 4
  d <- c(1, ny, 1, na, 2, 1)
  il <- build_at_age(n_yrs = ny, n_ages = na, n_sexes = 2,
                     ObsCatchAA = array(1e3, dim = d), UseCatchAA = array(1, dim = d),
                     CatchAA_Type = "spltRspltS", sigmaCAA_spec = "fix",
                     ln_sigmaCAA = array(log(0.3), dim = c(na, 2, 1)))
  model <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)

  out <- get_osa(model = model, data = il$data, index_source = "CatchAA")
  expect_equal(nrow(out$res), ny * na * 2)
  expect_true(all(c("age", "sex") %in% names(out$res)))
  expect_setequal(unique(out$res$sex), 1:2)
  expect_setequal(unique(out$res$age), il$data$ages)
  expect_true(all(is.finite(out$res$resid)))
})


# Everything at once ------------------------------------------------------

test_that("streams with different margins, densities and correlations coexist", {
  ny <- 15; na <- 5; nr <- 2; ns <- 2

  caa <- array(500, dim = c(nr, ny, 1, na, ns, 2))   # split region, split sex
  caa_use <- array(1, dim = dim(caa))
  saa <- array(0, dim = c(nr, ny, 1, na, ns, 1))     # split region, summed over sexes
  saa_use <- array(0, dim = dim(saa)); saa_use[, , , , 1, ] <- 1
  saa[saa_use == 1] <- 1e4

  il <- build_at_age(
    n_yrs = ny, n_ages = na, n_regions = nr, n_sexes = ns, n_fleets = 2,
    ObsCatchAA = caa, UseCatchAA = caa_use,
    CatchAA_Type = "spltRspltS",
    CatchAA_LikeType = c("lognormal", "normal"),
    AgeObsCorr_catch = c("1dar1", "iid"),
    srv_extra = list(ObsSrvIdxAA = saa, UseSrvIdxAA = saa_use,
                     SrvIdxAA_Type = "spltRaggS",
                     SrvIdxAA_sigma_form = "est_quadrature",
                     ObsSrvIdxAA_SE = array(0.2, dim = dim(saa)),
                     AgeObsCorr_srv_idx = "us"))

  # each stream stated its own margins, and the flags landed where they belong
  expect_equal(as.numeric(il$data$CatchAA_Type), c(3, 3))
  expect_equal(as.numeric(il$data$SrvIdxAA_Type), 1)
  expect_equal(as.numeric(il$data$CatchAA_LikeType), c(0, 1))
  expect_equal(as.numeric(il$data$AgeObsCorr_catch), c(1, 0))
  expect_equal(as.numeric(il$data$AgeObsCorr_srv_idx), 2)
  expect_equal(as.numeric(il$data$SrvIdxAA_sigma_form), 3)

  model <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)
  r <- model$rep
  expect_true(is.finite(r$jnLL))
  expect_true(all(is.finite(model$gr(model$par))))

  # the catch stream reads each region and sex separately
  expect_equal(as.numeric(r$PredCatchAA[2, 7, 1, , 2, 1]),
               as.numeric(r$CAA[1, 2, 7, 1, , 2, 1]))
  # the survey index keeps its regions apart while summing sexes
  expect_equal(as.numeric(r$PredSrvIdxAA[2, 7, 1, , 1, 1]),
               as.numeric(apply(r$SrvIAA[1, 2, 7, 1, , , 1], 1, sum)))
  expect_true(all(r$PredSrvIdxAA[, , , , 2, ] == 0))
})


# Plotting ----------------------------------------------------------------

test_that("the at-age fits plot carries the sex and the error the fit used", {
  ny <- 12; na <- 4
  d <- c(1, ny, 1, na, 2, 1)
  il <- build_at_age(n_yrs = ny, n_ages = na, n_sexes = 2,
                     ObsCatchAA = array(1e3, dim = d), UseCatchAA = array(1, dim = d),
                     CatchAA_Type = "spltRspltS", ObsCatchAA_SE = array(0.15, dim = d),
                     CatchAA_sigma_form = "est_quadrature")
  p <- get_at_age_fits_plot(list(il$data), list(at_age_rep(il)), "m1", stream = "CatchAA")

  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), ny * na * 2)
  expect_setequal(unique(p$data$Sex), 1:2)
  # the interval is drawn from the standard deviation the likelihood used, not
  # from the estimated component alone
  expect_equal(unique(p$data$sigma), sqrt(0.15^2 + 0.5^2))
})


# Input lists written before the sex margin --------------------------------

test_that("an input list carrying an older shape is refused, not reinterpreted", {
  ny <- 12; na <- 4
  d <- c(1, ny, 1, na, 1, 1)
  il <- build_at_age(n_yrs = ny, n_ages = na, ObsCatchAA = array(1e3, dim = d),
                     UseCatchAA = array(1, dim = d))

  # an array or parameter left at its old shape would be indexed by position and
  # read the wrong age or sex, so the objective refuses it by name
  for(nm in c("ObsCatchAA", "UseCatchAA", "ObsCatchAA_SE")) {
    old <- il
    old$data[[nm]] <- array(1, dim = c(1, ny, 1, na, 1))
    expect_error(fit_model(old$data, old$par, old$map, do_optim = FALSE, silent = TRUE),
                 "carry a sex margin")
  } # end nm loop

  old <- il
  old$par$ln_sigmaCAA <- array(log(0.5), dim = c(na, 1))
  old$map$ln_sigmaCAA <- factor(rep(1, na))
  expect_error(fit_model(old$data, old$par, old$map, do_optim = FALSE, silent = TRUE),
               "carry a sex margin")
})
