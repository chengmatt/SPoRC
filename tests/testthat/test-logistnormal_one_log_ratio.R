# A logistic normal composition left with exactly two positive bins has a single log ratio, so its covariance is 1 x 1.
# Subsetting a matrix down to 1 x 1 drops it to a number, and diag(rep(x, 1)) reads x as a size, so on plain numbers
# these cells used to stop (or take the wrong size) before the density was evaluated.
#
# Checked against a single normal on the one log ratio, written out here, in every logistic normal type, and in the
# simulation draw.

library(SPoRC)
library(testthat)

obs <- c(6, 0, 4, 0) # two positive bins, one log ratio
exp_n <- c(30, 20, 40, 10) # expected numbers
resid <- (log(6) - log(4)) - (log(30) - log(40)) # observed log ratio minus expected
obs_pos <- c(6, 4) / 10 # observed proportions over the positive bins
ISS <- 50 # input sample size
ln_theta <- log(0.6) # log standard deviation
trans_rho <- 0.4 # unconstrained correlation
rho_trans <- function(x) 2 / (1 + exp(-2 * x)) - 1 # types 2 - 4
ln_corr <- array(trans_rho, dim = c(1, 1, 3))

# one region, one sex, ages with an identity ageing error, plain numbers
comp_nll <- function(like_type, comp_type = 0, obs_arr = array(obs, dim = c(1, 4, 1)), exp_arr = array(exp_n, dim = c(1, 4, 1))) {
  n_sexes <- dim(obs_arr)[3]
  SPoRC:::Get_Comp_Likelihoods(Exp = exp_arr,
                               Obs = obs_arr,
                               ISS = array(ISS, dim = c(1, n_sexes)),
                               Wt_Mltnml = array(1, dim = c(1, n_sexes)),
                               ln_theta_agg = ln_theta,
                               ln_theta = array(ln_theta, dim = c(1, n_sexes)),
                               LN_corr_pars = array(trans_rho, dim = c(1, n_sexes, 3)),
                               LN_corr_pars_agg = trans_rho,
                               Comp_Type = comp_type,
                               Likelihood_Type = like_type,
                               n_regions = 1,
                               n_model_bins = 4,
                               n_obs_bins = 4,
                               n_sexes = n_sexes,
                               age_or_len = 0,
                               AgeingError = diag(4),
                               use = 1,
                               addtocomp = 1e-3,
                               comp_bins = NULL,
                               comp_const_obs = 0)
}

test_that("one log ratio gives the single normal in every logistic normal type", {

  theta <- exp(ln_theta)
  rho <- rho_trans(trans_rho)
  var_miss0 <- theta^2 / ISS # the zeros dropped forms scale by the sample size, and their correlation of a bin with itself is 1

  # zeros removed and renormalized, no change of variables
  expect_equal(sum(comp_nll(2)), -dnorm(resid, 0, theta, log = TRUE), tolerance = 1e-12)
  expect_equal(sum(comp_nll(3)), -dnorm(resid, 0, sqrt(theta^2 / (1 - rho^2)), log = TRUE), tolerance = 1e-12)

  # zeros dropped, with the change of variables
  for(like_type in 5:6) {
    expect_equal(sum(comp_nll(like_type)), -(dnorm(resid, 0, sqrt(var_miss0), log = TRUE) - sum(log(obs_pos))), tolerance = 1e-12)
  } # end like_type loop

  # joint across two sexes, one positive bin in each
  obs_joint <- array(c(6, 0, 0, 0, 4, 0, 0, 0), dim = c(1, 4, 2))
  exp_joint <- array(c(30, 20, 40, 10, 40, 30, 20, 10), dim = c(1, 4, 2))
  resid_joint <- (log(6) - log(4)) - (log(30) - log(40))
  expect_equal(sum(comp_nll(4, comp_type = 2, obs_arr = obs_joint, exp_arr = exp_joint)),
               -dnorm(resid_joint, 0, sqrt(theta^2 / (1 - rho^2)^2), log = TRUE), tolerance = 1e-12)
  expect_equal(sum(comp_nll(7, comp_type = 2, obs_arr = obs_joint, exp_arr = exp_joint)),
               -(dnorm(resid_joint, 0, sqrt(var_miss0), log = TRUE) - sum(log(obs_pos))), tolerance = 1e-12)

})

test_that("the simulation draws a two bin composition in every logistic normal type", {

  set.seed(1)
  for(comp_like in 2:7) {
    draw <- SPoRC:::rlogistnormal(exp = c(0.4, 0.6), pars = c(0.5, 0.3, 0.2), comp_like = comp_like, n_sexes = 1, ISS = ISS)
    expect_length(draw, 2)
    expect_equal(sum(draw), 1, tolerance = 1e-12)
    expect_true(all(draw > 0))
  } # end comp_like loop

})
