library(SPoRC)
library(testthat)

# The logistic normal with the zeros dropped, which renormalizes over the bins that are left,
# scales the variance by the input sample size and takes the change of variables off so the
# result is a density on the composition. Checked against an independent reference over every
# composition type, sex and region layout, ages and lengths, and a restricted bin range.

# Reference ------------------------------------------------------------------

# written from the definition rather than from the package, so a shared mistake cannot pass
ref_miss0 <- function(obs, pred, ln_sigma, ISS, ar1 = FALSE, trans_rho = 0, lag_bins = NULL) {

  pos <- which(obs > 1e-15)
  n_pos <- length(pos)
  if(n_pos < 2) return(0)

  obs_pos <- obs[pos] / sum(obs[pos])
  pred_pos <- pred[pos] / sum(pred[pos])

  y <- log(obs_pos[-n_pos]) - log(obs_pos[n_pos]) # additive log ratio
  mu <- log(pred_pos[-n_pos]) - log(pred_pos[n_pos])

  variance <- exp(2 * ln_sigma) / ISS
  n_tr <- n_pos - 1
  lags <- if(is.null(lag_bins)) pos else lag_bins[pos]

  Sigma <- matrix(0, n_tr, n_tr)
  for(i in 1:n_tr) for(j in 1:n_tr) {
    Sigma[i, j] <- if(ar1) variance * stats::plogis(trans_rho)^abs(lags[i] - lags[j]) else if(i == j) variance else 0
  }

  chol_S <- chol(Sigma)
  z <- backsolve(chol_S, y - mu, transpose = TRUE)
  loglik <- -0.5 * sum(z^2) - sum(log(diag(chol_S))) - 0.5 * n_tr * log(2 * pi)

  -1 * (loglik - sum(log(obs_pos))) # negative log likelihood, jacobian included
}

# a composition with genuine zeros in it
draw_comp <- function(n_bins, sample_size = 60, seed = 1) {
  set.seed(seed)
  prob <- stats::runif(n_bins, 0.05, 1)
  prob <- prob / sum(prob)
  obs <- as.numeric(stats::rmultinom(1, sample_size, prob))
  list(obs = obs / sum(obs), pred = prob)
}


test_that("the density matches the reference with and without zeros, iid and ar1", {

  for(seed in 1:5) {
    cell <- draw_comp(10, sample_size = 40, seed = seed)
    expect_equal(SPoRC:::get_logistnormal_miss0_nLL(cell$obs, cell$pred, ln_sigma = -0.3, ISS = 40),
                 ref_miss0(cell$obs, cell$pred, ln_sigma = -0.3, ISS = 40), tolerance = 1e-10)
    expect_equal(SPoRC:::get_logistnormal_miss0_nLL(cell$obs, cell$pred, ln_sigma = -0.3, ISS = 40,
                                                    corr_type = 1, trans_rho = 0.8),
                 ref_miss0(cell$obs, cell$pred, ln_sigma = -0.3, ISS = 40, ar1 = TRUE, trans_rho = 0.8),
                 tolerance = 1e-10)
  } # end seed loop

  # a cell with one positive bin has no log ratio to take
  expect_equal(SPoRC:::get_logistnormal_miss0_nLL(c(1, 0, 0), c(0.5, 0.3, 0.2), -0.3, 40), 0)
  expect_equal(SPoRC:::get_logistnormal_miss0_nLL(c(0, 0, 0), c(0.5, 0.3, 0.2), -0.3, 40), 0)
})


test_that("the sample size scales the variance and the jacobian is on the positive bins", {

  cell <- draw_comp(8, sample_size = 100, seed = 7)

  # doubling the sample size is the same as halving the variance
  a <- SPoRC:::get_logistnormal_miss0_nLL(cell$obs, cell$pred, ln_sigma = 0, ISS = 200)
  b <- SPoRC:::get_logistnormal_miss0_nLL(cell$obs, cell$pred, ln_sigma = -0.5 * log(2), ISS = 100)
  expect_equal(a, b, tolerance = 1e-10)

  # the form that keeps every bin has no jacobian, so the difference is exactly that term
  pos <- which(cell$obs > 1e-15)
  obs_pos <- cell$obs[pos] / sum(cell$obs[pos])
  pred_pos <- cell$pred[pos] / sum(cell$pred[pos])
  Sigma <- diag(rep(exp(2 * -0.3) / 50, length(pos) - 1))
  no_jac <- -SPoRC:::dlogistnormal(obs_pos, pred_pos, Sigma, TRUE)
  with_jac <- SPoRC:::get_logistnormal_miss0_nLL(cell$obs, cell$pred, ln_sigma = -0.3, ISS = 50)
  expect_equal(with_jac - no_jac, sum(log(obs_pos)), tolerance = 1e-10)
})


test_that("the ar1 is spaced by bin, so a gap in the observed bins is a longer lag", {

  obs <- c(0.4, 0, 0.35, 0.25) # bins 1, 3 and 4 are positive
  pred <- c(0.3, 0.2, 0.3, 0.2)

  got <- SPoRC:::get_logistnormal_miss0_nLL(obs, pred, ln_sigma = -0.2, ISS = 50,
                                            corr_type = 1, trans_rho = 0.9)
  expect_equal(got, ref_miss0(obs, pred, ln_sigma = -0.2, ISS = 50, ar1 = TRUE, trans_rho = 0.9),
               tolerance = 1e-10)

  # spacing by position instead would put bins 1 and 3 at lag 1 rather than lag 2
  by_position <- ref_miss0(obs[obs > 0], pred[obs > 0], ln_sigma = -0.2, ISS = 50,
                           ar1 = TRUE, trans_rho = 0.9)
  expect_false(isTRUE(all.equal(got, by_position)))
})


# Composition types ----------------------------------------------------------

# one cell through Get_Comp_Likelihoods, so the region, sex and bin handling is what is checked
call_comp <- function(Obs, Exp, ISS, ln_theta, corr, comp_type, like_type,
                      n_regions, n_bins, n_sexes, use, age_or_len = 0, comp_bins = NULL) {
  SPoRC:::Get_Comp_Likelihoods(
    Exp = Exp,
    Obs = Obs,
    ISS = ISS,
    Wt_Mltnml = array(1, dim = c(n_regions, n_sexes)),
    ln_theta_agg = ln_theta[1, 1],
    ln_theta = ln_theta,
    LN_corr_pars = corr,
    LN_corr_pars_agg = corr[1, 1, 1],
    Comp_Type = comp_type,
    Likelihood_Type = like_type,
    n_regions = n_regions,
    n_model_bins = n_bins,
    n_obs_bins = n_bins,
    n_sexes = n_sexes,
    age_or_len = age_or_len,
    AgeingError = diag(n_bins),
    use = use,
    addtocomp = 0,
    comp_bins = comp_bins
  )
}


test_that("split by region and sex gives each cell its own density", {

  n_regions <- 2; n_sexes <- 2; n_bins <- 8
  Obs <- array(0, dim = c(n_regions, n_bins, n_sexes))
  Exp <- array(0, dim = c(n_regions, n_bins, n_sexes))
  ISS <- array(0, dim = c(n_regions, n_sexes))
  seed <- 0

  for(r in 1:n_regions) for(s in 1:n_sexes) {
    seed <- seed + 1
    cell <- draw_comp(n_bins, sample_size = 50, seed = seed)
    Obs[r, , s] <- cell$obs
    Exp[r, , s] <- cell$pred
    ISS[r, s] <- 50 + 10 * seed # a different sample size in each cell
  }

  ln_theta <- matrix(c(-0.3, -0.1, 0.2, 0.05), n_regions, n_sexes)
  corr <- array(0.6, dim = c(n_regions, n_sexes, 3))
  use <- rep(1, n_regions)

  for(like_type in c(5, 6)) {
    got <- call_comp(Obs, Exp, ISS, ln_theta, corr, comp_type = 1, like_type = like_type,
                     n_regions = n_regions, n_bins = n_bins, n_sexes = n_sexes, use = use)

    for(r in 1:n_regions) for(s in 1:n_sexes) {
      expect_equal(got[r, s],
                   ref_miss0(Obs[r, , s], Exp[r, , s], ln_theta[r, s], ISS[r, s],
                             ar1 = (like_type == 6), trans_rho = corr[r, s, 1]),
                   tolerance = 1e-9)
    }
  } # end like_type loop
})


test_that("joint by sex puts the whole bin by sex stack in one density", {

  n_regions <- 2; n_sexes <- 2; n_bins <- 6
  Obs <- array(0, dim = c(n_regions, n_bins, n_sexes))
  Exp <- array(0, dim = c(n_regions, n_bins, n_sexes))
  seed <- 10

  for(r in 1:n_regions) {
    seed <- seed + 1
    stack <- draw_comp(n_bins * n_sexes, sample_size = 80, seed = seed)
    Obs[r, , ] <- matrix(stack$obs, n_bins, n_sexes)
    Exp[r, , ] <- matrix(stack$pred, n_bins, n_sexes)
  }

  ISS <- array(c(70, 90), dim = c(n_regions, n_sexes))
  ln_theta <- matrix(c(-0.25, 0.1), n_regions, n_sexes)
  corr <- array(0.4, dim = c(n_regions, n_sexes, 3))
  use <- rep(1, n_regions)

  for(like_type in c(5, 6)) {
    got <- call_comp(Obs, Exp, ISS, ln_theta, corr, comp_type = 2, like_type = like_type,
                     n_regions = n_regions, n_bins = n_bins, n_sexes = n_sexes, use = use)

    for(r in 1:n_regions) {
      expect_equal(got[r, 1],
                   ref_miss0(as.vector(Obs[r, , ]), as.vector(Exp[r, , ]), ln_theta[r, 1], ISS[r, 1],
                             ar1 = (like_type == 6), trans_rho = corr[r, 1, 1]),
                   tolerance = 1e-9)
      expect_equal(got[r, 2], 0) # the second sex slot holds nothing under a joint fit
    }
  } # end like_type loop
})


test_that("aggregated over regions and sexes is one density on the pooled composition", {

  n_regions <- 2; n_sexes <- 1; n_bins <- 7
  cell <- draw_comp(n_bins, sample_size = 120, seed = 21)

  Obs <- array(0, dim = c(n_regions, n_bins, n_sexes))
  Obs[1, , 1] <- cell$obs
  Exp <- array(0, dim = c(n_regions, n_bins, n_sexes))
  for(r in 1:n_regions) Exp[r, , 1] <- cell$pred # the same shape in both, so pooling returns it

  ISS <- array(120, dim = c(n_regions, n_sexes))
  ln_theta <- matrix(-0.35, n_regions, n_sexes)
  corr <- array(0.55, dim = c(n_regions, n_sexes, 3))
  use <- c(1, 0) # aggregated comps sit in the first slot

  for(like_type in c(5, 6)) {
    got <- call_comp(Obs, Exp, ISS, ln_theta, corr, comp_type = 0, like_type = like_type,
                     n_regions = n_regions, n_bins = n_bins, n_sexes = n_sexes, use = use)
    expect_equal(got[1, 1],
                 ref_miss0(cell$obs, cell$pred, ln_theta[1, 1], ISS[1, 1],
                           ar1 = (like_type == 6), trans_rho = corr[1, 1, 1]),
                 tolerance = 1e-9)
  } # end like_type loop
})


test_that("length compositions take the same route as ages", {

  n_regions <- 1; n_sexes <- 1; n_bins <- 9
  cell <- draw_comp(n_bins, sample_size = 45, seed = 31)

  Obs <- array(cell$obs, dim = c(n_regions, n_bins, n_sexes))
  Exp <- array(cell$pred, dim = c(n_regions, n_bins, n_sexes))
  ISS <- array(45, dim = c(n_regions, n_sexes))
  ln_theta <- matrix(-0.2, n_regions, n_sexes)
  corr <- array(0.3, dim = c(n_regions, n_sexes, 3))

  # age_or_len = 1 with no bin map skips the ageing error multiply
  got <- SPoRC:::Get_Comp_Likelihoods(
    Exp = Exp, Obs = Obs, ISS = ISS, Wt_Mltnml = array(1, dim = c(1, 1)),
    ln_theta_agg = ln_theta[1, 1], ln_theta = ln_theta,
    LN_corr_pars = corr, LN_corr_pars_agg = corr[1, 1, 1],
    Comp_Type = 1, Likelihood_Type = 5,
    n_regions = 1, n_model_bins = n_bins, n_obs_bins = n_bins, n_sexes = 1,
    age_or_len = 1, AgeingError = NA, use = 1, addtocomp = 0
  )

  expect_equal(got[1, 1], ref_miss0(cell$obs, cell$pred, ln_theta[1, 1], 45), tolerance = 1e-9)
})


test_that("a restricted bin range keeps the ar1 spaced over the original bins", {

  n_bins <- 10
  fit_bins <- c(3, 4, 5, 8, 9) # a gap between 5 and 8
  cell <- draw_comp(n_bins, sample_size = 200, seed = 41)

  Obs <- array(cell$obs, dim = c(1, n_bins, 1))
  Exp <- array(cell$pred, dim = c(1, n_bins, 1))
  ISS <- array(200, dim = c(1, 1))
  ln_theta <- matrix(-0.15, 1, 1)
  corr <- array(0.85, dim = c(1, 1, 3))

  got <- call_comp(Obs, Exp, ISS, ln_theta, corr, comp_type = 1, like_type = 6,
                   n_regions = 1, n_bins = n_bins, n_sexes = 1, use = 1, comp_bins = fit_bins)

  obs_fit <- cell$obs[fit_bins] / sum(cell$obs[fit_bins])
  pred_fit <- cell$pred[fit_bins] / sum(cell$pred[fit_bins])
  expect_equal(got[1, 1],
               ref_miss0(obs_fit, pred_fit, ln_theta[1, 1], 200, ar1 = TRUE,
                         trans_rho = corr[1, 1, 1], lag_bins = fit_bins),
               tolerance = 1e-9)
})


# Through the model ----------------------------------------------------------

# the seasonal test model, with its fishery compositions moved onto the new likelihood
miss0_model <- function(like_type = "iid-Logistic-Normal-miss0") {

  sim <- seasonal_M_sim()
  input_list <- seasonal_M_input(sim, list(M_spec = "fix", Fixed_natmort = seasonal_M_fixed(seasonal_M_cfg$M)))

  input_list$data$FishAgeComps_LikeType <- if(like_type == "iid-Logistic-Normal-miss0") 5L else 6L
  input_list$par$ln_FishAge_theta[] <- -0.3
  input_list$par$FishAge_theta_agg <- NULL
  input_list$par$ln_FishAge_theta_agg[] <- -0.3
  input_list$par$FishAge_corr_pars_agg[] <- 0.7

  input_list
}


test_that("the objective evaluates and differentiates on the new likelihood", {

  for(lt in c("iid-Logistic-Normal-miss0", "1d-Logistic-Normal-miss0")) {

    input_list <- miss0_model(lt)
    fit <- fit_model(data = input_list$data, parameters = input_list$par, mapping = input_list$map,
                     random = NULL, do_optim = FALSE, silent = TRUE)

    expect_true(is.finite(fit$rep$jnLL))
    expect_true(all(is.finite(fit$gr(fit$par))))
    expect_true(sum(fit$rep$FishAgeComps_nLL) != 0)
  } # end lt loop
})


test_that("the reported composition likelihood is the density the helper gives", {

  input_list <- miss0_model("1d-Logistic-Normal-miss0")
  fit <- fit_model(data = input_list$data, parameters = input_list$par, mapping = input_list$map,
                   random = NULL, do_optim = FALSE, silent = TRUE)

  # the model fits aggregated compositions, so each year is one density on the pooled vector
  n_yrs <- seasonal_M_cfg$n_yrs
  n_ages <- seasonal_M_cfg$n_ages
  ln_theta <- as.numeric(input_list$par$ln_FishAge_theta_agg[1])
  trans_rho <- as.numeric(input_list$par$FishAge_corr_pars_agg[1])

  for(y in c(1, 5, 12, n_yrs)) {
    seas <- which(apply(input_list$data$UseFishAgeComps[, y, , 1, drop = FALSE], 3, sum) > 0)[1]
    obs <- input_list$data$ObsFishAgeComps[1, y, seas, , 1, 1]
    obs <- obs / sum(obs)
    pred <- apply(fit$rep$CAA[, , y, seas, , , 1, drop = FALSE], 5, sum)
    pred <- pred / sum(pred)
    iss <- input_list$data$ISS_FishAgeComps[1, y, seas, 1, 1]

    expect_equal(sum(fit$rep$FishAgeComps_nLL[, y, seas, , 1]),
                 ref_miss0(obs, pred, ln_theta, iss, ar1 = TRUE, trans_rho = trans_rho),
                 tolerance = 1e-8)
  } # end y loop
})


test_that("the correlation parameter is estimated for the ar1 form and fixed for the iid one", {

  sim <- seasonal_M_sim()
  base <- seasonal_M_input(sim, list(M_spec = "fix", Fixed_natmort = seasonal_M_fixed(seasonal_M_cfg$M)))
  sim_data <- base$data

  build <- function(like_string) {
    input_list <- Setup_Mod_FishIdx_and_Comps(
      base,
      ObsFishIdx = sim_data$ObsFishIdx,
      ObsFishIdx_SE = sim_data$ObsFishIdx_SE,
      UseFishIdx = sim_data$UseFishIdx,
      ObsFishAgeComps = sim_data$ObsFishAgeComps,
      UseFishAgeComps = sim_data$UseFishAgeComps,
      ISS_FishAgeComps = sim_data$ISS_FishAgeComps,
      ObsFishLenComps = sim_data$ObsFishLenComps,
      UseFishLenComps = sim_data$UseFishLenComps,
      ISS_FishLenComps = sim_data$ISS_FishLenComps,
      fish_idx_type = "none",
      FishAgeComps_LikeType = like_string,
      FishLenComps_LikeType = "none",
      FishAgeComps_Type = "agg_Year_1-terminal_Fleet_1",
      FishLenComps_Type = "none_Year_1-terminal_Fleet_1"
    )
    input_list
  }

  ar1 <- build("1d-Logistic-Normal-miss0")
  iid <- build("iid-Logistic-Normal-miss0")

  expect_equal(ar1$data$FishAgeComps_LikeType, 6)
  expect_equal(iid$data$FishAgeComps_LikeType, 5)
  expect_true(any(!is.na(ar1$map$FishAge_corr_pars_agg)))
  expect_true(all(is.na(iid$map$FishAge_corr_pars_agg)))

  # the dispersion is estimated either way, the same as every non-multinomial form
  expect_true(any(!is.na(ar1$map$ln_FishAge_theta_agg)))
  expect_true(any(!is.na(iid$map$ln_FishAge_theta_agg)))
})


test_that("the internal composition residuals are refused rather than computed wrongly", {

  input_list <- list(data = list(do_internal_comp_osa = TRUE))
  expect_error(SPoRC:::check_miss0_osa(input_list, c(0, 5), "FishAgeComps_LikeType"),
               "zeros dropped")
  expect_error(SPoRC:::check_miss0_osa(input_list, c(6), "SrvAgeComps_LikeType"),
               "fixed length vector")

  # every other likelihood is unaffected, and so is a model not asking for the residuals
  expect_null(SPoRC:::check_miss0_osa(input_list, c(0, 1, 2, 3, 4), "FishAgeComps_LikeType"))
  expect_null(SPoRC:::check_miss0_osa(list(data = list(do_internal_comp_osa = FALSE)), c(5, 6), "x"))
})


test_that("the operating model draws on the same terms the estimation model fits on", {

  set.seed(99)
  expected <- c(0.05, 0.12, 0.2, 0.25, 0.18, 0.12, 0.08)

  for(comp_like in c(5, 6)) {
    draws <- replicate(400, SPoRC:::rlogistnormal(exp = expected,
                                                  pars = c(exp(-0.5), stats::plogis(0.6)),
                                                  comp_like = comp_like,
                                                  n_sexes = 1,
                                                  ISS = 100))
    expect_equal(dim(draws), c(length(expected), 400))
    expect_true(all(draws > 0))                       # a logistic normal draw has no zeros
    expect_equal(as.numeric(colSums(draws)), rep(1, 400), tolerance = 1e-10)
  } # end comp_like loop

  # a larger sample size draws closer to the expectation
  tight <- replicate(400, SPoRC:::rlogistnormal(expected, c(exp(-0.5), 0.6), 5, 1, ISS = 1000))
  loose <- replicate(400, SPoRC:::rlogistnormal(expected, c(exp(-0.5), 0.6), 5, 1, ISS = 25))
  expect_lt(mean(apply(tight, 2, function(x) sum(abs(x - expected)))),
            mean(apply(loose, 2, function(x) sum(abs(x - expected)))))

  # the correlation is put on the positive scale for these forms and the full range for the others
  expect_equal(SPoRC:::comp_corr_natural(0.6, 5), stats::plogis(0.6))
  expect_equal(SPoRC:::comp_corr_natural(0.6, 3), SPoRC:::rho_trans(0.6))
})


# The separable bin by sex form -----------------------------------------------

ref_miss0_2d <- function(obs, pred, ln_sigma, ISS, trans_rho_bin, trans_rho_sex, n_bins, n_sexes) {

  rho_bin <- stats::plogis(trans_rho_bin)
  rho_sex <- stats::plogis(trans_rho_sex)
  sex_corr <- matrix(rho_sex, n_sexes, n_sexes)
  diag(sex_corr) <- 1
  bin_corr <- outer(1:n_bins, 1:n_bins, function(i, j) rho_bin^abs(i - j))
  full <- kronecker(sex_corr, bin_corr)

  pos <- which(obs > 1e-15)
  n_pos <- length(pos)
  if(n_pos < 2) return(0)

  obs_pos <- obs[pos] / sum(obs[pos])
  pred_pos <- pred[pos] / sum(pred[pos])
  y <- log(obs_pos[-n_pos]) - log(obs_pos[n_pos])
  mu <- log(pred_pos[-n_pos]) - log(pred_pos[n_pos])

  keep <- pos[seq_len(n_pos - 1)]
  Sigma <- matrix(full[keep, keep], n_pos - 1, n_pos - 1) * exp(2 * ln_sigma) / ISS

  chol_S <- chol(Sigma)
  z <- backsolve(chol_S, y - mu, transpose = TRUE)
  -1 * (-0.5 * sum(z^2) - sum(log(diag(chol_S))) - 0.5 * (n_pos - 1) * log(2 * pi) - sum(log(obs_pos)))
}


test_that("the separable bin by sex form matches the reference", {

  n_bins <- 6; n_sexes <- 2
  ln_theta <- -0.35; rho_bin <- 0.7; rho_sex <- 0.4

  for(seed in 51:55) {
    set.seed(seed)
    prob <- stats::runif(n_bins * n_sexes, 0.02, 1)
    prob <- prob / sum(prob)
    obs <- as.numeric(stats::rmultinom(1, 45, prob))
    obs <- obs / sum(obs)

    Obs <- array(0, dim = c(1, n_bins, n_sexes)); Obs[1, , ] <- matrix(obs, n_bins, n_sexes)
    Exp <- array(0, dim = c(1, n_bins, n_sexes)); Exp[1, , ] <- matrix(prob, n_bins, n_sexes)
    corr <- array(0, dim = c(1, n_sexes, 3)); corr[1, 1, 1] <- rho_bin; corr[1, 1, 2] <- rho_sex

    got <- SPoRC:::Get_Comp_Likelihoods(
      Exp = Exp, Obs = Obs, ISS = array(45, c(1, n_sexes)),
      Wt_Mltnml = array(1, c(1, n_sexes)), ln_theta_agg = ln_theta,
      ln_theta = matrix(ln_theta, 1, n_sexes), LN_corr_pars = corr, LN_corr_pars_agg = 0,
      Comp_Type = 2, Likelihood_Type = 7, n_regions = 1, n_model_bins = n_bins,
      n_obs_bins = n_bins, n_sexes = n_sexes, age_or_len = 0, AgeingError = diag(n_bins),
      use = 1, addtocomp = 0)

    expect_equal(got[1, 1],
                 ref_miss0_2d(obs, prob, ln_theta, 45, rho_bin, rho_sex, n_bins, n_sexes),
                 tolerance = 1e-9)
  } # end seed loop
})


test_that("a zero sex correlation collapses the separable form onto the bin only form", {

  n_bins <- 7; n_sexes <- 2
  set.seed(61)
  prob <- stats::runif(n_bins * n_sexes, 0.05, 1)
  prob <- prob / sum(prob)
  obs <- as.numeric(stats::rmultinom(1, 300, prob))
  obs <- obs / sum(obs)

  # plogis(-30) is zero to machine precision, so the sexes are uncorrelated
  a <- ref_miss0_2d(obs, prob, -0.2, 300, trans_rho_bin = 0.9, trans_rho_sex = -30, n_bins, n_sexes)

  # the bin only form spaced along the stack differs, since it correlates across the sex boundary
  b <- SPoRC:::get_logistnormal_miss0_nLL(obs, prob, -0.2, 300, corr_type = 1, trans_rho = 0.9)
  expect_false(isTRUE(all.equal(a, b)))

  # but a single sex makes them the same statement
  one_sex <- ref_miss0_2d(obs[1:n_bins] / sum(obs[1:n_bins]), prob[1:n_bins] / sum(prob[1:n_bins]),
                          -0.2, 300, trans_rho_bin = 0.9, trans_rho_sex = -30, n_bins, n_sexes = 1)
  direct <- SPoRC:::get_logistnormal_miss0_nLL(obs[1:n_bins] / sum(obs[1:n_bins]),
                                               prob[1:n_bins] / sum(prob[1:n_bins]),
                                               -0.2, 300, corr_type = 1, trans_rho = 0.9)
  expect_equal(one_sex, direct, tolerance = 1e-9)
})


test_that("the separable form needs a composition joint across sexes", {

  sim <- seasonal_M_sim()
  base <- seasonal_M_input(sim, list(M_spec = "fix", Fixed_natmort = seasonal_M_fixed(seasonal_M_cfg$M)))
  sim_data <- base$data

  expect_error(
    Setup_Mod_FishIdx_and_Comps(
      base,
      ObsFishIdx = sim_data$ObsFishIdx, ObsFishIdx_SE = sim_data$ObsFishIdx_SE,
      UseFishIdx = sim_data$UseFishIdx,
      ObsFishAgeComps = sim_data$ObsFishAgeComps, UseFishAgeComps = sim_data$UseFishAgeComps,
      ISS_FishAgeComps = sim_data$ISS_FishAgeComps,
      ObsFishLenComps = sim_data$ObsFishLenComps, UseFishLenComps = sim_data$UseFishLenComps,
      ISS_FishLenComps = sim_data$ISS_FishLenComps,
      fish_idx_type = "none",
      FishAgeComps_LikeType = "2d-Logistic-Normal-miss0",
      FishLenComps_LikeType = "none",
      FishAgeComps_Type = "agg_Year_1-terminal_Fleet_1",
      FishLenComps_Type = "none_Year_1-terminal_Fleet_1"),
    "joint across sexes")
})


test_that("the aggregated iid logistic normal reads the same standard deviation as every other route", {

  n_bins <- 6
  obs <- c(0.05, 0.2, 0.3, 0.25, 0.15, 0.05)
  pred <- c(0.08, 0.18, 0.32, 0.22, 0.14, 0.06)
  ln_theta <- -0.4

  Obs <- array(0, dim = c(1, n_bins, 1)); Obs[1, , 1] <- obs
  Exp <- array(0, dim = c(1, n_bins, 1)); Exp[1, , 1] <- pred

  got <- SPoRC:::Get_Comp_Likelihoods(
    Exp = Exp, Obs = Obs, ISS = array(100, c(1, 1)), Wt_Mltnml = array(1, c(1, 1)),
    ln_theta_agg = ln_theta, ln_theta = matrix(ln_theta, 1, 1),
    LN_corr_pars = array(0, c(1, 1, 3)), LN_corr_pars_agg = 0,
    Comp_Type = 0, Likelihood_Type = 2, n_regions = 1, n_model_bins = n_bins,
    n_obs_bins = n_bins, n_sexes = 1, age_or_len = 0, AgeingError = diag(n_bins),
    use = 1, addtocomp = 0)

  # the variance is the square of the parameter, not its reciprocal
  expect_equal(got[1, 1],
               -SPoRC:::dlogistnormal(obs, pred, diag(rep(exp(ln_theta)^2, n_bins - 1)), TRUE),
               tolerance = 1e-10)

  # and the split by region and sex route says the same thing for one region and one sex
  split_got <- SPoRC:::Get_Comp_Likelihoods(
    Exp = Exp, Obs = Obs, ISS = array(100, c(1, 1)), Wt_Mltnml = array(1, c(1, 1)),
    ln_theta_agg = ln_theta, ln_theta = matrix(ln_theta, 1, 1),
    LN_corr_pars = array(0, c(1, 1, 3)), LN_corr_pars_agg = 0,
    Comp_Type = 1, Likelihood_Type = 2, n_regions = 1, n_model_bins = n_bins,
    n_obs_bins = n_bins, n_sexes = 1, age_or_len = 0, AgeingError = diag(n_bins),
    use = 1, addtocomp = 0)

  expect_equal(got[1, 1], split_got[1, 1], tolerance = 1e-10)
})
