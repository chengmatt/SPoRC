library(SPoRC)
library(testthat)

test_that("get_index_nLL dispatches the three index error structures", {

  obs <- c(10, 12, 9, 15)
  pred <- c(11, 11.5, 9.5, 14)
  sigma <- c(0.2, 0.25, 0.2, 0.3)

  test_that("lognormal matches a log-scale dnorm and honors the additive constant", {
    expect_equal(SPoRC:::get_index_nLL(obs, pred, sigma, 0),
                 -dnorm(log(obs), log(pred), sigma, TRUE), tolerance = 1e-12)
    expect_equal(SPoRC:::get_index_nLL(obs, pred, sigma, 0, const = 0.01),
                 -dnorm(log(obs + 0.01), log(pred + 0.01), sigma, TRUE), tolerance = 1e-12)
  })

  test_that("normal is on the arithmetic scale and ignores the constant", {
    expect_equal(SPoRC:::get_index_nLL(obs, pred, sigma, 1),
                 -dnorm(obs, pred, sigma, TRUE), tolerance = 1e-12)
    expect_equal(SPoRC:::get_index_nLL(obs, pred, sigma, 1, const = 0.01),
                 SPoRC:::get_index_nLL(obs, pred, sigma, 1), tolerance = 1e-12)
  })

  test_that("a diagonal covariance reduces the multivariate normal to independent normals", {
    Sigma <- diag(sigma^2)
    nLL <- SPoRC:::get_index_nLL(obs, pred, sigma, 2, Sigma)
    expect_equal(sum(nLL), sum(-dnorm(obs, pred, sigma, TRUE)), tolerance = 1e-10)
  })

  test_that("the multivariate normal total lands in the first element", {
    Sigma <- diag(sigma^2)
    nLL <- SPoRC:::get_index_nLL(obs, pred, sigma, 2, Sigma)
    expect_length(nLL, length(obs))
    expect_equal(nLL[-1], rep(0, length(obs) - 1))
  })

  test_that("correlation changes the multivariate normal answer", {
    Sigma <- diag(sigma^2)
    Sigma_corr <- Sigma
    Sigma_corr[1,2] <- Sigma_corr[2,1] <- 0.5 * sigma[1] * sigma[2]
    indep <- SPoRC:::get_index_nLL(obs, pred, sigma, 2, Sigma)[1]
    corr <- SPoRC:::get_index_nLL(obs, pred, sigma, 2, Sigma_corr)[1]
    expect_false(isTRUE(all.equal(indep, corr)))
  })

  test_that("the multivariate normal refuses to run without a covariance", {
    expect_error(SPoRC:::get_index_nLL(obs, pred, sigma, 2))
  })

})

test_that("parse_bin_subset builds the per-fleet bin selection array", {

  test_that("NULL selects every age for every fleet", {
    expect_equal(SPoRC:::parse_bin_subset(NULL, 5, 2, "x"), array(1, dim = c(5, 2)))
  })

  test_that("a list of age vectors selects per fleet, NULL meaning all ages", {
    out <- SPoRC:::parse_bin_subset(list(NULL, 1, c(2, 4)), 5, 3, "x")
    expect_equal(out[,1], rep(1, 5))
    expect_equal(out[,2], c(1, 0, 0, 0, 0))
    expect_equal(out[,3], c(0, 1, 0, 1, 0))
  })

  test_that("an array is accepted directly", {
    arr <- cbind(c(1,1,0), c(0,1,1))
    expect_equal(SPoRC:::parse_bin_subset(arr, 3, 2, "x"), array(as.numeric(arr), dim = c(3, 2)))
  })

  test_that("errors on out-of-range ages, wrong length, non-binary values, and empty fleets", {
    expect_error(SPoRC:::parse_bin_subset(list(9), 5, 1, "x"))
    expect_error(SPoRC:::parse_bin_subset(list(1, 2), 5, 3, "x"))
    expect_error(SPoRC:::parse_bin_subset(array(0.5, dim = c(3, 1)), 3, 1, "x"))
    expect_error(SPoRC:::parse_bin_subset(array(0, dim = c(3, 1)), 3, 1, "x"))
  })

})

test_that("parse_idx_cov validates only the fleets that need a covariance", {

  use_arr <- array(0, dim = c(1, 4, 1, 2))
  use_arr[1, 1:3, 1, 1] <- 1
  use_arr[1, 1:4, 1, 2] <- 1
  Sigma <- diag(c(1, 2, 3))

  test_that("the covariance is returned as given and a non-mvn fleet needs no matrix", {
    out <- SPoRC:::parse_idx_cov(list(Sigma, NULL), c(2, 0), use_arr, 2, "x")
    expect_equal(out[[1]], Sigma, tolerance = 1e-12)
    expect_null(out[[2]])
  })

  test_that("errors when the matrix is missing, the wrong size, asymmetric, or not positive definite", {
    expect_error(SPoRC:::parse_idx_cov(NULL, c(2, 0), use_arr, 2, "x"))
    expect_error(SPoRC:::parse_idx_cov(list(NULL, NULL), c(2, 0), use_arr, 2, "x"))
    expect_error(SPoRC:::parse_idx_cov(list(diag(4), NULL), c(2, 0), use_arr, 2, "x")) # fleet 1 fits 3 observations
    # RTMB::dmvnorm reads only the lower triangle, so an asymmetric matrix must not slip through
    asym <- Sigma; asym[1,2] <- 0.5
    expect_error(SPoRC:::parse_idx_cov(list(asym, NULL), c(2, 0), use_arr, 2, "x"))
    expect_error(SPoRC:::parse_idx_cov(list(matrix(c(1, 2, 2, 1), 2), NULL), c(2, 0), use_arr, 2, "x"))
  })

})

test_that("the multivariate normal index likelihood agrees with a hand-computed full density", {

  set.seed(31)
  n <- 6
  A <- matrix(rnorm(n * n), n)
  Sigma <- crossprod(A) + diag(n)
  obs <- rnorm(n, 10, 2)
  pred <- rnorm(n, 10, 2)

  test_that("RTMB's dmvnorm carries the 2*pi constant, so no offset is needed", {
    resid <- obs - pred
    full <- -0.5 * (n * log(2 * pi) + determinant(Sigma, logarithm = TRUE)$modulus[1] +
                      sum(resid * solve(Sigma, resid)))
    expect_equal(SPoRC:::get_index_nLL(obs, pred, rep(1, n), 2, Sigma)[1], -full, tolerance = 1e-10)
  })

  test_that("it is differentiable through an RTMB tape with a constant covariance", {
    f <- function(pars) {
      "c" <- RTMB::ADoverload("c")
      "[<-" <- RTMB::ADoverload("[<-")
      sum(SPoRC:::get_index_nLL(obs, pred + pars[1], rep(1, n), 2, Sigma))
    }
    obj <- RTMB::MakeADFun(f, c(0), silent = TRUE)
    expect_no_error(obj$fn(obj$par))
    expect_no_error(obj$gr(obj$par))
    fd <- (obj$fn(1e-5) - obj$fn(-1e-5)) / 2e-5
    expect_equal(as.numeric(obj$gr(0)), fd, tolerance = 1e-6)
  })

})

test_that("the survey observation model honors age subsets and analytic catchability", {

  n_pop <- 1; n_regions <- 1; n_yrs <- 3; n_seas <- 1; n_srv <- 2; n_sexes <- 1; n_ages <- 4

  NAA <- array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes))
  NAA[1,1,,1,,1] <- matrix(c(100, 80, 60, 40,
                             120, 90, 70, 50,
                             110, 85, 65, 45), nrow = n_yrs, byrow = TRUE)
  ZAA <- array(0.3, dim = dim(NAA))
  srv_sel <- array(1, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_srv))
  WAA_srv <- array(2, dim = dim(srv_sel))
  SrvIAA <- array(0, dim = dim(srv_sel))
  PredSrvIdx <- array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_srv))

  run <- function(srv_idx_ages = NULL, srv_q_type = NULL, ObsSrvIdx = NULL, UseSrvIdx = NULL, ln_q = 0) {
    get_survey_observation_model(
      n_pop = n_pop, n_regions = n_regions, n_yrs = n_yrs, n_seas = n_seas,
      n_srv_fleets = n_srv, n_sexes = n_sexes,
      srv_q_blocks = array(1, dim = c(n_regions, n_yrs, n_srv)),
      ln_srv_q = array(ln_q, dim = c(n_regions, 1, n_srv)),
      srv_q = array(0, dim = c(n_regions, n_yrs, n_srv)),
      do_srv_q_cov = 0, srv_q_cov = NULL, srv_q_coeff = NULL,
      srv_selex_type = 0, srv_sel = srv_sel, srv_sel_l = NULL, SizeAgeTrans = NULL,
      NAA = NAA, ZAA = ZAA, t_srv = array(0, dim = c(n_regions, n_seas, n_srv)),
      SrvIAA = SrvIAA, fit_lengths = 0, SrvIAL = NULL,
      srv_idx_type = c(0, 1), WAA_srv = WAA_srv, PredSrvIdx = PredSrvIdx,
      srv_idx_ages = srv_idx_ages, srv_q_type = srv_q_type,
      ObsSrvIdx = ObsSrvIdx, UseSrvIdx = UseSrvIdx
    )
  }

  test_that("with no age restriction the index sums every age", {
    out <- run()
    expect_equal(out$PredSrvIdx[1,1,,1,1], rowSums(NAA[1,1,,1,,1]))
    expect_equal(out$PredSrvIdx[1,1,,1,2], 2 * rowSums(NAA[1,1,,1,,1]))
  })

  test_that("restricting a fleet to one age makes it an index of that age alone", {
    ages <- array(1, dim = c(n_ages, n_srv))
    ages[,1] <- c(1, 0, 0, 0)
    out <- run(srv_idx_ages = ages)
    expect_equal(out$PredSrvIdx[1,1,,1,1], NAA[1,1,,1,1,1])
    # the unrestricted fleet is untouched
    expect_equal(out$PredSrvIdx[1,1,,1,2], 2 * rowSums(NAA[1,1,,1,,1]))
  })

  test_that("restricting to a contiguous run sums only those ages", {
    ages <- array(1, dim = c(n_ages, n_srv))
    ages[,1] <- c(0, 1, 1, 0)
    out <- run(srv_idx_ages = ages)
    expect_equal(out$PredSrvIdx[1,1,,1,1], rowSums(NAA[1,1,,1,2:3,1]))
  })

  test_that("the arithmetic solve makes the predicted and observed means agree", {
    Use <- array(0, dim = c(n_regions, n_yrs, n_seas, n_srv))
    Use[1,,1,1] <- 1
    Obs <- array(NA_real_, dim = dim(Use))
    Obs[1,,1,1] <- c(500, 700, 600)
    out <- run(srv_q_type = c(1, 0), ObsSrvIdx = Obs, UseSrvIdx = Use)
    expect_equal(mean(out$PredSrvIdx[1,1,,1,1]), mean(Obs[1,,1,1]), tolerance = 1e-10)
    expect_equal(as.numeric(unique(out$srv_q[1,,1])), mean(Obs[1,,1,1]) / mean(rowSums(NAA[1,1,,1,,1])), tolerance = 1e-10)
  })

  test_that("the log-scale solve makes the mean log residual zero", {
    Use <- array(0, dim = c(n_regions, n_yrs, n_seas, n_srv))
    Use[1,,1,1] <- 1
    Obs <- array(NA_real_, dim = dim(Use))
    Obs[1,,1,1] <- c(500, 700, 600)
    out <- run(srv_q_type = c(2, 0), ObsSrvIdx = Obs, UseSrvIdx = Use)
    expect_equal(mean(log(Obs[1,,1,1]) - log(out$PredSrvIdx[1,1,,1,1])), 0, tolerance = 1e-10)
  })

  test_that("an analytic solve ignores ln_srv_q, an estimated one does not", {
    Use <- array(0, dim = c(n_regions, n_yrs, n_seas, n_srv))
    Use[1,,1,1] <- 1
    Obs <- array(NA_real_, dim = dim(Use))
    Obs[1,,1,1] <- c(500, 700, 600)
    analytic_a <- run(srv_q_type = c(1, 0), ObsSrvIdx = Obs, UseSrvIdx = Use, ln_q = 0)
    analytic_b <- run(srv_q_type = c(1, 0), ObsSrvIdx = Obs, UseSrvIdx = Use, ln_q = 2)
    expect_equal(analytic_a$PredSrvIdx[1,1,,1,1], analytic_b$PredSrvIdx[1,1,,1,1], tolerance = 1e-10)
    expect_false(isTRUE(all.equal(analytic_a$PredSrvIdx[1,1,,1,2], analytic_b$PredSrvIdx[1,1,,1,2])))
  })

  test_that("only the years with observations enter the solve", {
    Use <- array(0, dim = c(n_regions, n_yrs, n_seas, n_srv))
    Use[1,1:2,1,1] <- 1
    Obs <- array(NA_real_, dim = dim(Use))
    Obs[1,1:2,1,1] <- c(500, 700)
    out <- run(srv_q_type = c(1, 0), ObsSrvIdx = Obs, UseSrvIdx = Use)
    expect_equal(mean(out$PredSrvIdx[1,1,1:2,1,1]), mean(Obs[1,1:2,1,1]), tolerance = 1e-10)
  })

})
