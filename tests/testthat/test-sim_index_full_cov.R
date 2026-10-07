# A multivariate normal index is drawn over its whole covariance for the cells it covers, so the draws
# have the correlation the estimation model's dmvnorm reads, not the one-factor approximation of it.

library(SPoRC)
library(testthat)

test_that("the errors inside the covariance are drawn with all of it, an ar1 included", {

  # an ar1 across ten years, which one common factor gets badly wrong at short lags
  n_cells <- 10
  phi <- 0.9
  sd_cell <- 0.2
  S <- sd_cell^2 * phi^abs(outer(seq_len(n_cells), seq_len(n_cells), "-"))
  use <- array(1, dim = c(1, n_cells, 1, 1)) # region, year, season, fleet

  sim_env <- new.env()
  sim_env$n_sims <- 20000
  sim_env$srv_idx_mvn <- build_idx_factor(list(S), 2, use, 1, "SrvIdx_Cov")
  set.seed(8)
  draw_sim_idx_mvn(sim_env)
  eps <- sim_env$srv_idx_mvn_eps[[1]] # covariance row by replicate

  empirical <- stats::cor(t(eps))
  expect_equal(empirical[1, 2], phi, tolerance = 0.02) # lag one
  expect_equal(empirical[1, 10], phi^9, tolerance = 0.05) # lag nine
  expect_equal(apply(eps, 1, stats::sd), rep(sd_cell, n_cells), tolerance = 0.03)

  # the factor approximation it replaces, for comparison: 0.62 at lag one where the covariance says 0.9
  fac <- cov_to_factor(S)
  factor_lag1 <- fac$lambda[1] * fac$lambda[2]
  expect_gt(abs(factor_lag1 - phi), 0.2)
})

test_that("a cell inside the covariance resolves to its row and one outside to none", {

  use <- array(0, dim = c(1, 4, 1, 1))
  use[1, c(1, 2, 4), 1, 1] <- 1
  fac <- build_idx_factor(list(diag(c(1, 4, 9))), 2, use, 1, "SrvIdx_Cov")[[1]]

  expect_equal(resolve_idx_factor(fac, 1, 4, 1)$row, 3) # the third observed cell
  expect_true(is.na(resolve_idx_factor(fac, 1, 3, 1)$row)) # a year the fleet did not observe
  expect_true(is.na(resolve_idx_factor(fac, 1, 9, 1)$row)) # a projection year
  expect_equal(fac$chol_lower %*% t(fac$chol_lower), diag(c(1, 4, 9)))
})

test_that("a fleet without a covariance draws no joint errors", {

  sim_env <- new.env()
  sim_env$n_sims <- 3
  sim_env$fish_idx_mvn <- build_idx_factor(list(NULL, diag(2)), c(0, 2), array(1, dim = c(1, 2, 1, 2)), 2, "FishIdx_Cov")
  draw_sim_idx_mvn(sim_env)

  expect_null(sim_env$fish_idx_mvn_eps[[1]])
  expect_equal(dim(sim_env$fish_idx_mvn_eps[[2]]), c(2, 3))
  expect_length(sim_env$srv_idx_mvn_eps, 0) # no survey covariance at all
})
