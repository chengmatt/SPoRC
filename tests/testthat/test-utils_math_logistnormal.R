library(testthat)
library(SPoRC)

test_that("AR(1) and logistic-normal covariance utilities produce correct output", {

  # get_AR1_CorrMat
  test_that("get_AR1_CorrMat returns correct structure and values", {
    mat <- SPoRC:::get_AR1_CorrMat(4, 0.5)
    expect_equal(dim(mat), c(4, 4))
    # diagonal should be 1 (rho^0)
    expect_equal(diag(mat), rep(1, 4))
    # off-diagonals: lag-1 = 0.5, lag-2 = 0.25, lag-3 = 0.125
    expect_equal(mat[1, 2], 0.5)
    expect_equal(mat[1, 3], 0.25)
    expect_equal(mat[1, 4], 0.125)
    # matrix should be symmetric
    expect_equal(mat, t(mat))
    # rho = 0 should return identity
    expect_equal(SPoRC:::get_AR1_CorrMat(3, 0), diag(3))
  })

  # get_Constant_CorrMat
  test_that("get_Constant_CorrMat returns correct structure and values", {
    mat <- SPoRC:::get_Constant_CorrMat(3, 0.6)
    expect_equal(dim(mat), c(3, 3))
    # diagonal should be 1
    expect_equal(diag(mat), rep(1, 3))
    # all off-diagonals should equal rho
    off_diag <- mat[lower.tri(mat)]
    expect_true(all(off_diag == 0.6))
    # symmetric
    expect_equal(mat, t(mat))
    # rho = 0 should return identity
    expect_equal(SPoRC:::get_Constant_CorrMat(3, 0), diag(3))
  })

  # rho_trans
  test_that("rho_trans maps to (-1, 1) correctly", {
    # output should be in (-1, 1)
    expect_true(rho_trans(0) > -1 && rho_trans(0) < 1)
    expect_true(rho_trans(10) < 1)
    expect_true(rho_trans(-10) > -1)
    # rho_trans(0) should equal 0
    expect_equal(rho_trans(0), 0)
    # should be monotonically increasing
    expect_true(rho_trans(1) > rho_trans(0))
    expect_true(rho_trans(-1) < rho_trans(0))
    # antisymmetric around 0
    expect_equal(rho_trans(1), -rho_trans(-1))
  })

  # get_logistN_Sigma
  test_that("get_logistN_Sigma comp_like = 2 (iid) returns scaled identity", {
    Sigma <- get_logistN_Sigma(comp_like = 2, n_bins = 5, n_sexes = NULL, theta = 0.5)
    expect_equal(dim(Sigma), c(5, 5))
    expect_equal(diag(Sigma), rep(0.25, 5))   # theta^2
    expect_equal(sum(Sigma[lower.tri(Sigma)]), 0)  # no off-diagonal
  })

  test_that("get_logistN_Sigma comp_like = 3 (AR1 bins) scales correctly", {
    Sigma <- get_logistN_Sigma(comp_like = 3, n_bins = 4, n_sexes = NULL,
                               theta = 2, corr_b = 0.5)
    expect_equal(dim(Sigma), c(4, 4))
    # diagonal should be theta^2
    expect_equal(diag(Sigma), rep(4, 4))
    # lag-1 off-diagonal should be theta^2 * rho
    expect_equal(Sigma[1, 2], 4 * 0.5)
    expect_equal(Sigma[1, 3], 4 * 0.25)
    expect_equal(Sigma, t(Sigma))
  })

  test_that("get_logistN_Sigma comp_like = 4 (AR1 bins x constant sexes) has correct dimensions and structure", {
    n_bins <- 3; n_sexes <- 2; theta <- 1; corr_b <- 0.4; corr_s <- 0.7
    Sigma <- get_logistN_Sigma(comp_like = 4, n_bins = n_bins, n_sexes = n_sexes,
                               theta = theta, corr_b = corr_b, corr_s = corr_s)
    expect_equal(dim(Sigma), c(n_bins * n_sexes, n_bins * n_sexes))
    expect_equal(Sigma, t(Sigma))   # symmetric
    # diagonal should all be theta^2
    expect_equal(diag(Sigma), rep(theta^2, n_bins * n_sexes))
    # corr_s = 0 and corr_b = 0 should give identity (scaled by theta^2)
    Sigma_indep <- get_logistN_Sigma(comp_like = 4, n_bins = 3, n_sexes = 2,
                                     theta = 1, corr_b = 0, corr_s = 0)
    expect_equal(Sigma_indep, diag(6))
  })

})
