library(SPoRC)
library(testthat)

test_that("Custom Random Samplers produce consistent results", {

  set.seed(42)

  # --- rlogistnormal ---
  test_that("rlogistnormal returns a vector summing to 1", {
    exp_comp <- c(0.2, 0.3, 0.5)
    result <- rlogistnormal(exp_comp, pars = c(0.2), comp_like = 2, n_sexes = 1)
    expect_length(result, 3)
    expect_equal(sum(result), 1, tolerance = 1e-10)
    expect_true(all(result > 0))
  })

  test_that("rlogistnormal iid (comp_like=2) mean converges to expected proportions", {
    exp_comp <- c(0.2, 0.3, 0.5)
    draws <- replicate(5000, rlogistnormal(exp_comp, pars = c(0.05), comp_like = 2, n_sexes = 1))
    # draws is 3 x 5000; row means should be close to exp_comp with tight sigma
    row_means <- rowMeans(draws)
    expect_equal(row_means, exp_comp, tolerance = 0.02)
  })

  test_that("rlogistnormal AR1 by bin (comp_like=3) returns valid composition", {
    exp_comp <- c(0.1, 0.2, 0.3, 0.4)
    result <- rlogistnormal(exp_comp, pars = c(0.2, 0.5), comp_like = 3, n_sexes = 1)
    expect_length(result, 4)
    expect_equal(sum(result), 1, tolerance = 1e-10)
    expect_true(all(result > 0))
  })

  test_that("rlogistnormal Kronecker (comp_like=4) returns valid composition", {
    # 2 sexes x 3 bins = 6 categories
    exp_comp <- rep(1/6, 6)
    result <- rlogistnormal(exp_comp, pars = c(0.2, 0.3, 0.4), comp_like = 4, n_sexes = 2)
    expect_length(result, 6)
    expect_equal(sum(result), 1, tolerance = 1e-10)
    expect_true(all(result > 0))
  })

  test_that("rlogistnormal variance increases with sigma", {
    exp_comp <- c(0.3, 0.3, 0.4)
    draws_tight <- replicate(2000, rlogistnormal(exp_comp, pars = c(0.05), comp_like = 2, n_sexes = 1))
    draws_wide  <- replicate(2000, rlogistnormal(exp_comp, pars = c(1.0),  comp_like = 2, n_sexes = 1))
    # variance of first category should be larger under wide sigma
    expect_gt(var(draws_wide[1, ]), var(draws_tight[1, ]))
  })

  # --- rdirM ---
  test_that("rdirM returns a matrix of correct dimensions", {
    alpha <- c(2, 5, 3)
    result <- rdirM(n = 10, N = 100, alpha = alpha)
    expect_equal(dim(result), c(3, 10))
  })

  test_that("rdirM columns each sum to N", {
    alpha <- c(2, 5, 3)
    result <- rdirM(n = 20, N = 50, alpha = alpha)
    expect_true(all(colSums(result) == 50))
  })

  test_that("rdirM all counts are non-negative integers", {
    result <- rdirM(n = 10, N = 100, alpha = c(3, 3, 4))
    expect_true(all(result >= 0))
    expect_true(all(result == floor(result)))
  })

  test_that("rdirM row means converge to alpha / sum(alpha)", {
    alpha <- c(2, 5, 3)
    expected_props <- alpha / sum(alpha)
    result <- rdirM(n = 5000, N = 200, alpha = alpha)
    observed_props <- rowMeans(result) / 200
    expect_equal(observed_props, expected_props, tolerance = 0.02)
  })

  test_that("rdirM overdispersion exceeds multinomial variance", {
    # DM variance > multinomial variance for the same expected proportions
    alpha <- c(1, 1, 1)   # high overdispersion (small alpha)
    N <- 100
    result_dm <- rdirM(n = 3000, N = N, alpha = alpha)
    # multinomial variance for p=1/3: N * p * (1-p) = 100 * 1/3 * 2/3 ~ 22.2
    multinom_var <- N * (1/3) * (2/3)
    dm_var <- var(result_dm[1, ])
    expect_gt(dm_var, multinom_var)
  })

  # --- rinvgauss_rec ---
  test_that("rinvgauss_rec returns a vector of the correct length", {
    rec <- c(100, 200, 150, 180, 130)
    result <- rinvgauss_rec(sims = 50, recruitment = rec)
    expect_length(result, 50)
  })

  test_that("rinvgauss_rec returns strictly positive values", {
    rec <- c(100, 200, 150, 180, 130)
    result <- rinvgauss_rec(sims = 500, recruitment = rec)
    expect_true(all(result > 0))
  })

  test_that("rinvgauss_rec mean converges to arithmetic mean of recruitment", {
    rec <- c(100, 200, 150, 180, 130, 160, 140)
    result <- rinvgauss_rec(sims = 10000, recruitment = rec)
    expect_equal(mean(result), mean(rec), tolerance = mean(rec) * 0.05)
  })

  test_that("rinvgauss_rec is right-skewed (mean > median)", {
    # inverse-Gaussian is right-skewed; mean should exceed median
    rec <- c(50, 300, 100, 400, 80, 250)
    result <- rinvgauss_rec(sims = 5000, recruitment = rec)
    expect_gt(mean(result), median(result))
  })

  test_that("rinvgauss_rec variance scales with input recruitment variability", {
    low_var_rec  <- rep(100, 10)  + rnorm(10, 0, 5)
    high_var_rec <- c(50, 500, 80, 400, 60, 300, 90, 450, 70, 350)
    # avoid degenerate case where all recruitment identical
    low_var_rec <- pmax(low_var_rec, 10)
    result_low  <- rinvgauss_rec(sims = 2000, recruitment = low_var_rec)
    result_high <- rinvgauss_rec(sims = 2000, recruitment = high_var_rec)
    expect_gt(var(result_high), var(result_low))
  })

})
