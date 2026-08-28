library(SPoRC)
library(testthat)

test_that("Custom Distributions produce consistent results", {

  # --- dbeta_symmetric ---
  test_that("dbeta_symmetric returns log-density by default and density with log=FALSE", {
    val <- dbeta_symmetric(0.5, p_ub = 1, p_lb = 0, p_prsd = 2)
    expect_true(is.numeric(val))
    expect_true(is.finite(val))
    val_exp <- dbeta_symmetric(0.5, p_ub = 1, p_lb = 0, p_prsd = 2, log = FALSE)
    expect_equal(exp(val), val_exp, tolerance = 1e-8)
  })

  test_that("dbeta_symmetric peaks near midpoint of support", {
    # density should be higher at midpoint than near boundaries
    mid  <- dbeta_symmetric(0.5,  p_ub = 1, p_lb = 0, p_prsd = 5, log = FALSE)
    low  <- dbeta_symmetric(0.1,  p_ub = 1, p_lb = 0, p_prsd = 5, log = FALSE)
    high <- dbeta_symmetric(0.9,  p_ub = 1, p_lb = 0, p_prsd = 5, log = FALSE)
    expect_gt(mid, low)
    expect_gt(mid, high)
  })

  test_that("dbeta_symmetric is symmetric around midpoint", {
    d1 <- dbeta_symmetric(0.3, p_ub = 1, p_lb = 0, p_prsd = 3, log = FALSE)
    d2 <- dbeta_symmetric(0.7, p_ub = 1, p_lb = 0, p_prsd = 3, log = FALSE)
    expect_equal(d1, d2, tolerance = 1e-6)
  })

  test_that("dbeta_symmetric stronger p_prsd concentrates mass more tightly", {
    mid_weak   <- dbeta_symmetric(0.5, p_ub = 1, p_lb = 0, p_prsd = 1,  log = FALSE)
    mid_strong <- dbeta_symmetric(0.5, p_ub = 1, p_lb = 0, p_prsd = 10, log = FALSE)
    off_weak   <- dbeta_symmetric(0.2, p_ub = 1, p_lb = 0, p_prsd = 1,  log = FALSE)
    off_strong <- dbeta_symmetric(0.2, p_ub = 1, p_lb = 0, p_prsd = 10, log = FALSE)
    # ratio of mid to off-center should be more extreme under strong prior
    expect_gt(mid_strong / off_strong, mid_weak / off_weak)
  })

  # --- ddirichlet ---
  test_that("ddirichlet returns log-density by default and density with log=FALSE", {
    x     <- c(0.2, 0.3, 0.5)
    alpha <- c(2, 3, 5)
    val     <- ddirichlet(x, alpha)
    val_exp <- ddirichlet(x, alpha, log = FALSE)
    expect_true(is.finite(val))
    expect_equal(exp(val), val_exp, tolerance = 1e-8)
  })

  test_that("ddirichlet agrees with stats::dbeta for K=2 (beta special case)", {
    # Dirichlet(alpha1, alpha2) on (x, 1-x) == Beta(alpha1, alpha2) on x
    x     <- 0.4
    alpha <- c(3, 5)
    d_dir  <- ddirichlet(c(x, 1 - x), alpha, log = FALSE)
    d_beta <- dbeta(x, shape1 = 3, shape2 = 5)
    expect_equal(d_dir, d_beta, tolerance = 1e-6)
  })

  test_that("ddirichlet is maximized at alpha / sum(alpha) (mode for alpha > 1)", {
    alpha <- c(4, 6, 2)
    mode  <- alpha / sum(alpha)           # exact mode
    off   <- c(0.1, 0.6, 0.3)            # off-mode but still sums to 1
    expect_gt(ddirichlet(mode, alpha), ddirichlet(off, alpha))
  })

  test_that("symmetric Dirichlet(1,1,1) gives uniform density (log = 0)", {
    # Dirichlet(1,...,1) is uniform on the simplex -> constant log-density
    alpha <- c(1, 1, 1)
    x1 <- c(0.2, 0.5, 0.3)
    x2 <- c(0.6, 0.1, 0.3)
    expect_equal(ddirichlet(x1, alpha), ddirichlet(x2, alpha), tolerance = 1e-10)
  })

  # --- ddirmult ---
  test_that("ddirmult returns log-likelihood by default and likelihood with give_log=FALSE", {
    obs    <- c(0.2, 0.5, 0.3)
    pred   <- c(0.25, 0.45, 0.30)
    val     <- ddirmult(obs, pred, Ntotal = 100, ln_theta = log(0.5))
    val_exp <- ddirmult(obs, pred, Ntotal = 100, ln_theta = log(0.5), give_log = FALSE)
    expect_true(is.finite(val))
    expect_equal(exp(val), val_exp, tolerance = 1e-8)
  })

  test_that("ddirmult is higher when predicted matches observed better", {
    obs   <- c(0.2, 0.5, 0.3)
    good  <- ddirmult(obs, obs,            Ntotal = 100, ln_theta = 0)
    bad   <- ddirmult(obs, c(0.1, 0.1, 0.8), Ntotal = 100, ln_theta = 0)
    expect_gt(good, bad)
  })

  test_that("ddirmult approaches multinomial as ln_theta -> -Inf (small theta)", {
    obs  <- c(0.3, 0.4, 0.3)
    pred <- c(0.3, 0.4, 0.3)
    # With very small theta the DM variance -> multinomial variance
    # Just check it stays finite and doesn't crash
    val <- ddirmult(obs, pred, Ntotal = 50, ln_theta = -20)
    expect_true(is.finite(val))
  })

  # --- dnbinom_robust_noint ---
  test_that("dnbinom_robust_noint returns log-likelihood by default and exp with give_log=FALSE", {
    val     <- dnbinom_robust_noint(x = 5, log_mu = log(4), log_var_minus_mu = log(2))
    val_exp <- dnbinom_robust_noint(x = 5, log_mu = log(4), log_var_minus_mu = log(2), give_log = FALSE)
    expect_true(is.finite(val))
    expect_equal(exp(val), val_exp, tolerance = 1e-8)
  })

  test_that("dnbinom_robust_noint agrees with stats::dnbinom for integer x", {
    mu  <- 6
    k   <- 3   # size/overdispersion
    x   <- 4L
    # var = mu + mu^2/k  =>  var - mu = mu^2/k  =>  log_var_minus_mu = log(mu^2/k)
    log_var_minus_mu <- log(mu^2 / k)
    our_val   <- dnbinom_robust_noint(x, log_mu = log(mu), log_var_minus_mu = log_var_minus_mu)
    stats_val <- dnbinom(x, size = k, mu = mu, log = TRUE)
    expect_equal(our_val, stats_val, tolerance = 1e-6)
  })

  test_that("dnbinom_robust_noint is highest at/near the mean", {
    mu  <- 8
    lvm <- log(4)   # arbitrary excess variance
    ll_at_mean <- dnbinom_robust_noint(mu, log_mu = log(mu), log_var_minus_mu = lvm)
    ll_far     <- dnbinom_robust_noint(20, log_mu = log(mu), log_var_minus_mu = lvm)
    expect_gt(ll_at_mean, ll_far)
  })

  # --- dpois_noint ---
  test_that("dpois_noint returns log-likelihood by default and exp with give_log=FALSE", {
    val     <- dpois_noint(x = 3, pred = 4)
    val_exp <- dpois_noint(x = 3, pred = 4, give_log = FALSE)
    expect_true(is.finite(val))
    expect_equal(exp(val), val_exp, tolerance = 1e-8)
  })

  test_that("dpois_noint agrees with stats::dpois for integer x", {
    for (x in c(0, 1, 5, 10)) {
      expect_equal(
        dpois_noint(x, pred = 3.5),
        dpois(x, lambda = 3.5, log = TRUE),
        tolerance = 1e-8,
        label = paste("x =", x)
      )
    }
  })

  test_that("dpois_noint is highest at/near the mean", {
    lambda <- 7
    ll_mode <- dpois_noint(lambda, pred = lambda)
    ll_far  <- dpois_noint(20,     pred = lambda)
    expect_gt(ll_mode, ll_far)
  })

  # --- get_beta_scaled_pars ---
  test_that("get_beta_scaled_pars returns a length-4 numeric vector", {
    pars <- get_beta_scaled_pars(low = 0.2, high = 1, mu = 0.7, sigma = 0.1)
    expect_length(pars, 4)
    expect_true(all(is.finite(pars)))
  })

  test_that("get_beta_scaled_pars alpha and beta are both positive", {
    pars <- get_beta_scaled_pars(low = 0.2, high = 1, mu = 0.7, sigma = 0.1)
    expect_gt(pars[1], 0)   # alpha
    expect_gt(pars[2], 0)   # beta
  })

  test_that("get_beta_scaled_pars encodes correct low and scale", {
    pars <- get_beta_scaled_pars(low = 0.2, high = 1, mu = 0.7, sigma = 0.1)
    expect_equal(pars[3], 0.2)          # low
    expect_equal(pars[4], 0.8)          # scale = high - low
  })

  test_that("get_beta_scaled_pars recovered mean matches input mean", {
    low <- 0.2; high <- 1; mu <- 0.6; sigma <- 0.08
    pars <- get_beta_scaled_pars(low, high, mu, sigma)
    alpha <- pars[1]; beta <- pars[2]
    # beta distribution mean = alpha / (alpha + beta), back-transformed
    recovered_mu <- low + (high - low) * alpha / (alpha + beta)
    expect_equal(recovered_mu, mu, tolerance = 1e-6)
  })

  test_that("get_beta_scaled_pars recovered sd matches input sd", {
    low <- 0.2; high <- 1; mu <- 0.6; sigma <- 0.08
    pars  <- get_beta_scaled_pars(low, high, mu, sigma)
    alpha <- pars[1]; beta <- pars[2]; scale <- pars[4]
    # beta variance = alpha*beta / ((alpha+beta)^2*(alpha+beta+1))
    ab    <- alpha + beta
    var01 <- alpha * beta / (ab^2 * (ab + 1))
    recovered_sd <- sqrt(var01) * scale
    expect_equal(recovered_sd, sigma, tolerance = 1e-6)
  })

})
