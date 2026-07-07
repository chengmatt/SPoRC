library(SPoRC)
library(testthat)

test_that("resolve_sel_pen_wts works", {

  zeros_all_terms <- c(smooth_bin_curve = 0, smooth_bin_diff = 0, smooth_yr_diff = 0, smooth_yr_curve = 0, smooth_dome = 0, smooth_mean_center = 0)

  test_that("NULL pen_wts defaults every term to 0", {
    expect_equal(resolve_sel_pen_wts(NULL), zeros_all_terms)
  })

  test_that("explicit pen_wts overrides named terms and zeroes unnamed ones", {
    out <- resolve_sel_pen_wts(list(smooth_bin_curve = 5))
    expect_equal(out, c(smooth_bin_curve = 5, smooth_bin_diff = 0, smooth_yr_diff = 0, smooth_yr_curve = 0, smooth_dome = 0, smooth_mean_center = 0))
  })

  test_that("explicit pen_wts can set multiple terms independently", {
    out <- resolve_sel_pen_wts(list(smooth_yr_diff = 1, smooth_dome = 30, smooth_mean_center = 10000))
    expect_equal(out, c(smooth_bin_curve = 0, smooth_bin_diff = 0,
                        smooth_yr_diff = 1, smooth_yr_curve = 0, smooth_dome = 30, smooth_mean_center = 10000))
  })

  test_that("explicit pen_wts can set the smooth_bin_diff term independently of the other smooth_* terms", {
    out <- resolve_sel_pen_wts(list(smooth_bin_diff = 7))
    expect_equal(out, c(smooth_bin_curve = 0, smooth_bin_diff = 7,
                        smooth_yr_diff = 0, smooth_yr_curve = 0, smooth_dome = 0, smooth_mean_center = 0))
  })

  test_that("errors on unrecognized term names", {
    expect_error(resolve_sel_pen_wts(list(bogus = 1)))
    expect_error(resolve_sel_pen_wts(list(bin_curve = 1))) # removed legacy name
  })

})

test_that("Get_Selex_Smoothness_Penalty works", {

  # A perfectly smooth (log-linear-in-age, constant-across-years) surface should have
  # exactly zero second-difference penalty in both directions, and no dome penalty
  # (strictly increasing across bins).
  n_yrs <- 4; n_bins <- 6
  log_sel_smooth <- outer(rep(1, n_yrs), seq(0.1, 1.5, length.out = n_bins))
  sel_vals_smooth <- array(exp(log_sel_smooth), dim = c(1, n_yrs, n_bins, 1, 1))

  test_that("zero penalty for a perfectly linear-in-bin, constant-in-year surface", {
    pen <- Get_Selex_Smoothness_Penalty(sel_vals_smooth, wt_bin_curve = 1, wt_yr_curve = 1, wt_dome = 1)
    expect_equal(pen, 0, tolerance = 1e-10)
  })

  test_that("weight of 0 disables a term even when the surface is not smooth", {
    set.seed(1)
    sel_vals_jagged <- array(exp(matrix(rnorm(n_yrs * n_bins), nrow = n_yrs)), dim = c(1, n_yrs, n_bins, 1, 1))
    pen_all_zero <- Get_Selex_Smoothness_Penalty(sel_vals_jagged, wt_bin_curve = 0, wt_yr_curve = 0, wt_dome = 0)
    expect_equal(pen_all_zero, 0)
  })

  test_that("bin_curve penalty scales linearly with its weight and is <= 0 (positive-scale convention, subtracted)", {
    set.seed(2)
    sel_vals_jagged <- array(exp(matrix(rnorm(n_yrs * n_bins), nrow = n_yrs)), dim = c(1, n_yrs, n_bins, 1, 1))
    pen_wt1 <- Get_Selex_Smoothness_Penalty(sel_vals_jagged, wt_bin_curve = 1)
    pen_wt2 <- Get_Selex_Smoothness_Penalty(sel_vals_jagged, wt_bin_curve = 2)
    expect_lt(pen_wt1, 0)
    expect_equal(pen_wt2, 2 * pen_wt1, tolerance = 1e-10)
  })

  test_that("yr_curve penalty only depends on the year dimension weight", {
    set.seed(3)
    sel_vals_jagged <- array(exp(matrix(rnorm(n_yrs * n_bins), nrow = n_yrs)), dim = c(1, n_yrs, n_bins, 1, 1))
    pen_age_only <- Get_Selex_Smoothness_Penalty(sel_vals_jagged, wt_bin_curve = 1, wt_yr_curve = 0)
    pen_yr_only  <- Get_Selex_Smoothness_Penalty(sel_vals_jagged, wt_bin_curve = 0, wt_yr_curve = 1)
    pen_both     <- Get_Selex_Smoothness_Penalty(sel_vals_jagged, wt_bin_curve = 1, wt_yr_curve = 1)
    expect_equal(pen_both, pen_age_only + pen_yr_only, tolerance = 1e-10)
  })

  test_that("dome penalty is zero for monotonically increasing selectivity, negative for decreasing", {
    increasing <- array(exp(matrix(rep(seq(0.1, 2, length.out = n_bins), n_yrs), nrow = n_yrs, byrow = TRUE)),
                        dim = c(1, n_yrs, n_bins, 1, 1))
    decreasing <- array(exp(matrix(rep(seq(2, 0.1, length.out = n_bins), n_yrs), nrow = n_yrs, byrow = TRUE)),
                        dim = c(1, n_yrs, n_bins, 1, 1))
    expect_equal(Get_Selex_Smoothness_Penalty(increasing, wt_dome = 1), 0)
    expect_lt(Get_Selex_Smoothness_Penalty(decreasing, wt_dome = 1), 0)
  })

  test_that("dome penalty is differentiable through RTMB (max-based hinge, not if())", {
    f <- function(pars) {
      "c" <- RTMB::ADoverload("c")
      "[<-" <- RTMB::ADoverload("[<-")
      sel_vals <- array(0, dim = c(1, 1, 3, 1, 1))
      sel_vals[1,1,1,1,1] <- exp(pars[1])
      sel_vals[1,1,2,1,1] <- exp(pars[2])
      sel_vals[1,1,3,1,1] <- exp(pars[3])
      -Get_Selex_Smoothness_Penalty(sel_vals, wt_dome = 1)
    }
    obj <- RTMB::MakeADFun(f, c(1, 0.5, 0.2), silent = TRUE) # decreasing across bins
    expect_no_error(obj$fn(obj$par))
    expect_no_error(obj$gr(obj$par))
  })

  test_that("bin_diff penalty is zero for a flat (constant-across-bins) selectivity curve", {
    flat <- array(exp(1.3), dim = c(1, n_yrs, n_bins, 1, 1))
    expect_equal(Get_Selex_Smoothness_Penalty(flat, wt_bin_diff = 1), 0, tolerance = 1e-10)
  })

  test_that("bin_diff penalty is unconditional (unlike wt_dome, an increasing curve also contributes)", {
    increasing <- array(exp(matrix(rep(seq(0.1, 2, length.out = n_bins), n_yrs), nrow = n_yrs, byrow = TRUE)),
                        dim = c(1, n_yrs, n_bins, 1, 1))
    decreasing <- array(exp(matrix(rep(seq(2, 0.1, length.out = n_bins), n_yrs), nrow = n_yrs, byrow = TRUE)),
                        dim = c(1, n_yrs, n_bins, 1, 1))
    # wt_dome is 0 for the increasing curve (only decreases penalized); wt_bin_diff penalizes both
    expect_equal(Get_Selex_Smoothness_Penalty(increasing, wt_dome = 1), 0)
    expect_lt(Get_Selex_Smoothness_Penalty(increasing, wt_bin_diff = 1), 0)
    # a monotonic curve with the same step size in both directions gives the same magnitude penalty either way
    expect_equal(Get_Selex_Smoothness_Penalty(increasing, wt_bin_diff = 1),
                Get_Selex_Smoothness_Penalty(decreasing, wt_bin_diff = 1), tolerance = 1e-10)
  })

  test_that("bin_diff penalty matches manual computation and is normalized by n_bins", {
    set.seed(6)
    sel_vals_jagged <- array(exp(matrix(rnorm(n_yrs * n_bins), nrow = n_yrs)), dim = c(1, n_yrs, n_bins, 1, 1))
    pen <- Get_Selex_Smoothness_Penalty(sel_vals_jagged, wt_bin_diff = 1)
    raw_ss <- sum(sapply(1:n_yrs, function(y) sum(diff(log(sel_vals_jagged[1,y,,1,1]))^2)))
    expect_equal(pen, -raw_ss / n_bins, tolerance = 1e-10)
  })

  test_that("bin_curve and yr_curve penalties are normalized by n_bins / n_yrs", {
    set.seed(4)
    sel_vals_jagged <- array(exp(matrix(rnorm(n_yrs * n_bins), nrow = n_yrs)), dim = c(1, n_yrs, n_bins, 1, 1))
    pen_bin <- Get_Selex_Smoothness_Penalty(sel_vals_jagged, wt_bin_curve = 1)
    pen_yr  <- Get_Selex_Smoothness_Penalty(sel_vals_jagged, wt_yr_curve = 1)

    # manually compute the unnormalized sum of squared second differences
    raw_bin_ss <- sum(sapply(1:n_yrs, function(y) sum(diff(diff(log(sel_vals_jagged[1,y,,1,1])))^2)))
    raw_yr_ss  <- sum(sapply(1:n_bins, function(b) sum(diff(diff(log(sel_vals_jagged[1,,b,1,1])))^2)))

    expect_equal(pen_bin, -raw_bin_ss / n_bins, tolerance = 1e-10)
    expect_equal(pen_yr, -raw_yr_ss / n_yrs, tolerance = 1e-10)
  })

  test_that("yr_diff (inter-annual first difference) is zero for a time-invariant surface, negative and normalized otherwise", {
    time_invariant <- array(exp(matrix(rep(seq(0.1, 1.5, length.out = n_bins), n_yrs), nrow = n_yrs, byrow = TRUE)),
                            dim = c(1, n_yrs, n_bins, 1, 1))
    expect_equal(Get_Selex_Smoothness_Penalty(time_invariant, wt_yr_diff = 1), 0)

    set.seed(5)
    sel_vals_jagged <- array(exp(matrix(rnorm(n_yrs * n_bins), nrow = n_yrs)), dim = c(1, n_yrs, n_bins, 1, 1))
    pen <- Get_Selex_Smoothness_Penalty(sel_vals_jagged, wt_yr_diff = 1)
    raw_ss <- sum(sapply(1:n_bins, function(b) sum(diff(log(sel_vals_jagged[1,,b,1,1]))^2)))
    expect_equal(pen, -raw_ss / n_yrs, tolerance = 1e-10)
  })

  test_that("mean_center penalty is zero when every year's log-selectivity averages to zero, negative otherwise", {
    n_bins_mc <- 3
    zero_mean <- array(0, dim = c(1, 2, n_bins_mc, 1, 1))
    zero_mean[1,1,,1,1] <- exp(c(-1, 0, 1))  # mean(log(.)) == 0
    zero_mean[1,2,,1,1] <- exp(c(-2, 0, 2))  # mean(log(.)) == 0
    expect_equal(Get_Selex_Smoothness_Penalty(zero_mean, wt_mean_center = 1), 0, tolerance = 1e-10)

    shifted <- array(0, dim = c(1, 2, n_bins_mc, 1, 1))
    shifted[1,1,,1,1] <- exp(c(-1, 0, 1) + 0.5)
    shifted[1,2,,1,1] <- exp(c(-2, 0, 2) + 0.5)
    pen <- Get_Selex_Smoothness_Penalty(shifted, wt_mean_center = 1)
    expect_equal(pen, -2 * 0.5^2, tolerance = 1e-10) # two years, each contributing -1*0.5^2 = -0.25
  })

  test_that("mean_center penalty scales with the ADMB default weight of 10000", {
    sel_vals <- array(0, dim = c(1, 1, 3, 1, 1))
    sel_vals[1,1,,1,1] <- exp(c(-1, 0, 1) + 0.1) # mean(log(.)) == 0.1
    pen <- Get_Selex_Smoothness_Penalty(sel_vals, wt_mean_center = 10000)
    expect_equal(pen, -10000 * 0.1^2, tolerance = 1e-8)
  })

})

test_that("Get_sel_PE_loglik modular penalty weights work", {

  test_that("bin_curve/yr_curve pen_wts are ignored by Get_sel_PE_loglik (caller's responsibility now)", {
    n_yrs <- 5; n_bins <- 4; n_sexes <- 1
    map_sel_devs <- array(0, dim = c(1, n_yrs, n_bins, n_sexes))
    map_sel_devs[1,,,1] <- matrix(rep(1:n_bins, each = n_yrs), nrow = n_yrs) # one unique dev per bin, shared across years
    set.seed(10)
    ln_devs <- array(rnorm(n_yrs * n_bins) * 0.1, dim = c(1, n_yrs, n_bins, n_sexes, 1))
    PE_pars <- array(0, dim = c(1, 4, n_sexes, 1))

    ll_zero <- Get_sel_PE_loglik(PE_model = 3, PE_pars = PE_pars, ln_devs = ln_devs, map_sel_devs = map_sel_devs,
                                 pen_wts = c(yr_devs = 0, bin_curve = 0, yr_curve = 0),
                                 min_sel_devs_shared_bins = 1:n_bins)
    ll_pen <- Get_sel_PE_loglik(PE_model = 3, PE_pars = PE_pars, ln_devs = ln_devs, map_sel_devs = map_sel_devs,
                                pen_wts = c(yr_devs = 1, bin_curve = 1, yr_curve = 1),
                                min_sel_devs_shared_bins = 1:n_bins)

    # PE_model = 3: yr_devs only applies to PE_model 1-2, and bin_curve/yr_curve are no
    # longer read here at all -> all three weights are inert, so results must be identical
    expect_equal(ll_pen, ll_zero, tolerance = 1e-10)
  })

  test_that("bin_curve and yr_curve are now applied by the caller via Get_Selex_Smoothness_Penalty, normalized", {
    n_yrs <- 5; n_bins <- 4; n_sexes <- 1
    set.seed(11)
    sel_vals <- array(exp(matrix(rnorm(n_yrs * n_bins), nrow = n_yrs)), dim = c(1, n_yrs, n_bins, n_sexes, 1))

    pen_base <- Get_Selex_Smoothness_Penalty(sel_vals, wt_bin_curve = 0, wt_yr_curve = 0, normalize = TRUE)
    pen_age  <- Get_Selex_Smoothness_Penalty(sel_vals, wt_bin_curve = 1, wt_yr_curve = 0, normalize = TRUE)
    pen_yr   <- Get_Selex_Smoothness_Penalty(sel_vals, wt_bin_curve = 0, wt_yr_curve = 1, normalize = TRUE)

    expect_true(pen_age != pen_base)
    expect_true(pen_yr != pen_base)
    expect_true(pen_age != pen_yr)
  })

  test_that("yr_devs weight only affects PE_model 1-2, not 3-5", {
    n_yrs <- 5; n_bins <- 4; n_sexes <- 1
    map_sel_devs <- array(0, dim = c(1, n_yrs, n_bins, n_sexes))
    map_sel_devs[1,,,1] <- matrix(rep(1:n_bins, each = n_yrs), nrow = n_yrs) # one unique dev per bin, shared across years
    set.seed(12)
    ln_devs <- array(rnorm(n_yrs * n_bins) * 0.1, dim = c(1, n_yrs, n_bins, n_sexes, 1))
    PE_pars <- array(0, dim = c(1, 4, n_sexes, 1))

    ll_no_yrdevs <- Get_sel_PE_loglik(1, PE_pars, ln_devs, map_sel_devs,
                                      pen_wts = c(yr_devs = 0, bin_curve = 0, yr_curve = 0), 1:n_bins)
    ll_yrdevs <- Get_sel_PE_loglik(1, PE_pars, ln_devs, map_sel_devs,
                                   pen_wts = c(yr_devs = 1, bin_curve = 0, yr_curve = 0), 1:n_bins)
    expect_true(ll_yrdevs != ll_no_yrdevs)
  })

})
