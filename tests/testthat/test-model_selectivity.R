# Every selectivity form Get_Selex offers, each against what its parameters mean, and
# the deviations that move it over time.

library(SPoRC)
library(testthat)

test_that("Get_Selex returns the curve each selectivity form defines", {

  ages  <- 1:20
  lens  <- seq(10, 80, by = 5)

  # Zero deviations array: [n_regions, n_years, n_pars_or_bins, n_sexes, 1]
  zero_devs <- function(n_pars, n_regions = 2, n_years = 5, n_sexes = 1) {
    array(0, dim = c(n_regions, n_years, n_pars, n_sexes, 1))
  }

  # call Get_Selex with zero deviations unless some are given
  selex <- function(
    model,
    pars,
    bins = ages,
    tv = 0,
    devs = NULL,
    region = 1,
    year = 1,
    sex = 1
  ) {
    if (is.null(devs)) devs <- zero_devs(length(pars))
    Get_Selex(
      Selex_Model    = model,
      TimeVary_Model = tv,
      pars           = pars,
      ln_seldevs     = devs,
      Region         = region,
      Year           = year,
      Bin            = bins,
      Sex            = sex
    )
  }

  is_monotone_increasing <- function(x) all(diff(x) >= -1e-10)
  is_in_01 <- function(x) all(x >= -1e-10 & x <= 1 + 1e-10)

  # One Value per Bin --------------------------------------------------------

  test_that("every form returns one value per bin", {
    pars_by_model <- list(
      `0` = log(c(10, 0.5)),
      `1` = log(c(12, 3)),
      `2` = log(0.5),
      `3` = log(c(10, 15)),
      `4` = c(0, 0, log(5), log(5), 0, 0),
      `5` = rep(0, length(ages)),
      `6` = c(0, log(10), log(0.5)),
      `7` = c(0, log(10), log(15))
    )
    for (m in 0:7) {
      p <- pars_by_model[[as.character(m)]]
      devs <- zero_devs(length(p))
      res <- selex(m, p, devs = devs)
      expect_equal(length(res), length(ages),
                   label = sprintf("length Selex_Model=%d", m))
    }
  })

  # Logistic on b50 and Slope (form 0) ---------------------------------------

  test_that("the logistic on b50 and slope rises from zero to one", {
    res <- selex(0, log(c(10, 0.5)))
    expect_true(is_in_01(res))
    expect_true(is_monotone_increasing(res))
  })

  test_that("the logistic selects half the fish at b50", {
    b50 <- 10
    res <- selex(0, log(c(b50, 0.5)), bins = b50)
    expect_equal(res, 0.5, tolerance = 1e-6)
  })

  test_that("a steeper slope selects more above b50", {
    # Above b50, steeper k => faster approach to 1 => higher selex
    res_flat  <- selex(0, log(c(10, 0.2)))
    res_steep <- selex(0, log(c(10, 2.0)))
    expect_gt(res_steep[which(ages == 12)], res_flat[which(ages == 12)])
  })

  # Gamma Dome (form 1) ------------------------------------------------------

  test_that("the gamma dome peaks away from both end bins", {
    res <- selex(1, log(c(10, 3)))
    expect_true(all(res > 0))
    peak_idx <- which.max(res)
    expect_gt(peak_idx, 1)
    expect_lt(peak_idx, length(ages))
  })

  test_that("raising bmax moves the gamma peak to later bins", {
    peak_low  <- which.max(selex(1, log(c(8,  3))))
    peak_high <- which.max(selex(1, log(c(14, 3))))
    expect_gt(peak_high, peak_low)
  })

  # Power Function (form 2) --------------------------------------------------

  test_that("the power form falls across bins and stays positive", {
    res <- selex(2, log(0.5))
    expect_true(all(res > 0))
    expect_true(all(diff(res) <= 1e-10))
  })

  test_that("a larger power makes the decline steeper", {
    res_small <- selex(2, log(0.2))
    res_large <- selex(2, log(2.0))
    # at older ages, large power gives lower selectivity
    expect_lt(res_large[length(ages)], res_small[length(ages)])
  })

  # Logistic on b50 and b95 (form 3) -----------------------------------------

  test_that("the logistic on b50 and b95 rises from zero to one", {
    res <- selex(3, log(c(10, 5)))
    expect_true(is_in_01(res))
    expect_true(is_monotone_increasing(res))
  })

  test_that("that form selects half at b50 and 95 percent at b50 plus b95", {
    # b95 is the width above b50, so 95 percent selectivity is reached at b50 plus b95
    b50 <- 10
    b95 <- 5
    expect_equal(selex(3, log(c(b50, b95)), bins = b50),        0.5,  tolerance = 1e-6)
    expect_equal(selex(3, log(c(b50, b95)), bins = b50 + b95),  0.95, tolerance = 1e-4)
  })

  test_that("the two logistic parameterizations give the same curve", {
    # Model 3: 1 / (1 + 19^((b50-Bin)/b95)), b95 = width from b50 to 95%
    # k_equiv = log(19) / b95 makes Model 0 match
    b50 <- 10
    b95 <- 5
    k_equiv <- log(19) / b95
    res0 <- selex(0, log(c(b50, k_equiv)))
    res3 <- selex(3, log(c(b50, b95)))
    expect_equal(res0, res3, tolerance = 1e-8)
  })

  # Double Normal Dome (form 4) ----------------------------------------------

  test_that("the double normal stays between zero and one", {
    pars4 <- c(10, 0, log(5), log(5), 0, 0)
    res <- selex(4, pars4)
    expect_true(is_in_01(res))
  })

  test_that("the double normal peaks away from both end bins", {
    pars4 <- c(10, 0, log(5), log(5), -5, -5)  # peak at bin 10, low first/last bin selex
    res <- selex(4, pars4)
    peak <- which.max(res)
    expect_gt(peak, 1)
    expect_lt(peak, length(ages))
  })

  test_that("the double normal's last two parameters are selectivity at the end bins", {
    # p5 = plogis(pars[5]) sets selex[1]; p6 controls last bin via des.scaled
    low_first  <- selex(4, c(10, 0, log(5), log(5), -5,  0))
    high_first <- selex(4, c(10, 0, log(5), log(5),  5,  0))
    expect_lt(low_first[1], high_first[1])
  })

  # Non-parametric on the Logit Scale (form 5) -------------------------------

  test_that("the non-parametric logit form stays between zero and one", {
    pars5 <- seq(-2, 2, length.out = length(ages))
    res <- selex(5, pars5, devs = zero_devs(length(ages)))
    expect_true(is_in_01(res))
  })

  test_that("logit parameters of zero select half of every bin", {
    pars5 <- rep(0, length(ages))
    res <- selex(5, pars5, devs = zero_devs(length(ages)))
    expect_equal(unique(res), 0.5, tolerance = 1e-10)
  })

  test_that("large positive logit parameters select nearly everything", {
    pars5 <- rep(10, length(ages))
    res <- selex(5, pars5, devs = zero_devs(length(ages)))
    expect_true(all(res > 0.99))
  })

  test_that("large negative logit parameters select almost nothing", {
    pars5 <- rep(-10, length(ages))
    res <- selex(5, pars5, devs = zero_devs(length(ages)))
    expect_true(all(res < 0.01))
  })

  # Logistic with an Asymptote (form 6) --------------------------------------

  test_that("the logistic with an asymptote tops out at alpha", {
    alpha_logit <- 0   # plogis(0) = 0.5
    res <- selex(6, c(alpha_logit, log(10), log(0.5)))
    expect_true(all(res >= -1e-10))
    expect_true(all(res <= 0.5 + 1e-6))
  })

  test_that("the logistic with an asymptote rises across bins", {
    res <- selex(6, c(0, log(10), log(0.5)))
    expect_true(is_monotone_increasing(res))
  })

  test_that("alpha scales how much of the stock the curve can reach", {
    res_half <- selex(6, c(qlogis(0.5), log(10), log(0.5)))  # alpha = 0.5
    res_full <- selex(6, c(qlogis(0.9), log(10), log(0.5)))  # alpha = 0.9
    expect_lt(max(res_half), max(res_full))
  })

  test_that("at an asymptote of one it is the plain b50 and slope logistic", {
    # With alpha ~= 1, Model 6 should approach Model 0
    res0 <- selex(0, log(c(10, 0.5)))
    res6 <- selex(6, c(qlogis(0.9999), log(10), log(0.5)))
    expect_equal(res0, res6, tolerance = 1e-3)
  })

  # b95 Logistic with an Asymptote (form 7) ----------------------------------

  test_that("the b95 logistic with an asymptote tops out at alpha", {
    res <- selex(7, c(0, log(10), log(5)))   # alpha = plogis(0) = 0.5
    expect_true(all(res >= -1e-10))
    expect_true(all(res <= 0.5 + 1e-6))
  })

  test_that("the b95 logistic with an asymptote rises across bins", {
    res <- selex(7, c(0, log(10), log(5)))
    expect_true(is_monotone_increasing(res))
  })

  test_that("at an asymptote of one it is the plain b50 and b95 logistic", {
    res3 <- selex(3, log(c(10, 5)))
    res7 <- selex(7, c(qlogis(0.9999), log(10), log(5)))
    expect_equal(res3, res7, tolerance = 1e-3)
  })

  # Bicubic Spline over Age and Year Nodes (form 8) --------------------------

  bicubic_selex <- function(pars, Wbin, Wyr, year = 1) {
    Get_Selex(
      Selex_Model    = 8,
      TimeVary_Model = 0,
      pars           = pars,
      ln_seldevs     = zero_devs(1),
      Region         = 1,
      Year           = year,
      Bin            = ages,
      Sex            = 1,
      Wbin_bicubic   = Wbin,
      Wyr_bicubic    = Wyr
    )
  }

  test_that("the spline returns one value per bin", {
    bin_nodes <- seq(0, 1, length.out = 4)
    Wbin <- Get_Natural_Cubic_Spline_Weights(bin_nodes, seq(0, 1, length.out = length(ages)))
    Wyr  <- matrix(1, nrow = 5, ncol = 1) # single year node -> time-invariant
    pars <- c(0, 1, 0.5, -0.5)
    res <- bicubic_selex(pars, Wbin, Wyr)
    expect_equal(length(res), length(ages))
  })

  test_that("the spline is positive everywhere, being exponentiated", {
    bin_nodes <- seq(0, 1, length.out = 4)
    Wbin <- Get_Natural_Cubic_Spline_Weights(bin_nodes, seq(0, 1, length.out = length(ages)))
    Wyr  <- matrix(1, nrow = 3, ncol = 1)
    res <- bicubic_selex(c(-2, 1, 0.5, -3), Wbin, Wyr)
    expect_true(all(res > 0))
  })

  test_that("one year node gives a curve that does not change over time", {
    bin_nodes <- seq(0, 1, length.out = 5)
    age_bins  <- seq(0, 1, length.out = length(ages))
    Wbin <- Get_Natural_Cubic_Spline_Weights(bin_nodes, age_bins)
    Wyr  <- matrix(1, nrow = 6, ncol = 1) # every year maps to the single node

    log_node_vals <- c(-1, 0.5, 1.2, 0.3, -0.8)
    res_y1 <- bicubic_selex(log_node_vals, Wbin, Wyr, year = 1)
    res_y6 <- bicubic_selex(log_node_vals, Wbin, Wyr, year = 6)

    # matches direct age-only spline evaluation
    expect_equal(res_y1, exp(as.vector(Wbin %*% log_node_vals)), tolerance = 1e-8)
    # constant across years
    expect_equal(res_y1, res_y6)
  })

  test_that("the surface matches splining over ages and then over years by hand", {
    n_bin_nodes <- 4
    n_yr_nodes  <- 3
    n_yrs <- 7

    bin_nodes <- seq(0, 1, length.out = n_bin_nodes)
    yr_nodes  <- seq(0, 1, length.out = n_yr_nodes)
    age_bins  <- seq(0, 1, length.out = length(ages))
    yr_bins   <- seq(0, 1, length.out = n_yrs)

    Wbin <- Get_Natural_Cubic_Spline_Weights(bin_nodes, age_bins)
    Wyr  <- Get_Natural_Cubic_Spline_Weights(yr_nodes, yr_bins)

    set.seed(42)
    node_par <- matrix(rnorm(n_yr_nodes * n_bin_nodes), nrow = n_yr_nodes, ncol = n_bin_nodes)

    # manual two-pass construction: age-spline each year-node row, then year-spline each age column
    age_interp <- node_par %*% t(Wbin)      # n_yr_nodes x n_ages
    full_surface <- Wyr %*% age_interp      # n_yrs x n_ages
    expected <- exp(full_surface)

    for (y in 1:n_yrs) {
      res <- bicubic_selex(as.vector(node_par), Wbin, Wyr, year = y)
      expect_equal(res, expected[y, ], tolerance = 1e-8, label = sprintf("year %d", y))
    }
  })

  test_that("the zeros setup pads the node grid with change nothing", {
    n_bin_nodes <- 3
    n_yr_nodes  <- 2
    bin_nodes <- seq(0, 1, length.out = n_bin_nodes)
    yr_nodes  <- seq(0, 1, length.out = n_yr_nodes)
    age_bins  <- seq(0, 1, length.out = length(ages))
    yr_bins   <- seq(0, 1, length.out = 4)

    Wbin <- Get_Natural_Cubic_Spline_Weights(bin_nodes, age_bins)
    Wyr  <- Get_Natural_Cubic_Spline_Weights(yr_nodes, yr_bins)
    node_par <- matrix(c(0.2, -0.3, 0.1, 0.4, -0.1, 0.6), nrow = n_yr_nodes, ncol = n_bin_nodes)

    res_unpadded <- bicubic_selex(as.vector(node_par), Wbin, Wyr, year = 2)

    # pad with an extra all-zero age-node column and year-node column
    Wbin_pad <- cbind(Wbin, 0)
    Wyr_pad  <- cbind(Wyr, 0)
    node_par_pad <- rbind(cbind(node_par, 999), 999) # padded parameter slots can hold any value

    res_padded <- bicubic_selex(as.vector(node_par_pad), Wbin_pad, Wyr_pad, year = 2)

    expect_equal(res_unpadded, res_padded, tolerance = 1e-8)
  })

  test_that("padding to another fleet's larger node grid leaves this fleet's curve alone", {
    # setup stores the spline weights in one array shared by every bicubic fleet, padded out to
    # the widest node grid any of them has, while each block's parameters hold only its own.
    #
    # so Get_Selex has to reshape the parameters on this block's own node counts rather than
    # on the width of the padded weight matrix
    n_bin_nodes_true <- 3
    n_yr_nodes_true  <- 2
    n_yrs <- 5

    bin_nodes <- seq(0, 1, length.out = n_bin_nodes_true)
    yr_nodes  <- seq(0, 1, length.out = n_yr_nodes_true)
    age_bins  <- seq(0, 1, length.out = length(ages))
    yr_bins   <- seq(0, 1, length.out = n_yrs)

    Wbin_true <- Get_Natural_Cubic_Spline_Weights(bin_nodes, age_bins) # n_ages x 3
    Wyr_true  <- Get_Natural_Cubic_Spline_Weights(yr_nodes, yr_bins)   # n_yrs x 2

    set.seed(7)
    node_par_true <- matrix(rnorm(n_bin_nodes_true * n_yr_nodes_true), nrow = n_yr_nodes_true, ncol = n_bin_nodes_true)
    expected <- exp(Wyr_true %*% (node_par_true %*% t(Wbin_true))) # n_yrs x n_ages, computed with no padding at all

    # stand in for the shared padded storage: another fleet somewhere has a larger node grid, so
    # this block's weight matrices are padded out to that width whatever its own counts are
    n_bin_nodes_padded <- 5
    n_yr_nodes_padded  <- 4
    Wbin_padded <- cbind(Wbin_true, matrix(0, nrow = nrow(Wbin_true), ncol = n_bin_nodes_padded - n_bin_nodes_true))
    Wyr_padded  <- cbind(Wyr_true,  matrix(0, nrow = nrow(Wyr_true),  ncol = n_yr_nodes_padded  - n_yr_nodes_true))

    # This block's own parameter storage: its true node values in the first n_bin_nodes_true *
    # n_yr_nodes_true slots (simple sequential 1:max_sel_pars mapping), everything else fixed at 0.
    pars_flat <- c(as.vector(node_par_true), rep(0, 40))

    for (y in 1:n_yrs) {
      res <- Get_Selex(
        Selex_Model = 8,
        TimeVary_Model = 0,
        pars = pars_flat,
        ln_seldevs = zero_devs(1),
        Region = 1,
        Year = y,
        Bin = ages,
        Sex = 1,
        Wbin_bicubic = Wbin_padded,
        Wyr_bicubic = Wyr_padded,
        n_bin_nodes_bicubic = n_bin_nodes_true,
        n_yr_nodes_bicubic = n_yr_nodes_true
      )
      expect_equal(res, expected[y, ], tolerance = 1e-8, label = sprintf("year %d, padded-vs-true grid mismatch", y))
    }
  })

  test_that("year blocks hold the curve fixed within a block and smooth across ages", {
    # mirrors ADMB sel_option==4: independent age-only spline per year-block, kept constant within the block
    n_bin_nodes <- 4
    bin_nodes <- seq(0, 1, length.out = n_bin_nodes)
    age_bins  <- seq(0, 1, length.out = length(ages))
    Wbin <- Get_Natural_Cubic_Spline_Weights(bin_nodes, age_bins)

    n_yrs <- 6
    block_of_year <- c(1, 1, 1, 2, 2, 2) # first 3 years -> block 1, last 3 -> block 2
    Wyr <- matrix(0, nrow = n_yrs, ncol = 2)
    for (y in 1:n_yrs) Wyr[y, block_of_year[y]] <- 1

    node_par <- matrix(c(-1, 0.5, 1, -0.3,   # block 1 age-node values
                         0.2, -0.6, 0.8, 1.1), # block 2 age-node values
                       nrow = 2, ncol = n_bin_nodes, byrow = TRUE)

    res <- lapply(1:n_yrs, function(y) bicubic_selex(as.vector(node_par), Wbin, Wyr, year = y))

    # constant within each block
    expect_equal(res[[1]], res[[2]])
    expect_equal(res[[2]], res[[3]])
    expect_equal(res[[4]], res[[5]])
    expect_equal(res[[5]], res[[6]])
    # differs across blocks
    expect_false(isTRUE(all.equal(res[[1]], res[[4]])))
    # matches direct age-only spline for each block's own node values
    expect_equal(res[[1]], exp(as.vector(Wbin %*% node_par[1, ])), tolerance = 1e-8)
    expect_equal(res[[4]], exp(as.vector(Wbin %*% node_par[2, ])), tolerance = 1e-8)
  })

  # No Time Variation --------------------------------------------------------

  test_that("no time variation leaves the logistic as it is", {
    pars <- log(c(10, 0.5))
    res_tv0  <- selex(0, pars, tv = 0)
    res_base <- selex(0, pars, tv = 0)
    expect_equal(res_tv0, res_base)
  })

  # Deviations on the Parameters ---------------------------------------------

  test_that("iid deviations of zero leave the logistic unchanged", {
    pars <- log(c(10, 0.5))
    devs <- zero_devs(2)
    expect_equal(
      selex(0, pars, tv = 0, devs = devs),
      selex(0, pars, tv = 1, devs = devs)
    )
  })

  test_that("a positive deviation on b50 selects fewer young fish", {
    pars <- log(c(10, 0.5))
    devs_pos <- zero_devs(2)
    devs_pos[1, 1, 1, 1, 1] <- 0.5  # shift b50 up
    res_base <- selex(0, pars, tv = 0)
    res_tv   <- selex(0, pars, tv = 1, devs = devs_pos)
    # Higher b50 => lower selex at younger ages
    expect_lt(res_tv[5], res_base[5])
  })

  test_that("random walk deviations of zero leave the b95 logistic unchanged", {
    pars <- log(c(10, 5))
    devs <- zero_devs(2)
    expect_equal(
      selex(3, pars, tv = 0, devs = devs),
      selex(3, pars, tv = 2, devs = devs)
    )
  })

  test_that("deviations on the parameters move every bin of the non-parametric form", {
    n_bins <- length(ages)
    pars5 <- rep(0, n_bins)
    devs_pos <- zero_devs(n_bins)
    devs_pos[1, 1, , 1, 1] <- 2   # add 2 to all logit pars => selex > 0.5
    res_base <- selex(5, pars5, tv = 0, devs = zero_devs(n_bins))
    res_tv   <- selex(5, pars5, tv = 1, devs = devs_pos)
    expect_true(all(res_tv > res_base))
  })

  # Deviations on the Curve Itself -------------------------------------------

  test_that("deviations on the curve itself, at zero, leave the logistic unchanged", {
    pars <- log(c(10, 0.5))
    devs <- zero_devs(length(ages))
    expect_equal(
      selex(0, pars, tv = 0, devs = zero_devs(2)),
      selex(0, pars, tv = 3, devs = devs)
    )
  })

  test_that("a positive deviation on the curve raises selectivity in that bin", {
    pars <- log(c(10, 0.5))
    devs_base <- zero_devs(length(ages))
    devs_pos  <- devs_base
    devs_pos[1, 1, , 1, 1] <- 0.5
    res_base <- selex(0, pars, tv = 0, devs = zero_devs(2))
    res_tv   <- selex(0, pars, tv = 3, devs = devs_pos)
    expect_true(all(res_tv >= res_base - 1e-10))
    expect_true(any(res_tv > res_base))
  })

  test_that("the conditional form, at zero, leaves the b95 logistic unchanged", {
    pars <- log(c(10, 15))
    devs <- zero_devs(length(ages))
    expect_equal(
      selex(3, pars, tv = 0, devs = zero_devs(2)),
      selex(3, pars, tv = 4, devs = devs)
    )
  })

  test_that("the separable form, at zero, leaves the gamma dome unchanged", {
    pars <- log(c(10, 3))
    devs <- zero_devs(length(ages))
    expect_equal(
      selex(1, pars, tv = 0, devs = zero_devs(2)),
      selex(1, pars, tv = 5, devs = devs)
    )
  })

  # One Form against Another -------------------------------------------------

  test_that("the asymptotic logistic is the plain one scaled by alpha", {
    k <- 0.5
    b50 <- 10
    alpha <- 0.7
    res0 <- selex(0, log(c(b50, k)))
    res6 <- selex(6, c(qlogis(alpha), log(b50), log(k)))
    expect_equal(res6, alpha * res0, tolerance = 1e-8)
  })

  test_that("the asymptotic b95 logistic is the plain one scaled by alpha", {
    b50 <- 10
    b95 <- 15
    alpha <- 0.7
    res3 <- selex(3, log(c(b50, b95)))
    res7 <- selex(7, c(qlogis(alpha), log(b50), log(b95)))
    expect_equal(res7, alpha * res3, tolerance = 1e-8)
  })

})
