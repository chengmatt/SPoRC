library(SPoRC)
library(testthat)

test_that("Get_Selex works", {

  ages  <- 1:20
  lens  <- seq(10, 80, by = 5)

  # Zero deviations array: [n_regions, n_years, n_pars_or_bins, n_sexes, 1]
  zero_devs <- function(n_pars, n_regions = 2, n_years = 5, n_sexes = 1) {
    array(0, dim = c(n_regions, n_years, n_pars, n_sexes, 1))
  }

  # Convenience wrapper
  selex <- function(model, pars, bins = ages,
                    tv = 0, devs = NULL,
                    region = 1, year = 1, sex = 1) {
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

  # ── output length ─────────────────────────────────────────────────────────────

  test_that("output length equals length(Bin) for all Selex_Models", {
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

  # ── Model 0: logistic (b50, slope) ───────────────────────────────────────────

  test_that("Model 0: values in [0,1] and monotone increasing", {
    res <- selex(0, log(c(10, 0.5)))
    expect_true(is_in_01(res))
    expect_true(is_monotone_increasing(res))
  })

  test_that("Model 0: value at b50 is ~0.5", {
    b50 <- 10
    res <- selex(0, log(c(b50, 0.5)), bins = b50)
    expect_equal(res, 0.5, tolerance = 1e-6)
  })

  test_that("Model 0: steeper slope gives higher selex at bins above b50", {
    # Above b50, steeper k => faster approach to 1 => higher selex
    res_flat  <- selex(0, log(c(10, 0.2)))
    res_steep <- selex(0, log(c(10, 2.0)))
    expect_gt(res_steep[which(ages == 12)], res_flat[which(ages == 12)])
  })

  # ── Model 1: gamma dome ───────────────────────────────────────────────────────

  test_that("Model 1: values > 0, dome shape (peak not at first or last bin)", {
    res <- selex(1, log(c(10, 3)))
    expect_true(all(res > 0))
    peak_idx <- which.max(res)
    expect_gt(peak_idx, 1)
    expect_lt(peak_idx, length(ages))
  })

  test_that("Model 1: peak shifts right when bmax increases", {
    peak_low  <- which.max(selex(1, log(c(8,  3))))
    peak_high <- which.max(selex(1, log(c(14, 3))))
    expect_gt(peak_high, peak_low)
  })

  # ── Model 2: power function ───────────────────────────────────────────────────

  test_that("Model 2: monotone decreasing, all positive", {
    res <- selex(2, log(0.5))
    expect_true(all(res > 0))
    expect_true(all(diff(res) <= 1e-10))
  })

  test_that("Model 2: larger power => steeper decline", {
    res_small <- selex(2, log(0.2))
    res_large <- selex(2, log(2.0))
    # at older ages, large power gives lower selectivity
    expect_lt(res_large[length(ages)], res_small[length(ages)])
  })

  # ── Model 3: logistic (b50, b95) ─────────────────────────────────────────────

  test_that("Model 3: values in [0,1] and monotone increasing", {
    res <- selex(3, log(c(10, 5)))
    expect_true(is_in_01(res))
    expect_true(is_monotone_increasing(res))
  })

  test_that("Model 3: value at b50 is ~0.5 and at b50+b95 is ~0.95", {
    # b95 is the WIDTH from b50, so 95% selectivity is reached at bin = b50 + b95
    b50 <- 10; b95 <- 5
    expect_equal(selex(3, log(c(b50, b95)), bins = b50),        0.5,  tolerance = 1e-6)
    expect_equal(selex(3, log(c(b50, b95)), bins = b50 + b95),  0.95, tolerance = 1e-4)
  })

  test_that("Models 0 and 3 agree when slope parameterizations are equivalent", {
    # Model 3: 1 / (1 + 19^((b50-Bin)/b95)), b95 = width from b50 to 95%
    # k_equiv = log(19) / b95 makes Model 0 match
    b50 <- 10; b95 <- 5
    k_equiv <- log(19) / b95
    res0 <- selex(0, log(c(b50, k_equiv)))
    res3 <- selex(3, log(c(b50, b95)))
    expect_equal(res0, res3, tolerance = 1e-8)
  })

  # ── Model 4: double-normal dome ──────────────────────────────────────────────

  test_that("Model 4: values in [0,1]", {
    pars4 <- c(0, 0, log(5), log(5), 0, 0)
    res <- selex(4, pars4)
    expect_true(is_in_01(res))
  })

  test_that("Model 4: dome-shaped (max not at bin 1 or last bin)", {
    pars4 <- c(0, 0, log(5), log(5), -5, -5)  # low first/last bin selex
    res <- selex(4, pars4)
    peak <- which.max(res)
    expect_gt(peak, 1)
    expect_lt(peak, length(ages))
  })

  test_that("Model 4: p5/p6 control first/last bin selectivity", {
    # p5 = plogis(pars[5]) sets selex[1]; p6 controls last bin via des.scaled
    low_first  <- selex(4, c(0, 0, log(5), log(5), -5,  0))
    high_first <- selex(4, c(0, 0, log(5), log(5),  5,  0))
    expect_lt(low_first[1], high_first[1])
  })

  # ── Model 5: non-parametric ───────────────────────────────────────────────────

  test_that("Model 5: values in [0,1]", {
    pars5 <- seq(-2, 2, length.out = length(ages))
    res <- selex(5, pars5, devs = zero_devs(length(ages)))
    expect_true(is_in_01(res))
  })

  test_that("Model 5: all-zero logit pars give selex = 0.5 everywhere", {
    pars5 <- rep(0, length(ages))
    res <- selex(5, pars5, devs = zero_devs(length(ages)))
    expect_equal(unique(res), 0.5, tolerance = 1e-10)
  })

  test_that("Model 5: large positive logit pars give selex near 1", {
    pars5 <- rep(10, length(ages))
    res <- selex(5, pars5, devs = zero_devs(length(ages)))
    expect_true(all(res > 0.99))
  })

  test_that("Model 5: large negative logit pars give selex near 0", {
    pars5 <- rep(-10, length(ages))
    res <- selex(5, pars5, devs = zero_devs(length(ages)))
    expect_true(all(res < 0.01))
  })

  # ── Model 6: logistic with asymptote (b50, k) ────────────────────────────────

  test_that("Model 6: values in [0, alpha]", {
    alpha_logit <- 0   # plogis(0) = 0.5
    res <- selex(6, c(alpha_logit, log(10), log(0.5)))
    expect_true(all(res >= -1e-10))
    expect_true(all(res <= 0.5 + 1e-6))
  })

  test_that("Model 6: monotone increasing", {
    res <- selex(6, c(0, log(10), log(0.5)))
    expect_true(is_monotone_increasing(res))
  })

  test_that("Model 6: asymptote parameter scales max selectivity", {
    res_half <- selex(6, c(qlogis(0.5), log(10), log(0.5)))  # alpha = 0.5
    res_full <- selex(6, c(qlogis(0.9), log(10), log(0.5)))  # alpha = 0.9
    expect_lt(max(res_half), max(res_full))
  })

  test_that("Model 6 and Model 0 agree when alpha -> 1", {
    # With alpha ~= 1, Model 6 should approach Model 0
    res0 <- selex(0, log(c(10, 0.5)))
    res6 <- selex(6, c(qlogis(0.9999), log(10), log(0.5)))
    expect_equal(res0, res6, tolerance = 1e-3)
  })

  # ── Model 7: logistic with asymptote (b50, b95) ──────────────────────────────

  test_that("Model 7: values in [0, alpha]", {
    res <- selex(7, c(0, log(10), log(5)))   # alpha = plogis(0) = 0.5
    expect_true(all(res >= -1e-10))
    expect_true(all(res <= 0.5 + 1e-6))
  })

  test_that("Model 7: monotone increasing", {
    res <- selex(7, c(0, log(10), log(5)))
    expect_true(is_monotone_increasing(res))
  })

  test_that("Model 7 and Model 3 agree when alpha -> 1", {
    res3 <- selex(3, log(c(10, 5)))
    res7 <- selex(7, c(qlogis(0.9999), log(10), log(5)))
    expect_equal(res3, res7, tolerance = 1e-3)
  })

  # ── Model 8: bicubic spline (age node x year node grid) ──────────────────────

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

  test_that("output length equals length(Bin)", {
    bin_nodes <- seq(0, 1, length.out = 4)
    Wbin <- Get_Natural_Cubic_Spline_Weights(bin_nodes, seq(0, 1, length.out = length(ages)))
    Wyr  <- matrix(1, nrow = 5, ncol = 1) # single year node -> time-invariant
    pars <- c(0, 1, 0.5, -0.5)
    res <- bicubic_selex(pars, Wbin, Wyr)
    expect_equal(length(res), length(ages))
  })

  test_that("values are strictly positive (log-scale spline, exponentiated)", {
    bin_nodes <- seq(0, 1, length.out = 4)
    Wbin <- Get_Natural_Cubic_Spline_Weights(bin_nodes, seq(0, 1, length.out = length(ages)))
    Wyr  <- matrix(1, nrow = 3, ncol = 1)
    res <- bicubic_selex(c(-2, 1, 0.5, -3), Wbin, Wyr)
    expect_true(all(res > 0))
  })

  test_that("n_yr_nodes == 1 gives a time-invariant age-only natural cubic spline", {
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

  test_that("full bicubic surface matches manual two-pass (age-then-year) spline construction", {
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

  test_that("zero-padded extra node columns/rows are harmless (padding convention used by setup code)", {
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

  test_that("Wbin/Wyr padded to a *different* fleet/block's larger grid does not scramble this block's own node values", {
    # Regression test: in Setup_Fishery.R/Setup_Survey.R, Wbin_bicubic/Wyr_bicubic are stored in
    # arrays shared across all fleets/blocks using Selex_Model == 8, padded to the *global* max
    # bin-node/year-node width across all of them independently. A block whose own true grid is
    # smaller than that global max (in either dimension) receives Wbin_bicubic/Wyr_bicubic wider
    # than its own true n_bin_nodes/n_yr_nodes, while its own flattened parameter vector only ever
    # holds its own true n_bin_nodes * n_yr_nodes estimated values (see
    # do_fish_fixed_sel_pars_mapping's "max_sel_pars <- n_bin_nodes_this * n_yr_nodes_this").
    # Get_Selex must use this block's own true node counts (n_bin_nodes_bicubic/n_yr_nodes_bicubic)
    # to reshape pars, not ncol(Wbin_bicubic)/ncol(Wyr_bicubic), which reflect only the padded width.
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

    # Simulate the shared, globally-padded storage: some *other* fleet/block elsewhere uses a
    # larger bicubic grid (e.g. Bin_5_Yr_4), so this block's own Wbin/Wyr get zero-padded out to
    # that wider shared width, independent of this block's own true node counts.
    n_bin_nodes_padded <- 5
    n_yr_nodes_padded  <- 4
    Wbin_padded <- cbind(Wbin_true, matrix(0, nrow = nrow(Wbin_true), ncol = n_bin_nodes_padded - n_bin_nodes_true))
    Wyr_padded  <- cbind(Wyr_true,  matrix(0, nrow = nrow(Wyr_true),  ncol = n_yr_nodes_padded  - n_yr_nodes_true))

    # This block's own parameter storage: its true node values in the first n_bin_nodes_true *
    # n_yr_nodes_true slots (simple sequential 1:max_sel_pars mapping), everything else fixed at 0.
    pars_flat <- c(as.vector(node_par_true), rep(0, 40))

    for (y in 1:n_yrs) {
      res <- Get_Selex(
        Selex_Model = 8, TimeVary_Model = 0, pars = pars_flat, ln_seldevs = zero_devs(1),
        Region = 1, Year = y, Bin = ages, Sex = 1,
        Wbin_bicubic = Wbin_padded, Wyr_bicubic = Wyr_padded,
        n_bin_nodes_bicubic = n_bin_nodes_true, n_yr_nodes_bicubic = n_yr_nodes_true
      )
      expect_equal(res, expected[y, ], tolerance = 1e-8, label = sprintf("year %d, padded-vs-true grid mismatch", y))
    }
  })

  test_that("year-block (nearest-node indicator) Wyr gives piecewise-constant-in-year, smooth-in-age selectivity", {
    # mirrors ADMB sel_option==4: independent age-only spline per year-block, held constant within the block
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

  # ── TimeVary_Model = 0: no devs → same as base ───────────────────────────────

  test_that("TimeVary_Model=0 with zero devs equals base model (Model 0)", {
    pars <- log(c(10, 0.5))
    res_tv0  <- selex(0, pars, tv = 0)
    res_base <- selex(0, pars, tv = 0)
    expect_equal(res_tv0, res_base)
  })

  # ── TimeVary_Model = 1/2: parameter-level deviations ─────────────────────────

  test_that("TimeVary_Model=1 with zero devs equals TimeVary_Model=0 (Model 0)", {
    pars <- log(c(10, 0.5))
    devs <- zero_devs(2)
    expect_equal(
      selex(0, pars, tv = 0, devs = devs),
      selex(0, pars, tv = 1, devs = devs)
    )
  })

  test_that("TimeVary_Model=1 positive b50 dev shifts logistic right (Model 0)", {
    pars <- log(c(10, 0.5))
    devs_pos <- zero_devs(2); devs_pos[1, 1, 1, 1, 1] <- 0.5  # shift b50 up
    res_base <- selex(0, pars, tv = 0)
    res_tv   <- selex(0, pars, tv = 1, devs = devs_pos)
    # Higher b50 => lower selex at younger ages
    expect_lt(res_tv[5], res_base[5])
  })

  test_that("TimeVary_Model=2 with zero devs equals TimeVary_Model=0 (Model 3)", {
    pars <- log(c(10, 5))
    devs <- zero_devs(2)
    expect_equal(
      selex(3, pars, tv = 0, devs = devs),
      selex(3, pars, tv = 2, devs = devs)
    )
  })

  test_that("TimeVary_Model=1 modifies non-parametric pars (Model 5)", {
    n_bins <- length(ages)
    pars5 <- rep(0, n_bins)
    devs_pos <- zero_devs(n_bins)
    devs_pos[1, 1, , 1, 1] <- 2   # add 2 to all logit pars => selex > 0.5
    res_base <- selex(5, pars5, tv = 0, devs = zero_devs(n_bins))
    res_tv   <- selex(5, pars5, tv = 1, devs = devs_pos)
    expect_true(all(res_tv > res_base))
  })

  # ── TimeVary_Model = 3–5: bin-level multiplicative deviations ─────────────────

  test_that("TimeVary_Model=3 with zero devs equals TimeVary_Model=0 (Model 0)", {
    pars <- log(c(10, 0.5))
    devs <- zero_devs(length(ages))
    expect_equal(
      selex(0, pars, tv = 0, devs = zero_devs(2)),
      selex(0, pars, tv = 3, devs = devs)
    )
  })

  test_that("TimeVary_Model=3 positive bin devs increase selectivity (Model 0)", {
    pars <- log(c(10, 0.5))
    devs_base <- zero_devs(length(ages))
    devs_pos  <- devs_base; devs_pos[1, 1, , 1, 1] <- 0.5
    res_base <- selex(0, pars, tv = 0, devs = zero_devs(2))
    res_tv   <- selex(0, pars, tv = 3, devs = devs_pos)
    expect_true(all(res_tv >= res_base - 1e-10))
    expect_true(any(res_tv > res_base))
  })

  test_that("TimeVary_Model=4 with zero devs equals TimeVary_Model=0 (Model 3)", {
    pars <- log(c(10, 15))
    devs <- zero_devs(length(ages))
    expect_equal(
      selex(3, pars, tv = 0, devs = zero_devs(2)),
      selex(3, pars, tv = 4, devs = devs)
    )
  })

  test_that("TimeVary_Model=5 with zero devs equals TimeVary_Model=0 (Model 1)", {
    pars <- log(c(10, 3))
    devs <- zero_devs(length(ages))
    expect_equal(
      selex(1, pars, tv = 0, devs = zero_devs(2)),
      selex(1, pars, tv = 5, devs = devs)
    )
  })

  # ── cross-model sanity ────────────────────────────────────────────────────────

  test_that("Models 0 and 6 produce same shape; Model 6 is scaled by alpha", {
    k <- 0.5; b50 <- 10; alpha <- 0.7
    res0 <- selex(0, log(c(b50, k)))
    res6 <- selex(6, c(qlogis(alpha), log(b50), log(k)))
    expect_equal(res6, alpha * res0, tolerance = 1e-8)
  })

  test_that("Models 3 and 7 produce same shape; Model 7 is scaled by alpha", {
    b50 <- 10; b95 <- 15; alpha <- 0.7
    res3 <- selex(3, log(c(b50, b95)))
    res7 <- selex(7, c(qlogis(alpha), log(b50), log(b95)))
    expect_equal(res7, alpha * res3, tolerance = 1e-8)
  })

})
