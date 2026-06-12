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
