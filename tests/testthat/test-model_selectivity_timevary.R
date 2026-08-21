library(SPoRC)
library(testthat)

# Get_Selex applies annual deviations to the transformed parameters when
# TimeVary_Model is 1 (iid) or 2 (random walk). The existing suite covers the
# deviation path only for the logistic and non-parametric forms, so every
# remaining Selex_Model is driven here.

bins <- 1:20
n_regions <- 2
n_years <- 5
n_sexes <- 1

# Deviations array: [region, year, par, sex, 1]
make_devs <- function(vals, region = 1, year = 2, sex = 1) {
  d <- array(0, dim = c(n_regions, n_years, length(vals), n_sexes, 1))
  d[region, year, , sex, 1] <- vals
  return(d)
}

selex <- function(model, pars, devs, tv = 1, year = 2, region = 1, sex = 1) {
  SPoRC:::Get_Selex(Selex_Model = model, TimeVary_Model = tv, pars = pars,
                    ln_seldevs = devs, Region = region, Year = year,
                    Bin = bins, Sex = sex)
}

test_that("deviations on log-scale parameters shift the parameters directly", {
  # For forms whose parameters are all exp()-transformed, multiplying by
  # exp(dev) is identical to adding the deviation on the log scale. That
  # equivalence pins the deviation wiring without restating each formula.
  cases <- list(
    list(model = 0, pars = log(c(10, 0.5)), devs = c(0.15, -0.2)),  # logistic b50/slope
    list(model = 1, pars = log(c(12, 3)),   devs = c(0.1, 0.25)),   # gamma dome
    list(model = 2, pars = log(0.5),        devs = 0.3),            # power
    list(model = 3, pars = log(c(10, 15)),  devs = c(-0.1, 0.2))    # logistic b50/b95
  )

  for(case in cases) {
    devs <- make_devs(case$devs)
    tv_sel <- selex(case$model, case$pars, devs, tv = 1)
    eq_sel <- selex(case$model, case$pars + case$devs, devs * 0, tv = 0)
    expect_equal(tv_sel, eq_sel, tolerance = 1e-12,
                 info = paste("Selex_Model", case$model))
    # The deviation must actually change the curve, or the test proves nothing.
    expect_false(isTRUE(all.equal(tv_sel, selex(case$model, case$pars, devs * 0, tv = 0))),
                 info = paste("Selex_Model", case$model))
  } # end case loop
})

test_that("the gamma dome peaks at the deviated bmax", {
  # bmax is the bin at maximum selectivity, so a positive deviation must push
  # the mode to a higher bin.
  pars <- log(c(8, 3))
  base <- selex(1, pars, make_devs(c(0, 0)), tv = 0)
  up <- selex(1, pars, make_devs(c(0.4, 0)), tv = 1)
  expect_gt(bins[which.max(up)], bins[which.max(base)])
})

test_that("the power form responds to its single deviation", {
  pars <- log(0.5)
  base <- selex(2, pars, make_devs(0), tv = 0)
  up <- selex(2, pars, make_devs(0.3), tv = 1)
  # selex = 1 / Bin^power, so a larger power lowers selectivity above bin 1.
  expect_true(all(up[-1] < base[-1]))
  expect_equal(up[1], base[1], tolerance = 1e-12)
})

test_that("the double normal applies a deviation to all six parameters", {
  pars <- c(10, 0, log(5), log(5), -1, -1)
  devs <- c(0.1, 0.2, 0.3, -0.1, 0.25, -0.15)
  d <- make_devs(devs)

  tv_sel <- selex(4, pars, d, tv = 1)
  base <- selex(4, pars, d * 0, tv = 0)

  expect_length(tv_sel, length(bins))
  expect_true(all(is.finite(tv_sel)))
  expect_false(isTRUE(all.equal(tv_sel, base)))

  # The ascending limb is rescaled so that it passes through the fifth
  # parameter at the first bin, which makes this a check on that parameter's
  # deviation up to the joiner's contribution there. The joiner is around
  # 1e-8 this far below the peak, so the comparison is to that.
  expect_equal(unname(tv_sel[1]), stats::plogis(pars[5]) * exp(devs[5]),
               tolerance = 1e-6)
  expect_equal(unname(base[1]), stats::plogis(pars[5]), tolerance = 1e-6)
})

test_that("the asymptotic logistic (b50/k) matches a hand calculation", {
  pars <- c(-0.5, log(10), log(0.5))
  devs <- c(-0.2, 0.15, 0.1)
  tv_sel <- selex(6, pars, make_devs(devs), tv = 1)

  alpha <- stats::plogis(pars[1]) * exp(devs[1])
  b50 <- exp(pars[2] + devs[2])
  k <- exp(pars[3] + devs[3])
  expect_equal(tv_sel, alpha / (1 + exp(-k * (bins - b50))), tolerance = 1e-12)
})

test_that("the asymptotic logistic (b50/b95) matches a hand calculation", {
  pars <- c(-0.5, log(10), log(15))
  devs <- c(-0.2, 0.15, 0.1)
  tv_sel <- selex(7, pars, make_devs(devs), tv = 1)

  alpha <- stats::plogis(pars[1]) * exp(devs[1])
  b50 <- exp(pars[2] + devs[2])
  b95 <- exp(pars[3] + devs[3])
  expect_equal(tv_sel, alpha / (1 + 19^((b50 - bins) / b95)), tolerance = 1e-12)
})

test_that("the asymptote deviation is multiplicative, not on the logit scale", {
  # alpha is plogis()-transformed but the deviation multiplies the result, so
  # it is not equivalent to shifting the parameter before the transform.
  pars <- c(-0.5, log(10), log(0.5))
  devs <- c(-0.2, 0, 0)
  tv_sel <- selex(6, pars, make_devs(devs), tv = 1)
  logit_shift <- selex(6, pars + devs, make_devs(c(0, 0, 0)), tv = 0)
  expect_false(isTRUE(all.equal(tv_sel, logit_shift)))
  expect_equal(max(tv_sel) / max(selex(6, pars, make_devs(c(0, 0, 0)), tv = 0)),
               exp(devs[1]), tolerance = 1e-6)
})

test_that("the non-parametric form adds deviations on the logit scale", {
  pars <- rep(0, length(bins))
  devs <- seq(-0.5, 0.5, length.out = length(bins))
  tv_sel <- selex(5, pars, make_devs(devs), tv = 1)
  expect_equal(tv_sel, stats::plogis(pars + devs), tolerance = 1e-12)
})

test_that("random walk deviations behave identically to iid deviations", {
  # TimeVary_Model 1 and 2 share a single branch; only the penalty differs.
  cases <- list(
    list(model = 0, pars = log(c(10, 0.5)), devs = c(0.15, -0.2)),
    list(model = 1, pars = log(c(12, 3)),   devs = c(0.1, 0.25)),
    list(model = 2, pars = log(0.5),        devs = 0.3),
    list(model = 4, pars = c(0, 0, log(5), log(5), -1, -1),
         devs = c(0.1, 0.2, 0.3, -0.1, 0.25, -0.15)),
    list(model = 6, pars = c(-0.5, log(10), log(0.5)), devs = c(-0.2, 0.15, 0.1)),
    list(model = 7, pars = c(-0.5, log(10), log(15)),  devs = c(-0.2, 0.15, 0.1))
  )

  for(case in cases) {
    d <- make_devs(case$devs)
    expect_equal(selex(case$model, case$pars, d, tv = 1),
                 selex(case$model, case$pars, d, tv = 2),
                 info = paste("Selex_Model", case$model))
  } # end case loop
})

test_that("deviations are read from the requested region, year and sex", {
  # A deviation stored against one year must not leak into another.
  pars <- log(c(10, 0.5))
  d <- make_devs(c(0.3, 0.3), region = 2, year = 4)

  hit <- selex(0, pars, d, tv = 1, region = 2, year = 4)
  miss_year <- selex(0, pars, d, tv = 1, region = 2, year = 1)
  miss_region <- selex(0, pars, d, tv = 1, region = 1, year = 4)
  none <- selex(0, pars, d * 0, tv = 0)

  expect_false(isTRUE(all.equal(hit, none)))
  expect_equal(miss_year, none, tolerance = 1e-12)
  expect_equal(miss_region, none, tolerance = 1e-12)
})

test_that("zero deviations leave every form unchanged", {
  cases <- list(
    list(model = 0, pars = log(c(10, 0.5))),
    list(model = 1, pars = log(c(12, 3))),
    list(model = 2, pars = log(0.5)),
    list(model = 3, pars = log(c(10, 15))),
    list(model = 4, pars = c(0, 0, log(5), log(5), -1, -1)),
    list(model = 5, pars = rep(0, length(bins))),
    list(model = 6, pars = c(-0.5, log(10), log(0.5))),
    list(model = 7, pars = c(-0.5, log(10), log(15)))
  )

  for(case in cases) {
    d <- make_devs(rep(0, length(case$pars)))
    expect_equal(selex(case$model, case$pars, d, tv = 1),
                 selex(case$model, case$pars, d, tv = 0),
                 tolerance = 1e-12, info = paste("Selex_Model", case$model))
  } # end case loop
})
