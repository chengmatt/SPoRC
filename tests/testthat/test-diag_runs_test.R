library(SPoRC)
library(testthat)

# do_runs_test() returns three-sigma control limits and a runs-test p-value for
# a residual series. The limits come from the average moving range, so they are
# checkable against a hand calculation.

amr_limits <- function(x, mu) {
  # Mirrors the Nelson (1982) moving-range calculation the function uses.
  mr <- abs(diff(x - mu))
  mr <- mr[mr < 3.267 * mean(mr, na.rm = TRUE)]
  stdev <- mean(mr, na.rm = TRUE) / 1.128
  c(mu - 3 * stdev, mu + 3 * stdev)
}

test_that("control limits are centered on zero for residual-type input", {
  set.seed(42)
  x <- rnorm(50)
  out <- SPoRC::do_runs_test(x)
  expect_named(out, c("sig3lim", "p.runs"))
  expect_length(out$sig3lim, 2)
  expect_lt(out$sig3lim[1], 0)
  expect_gt(out$sig3lim[2], 0)
  # Centered on mu = 0, so the limits are symmetric.
  expect_equal(out$sig3lim[1], -out$sig3lim[2])
  expect_equal(out$sig3lim, amr_limits(x, 0))
})

test_that("type defaults to residual when NULL", {
  set.seed(7)
  x <- rnorm(30)
  expect_equal(SPoRC::do_runs_test(x, type = NULL), SPoRC::do_runs_test(x, type = "resid"))
})

test_that("a non-residual type centers the limits on the series mean", {
  set.seed(11)
  x <- rnorm(40, mean = 5)
  out <- SPoRC::do_runs_test(x, type = "mean")
  expect_equal(out$sig3lim, amr_limits(x, mean(x)))
  # The limits bracket the sample mean rather than zero.
  expect_lt(out$sig3lim[1], mean(x))
  expect_gt(out$sig3lim[2], mean(x))
  expect_gt(out$sig3lim[1], 0)
})

test_that("an all-positive series short-circuits to a significant p-value", {
  # With no sign changes there are no runs to test, so the function reports the
  # floor value of 0.001 rather than calling randtests.
  out <- SPoRC::do_runs_test(c(1, 2, 3, 4, 5))
  expect_equal(out$p.runs, 0.001)
})

test_that("an all-negative series short-circuits the same way", {
  out <- SPoRC::do_runs_test(c(-1, -2, -3, -4, -5))
  expect_equal(out$p.runs, 0.001)
})

test_that("a strongly alternating series is flagged as non-random", {
  x <- rep(c(1, -1), 15)
  out <- SPoRC::do_runs_test(x)
  expect_lt(out$p.runs, 0.05)
})

test_that("independent noise is not flagged as non-random", {
  set.seed(123)
  out <- SPoRC::do_runs_test(rnorm(100))
  expect_gt(out$p.runs, 0.05)
})

test_that("a one-sided test is available for positive autocorrelation", {
  set.seed(99)
  x <- as.numeric(arima.sim(list(ar = 0.8), n = 100))
  two <- SPoRC::do_runs_test(x, mixing = "two.sided")
  less <- SPoRC::do_runs_test(x, mixing = "less")
  # Both flag the series, and the limits do not depend on the alternative.
  expect_lt(less$p.runs, 0.05)
  expect_equal(two$sig3lim, less$sig3lim)
})

test_that("the p-value is rounded to three decimals", {
  set.seed(5)
  out <- SPoRC::do_runs_test(rnorm(60))
  expect_equal(out$p.runs, round(out$p.runs, 3))
})

test_that("a two-element series returns a NaN p-value", {
  # randtests::runs.test cannot compute a variance from a single run, and the
  # NA guard inside do_runs_test assigns to `p.value` while the returned value
  # is built from `pvalue`, so the guard does not take effect. Pinned as-is:
  # callers must treat p.runs as possibly non-finite for very short series.
  out <- SPoRC::do_runs_test(c(1, -1))
  expect_true(is.nan(out$p.runs))
})
