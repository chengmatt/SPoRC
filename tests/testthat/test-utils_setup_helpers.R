library(SPoRC)
library(testthat)

# extend_years() and convert_to_numeric() are pure setup helpers. Their fill
# strategies and coercion paths are reachable without building a model, so the
# edge cases are tested directly here.

test_that("extend_years appends n_years slices to the year dimension", {
  # n_years is the number of slices appended, not the resulting total.
  a <- array(1:8, dim = c(2, 4))
  out <- SPoRC:::extend_years(a, n_years = 3, yr_dim = 2, fill = "zeros")
  expect_equal(dim(out), c(2, 7))
  # abind attaches dimnames to the result, so compare values only.
  expect_equal(out[, 1:4], a, ignore_attr = TRUE)
  expect_true(all(out[, 5:7] == 0))
})

test_that("extend_years fills with zeros for both zeros and F_pattern", {
  # F_pattern defers the actual sample sizes to the closed-loop simulation, so
  # at setup it is indistinguishable from a zero fill.
  a <- array(1:12, dim = c(2, 3, 2))
  zeros <- SPoRC:::extend_years(a, n_years = 2, yr_dim = 2, fill = "zeros")
  fpat <- SPoRC:::extend_years(a, n_years = 2, yr_dim = 2, fill = "F_pattern")
  expect_equal(zeros, fpat)
  expect_equal(dim(zeros), c(2, 5, 2))
  expect_true(all(zeros[, 4:5, ] == 0))
})

test_that("extend_years 'last' repeats the final populated year slice", {
  a <- array(0, dim = c(2, 3))
  a[, 1] <- c(1, 2)
  a[, 2] <- c(3, 4)
  a[, 3] <- c(5, 6)
  out <- SPoRC:::extend_years(a, n_years = 2, yr_dim = 2, fill = "last")
  expect_equal(dim(out), c(2, 5))
  expect_equal(out[, 4], c(5, 6))
  expect_equal(out[, 5], c(5, 6))
})

test_that("extend_years 'last' skips trailing all-NA years", {
  # Biological arrays are often padded with NA past the last data year, and the
  # fill must reach back to the last year that actually holds values.
  a <- array(NA_real_, dim = c(2, 4))
  a[, 1] <- c(1, 2)
  a[, 2] <- c(3, 4)
  out <- SPoRC:::extend_years(a, n_years = 2, yr_dim = 2, fill = "last")
  expect_equal(out[, 5], c(3, 4))
  expect_equal(out[, 6], c(3, 4))
})

test_that("extend_years 'last' yields NA when no year holds a value", {
  a <- array(NA_real_, dim = c(2, 3))
  out <- SPoRC:::extend_years(a, n_years = 2, yr_dim = 2, fill = "last")
  expect_equal(dim(out), c(2, 5))
  expect_true(all(is.na(out[, 4:5])))
})

test_that("extend_years 'mean' averages over years excluding zeros and NA", {
  a <- array(0, dim = c(2, 4))
  a[1, ] <- c(2, 4, 0, NA)   # mean of 2 and 4
  a[2, ] <- c(10, NA, 0, 20) # mean of 10 and 20
  out <- SPoRC:::extend_years(a, n_years = 2, yr_dim = 2, fill = "mean")
  expect_equal(out[1, 5], 3)
  expect_equal(out[2, 5], 15)
  expect_equal(out[, 5], out[, 6])
})

test_that("extend_years 'mean' returns zero for an element with no valid years", {
  a <- array(0, dim = c(2, 3))
  a[1, ] <- c(1, 3, 0)
  a[2, ] <- c(0, NA, 0)
  out <- SPoRC:::extend_years(a, n_years = 1, yr_dim = 2, fill = "mean")
  expect_equal(out[1, 4], 2)
  expect_equal(out[2, 4], 0)
})

test_that("extend_years accepts a numeric constant fill", {
  a <- array(1:8, dim = c(2, 4))
  out <- SPoRC:::extend_years(a, n_years = 2, yr_dim = 2, fill = 7)
  expect_equal(dim(out), c(2, 6))
  expect_equal(out[, 1:4], a, ignore_attr = TRUE)
  expect_true(all(out[, 5:6] == 7))
})

test_that("extend_years extends along a non-trailing year dimension", {
  # Year is the third axis for most model arrays, so yr_dim is rarely the last.
  a <- array(1, dim = c(2, 3, 4, 2))
  out <- SPoRC:::extend_years(a, n_years = 2, yr_dim = 3, fill = "zeros")
  expect_equal(dim(out), c(2, 3, 6, 2))
  expect_true(all(out[, , 1:4, ] == 1))
  expect_true(all(out[, , 5:6, ] == 0))
})

test_that("convert_to_numeric passes numeric input through unchanged", {
  lookup <- list("none" = 999, "multinomial" = 0, "dirichlet" = 1)
  expect_equal(SPoRC:::convert_to_numeric(c(0, 1, 999), lookup), c(0, 1, 999))
  a <- array(1:4, dim = c(2, 2))
  expect_equal(SPoRC:::convert_to_numeric(a, lookup), a)
})

test_that("convert_to_numeric maps character labels to their codes", {
  lookup <- list("none" = 999, "multinomial" = 0, "dirichlet" = 1)
  out <- SPoRC:::convert_to_numeric(c("multinomial", "none", "dirichlet"), lookup)
  expect_equal(unname(out), c(0, 999, 1))
})

test_that("convert_to_numeric rejects unrecognised labels for a named-vector lookup", {
  lookup <- c("none" = 999, "multinomial" = 0, "dirichlet" = 1)
  expect_error(SPoRC:::convert_to_numeric("multinomal", lookup), "multinomal")
  expect_error(SPoRC:::convert_to_numeric("multinomal", lookup), "Valid options")
  # Every invalid entry is reported, not just the first.
  err <- expect_error(SPoRC:::convert_to_numeric(c("none", "bad1", "bad2"), lookup))
  expect_match(conditionMessage(err), "bad1")
  expect_match(conditionMessage(err), "bad2")
})

test_that("convert_to_numeric rejects unrecognised labels for a list lookup", {
  # Every caller in R/setup_*.R passes a list, so this is the path that matters.
  lookup <- list("none" = 999, "multinomial" = 0, "dirichlet" = 1)

  expect_error(SPoRC:::convert_to_numeric("multinomal", lookup), "multinomal")
  expect_error(SPoRC:::convert_to_numeric("multinomal", lookup), "Valid options")

  # A mixed vector errors rather than silently shrinking.
  err <- expect_error(SPoRC:::convert_to_numeric(c("none", "typo"), lookup))
  expect_match(conditionMessage(err), "typo")

  # Valid labels are unaffected.
  expect_equal(SPoRC:::convert_to_numeric(c("none", "dirichlet"), lookup), c(999, 1))
})

test_that("convert_to_numeric rejects input that is neither numeric nor character", {
  lookup <- list("none" = 999)
  expect_error(SPoRC:::convert_to_numeric(TRUE, lookup), "numeric or character")
  expect_error(SPoRC:::convert_to_numeric(list(1, 2), lookup), "numeric or character")
})

test_that("convert_to_numeric keeps the shape of a character array", {
  # Setup_Sim_Survey passes [n_yrs x n_srv_fleets] comp-type arrays straight to
  # check_sim_dimensions(), so the dim attribute has to survive the conversion.
  lookup <- list("none" = 999, "multinomial" = 0)
  x <- array(c("none", "multinomial", "none", "multinomial"), dim = c(2, 2))
  out <- SPoRC:::convert_to_numeric(x, lookup)
  expect_equal(dim(out), c(2, 2))
  expect_equal(as.vector(out), c(999, 0, 999, 0))
})

test_that("convert_to_numeric rejects a bad label inside a character array", {
  lookup <- list("none" = 999, "multinomial" = 0)
  x <- array(c("none", "multinomial", "none", "typo"), dim = c(2, 2))
  expect_error(SPoRC:::convert_to_numeric(x, lookup), "typo")
})

test_that("resolve_sel_pen_wts defaults every term to zero", {
  out <- SPoRC:::resolve_sel_pen_wts(NULL)
  expect_length(out, 6)
  expect_true(all(out == 0))
  expect_setequal(names(out),
                  c("smooth_bin_curve", "smooth_bin_diff", "smooth_yr_diff",
                    "smooth_yr_curve", "smooth_dome", "smooth_mean_center"))
})

test_that("resolve_sel_pen_wts fills supplied terms and zeroes the rest", {
  out <- SPoRC:::resolve_sel_pen_wts(c(smooth_dome = 2.5, smooth_yr_diff = 1))
  expect_equal(unname(out[["smooth_dome"]]), 2.5)
  expect_equal(unname(out[["smooth_yr_diff"]]), 1)
  expect_equal(unname(out[["smooth_bin_curve"]]), 0)
  expect_length(out, 6)
})

test_that("resolve_sel_pen_wts rejects unnamed or misspelled terms", {
  expect_error(SPoRC:::resolve_sel_pen_wts(c(1, 2)), "named numeric")
  expect_error(SPoRC:::resolve_sel_pen_wts(c(smooth_dom = 1)), "named numeric")
})

test_that("safe_extract returns zero for absent or NULL fields", {
  obj <- list(a = 5, b = NULL)
  expect_equal(SPoRC:::safe_extract(obj, "a"), 5)
  expect_equal(SPoRC:::safe_extract(obj, "b"), 0)
  expect_equal(SPoRC:::safe_extract(obj, "missing"), 0)
})
