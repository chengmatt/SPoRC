library(SPoRC)
library(testthat)

# Retention selectivity estimated over length bins rather than age bins. No other test
# switches ret_selex_type to length, so the branch that picks the length bins and the one
# that maps the length-based curve onto ages through the size-age transition are otherwise
# never run.

n_lens_test <- 8

build <- function(...) suppressWarnings(suppressMessages(objective_setup_input(...)))

length_based <- function() build(n_lens = n_lens_test, fishsel = list(ret_selex_type = "length"))
age_based <- function() build(n_lens = n_lens_test)


test_that("length based retention selectivity is estimated over length bins", {
  input <- length_based()
  model <- evaluate_input(input)

  expect_equal(input$data$ret_selex_type, 1)

  # the length-based curve is reported only under this branch, and it is indexed by length
  expect_false(is.null(model$rep$ret_sel_l))
  expect_equal(dim(model$rep$ret_sel_l)[3], n_lens_test)

  # the same test setup with age-based retention reports no length-based retention curve
  age_model <- evaluate_input(age_based())
  expect_equal(age_model$data$ret_selex_type, 0)
  expect_null(age_model$rep$ret_sel_l)
})


test_that("length based retention is mapped onto ages through the size-age transition", {
  input <- length_based()
  model <- evaluate_input(input)

  size_age <- input$data$SizeAgeTrans
  ret_sel_l <- model$rep$ret_sel_l
  ret_sel <- model$rep$ret_sel

  for(y in seq_along(input$data$years)) {
    expected <- as.vector(ret_sel_l[1, y, , 1, 1] %*% size_age[1, 1, y, 1, , , 1])
    expect_equal(as.vector(ret_sel[1, 1, y, 1, , 1, 1]), expected, tolerance = 1e-10,
                 info = paste("year", y))
  }

  # retention varies across ages here, so the mapping is doing something rather than
  # returning a flat curve that would match trivially
  expect_gt(stats::sd(ret_sel[1, 1, 1, 1, , 1, 1]), 0)

  # every retention value stays a proportion, as a weighted average of a bounded curve must
  expect_true(all(ret_sel >= 0))
  expect_true(all(ret_sel <= 1))
})


test_that("length based retention gives a different jnLL from age based, and still decomposes", {
  model <- evaluate_input(length_based())
  age_model <- evaluate_input(age_based())

  expect_false(isTRUE(all.equal(as.numeric(model$rep$jnLL), as.numeric(age_model$rep$jnLL))))
  expect_jnLL_decomposes(model, label = "length based retention")
  expect_jnLL_decomposes(age_model, label = "age based retention with lengths fitted")
})
