library(SPoRC)
library(testthat)

# Selectivity supplied as a fixed input array and defined over length bins, for the total
# fishery, retention, and the survey at once. Estimated length-based selectivity is covered
# elsewhere, and fixed age-based selectivity is covered by the simulation tests, but the
# combination of the two is not: no other test sets use_fixed_*_sel = 1 alongside a length
# based selex type, so the branch that reads sel_input into the length-based array never
# runs. That branch previously indexed a seven dimensional array with six subscripts and
# would have errored on the first call.

n_lens_test <- 8

build <- function(...) suppressWarnings(suppressMessages(objective_setup_input(...)))

# A separate curve per class, and a shift per year within each class, so a slice taken from
# the wrong class or the wrong year cannot match the expected values by coincidence. This is
# the property worth defending here: the three classes share one code path, so a mix-up
# between them would otherwise be silent.
fixed_len_sel <- function(d, offset) {
  arr <- array(0, dim = c(d$n_pop, d$n_regions, length(d$years), d$n_seas,
                          length(d$lens), d$n_sexes, 1))
  for(y in seq_along(d$years)) {
    arr[1, 1, y, 1, , 1, 1] <- stats::plogis(seq_along(d$lens) - offset - y / length(d$years))
  }
  arr
}

# The default test setup is built once to read its dimensions, then rebuilt with the fixed
# length-based inputs. The underlying simulation is cached, so the second build is cheap.
fixed_length_input <- function() {
  d <- build(n_lens = n_lens_test)$data

  fish_in <- fixed_len_sel(d, offset = 2)
  ret_in <- fixed_len_sel(d, offset = 5)
  srv_in <- fixed_len_sel(d, offset = 3.5)

  input <- build(
    n_lens = n_lens_test,
    fishsel = list(
      fish_selex_type = "length",
      use_fixed_fish_sel = 1,
      fish_fixed_sel_pars_spec = "fix_fish_sel_input",
      fish_sel_input = fish_in,
      ret_selex_type = "length",
      use_fixed_ret_sel = 1,
      ret_fixed_sel_pars_spec = "fix_ret_sel_input",
      ret_sel_input = ret_in
    ),
    srvsel = list(
      srv_selex_type = "length",
      use_fixed_srv_sel = 1,
      srv_fixed_sel_pars_spec = "fix_srv_sel_input",
      srv_sel_input = srv_in
    )
  )

  list(input = input, fish = fish_in, ret = ret_in, srv = srv_in)
}


test_that("fixed length based selectivity is read into the length based arrays unchanged", {
  fx <- fixed_length_input()
  model <- evaluate_input(fx$input)

  expect_equal(fx$input$data$fish_selex_type, 1)
  expect_equal(fx$input$data$ret_selex_type, 1)
  expect_equal(fx$input$data$srv_selex_type, 1)

  years <- seq_along(fx$input$data$years)

  for(y in years) {
    expect_equal(as.vector(model$rep$fish_sel_l[1, y, , 1, 1]),
                 as.vector(fx$fish[1, 1, y, 1, , 1, 1]), tolerance = 1e-12,
                 info = paste("fishery, year", y))
    expect_equal(as.vector(model$rep$ret_sel_l[1, y, , 1, 1]),
                 as.vector(fx$ret[1, 1, y, 1, , 1, 1]), tolerance = 1e-12,
                 info = paste("retention, year", y))
    expect_equal(as.vector(model$rep$srv_sel_l[1, y, , 1, 1]),
                 as.vector(fx$srv[1, 1, y, 1, , 1, 1]), tolerance = 1e-12,
                 info = paste("survey, year", y))
  }
})


test_that("each selectivity class reads its own fixed input rather than another class's", {
  fx <- fixed_length_input()
  model <- evaluate_input(fx$input)

  fish_l <- model$rep$fish_sel_l
  ret_l <- model$rep$ret_sel_l
  srv_l <- model$rep$srv_sel_l

  expect_false(isTRUE(all.equal(fish_l, ret_l)))
  expect_false(isTRUE(all.equal(fish_l, srv_l)))
  expect_false(isTRUE(all.equal(ret_l, srv_l)))

  # the curves also have to differ across years, or the year index could be wrong and the
  # per-year checks above would still pass
  expect_gt(stats::sd(apply(fish_l[1, , , 1, 1], 1, mean)), 0)
})


test_that("fixed length based selectivity maps onto ages through the size-age transition", {
  fx <- fixed_length_input()
  model <- evaluate_input(fx$input)

  size_age <- fx$input$data$SizeAgeTrans

  for(y in seq_along(fx$input$data$years)) {
    expect_equal(as.vector(model$rep$fish_sel[1, 1, y, 1, , 1, 1]),
                 as.vector(model$rep$fish_sel_l[1, y, , 1, 1] %*% size_age[1, 1, y, 1, , , 1]),
                 tolerance = 1e-10, info = paste("fishery, year", y))
    expect_equal(as.vector(model$rep$ret_sel[1, 1, y, 1, , 1, 1]),
                 as.vector(model$rep$ret_sel_l[1, y, , 1, 1] %*% size_age[1, 1, y, 1, , , 1]),
                 tolerance = 1e-10, info = paste("retention, year", y))
    expect_equal(as.vector(model$rep$srv_sel[1, 1, y, 1, , 1, 1]),
                 as.vector(model$rep$srv_sel_l[1, y, , 1, 1] %*% size_age[1, 1, y, 1, , , 1]),
                 tolerance = 1e-10, info = paste("survey, year", y))
  }

  # selectivity at age varies, so the mapping is doing work rather than returning a flat
  # curve that any transition matrix would reproduce
  expect_gt(stats::sd(model$rep$fish_sel[1, 1, 1, 1, , 1, 1]), 0)

  # a weighted average of a bounded curve stays bounded
  for(quant_name in c("fish_sel", "ret_sel", "srv_sel")) {
    expect_true(all(model$rep[[quant_name]] >= 0), info = quant_name)
    expect_true(all(model$rep[[quant_name]] <= 1), info = quant_name)
  }
})


test_that("fixed length based selectivity gives a finite jnLL that still decomposes", {
  model <- evaluate_input(fixed_length_input()$input)

  expect_true(is.finite(as.numeric(model$rep$jnLL)))
  expect_jnLL_decomposes(model, label = "fixed length based selectivity")
})
