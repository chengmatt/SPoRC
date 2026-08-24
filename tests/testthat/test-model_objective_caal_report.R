library(SPoRC)
library(testthat)

# Joint arrays at length and age, reported when do_caal == 1. Nothing is fit to them, so
# the tests here are about the arrays being the joint distribution they claim to be: their
# length margin has to be the catch or index at age, and their age margin the catch or
# index at length. The flag also has to leave the objective untouched, since a reporting
# switch that moved the likelihood would change every model that turns it on.
#
# The arrays are dimensioned like CAA and CAL, that is population specific, and no region
# level version is reported. Summing the population margin away is left to the caller, the
# same way the composition likelihoods already do it on CAA and CAL.

n_lens_test <- 12

build <- function(...) suppressWarnings(suppressMessages(objective_fixture_input(...)))
run <- function(input) suppressWarnings(suppressMessages(evaluate_input(input)))

with_caal <- function() {
  input <- build(n_lens = n_lens_test)
  input$data$do_caal <- 1
  input
}

# the joint arrays and the marginal arrays they have to agree with, by fleet type
caal_pairs <- list(
  list(joint = "Fish_caal",    at_age = "CAA",    at_len = "CAL"),
  list(joint = "Fish_caal_discard",    at_age = "DAA",    at_len = "DAL"),
  list(joint = "Srv_caal", at_age = "SrvIAA", at_len = "SrvIAL")
)


test_that("the joint arrays are reported only when asked for", {
  off <- run(build(n_lens = n_lens_test))
  expect_equal(off$data$do_caal, 0)
  for(pair in caal_pairs) expect_null(off$rep[[pair$joint]])

  on <- run(with_caal())
  for(pair in caal_pairs) expect_false(is.null(on$rep[[pair$joint]]))
})


test_that("the joint arrays are dimensioned like the marginals they extend", {
  model <- run(with_caal())
  for(pair in caal_pairs) {
    joint <- model$rep[[pair$joint]]
    at_age <- model$rep[[pair$at_age]]
    at_len <- model$rep[[pair$at_len]]
    # population, region, year and season come first, exactly as in CAA and CAL
    expect_equal(dim(joint)[1:4], dim(at_age)[1:4], info = pair$joint)
    # then a length and an age dimension at once, then sex and fleet
    expect_equal(dim(joint)[5], dim(at_len)[5], info = pair$joint)
    expect_equal(dim(joint)[6], dim(at_age)[5], info = pair$joint)
    expect_equal(dim(joint)[7:8], dim(at_age)[6:7], info = pair$joint)
  } # end pair loop
})


test_that("turning the joint arrays on leaves the objective unchanged", {
  off <- run(build(n_lens = n_lens_test))
  on <- run(with_caal())
  expect_equal(on$fn(on$par), off$fn(off$par), tolerance = 0)
})


test_that("the length margin of each joint array is the marginal at age", {
  model <- run(with_caal())
  for(pair in caal_pairs) {
    joint <- model$rep[[pair$joint]]
    expect_equal(apply(joint, c(1,2,3,4,6,7,8), sum), model$rep[[pair$at_age]],
                 tolerance = 1e-12, info = pair$joint)
  } # end pair loop
})


test_that("the age margin of each joint array is the marginal at length", {
  model <- run(with_caal())
  for(pair in caal_pairs) {
    joint <- model$rep[[pair$joint]]
    expect_equal(apply(joint, c(1,2,3,4,5,7,8), sum), model$rep[[pair$at_len]],
                 tolerance = 1e-12, info = pair$joint)
  } # end pair loop
})


test_that("each joint cell is the size-age transition scaled by the marginal at age", {
  input <- with_caal()
  model <- run(input)

  size_age <- input$data$SizeAgeTrans
  n_ages <- length(input$data$ages)

  # SizeAgeTrans holds P(len | age), so scaling age column a by the catch at age a is the
  # joint array by definition. Checked cell by cell on one stratum rather than through the
  # same vectorized expression the model uses.
  p <- 1; r <- 1; y <- 7; seas <- 1; s <- 1; f <- 1
  expected <- matrix(0, nrow = n_lens_test, ncol = n_ages)
  for(a in 1:n_ages) expected[,a] <- size_age[p,r,y,seas,,a,s] * model$rep$CAA[p,r,y,seas,a,s,f]

  expect_equal(model$rep$Fish_caal[p,r,y,seas,,,s,f], expected, tolerance = 1e-12)
})


test_that("conditioning a row on its own row sum gives the age composition at that length", {
  model <- run(with_caal())

  p <- 1; r <- 1; y <- 7; seas <- 1; l <- 5; s <- 1; f <- 1
  row <- model$rep$Fish_caal[p,r,y,seas,l,,s,f]
  expect_equal(sum(row), model$rep$CAL[p,r,y,seas,l,s,f], tolerance = 1e-12) # denominator is CAL
  expect_equal(sum(row / sum(row)), 1, tolerance = 1e-12)
})


test_that("do_caal is rejected without length compositions", {
  input <- build()
  expect_equal(input$data$fit_lengths, 0)
  expect_error(
    suppressMessages(Setup_Mod_Biologicals(
      input_list = input, WAA = input$data$WAA, MatAA = input$data$MatAA,
      fit_lengths = 0, SizeAgeTrans = NA, do_caal = 1,
      M_spec = "fix", Fixed_natmort = input$data$Fixed_natmort
    )),
    "do_caal"
  )
})
