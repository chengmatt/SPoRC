library(SPoRC)
library(testthat)

# Three hooks that let the projection be driven as a management strategy loop:
# recruitment deviations supplied from outside, a control rule that reads the
# population rather than spawning biomass alone, and catch at age handed back so
# a utility can be written on the size of the fish taken.
#
# The first two matter for differentiating through the loop. Deviations drawn
# outside and kept fixed are what makes an expectation over replicates a
# deterministic function of a rule's parameters; a rule that only sees one
# number cannot express a policy on more than one observation.

test_that("rec_devs of one leaves the projection unchanged", {

  base <- project_at_F(f = 0.06, n_proj_yrs = 8)
  ones <- project_at_F(f = 0.06, n_proj_yrs = 8,
                       rec_devs = array(1, dim = c(sgl_rg_sable_data$n_pop,
                                                   sgl_rg_sable_data$n_regions, 8)))

  expect_equal(base$proj_SSB, ones$proj_SSB)
  expect_equal(base$proj_Catch, ones$proj_Catch)
})


test_that("rec_devs scales the recruitment the option produced", {

  base <- project_at_F(f = 0.06, n_proj_yrs = 8)

  devs <- array(1, dim = c(sgl_rg_sable_data$n_pop, sgl_rg_sable_data$n_regions, 8))
  devs[,,4] <- 2
  doubled <- project_at_F(f = 0.06, n_proj_yrs = 8, rec_devs = devs)

  # recruits enter at age one, so year four's age-one numbers double and no
  # earlier year moves
  age1 <- function(o, y) sum(o$proj_NAA[,,y,1,1,])
  expect_equal(age1(doubled, 4), 2 * age1(base, 4))
  expect_equal(age1(doubled, 3), age1(base, 3))
})


test_that("rec_devs is checked for shape and sign", {

  n_pop <- sgl_rg_sable_data$n_pop
  n_regions <- sgl_rg_sable_data$n_regions

  expect_error(project_at_F(f = 0.06, n_proj_yrs = 8,
                            rec_devs = array(1, dim = c(n_pop, n_regions, 3))),
               "rec_devs should be dimensioned")
  expect_error(project_at_F(f = 0.06, n_proj_yrs = 8,
                            rec_devs = array(-1, dim = c(n_pop, n_regions, 8))),
               "negative or non-finite")
})


test_that("a control rule that asks for the state gets it, and one that does not is unaffected", {

  seen <- new.env()
  seen$n <- 0

  plain <- function(x, frp, brp) frp * min(1, x / brp)
  aware <- function(x, frp, brp, state) {
    seen$n <- seen$n + 1
    seen$names <- names(state)
    seen$ssb <- sum(state$SSB)
    seen$x <- x
    frp * min(1, x / brp)
  }

  a <- project_with_HCR(plain)
  b <- project_with_HCR(aware)

  expect_equal(a$proj_F, b$proj_F)
  expect_equal(a$proj_Catch, b$proj_Catch)
  expect_gt(seen$n, 0)
  expect_setequal(seen$names, c("y", "NAA", "SSB", "Total_Biom", "Catch", "r"))
  # the state's spawning biomass is the same quantity the rule is handed as x
  expect_equal(seen$ssb, seen$x)
})


test_that("catch at age is returned and weighs up to the catch biomass", {

  out <- project_at_F(f = 0.06, n_proj_yrs = 6)
  W <- sgl_rg_sable_data$WAA[1, 1, length(sgl_rg_sable_data$years), 1, , ]

  expect_false(is.null(out$proj_CAA))
  biom <- apply(out$proj_CAA[1,1,,1,,,], c(1, 4), function(m) sum(m * W))
  expect_equal(biom, out$proj_Catch[1,1,,1,], tolerance = 1e-10)
})
