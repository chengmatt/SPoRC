# Checks map_move_devs takes cells out of the movement penalty: one at a time under the iid forms
# (1e-12), the way a dsem takes a series over, and a whole block at a time under ar1.

move_pen_setup <- function(seed = 9, n_yrs = 8) {
  set.seed(seed)
  d <- c(1, 2, 2, n_yrs, 1, 3, 1) # pop, from, to, year, season, age, sex
  pe <- array(0, dim = c(2, 2, 3)) # from, to, [log sd, rho age, rho year]
  pe[,,1] <- log(0.25)
  pe[,,3] <- 0.7
  list(d = d, n_yrs = n_yrs,
       move_devs = array(rnorm(prod(d), 0, 0.3), dim = d),
       PE_pars = pe,
       pairs = cbind(region_from = rep(1:2, each = 2), region_to = rep(1:2, 2)))
}

# a year series per pair, every other dim shared and read at its first level
move_pen <- function(s, mirror, year_re = 1) Get_move_PE_loglik(move_year_re = year_re, move_age_re = 0, move_pop_re = 0, move_seas_re = 0, move_sex_re = 0,
                                                                PE_pars = s$PE_pars, move_pop_corr_pars = 0, move_seas_corr_pars = 0, move_sex_corr_pars = 0,
                                                                move_devs = s$move_devs, map_move_devs = mirror, move_pairs = s$pairs,
                                                                move_pe_block = 1:4, move_pop_block = 1, move_year_block = 1:s$n_yrs,
                                                                move_seas_block = 1, move_age_block = rep(1, 3), move_sex_block = 1, move_dsem = 0)

mvn_logpdf <- function(x, S) {
  L <- chol(S)
  z <- backsolve(L, x, transpose = TRUE)
  -0.5 * (length(x) * log(2 * pi) + 2 * sum(log(diag(L))) + sum(z^2))
}

test_that("a mirror keeping every cell matches the penalty by hand", {

  s <- move_pen_setup()
  full <- array(1, dim = s$d)
  # a year only form reads pop 1, season 1, age 1, sex 1 across from, to and year
  by_hand <- sum(dnorm(as.vector(s$move_devs[1,,,,1,1,1]), 0, 0.25, log = TRUE))
  expect_equal(move_pen(s, full), by_hand, tolerance = 1e-12)

})

test_that("blanked cells drop exactly their own contribution under the iid form", {

  s <- move_pen_setup()
  full <- array(1, dim = s$d)
  mirror <- full
  mirror[1,1,2,3:5,1,1,1] <- NA # region 1 to region 2, years 3 to 5, the cells the penalty reads

  dropped <- sum(dnorm(s$move_devs[1,1,2,3:5,1,1,1], 0, 0.25, log = TRUE))
  expect_equal(move_pen(s, mirror), move_pen(s, full) - dropped, tolerance = 1e-12)

  # a blank on a cell the penalty never reads (age 2 under a year-only model) changes nothing
  elsewhere <- full; elsewhere[1,1,2,3,1,2,1] <- NA
  expect_equal(move_pen(s, elsewhere), move_pen(s, full), tolerance = 1e-12)

  expect_equal(move_pen(s, array(NA, dim = s$d)), 0, tolerance = 1e-12) # everything out

})

test_that("an ar1 block leaves whole or not at all", {

  s <- move_pen_setup()
  full <- array(1, dim = s$d)
  rho <- SPoRC:::rho_trans(0.7)
  S <- 0.25^2 / (1 - rho^2) * rho^abs(outer(1:s$n_yrs, 1:s$n_yrs, "-")) # the parameter is the conditional sd
  by_hand <- 0
  for(r in 1:2) for(rr in 1:2) by_hand <- by_hand + mvn_logpdf(s$move_devs[1,r,rr,,1,1,1], S)
  expect_equal(move_pen(s, full, year_re = 2), by_hand, tolerance = 1e-12)

  # a block blanked entirely drops its own density
  mirror <- full; mirror[1,1,2,,1,,1] <- NA
  expect_equal(move_pen(s, mirror, year_re = 2), by_hand - mvn_logpdf(s$move_devs[1,1,2,,1,1,1], S), tolerance = 1e-12)

  # a block blanked in part cannot leave, since its cells are correlated
  partial <- full; partial[1,1,2,3:5,1,1,1] <- NA
  expect_error(move_pen(s, partial, year_re = 2), "cannot be dropped one at a time")

})
