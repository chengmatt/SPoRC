# Natural mortality and NAA process error blocks: a cell reads ln_M[popblk, regionblk, yearblk, seasblk,
# ageblk, sexblk] for the blocks it sits in, the M prior reads the same slot, and a retrospective peel keeps
# the arrays whole and fixes the blocks it leaves with no cell.

library(SPoRC)
library(testthat)


block_dims <- list(n_regions = 1, n_sexes = 2, n_fish_fleets = 1, n_seas = 1, n_yrs = 22, n_ages = 8)

# which block of a spec list holds position i
block_of <- function(spec, i) which(vapply(spec, function(v) i %in% v, logical(1)))

test_that("a cell reads the ln_M and ln_sigmaNAA slot of its own age and sex blocks", {

  age_spec <- list(1:2, 3:5, 6:8)
  sex_spec <- list(1, 2)
  il <- sweep_input(dims = block_dims, biol = list(M_spec = "est_ln_M", M_ageblk_spec = age_spec, M_sexblk_spec = sex_spec,
                                                   NAA_re = "iid", NAA_sigma_ageblk_spec = age_spec,
                                                   NAA_sigma_sexblk_spec = sex_spec))
  expect_equal(dim(il$par$ln_M), c(1, 1, 1, 1, 3, 2))
  expect_equal(dim(il$par$ln_sigmaNAA), c(1, 1, 1, 1, 3, 2))

  # with two blocked dims, the stored order and the loop order differ, so check every age and sex
  for(a in 1:8) for(s in 1:2) {
    want <- c(1, 1, 1, 1, block_of(age_spec, a), block_of(sex_spec, s))
    expect_equal(as.vector(arrayInd(il$data$M_blocks[1, 1, 5, 1, a, s], dim(il$par$ln_M))), want)
    expect_equal(as.vector(arrayInd(il$data$naa_sigma_blocks[1, 1, 5, 1, a, s], dim(il$par$ln_sigmaNAA))), want)
  } # end a, s loop
})

test_that("the M prior reads the block it names, and refuses a block the spec does not make", {

  # ages 3 to 5 are the second block, so a prior on age block 2 must land on ln_M[,,,,2,]
  ln_M <- array(log(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)), dim = c(1, 1, 1, 1, 3, 2))
  prior <- data.frame(popblk = 1, regionblk = 1, yearblk = 1, ageblk = 2, sexblk = 2, mu = 0.5, sd = 0.1)
  expect_equal(SPoRC:::get_natmort_prior(prior, ln_M, M_blocks = NULL),
               -dnorm(ln_M[1, 1, 1, 1, 2, 2], log(0.5), 0.1, log = TRUE))

  # an age rather than a block number is refused by name
  expect_error(sweep_input(dims = block_dims, biol = list(M_spec = "est_ln_M", M_ageblk_spec = list(1:2, 3:5, 6:8),
                                                          Use_M_prior = 1, M_prior = transform(prior, ageblk = 6))),
               "M_prior\\$ageblk has 6, but M_ageblk_spec makes 3 blocks")
})

test_that("a retrospective peel keeps the block arrays whole and fixes the blocks it empties", {

  il <- sweep_input(dims = block_dims, biol = list(M_spec = "est_ln_M", M_yearblk_spec = list(1:16, 17:22),
                                                   M_ageblk_spec = list(1:4, 5:8), NAA_re = "iid",
                                                   NAA_sigma_yearblk_spec = list(1:16, 17:22),
                                                   NAA_sigma_ageblk_spec = list(1:4, 5:8)))
  peel <- SPoRC:::truncate_yr(6, il$data, il$par, il$map)

  # the second year block covers 2007 to 2012 only, so peeling six years leaves it with no cell
  expect_equal(dim(peel$retro_parameters$ln_M), dim(il$par$ln_M))
  map_M <- array(peel$retro_mapping$ln_M, dim = dim(il$par$ln_M))
  expect_true(all(is.na(map_M[1, 1, 2, 1, , ])))
  expect_false(anyNA(map_M[1, 1, 1, 1, , ]))

  expect_equal(dim(peel$retro_parameters$ln_sigmaNAA), dim(il$par$ln_sigmaNAA))
  map_sigma <- array(peel$retro_mapping$ln_sigmaNAA, dim = dim(il$par$ln_sigmaNAA))
  expect_true(all(is.na(map_sigma[1, 1, 2, 1, , ])))
  expect_false(anyNA(map_sigma[1, 1, 1, 1, , ]))
})
