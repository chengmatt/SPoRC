library(SPoRC)
library(testthat)

# ── do_dmr_dev_mapping keys on fished cells, and matches get_dmr_penalty ─────
#
# dmr enters the likelihood through total mortality: ZAA has disc_FAA, and
# CAA = ret_FAA/ZAA * NAA * (1 - exp(-ZAA)), so dmr moves predicted retained
# catch, indices, and compositions in any fished cell where retention is less
# than one. Discard observations are not the boundary of identifiability -- dmr
# even cancels out of PredDiscard, which divides DAA back through by it.
#
# So deviations are estimated wherever the objective computes a non-zero dmr,
# i.e. every cell that is not a true closure, and get_dmr_penalty() must key on
# that same set so that every estimated deviation is penalized and no deviation
# fixed at zero contributes to the penalty.

n_pop <- 2; n_regions <- 2; n_yrs <- 6; n_seas <- 2; n_fleets <- 2
agg_dims <- c(n_regions, n_yrs, n_seas, n_fleets)
pop_dims <- c(n_pop, n_regions, n_yrs, n_seas, n_fleets)

make_catch_input_list <- function() {
  list(
    data = list(
      n_pop = n_pop,
      n_regions = n_regions,
      years = 1:n_yrs,
      n_seas = n_seas,
      n_fish_fleets = n_fleets
    ),
    par = list(),
    map = list(),
    verbose = FALSE,
    store_config = FALSE
  )
}

run_setup <- function(
  UseCatch,
  UseDiscard = array(0, dim = agg_dims),
  UseCatch_pop = array(0, dim = pop_dims),
  UseDiscard_pop = array(0, dim = pop_dims),
  ObsCatch = array(100, dim = agg_dims),
  dmr_dev_spec = "est_all",
  Use_dmr_pen = 1
) {
  suppressWarnings(Setup_Mod_Catch_and_F(
    input_list = make_catch_input_list(),
    ObsCatch = ObsCatch,
    UseCatch = UseCatch,
    ObsCatch_pop = array(50, dim = pop_dims),
    UseCatch_pop = UseCatch_pop,
    ObsDiscard = array(0.1, dim = agg_dims),
    UseDiscard = UseDiscard,
    ObsDiscard_pop = array(0.05, dim = pop_dims),
    UseDiscard_pop = UseDiscard_pop,
    Use_dmr_pen = Use_dmr_pen,
    dmr_dev_spec = dmr_dev_spec
  ))
}

# TRUE where a deviation is estimated (map index not NA)
estimated_cells <- function(il) array(!is.na(il$map$logit_dmr_devs), dim = agg_dims)

# TRUE where the penalty contributes, evaluated at a non-zero deviation so any
# penalized cell is non-zero. The penalty reads the map mirror from the data
# list, which is what keeps it in step with what is actually estimated
penalized_cells <- function(il) {
  pen <- SPoRC:::get_dmr_penalty(
    logit_dmr_devs = array(0.5, dim = agg_dims),
    ln_sigma_dmr = array(log(1), dim = c(n_regions, n_seas, n_fleets)),
    map_logit_dmr_devs = il$data$map_logit_dmr_devs,
    n_fish_fleets = n_fleets,
    n_yrs = n_yrs,
    n_regions = n_regions,
    n_seas = n_seas
  )
  pen != 0
}

test_that("deviations are estimated in fished cells and fixed in true closures", {
  UseCatch <- array(0, dim = agg_dims)
  UseCatch[, 2:5, , ] <- 1 # fished in the middle four years only

  il <- run_setup(UseCatch)
  est <- estimated_cells(il)

  expect_true(all(est[, 2:5, , ]))
  expect_false(any(est[, c(1, 6), , ])) # true closures: dmr is pinned to 0 there
})

test_that("a fished cell gets a deviation whether or not discard is observed", {
  UseCatch <- array(1, dim = agg_dims)
  UseDiscard <- array(0, dim = agg_dims)
  UseDiscard[, 4:6, , ] <- 1 # discard observed only in the last three years

  il <- run_setup(UseCatch, UseDiscard = UseDiscard)

  # dmr is informed through ZAA in years 1-3 too, so those deviations are
  # estimated -- and penalized, per the matching-set test below
  expect_true(all(estimated_cells(il)))
})

test_that("a true closure gets no deviation even where discard is observed", {
  UseCatch <- array(1, dim = agg_dims)
  UseCatch[1, 3, 1, 1] <- 0 # closed cell, recorded zero catch
  ObsCatch <- array(100, dim = agg_dims)
  ObsCatch[1, 3, 1, 1] <- 0
  UseDiscard <- array(1, dim = agg_dims) # discard nominally observed everywhere

  il <- run_setup(UseCatch, UseDiscard = UseDiscard, ObsCatch = ObsCatch)
  est <- estimated_cells(il)

  expect_false(est[1, 3, 1, 1]) # dmr is 0 in the objective here, so no deviation
  expect_equal(sum(!est), 1)
})

test_that("a missing (NA) catch observation activates a deviation", {
  UseCatch <- array(1, dim = agg_dims)
  UseCatch[2, 4, 2, 1] <- 0 # not fit
  ObsCatch <- array(100, dim = agg_dims)
  ObsCatch[2, 4, 2, 1] <- NA # but missing rather than a recorded zero

  il <- run_setup(UseCatch, ObsCatch = ObsCatch)

  # fishing is assumed to have continued through a missing observation, so dmr
  # is non-zero there and the deviation is estimated
  expect_true(all(estimated_cells(il)))
})

test_that("population-specific catch alone activates a deviation", {
  UseCatch <- array(0, dim = agg_dims)
  UseCatch_pop <- array(0, dim = pop_dims)
  UseCatch_pop[2, 1, 5, 1, 2] <- 1

  il <- run_setup(UseCatch, UseCatch_pop = UseCatch_pop)
  est <- estimated_cells(il)

  expect_true(est[1, 5, 1, 2])
  expect_equal(sum(est), 1)
})

test_that("the estimated set matches the penalized set exactly", {
  set.seed(11)
  UseCatch <- array(rbinom(prod(agg_dims), 1, 0.7), dim = agg_dims)
  UseCatch_pop <- array(rbinom(prod(pop_dims), 1, 0.15), dim = pop_dims)
  UseDiscard <- array(rbinom(prod(agg_dims), 1, 0.4), dim = agg_dims)
  ObsCatch <- array(100, dim = agg_dims)
  ObsCatch[1, 2, 1, 1] <- NA # exercise the missing-observation branch on both sides

  # the discard indicator disagrees with the fished set in both
  # directions, so the match below is a property of the fished set, not a
  # coincidence of the test setup
  fished <- UseCatch == 1 | apply(UseCatch_pop == 1, c(2,3,4,5), any) | is.na(ObsCatch)
  expect_true(any(fished & UseDiscard == 0))
  expect_true(any(!fished & UseDiscard == 1))

  il <- run_setup(UseCatch, UseDiscard = UseDiscard, UseCatch_pop = UseCatch_pop, ObsCatch = ObsCatch)
  expect_equal(estimated_cells(il), penalized_cells(il))
})

test_that("a deviation mapped off by hand is dropped from the penalty", {
  il <- run_setup(array(1, dim = agg_dims))
  expect_true(all(penalized_cells(il))) # every cell penalized to begin with

  # a user turns deviations off for years 3:5 after setup, then builds the model
  mp <- array(il$map$logit_dmr_devs, dim = agg_dims)
  mp[, 3:5, , ] <- NA
  il$map$logit_dmr_devs <- factor(mp)

  il$data <- SPoRC:::sync_dev_map_data(il$data, il$map)
  pen <- penalized_cells(il)

  expect_false(any(pen[, 3:5, , ])) # no penalty on the deviations they pinned
  expect_true(all(pen[, c(1, 2, 6), , ])) # and the rest are untouched
  expect_equal(pen, estimated_cells(il))
})

test_that("each estimated deviation gets its own index", {
  UseCatch <- array(0, dim = agg_dims)
  UseCatch[1, 2:4, , ] <- 1

  il <- run_setup(UseCatch)
  idx <- il$map$logit_dmr_devs[!is.na(il$map$logit_dmr_devs)]

  expect_equal(length(unique(idx)), sum(UseCatch)) # no unintended sharing
})

test_that("dmr_dev_spec = 'fix' maps every deviation to NA regardless of catch", {
  il <- run_setup(array(1, dim = agg_dims), dmr_dev_spec = "fix", Use_dmr_pen = 0)
  expect_true(all(is.na(il$map$logit_dmr_devs)))
  expect_equal(length(il$map$logit_dmr_devs), prod(agg_dims))
})
