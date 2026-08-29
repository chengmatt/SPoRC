# The same population, described at two resolutions.
#
# A model split across regions, sexes, seasons or fleets that are all identical
# describes exactly the population a model with one of each describes. Where the
# two disagree, an index is walking a margin it should not, or a quantity is
# being summed on the wrong side of an operation.
#
# This is the check a stored value cannot make. A stride that reads the wrong
# margin still returns a stable number, and a regression test pinned to that
# number passes forever. The relations below are true of the model regardless of
# what it currently computes, so nothing about them can be pinned wrong.
#
# The fixture dimensions are deliberately unequal (see helper-sweep_fixture.R):
# when two dimensions have the same extent, a transposition between them is
# invisible.


test_that("regions that mix completely describe one region", {
  # Under full mixing with identical biology in every region, a fish's fate does
  # not depend on which region it is in, so the summed population is the
  # single-region population.
  one <- collapse_rep(nr = 1)

  for(nr in c(2, 3)) {
    expect_collapses(one, collapse_rep(nr = nr), sprintf("%d regions vs 1", nr))
  }
})


test_that("sexes with identical biology describe one sex", {
  # Weight, maturity, mortality and selectivity are the same for both sexes here,
  # so splitting the population by sex changes how it is stored and nothing about
  # how it develops.
  expect_collapses(collapse_rep(nx = 1), collapse_rep(nx = 2), "2 sexes vs 1")
})


test_that("fleets sharing a selectivity describe one fleet", {
  # Two fleets each taking half the catch at half the fishing mortality remove
  # exactly what one fleet taking all of it removes. The halving is explicit
  # because these models are evaluated rather than fitted, so each fleet's F
  # comes from its starting value rather than from the catch it is given.
  one <- collapse_rep(nf = 1)

  for(nf in c(2, 5)) {
    expect_collapses(one, collapse_rep(nf = nf, f_scale = 1 / nf),
                     sprintf("%d fleets vs 1", nf))
  }
})


test_that("the collapse fixture is actually sensitive to the dynamics", {
  # A relation that holds because both sides are trivially equal proves nothing.
  # Changing the fishing mortality has to move the very quantities the collapse
  # tests compare, or those tests would pass against a broken model.
  base <- collapse_rep(nr = 1)
  harder <- collapse_rep(nr = 1, f_scale = 4)

  moved <- vapply(c("NAA", "SSB", "Total_Biom"), function(nm)
    max(abs(apply(base[[nm]], 3, sum) - apply(harder[[nm]], 3, sum))) /
      max(abs(apply(base[[nm]], 3, sum))), numeric(1))

  expect_true(all(moved > 0.05))
})


test_that("splitting a region does not change what is predicted for the fishery", {
  # The population collapsing is one claim; the observation layer reading that
  # population on the right margins is another. Predicted catch at age is where
  # the two meet, so it is compared across the same splits.
  one <- collapse_rep(nr = 1)

  for(nr in c(2, 3)) {
    fine <- collapse_rep(nr = nr)
    expect_collapses(one, fine, sprintf("predicted catch, %d regions vs 1", nr),
                     what = "CAA")
  }
})


test_that("seasons that share the year's fishing describe one season", {
  # A year cut into k seasons, each taking 1/k of the fishing mortality, removes
  # over the year exactly what a single season taking all of it removes.
  #
  # Numbers at age carry a season margin and are recorded within each season, so
  # the same fish appear once per season and summing over that margin counts them
  # k times. The comparison is made in the first season, where both models are at
  # the same point in the year. The annual quantities have no season margin and
  # are compared whole.
  one <- collapse_rep(ns = 1)

  for(ns in c(2, 3)) {
    fine <- collapse_rep(ns = ns, f_scale = 1 / ns)
    expect_collapses(one, fine, sprintf("%d seasons vs 1", ns),
                     what = c("SSB", "Total_Biom", "Rec"))

    season1 <- function(r) apply(r$NAA[, , , 1, , , drop = FALSE], 3, sum)
    expect_equal(season1(fine), season1(one), tolerance = 1e-10,
                 label = sprintf("%d seasons vs 1: NAA in season 1", ns))
  }
})
