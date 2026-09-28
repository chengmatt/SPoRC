# The same population, described at two resolutions.
#
# A model split across identical regions, sexes, seasons or fleets describes exactly the population a
# model with one of each describes. Where they disagree, an index is walking a dim it should not.
#
# A stored value cannot make this check, since a read off the wrong dimension still returns a
# stable number forever. No two dimensions are the same size, since equal ones hide a transpose.


test_that("regions that mix completely describe one region", {
  # under full mixing with the same biology everywhere, a fish's fate does not depend on which
  # region it is in, so the summed population is the single-region one
  one <- collapse_rep(nr = 1)

  for(nr in c(2, 3)) {
    expect_collapses(one, collapse_rep(nr = nr), sprintf("%d regions vs 1", nr))
  }
})


test_that("sexes with identical biology describe one sex", {
  # weight, maturity, mortality and selectivity are the same for both sexes here, so splitting
  # by sex changes how the population is stored and nothing about how it develops
  expect_collapses(collapse_rep(nx = 1), collapse_rep(nx = 2), "2 sexes vs 1")
})


test_that("fleets sharing a selectivity describe one fleet", {
  # two fleets each taking half the catch at half the fishing mortality remove exactly what one
  # fleet taking all of it removes. the halving is explicit, each fleet's F being its start
  one <- collapse_rep(nf = 1)

  for(nf in c(2, 5)) {
    expect_collapses(one, collapse_rep(nf = nf, f_scale = 1 / nf),
                     sprintf("%d fleets vs 1", nf))
  }
})


test_that("the collapse test setup is actually sensitive to the dynamics", {
  # changing the fishing mortality has to move the very quantities the collapse tests compare,
  # or those tests would hold because both sides are trivially equal
  base <- collapse_rep(nr = 1)
  harder <- collapse_rep(nr = 1, f_scale = 4)

  moved <- vapply(c("NAA", "SSB", "Total_Biom"), function(quant_name)
    max(abs(apply(base[[quant_name]], 3, sum) - apply(harder[[quant_name]], 3, sum))) /
      max(abs(apply(base[[quant_name]], 3, sum))), numeric(1))

  expect_true(all(moved > 0.05))
})


test_that("splitting a region does not change what is predicted for the fishery", {
  # the population collapsing is one thing, the observations reading it on the right dimensions
  # another, and predicted catch at age is where the two meet
  one <- collapse_rep(nr = 1)

  for(nr in c(2, 3)) {
    fine <- collapse_rep(nr = nr)
    expect_collapses(one, fine, sprintf("predicted catch, %d regions vs 1", nr),
                     what = "CAA")
  }
})


test_that("seasons that share the year's fishing describe one season", {
  # a year cut into k seasons, each taking a kth of the fishing mortality, removes over the year
  # exactly what one season taking all of it removes.
  #
  # numbers at age are recorded within each season, so the same fish appear k times over the
  # season dimension and the comparison is made in the first, where both models agree
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
