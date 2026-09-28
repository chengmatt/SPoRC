# Reference points against what they are defined to be, rather than against a stored baseline that came
# from the same code and so cannot catch one that was wrong when it was taken.
#
# An SPR reference point is the F at which spawning biomass per recruit is a stated fraction of its unfished
# value, so computing that fraction back from b_ref_pt / virgin_b_ref_pt must return what was asked for.

data("sgl_rg_sable_rep")
data("sgl_rg_sable_data")

refpt_spr <- function(spr_x) {
  Get_Reference_Points(
    data = sgl_rg_sable_data,
    rep = sgl_rg_sable_rep,
    SPR_x = spr_x,
    type = "single_region",
    what = "SPR",
    calc_rec_st_yr = 20,
    rec_age = 2
  )
}

#' Spawning potential ratio actually achieved at a reference point
#'
#' @keywords internal
realized_spr <- function(ref) as.numeric(ref$b_ref_pt[1] / ref$virgin_b_ref_pt[1])


test_that("the F returned for a target SPR achieves that SPR", {
  # the definition read back to the solver, so a reference point that solved the wrong
  # equation or reported a different quantity fails whatever number it returns
  for(spr_x in c(0.05, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 0.95)) {
    ref <- refpt_spr(spr_x)
    expect_equal(realized_spr(ref), spr_x, tolerance = 1e-6,
                 label = sprintf("realized SPR at SPR_x = %g", spr_x))
  }
})


test_that("fishing harder leaves a smaller share of the unfished stock", {
  # spawning biomass per recruit falls with F, so a lower target needs a higher F. this is
  # what says the solver walks the curve the right way
  targets <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  f <- vapply(targets, function(x) as.numeric(refpt_spr(x)$f_ref_pt[1]), numeric(1))

  expect_true(all(diff(f) < 0),
              label = paste("F decreasing in SPR target; got", paste(signif(f, 4), collapse = " ")))
})


test_that("an SPR target of nearly one leaves the stock nearly unfished", {
  # the endpoint the curve has to pass through. an offset or a scaling error can still be
  # monotone and still hit targets in the middle of the range
  expect_lt(as.numeric(refpt_spr(0.99)$f_ref_pt[1]), 0.005)
  expect_gt(as.numeric(refpt_spr(0.05)$f_ref_pt[1]),
            as.numeric(refpt_spr(0.5)$f_ref_pt[1]))
})


test_that("the unfished reference does not depend on the target asked for", {
  # unfished spawning biomass belongs to the stock, not to the target, so if it moves with
  # SPR_x the unfished calculation is reading fishing mortality somewhere
  virgin <- vapply(c(0.2, 0.4, 0.6, 0.8),
                   function(x) as.numeric(refpt_spr(x)$virgin_b_ref_pt[1]), numeric(1))

  expect_equal(virgin, rep(virgin[1], length(virgin)), tolerance = 1e-10)
})


test_that("the biomass reference point scales with the target as the ratio says", {
  # b_ref_pt is the unfished value times the target, so the two returned numbers have to
  # agree with the target that produced them
  for(spr_x in c(0.2, 0.4, 0.6)) {
    ref <- refpt_spr(spr_x)
    expect_equal(as.numeric(ref$b_ref_pt[1]),
                 spr_x * as.numeric(ref$virgin_b_ref_pt[1]), tolerance = 1e-8,
                 label = sprintf("b_ref_pt at SPR_x = %g", spr_x))
  }
})


test_that("these checks would notice a reference point that ignored its target", {
  # every test above compares the solver against its own target, so asserting that F
  # actually varies is what stops them passing on a solver that returns a constant
  f <- vapply(c(0.2, 0.4, 0.6), function(x) as.numeric(refpt_spr(x)$f_ref_pt[1]), numeric(1))

  expect_gt(diff(range(f)) / max(f), 0.5)
})


# the checks above stay inside the reference point solver. the ones below put it against
# the projection, which reaches the same equilibrium by stepping the population forward

test_that("projecting at the reference F reaches the reference biomass", {
  # b_ref_pt says what spawning biomass a stock fished at f_ref_pt settles at, and the
  # projection reaches that number the long way
  for(spr_x in c(0.3, 0.4, 0.5)) {
    ref <- refpt_spr(spr_x)
    out <- project_at_F(as.numeric(ref$f_ref_pt[1]), n_proj_yrs = 400)

    expect_equal(equilibrium_ssb(out), as.numeric(ref$b_ref_pt[1]), tolerance = 1e-8,
                 label = sprintf("projected equilibrium SSB at F_SPR%g", spr_x * 100))
  }
})


test_that("an unfished projection reaches the virgin biomass", {
  # The same statement at F = 0, where the reference point is the unfished
  # spawning biomass and the projection should take no catch at all.
  ref <- refpt_spr(0.4)
  out <- project_at_F(0, n_proj_yrs = 400)

  expect_equal(equilibrium_ssb(out), as.numeric(ref$virgin_b_ref_pt[1]), tolerance = 1e-8)
  expect_equal(equilibrium_catch(out), 0, tolerance = 1e-12)
})


test_that("equilibrium biomass falls as fishing mortality rises", {
  # the equilibrium the projection settles at falls with F, which is the property the
  # solver inverts. read on the projection side, so both would have to be wrong together
  f <- c(0, 0.02, 0.05, 0.1, 0.2, 0.4)
  ssb <- vapply(f, function(x) equilibrium_ssb(project_at_F(x, n_proj_yrs = 250)), numeric(1))

  expect_true(all(diff(ssb) < 0),
              label = paste("equilibrium SSB decreasing in F; got", paste(signif(ssb, 4), collapse = " ")))
})


test_that("the projection has actually equilibrated where these tests read it", {
  # every comparison above reads one year of a long projection as the equilibrium, so if
  # the run were still moving the agreement would depend on the year chosen
  out <- project_at_F(0.0862541, n_proj_yrs = 400)
  last <- equilibrium_ssb(out)
  earlier <- proj_year_total(out$proj_SSB, offset = 50)

  expect_equal(last, earlier, tolerance = 1e-8)
})
