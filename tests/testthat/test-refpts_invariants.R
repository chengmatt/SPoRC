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
  # The definition, stated back to the solver. A reference point that solved the
  # wrong equation, or reported a quantity other than the one it solved for,
  # fails here whatever number it happens to return.
  for(spr_x in c(0.05, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 0.95)) {
    ref <- refpt_spr(spr_x)
    expect_equal(realized_spr(ref), spr_x, tolerance = 1e-6,
                 label = sprintf("realized SPR at SPR_x = %g", spr_x))
  }
})


test_that("fishing harder leaves a smaller share of the unfished stock", {
  # SPR is decreasing in F, so the F that achieves a lower SPR target must be
  # higher. This is what says the solver is walking the curve in the right
  # direction rather than landing on the target from an arbitrary place.
  targets <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  f <- vapply(targets, function(x) as.numeric(refpt_spr(x)$f_ref_pt[1]), numeric(1))

  expect_true(all(diff(f) < 0),
              label = paste("F decreasing in SPR target; got", paste(signif(f, 4), collapse = " ")))
})


test_that("an SPR target of nearly one leaves the stock nearly unfished", {
  # The endpoint the curve has to pass through. A solver with an offset or a
  # scaling error can still be monotone and still hit its targets in the middle
  # of the range while missing here.
  expect_lt(as.numeric(refpt_spr(0.99)$f_ref_pt[1]), 0.005)
  expect_gt(as.numeric(refpt_spr(0.05)$f_ref_pt[1]),
            as.numeric(refpt_spr(0.5)$f_ref_pt[1]))
})


test_that("the unfished reference does not depend on the target asked for", {
  # Virgin spawning biomass is a property of the stock, not of the target. If it
  # moves with SPR_x then the unfished calculation is picking up the fished
  # mortality somewhere.
  virgin <- vapply(c(0.2, 0.4, 0.6, 0.8),
                   function(x) as.numeric(refpt_spr(x)$virgin_b_ref_pt[1]), numeric(1))

  expect_equal(virgin, rep(virgin[1], length(virgin)), tolerance = 1e-10)
})


test_that("the biomass reference point scales with the target as the ratio says", {
  # b_ref_pt is the unfished value times the target, so the two returned
  # quantities have to agree with the target that produced them. This is the same
  # statement as the first test read the other way round, and it fails separately
  # if b_ref_pt is computed from something other than the SPR that was solved.
  for(spr_x in c(0.2, 0.4, 0.6)) {
    ref <- refpt_spr(spr_x)
    expect_equal(as.numeric(ref$b_ref_pt[1]),
                 spr_x * as.numeric(ref$virgin_b_ref_pt[1]), tolerance = 1e-8,
                 label = sprintf("b_ref_pt at SPR_x = %g", spr_x))
  }
})


test_that("these checks would notice a reference point that ignored its target", {
  # Every test above compares the solver against its own target. If the solver
  # returned the same F whatever was asked, the monotonicity test would fail, but
  # the round trip would too only if the reported ratio moved with it. Asserting
  # the F actually varies keeps the suite from passing on a constant.
  f <- vapply(c(0.2, 0.4, 0.6), function(x) as.numeric(refpt_spr(x)$f_ref_pt[1]), numeric(1))

  expect_gt(diff(range(f)) / max(f), 0.5)
})


# The checks above are internal to the reference point solver. The ones below put
# it against the projection engine, which computes the same equilibrium by
# stepping the population forward rather than by solving for it. The two share no
# code, so agreement between them is evidence about both.

test_that("projecting at the reference F reaches the reference biomass", {
  # b_ref_pt says what spawning biomass a stock fished at f_ref_pt settles at.
  # Do_Population_Projection reaches that number the long way. A solver that
  # returns a consistent pair of its own quantities but describes a stock the
  # dynamics never produce fails here and nowhere else.
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
  # Monotonicity of the equilibrium the projection converges to, which is the
  # property the reference point solver is inverting. Checked on the projection
  # side so a solver and a dynamics that agreed on a non-monotone curve would
  # still be caught.
  f <- c(0, 0.02, 0.05, 0.1, 0.2, 0.4)
  ssb <- vapply(f, function(x) equilibrium_ssb(project_at_F(x, n_proj_yrs = 250)), numeric(1))

  expect_true(all(diff(ssb) < 0),
              label = paste("equilibrium SSB decreasing in F; got", paste(signif(ssb, 4), collapse = " ")))
})


test_that("the projection has actually equilibrated where these tests read it", {
  # Every comparison above reads one year of a long projection and calls it the
  # equilibrium. If the run were still moving there, the agreement would be a
  # coincidence of the year chosen.
  out <- project_at_F(0.0862541, n_proj_yrs = 400)
  last <- equilibrium_ssb(out)
  earlier <- proj_year_total(out$proj_SSB, offset = 50)

  expect_equal(last, earlier, tolerance = 1e-8)
})
