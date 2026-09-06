library(SPoRC)
library(testthat)
data("sgl_rg_sable_rep")
data("sgl_rg_sable_data")
data("mlt_rg_sable_rep")
data("mlt_rg_sable_data")

# MSY is the maximum of equilibrium yield over a stock-recruit curve. A fit with
# rec_model = "mean_rec" never estimated one, and the equilibrium recruitment
# helpers branch on Ricker against everything else, so rec_model 0 used to fall
# through to Beverton-Holt at the default steepness of 0.6 without complaint.
# These tests pin the error, and pin that SPR is unaffected.

test_that("MSY reference points error rather than defaulting to Beverton-Holt under mean recruitment", {

  mean_rec_data <- sgl_rg_sable_data
  mean_rec_data$rec_model <- 0 # mean_rec

  expect_error(
    Get_Reference_Points(data = mean_rec_data,
                         rep = sgl_rg_sable_rep,
                         SPR_x = 0.4,
                         type = 'single_region',
                         what = 'MSY'),
    "requires a stock-recruit curve"
  )

  # The message has to name the escape route, since the whole point is that the
  # user asked for something undefined and needs to be told what to ask for.
  expect_error(
    Get_Reference_Points(data = mean_rec_data,
                         rep = sgl_rg_sable_rep,
                         SPR_x = 0.4,
                         type = 'single_region',
                         what = 'MSY'),
    "what = 'SPR'"
  )

})

test_that("The deprecated BH_MSY alias reaches the same guard", {

  mean_rec_data <- sgl_rg_sable_data
  mean_rec_data$rec_model <- 0 # mean_rec

  # BH_MSY warns and maps to MSY, so the guard has to fire on the mapped name.
  expect_error(
    suppressWarnings(
      Get_Reference_Points(data = mean_rec_data,
                           rep = sgl_rg_sable_rep,
                           SPR_x = 0.4,
                           type = 'single_region',
                           what = 'BH_MSY')
    ),
    "requires a stock-recruit curve"
  )

})

test_that("Every multi region MSY variant is guarded", {

  mean_rec_data <- mlt_rg_sable_data
  mean_rec_data$rec_model <- 0 # mean_rec

  for(w in c("independent_MSY", "global_MSY", "local_MSY")) {
    expect_error(
      Get_Reference_Points(data = mean_rec_data,
                           rep = mlt_rg_sable_rep,
                           SPR_x = 0.4,
                           type = 'multi_region',
                           what = w),
      "requires a stock-recruit curve",
      info = w
    )
  } # end w loop

})

test_that("A stock-recruit penalty under mean recruitment is still rejected, with its own note", {

  # sr_penalty fits a curve, but only as a penalty against the recruitment
  # deviations, and its scale sr_R0 is not reported, so the solvers cannot use it.
  mean_rec_data <- sgl_rg_sable_data
  mean_rec_data$rec_model <- 0 # mean_rec
  mean_rec_data$sr_penalty <- 1 # beverton-holt penalty on the residual

  expect_error(
    Get_Reference_Points(data = mean_rec_data,
                         rep = sgl_rg_sable_rep,
                         SPR_x = 0.4,
                         type = 'single_region',
                         what = 'MSY'),
    "sr_penalty = 1"
  )

})

test_that("SPR reference points are unaffected by rec_model, since they never touch the curve", {

  # SPR scales spawning biomass per recruit by mean recruitment from rep$Rec, so
  # it is well defined under every rec_model and must return the same numbers.
  bh_spr <- Get_Reference_Points(data = sgl_rg_sable_data,
                                 rep = sgl_rg_sable_rep,
                                 SPR_x = 0.4,
                                 type = 'single_region',
                                 what = 'SPR')

  mean_rec_data <- sgl_rg_sable_data
  mean_rec_data$rec_model <- 0 # mean_rec

  mean_rec_spr <- expect_no_error(
    Get_Reference_Points(data = mean_rec_data,
                         rep = sgl_rg_sable_rep,
                         SPR_x = 0.4,
                         type = 'single_region',
                         what = 'SPR')
  )

  expect_equal(mean_rec_spr$f_ref_pt, bh_spr$f_ref_pt)
  expect_equal(mean_rec_spr$b_ref_pt, bh_spr$b_ref_pt)
  expect_equal(mean_rec_spr$virgin_b_ref_pt, bh_spr$virgin_b_ref_pt)

})

test_that("The solvers themselves reject rec_model 0, for data lists built by hand", {

  # Get_Reference_Points is the guarded entry point, but the solvers are reachable
  # directly, so they have their own check rather than trusting the caller.
  for(solver in c("single_region_Fmsy", "global_Fmsy", "local_Fmsy_sglpop", "local_Fmsy_multipop")) {
    expect_error(
      do.call(get(solver, envir = asNamespace("SPoRC")),
              list(pars = list(log_Fmsy = log(0.1)), data = list(rec_model = 0))),
      "no stock-recruit curve",
      info = solver
    )
  } # end solver loop

})

test_that("Beverton-Holt and Ricker MSY still run, and an absent rec_model still means Beverton-Holt", {

  bh_data <- sgl_rg_sable_data
  bh_data$rec_model <- 1 # beverton-holt

  bh <- Get_Reference_Points(data = bh_data,
                             rep = sgl_rg_sable_rep,
                             SPR_x = 0.4,
                             type = 'single_region',
                             what = 'MSY')

  # The packaged data lists predate rec_model, so NULL has to keep meaning 1.
  expect_null(sgl_rg_sable_data$rec_model)
  absent <- Get_Reference_Points(data = sgl_rg_sable_data,
                                 rep = sgl_rg_sable_rep,
                                 SPR_x = 0.4,
                                 type = 'single_region',
                                 what = 'MSY')

  expect_equal(absent$f_ref_pt, bh$f_ref_pt)
  expect_equal(absent$b_ref_pt, bh$b_ref_pt)

  ricker_data <- sgl_rg_sable_data
  ricker_data$rec_model <- 2 # ricker

  ricker <- Get_Reference_Points(data = ricker_data,
                                 rep = sgl_rg_sable_rep,
                                 SPR_x = 0.4,
                                 type = 'single_region',
                                 what = 'MSY')

  # The two curves agree only at (S0, R0), so Fmsy has to differ between them.
  expect_true(all(is.finite(ricker$f_ref_pt)))
  expect_false(isTRUE(all.equal(ricker$f_ref_pt, bh$f_ref_pt)))

})
