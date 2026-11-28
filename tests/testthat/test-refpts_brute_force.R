# SPR and MSY reference points for one region against a year class followed for 1500 years and a yield curve
# maximized numerically (1e-6). Spawning is part way through the year and natural mortality differs at every
# age, so spawn timing and the oldest ages both enter the answer. The fishery discards, with retention rising by age.

# Setup ----------------------------------------------------------------------

# one region, sex and fleet under Beverton-Holt recruitment, evaluated at its starting values
refpt_hand_model <- local({
  cache <- NULL
  function() {
    if(!is.null(cache)) return(cache)

    n_ages <- 6
    natmort_at_age <- c(0.45, 0.35, 0.28, 0.22, 0.17, 0.12) # per year, a different value at every age

    il <- sweep_input(
      dims = list(n_regions = 1, n_sexes = 1, n_fish_fleets = 1, n_yrs = 10, n_ages = n_ages),
      rec = list(rec_model = "bh_rec", h_spec = "fix", t_spawn = 0.4) # spawning 40 pct of the way into the year
    )
    for(a in 1:n_ages) il$data$Fixed_natmort[,,,,a,] <- natmort_at_age[a]
    fit <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)

    # the report read as a fishery that discards the young fish, 40 pct of them dying, so retention
    # differs at every age up to the plus group
    retention_at_age <- c(0.30, 0.50, 0.70, 0.85, 0.95, 1)
    for(a in 1:n_ages) fit$rep$ret_sel[,,,,a,,] <- retention_at_age[a]
    fit$rep$dmr[] <- 0.4

    cache <<- list(data = il$data, rep = fit$rep)
    cache
  }
})

# spawning biomass and retained yield for one recruit fished at F every year of its life
refpt_hand_per_recruit <- function(Fmort, data, rep) {

  n_yrs <- length(data$years)
  n_ages <- length(data$ages)
  numbers <- 1
  sbpr <- ypr <- 0

  for(year in 1:1500) {

    a <- min(year, n_ages) # the plus group keeps its own values every year after
    ret_F <- Fmort * rep$fish_sel[1,1,n_yrs,1,a,1,1] * rep$ret_sel[1,1,n_yrs,1,a,1,1] # retained
    disc_F <- Fmort * rep$fish_sel[1,1,n_yrs,1,a,1,1] * (1 - rep$ret_sel[1,1,n_yrs,1,a,1,1]) * rep$dmr[1,n_yrs,1,1] # discarded and dead
    Z <- rep$natmort[1,1,n_yrs,1,a,1] + ret_F + disc_F

    # one sex holds both, so half the mature biomass at spawning is female
    sbpr <- sbpr + 0.5 * numbers * exp(-data$t_spawn * Z) * data$WAA[1,1,n_yrs,1,a,1] * data$MatAA[1,1,n_yrs,1,a,1]
    ypr <- ypr + ret_F / Z * numbers * (1 - exp(-Z)) * data$WAA_fish[1,1,n_yrs,1,a,1,1]
    numbers <- numbers * exp(-Z)

  } # end year loop

  list(sbpr = sbpr, ypr = ypr)
}

# equilibrium recruitment, spawning biomass and yield at F under the Beverton-Holt curve
refpt_hand_equilibrium <- function(Fmort, data, rep) {

  h <- as.numeric(rep$h_trans[1]) # steepness
  R0 <- as.numeric(rep$R0[1]) # unfished recruitment
  fished <- refpt_hand_per_recruit(Fmort, data, rep)
  unfished <- refpt_hand_per_recruit(0, data, rep)

  rec <- R0 * (4 * h * fished$sbpr - (1 - h) * unfished$sbpr) / ((5 * h - 1) * fished$sbpr)
  list(rec = rec, ssb = rec * fished$sbpr, yield = rec * fished$ypr)
}

# Tests ----------------------------------------------------------------------

test_that("the F at a target SPR leaves that share of unfished spawning biomass per recruit", {

  model <- refpt_hand_model()
  unfished <- refpt_hand_per_recruit(0, model$data, model$rep)

  for(SPR_x in c(0.2, 0.4, 0.6)) {

    ref <- suppressWarnings(Get_Reference_Points(data = model$data, rep = model$rep, SPR_x = SPR_x, t_spawn = model$data$t_spawn,
                                                 type = "single_region", what = "SPR"))
    fished <- refpt_hand_per_recruit(as.numeric(ref$f_ref_pt[1]), model$data, model$rep)

    expect_equal(fished$sbpr / unfished$sbpr, SPR_x, tolerance = 1e-6, label = sprintf("realized SPR at a target of %g", SPR_x))

    # the biomass reference points are per recruit values times mean recruitment
    expect_equal(as.numeric(ref$b_ref_pt[1]), fished$sbpr * as.numeric(ref$mean_rec[1]), tolerance = 1e-6)
    expect_equal(as.numeric(ref$virgin_b_ref_pt[1]), unfished$sbpr * as.numeric(ref$mean_rec[1]), tolerance = 1e-6)

  } # end SPR_x loop

})

test_that("Fmsy maximizes equilibrium yield and Bmsy is the spawning biomass there", {

  model <- refpt_hand_model()
  ref <- suppressWarnings(Get_Reference_Points(data = model$data, rep = model$rep, t_spawn = model$data$t_spawn,
                                               type = "single_region", what = "MSY"))

  # the peak of the yield curve, found numerically
  peak <- stats::optimize(function(Fmort) refpt_hand_equilibrium(Fmort, model$data, model$rep)$yield,
                          interval = c(0.001, 2), maximum = TRUE, tol = 1e-10)
  at_peak <- refpt_hand_equilibrium(peak$maximum, model$data, model$rep)

  expect_equal(as.numeric(ref$f_ref_pt[1]), peak$maximum, tolerance = 1e-5)
  expect_equal(as.numeric(ref$b_ref_pt[1]), at_peak$ssb, tolerance = 1e-5)
  expect_equal(as.numeric(ref$virgin_b_ref_pt[1]), refpt_hand_equilibrium(0, model$data, model$rep)$ssb, tolerance = 1e-6)

})

test_that("spawn timing left unset is read from the model", {

  model <- refpt_hand_model()
  unset <- suppressWarnings(Get_Reference_Points(data = model$data, rep = model$rep, SPR_x = 0.4, type = "single_region", what = "SPR"))
  given <- suppressWarnings(Get_Reference_Points(data = model$data, rep = model$rep, SPR_x = 0.4, t_spawn = model$data$t_spawn,
                                                 type = "single_region", what = "SPR"))
  expect_equal(unset$f_ref_pt, given$f_ref_pt, tolerance = 1e-12)
  expect_equal(unset$b_ref_pt, given$b_ref_pt, tolerance = 1e-12)

  # a data list with no spawn timing is read as spawning at the start of the year
  data_without <- model$data
  data_without$t_spawn <- NULL
  unset_without <- suppressWarnings(Get_Reference_Points(data = data_without, rep = model$rep, SPR_x = 0.4, type = "single_region", what = "SPR"))
  year_start <- suppressWarnings(Get_Reference_Points(data = data_without, rep = model$rep, SPR_x = 0.4, t_spawn = 0,
                                                      type = "single_region", what = "SPR"))
  expect_equal(unset_without$f_ref_pt, year_start$f_ref_pt, tolerance = 1e-12)

})

test_that("spawn timing changes the reference points in this setup", {

  # if it did not, the tests above could not tell a calculation that ignored spawn timing from one that used it
  model <- refpt_hand_model()
  at_spawn_time <- suppressWarnings(Get_Reference_Points(data = model$data, rep = model$rep, SPR_x = 0.4, t_spawn = model$data$t_spawn,
                                                         type = "single_region", what = "SPR"))
  at_year_start <- suppressWarnings(Get_Reference_Points(data = model$data, rep = model$rep, SPR_x = 0.4, t_spawn = 0,
                                                         type = "single_region", what = "SPR"))

  expect_gt(abs(as.numeric(at_year_start$f_ref_pt[1]) / as.numeric(at_spawn_time$f_ref_pt[1]) - 1), 0.02)

})
