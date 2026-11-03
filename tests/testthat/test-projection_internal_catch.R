# Catch given by year, inside the objective and in the post hoc projection. A catch target is solved
# numerically, then refined with two Newton steps the model differentiates through, as a catch cap is. On
# dusky the in-model projection takes the catch and equals the post hoc one (1e-6), the catch has a standard
# error of zero, F's derivative matches central differences of the re-solved projection (1e-4), and years
# left NA are fished under f_opt with any constraint. Annual and seasonal targets on a two-region,
# two-season stock match numerical gradients (1e-5), including the gradient in the catch itself.

library(SPoRC)
library(testthat)
data("sgl_rg_dusky_data")

n_proj_yrs <- 6
build_ext <- build_goa_dusky_input
body(build_ext) <- do.call(substitute, list(body(build_ext), list(Setup_Mod_Dim = quote(function(...) Setup_Mod_Dim(..., n_proj_yrs_devs = 5)))))
input_base <- build_ext(sgl_rg_dusky_data)
n_yrs <- length(input_base$data$years)
rec_yrs <- 1:(n_yrs - 1) # the setup's default window, calc_rec_st_yr 1 and rec_age 1

# the post hoc projection at a fit's own values: terminal year demographics and the fit's projected deviations,
# built by hand since dusky is one population, region, season and sex
post_hoc_projection <- function(fit, fmort_opt, catch_input = NULL, catch_fallback_opt = "Input", f = 0, f_ref_pt = 0, b_ref_pt = 1,
                                HCR_function = NULL, F_constraint = NULL) {

  d <- fit$data
  rp <- fit$rep
  n_ages <- length(d$ages)
  n_sexes <- d$n_sexes
  n_fleets <- d$n_fish_fleets
  biol_d <- c(1, 1, n_proj_yrs, 1, n_ages, n_sexes)

  WAA <- MatAA <- natmort <- array(0, dim = biol_d)
  WAA_fish <- fish_sel <- ret_sel <- array(0, dim = c(biol_d, n_fleets))
  for(y in 1:n_proj_yrs) {
    WAA[1,1,y,1,,] <- d$WAA[1,1,n_yrs,1,,] # weight and maturity are data for dusky, so the report has neither
    MatAA[1,1,y,1,,] <- d$MatAA[1,1,n_yrs,1,,]
    natmort[1,1,y,1,,] <- rp$natmort[1,1,n_yrs,1,,]
    WAA_fish[1,1,y,1,,,] <- d$WAA_fish[1,1,n_yrs,1,,,]
    fish_sel[1,1,y,1,,,] <- rp$fish_sel[1,1,n_yrs,1,,,]
    ret_sel[1,1,y,1,,,] <- rp$ret_sel[1,1,n_yrs,1,,,]
  } # end y loop

  rec_devs <- array(1, dim = c(1, 1, n_proj_yrs)) # the fit's projected deviations
  rec_devs[1,1,2:n_proj_yrs] <- exp(fit$env$parList()$ln_RecDevs[1,1,n_yrs + 1:(n_proj_yrs - 1)])

  Do_Population_Projection(n_proj_yrs = n_proj_yrs, n_pop = 1, n_regions = 1, n_ages = n_ages, n_sexes = n_sexes, n_fish_fleets = n_fleets,
                           sexratio = array(d$sexratio[1,1,n_yrs], dim = c(1, 1, n_proj_yrs, n_sexes)),
                           do_recruits_move = 0,
                           recruitment = array(rp$Rec[1,1,rec_yrs], dim = c(1, 1, length(rec_yrs))),
                           terminal_NAA = array(rp$NAA[1,1,n_yrs,,,], dim = c(1, 1, 1, n_ages, n_sexes)),
                           terminal_NAA0 = array(rp$NAA0[1,1,n_yrs,,,], dim = c(1, 1, 1, n_ages, n_sexes)),
                           terminal_F = array(rp$Fmort[1,n_yrs,,], dim = c(1, 1, n_fleets)),
                           dmr = array(rp$dmr[1,n_yrs,,], dim = c(1, 1, n_fleets)),
                           natmort = natmort, WAA = WAA, WAA_fish = WAA_fish, MatAA = MatAA, fish_sel = fish_sel, ret_sel = ret_sel,
                           Movement = array(1, dim = c(1, 1, 1, n_proj_yrs, 1, n_ages, n_sexes)),
                           f_ref_pt = array(if(catch_fallback_opt == "Input") f else f_ref_pt, dim = c(1, n_proj_yrs)),
                           b_ref_pt = array(b_ref_pt, dim = c(1, 1, n_proj_yrs)),
                           HCR_function = HCR_function, recruitment_opt = "mean_rec", fmort_opt = fmort_opt,
                           catch_input = catch_input, catch_fallback_opt = catch_fallback_opt,
                           rec_devs = rec_devs, F_constraint = F_constraint, t_spawn = 0)
}

annual_catch <- function(rep) {
  catch <- rep(0, n_proj_yrs)
  for(y in 1:n_proj_yrs) catch[y] <- sum(rep$proj_Catch[,,y,,])
  catch
}

se_of <- function(sdr, name) as.numeric(sdr$sd[names(sdr$value) == name])

# the plain fit, and a catch below the terminal year's
fit_plain <- fit_model(input_base$data, input_base$par, input_base$map, newton_loops = 2, silent = TRUE)
catch_amt <- 0.8 * sum(fit_plain$rep$PredCatch[,,n_yrs,,])
catch_all <- array(NA, dim = c(1, n_proj_yrs))
catch_all[1, 2:n_proj_yrs] <- catch_amt

test_that("Setup_Mod_Projection checks catch_input and stores it with year 1 left out", {

  expect_error(Setup_Mod_Projection(input_base, n_proj_yrs = n_proj_yrs, catch_input = catch_amt), "array")
  expect_error(Setup_Mod_Projection(input_base, n_proj_yrs = n_proj_yrs, catch_input = array(catch_amt, dim = c(1, 5))), "dimensioned")
  bad <- catch_all
  bad[1, 3] <- -1
  expect_error(Setup_Mod_Projection(input_base, n_proj_yrs = n_proj_yrs, catch_input = bad), "negative")
  only_one <- array(NA, dim = c(1, n_proj_yrs))
  only_one[1, 1] <- catch_amt
  expect_error(Setup_Mod_Projection(input_base, n_proj_yrs = n_proj_yrs, catch_input = only_one), "projection year 2")

  # a year is given in every region or in none
  two <- input_base
  two$data$n_regions <- 2
  part <- array(NA, dim = c(2, n_proj_yrs))
  part[1, 2] <- catch_amt
  part[, 3:n_proj_yrs] <- catch_amt
  expect_error(Setup_Mod_Projection(two, n_proj_yrs = n_proj_yrs, catch_input = part), "partly given in projection year\\(s\\) 2")

  # stored with year 1 blank, since year 1 replays the terminal year, and stored empty when not given
  given <- catch_all
  given[1, 1] <- 1
  ip <- Setup_Mod_Projection(input_base, n_proj_yrs = n_proj_yrs, catch_input = given)
  expect_true(is.na(ip$data$proj_catch_input[1, 1]))
  expect_equal(ip$data$proj_catch_input[1, 2:n_proj_yrs], rep(catch_amt, n_proj_yrs - 1))
  ip_none <- Setup_Mod_Projection(input_base, n_proj_yrs = n_proj_yrs)
  expect_true("proj_catch_input" %in% names(ip_none$data))
  expect_null(ip_none$data$proj_catch_input)

})

test_that("catch given in every year is taken, equals the post hoc projection, and has no standard error", {

  ip <- Setup_Mod_Projection(input_base, n_proj_yrs = n_proj_yrs, catch_input = catch_all)
  fit <- fit_model(ip$data, input_base$par, ip$map, newton_loops = 2, silent = TRUE)

  # the catch is taken, and the projection is the post hoc one
  expect_equal(annual_catch(fit$rep)[2:n_proj_yrs], rep(catch_amt, n_proj_yrs - 1), tolerance = 1e-6)
  ph <- post_hoc_projection(fit, "Catch", catch_input = catch_all)
  expect_equal(fit$rep$proj_SSB, ph$proj_SSB, tolerance = 1e-6)
  expect_equal(fit$rep$proj_F, ph$proj_F, tolerance = 1e-6)

  # the catch is fixed, so its uncertainty moves into F and spawning biomass
  sdr <- RTMB::sdreport(fit)
  expect_true(all(se_of(sdr, "log_proj_Catch")[2:n_proj_yrs] < 1e-8))
  expect_true(all(se_of(sdr, "log_proj_F")[2:n_proj_yrs] > 0.05))
  expect_gt(se_of(sdr, "log_proj_SSB")[n_proj_yrs], se_of(sdr, "log_proj_SSB")[2])

})

test_that("F's derivative through the catch solve matches central differences of the re-solved projection", {

  ip <- Setup_Mod_Projection(input_base, n_proj_yrs = n_proj_yrs, catch_input = catch_all)
  fit <- fit_model(ip$data, input_base$par, ip$map, newton_loops = 2, silent = TRUE)
  sdr <- RTMB::sdreport(fit)
  rows <- which(names(sdr$value) == "log_proj_F")[2:n_proj_yrs]

  # the model's own derivative of every reported quantity, in the order sdreport reports them
  obj_rep <- RTMB::MakeADFun(cmb(SPoRC_rtmb, fit$data), parameters = fit$env$parList(), map = ip$map, ADreport = TRUE, silent = TRUE)
  p0 <- fit$env$last.par.best
  expect_equal(as.numeric(obj_rep$fn(p0))[rows], as.numeric(sdr$value[rows]), tolerance = 1e-10)
  jac <- obj_rep$gr(p0)[rows, , drop = FALSE]

  # along three directions, against the plain projection solved again at each side
  log_F_at <- function(p) log(fit$report(p)$proj_F[1, 2:n_proj_yrs])
  set.seed(1)
  h <- 1e-4
  for(k in 1:3) {
    v <- stats::rnorm(length(p0))
    v <- v / sqrt(sum(v^2))
    central <- (log_F_at(p0 + h * v) - log_F_at(p0 - h * v)) / (2 * h)
    expect_equal(as.numeric(jac %*% v), central, tolerance = 1e-4)
  } # end k loop

})

test_that("years left NA are fished under f_opt, and constraints apply only in those years", {

  rp <- Get_Reference_Points(data = fit_plain$data, rep = fit_plain$rep, SPR_x = 0.4, t_spawn = 0, calc_rec_st_yr = 1, rec_age = 1,
                             sex_ratio_f = array(0.5, dim = c(1, 1)), type = "single_region", what = "SPR")
  catch_first <- array(NA, dim = c(1, n_proj_yrs))
  catch_first[1, 2:3] <- catch_amt
  cap <- 0.9 * catch_amt # below the given catch, so it would bind there if it applied

  ip <- Setup_Mod_Projection(input_base, n_proj_yrs = n_proj_yrs, f_opt = "hcr", f_ref_pt = rp$f_ref_pt, b_ref_pt = rp$b_ref_pt,
                             catch_input = catch_first, F_constraint = catch_cap(cap))
  fit <- fit_model(ip$data, input_base$par, ip$map, newton_loops = 2, silent = TRUE)
  catch <- annual_catch(fit$rep)

  expect_equal(catch[2:3], rep(catch_amt, 2), tolerance = 1e-6) # the given catch, above the cap
  expect_true(all(catch[4:n_proj_yrs] <= cap * (1 + 1e-4))) # the rule's years, held under it
  expect_true(any(fit$rep$proj_F_rule[1, 4:n_proj_yrs] > fit$rep$proj_F[1, 4:n_proj_yrs] * (1 + 1e-6))) # and the cap bound there

  ph <- post_hoc_projection(fit, "Catch", catch_input = catch_first, catch_fallback_opt = "HCR",
                            f_ref_pt = as.numeric(rp$f_ref_pt), b_ref_pt = as.numeric(rp$b_ref_pt),
                            HCR_function = function(x, frp, brp) hcr_threshold(x, frp, brp, 0.05), F_constraint = catch_cap(cap))
  expect_equal(fit$rep$proj_F, ph$proj_F, tolerance = 1e-6)
  expect_equal(fit$rep$proj_SSB, ph$proj_SSB, tolerance = 1e-6)

})

# a two-region, two-season, one-sex, one-fleet stock with mixing, so the joint and the seasonal solves are both live
two_season_inputs <- function(n_proj_yrs = 4, n_ages = 10) {
  n_regions <- 2
  n_seas <- 2
  ages <- 1:n_ages
  biol_d <- c(1, n_regions, n_proj_yrs, n_seas, n_ages, 1)

  natmort <- array(0.15, dim = biol_d)
  WAA <- MatAA <- array(0, dim = biol_d)
  fish_sel <- array(0, dim = c(biol_d, 1))
  for(r in 1:n_regions) for(y in 1:n_proj_yrs) for(seas in 1:n_seas) {
    WAA[1,r,y,seas,,1] <- (1 - exp(-0.2 * ages))^3 * 5
    MatAA[1,r,y,seas,,1] <- 1 / (1 + exp(-0.8 * (ages - 4)))
    fish_sel[1,r,y,seas,,1,1] <- 1 / (1 + exp(-1.1 * (ages - 3 - 0.3 * r)))
  } # end r, y, seas loops

  Movement <- array(0, dim = c(1, n_regions, n_regions, n_proj_yrs, n_seas, n_ages, 1))
  for(y in 1:n_proj_yrs) for(seas in 1:n_seas) for(a in 1:n_ages) Movement[1,,,y,seas,a,1] <- matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2)
  sgl_move <- array(0, dim = c(1, n_regions, n_regions, n_proj_yrs, n_ages, 1))
  for(y in 1:n_proj_yrs) for(a in 1:n_ages) sgl_move[1,,,y,a,1] <- diag(n_regions)

  terminal_NAA <- array(0, dim = c(1, n_regions, n_seas, n_ages, 1))
  for(r in 1:n_regions) for(seas in 1:n_seas) terminal_NAA[1,r,seas,,1] <- 1e6 * exp(-0.3 * (ages - 1)) * (1 + 0.3 * r)
  terminal_F <- array(0, dim = c(n_regions, n_seas, 1))
  terminal_F[,1,1] <- c(0.06, 0.05)
  terminal_F[,2,1] <- c(0.03, 0.04)

  list(n_proj_yrs = n_proj_yrs, n_pop = 1, n_regions = n_regions, n_ages = n_ages, n_sexes = 1, n_fish_fleets = 1, n_seas = n_seas,
       seasdur = c(0.5, 0.5), spawn_seas = 1, t_spawn = 0, natal_region = 1, do_recruits_move = 0,
       rec_seas_prop = array(c(1, 0), dim = c(1, n_seas)),
       stray_rate = array(0, dim = c(1, n_proj_yrs)),
       sexratio = array(1, dim = c(1, n_regions, n_proj_yrs, 1)),
       recruitment = array(3e5, dim = c(1, n_regions, 10)),
       terminal_NAA = terminal_NAA, terminal_NAA0 = terminal_NAA * 1.5, terminal_F = terminal_F,
       natmort = natmort, WAA = WAA, WAA_fish = array(WAA, dim = c(biol_d, 1)), MatAA = MatAA, fish_sel = fish_sel,
       Movement = Movement, sgl_seas_spawning_movement = sgl_move, recruitment_opt = "mean_rec",
       f_ref_pt = array(0.08, dim = c(n_regions, n_proj_yrs)))
}

test_that("annual and seasonal catch targets on two regions take the catch, and their gradients match finite differences", {

  inp <- two_season_inputs()
  annual <- array(NA, dim = c(2, 4))
  annual[, 2:4] <- c(9e4, 7e4)
  seasonal <- array(NA, dim = c(2, 4, 2))
  seasonal[, 2:4, 1] <- c(5e4, 4e4)
  seasonal[, 2:4, 2] <- c(3e4, 3.5e4)

  for(targets in list(annual, seasonal)) {

    out <- do.call(Do_Population_Projection, c(inp, list(fmort_opt = "Catch", catch_input = targets, catch_fallback_opt = "Input")))
    for(r in 1:2) for(y in 2:4) {
      if(length(dim(targets)) == 2) expect_equal(sum(out$proj_Catch[,r,y,,]), targets[r,y], tolerance = 1e-6)
      else for(seas in 1:2) expect_equal(sum(out$proj_Catch[,r,y,seas,]), targets[r,y,seas], tolerance = 1e-6)
    } # end r, y loops

    # a log scale on the terminal numbers and on natural mortality, so the solved F moves with both
    f <- function(p, what) {
      inp_p <- inp
      inp_p$terminal_NAA <- inp$terminal_NAA * exp(p[1])
      inp_p$natmort <- inp$natmort * exp(p[2])
      out <- do.call(Do_Population_Projection, c(inp_p, list(fmort_opt = "Catch", catch_input = targets, catch_fallback_opt = "Input")))
      if(what == "proj_F") sum(out$proj_F[, 2:4]) else sum(out$proj_SSB[,, 2:4])
    }
    p0 <- c(0.1, -0.05)
    for(what in c("proj_F", "proj_SSB")) {
      f_what <- function(p) f(p, what)
      tape <- RTMB::MakeTape(f_what, p0)
      expect_equal(tape(p0), f_what(p0), tolerance = 1e-10)
      expect_equal(as.vector(tape$jacobian(p0)), numDeriv::grad(f_what, p0), tolerance = 1e-5)
    } # end what loop

    # and in the catch itself, which may be a quantity the model differentiates. year 1 is filled, since it is
    # not read and an NA scaled by a differentiated number is no longer NA
    targets_full <- targets
    targets_full[is.na(targets_full)] <- 1
    g <- function(log_scale) {
      out <- do.call(Do_Population_Projection, c(inp, list(fmort_opt = "Catch", catch_input = targets_full * exp(log_scale), catch_fallback_opt = "Input")))
      sum(out$proj_SSB[,,2:4])
    }
    tape_g <- RTMB::MakeTape(g, 0)
    expect_equal(as.vector(tape_g$jacobian(0)), numDeriv::grad(g, 0), tolerance = 1e-5)

  } # end targets loop

})
