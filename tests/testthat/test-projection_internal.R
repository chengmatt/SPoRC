# The projection inside the objective on dusky: off leaves the fit alone; terminal F, input F and the
# threshold rule with a cap each reproduce the post hoc projection (1e-6); sdreport gives the terminal
# SE in year one (1e-3), a catch SE of zero where the cap binds, and the tape branches as the report does.

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
post_hoc_projection <- function(fit, fmort_opt, f = 0, f_ref_pt = 0, b_ref_pt = 1, HCR_function = NULL, F_constraint = NULL, rec_devs = NULL) {

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

  if(is.null(rec_devs)) { # the fit's projected deviations unless given
    rec_devs <- array(1, dim = c(1, 1, n_proj_yrs))
    rec_devs[1,1,2:n_proj_yrs] <- exp(fit$env$parList()$ln_RecDevs[1,1,n_yrs + 1:(n_proj_yrs - 1)])
  }

  Do_Population_Projection(n_proj_yrs = n_proj_yrs, n_pop = 1, n_regions = 1, n_ages = n_ages, n_sexes = n_sexes, n_fish_fleets = n_fleets,
                           sexratio = array(d$sexratio[1,1,n_yrs], dim = c(1, 1, n_proj_yrs, n_sexes)), # one sex, so the terminal year's value repeated
                           do_recruits_move = 0,
                           recruitment = array(rp$Rec[1,1,rec_yrs], dim = c(1, 1, length(rec_yrs))),
                           terminal_NAA = array(rp$NAA[1,1,n_yrs,,,], dim = c(1, 1, 1, n_ages, n_sexes)),
                           terminal_NAA0 = array(rp$NAA0[1,1,n_yrs,,,], dim = c(1, 1, 1, n_ages, n_sexes)),
                           terminal_F = array(rp$Fmort[1,n_yrs,,], dim = c(1, 1, n_fleets)),
                           dmr = array(rp$dmr[1,n_yrs,,], dim = c(1, 1, n_fleets)),
                           natmort = natmort, WAA = WAA, WAA_fish = WAA_fish, MatAA = MatAA, fish_sel = fish_sel, ret_sel = ret_sel,
                           Movement = array(1, dim = c(1, 1, 1, n_proj_yrs, 1, n_ages, n_sexes)),
                           f_ref_pt = array(if(fmort_opt == "Input") f else f_ref_pt, dim = c(1, n_proj_yrs)),
                           b_ref_pt = array(b_ref_pt, dim = c(1, 1, n_proj_yrs)),
                           HCR_function = HCR_function, recruitment_opt = "mean_rec", fmort_opt = fmort_opt,
                           rec_devs = rec_devs, F_constraint = F_constraint, t_spawn = 0)
}

annual_catch <- function(out) {
  catch <- rep(0, n_proj_yrs)
  for(y in 1:n_proj_yrs) catch[y] <- sum(out$proj_Catch[,,y,,])
  catch
}

se_of <- function(sdr, name) as.numeric(sdr$sd[names(sdr$value) == name])

# the plain fit every comparison starts from
fit_plain <- fit_model(input_base$data, input_base$par, input_base$map, newton_loops = 2, silent = TRUE)
sdr_plain <- RTMB::sdreport(fit_plain)

test_that("Setup_Mod_Projection refuses what the fit cannot supply and stores its fields", {

  expect_error(Setup_Mod_Projection(input_base, n_proj_yrs = 8), "n_proj_yrs_devs")
  expect_error(Setup_Mod_Projection(input_base, n_proj_yrs = n_proj_yrs, f_opt = "input"), "f_input")
  expect_error(Setup_Mod_Projection(input_base, n_proj_yrs = n_proj_yrs, f_opt = "hcr"), "f_ref_pt")
  expect_error(Setup_Mod_Projection(input_base, n_proj_yrs = n_proj_yrs, f_opt = "hcr", f_ref_pt = 0.1, b_ref_pt = 100, hcr_function = function(x) x), "x, frp and brp")
  expect_error(Setup_Mod_Projection(input_base, n_proj_yrs = n_proj_yrs, F_constraint = "cap"), "function")
  expect_error(Setup_Mod_Projection(input_base, n_proj_yrs = n_proj_yrs, calc_rec_yrs = 2000:2010), "model years")

  # the optional fields are stored even when empty, since the objective reads them by name
  in_off <- Setup_Mod_Projection(input_base, n_proj_yrs = 0)
  expect_true(all(c("proj_hcr_function", "proj_F_constraint") %in% names(in_off$data)))
  expect_equal(in_off$data$proj_n_yrs, 0)

  # one reference point per region is repeated over the projection years
  in_hcr <- Setup_Mod_Projection(input_base, n_proj_yrs = n_proj_yrs, f_opt = "hcr", f_ref_pt = 0.1, b_ref_pt = 100)
  expect_equal(dim(in_hcr$data$proj_f_ref_pt), c(1, n_proj_yrs))
  expect_equal(dim(in_hcr$data$proj_b_ref_pt), c(1, 1, n_proj_yrs))
  expect_equal(in_hcr$data$proj_calc_rec_yrs, rec_yrs)
})

test_that("with the projection off the fit is what it always was", {

  in_off <- Setup_Mod_Projection(input_base, n_proj_yrs = 0)
  fit_off <- fit_model(in_off$data, in_off$par, in_off$map, newton_loops = 2, silent = TRUE)

  expect_equal(fit_off$optim$objective, fit_plain$optim$objective, tolerance = 1e-10)
  expect_equal(fit_off$env$last.par.best, fit_plain$env$last.par.best, tolerance = 1e-8)
  expect_null(fit_off$rep$proj_SSB)
})

test_that("terminal F and input F reproduce the post hoc projection, and the likelihood is untouched", {

  in_term <- Setup_Mod_Projection(input_base, n_proj_yrs = n_proj_yrs, f_opt = "terminal")
  fit_term <- fit_model(in_term$data, in_term$par, in_term$map, newton_loops = 2, silent = TRUE)
  expect_equal(fit_term$optim$objective, fit_plain$optim$objective, tolerance = 1e-10)
  expect_equal(fit_term$fn(fit_term$par), fit_plain$optim$objective, tolerance = 1e-8) # the rebuilt tape at the estimates

  terminal_f <- sum(fit_term$rep$Fmort[1,n_yrs,,])
  ext <- post_hoc_projection(fit_term, fmort_opt = "Input", f = terminal_f)
  expect_equal(as.numeric(fit_term$rep$proj_SSB[1,1,]), as.numeric(ext$proj_SSB[1,1,]), tolerance = 1e-6)
  expect_equal(annual_catch(fit_term$rep), annual_catch(ext), tolerance = 1e-6)
  expect_equal(as.numeric(fit_term$rep$proj_Total_Biom[1,1,]), as.numeric(ext$proj_Total_Biom[1,1,]), tolerance = 1e-6)
  expect_equal(fit_term$rep$proj_F_rule, fit_term$rep$proj_F)

  f_input <- array(0.05, dim = c(1, n_proj_yrs))
  in_input <- Setup_Mod_Projection(input_base, n_proj_yrs = n_proj_yrs, f_opt = "input", f_input = f_input)
  fit_input <- fit_model(in_input$data, in_input$par, in_input$map, newton_loops = 2, silent = TRUE)
  ext <- post_hoc_projection(fit_input, fmort_opt = "Input", f = 0.05)
  expect_equal(as.numeric(fit_input$rep$proj_SSB[1,1,]), as.numeric(ext$proj_SSB[1,1,]), tolerance = 1e-6)
  expect_equal(as.numeric(fit_input$rep$proj_F[1,2:n_proj_yrs]), rep(0.05, n_proj_yrs - 1))
})

test_that("the threshold rule on given reference points with a cap reproduces the post hoc projection and its uncertainty", {

  rp <- Get_Reference_Points(data = fit_plain$data, rep = fit_plain$rep, SPR_x = 0.4, type = "single_region", what = "SPR")
  free <- post_hoc_projection(fit_plain, fmort_opt = "HCR", f_ref_pt = rp$f_ref_pt, b_ref_pt = rp$b_ref_pt, HCR_function = hcr_threshold)
  cap <- 0.85 * annual_catch(free)[2] # binds in the first projected year

  in_hcr <- Setup_Mod_Projection(input_base, n_proj_yrs = n_proj_yrs, f_opt = "hcr", f_ref_pt = rp$f_ref_pt, b_ref_pt = rp$b_ref_pt,
                                 F_constraint = catch_cap(cap))
  fit <- fit_model(in_hcr$data, in_hcr$par, in_hcr$map, newton_loops = 2, silent = TRUE, do_sdrep = TRUE)
  sdr <- fit$sd_rep
  ext <- post_hoc_projection(fit, fmort_opt = "HCR", f_ref_pt = rp$f_ref_pt, b_ref_pt = rp$b_ref_pt, HCR_function = hcr_threshold,
                             F_constraint = catch_cap(cap))

  # the same projection as the post hoc route
  expect_equal(fit$optim$objective, fit_plain$optim$objective, tolerance = 1e-10)
  expect_equal(as.numeric(fit$rep$proj_SSB[1,1,]), as.numeric(ext$proj_SSB[1,1,]), tolerance = 1e-6)
  expect_equal(annual_catch(fit$rep), annual_catch(ext), tolerance = 1e-6)
  expect_equal(as.numeric(fit$rep$proj_F[1,1:n_proj_yrs]), as.numeric(ext$proj_F[1,1:n_proj_yrs]), tolerance = 1e-6)
  expect_equal(annual_catch(fit$rep)[2], cap, tolerance = 1e-4)
  expect_lt(fit$rep$proj_F[1,2], fit$rep$proj_F_rule[1,2])

  # the tape was built at the estimates, so its branch is the report's
  expect_equal(as.numeric(sdr$value[names(sdr$value) == "log_proj_SSB"]), as.numeric(log(fit$rep$proj_SSB[1,1,] + 1e-5)), tolerance = 1e-10)

  # year one replays the terminal year, so its SE is the terminal SE; the capped year's catch has none;
  # the slack years fish the reference point given as data, so their F has none either
  expect_equal(se_of(sdr, "log_proj_SSB")[1], se_of(sdr_plain, "log_SSB")[n_yrs], tolerance = 1e-3)
  expect_lt(se_of(sdr, "log_proj_Catch")[2], 1e-8)
  expect_gt(se_of(sdr, "log_proj_F")[2], 0.05) # the capped F moves with the biomass
  slack <- which(fit$rep$proj_F[1,1:n_proj_yrs] == fit$rep$proj_F_rule[1,1:n_proj_yrs])[-1]
  expect_true(all(se_of(sdr, "log_proj_F")[slack] < 1e-8))
  expect_true(fit$post_optim_sanity_checks)
})

test_that("the threshold rule on F40 and B40 solved on the tape reproduces the post hoc projection and its uncertainty", {

  in_ref <- Setup_Mod_Projection(input_base, n_proj_yrs = n_proj_yrs, f_opt = "refpt", refpt_what = "SPR", SPR_x = 0.4)
  fit <- fit_model(in_ref$data, in_ref$par, in_ref$map, newton_loops = 2, silent = TRUE, do_sdrep = TRUE)
  sdr <- fit$sd_rep
  rp <- Get_Reference_Points(data = fit$data, rep = fit$rep, SPR_x = 0.4, type = "single_region", what = "SPR")
  ext <- post_hoc_projection(fit, fmort_opt = "HCR", f_ref_pt = rp$f_ref_pt, b_ref_pt = rp$b_ref_pt, HCR_function = hcr_threshold)

  # the same reference points and projection as the post hoc route
  expect_equal(as.numeric(fit$rep$proj_f_ref_pt), as.numeric(rp$f_ref_pt), tolerance = 1e-8)
  expect_equal(as.numeric(fit$rep$proj_b_ref_pt), as.numeric(rp$b_ref_pt), tolerance = 1e-8)
  expect_equal(as.numeric(fit$rep$proj_virgin_b_ref_pt), as.numeric(rp$virgin_b_ref_pt), tolerance = 1e-8)
  expect_equal(as.numeric(fit$rep$proj_SSB[1,1,]), as.numeric(ext$proj_SSB[1,1,]), tolerance = 1e-6)
  expect_equal(annual_catch(fit$rep), annual_catch(ext), tolerance = 1e-6)
  expect_equal(as.numeric(fit$rep$proj_status[1,]), as.numeric(fit$rep$proj_SSB[1,1,]) / sum(rp$b_ref_pt), tolerance = 1e-8)

  # their standard errors are the delta method's, and a year fished at F40 has F40's
  rpu <- Get_Reference_Point_Uncertainty(fit_plain, SPR_x = 0.4, type = "single_region", what = "SPR", sd_rep = sdr_plain)
  expect_equal(se_of(sdr, "log_proj_f_ref_pt"), rpu$refpts$log_se[rpu$refpts$quantity == "f_ref_pt"], tolerance = 2e-3)
  expect_equal(se_of(sdr, "log_proj_b_ref_pt"), rpu$refpts$log_se[rpu$refpts$quantity == "b_ref_pt"], tolerance = 2e-3)
  expect_equal(se_of(sdr, "log_proj_F")[2], se_of(sdr, "log_proj_f_ref_pt"), tolerance = 1e-6)
  expect_equal(se_of(sdr, "log_proj_SSB")[1], se_of(sdr_plain, "log_SSB")[n_yrs], tolerance = 1e-3)
  expect_gt(se_of(sdr, "log_proj_status")[1], 0)
  expect_lt(se_of(sdr, "log_proj_status")[1], se_of(sdr, "log_proj_SSB")[1]) # status and biomass move together, so their ratio is tighter
})

test_that("recruitment runs at the mean or the terminal deviation when asked, and propagating needs the fit's deviations", {

  # a fit without projected deviations can only project at the mean or the terminal deviation
  input_plain <- build_goa_dusky_input(sgl_rg_dusky_data)
  expect_error(Setup_Mod_Projection(input_plain, n_proj_yrs = n_proj_yrs), "rec_opt")
  expect_error(Setup_Mod_Projection(input_plain, n_proj_yrs = n_proj_yrs, rec_opt = "mean", sel_opt = "propagate"), "sel_opt")

  in_mean <- Setup_Mod_Projection(input_plain, n_proj_yrs = n_proj_yrs, rec_opt = "mean")
  fit_mean <- fit_model(in_mean$data, in_mean$par, in_mean$map, newton_loops = 2, silent = TRUE)
  ext_mean <- post_hoc_projection(fit_mean, fmort_opt = "Input", f = sum(fit_mean$rep$Fmort[1,n_yrs,,]), rec_devs = array(1, dim = c(1, 1, n_proj_yrs)))
  expect_equal(as.numeric(fit_mean$rep$proj_SSB[1,1,]), as.numeric(ext_mean$proj_SSB[1,1,]), tolerance = 1e-6)

  # the terminal deviation kept in every projected year
  in_term <- Setup_Mod_Projection(input_base, n_proj_yrs = n_proj_yrs, rec_opt = "terminal")
  fit_term <- fit_model(in_term$data, in_term$par, in_term$map, newton_loops = 2, silent = TRUE)
  dev_term <- exp(fit_term$env$parList()$ln_RecDevs[1,1,n_yrs])
  rp <- fit_term$rep
  d <- fit_term$data
  ext_term <- post_hoc_projection(fit_term, fmort_opt = "Input", f = sum(rp$Fmort[1,n_yrs,,]), rec_devs = array(dev_term, dim = c(1, 1, n_proj_yrs)))
  expect_equal(as.numeric(fit_term$rep$proj_SSB[1,1,]), as.numeric(ext_term$proj_SSB[1,1,]), tolerance = 1e-6)
  expect_equal(as.numeric(fit_term$rep$proj_NAA[1,1,3,1,1,1]), dev_term * mean(rp$Rec[1,1,rec_yrs]) * d$sexratio[1,1,n_yrs], tolerance = 1e-6) # year three's recruits
})
