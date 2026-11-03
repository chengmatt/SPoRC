# get_key_quants under the same choices as Setup_Mod_Projection gives the in-model projection's catch advice
# (1e-6): a catch cap on dusky, nonzero projected state innovations and recruitment deviations on a state-space
# toy, propagated growth and selectivity on the Pacific cod bridge. the deviations are refused without the fit's
# parameters, and the defaults give what the function always gave.

library(SPoRC)
library(testthat)
data("sgl_rg_dusky_data")
data("sgl_rg_ebs_pcod_data")

n_proj_yrs <- 4
with_devs <- function(build) {
  body(build) <- do.call(substitute, list(body(build), list(Setup_Mod_Dim = bquote(function(...) Setup_Mod_Dim(..., n_proj_yrs_devs = .(n_proj_yrs - 1))))))
  build
}
rule <- hcr_threshold

test_that("a catch cap gives the in-model catch advice on dusky", {

  il <- with_devs(build_goa_dusky_input)(sgl_rg_dusky_data)
  fit0 <- fit_model(il$data, il$par, il$map, newton_loops = 2, silent = TRUE)
  refpt_opt <- list(SPR_x = 0.4, t_spawn = 0, sex_ratio_f = array(0.5, dim = c(1, 1)), calc_rec_st_yr = 1, rec_age = 1, type = "single_region", what = "SPR")
  rp <- Get_Reference_Points(data = fit0$data, rep = fit0$rep, SPR_x = 0.4, t_spawn = 0, sex_ratio_f = refpt_opt$sex_ratio_f, type = "single_region", what = "SPR")

  # the same projection inside the model, at the estimates
  free <- Setup_Mod_Projection(il, n_proj_yrs = n_proj_yrs, f_opt = "hcr", f_ref_pt = rp$f_ref_pt, b_ref_pt = rp$b_ref_pt, rec_opt = "propagate")
  fit_free <- fit_model(free$data, fit0$env$parList(), free$map, do_optim = FALSE, silent = TRUE)
  cap <- 0.85 * sum(fit_free$rep$proj_Catch[,,2,,])
  capped <- Setup_Mod_Projection(il, n_proj_yrs = n_proj_yrs, f_opt = "hcr", f_ref_pt = rp$f_ref_pt, b_ref_pt = rp$b_ref_pt, rec_opt = "propagate", F_constraint = catch_cap(cap))
  fit_cap <- fit_model(capped$data, fit0$env$parList(), capped$map, do_optim = FALSE, silent = TRUE)

  proj_opt <- list(n_proj_yrs = n_proj_yrs, n_avg_yrs = 1, HCR_function = rule, recruitment_opt = "mean_rec", fmort_opt = "HCR", rec_opt = "propagate")
  kq_free <- get_key_quants(data = list(fit_free$data), rep = list(fit_free$rep), reference_points_opt = refpt_opt, proj_model_opt = proj_opt,
                            model_names = "dusky", pars = list(fit_free$env$parList()))
  proj_opt$F_constraint <- catch_cap(cap)
  kq_cap <- get_key_quants(data = list(fit_cap$data), rep = list(fit_cap$rep), reference_points_opt = refpt_opt, proj_model_opt = proj_opt,
                           model_names = "dusky", pars = list(fit_cap$env$parList()))

  expect_equal(kq_free[[1]]$Catch_Advice, sum(fit_free$rep$proj_Catch[,,2,,]), tolerance = 1e-6)
  expect_equal(kq_free[[1]]$Terminal_SSB, sum(fit_free$rep$proj_SSB[,,1]), tolerance = 1e-6)
  expect_equal(kq_free[[1]]$F_Ref_Pt, round(as.numeric(rp$f_ref_pt), 5))
  expect_equal(kq_cap[[1]]$Catch_Advice, cap, tolerance = 1e-4)
  expect_equal(kq_cap[[1]]$Catch_Advice, sum(fit_cap$rep$proj_Catch[,,2,,]), tolerance = 1e-6)

  # the deviations come from the fit, so they are refused without it, and the defaults are the old behavior
  expect_error(get_key_quants(data = list(fit_free$data), rep = list(fit_free$rep), reference_points_opt = refpt_opt, proj_model_opt = proj_opt, model_names = "dusky"), "pars")
  proj_opt$rec_opt <- NULL
  proj_opt$F_constraint <- NULL
  kq_old <- get_key_quants(data = list(fit_free$data), rep = list(fit_free$rep), reference_points_opt = refpt_opt, proj_model_opt = proj_opt, model_names = "dusky")
  expect_equal(kq_old[[1]]$Terminal_SSB, kq_free[[1]]$Terminal_SSB, tolerance = 1e-8)

})

test_that("propagated growth and selectivity give the in-model catch advice on the Pacific cod bridge", {

  il <- suppressWarnings(seed_ebs_pcod_mle(suppressMessages(with_devs(build_ebs_pcod_input)(sgl_rg_ebs_pcod_data)), sgl_rg_ebs_pcod_data))
  il$data$growth_tv_type <- 0 # curve mode builds every year up front, which is what projected years need
  obj0 <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)
  n_pop <- il$data$n_pop
  n_regions <- il$data$n_regions
  t_spawn <- if(is.null(il$data$t_spawn)) 0 else il$data$t_spawn
  refpt_opt <- list(SPR_x = 0.4, t_spawn = t_spawn, sex_ratio_f = array(0.5, dim = c(n_pop, n_regions)), calc_rec_st_yr = 1, rec_age = 1, type = "single_region", what = "SPR")
  rp <- Get_Reference_Points(data = obj0$data, rep = obj0$rep, SPR_x = 0.4, t_spawn = t_spawn, sex_ratio_f = refpt_opt$sex_ratio_f, type = "single_region", what = "SPR")

  # the same projection inside the model, fishing at F40 every year with growth and selectivity propagated
  ip <- Setup_Mod_Projection(il, n_proj_yrs = n_proj_yrs, f_opt = "input", f_input = array(as.numeric(rp$f_ref_pt), dim = c(n_regions, n_proj_yrs)),
                             rec_opt = "mean", growth_opt = "propagate", sel_opt = "propagate")
  obj <- fit_model(ip$data, ip$par, ip$map, do_optim = FALSE, silent = TRUE)

  proj_opt <- list(n_proj_yrs = n_proj_yrs, n_avg_yrs = 1, HCR_function = rule, recruitment_opt = "mean_rec", fmort_opt = "Input", growth_opt = "propagate", sel_opt = "propagate")
  kq <- get_key_quants(data = list(obj$data), rep = list(obj$rep), reference_points_opt = refpt_opt, proj_model_opt = proj_opt, model_names = "pcod")
  expect_equal(kq[[1]]$Catch_Advice, sum(obj$rep$proj_Catch[,,2,,]), tolerance = 1e-6)
  expect_equal(kq[[1]]$Terminal_SSB, sum(obj$rep$proj_SSB[,,1]), tolerance = 1e-6)

  # averaging instead moves the advice, since the projected growth sits below the terminal year
  proj_opt$growth_opt <- "average"
  proj_opt$sel_opt <- "average"
  kq_avg <- get_key_quants(data = list(obj$data), rep = list(obj$rep), reference_points_opt = refpt_opt, proj_model_opt = proj_opt, model_names = "pcod")
  expect_false(isTRUE(all.equal(kq_avg[[1]]$Catch_Advice, kq[[1]]$Catch_Advice, tolerance = 1e-4)))

})

test_that("nonzero projected state innovations and recruitment deviations give the in-model catch advice on a state-space toy", {

  naa_dims <- list(n_regions = 1, n_sexes = 1, n_fish_fleets = 1, n_seas = 1, n_yrs = 22, n_ages = 8)
  n_yrs <- naa_dims$n_yrs
  build <- sweep_input
  body(build) <- do.call(substitute, list(body(build), list(Setup_Mod_Dim = bquote(function(...) Setup_Mod_Dim(..., n_proj_yrs_devs = .(n_proj_yrs - 1))))))

  # a deterministic pass seeds the state over the fitted years
  il0 <- build(dims = naa_dims)
  obj0 <- fit_model(il0$data, il0$par, il0$map, do_optim = FALSE, silent = TRUE)
  seed <- log(array(obj0$report(obj0$par)$NAA[,,1:n_yrs,1,,], dim = c(1, 1, n_yrs, 1, naa_dims$n_ages, 1)))
  il <- build(dims = naa_dims, biol = list(NAA_re = "iid", NAA_sigma_spec = "fix", M_spec = "fix", ln_NAA = seed))
  refpt_opt <- list(SPR_x = 0.4, t_spawn = 0, sex_ratio_f = array(1, dim = c(1, 1)), calc_rec_st_yr = 1, rec_age = 1, type = "single_region", what = "SPR")
  obj_base <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)
  rp <- Get_Reference_Points(data = obj_base$data, rep = obj_base$rep, SPR_x = 0.4, t_spawn = 0, sex_ratio_f = refpt_opt$sex_ratio_f, type = "single_region", what = "SPR")

  # a fit's projected deviations sit at zero, so set them away from it to see both routes read them
  pert <- il$par
  pert$ln_NAA[1,1,n_yrs + 1:(n_proj_yrs - 1),1,il$data$naa_re_ages,1] <- 0.2
  pert$ln_RecDevs[1,1,n_yrs + 1:(n_proj_yrs - 1)] <- c(0.6, -0.4, 0.3)
  ip <- Setup_Mod_Projection(il, n_proj_yrs = n_proj_yrs, f_opt = "hcr", f_ref_pt = rp$f_ref_pt, b_ref_pt = rp$b_ref_pt, rec_opt = "propagate", naa_opt = "propagate")
  obj <- fit_model(ip$data, pert, ip$map, do_optim = FALSE, silent = TRUE)

  proj_opt <- list(n_proj_yrs = n_proj_yrs, n_avg_yrs = 1, HCR_function = rule, recruitment_opt = "mean_rec", fmort_opt = "HCR", rec_opt = "propagate", naa_opt = "propagate")
  kq <- get_key_quants(data = list(obj$data), rep = list(obj$rep), reference_points_opt = refpt_opt, proj_model_opt = proj_opt, model_names = "toy", pars = list(pert))
  expect_equal(kq[[1]]$Catch_Advice, sum(obj$rep$proj_Catch[,,2,,]), tolerance = 1e-6)
  expect_equal(kq[[1]]$Terminal_SSB, sum(obj$rep$proj_SSB[,,1]), tolerance = 1e-6)

  # each set of deviations moves the advice: the state alone, then neither
  proj_opt$rec_opt <- "mean"
  kq_state <- get_key_quants(data = list(obj$data), rep = list(obj$rep), reference_points_opt = refpt_opt, proj_model_opt = proj_opt, model_names = "toy", pars = list(pert))
  proj_opt$naa_opt <- "terminal"
  kq_none <- get_key_quants(data = list(obj$data), rep = list(obj$rep), reference_points_opt = refpt_opt, proj_model_opt = proj_opt, model_names = "toy")
  expect_gt(abs(kq[[1]]$Catch_Advice / kq_state[[1]]$Catch_Advice - 1), 0.01)
  expect_gt(abs(kq_state[[1]]$Catch_Advice / kq_none[[1]]$Catch_Advice - 1), 0.1)

})
