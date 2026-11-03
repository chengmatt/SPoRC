# The female share of recruits behind the reference points solved inside the model. By default a two-sex model
# reads its own recruitment sex ratio, so B40 equals Get_Reference_Points at that share (1e-8) and an estimated
# sex ratio moves B40 with the derivative central differences give (1e-5); a share given to the setup is used
# instead. The threshold rule is the default rule, with its cutoff set inside it rather than by hcr_alpha.

library(SPoRC)
library(testthat)

n_proj_yrs <- 4
two_sex_dims <- list(n_regions = 1, n_sexes = 2, n_fish_fleets = 1, n_seas = 1, n_yrs = 22, n_ages = 8)
female_share <- 0.6

# a two-sex toy with its female share of recruits set away from one half
il <- sweep_input(dims = two_sex_dims)
il$par$sexratio_pars[] <- stats::qlogis(female_share)

post_hoc_refpts <- function(obj, share) {
  Get_Reference_Points(data = obj$data, rep = obj$rep, SPR_x = 0.4, type = "single_region", what = "SPR",
                       sex_ratio_f = array(share, dim = c(1, 1)))
}

test_that("a two-sex model's reference points start from its own female share unless one is given", {

  ip <- Setup_Mod_Projection(il, n_proj_yrs = n_proj_yrs, f_opt = "refpt", refpt_what = "SPR", SPR_x = 0.4, rec_opt = "mean")
  expect_null(ip$data$proj_sex_ratio_f)
  obj <- fit_model(ip$data, il$par, ip$map, do_optim = FALSE, silent = TRUE)
  expect_equal(obj$rep$sexratio[1,1,two_sex_dims$n_yrs,1], female_share)

  # B40 is the post hoc one at the model's share, not at one half; F40 is a ratio, so the share leaves it alone
  at_model <- post_hoc_refpts(obj, female_share)
  at_half <- post_hoc_refpts(obj, 0.5)
  expect_equal(as.numeric(obj$rep$proj_b_ref_pt), as.numeric(at_model$b_ref_pt), tolerance = 1e-8)
  expect_equal(as.numeric(obj$rep$proj_f_ref_pt), as.numeric(at_model$f_ref_pt), tolerance = 1e-8)
  expect_gt(abs(as.numeric(obj$rep$proj_b_ref_pt) / as.numeric(at_half$b_ref_pt) - 1), 0.1)

  # a share given to the setup is used instead
  ip_half <- Setup_Mod_Projection(il, n_proj_yrs = n_proj_yrs, f_opt = "refpt", refpt_what = "SPR", SPR_x = 0.4, rec_opt = "mean",
                                  sex_ratio_f = array(0.5, dim = c(1, 1)))
  obj_half <- fit_model(ip_half$data, il$par, ip_half$map, do_optim = FALSE, silent = TRUE)
  expect_equal(as.numeric(obj_half$rep$proj_b_ref_pt), as.numeric(at_half$b_ref_pt), tolerance = 1e-8)

})

test_that("an estimated sex ratio moves B40 with the derivative central differences give", {

  ip <- Setup_Mod_Projection(il, n_proj_yrs = n_proj_yrs, f_opt = "refpt", refpt_what = "SPR", SPR_x = 0.4, rec_opt = "mean")
  ip$map$sexratio_pars <- factor(1) # estimated, so it is a parameter the reference points can be differentiated in
  obj <- RTMB::MakeADFun(cmb(SPoRC_rtmb, ip$data), parameters = il$par, map = ip$map, ADreport = TRUE, silent = TRUE)
  p0 <- obj$par
  pos <- which(names(p0) == "sexratio_pars")
  expect_length(pos, 1)

  # the row of log_proj_b_ref_pt among the reported quantities, read off a plain sdreport-style build
  obj_fit <- RTMB::MakeADFun(cmb(SPoRC_rtmb, ip$data), parameters = il$par, map = ip$map, silent = TRUE)
  sdr <- RTMB::sdreport(obj_fit, par.fixed = p0, hessian.fixed = diag(length(p0)))
  row <- which(names(sdr$value) == "log_proj_b_ref_pt")
  slope_model <- obj$gr(p0)[row, pos]

  log_b40_at <- function(p) log(as.numeric(obj_fit$report(p)$proj_b_ref_pt))
  h <- 1e-4
  up <- p0
  up[pos] <- up[pos] + h
  down <- p0
  down[pos] <- down[pos] - h
  slope_central <- (log_b40_at(up) - log_b40_at(down)) / (2 * h)
  expect_equal(slope_model, slope_central, tolerance = 1e-5)
  expect_gt(abs(slope_model), 0.1)

})

test_that("the threshold rule is the default rule and its cutoff is set inside it", {

  ip <- Setup_Mod_Projection(il, n_proj_yrs = n_proj_yrs, f_opt = "refpt", refpt_what = "SPR", SPR_x = 0.4, rec_opt = "mean")
  expect_identical(ip$data$proj_hcr_function, hcr_threshold)
  expect_false("proj_hcr_alpha" %in% names(ip$data))
  expect_error(Setup_Mod_Projection(il, n_proj_yrs = n_proj_yrs, hcr_alpha = 0.1), "unused argument")

  # a cutoff of one half is a different rule, which shows in the stored function
  cutoff <- function(x, frp, brp) hcr_threshold(x, frp, brp, alpha = 0.5)
  ip_cut <- Setup_Mod_Projection(il, n_proj_yrs = n_proj_yrs, f_opt = "refpt", refpt_what = "SPR", SPR_x = 0.4, rec_opt = "mean",
                                 hcr_function = cutoff)
  expect_identical(ip_cut$data$proj_hcr_function, cutoff)

})
