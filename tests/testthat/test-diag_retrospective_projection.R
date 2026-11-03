# A retrospective peel of a fit with projected deviation years: every extended parameter and data array keeps its
# projected cells after the kept fitted years, the projection's mean recruitment window follows the new terminal
# year, and the peeled model fits with the projection on. Dusky covers recruitment, selectivity and catchability
# deviations, a two-region toy movement, and the Pacific cod bridge growth.

library(SPoRC)
library(testthat)
data("sgl_rg_dusky_data")
data("sgl_rg_ebs_pcod_data")

n_proj_yrs <- 6
n_ext <- n_proj_yrs - 1
with_devs <- function(build) {
  body(build) <- do.call(substitute, list(body(build), list(Setup_Mod_Dim = bquote(function(...) Setup_Mod_Dim(..., n_proj_yrs_devs = .(n_ext))))))
  build
}

# every parameter array whose year dim runs past the fitted years, with that dim, and the data arrays beside them
par_year_dims <- list(ln_RecDevs = 3, ln_fishsel_devs = 2, ln_fishsel_bin_devs = 2, ln_retsel_devs = 2, ln_retsel_bin_devs = 2,
                      ln_srvsel_devs = 2, ln_srvsel_bin_devs = 2, ln_fish_q_devs = 2, ln_srv_q_devs = 2,
                      ln_growth_devs = 3, ln_growth_semipar_devs = 3, move_devs = 4, ln_NAA = 3)
data_year_dims <- list(map_ln_RecDevs = 3, map_ln_fishsel_devs = 2, map_ln_fishsel_bin_devs = 2, map_ln_retsel_devs = 2, map_ln_retsel_bin_devs = 2,
                       map_ln_srvsel_devs = 2, map_ln_srvsel_bin_devs = 2, map_ln_fish_q_devs = 2, map_ln_srv_q_devs = 2,
                       map_ln_growth_devs = 3, map_ln_growth_semipar_devs = 3, map_move_devs = 4,
                       fish_sel_bicubic_Wyr = 2, ret_sel_bicubic_Wyr = 2, srv_sel_bicubic_Wyr = 2, move_year_block = 1)

# the slice of an array along one dim
along <- function(x, yr_dim, idx) {
  if(is.null(dim(x))) return(x[idx])
  i <- rep(list(TRUE), length(dim(x)))
  i[[yr_dim]] <- idx
  do.call("[", c(list(x), i, list(drop = FALSE)))
}

# a peel keeps the first n_keep fitted years of every extended array and every projected cell unchanged
expect_projected_cells_kept <- function(il, out, j) {
  n_yrs <- length(il$data$years)
  n_keep <- n_yrs - j
  for(nm in names(par_year_dims)) {
    par <- il$par[[nm]]
    if(is.null(par) || length(par) == 0) next
    yd <- par_year_dims[[nm]]
    n_past <- dim(par)[yd] - n_yrs
    if(n_past <= 0) next
    peeled <- out$retro_parameters[[nm]]
    expect_equal(dim(peeled)[yd], n_keep + n_past, info = nm)
    expect_equal(along(peeled, yd, n_keep + seq_len(n_past)), along(par, yd, n_yrs + seq_len(n_past)), info = nm)
    expect_equal(along(peeled, yd, seq_len(n_keep)), along(par, yd, seq_len(n_keep)), info = nm)
    if(!is.null(il$map[[nm]])) expect_equal(length(out$retro_mapping[[nm]]), length(peeled), info = nm)
  } # end nm loop
  for(nm in names(data_year_dims)) {
    x <- il$data[[nm]]
    if(is.null(x) || length(x) == 0) next
    yd <- data_year_dims[[nm]]
    n_total <- if(is.null(dim(x))) length(x) else dim(x)[yd]
    n_past <- n_total - n_yrs
    if(n_past <= 0) next
    peeled <- out$retro_data[[nm]]
    n_peeled <- if(is.null(dim(peeled))) length(peeled) else dim(peeled)[yd]
    expect_equal(n_peeled, n_keep + n_past, info = nm)
    expect_equal(along(peeled, yd, n_keep + seq_len(n_past)), along(x, yd, n_yrs + seq_len(n_past)), info = nm)
  } # end nm loop
}

test_that("a dusky peel keeps the projected cells, moves the recruitment window, and fits with the projection on", {

  il <- Setup_Mod_Projection(with_devs(build_goa_dusky_input)(sgl_rg_dusky_data), n_proj_yrs = n_proj_yrs, f_opt = "terminal")
  n_yrs <- length(il$data$years)
  j <- 2
  out <- truncate_yr(j = j, data = il$data, parameters = il$par, mapping = il$map)
  expect_projected_cells_kept(il, out, j)

  # the mean recruitment window ends one year before the new terminal year, as it did before the peel
  expect_equal(max(il$data$proj_calc_rec_yrs), n_yrs - 1)
  expect_equal(max(out$retro_data$proj_calc_rec_yrs), n_yrs - j - 1)
  expect_equal(out$retro_data$n_proj_yrs_devs, n_ext)

  fit <- fit_model(out$retro_data, out$retro_parameters, out$retro_mapping, newton_loops = 1, silent = TRUE)
  expect_equal(length(fit$rep$proj_SSB[1,1,]), n_proj_yrs)
  expect_equal(fit$rep$proj_SSB[1,1,1], fit$rep$SSB[1,1,n_yrs - j], tolerance = 1e-8) # year one replays the new terminal year
  expect_true(all(is.finite(fit$rep$proj_SSB)))

})

test_that("a two-region toy with movement deviations peels and fits with the projection on", {

  dims2 <- list(n_regions = 2, n_sexes = 1, n_fish_fleets = 1, n_seas = 1, n_yrs = 22, n_ages = 8)
  il <- Setup_Mod_Projection(with_devs(sweep_input)(dims = dims2), n_proj_yrs = n_proj_yrs, f_opt = "terminal")
  j <- 2
  out <- truncate_yr(j = j, data = il$data, parameters = il$par, mapping = il$map)
  expect_projected_cells_kept(il, out, j)
  expect_equal(dim(out$retro_parameters$move_devs)[4], dims2$n_yrs - j + n_ext)

  fit <- fit_model(out$retro_data, out$retro_parameters, out$retro_mapping, newton_loops = 1, silent = TRUE)
  expect_equal(dim(fit$rep$proj_SSB)[3], n_proj_yrs)
  expect_true(all(is.finite(fit$rep$proj_SSB)))

})

test_that("the Pacific cod bridge with growth deviations peels and its objective evaluates with the projection on", {

  il <- suppressWarnings(seed_ebs_pcod_mle(suppressMessages(with_devs(build_ebs_pcod_input)(sgl_rg_ebs_pcod_data)), sgl_rg_ebs_pcod_data))
  il$data$growth_tv_type <- 0 # curve mode builds every year up front, which is what projected years need
  il <- Setup_Mod_Projection(il, n_proj_yrs = n_proj_yrs, f_opt = "terminal", growth_opt = "propagate", sel_opt = "propagate")
  j <- 3
  out <- truncate_yr(j = j, data = il$data, parameters = il$par, mapping = il$map)
  expect_projected_cells_kept(il, out, j)

  obj <- fit_model(out$retro_data, out$retro_parameters, out$retro_mapping, do_optim = FALSE, silent = TRUE)
  expect_equal(dim(obj$rep$WAA)[3], length(il$data$years) - j + n_ext)
  expect_equal(dim(obj$rep$proj_SSB)[3], n_proj_yrs)
  expect_true(all(is.finite(obj$rep$proj_SSB)))

})
