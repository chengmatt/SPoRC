# plot_resids draws its second panel against the residual's bin labels, from get_osa's bins and bin_label
# arguments. Left out, the label columns never reached the frame; they are filled from the data now.

library(SPoRC)
library(testthat)
data("mlt_rg_goa_rex_data")

test_that("get_osa labels its bins without being told, and both residual plots draw", {

  dat <- mlt_rg_goa_rex_data
  input_list <- suppressWarnings(suppressMessages(build_goa_rex_input(dat)))
  input_list$data$do_internal_comp_osa <- TRUE
  fit <- fit_model(input_list$data, input_list$par, input_list$map, do_optim = FALSE, silent = TRUE)

  # a marginal length composition and a conditional age-at-length data source, neither
  # given bins
  for(src in c("SrvLen", "Srv_caal")) {

    osa <- get_osa(model = fit, data = input_list$data, comp_source = src)
    expect_true(all(c("index", "index_label") %in% names(osa$res)))
    expect_equal(unique(osa$res$index_label), if(src == "SrvLen") "Length" else "Age")

    # every plot has to survive being drawn, not just built
    plots <- plot_resids(osa)
    expect_gt(length(plots), 1)
    # a handful of non-finite residuals is ordinary, so only the draw itself is checked
    for(pl in plots) expect_no_error(suppressWarnings(ggplot2::ggplot_build(pl)))

  } # end src loop

  # lengths take the model's length bins, ages the bin index when ageing error
  # leaves fewer observed bins than the model has
  len_bins <- osa_default_bins(input_list$data, "SrvLen")
  expect_equal(len_bins$bins, input_list$data$lens)
  age_bins <- osa_default_bins(input_list$data, "Srv_caal")
  expect_equal(length(age_bins$bins), dim(input_list$data$ObsSrv_caal)[5])
  expect_lt(length(age_bins$bins), length(input_list$data$ages))

  # an explicit choice still wins
  osa_explicit <- get_osa(
    model = fit,
    data = input_list$data,
    comp_source = "SrvLen",
    bins = seq_along(input_list$data$lens),
    bin_label = "Bin"
  )
  expect_equal(unique(osa_explicit$res$index_label), "Bin")
})


test_that("the external OSA path takes a plain year vector and skips regions with no data", {

  skip_if_not_installed("compResidual")

  dat <- mlt_rg_goa_rex_data
  input_list <- suppressWarnings(suppressMessages(build_goa_rex_input(dat)))
  fit <- fit_model(input_list$data, input_list$par, input_list$map, do_optim = FALSE, silent = TRUE)
  cp <- get_comp_prop(
    input_list$data,
    fit$rep,
    age_labels = 1:dim(input_list$data$ObsSrvAgeComps)[4],
    len_labels = input_list$data$lens,
    year_labels = input_list$data$years
  )

  n_reg <- input_list$data$n_regions
  n_sex <- input_list$data$n_sexes
  n_len <- length(input_list$data$lens)
  ISS <- input_list$data$ISS_SrvLenComps
  yrs <- which(input_list$data$UseSrvLenComps[1, , 1, 1] == 1)

  # survey fleet one samples region one only, so the split types have a region
  # with no years at all, which used to collapse the slice
  expect_equal(sum(input_list$data$UseSrvLenComps[2, , 1, 1]), 0)

  # every composition type off the same plain year vector, and N at the model's
  # full year dimension in every case: comp_type 0 used to want N pre-sliced to
  # `years` while 1 and 2 wanted the full array, so a caller had to know which
  # branch it was calling into. All three now index N by years/years_by_region
  # themselves, the same way they already index obs_mat and exp_mat.
  extra <- list(
    list(N = ISS[1, , 1, 1, 1]),
    list(
      N = ISS[, , 1, , 1],
      DM_theta = array(0, c(n_reg, n_sex)),
      LN_Sigma = array(0.1, c(n_reg, n_len, n_len, n_sex))
    ),
    list(
      N = ISS[, , 1, 1, 1],
      DM_theta = rep(0, n_reg),
      LN_Sigma = array(0.1, c(n_reg, n_len, n_len))
    )
  )

  for(ct in c(0, 1, 2)) {
    # compResidual reports NaNs of its own on cells it cannot evaluate
    out <- suppressWarnings(do.call(get_osa, c(list(
      obs_mat = cp$Obs_SrvLen_mat,
      exp_mat = cp$Pred_SrvLen_mat,
      years = yrs,
      seas = 1,
      fleet = 1,
      bins = input_list$data$lens,
      comp_type = ct,
      comp_like = 0,
      bin_label = "Length"
    ), extra[[ct + 1]])))
    expect_gt(nrow(out$res), 0)
    expect_true(all(c("year", "index", "resid") %in% names(out$res)))
    # nothing from the region that was never sampled
    if(ct > 0) expect_false(2 %in% unique(out$res$region))
  } # end ct loop

  # N is read at the years actually used, not the position in a pre-sliced
  # vector: feeding a vector already sliced to length(yrs) would misalign against
  # the aggregated years index and give the wrong sample size to every year
  out_agg <- suppressWarnings(get_osa(
    obs_mat = cp$Obs_SrvLen_mat,
    exp_mat = cp$Pred_SrvLen_mat,
    N = ISS[1, , 1, 1, 1],
    years = yrs,
    seas = 1,
    fleet = 1,
    bins = input_list$data$lens,
    comp_type = 0,
    comp_like = 0,
    bin_label = "Length"
  ))
  # the residual frame's year column holds the year INDICES passed in, not
  # calendar years, since run_external_comp_osa names its dimnames from `years`
  # verbatim
  expect_equal(sort(unique(out_agg$res$year)), sort(yrs))

  # and a per-region list still works, which is what the split types took before
  yr_list <- lapply(1:n_reg, function(r) which(input_list$data$UseSrvLenComps[r, , 1, 1] == 1))
  out_list <- suppressWarnings(get_osa(
    obs_mat = cp$Obs_SrvLen_mat,
    exp_mat = cp$Pred_SrvLen_mat,
    N = ISS[, , 1, 1, 1],
    DM_theta = rep(0, n_reg),
    LN_Sigma = array(0.1, c(n_reg, n_len, n_len)),
    years = yr_list,
    seas = 1,
    fleet = 1,
    bins = input_list$data$lens,
    comp_type = 2,
    comp_like = 0,
    bin_label = "Length"
  ))
  expect_gt(nrow(out_list$res), 0)
})


test_that("CAAL residual bubbles separate by the length bin the fish were aged from", {

  dat <- mlt_rg_goa_rex_data
  input_list <- suppressWarnings(suppressMessages(build_goa_rex_input(dat)))
  input_list$data$do_internal_comp_osa <- TRUE
  fit <- fit_model(input_list$data, input_list$par, input_list$map, do_optim = FALSE, silent = TRUE)

  osa <- get_osa(model = fit, data = input_list$data, comp_source = "Srv_caal")

  # without the length bin in the facets every bin lands on the same year and age
  stacked <- table(osa$res$year, osa$res$index)
  expect_gt(stats::median(stacked[stacked > 0]), 20)

  # the bubble plot is the second element, and it now facets on len
  plots <- plot_resids(osa)
  facet_vars <- unlist(lapply(plots[[2]]$facet$params[c("rows", "cols")], names))
  expect_true("len" %in% facet_vars)

  # a marginal composition has no length bin to facet on, so it is untouched
  osa_len <- get_osa(model = fit, data = input_list$data, comp_source = "SrvLen")
  facet_len <- unlist(lapply(plot_resids(osa_len)[[2]]$facet$params[c("rows", "cols")], names))
  expect_false("len" %in% facet_len)
})
