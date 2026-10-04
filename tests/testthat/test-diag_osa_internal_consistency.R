# Internal one-step-ahead compositions evaluate the same likelihood as the ordinary path, so on whole-fish data
# the two objectives differ by a constant and share their gradient, for every composition source at once. Also
# checks residuals on recorded length bins and age-at-length rows are standard normal and labeled by length.

library(SPoRC)
library(testthat)
data("sgl_rg_ebs_pcod_data")
data("mlt_rg_goa_rex_data")

test_that("the internal OSA path is the ordinary composition likelihood up to a constant", {

  # pcod: lengths on recorded bins through LenBinMap and survey ages through ageing error, joint, aggregated
  # and Dirichlet-multinomial
  pcod <- seed_ebs_pcod_mle(suppressWarnings(suppressMessages(build_ebs_pcod_input(sgl_rg_ebs_pcod_data))), sgl_rg_ebs_pcod_data)
  aggregated <- pcod$data
  aggregated$FishLenComps_Type[] <- 0
  dirichlet <- pcod$data
  dirichlet$SrvLenComps_LikeType[] <- 1
  for(data in list(pcod$data, aggregated, dirichlet)) {
    gaps <- osa_path_gaps(data, pcod$par, pcod$map)
    expect_equal(gaps[[1]][["objective"]], gaps[[2]][["objective"]], tolerance = 1e-8) # a constant, whatever the parameters
    expect_lt(max(gaps[[1]][["gradient"]], gaps[[2]][["gradient"]]), 1e-7)
  } # end data loop

  # the iid and AR1 logistic normals, on survey lengths with fish in every cell. with empty cells the two paths
  # differ, since the ordinary path drops those bins and the OSA path keeps them at the added constant
  no_empty <- pcod$data
  for(y in which(no_empty$UseSrvLenComps[1,,1,1] == 1)) {
    p <- no_empty$ObsSrvLenComps[1,y,1,,1,1] / sum(no_empty$ObsSrvLenComps[1,y,1,,1,1])
    no_empty$ObsSrvLenComps[1,y,1,,1,1] <- (p + 0.01) / sum(p + 0.01)
  } # end y loop
  for(like_type in c(2, 3)) {
    logistic_normal <- no_empty
    logistic_normal$SrvLenComps_LikeType[] <- like_type
    # the OSA path adds the tiny constant to the expected proportions and the ordinary path does not, which moves
    # the gap a little in bins expected near empty
    gaps <- osa_path_gaps(logistic_normal, pcod$par, pcod$map)
    expect_equal(gaps[[1]][["objective"]], gaps[[2]][["objective"]], tolerance = 1e-5)
    expect_lt(max(gaps[[1]][["gradient"]], gaps[[2]][["gradient"]]), 1e-5)
  } # end like_type loop

  # rex: two regions and two sexes, ages and lengths joint across sexes, and age-at-length split by sex on
  # rows merging pairs of model length bins
  rex <- rex_caal_input()
  merged <- caal_on_rows(rex$data, bin_pairs(length(rex$data$lens)))
  gaps <- osa_path_gaps(merged, rex$par, rex$map)
  expect_equal(gaps[[1]][["objective"]], gaps[[2]][["objective"]], tolerance = 1e-8)
  expect_lt(max(gaps[[1]][["gradient"]], gaps[[2]][["gradient"]]), 1e-7)

})

test_that("recorded length bins and age-at-length rows are labeled by the model lengths they cover", {

  data <- list(lens = c(12.5, 17.5, 22.5, 27.5, 32.5), LenBinMap = bin_pairs(5), CAAL_LenBinMap = bin_pairs(5))
  expect_equal(obs_len_labels(data), c(12.5, 22.5, 32.5))
  expect_equal(caal_row_lens(data), c(12.5, 22.5, 32.5))
  data$LenBinMap <- data$CAAL_LenBinMap <- NULL
  expect_equal(obs_len_labels(data), data$lens)
  expect_equal(caal_row_lens(data), data$lens)

})

test_that("residuals on recorded length bins and age-at-length rows are standard normal under the generating model", {

  skip_on_cran()

  # the CAAL self test's operating model with lengths recorded on pairs of bins and ages read on pairs of bins
  pairs <- bin_pairs(caal_cfg$n_lens)
  om <- caal_make_om(len_bin_map = pairs, caal_len_bin_map = pairs)
  sim_data <- simulation_data_to_SPoRC(sim_env = om, y = caal_cfg$n_yrs, sim = 1)
  expect_equal(dim(sim_data$ObsFishLenComps)[4], ncol(pairs))
  expect_equal(dim(sim_data$ObsFish_caal)[4], ncol(pairs))
  input <- caal_build_input(sim_data, osa = TRUE, len_bin_map = pairs, caal_len_bin_map = pairs)
  fit <- suppressWarnings(fit_model(input$data, input$par, input$map, random = NULL, silent = TRUE))

  # lengths, labeled by the first model length in each recorded bin without being told. the last bin of a
  # composition is set by the others, so it has no residual
  pair_lens <- (caal_cfg$len_lower + 2.5)[seq(1, caal_cfg$n_lens, by = 2)]
  len_res <- get_osa(model = fit, data = fit$data, comp_source = "FishLen")$res
  expect_setequal(unique(len_res$index), head(pair_lens, -1))
  r_len <- len_res$resid[is.finite(len_res$resid)]
  expect_gt(length(r_len), 100)
  expect_lt(abs(stats::sd(r_len) - 1), 0.15)

  # age-at-length rows, each labeled by its first model length, with the same spread as the lengths
  caal_res <- get_osa(model = fit, data = fit$data, comp_source = "Fish_caal", bins = 1:caal_cfg$n_ages, bin_label = "Age")$res
  expect_setequal(unique(caal_res$len), pair_lens)
  r_caal <- caal_res$resid[is.finite(caal_res$resid)]
  expect_gt(length(r_caal), 500)
  expect_lt(abs(stats::sd(r_caal) - 1), 0.15)
  expect_lt(abs(mean(r_caal)), abs(mean(r_len)) + 0.25) # discreteness leaves the same small offset as in the CAAL self test

})
