# Conditional age-at-length on rows coarser than the model's length bins (CAAL_LenBinMap), on GOA rex sole:
# a row's expected ages sum the bins it covers in the fit, its diagnostics and the operating model, and a
# perfect-data self test on rows merging pairs of bins recovers the operating model.

library(SPoRC)
library(testthat)
data("mlt_rg_goa_rex_data")

test_that("an identity map is the same model as none, and the map is checked", {

  input_list <- rex_caal_input()
  n_lens <- length(input_list$data$lens)
  base <- fit_model(input_list$data, input_list$par, input_list$map, do_optim = FALSE, silent = TRUE)
  ident_data <- input_list$data
  ident_data$CAAL_LenBinMap <- diag(n_lens)
  ident <- fit_model(ident_data, input_list$par, input_list$map, do_optim = FALSE, silent = TRUE)
  expect_equal(ident$fn(ident$par), base$fn(base$par))
  expect_equal(ident$gr(ident$par), base$gr(base$par))
  expect_equal(caal_row_lens(ident_data), input_list$data$lens)

  expect_error(check_caal_len_bin_map(diag(3), n_lens, "CAAL_LenBinMap"), "one row per model length bin")
  expect_error(check_caal_len_bin_map(0.5 * diag(n_lens), n_lens, "CAAL_LenBinMap"), "only 0 and 1")
  expect_error(check_caal_len_bin_map(cbind(diag(n_lens), 0), n_lens, "CAAL_LenBinMap"), "covering no model length bin")

})

test_that("a row's expected ages sum the bins it covers, and rows need not cover every bin", {

  input_list <- rex_caal_input()
  n_lens <- length(input_list$data$lens)

  # each row covers only the first bin of a pair, as rows given by their own lower edges can
  first_of_pair <- matrix(0, n_lens, ceiling(n_lens / 2))
  for(k in 1:ncol(first_of_pair)) first_of_pair[2 * k - 1, k] <- 1
  sparse_data <- caal_on_rows(input_list$data, first_of_pair)
  obj <- fit_model(sparse_data, input_list$par, input_list$map, do_optim = FALSE, silent = TRUE)
  expect_true(is.finite(obj$fn(obj$par)))
  expect_equal(dim(obj$rep$Srv_caal_nLL)[4], ncol(first_of_pair))
  expect_equal(caal_row_lens(sparse_data), input_list$data$lens[seq(1, n_lens, by = 2)])

  # rows merging pairs of bins: the expected row is the pair's numbers at length and age summed, then aged
  pairs <- matrix(0, n_lens, ceiling(n_lens / 2))
  for(l in 1:n_lens) pairs[l, ceiling(l / 2)] <- 1
  merged <- caal_on_rows(input_list$data, pairs)
  obj <- fit_model(merged, input_list$par, input_list$map, do_optim = FALSE, silent = TRUE)
  prop <- get_caal_prop(obj$data, obj$rep)
  y <- which(apply(merged$UseSrv_caal[1,,1,,1], 1, sum) > 0)[1]
  k <- which(merged$UseSrv_caal[1,y,1,,1] == 1)[1]
  by_hand <- as.vector(colSums(obj$rep$Srv_caal[1,1,y,1,pairs[,k] == 1,,1,1]) %*% obj$data$AgeingError_srv[y,,,1])
  expect_equal(prop$Pred_Srv_caal[1,y,1,k,,1,1], by_hand / sum(by_hand), tolerance = 1e-12)

})

test_that("a perfect-data self test on rows merging pairs of bins recovers the operating model", {

  skip_on_cran()
  input_list <- rex_caal_input()
  n_lens <- length(input_list$data$lens)
  pairs <- matrix(0, n_lens, ceiling(n_lens / 2))
  for(l in 1:n_lens) pairs[l, ceiling(l / 2)] <- 1
  merged <- caal_on_rows(input_list$data, pairs)
  fit <- suppressWarnings(fit_model(merged, input_list$par, input_list$map, do_optim = TRUE, newton_loops = 2, silent = TRUE))

  set.seed(21)
  sim_path <- tempfile(fileext = ".rds")
  st <- simulation_self_test(data = fit$data, parameters = fit$parameters, mapping = fit$mapping, random = NULL,
                             rep = fit$rep, sd_rep = list(par.fixed = fit$env$last.par.best), n_sims = 1, newton_loops = 2,
                             what = c("SSB", "Rec", "srv_q"), what_par = "ln_growth_pars", perfect_data = TRUE, output_path = sim_path)
  sim <- readRDS(sim_path)

  # the operating model draws on the same rows, and on perfect data its rows are the fit's expected ones
  expect_equal(dim(sim$ObsSrv_caal)[4], ncol(pairs))
  expected <- get_caal_prop(fit$data, fit$rep)$Pred_Srv_caal
  worst <- 0
  for(r in 1:2) for(y in seq_along(fit$data$years)) for(k in 1:ncol(pairs)) for(f in 1:2) for(s in 1:2) {
    simulated <- sim$ObsSrv_caal[r,y,1,k,,s,f,1]
    if(fit$data$UseSrv_caal[r,y,1,k,f] != 1 || sum(simulated) == 0) next
    worst <- max(worst, abs(simulated / sum(simulated) - expected[r,y,1,k,,s,f]))
  } # end r, y, k, f, s loops
  expect_lt(worst, 0.005)

  # and the refit recovers it. the first year and the last year classes are seen in few observations, so they
  # get a looser bound
  rel_err <- function(est, truth) max(abs(est / truth - 1))
  expect_lt(rel_err(st$SSB, st$truth$SSB), 0.01)
  expect_lt(rel_err(st$srv_q, st$truth$srv_q), 0.005)
  expect_lt(max(abs(st$ln_growth_pars - st$truth$ln_growth_pars)), 0.005)
  rec_err <- abs(as.vector(st$Rec[1,1,,1]) / as.vector(st$truth$Rec[1,1,,1]) - 1)
  n_yrs <- length(rec_err)
  expect_lt(max(rec_err[2:(n_yrs - 4)]), 0.01)
  expect_lt(max(rec_err[c(1, (n_yrs - 3):n_yrs)]), 0.05)

})
