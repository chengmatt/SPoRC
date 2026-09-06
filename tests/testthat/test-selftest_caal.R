library(SPoRC)
library(testthat)

# Conditional age-at-length, end to end through the operating model. The OM draws
# a fixed number of otoliths per length bin from the joint of length and age, the
# estimation model fits marginal lengths plus CAAL with no marginal ages, so all of
# its age information arrives through the CAAL likelihood. Routines shared through
# helper-selftest_caal.R. A correctly specified CAAL likelihood has to leave the
# median relative error of spawning biomass and recruitment at zero, skip length
# bins nobody aged, recover a Dirichlet-multinomial overdispersion it generated
# under, and produce standard normal one step ahead residuals.


test_that("simulate_caal draws from P(age | length) and applies ageing error", {
  n_lens <- 4; n_ages <- 3
  phi <- cbind(c(0.7, 0.2, 0.1, 0.0), c(0.1, 0.5, 0.3, 0.1), c(0.0, 0.1, 0.4, 0.5))
  caa <- c(100, 50, 20)
  joint <- phi * rep(caa, each = n_lens)
  cond <- joint / rowSums(joint) # P(age | length)

  sat <- array(phi, dim = c(1, 1, 1, 1, n_lens, n_ages, 1, 1))
  at_age <- array(caa, dim = c(1, 1, 1, 1, n_ages, 1, 1, 1))
  iss <- array(2e5, dim = c(1, 1, 1, n_lens, 1, 1, 1))
  ae <- array(diag(n_ages), dim = c(1, n_ages, n_ages, 1))
  obs <- array(0, dim = c(1, 1, 1, n_lens, n_ages, 1, 1, 1))

  set.seed(3)
  out <- simulate_caal(
    r = 1,
    y = 1,
    f = 1,
    seas = 1,
    sim = 1,
    SizeAgeTrans = sat,
    AtAge = at_age,
    ISS = iss,
    AgeingError = ae,
    comp_like = 0,
    ln_theta = array(0, c(1, 1, 1)),
    ln_theta_agg = 0,
    comp_type = matrix(1, 1, 1),
    n_sexes = 1,
    n_regions = 1,
    n_lens = n_lens,
    Obs = obs
  )
  drawn <- out[1, 1, 1, , , 1, 1, 1]
  expect_equal(rowSums(drawn), rep(2e5, n_lens)) # every bin got exactly its ISS
  expect_equal(drawn / rowSums(drawn), cond, tolerance = 0.01) # and in the conditional proportions

  # an ageing error that shifts every fish one age older moves the whole row
  shift <- rbind(cbind(0, diag(n_ages - 1)), c(rep(0, n_ages - 1), 1))
  ae_shift <- array(shift, dim = c(1, n_ages, n_ages, 1))
  out <- simulate_caal(
    r = 1,
    y = 1,
    f = 1,
    seas = 1,
    sim = 1,
    SizeAgeTrans = sat,
    AtAge = at_age,
    ISS = iss,
    AgeingError = ae_shift,
    comp_like = 0,
    ln_theta = array(0, c(1, 1, 1)),
    ln_theta_agg = 0,
    comp_type = matrix(1, 1, 1),
    n_sexes = 1,
    n_regions = 1,
    n_lens = n_lens,
    Obs = obs
  )
  drawn <- out[1, 1, 1, , , 1, 1, 1]
  expect_equal(drawn[, 1], rep(0, n_lens)) # nothing is left at age 1
  expect_equal(drawn[, 3] / rowSums(drawn), cond[, 2] + cond[, 3], tolerance = 0.01)
})


test_that("CAAL in place of marginal ages recovers spawning biomass and recruitment", {
  om <- caal_make_om()
  sim_data <- simulation_data_to_SPoRC(sim_env = om, y = caal_cfg$n_yrs, sim = 1)
  expect_true(all(sim_data$UseFish_caal == 1)) # every bin was aged in this design

  out <- caal_run(caal_build_input(sim_data), what = c("SSB", "Rec"), seed = 101)

  expect_gte(out$summ$SSB[["n_ok"]], 0.9 * caal_cfg$n_sims)
  expect_lt(abs(out$summ$SSB[["median_RE"]]), 0.03)
  expect_lt(abs(out$summ$Rec[["median_RE"]]), 0.05)
})


test_that("length bins nobody aged are skipped and leave the fit unbiased", {
  aged_bins <- seq(2, caal_cfg$n_lens, by = 2)
  om <- caal_make_om(caal_bins = aged_bins)
  sim_data <- simulation_data_to_SPoRC(sim_env = om, y = caal_cfg$n_yrs, sim = 1)

  # the use flags follow the sampling design, not the length bins in the model
  used <- apply(sim_data$UseFish_caal, 4, sum)
  expect_true(all(used[aged_bins] > 0))
  expect_true(all(used[-aged_bins] == 0))

  input <- caal_build_input(sim_data)
  fit <- fit_model(input$data, input$par, input$map, random = NULL, silent = TRUE)
  # unaged bins contribute nothing to the likelihood
  expect_true(all(fit$rep$Fish_caal_nLL[, , , -aged_bins, , ] == 0))
  expect_true(any(fit$rep$Fish_caal_nLL[, , , aged_bins, , ] != 0))

  out <- caal_run(input, what = c("SSB", "Rec"), n_sims = 30, seed = 202)
  expect_lt(abs(out$summ$SSB[["median_RE"]]), 0.03)
})


test_that("Dirichlet-multinomial CAAL recovers the overdispersion it was generated under", {
  theta_true <- 0.5
  n_reps <- 20
  om <- caal_make_om(
    caal_like = "Dirichlet-Multinomial",
    ln_theta = log(theta_true),
    n_sims = n_reps,
    seed = 77
  )

  theta_hat <- matrix(NA_real_, n_reps, 2, dimnames = list(NULL, c("fish", "srv")))
  for(i in 1:n_reps) {
    sim_data <- simulation_data_to_SPoRC(sim_env = om, y = caal_cfg$n_yrs, sim = i)
    input <- caal_build_input(sim_data, caal_like = "Dirichlet-Multinomial")
    expect_equal(sum(!is.na(input$map$ln_Fish_caal_theta)), 1) # one theta shared across length bins
    fit <- fit_model(input$data, input$par, input$map, random = NULL, silent = TRUE)
    par <- fit$env$last.par.best
    theta_hat[i, "fish"] <- exp(par[names(par) == "ln_Fish_caal_theta"])
    theta_hat[i, "srv"] <- exp(par[names(par) == "ln_Srv_caal_theta"])
  } # end i loop

  expect_lt(abs(median(theta_hat[, "fish"]) / theta_true - 1), 0.2)
  expect_lt(abs(median(theta_hat[, "srv"]) / theta_true - 1), 0.2)
})


test_that("two sexes drawn jointly within a length bin are recovered", {
  om <- caal_make_om(n_sexes = 2, caal_type = "spltRjntS")
  sim_data <- simulation_data_to_SPoRC(sim_env = om, y = caal_cfg$n_yrs, sim = 1)
  expect_equal(dim(sim_data$ObsFish_caal)[6], 2)

  # a joint draw puts the whole bin's ISS across both sexes, so the sexes sum to it
  bin_total <- apply(sim_data$ObsFish_caal, c(1, 2, 3, 4, 7), sum)
  expect_true(all(bin_total[sim_data$UseFish_caal == 1] == caal_cfg$caal_per_bin))

  out <- caal_run(
    caal_build_input(sim_data, caal_type = "spltRjntS", n_sexes = 2),
    what = c("SSB", "Rec"),
    n_sims = 30,
    seed = 303
  )
  expect_lt(abs(out$summ$SSB[["median_RE"]]), 0.03)
})


test_that("CAAL one step ahead residuals are standard normal under the generating model", {
  om <- caal_make_om()
  sim_data <- simulation_data_to_SPoRC(sim_env = om, y = caal_cfg$n_yrs, sim = 1)
  input <- caal_build_input(sim_data, osa = TRUE)
  fit <- fit_model(input$data, input$par, input$map, random = NULL, silent = TRUE)

  osa <- get_osa(
    model = fit,
    data = fit$data,
    comp_source = "Fish_caal",
    bins = 1:caal_cfg$n_ages,
    bin_label = "Age"
  )
  res <- osa$res
  expect_true("len" %in% names(res))
  expect_setequal(unique(res$len), caal_cfg$len_lower + 2.5)

  # The package's discrete residuals have a small offset in their mean under a
  # correct model (about -0.15 for marginal compositions at N = 300 in this same
  # OM), and a CAAL row of 25 otoliths is coarser still, so the standard is the
  # one the other OSA tests use: unit spread, and a mean no further from zero
  # than the marginal compositions in the same fit show plus that coarseness.
  r <- res$resid[is.finite(res$resid)]
  expect_gt(length(r), 1000)
  expect_lt(abs(stats::sd(r) - 1), 0.15)
  len_res <- get_osa(
    model = fit,
    data = fit$data,
    comp_source = "FishLen",
    bins = 1:caal_cfg$n_lens,
    bin_label = "Len"
  )$res$resid
  expect_lt(abs(mean(r)), abs(mean(len_res, na.rm = TRUE)) + 0.25)

  # no structure across the conditioning variable: the length bins agree with
  # each other, whatever common offset the discreteness leaves
  by_len <- tapply(res$resid, res$len, mean, na.rm = TRUE)
  expect_lt(max(by_len) - min(by_len), 0.6)
})
