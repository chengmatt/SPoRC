# The written out multinomial with nothing added to the proportions: a bin observed empty adds zero rather than
# 0 * log(0), so addtocomp = 0 gives the multinomial's likelihood instead of NaN, and nothing else changes.

test_that("an empty bin adds zero, and a bin with fish is unchanged", {
  obs <- c(0.5, 0.5, 0)
  pred <- c(0.4, 0.6, 0)
  expect_equal(comp_mltnml_term(obs, pred, 0), sum(obs[1:2] * log(pred[1:2])))
  expect_true(is.finite(comp_mltnml_term(obs, obs, 0))) # the offset term, which was 0 * log(0)
  expect_identical(comp_mltnml_term(obs + 1e-3, pred, 1e-3), sum((obs + 1e-3) * log(pred + 1e-3))) # with the constant, as before
})

test_that("with addtocomp = 0 the composition likelihood is the multinomial's, empty bins included", {
  counts <- c(30, 0, 50, 20, 0)
  pred <- c(0.3, 0.05, 0.4, 0.25, 0) # the last bin cannot hold fish
  for(comp_const_obs in c(0, 1)) {
    nLL <- Get_Comp_Likelihoods(comp_const_obs = comp_const_obs, Exp = pred, Obs = counts, ISS = 100, Wt_Mltnml = 1,
                                Comp_Type = 0, Likelihood_Type = 0, ln_theta = 0, ln_theta_agg = 0, LN_corr_pars = array(0, dim = c(1, 1, 3)),
                                LN_corr_pars_agg = 0, n_regions = 1, n_sexes = 1, age_or_len = 1, AgeingError = NA, use = 1,
                                n_model_bins = 5, n_obs_bins = 5, addtocomp = 0, comp_bins = NULL)
    props <- counts / sum(counts)
    seen <- counts > 0
    expect_equal(sum(nLL), -100 * sum(props[seen] * log(pred[seen] / props[seen])), info = comp_const_obs)
  } # end comp_const_obs loop
})

test_that("a fit with addtocomp = 0 and empty age bins has a finite objective and gradient", {
  set.seed(11)
  il <- naaom_build_em(naaom_om_data(naaom_make_om(NAA_re = "none")), NAA_re = "none")
  il$data$addtocomp <- 0
  expect_true(any(il$data$ObsFishAgeComps == 0) || any(il$data$ObsSrvAgeComps == 0)) # there are empty bins to test
  obj <- suppressWarnings(fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE))
  expect_true(is.finite(obj$fn(obj$par)))
  expect_true(all(is.finite(obj$gr(obj$par))))
})
