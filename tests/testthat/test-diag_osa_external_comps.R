library(SPoRC)
library(testthat)

test_that("OSA residuals are calibrated (SDNR ≈ 1) for all composition likelihoods", {
  skip_if_not_installed("compResidual")
  skip_if_not_installed("TMB")

  set.seed(42)
  n_years <- 50
  n_bins  <- 30
  true_p  <- rep(1/n_bins, n_bins)
  tol     <- 0.1  # SDNR within [0.9, 1.1]

  # --- Multinomial ---
  N_mn   <- 500
  obs_mn <- t(replicate(n_years, as.vector(rmultinom(1, N_mn, true_p))))
  res_mn <- run_external_comp_osa(
    obs         = obs_mn,
    exp         = matrix(true_p, nrow = n_years, ncol = n_bins, byrow = TRUE),
    N           = N_mn,
    fleet       = "test",
    index       = 1:n_bins,
    years       = 1:n_years,
    index_label = "age",
    comp_like   = 0
  )
  sdnr_mn <- sd(res_mn$res$resid, na.rm = TRUE)
  expect_equal(sdnr_mn, 1, tolerance = tol,
               label = paste("Multinomial SDNR:", round(sdnr_mn, 3)))

  # --- Dirichlet-multinomial ---
  N_dm   <- 500
  theta  <- 0.5
  obs_dm <- t(rdirM(n_years, N_dm, theta * N_dm * true_p))
  res_dm <- run_external_comp_osa(
    obs         = obs_dm,
    exp         = matrix(true_p, nrow = n_years, ncol = n_bins, byrow = TRUE),
    N           = N_dm,
    DM_theta    = theta,
    fleet       = "test",
    index       = 1:n_bins,
    years       = 1:n_years,
    index_label = "age",
    comp_like   = 1
  )
  sdnr_dm <- sd(res_dm$res$resid, na.rm = TRUE)
  expect_equal(sdnr_dm, 1, tolerance = tol,
               label = paste("Dir-Multinomial SDNR:", round(sdnr_dm, 3)))

  # --- Logistic-normal iid ---
  sigma_ln <- 0.1
  obs_ln2  <- t(replicate(n_years, rlogistnormal(true_p, sigma_ln, comp_like = 2, n_sexes = 1)))
  res_ln2  <- run_external_comp_osa(
    obs         = obs_ln2,
    exp         = matrix(true_p, nrow = n_years, ncol = n_bins, byrow = TRUE),
    LN_Sigma    = diag(n_bins) * sigma_ln^2,  # K x K; run_osa strips last row/col
    fleet       = "test",
    index       = 1:n_bins,
    years       = 1:n_years,
    index_label = "age",
    comp_like   = 2
  )
  sdnr_ln2 <- sd(res_ln2$res$resid, na.rm = TRUE)
  expect_equal(sdnr_ln2, 1, tolerance = tol,
               label = paste("LN-iid SDNR:", round(sdnr_ln2, 3)))

  # --- Logistic-normal AR1 ---
  rho_ln3  <- 0.7
  obs_ln3  <- t(replicate(n_years, rlogistnormal(true_p, c(sigma_ln, rho_ln3), comp_like = 3, n_sexes = 1)))
  res_ln3  <- run_external_comp_osa(
    obs         = obs_ln3,
    exp         = matrix(true_p, nrow = n_years, ncol = n_bins, byrow = TRUE),
    LN_Sigma    = get_AR1_CorrMat(n_bins, rho_ln3) * (sigma_ln^2 / (1 - rho_ln3^2)),  # K x K
    fleet       = "test",
    index       = 1:n_bins,
    years       = 1:n_years,
    index_label = "age",
    comp_like   = 3
  )
  sdnr_ln3 <- sd(res_ln3$res$resid, na.rm = TRUE)
  expect_equal(sdnr_ln3, 1, tolerance = tol,
               label = paste("LN-AR1 SDNR:", round(sdnr_ln3, 3)))

  # --- Logistic-normal AR1 x constant sex correlation (comp_like = 4) ---
  rho_sex  <- 0.5
  n_sexes  <- 2
  n_bins_sex <- n_bins * n_sexes  # total bins = bins per sex * n_sexes
  true_p_sex <- rep(1/n_bins_sex, n_bins_sex)

  obs_ln4 <- t(replicate(n_years, rlogistnormal(true_p_sex, c(sigma_ln, rho_ln3, rho_sex),
                                                comp_like = 4, n_sexes = n_sexes)))

  Sigma_ln4 <- kronecker(
    get_Constant_CorrMat(n_sexes, rho_sex),
    get_AR1_CorrMat(n_bins, rho_ln3)
  ) * (sigma_ln^2 / (1 - rho_ln3^2) / (1 - rho_sex^2))  # K*n_sexes x K*n_sexes

  res_ln4 <- run_external_comp_osa(
    obs         = obs_ln4,
    exp         = matrix(true_p_sex, nrow = n_years, ncol = n_bins_sex, byrow = TRUE),
    LN_Sigma    = Sigma_ln4,  # K x K (full); run_osa strips last row/col
    fleet       = "test",
    index       = 1:n_bins_sex,
    years       = 1:n_years,
    index_label = "age",
    comp_like   = 4
  )
  sdnr_ln4 <- sd(res_ln4$res$resid, na.rm = TRUE)
  expect_equal(sdnr_ln4, 1, tolerance = tol, label = paste("LN-AR1xSex SDNR:", round(sdnr_ln4, 3)))

})
