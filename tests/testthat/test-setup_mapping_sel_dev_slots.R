library(SPoRC)
library(testthat)

# Under iid/rw time variation a parametric fleet's deviation slots are its
# parameters, so slots beyond what the form reads must stay unmapped: a
# time-varying logistic fleet carries two deviation parameters per year, not
# n_ages. Dead slots would never be read by the model, leaving zero-gradient
# parameters that pollute the count and make the Hessian singular.

test_that("deviation slots beyond a fleet's selectivity form stay unmapped", {

  n_yrs <- 6; n_ages <- 8

  input_list <- Setup_Mod_Dim(years = 1:n_yrs, ages = 1:n_ages, lens = NULL,
                              n_regions = 1, n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = 2,
                              n_pop = 1, natal_region = 1, verbose = FALSE)
  input_list <- Setup_Mod_Rec(input_list = input_list, do_rec_bias_ramp = 0, sigmaR_switch = 1,
                              ln_sigmaR = array(log(0.4), c(2, 1, 1)), rec_model = "mean_rec",
                              sigmaR_spec = "fix", init_age_strc = 1, equil_init_age_strc = 2,
                              ln_global_R0 = log(5))
  waa <- array(1, dim = c(1, 1, n_yrs, 1, n_ages, 1))
  suppressWarnings(input_list <- Setup_Mod_Biologicals(
    input_list = input_list, WAA = waa, MatAA = waa,
    WAA_fish = array(1, dim = c(1, 1, n_yrs, 1, n_ages, 1, 1)),
    WAA_srv = array(1, dim = c(1, 1, n_yrs, 1, n_ages, 1, 2)),
    fit_lengths = 0, AgeingError = array(rep(diag(n_ages), each = n_yrs), dim = c(n_yrs, n_ages, n_ages)),
    M_spec = "fix", Fixed_natmort = array(0.3, dim = c(1, 1, n_yrs, n_ages, 1))))
  input_list <- Setup_Mod_Tagging(input_list = input_list, use_conv_fish_tagging = 0)
  input_list <- Setup_Mod_Movement(input_list = input_list, use_fixed_movement = 1,
                                   Fixed_Movement = NA, do_recruits_move = 0)
  suppressWarnings(input_list <- Setup_Mod_Catch_and_F(
    input_list = input_list, ObsCatch = array(5, dim = c(1, n_yrs, 1, 1)),
    UseCatch = array(1, dim = c(1, n_yrs, 1, 1)), Use_F_pen = 1,
    sigmaC_spec = "fix", ln_sigmaC = array(log(0.05), dim = c(1, n_yrs, 1, 1)),
    ln_sigmaF = array(log(1), dim = c(1, 1, 1))))
  suppressWarnings(input_list <- Setup_Mod_FishIdx_and_Comps(
    input_list = input_list,
    ObsFishIdx = array(NA_real_, dim = c(1, n_yrs, 1, 1)),
    ObsFishIdx_SE = array(0.2, dim = c(1, n_yrs, 1, 1)),
    UseFishIdx = array(0, dim = c(1, n_yrs, 1, 1)),
    ObsFishAgeComps = array(10, dim = c(1, n_yrs, 1, n_ages, 1, 1)),
    UseFishAgeComps = array(1, dim = c(1, n_yrs, 1, 1)),
    ISS_FishAgeComps = array(50, dim = c(1, n_yrs, 1, 1, 1)),
    ObsFishLenComps = NULL, UseFishLenComps = array(0, dim = c(1, n_yrs, 1, 1)),
    ISS_FishLenComps = NULL,
    fish_idx_type = "biom", FishAgeComps_LikeType = "Multinomial",
    FishLenComps_LikeType = "none", FishAgeComps_Type = "agg_Year_1-terminal_Fleet_1",
    FishLenComps_Type = "none_Year_1-terminal_Fleet_1"))
  suppressWarnings(input_list <- Setup_Mod_SrvIdx_and_Comps(
    input_list = input_list,
    ObsSrvIdx = array(100, dim = c(1, n_yrs, 1, 2)),
    ObsSrvIdx_SE = array(0.2, dim = c(1, n_yrs, 1, 2)),
    UseSrvIdx = array(1, dim = c(1, n_yrs, 1, 2)),
    ObsSrvAgeComps = array(10, dim = c(1, n_yrs, 1, n_ages, 1, 2)),
    UseSrvAgeComps = array(1, dim = c(1, n_yrs, 1, 2)),
    ISS_SrvAgeComps = array(50, dim = c(1, n_yrs, 1, 1, 2)),
    ObsSrvLenComps = NULL, UseSrvLenComps = array(0, dim = c(1, n_yrs, 1, 2)),
    ISS_SrvLenComps = NULL,
    srv_idx_type = c("biom", "biom"), SrvAgeComps_LikeType = c("Multinomial", "Multinomial"),
    SrvLenComps_LikeType = c("none", "none"),
    SrvAgeComps_Type = c("agg_Year_1-terminal_Fleet_1", "agg_Year_1-terminal_Fleet_2"),
    SrvLenComps_Type = c("none_Year_1-terminal_Fleet_1", "none_Year_1-terminal_Fleet_2")))
  input_list <- Setup_Mod_Fishsel_and_Q(input_list = input_list, fish_sel_model = "logist1_Fleet_1",
                                        fish_fixed_sel_pars_spec = "est_all", fish_q_spec = "est_all",
                                        use_fixed_ret_sel = 1)
  input_list <- Setup_Mod_Srvsel_and_Q(input_list = input_list,
                                       srv_sel_model = c("logist1_Fleet_1", "nonparlog_Fleet_2"),
                                       cont_tv_srv_sel = c("iid_Fleet_1", "rw_Fleet_2"),
                                       srv_sel_devs_spec = c("est_all", "est_all"),
                                       srvsel_pe_pars_spec = c("est_all", "est_all"),
                                       srv_sel_nonpar_est_bins = list(NULL, 1:n_ages),
                                       srv_fixed_sel_pars_spec = c("est_all", "est_all"),
                                       srv_q_spec = c("est_all", "est_all"))

  m <- array(as.numeric(input_list$map$ln_srvsel_devs), dim = dim(input_list$par$ln_srvsel_devs))

  test_that("a time-varying logistic fleet carries exactly its two parameter slots", {
    active <- apply(!is.na(m[1, , , 1, 1]), 2, any)
    expect_equal(active, c(TRUE, TRUE, rep(FALSE, n_ages - 2)))
  })

  test_that("a non-parametric fleet's deviations are indexed by bin, all active", {
    expect_true(all(!is.na(m[1, , , 1, 2])))
  })

  test_that("the data copy the penalty reads matches the factor map", {
    expect_identical(is.na(input_list$data$map_ln_srvsel_devs), is.na(m))
  })

})
