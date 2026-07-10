library(SPoRC)
library(testthat)

test_that("plot_resids() facets by fleet for every comp_type when multiple fleets are present", {

  mk_comp <- function(comp_type, n_fleet) {
    df <- expand.grid(fleet = 1:n_fleet, region = 1:2, sex = 1:2, year = 1:10, index = 1:5)
    df$resid <- rnorm(nrow(df))
    df$index_label <- "Age"
    df$comp_type <- comp_type
    df
  }
  mk_tag <- function(n_fleet) {
    df <- expand.grid(fleet = 1:n_fleet, recovery_year = 1:10, years_at_liberty = 1:3)
    df$resid <- rnorm(nrow(df))
    df$comp_type <- "Tag"
    df
  }

  for (ct in c("Aggregated", "SpltR_SpltS", "SpltR_JntS")) {
    for (n_fleet in c(1, 2)) {
      res <- list(res = mk_comp(ct, n_fleet))
      p <- plot_resids(res)
      expect_type(p, "list")
      expect_length(p, 2)
      expect_s3_class(p[[1]], "ggplot")
      expect_s3_class(p[[2]], "ggplot")
    }
  }

  for (n_fleet in c(1, 2)) {
    res <- list(res = mk_tag(n_fleet))
    p <- plot_resids(res)
    expect_s3_class(p[[1]], "ggplot")
  }
})

# --- Integration test: 2 fish fleets, 2 srv fleets, 2 regions, 2 sexes -------
test_that("get_osa(model = ...) + plot_resids() work for a multi-fleet, multi-sex, multi-region model", {

  set.seed(42)
  n_yrs <- 40; n_regions <- 2; n_sexes <- 2; n_fish_fleets <- 2; n_srv_fleets <- 2
  n_ages <- 8; n_pop <- 1; n_seas <- 1

  sim_list <- Setup_Sim_Dim(
    n_sims = 1, n_yrs = n_yrs, n_regions = n_regions, n_ages = n_ages, n_lens = NULL,
    n_sexes = n_sexes, n_fish_fleets = n_fish_fleets, n_srv_fleets = n_srv_fleets,
    n_seas = n_seas, n_pop = n_pop
  )
  sim_list <- Setup_Sim_Containers(sim_list)

  rep_all <- function(vec, sim_list, n_fleets = NULL) {
    base <- sim_list$n_pop * sim_list$n_regions * sim_list$n_yrs * sim_list$n_seas
    times <- sim_list$n_sexes * (if (is.null(n_fleets)) 1 else n_fleets)
    rep(vec, each = base, times = times)
  }

  fish_sel_vec <- 1 / (1 + exp(-3 * ((1:n_ages) - 2)))
  ret_sel_vec  <- 1 / (1 + exp(-3 * ((1:n_ages) - 4)))
  srv_sel_vec  <- 1 / (1 + exp(-1 * ((1:n_ages) - 3)))

  fish_sel_input <- array(rep_all(fish_sel_vec, sim_list, n_fish_fleets),
                          dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_fish_fleets))
  ret_sel_input <- array(rep_all(ret_sel_vec, sim_list, n_fish_fleets),
                         dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_fish_fleets))
  srv_sel_input <- array(rep_all(srv_sel_vec, sim_list, n_srv_fleets),
                         dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_srv_fleets))

  dmr_input <- array(0.4, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets, sim_list$n_sims))

  # Ramp F over time so the data have real depletion contrast (a flat F gives
  # zero contrast to separate R0 from selectivity x catchability
  Fmort_vec <- c(seq(0.04, 0.16, length.out = round(n_yrs * 0.7)),
                seq(0.16, 0.10, length.out = n_yrs - round(n_yrs * 0.7)))
  Fmort_input <- array(rep(Fmort_vec, each = n_regions, times = n_seas * n_fish_fleets * sim_list$n_sims),
                       dim = c(n_regions, n_yrs, n_seas, n_fish_fleets, sim_list$n_sims))
  FishAgeComps_Type_sim <- array(1, dim = c(n_yrs, n_fish_fleets))
  FishAgeComps_discard_Type_sim <- array(0, dim = c(n_yrs, n_fish_fleets))

  sim_list <- Setup_Sim_Fishing(
    sim_list = sim_list,
    fish_sel_input = replicate(sim_list$n_sims, fish_sel_input),
    ret_sel_input = replicate(sim_list$n_sims, ret_sel_input),
    dmr_input = dmr_input,
    Fmort_input = Fmort_input,
    FishAgeComps_Type = FishAgeComps_Type_sim,
    FishAgeComps_discard_Type = FishAgeComps_discard_Type_sim
  )

  sim_list <- Setup_Sim_Survey(
    sim_list = sim_list,
    srv_sel_input = replicate(sim_list$n_sims, srv_sel_input)
  )

  waa_vec <- 5 / (1 + exp(-3 * ((1:n_ages) - 3)))
  mat_vec <- 1 / (1 + exp(-3 * ((1:n_ages) - 3)))

  sim_list <- suppressWarnings(Setup_Sim_Biologicals(
    sim_list = sim_list,
    natmort_input = array(0.25, dim = c(n_pop, n_regions, n_yrs, n_ages, n_sexes, sim_list$n_sims)),
    WAA_input = array(rep_all(waa_vec, sim_list), dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, sim_list$n_sims)),
    WAA_fish_input = array(rep_all(waa_vec, sim_list, n_fish_fleets), dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_fish_fleets, sim_list$n_sims)),
    WAA_srv_input = array(rep_all(waa_vec, sim_list, n_srv_fleets), dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_srv_fleets, sim_list$n_sims)),
    MatAA_input = array(rep_all(mat_vec, sim_list), dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, sim_list$n_sims))
  ))

  sim_list <- Setup_Sim_Tagging(sim_list = sim_list, use_conv_fish_tagging = 0)

  stay_prob <- 0.7
  sim_list$Movement <- array(0, dim = c(n_pop, n_regions, n_regions, n_yrs, n_seas, n_ages, n_sexes, sim_list$n_sims))
  for (r_from in 1:n_regions) for (r_to in 1:n_regions) {
    prob <- if (r_to == r_from) stay_prob else (1 - stay_prob) / (n_regions - 1)
    sim_list$Movement[, r_from, r_to, , , , , ] <- prob
  }

  sim_list <- Setup_Sim_Rec(
    sim_list = sim_list,
    R0_input = array(6, dim = c(n_pop, n_regions, n_yrs, sim_list$n_sims)),
    ln_sigmaR = array(log(0.5), dim = c(2, n_pop, n_regions)),
    recruitment_opt = "mean_rec",
    init_age_strc = 1
  )

  sim_obj <- Simulate_Pop_Static(sim_list = sim_list, output_path = NULL)

  sim_data <- simulation_data_to_SPoRC(sim_env = sim_obj, y = sim_obj$n_years, sim = 1)

  input_list <- Setup_Mod_Dim(
    years = 1:sim_obj$n_years, ages = 1:sim_obj$n_ages, lens = sim_obj$n_lens,
    n_regions = sim_obj$n_regions, n_sexes = sim_obj$n_sexes,
    n_fish_fleets = sim_obj$n_fish_fleets, n_srv_fleets = sim_obj$n_srv_fleets,
    n_pop = sim_obj$n_pop, natal_region = sim_obj$natal_region,
    verbose = FALSE, do_internal_comp_osa = TRUE
  )

  input_list <- Setup_Mod_Rec(
    input_list = input_list, do_rec_bias_ramp = 0, sigmaR_switch = 1,
    ln_sigmaR = array(log(0.5), c(2, input_list$data$n_pop, input_list$data$n_regions)),
    rec_model = "mean_rec", sigmaR_spec = "fix", init_age_strc = 1, equil_init_age_strc = 2,
    ln_global_R0 = log(6)
  )

  input_list <- Setup_Mod_Biologicals(
    input_list = input_list, WAA = sim_data$WAA, MatAA = sim_data$MatAA,
    WAA_fish = sim_data$WAA_fish, WAA_srv = sim_data$WAA_srv,
    fit_lengths = 0, AgeingError = sim_data$AgeingError, M_spec = "fix",
    Fixed_natmort = array(0.25, dim = c(input_list$data$n_pop, input_list$data$n_regions,
                                        length(input_list$data$years), length(input_list$data$ages),
                                        input_list$data$n_sexes))
  )

  input_list <- Setup_Mod_Tagging(input_list = input_list, use_conv_fish_tagging = 0)
  input_list <- Setup_Mod_Movement(
    input_list = input_list, use_fixed_movement = 1, Fixed_Movement = NA, do_recruits_move = 0
  )

  input_list <- Setup_Mod_Catch_and_F(
    input_list = input_list,
    ObsCatch = sim_data$ObsCatch, UseCatch = sim_data$UseCatch, Use_F_pen = 1,
    sigmaC_spec = "fix", ln_sigmaC = sim_data$ln_sigmaC,
    ln_sigmaF = array(log(1), dim = c(input_list$data$n_regions, input_list$data$n_seas, input_list$data$n_fish_fleets)),
    ObsDiscard = sim_data$ObsDiscard, UseDiscard = sim_data$UseDiscard,
    sigma_dmr_spec = "fix", dmr_mean_spec = "est_all", ln_sigmaD = sim_data$ln_sigmaD
  )

  input_list <- Setup_Mod_FishIdx_and_Comps(
    input_list = input_list,
    ObsFishIdx = sim_data$ObsFishIdx, ObsFishIdx_SE = sim_data$ObsFishIdx_SE, UseFishIdx = sim_data$UseFishIdx,
    ObsFishAgeComps = sim_data$ObsFishAgeComps, ObsFishLenComps = sim_data$ObsFishLenComps,
    UseFishAgeComps = sim_data$UseFishAgeComps, UseFishLenComps = sim_data$UseFishLenComps,
    ISS_FishAgeComps = sim_data$ISS_FishAgeComps, ISS_FishLenComps = sim_data$ISS_FishLenComps,
    fish_idx_type = c("biom", "biom"),
    FishAgeComps_LikeType = c("Multinomial", "Multinomial"), FishLenComps_LikeType = c("none", "none"),
    FishAgeComps_Type = c("spltRspltS_Year_1-terminal_Fleet_1", "spltRspltS_Year_1-terminal_Fleet_2"),
    FishLenComps_Type = c("none_Year_1-terminal_Fleet_1", "none_Year_1-terminal_Fleet_2"),
    ObsFishAgeComps_discard = sim_data$ObsFishAgeComps_discard, UseFishAgeComps_discard = sim_data$UseFishAgeComps_discard,
    ISS_FishAgeComps_discard = sim_data$ISS_FishAgeComps_discard,
    FishAgeComps_discard_LikeType = c("Multinomial", "Multinomial"),
    FishAgeComps_discard_Type = c("agg_Year_1-terminal_Fleet_1", "agg_Year_1-terminal_Fleet_2")
  )

  input_list <- Setup_Mod_SrvIdx_and_Comps(
    input_list = input_list,
    ObsSrvIdx = sim_data$ObsSrvIdx, ObsSrvIdx_SE = sim_data$ObsSrvIdx_SE, UseSrvIdx = sim_data$UseSrvIdx,
    ObsSrvAgeComps = sim_data$ObsSrvAgeComps, ObsSrvLenComps = sim_data$ObsSrvLenComps,
    UseSrvAgeComps = sim_data$UseSrvAgeComps, UseSrvLenComps = sim_data$UseSrvLenComps,
    ISS_SrvAgeComps = sim_data$ISS_SrvAgeComps, ISS_SrvLenComps = sim_data$ISS_SrvLenComps,
    srv_idx_type = c("biom", "biom"),
    SrvAgeComps_LikeType = c("Multinomial", "Multinomial"), SrvLenComps_LikeType = c("none", "none"),
    SrvAgeComps_Type = c("spltRjntS_Year_1-terminal_Fleet_1", "spltRjntS_Year_1-terminal_Fleet_2"),
    SrvLenComps_Type = c("none_Year_1-terminal_Fleet_1", "none_Year_1-terminal_Fleet_2")
  )

  input_list <- Setup_Mod_Fishsel_and_Q(
    input_list = input_list,
    fish_sel_model = c("logist1_Fleet_1", "logist1_Fleet_2"),
    fish_fixed_sel_pars_spec = c("est_all", "est_all"), fish_q_spec = c("est_all", "est_all"),
    ret_sel_model = c("logist1_Fleet_1", "logist1_Fleet_2"),
    ret_fixed_sel_pars_spec = c("est_all", "est_all"), use_fixed_ret_sel = c(0, 0)
  )

  input_list <- Setup_Mod_Srvsel_and_Q(
    input_list = input_list,
    srv_sel_model = c("logist1_Fleet_1", "logist1_Fleet_2"),
    srv_fixed_sel_pars_spec = c("est_all", "est_all"), srv_q_spec = c("est_all", "est_all")
  )

  input_list <- Setup_Mod_Weighting(
    input_list = input_list, Wt_Catch = 1, Wt_FishIdx = 1, Wt_SrvIdx = 1, Wt_Rec = 1, Wt_F = 1,
    Wt_Discard = 1, Wt_D = 1,
    Wt_FishAgeComps = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years),
                                       input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_fish_fleets)),
    Wt_FishLenComps = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years),
                                       input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_fish_fleets)),
    Wt_SrvAgeComps = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years),
                                      input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_srv_fleets)),
    Wt_FishAgeComps_discard = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years),
                                               input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_fish_fleets))
  )

  model <- fit_model(
    input_list$data, input_list$par, input_list$map,
    random = NULL, silent = TRUE, do_optim = TRUE, newton_loops = 5
  )

  # Retained fishery age comps: split by region & sex, 2 fleets -> SpltR_SpltS
  osa_ret <- get_osa(model = model, data = input_list$data, comp_source = "FishAge",
                     family = "discrete", bins = input_list$data$ages, bin_label = "Age")
  expect_setequal(unique(osa_ret$res$fleet), c("1", "2"))
  expect_setequal(unique(osa_ret$res$region), 1:2)
  expect_setequal(unique(osa_ret$res$sex), 1:2)
  expect_equal(unique(osa_ret$res$comp_type), "SpltR_SpltS")
  expect_true(is.finite(sd(osa_ret$res$resid)))
  p_ret <- plot_resids(osa_ret)
  expect_s3_class(p_ret[[1]], "ggplot")
  expect_s3_class(p_ret[[2]], "ggplot")

  # Discard fishery age comps: aggregated, 2 fleets -> Aggregated
  osa_disc <- get_osa(model = model, data = input_list$data, comp_source = "FishAge", discard = TRUE,
                      family = "discrete", bins = input_list$data$ages, bin_label = "Age")
  expect_setequal(unique(osa_disc$res$fleet), c("1", "2"))
  expect_equal(unique(osa_disc$res$comp_type), "Aggregated")
  expect_true(is.finite(sd(osa_disc$res$resid)))
  p_disc <- plot_resids(osa_disc)
  expect_s3_class(p_disc[[1]], "ggplot")
  expect_s3_class(p_disc[[2]], "ggplot")

  # Survey age comps: split by region, joint sex, 2 fleets -> SpltR_JntS
  osa_srv <- get_osa(model = model, data = input_list$data, comp_source = "SrvAge", parallel = TRUE,
                     family = "discrete", bins = input_list$data$ages, bin_label = "Age")
  expect_setequal(unique(osa_srv$res$fleet), c("1", "2"))
  expect_setequal(unique(osa_srv$res$region), 1:2)
  expect_equal(unique(osa_srv$res$comp_type), "SpltR_JntS")
  expect_true(is.finite(sd(osa_srv$res$resid)))
  p_srv <- plot_resids(osa_srv)
  expect_s3_class(p_srv[[1]], "ggplot")
  expect_s3_class(p_srv[[2]], "ggplot")
})
