library(SPoRC)
library(testthat)

test_that("get_osa()/plot_resids() support index-type sources (Catch, Discard, FishIdx, SrvIdx)", {

  set.seed(123)
  sim_list <- Setup_Sim_Dim(
    n_sims        = 1,
    n_yrs         = 50,
    n_regions     = 1,
    n_ages        = 10,
    n_lens        = NULL,
    n_sexes       = 1,
    n_fish_fleets = 1,
    n_srv_fleets  = 1,
    n_pop         = 1
  )

  sim_list <- Setup_Sim_Containers(sim_list)

  sim_list <- Setup_Sim_Fishing(
    sim_list = sim_list,
    fish_sel_input = replicate(
      n = sim_list$n_sims,
      array(
        rep(1 / (1 + exp(-3 * ((1:sim_list$n_ages) - 2))), each = sim_list$n_yrs),
        dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs,
                sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes,
                sim_list$n_fish_fleets)
      )
    ),
    ret_sel_input = replicate(
      n = sim_list$n_sims,
      array(
        rep(0.5 / (1 + exp(-3 * ((1:sim_list$n_ages) - 5))), each = sim_list$n_yrs),
        dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs,
                sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes,
                sim_list$n_fish_fleets)
      )
    ),
    dmr_input = array(0.5, dim = c(sim_list$n_regions, sim_list$n_yrs,
                                   sim_list$n_seas, sim_list$n_fish_fleets,
                                   sim_list$n_sims))
  )

  sim_list <- Setup_Sim_Survey(
    sim_list = sim_list,
    srv_sel_input = replicate(
      n = sim_list$n_sims,
      array(
        rep(1 / (1 + exp(-1 * ((1:sim_list$n_ages) - 3))), each = sim_list$n_yrs),
        dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs,
                sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes,
                sim_list$n_srv)
      )
    )
  )

  sim_list <- suppressWarnings(
    Setup_Sim_Biologicals(
      sim_list = sim_list,
      natmort_input = replicate(
        n = sim_list$n_sims,
        array(0.3, dim = c(sim_list$n_pop, sim_list$n_regions,
                           sim_list$n_yrs, sim_list$n_ages, sim_list$n_sexes))
      ),
      WAA_input = replicate(
        n = sim_list$n_sims,
        array(
          rep(5 / (1 + exp(-3 * ((1:sim_list$n_ages) - 3))), each = sim_list$n_yrs),
          dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs,
                  sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes)
        )
      ),
      WAA_fish_input = replicate(
        n = sim_list$n_sims,
        array(
          rep(5 / (1 + exp(-3 * ((1:sim_list$n_ages) - 3))), each = sim_list$n_yrs),
          dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs,
                  sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes,
                  sim_list$n_fish_fleets)
        )
      ),
      WAA_srv_input = replicate(
        n = sim_list$n_sims,
        array(
          rep(5 / (1 + exp(-3 * ((1:sim_list$n_ages) - 3))), each = sim_list$n_yrs),
          dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs,
                  sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes,
                  sim_list$n_srv_fleets)
        )
      ),
      MatAA_input = replicate(
        n = sim_list$n_sims,
        array(
          rep(1 / (1 + exp(-3 * ((1:sim_list$n_ages) - 3))), each = sim_list$n_yrs),
          dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs,
                  sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes)
        )
      )
    )
  )

  sim_list <- Setup_Sim_Tagging(sim_list = sim_list, use_conv_fish_tagging = 0)

  sim_list$Movement <- array(
    1,
    dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_regions,
            sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages,
            sim_list$n_sexes, sim_list$n_sims)
  )

  sim_list <- Setup_Sim_Rec(
    sim_list       = sim_list,
    R0_input       = replicate(
      n    = sim_list$n_sims,
      expr = array(5, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs))
    ),
    ln_sigmaR      = array(log(0.5), dim = c(2, sim_list$n_pop, sim_list$n_region)),
    recruitment_opt = "mean_rec",
    init_age_strc  = 1
  )

  set.seed(777)
  sim_obj <- Simulate_Pop_Static(sim_list = sim_list, output_path = NULL)

  sim_data <- simulation_data_to_SPoRC(sim_env = sim_obj, y = sim_obj$n_years, sim = 1)

  input_list <- Setup_Mod_Dim(
    years          = 1:sim_obj$n_years,
    ages           = 1:sim_obj$n_ages,
    lens           = sim_obj$n_lens,
    n_regions      = sim_obj$n_regions,
    n_sexes        = sim_obj$n_sexes,
    n_fish_fleets  = sim_obj$n_fish_fleets,
    n_srv_fleets   = sim_obj$n_srv_fleets,
    n_pop          = sim_obj$n_pop,
    natal_region   = sim_obj$natal_region,
    verbose        = FALSE,
    do_internal_comp_osa = TRUE
  )

  input_list <- Setup_Mod_Rec(
    input_list        = input_list,
    do_rec_bias_ramp  = 0,
    sigmaR_switch     = 1,
    ln_sigmaR         = array(log(0.5), c(2, input_list$data$n_pop, input_list$data$n_regions)),
    rec_model         = "mean_rec",
    sigmaR_spec       = "fix",
    init_age_strc     = 1,
    equil_init_age_strc = 2,
    ln_global_R0      = log(5)
  )

  input_list <- Setup_Mod_Biologicals(
    input_list    = input_list,
    WAA           = sim_data$WAA,
    MatAA         = sim_data$MatAA,
    WAA_fish      = sim_data$WAA_fish,
    WAA_srv       = sim_data$WAA_srv,
    fit_lengths   = 0,
    AgeingError   = sim_data$AgeingError,
    M_spec        = "fix",
    Fixed_natmort = array(0.3, dim = c(input_list$data$n_pop,
                                       input_list$data$n_regions,
                                       length(input_list$data$years),
                                       length(input_list$data$ages),
                                       input_list$data$n_sexes))
  )

  input_list <- Setup_Mod_Tagging(input_list = input_list, use_conv_fish_tagging = 0)
  input_list <- Setup_Mod_Movement(
    input_list         = input_list,
    use_fixed_movement = 1,
    Fixed_Movement     = NA,
    do_recruits_move   = 0
  )

  input_list <- Setup_Mod_Catch_and_F(
    input_list     = input_list,
    ObsCatch       = sim_data$ObsCatch,
    UseCatch       = sim_data$UseCatch,
    Use_F_pen      = 1,
    sigmaC_spec    = "fix",
    ln_sigmaC      = sim_data$ln_sigmaC,
    ln_sigmaF      = array(log(1), dim = c(input_list$data$n_regions,
                                           input_list$data$n_seas,
                                           input_list$data$n_fish_fleets)),
    ObsDiscard     = sim_data$ObsDiscard,
    UseDiscard     = sim_data$UseDiscard,
    sigma_dmr_spec = "fix",
    dmr_mean_spec  = "est_all",
    ln_sigmaD      = sim_data$ln_sigmaD
  )

  input_list <- Setup_Mod_FishIdx_and_Comps(
    input_list = input_list,
    ObsFishIdx    = sim_data$ObsFishIdx,
    ObsFishIdx_SE = sim_data$ObsFishIdx_SE,
    UseFishIdx    = sim_data$UseFishIdx,
    ObsFishAgeComps  = sim_data$ObsFishAgeComps,
    ObsFishLenComps  = sim_data$ObsFishLenComps,
    UseFishAgeComps  = sim_data$UseFishAgeComps,
    UseFishLenComps  = sim_data$UseFishLenComps,
    ISS_FishAgeComps = sim_data$ISS_FishAgeComps,
    ISS_FishLenComps = sim_data$ISS_FishLenComps,
    fish_idx_type        = "biom",
    FishAgeComps_LikeType = "Multinomial",
    FishLenComps_LikeType = "none",
    FishAgeComps_Type     = "agg_Year_1-terminal_Fleet_1",
    FishLenComps_Type     = "none_Year_1-terminal_Fleet_1",
    ObsFishAgeComps_discard      = sim_data$ObsFishAgeComps_discard,
    UseFishAgeComps_discard      = sim_data$UseFishAgeComps_discard,
    ISS_FishAgeComps_discard     = sim_data$ISS_FishAgeComps_discard,
    FishAgeComps_discard_LikeType = rep("Multinomial", input_list$data$n_fish_fleets),
    FishAgeComps_discard_Type     = "agg_Year_1-terminal_Fleet_1"
  )

  input_list <- Setup_Mod_SrvIdx_and_Comps(
    input_list = input_list,
    ObsSrvIdx        = sim_data$ObsSrvIdx,
    ObsSrvIdx_SE     = sim_data$ObsSrvIdx_SE,
    UseSrvIdx        = sim_data$UseSrvIdx,
    ObsSrvAgeComps   = sim_data$ObsSrvAgeComps,
    ObsSrvLenComps   = sim_data$ObsSrvLenComps,
    UseSrvAgeComps   = sim_data$UseSrvAgeComps,
    UseSrvLenComps   = sim_data$UseSrvLenComps,
    ISS_SrvAgeComps  = sim_data$ISS_SrvAgeComps,
    ISS_SrvLenComps  = sim_data$ISS_SrvLenComps,
    srv_idx_type          = "biom",
    SrvAgeComps_LikeType  = "Multinomial",
    SrvLenComps_LikeType  = "none",
    SrvAgeComps_Type      = "agg_Year_1-terminal_Fleet_1",
    SrvLenComps_Type      = "none_Year_1-terminal_Fleet_1"
  )

  input_list <- Setup_Mod_Fishsel_and_Q(
    input_list              = input_list,
    fish_sel_model          = "logist1_Fleet_1",
    fish_fixed_sel_pars_spec = "est_all",
    fish_q_spec             = "est_all",
    ret_sel_model           = "asymplogist1_Fleet_1",
    ret_fixed_sel_pars_spec = "est_all",
    use_fixed_ret_sel       = 0
  )

  input_list <- Setup_Mod_Srvsel_and_Q(
    input_list              = input_list,
    srv_sel_model           = "logist1_Fleet_1",
    srv_fixed_sel_pars_spec = "est_all",
    srv_q_spec              = "est_all"
  )

  input_list <- Setup_Mod_Weighting(
    input_list     = input_list,
    Wt_Catch       = 1,
    Wt_FishIdx     = 1,
    Wt_SrvIdx      = 1,
    Wt_Rec         = 1,
    Wt_F           = 1,
    Wt_Discard     = 1,
    Wt_D           = 1,
    Wt_FishAgeComps = array(1, dim = c(input_list$data$n_regions,
                                       length(input_list$data$years),
                                       input_list$data$n_seas,
                                       input_list$data$n_sexes,
                                       input_list$data$n_fish_fleets)),
    Wt_FishLenComps = array(1, dim = c(input_list$data$n_regions,
                                       length(input_list$data$years),
                                       input_list$data$n_seas,
                                       input_list$data$n_sexes,
                                       input_list$data$n_fish_fleets)),
    Wt_SrvAgeComps = array(1, dim = c(input_list$data$n_regions,
                                      length(input_list$data$years),
                                      input_list$data$n_seas,
                                      input_list$data$n_sexes,
                                      input_list$data$n_srv_fleets)),
    Wt_FishAgeComps_discard = array(1, dim = c(input_list$data$n_regions,
                                               length(input_list$data$years),
                                               input_list$data$n_seas,
                                               input_list$data$n_sexes,
                                               input_list$data$n_fish_fleets))
  )

  model <- fit_model(
    input_list$data,
    input_list$par,
    input_list$map,
    random       = NULL,
    silent       = TRUE,
    do_optim     = TRUE,
    newton_loops = 3
  )

  expect_jnLL_decomposes(model)

  # test internals
  expect_error(index_osa_field_map("naw"), "must be one of")
  fm <- index_osa_field_map("SrvIdx", pop = TRUE)
  expect_equal(fm$Obs, "ObsSrvIdx_pop")
  expect_equal(fm$Use, "UseSrvIdx_pop")

  # Test if actuall runs
  for(src in c("Catch", "Discard", "FishIdx", "SrvIdx")) {
    out <- get_osa(model = model, data = input_list$data, index_source = src)
    expect_false(is.null(out))
    res <- out$res
    expect_true(all(c("fleet","region","year","season","pop","resid","idx_type") %in% names(res)))
    expect_false("comp_type" %in% names(res))  # index residuals use idx_type, not comp_type
    expect_equal(unique(res$idx_type), src)
    expect_true(all(is.finite(res$resid)))
    direct <- RTMB::oneStepPredict(
      model,
      observation.name = paste0("Obs", src),
      method = "oneStepGeneric",
      discrete = FALSE,
      parallel = FALSE,
      trace = FALSE
    )
    expect_equal(res$resid, direct$residual, tolerance = 1e-3)
    expect_equal(nrow(res), sum(input_list$data[[paste0("Use", src)]] == 1))
    plots <- plot_resids(out)
    expect_length(plots, 2)
    expect_s3_class(plots[[1]], "ggplot")
    expect_s3_class(plots[[2]], "ggplot")
    ggplot2::ggplot_build(plots[[1]])
    ggplot2::ggplot_build(plots[[2]])
  }

  # --- pop-specific sources absent in this single-pop setup -> NULL + warning
  expect_warning(
    out_pop <- get_osa(model = model, data = input_list$data, index_source = "SrvIdx", pop = TRUE),
    "No data found"
  )
  expect_null(out_pop)

  # --- osa_method = "cdf" is disallowed for every internal-OSA branch --------
  expect_error(
    get_osa(model = model, data = input_list$data, index_source = "SrvIdx", osa_method = "cdf"),
    "not permitted"
  )
  expect_error(
    get_osa(
      model = model,
      data = input_list$data,
      comp_source = "FishAge",
      family = "discrete",
      bins = input_list$data$ages,
      bin_label = "Age",
      osa_method = "cdf"
    ),
    "not permitted"
  )
  # a valid oneStep* method still works
  out_valid <- get_osa(
    model = model,
    data = input_list$data,
    index_source = "SrvIdx",
    osa_method = "oneStepGaussian"
  )
  expect_false(is.null(out_valid))
  expect_true(all(is.finite(out_valid$res$resid)))

})
