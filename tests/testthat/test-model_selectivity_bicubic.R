library(SPoRC)
library(testthat)

test_that("bicubic fishery selectivity (Selex_Model == 8) wires correctly through Setup_Mod_Fishsel_and_Q and SPoRC_rtmb", {

  # Build a minimal single-region/single-sex operating model + simulated data, purely as scaffolding
  # to get a valid input_list$data through the full setup pipeline. The *true* simulated selectivity
  # is an ordinary logistic curve -- what's under test is the *estimation* model's bicubic wiring
  # (Setup_Mod_Fishsel_and_Q's array construction and SPoRC_rtmb's Get_Selex call), not parameter
  # recovery, so no optimization is performed here.
  sim_list <- Setup_Sim_Dim(
    n_sims = 1,
    n_yrs = 12,
    n_regions = 1,
    n_ages = 8,
    n_lens = NULL,
    n_sexes = 1,
    n_fish_fleets = 1,
    n_srv_fleets = 1,
    n_pop = 1
  )
  sim_list <- Setup_Sim_Containers(sim_list)

  sim_list <- Setup_Sim_Fishing(
    sim_list = sim_list,
    fish_sel_input = replicate(
      n = sim_list$n_sims,
      array(rep(1 / (1 + exp(-3 * ((1:sim_list$n_ages) - 4))), each = sim_list$n_yrs),
            dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages,
                    sim_list$n_sexes, sim_list$n_fish_fleets))
    )
  )

  sim_list <- Setup_Sim_Survey(
    sim_list = sim_list,
    srv_sel_input = replicate(
      n = sim_list$n_sims,
      array(rep(1 / (1 + exp(-1 * ((1:sim_list$n_ages) - 2))), each = sim_list$n_yrs),
            dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages,
                    sim_list$n_sexes, sim_list$n_srv))
    )
  )

  sim_list <- suppressWarnings(
    Setup_Sim_Biologicals(
      sim_list = sim_list,
      natmort_input = replicate(n = sim_list$n_sims, array(0.3, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs,
                                                                        sim_list$n_ages, sim_list$n_sexes))),
      WAA_input = replicate(n = sim_list$n_sims, array(rep(5 / (1 + exp(-3 * ((1:sim_list$n_ages) - 3))), each = sim_list$n_yrs),
                                                       dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes))),
      WAA_fish_input = replicate(n = sim_list$n_sims, array(rep(5 / (1 + exp(-3 * ((1:sim_list$n_ages) - 3))), each = sim_list$n_yrs),
                                                            dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes, sim_list$n_fish_fleets))),
      WAA_srv_input = replicate(n = sim_list$n_sims, array(rep(5 / (1 + exp(-3 * ((1:sim_list$n_ages) - 3))), each = sim_list$n_yrs),
                                                           dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes, sim_list$n_srv_fleets))),
      MatAA_input = replicate(n = sim_list$n_sims, array(rep(1 / (1 + exp(-3 * ((1:sim_list$n_ages) - 3))), each = sim_list$n_yrs),
                                                         dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes)))
    )
  )

  sim_list <- Setup_Sim_Tagging(sim_list = sim_list, use_conv_fish_tagging = 0)
  sim_list$Movement <- array(1, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes, sim_list$n_sims))

  sim_list <- Setup_Sim_Rec(
    sim_list = sim_list,
    R0_input = replicate(n = sim_list$n_sims, expr = array(5, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs))),
    rinit_input = array(2, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_sims)),
    use_rinit = 1,
    ln_sigmaR = array(log(1), dim = c(2, sim_list$n_pop, sim_list$n_region)),
    recruitment_opt = 'mean_rec',
    init_age_strc = 1
  )

  set.seed(123)
  sim_obj <- Simulate_Pop_Static(sim_list = sim_list, output_path = NULL)
  sim_data <- simulation_data_to_SPoRC(sim_env = sim_obj, y = sim_obj$n_years, sim = 1)

  build_input_list <- function(
    fish_sel_model,
    fish_fixed_sel_pars_spec = "fix",
    srv_sel_model = c("logist2_Fleet_1"),
    srv_fixed_sel_pars_spec = c("est_all"),
    srv_extra_args = list(),
    fish_sel_pen_wts = NULL,
    srv_sel_pen_wts = NULL,
    ret_sel_pen_wts = NULL,
    ...
  ) {

    input_list <- Setup_Mod_Dim(
      years = 1:sim_obj$n_years,
      ages = 1:sim_obj$n_ages,
      lens = sim_obj$n_lens,
      n_regions = sim_obj$n_regions,
      n_sexes = sim_obj$n_sexes,
      n_fish_fleets = sim_obj$n_fish_fleets,
      n_srv_fleets = sim_obj$n_srv_fleets,
      n_pop = sim_obj$n_pop,
      natal_region = sim_obj$natal_region,
      verbose = FALSE
    )

    input_list <- Setup_Mod_Rec(
      input_list = input_list,
      do_rec_bias_ramp = 0,
      sigmaR_switch = 1,
      ln_sigmaR = array(log(1), c(2, input_list$data$n_pop, input_list$data$n_regions)),
      rec_model = "mean_rec",
      use_rinit = 1,
      sigmaR_spec = "fix",
      init_age_strc = 1,
      equil_init_age_strc = 2,
      ln_global_R0 = log(5),
      ln_rinit = log(2)
    )

    input_list <- Setup_Mod_Biologicals(
      input_list = input_list,
      WAA = sim_data$WAA,
      MatAA = sim_data$MatAA,
      WAA_fish = sim_data$WAA_fish,
      WAA_srv = sim_data$WAA_srv,
      fit_lengths = 0,
      AgeingError = sim_data$AgeingError,
      M_spec = "fix",
      Fixed_natmort = array(0.3, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years),
                                         length(input_list$data$ages), input_list$data$n_sexes))
    )

    input_list <- Setup_Mod_Tagging(input_list = input_list, use_conv_fish_tagging = 0)
    input_list <- Setup_Mod_Movement(
      input_list = input_list,
      use_fixed_movement = 1,
      Fixed_Movement = NA,
      do_recruits_move = 0
    )

    input_list <- suppressWarnings(
      Setup_Mod_Catch_and_F(
        input_list = input_list,
        ObsCatch = sim_data$ObsCatch,
        UseCatch = sim_data$UseCatch,
        Use_F_pen = 1,
        sigmaC_spec = "fix",
        ln_sigmaC = sim_data$ln_sigmaC,
        ln_sigmaF = array(log(1), dim = c(input_list$data$n_regions, input_list$data$n_seas, input_list$data$n_fish_fleets))
      )
    )

    input_list <- Setup_Mod_FishIdx_and_Comps(
      input_list = input_list,
      ObsFishIdx = sim_data$ObsFishIdx,
      ObsFishIdx_SE = sim_data$ObsFishIdx_SE,
      UseFishIdx = sim_data$UseFishIdx,
      ObsFishAgeComps = sim_data$ObsFishAgeComps,
      ObsFishLenComps = sim_data$ObsFishLenComps,
      UseFishAgeComps = sim_data$UseFishAgeComps,
      UseFishLenComps = sim_data$UseFishLenComps,
      ISS_FishAgeComps = sim_data$ISS_FishAgeComps,
      ISS_FishLenComps = sim_data$ISS_FishLenComps,
      fish_idx_type = c("biom"),
      FishAgeComps_LikeType = c("Multinomial"),
      FishLenComps_LikeType = c("none"),
      FishAgeComps_Type = c("agg_Year_1-terminal_Fleet_1"),
      FishLenComps_Type = c("none_Year_1-terminal_Fleet_1")
    )

    input_list <- Setup_Mod_SrvIdx_and_Comps(
      input_list = input_list,
      ObsSrvIdx = sim_data$ObsSrvIdx,
      ObsSrvIdx_SE = sim_data$ObsSrvIdx_SE,
      UseSrvIdx = sim_data$UseSrvIdx,
      ObsSrvAgeComps = sim_data$ObsSrvAgeComps,
      ObsSrvLenComps = sim_data$ObsSrvLenComps,
      UseSrvAgeComps = sim_data$UseSrvAgeComps,
      UseSrvLenComps = sim_data$UseSrvLenComps,
      ISS_SrvAgeComps = sim_data$ISS_SrvAgeComps,
      ISS_SrvLenComps = sim_data$ISS_SrvLenComps,
      srv_idx_type = c("biom"),
      SrvAgeComps_LikeType = c("Multinomial"),
      SrvLenComps_LikeType = c("none"),
      SrvAgeComps_Type = c("agg_Year_1-terminal_Fleet_1"),
      SrvLenComps_Type = c("none_Year_1-terminal_Fleet_1")
    )

    input_list <- Setup_Mod_Fishsel_and_Q(
      input_list = input_list,
      fish_sel_model = fish_sel_model,
      fish_fixed_sel_pars_spec = fish_fixed_sel_pars_spec,
      fish_q_spec = "est_all",
      ...
    )

    input_list <- do.call(Setup_Mod_Srvsel_and_Q, c(
      list(
        input_list = input_list,
        srv_sel_model = srv_sel_model,
        srv_fixed_sel_pars_spec = srv_fixed_sel_pars_spec,
        srv_q_spec = c("est_all")
      ),
      srv_extra_args
    ))

    input_list <- Setup_Mod_Weighting(
      input_list = input_list,
      Wt_Catch = 1,
      Wt_FishIdx = 1,
      Wt_SrvIdx = 1,
      Wt_Rec = 1,
      Wt_F = 1,
      Wt_Tagging = 0,
      Wt_FishAgeComps = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_fish_fleets)),
      Wt_FishLenComps = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_fish_fleets)),
      Wt_SrvAgeComps = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_srv_fleets)),
      Wt_SrvLenComps = array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_srv_fleets)),
      fish_sel_pen_wts = fish_sel_pen_wts,
      srv_sel_pen_wts = srv_sel_pen_wts,
      ret_sel_pen_wts = ret_sel_pen_wts
    )

    return(input_list)
  }

  test_that("model builds, evaluates, and differentiates without error with a single bicubic block", {
    n_bin_nodes <- 4; n_yr_nodes <- 3
    node_par <- matrix(
      c(-1.5, -0.2, 0.6, 1.2,
                         -0.8,  0.3, 0.9, 0.5,
                         -0.3,  0.1, 0.4, -0.2),
      nrow = n_yr_nodes,
      ncol = n_bin_nodes,
      byrow = TRUE
    )
    fixed_pars_init <- array(0, dim = c(1, n_bin_nodes * n_yr_nodes, 1, 1, 1))
    fixed_pars_init[1, , 1, 1, 1] <- as.vector(node_par)

    input_list <- build_input_list(
      fish_sel_model = paste0("bicubic_Bin_", n_bin_nodes, "_Yr_", n_yr_nodes, "_Fleet_1"),
      fish_fixed_sel_pars_spec = "est_all",
      fish_fixed_sel_pars = fixed_pars_init
    )

    obj <- RTMB::MakeADFun(
      cmb(SPoRC_rtmb, input_list$data),
      parameters = input_list$par,
      map = input_list$map,
      silent = TRUE
    )

    expect_no_error(obj$fn(obj$par))
    expect_no_error(obj$gr(obj$par))

    rep <- obj$report(obj$par)

    # manual computation using the exact weight matrices Setup_Mod_Fishsel_and_Q constructed
    Wbin <- input_list$data$fish_sel_bicubic_Wbin[1, , 1:n_bin_nodes, 1, 1]
    Wyr  <- input_list$data$fish_sel_bicubic_Wyr[1, , 1:n_yr_nodes, 1, 1]
    expected_surface <- exp(Wyr %*% (node_par %*% t(Wbin))) # n_yrs_total x n_ages

    n_yrs <- length(input_list$data$years)
    for (y in 1:n_yrs) {
      expect_equal(as.vector(rep$fish_sel[1, 1, y, 1, , 1, 1]), expected_surface[y, ],
                  tolerance = 1e-8, label = sprintf("fish_sel year %d", y))
    }
  })

  test_that("n_yr_nodes == 1 gives a time-invariant fishery selectivity surface end-to-end", {
    n_bin_nodes <- 5
    node_par <- c(-1, 0.2, 1.1, 0.4, -0.6)
    fixed_pars_init <- array(0, dim = c(1, n_bin_nodes, 1, 1, 1))
    fixed_pars_init[1, , 1, 1, 1] <- node_par

    input_list <- build_input_list(
      fish_sel_model = paste0("bicubic_Bin_", n_bin_nodes, "_Yr_1_Fleet_1"),
      fish_fixed_sel_pars_spec = "est_all",
      fish_fixed_sel_pars = fixed_pars_init
    )

    obj <- RTMB::MakeADFun(
      cmb(SPoRC_rtmb, input_list$data),
      parameters = input_list$par,
      map = input_list$map,
      silent = TRUE
    )
    rep <- obj$report(obj$par)

    n_yrs <- length(input_list$data$years)
    sel_y1 <- as.vector(rep$fish_sel[1, 1, 1, 1, , 1, 1])
    for (y in 2:n_yrs) {
      expect_equal(as.vector(rep$fish_sel[1, 1, y, 1, , 1, 1]), sel_y1, tolerance = 1e-8)
    }

    Wbin <- input_list$data$fish_sel_bicubic_Wbin[1, , 1:n_bin_nodes, 1, 1]
    expect_equal(sel_y1, exp(as.vector(Wbin %*% node_par)), tolerance = 1e-8)
  })

  test_that("multiple bicubic year-blocks give a piecewise-constant-in-year, smooth-in-age surface end-to-end", {
    n_bin_nodes <- 4
    half <- floor(sim_obj$n_years / 2)
    block1_yrs <- paste0(1, "-", half)
    block2_yrs <- paste0(half + 1, "-terminal")

    node_par_1 <- c(-1.2, 0.4, 0.8, -0.1)
    node_par_2 <- c(0.5, -0.7, 1.0, 0.9)
    fixed_pars_init <- array(0, dim = c(1, n_bin_nodes, 2, 1, 1))
    fixed_pars_init[1, , 1, 1, 1] <- node_par_1
    fixed_pars_init[1, , 2, 1, 1] <- node_par_2

    input_list <- build_input_list(
      fish_sel_model = c(
        paste0("bicubic_Bin_", n_bin_nodes, "_Yr_1_Block_1_Fleet_1"),
        paste0("bicubic_Bin_", n_bin_nodes, "_Yr_1_Block_2_Fleet_1")
      ),
      fish_sel_blocks = c(
        paste0("Block_1_Year_", block1_yrs, "_Fleet_1"),
        paste0("Block_2_Year_", block2_yrs, "_Fleet_1")
      ),
      fish_fixed_sel_pars_spec = "est_all",
      fish_fixed_sel_pars = fixed_pars_init
    )

    obj <- RTMB::MakeADFun(
      cmb(SPoRC_rtmb, input_list$data),
      parameters = input_list$par,
      map = input_list$map,
      silent = TRUE
    )
    rep <- obj$report(obj$par)

    Wbin <- input_list$data$fish_sel_bicubic_Wbin[1, , 1:n_bin_nodes, 1, 1] # same age nodes/bins for both blocks
    expected_1 <- exp(as.vector(Wbin %*% node_par_1))
    expected_2 <- exp(as.vector(Wbin %*% node_par_2))

    expect_equal(as.vector(rep$fish_sel[1, 1, 1, 1, , 1, 1]), expected_1, tolerance = 1e-8)
    expect_equal(as.vector(rep$fish_sel[1, 1, half, 1, , 1, 1]), expected_1, tolerance = 1e-8)
    expect_equal(as.vector(rep$fish_sel[1, 1, half + 1, 1, , 1, 1]), expected_2, tolerance = 1e-8)
    expect_equal(as.vector(rep$fish_sel[1, 1, sim_obj$n_years, 1, , 1, 1]), expected_2, tolerance = 1e-8)
  })

  # ── SelStyr: sub-range bicubic fit start, matching ADMB's fsh_sel_styr ───────

  test_that("SelStyr holds pre-fit years constant at the SelStyr year's curve and fits the spline only from SelStyr onward", {
    n_bin_nodes <- 4; n_yr_nodes <- 3
    selstyr_year <- 5 # years 1:4 should be edge-kept; the spline is fit only over years 5:n_years
    set.seed(202)
    node_par <- matrix(rnorm(n_bin_nodes * n_yr_nodes, sd = 0.5), nrow = n_yr_nodes, ncol = n_bin_nodes)
    fixed_pars_init <- array(0, dim = c(1, n_bin_nodes * n_yr_nodes, 1, 1, 1))
    fixed_pars_init[1, , 1, 1, 1] <- as.vector(node_par)

    input_list <- build_input_list(
      fish_sel_model = paste0("bicubic_Bin_", n_bin_nodes, "_Yr_", n_yr_nodes, "_SelStyr_", selstyr_year, "_Fleet_1"),
      fish_fixed_sel_pars_spec = "est_all",
      fish_fixed_sel_pars = fixed_pars_init
    )

    expect_equal(input_list$data$fish_sel_bicubic_selstyr[1, 1, 1], selstyr_year)

    obj <- RTMB::MakeADFun(
      cmb(SPoRC_rtmb, input_list$data),
      parameters = input_list$par,
      map = input_list$map,
      silent = TRUE
    )
    rep <- obj$report(obj$par)

    n_yrs <- length(input_list$data$years)
    fit_years <- selstyr_year:n_yrs

    Wbin <- input_list$data$fish_sel_bicubic_Wbin[1, , 1:n_bin_nodes, 1, 1]
    Wyr_fit <- Get_Natural_Cubic_Spline_Weights(seq(0, 1, length.out = n_yr_nodes), seq(0, 1, length.out = length(fit_years)))
    expected_fit_surface <- exp(Wyr_fit %*% (node_par %*% t(Wbin)))

    # years before SelStyr are all identical to the SelStyr year's fitted curve
    sel_at_selstyr <- as.vector(rep$fish_sel[1, 1, selstyr_year, 1, , 1, 1])
    for (y in 1:(selstyr_year - 1)) {
      expect_equal(as.vector(rep$fish_sel[1, 1, y, 1, , 1, 1]), sel_at_selstyr, tolerance = 1e-8, label = sprintf("pre-SelStyr year %d", y))
    }
    expect_equal(sel_at_selstyr, expected_fit_surface[1, ], tolerance = 1e-8)

    # years from SelStyr onward follow the spline fit over the sub-range only
    for (y in fit_years) {
      expect_equal(as.vector(rep$fish_sel[1, 1, y, 1, , 1, 1]), expected_fit_surface[y - selstyr_year + 1, ],
                  tolerance = 1e-8, label = sprintf("fit year %d", y))
    }
  })

  # ── Survey selectivity: same bicubic wiring, mirrored for srv_sel_model ──────

  test_that("survey bicubic model builds, evaluates, and differentiates without error with a single block", {
    n_bin_nodes <- 4; n_yr_nodes <- 3
    node_par <- matrix(
      c(-1.0, 0.3, 0.7, -0.2,
                         0.4, -0.5, 0.2, 0.9,
                         -0.6, 0.8, -0.1, 0.3),
      nrow = n_yr_nodes,
      ncol = n_bin_nodes,
      byrow = TRUE
    )
    fixed_pars_init <- array(0, dim = c(1, n_bin_nodes * n_yr_nodes, 1, 1, 1))
    fixed_pars_init[1, , 1, 1, 1] <- as.vector(node_par)

    input_list <- build_input_list(
      fish_sel_model = "logist2_Fleet_1",
      srv_sel_model = paste0("bicubic_Bin_", n_bin_nodes, "_Yr_", n_yr_nodes, "_Fleet_1"),
      srv_fixed_sel_pars_spec = "est_all",
      srv_extra_args = list(srv_fixed_sel_pars = fixed_pars_init)
    )

    obj <- RTMB::MakeADFun(
      cmb(SPoRC_rtmb, input_list$data),
      parameters = input_list$par,
      map = input_list$map,
      silent = TRUE
    )

    expect_no_error(obj$fn(obj$par))
    expect_no_error(obj$gr(obj$par))

    rep <- obj$report(obj$par)

    Wbin <- input_list$data$srv_sel_bicubic_Wbin[1, , 1:n_bin_nodes, 1, 1]
    Wyr  <- input_list$data$srv_sel_bicubic_Wyr[1, , 1:n_yr_nodes, 1, 1]
    expected_surface <- exp(Wyr %*% (node_par %*% t(Wbin)))

    n_yrs <- length(input_list$data$years)
    for (y in 1:n_yrs) {
      expect_equal(as.vector(rep$srv_sel[1, 1, y, 1, , 1, 1]), expected_surface[y, ],
                  tolerance = 1e-8, label = sprintf("srv_sel year %d", y))
    }
  })

  test_that("survey bicubic n_yr_nodes == 1 gives a time-invariant selectivity surface end-to-end", {
    n_bin_nodes <- 5
    node_par <- c(0.3, -0.4, 1.0, 0.1, -0.9)
    fixed_pars_init <- array(0, dim = c(1, n_bin_nodes, 1, 1, 1))
    fixed_pars_init[1, , 1, 1, 1] <- node_par

    input_list <- build_input_list(
      fish_sel_model = "logist2_Fleet_1",
      srv_sel_model = paste0("bicubic_Bin_", n_bin_nodes, "_Yr_1_Fleet_1"),
      srv_fixed_sel_pars_spec = "est_all",
      srv_extra_args = list(srv_fixed_sel_pars = fixed_pars_init)
    )

    obj <- RTMB::MakeADFun(
      cmb(SPoRC_rtmb, input_list$data),
      parameters = input_list$par,
      map = input_list$map,
      silent = TRUE
    )
    rep <- obj$report(obj$par)

    n_yrs <- length(input_list$data$years)
    sel_y1 <- as.vector(rep$srv_sel[1, 1, 1, 1, , 1, 1])
    for (y in 2:n_yrs) {
      expect_equal(as.vector(rep$srv_sel[1, 1, y, 1, , 1, 1]), sel_y1, tolerance = 1e-8)
    }

    Wbin <- input_list$data$srv_sel_bicubic_Wbin[1, , 1:n_bin_nodes, 1, 1]
    expect_equal(sel_y1, exp(as.vector(Wbin %*% node_par)), tolerance = 1e-8)
  })

  test_that("survey and fishery fleets can each independently use bicubic selectivity at once", {
    n_bin_nodes <- 4
    fish_node_par <- c(-0.5, 0.5, 0.9, -0.3)
    srv_node_par  <- c(0.6, -0.8, 0.2, 1.0)

    fish_fixed_pars_init <- array(0, dim = c(1, n_bin_nodes, 1, 1, 1))
    fish_fixed_pars_init[1, , 1, 1, 1] <- fish_node_par
    srv_fixed_pars_init <- array(0, dim = c(1, n_bin_nodes, 1, 1, 1))
    srv_fixed_pars_init[1, , 1, 1, 1] <- srv_node_par

    input_list <- build_input_list(
      fish_sel_model = paste0("bicubic_Bin_", n_bin_nodes, "_Yr_1_Fleet_1"),
      fish_fixed_sel_pars_spec = "est_all",
      srv_sel_model = paste0("bicubic_Bin_", n_bin_nodes, "_Yr_1_Fleet_1"),
      srv_fixed_sel_pars_spec = "est_all",
      srv_extra_args = list(srv_fixed_sel_pars = srv_fixed_pars_init),
      fish_fixed_sel_pars = fish_fixed_pars_init
    )

    obj <- RTMB::MakeADFun(
      cmb(SPoRC_rtmb, input_list$data),
      parameters = input_list$par,
      map = input_list$map,
      silent = TRUE
    )
    expect_no_error(obj$fn(obj$par))
    rep <- obj$report(obj$par)

    Wbin_fish <- input_list$data$fish_sel_bicubic_Wbin[1, , 1:n_bin_nodes, 1, 1]
    Wbin_srv  <- input_list$data$srv_sel_bicubic_Wbin[1, , 1:n_bin_nodes, 1, 1]

    expect_equal(as.vector(rep$fish_sel[1, 1, 1, 1, , 1, 1]), exp(as.vector(Wbin_fish %*% fish_node_par)), tolerance = 1e-8)
    expect_equal(as.vector(rep$srv_sel[1, 1, 1, 1, , 1, 1]), exp(as.vector(Wbin_srv %*% srv_node_par)), tolerance = 1e-8)
  })

  # ── ADMB-aligned smoothness penalty package (smooth_dome, smooth_bin_curve, smooth_yr_diff, smooth_yr_curve, smooth_mean_center), applied here to bicubic fleets ──

  test_that("bicubic penalty terms are zero by default (no behavior change unless opted into)", {
    n_bin_nodes <- 4; n_yr_nodes <- 3
    set.seed(99)
    node_par <- matrix(rnorm(n_bin_nodes * n_yr_nodes, sd = 0.5), nrow = n_yr_nodes, ncol = n_bin_nodes)
    fixed_pars_init <- array(0, dim = c(1, n_bin_nodes * n_yr_nodes, 1, 1, 1))
    fixed_pars_init[1, , 1, 1, 1] <- as.vector(node_par)

    input_list <- build_input_list(
      fish_sel_model = paste0("bicubic_Bin_", n_bin_nodes, "_Yr_", n_yr_nodes, "_Fleet_1"),
      fish_fixed_sel_pars_spec = "est_all",
      fish_fixed_sel_pars = fixed_pars_init
    )

    obj <- RTMB::MakeADFun(
      cmb(SPoRC_rtmb, input_list$data),
      parameters = input_list$par,
      map = input_list$map,
      silent = TRUE
    )
    rep_before <- obj$report(obj$par)

    input_list_pen <- build_input_list(
      fish_sel_model = paste0("bicubic_Bin_", n_bin_nodes, "_Yr_", n_yr_nodes, "_Fleet_1"),
      fish_fixed_sel_pars_spec = "est_all",
      fish_fixed_sel_pars = fixed_pars_init,
      fish_sel_pen_wts = list(
        smooth_bin_curve = 0,
        smooth_yr_diff = 0,
        smooth_yr_curve = 0,
        smooth_dome = 0,
        smooth_mean_center = 0
      )
    )
    obj_pen <- RTMB::MakeADFun(
      cmb(SPoRC_rtmb, input_list_pen$data),
      parameters = input_list_pen$par,
      map = input_list_pen$map,
      silent = TRUE
    )
    rep_after <- obj_pen$report(obj_pen$par)

    expect_equal(rep_before$jnLL, rep_after$jnLL, tolerance = 1e-10)
  })

  test_that("setting bicubic penalty weights matches manual Get_Selex_Smoothness_Penalty computation", {
    n_bin_nodes <- 4; n_yr_nodes <- 3
    set.seed(100)
    node_par <- matrix(rnorm(n_bin_nodes * n_yr_nodes, sd = 0.5), nrow = n_yr_nodes, ncol = n_bin_nodes)
    fixed_pars_init <- array(0, dim = c(1, n_bin_nodes * n_yr_nodes, 1, 1, 1))
    fixed_pars_init[1, , 1, 1, 1] <- as.vector(node_par)

    pen_wts <- list(
      smooth_bin_curve = 10,
      smooth_yr_diff = 1,
      smooth_yr_curve = 300,
      smooth_dome = 30,
      smooth_mean_center = 10000
    )

    build_with_pen <- function(wts) {
      build_input_list(
        fish_sel_model = paste0("bicubic_Bin_", n_bin_nodes, "_Yr_", n_yr_nodes, "_Fleet_1"),
        fish_fixed_sel_pars_spec = "est_all",
        fish_fixed_sel_pars = fixed_pars_init,
        fish_sel_pen_wts = wts
      )
    }

    input_list_zero <- build_with_pen(list(
      smooth_bin_curve = 0,
      smooth_yr_diff = 0,
      smooth_yr_curve = 0,
      smooth_dome = 0,
      smooth_mean_center = 0
    ))
    obj_zero <- RTMB::MakeADFun(
      cmb(SPoRC_rtmb, input_list_zero$data),
      parameters = input_list_zero$par,
      map = input_list_zero$map,
      silent = TRUE
    )
    jnLL_zero <- obj_zero$report(obj_zero$par)$jnLL

    input_list_pen <- build_with_pen(pen_wts)
    obj_pen <- RTMB::MakeADFun(
      cmb(SPoRC_rtmb, input_list_pen$data),
      parameters = input_list_pen$par,
      map = input_list_pen$map,
      silent = TRUE
    )
    jnLL_pen <- obj_pen$report(obj_pen$par)$jnLL

    # manual computation using the exact realized fish_sel surface and the same weights
    n_yrs <- length(input_list_pen$data$years)
    rep_pen <- obj_pen$report(obj_pen$par)
    n_ages_actual <- dim(rep_pen$fish_sel)[5]
    sel_vals <- array(rep_pen$fish_sel[1, 1, 1:n_yrs, 1, , 1, 1], dim = c(1, n_yrs, n_ages_actual, 1, 1))

    expected_penalty <- Get_Selex_Smoothness_Penalty(
      sel_vals,
      wt_bin_curve = pen_wts$smooth_bin_curve,
      wt_yr_diff = pen_wts$smooth_yr_diff,
      wt_yr_curve = pen_wts$smooth_yr_curve,
      wt_dome = pen_wts$smooth_dome,
      wt_mean_center = pen_wts$smooth_mean_center
    )

    expect_equal(jnLL_pen - jnLL_zero, -expected_penalty, tolerance = 1e-6)
  })

  test_that("bicubic penalty is applied unconditionally regardless of cont_tv_fish_sel (no time-varying-deviation gate)", {
    n_bin_nodes <- 4; n_yr_nodes <- 3
    set.seed(101)
    node_par <- matrix(rnorm(n_bin_nodes * n_yr_nodes, sd = 0.5), nrow = n_yr_nodes, ncol = n_bin_nodes)
    fixed_pars_init <- array(0, dim = c(1, n_bin_nodes * n_yr_nodes, 1, 1, 1))
    fixed_pars_init[1, , 1, 1, 1] <- as.vector(node_par)

    # cont_tv_fish_sel defaults to "none" for this fleet (no continuous time-varying deviations at all);
    # the smoothness penalty section is entirely independent of that mechanism.
    input_list <- build_input_list(
      fish_sel_model = paste0("bicubic_Bin_", n_bin_nodes, "_Yr_", n_yr_nodes, "_Fleet_1"),
      fish_fixed_sel_pars_spec = "est_all",
      fish_fixed_sel_pars = fixed_pars_init,
      fish_sel_pen_wts = list(smooth_dome = 30)
    )
    expect_equal(input_list$data$cont_tv_fish_sel[1,1], 0) # confirm "none"

    obj <- RTMB::MakeADFun(
      cmb(SPoRC_rtmb, input_list$data),
      parameters = input_list$par,
      map = input_list$map,
      silent = TRUE
    )
    expect_no_error(obj$fn(obj$par))
    expect_no_error(obj$gr(obj$par))
  })

  test_that("survey SelStyr holds pre-fit years constant and fits the spline only from SelStyr onward", {
    n_bin_nodes <- 4; n_yr_nodes <- 3
    selstyr_year <- 5
    set.seed(203)
    node_par <- matrix(rnorm(n_bin_nodes * n_yr_nodes, sd = 0.5), nrow = n_yr_nodes, ncol = n_bin_nodes)
    fixed_pars_init <- array(0, dim = c(1, n_bin_nodes * n_yr_nodes, 1, 1, 1))
    fixed_pars_init[1, , 1, 1, 1] <- as.vector(node_par)

    input_list <- build_input_list(
      fish_sel_model = "logist2_Fleet_1",
      srv_sel_model = paste0("bicubic_Bin_", n_bin_nodes, "_Yr_", n_yr_nodes, "_SelStyr_", selstyr_year, "_Fleet_1"),
      srv_fixed_sel_pars_spec = "est_all",
      srv_extra_args = list(srv_fixed_sel_pars = fixed_pars_init)
    )

    expect_equal(input_list$data$srv_sel_bicubic_selstyr[1, 1, 1], selstyr_year)

    obj <- RTMB::MakeADFun(
      cmb(SPoRC_rtmb, input_list$data),
      parameters = input_list$par,
      map = input_list$map,
      silent = TRUE
    )
    rep <- obj$report(obj$par)

    n_yrs <- length(input_list$data$years)
    fit_years <- selstyr_year:n_yrs

    Wbin <- input_list$data$srv_sel_bicubic_Wbin[1, , 1:n_bin_nodes, 1, 1]
    Wyr_fit <- Get_Natural_Cubic_Spline_Weights(seq(0, 1, length.out = n_yr_nodes), seq(0, 1, length.out = length(fit_years)))
    expected_fit_surface <- exp(Wyr_fit %*% (node_par %*% t(Wbin)))

    sel_at_selstyr <- as.vector(rep$srv_sel[1, 1, selstyr_year, 1, , 1, 1])
    for (y in 1:(selstyr_year - 1)) {
      expect_equal(as.vector(rep$srv_sel[1, 1, y, 1, , 1, 1]), sel_at_selstyr, tolerance = 1e-8, label = sprintf("pre-SelStyr year %d", y))
    }
    expect_equal(sel_at_selstyr, expected_fit_surface[1, ], tolerance = 1e-8)

    for (y in fit_years) {
      expect_equal(as.vector(rep$srv_sel[1, 1, y, 1, , 1, 1]), expected_fit_surface[y - selstyr_year + 1, ],
                  tolerance = 1e-8, label = sprintf("fit year %d", y))
    }
  })

  # ── Modular smoothness penalty applied to a non-bicubic (nonpar/blocked) fleet ──
  # Bridges an ADMB assessment's "selectivity kept constant between change years" coefficient
  # selectivity (e.g. EBS pollock's compute_fsh_selectivity / sel_devs_fsh) using SPoRC's existing
  # block system + nonpar (Selex_Model == 5) selectivity, with the same modular smoothness-penalty
  # terms used by the bicubic spline now generalized to apply here too (see model_objective.R's
  # "Modular Selectivity Smoothness Penalty" section).

  test_that("nonpar fishery selectivity with discrete blocks can opt into the modular smoothness penalty (pollock-style bridge)", {
    n_ages <- sim_obj$n_ages
    half <- floor(sim_obj$n_years / 2)
    block1_yrs <- paste0(1, "-", half)
    block2_yrs <- paste0(half + 1, "-terminal")

    set.seed(300)
    node_par_1 <- rnorm(n_ages, sd = 0.5)
    node_par_2 <- rnorm(n_ages, sd = 0.5)
    fixed_pars_init <- array(0, dim = c(1, n_ages, 2, 1, 1))
    fixed_pars_init[1, , 1, 1, 1] <- node_par_1
    fixed_pars_init[1, , 2, 1, 1] <- node_par_2

    pen_wts <- list(smooth_bin_curve = 5, smooth_bin_diff = 20, smooth_yr_diff = 2, smooth_dome = 15)

    build_nonpar <- function(wts) {
      build_input_list(
        fish_sel_model = c("nonpar_Block_1_Fleet_1", "nonpar_Block_2_Fleet_1"),
        fish_sel_blocks = c(
          paste0("Block_1_Year_", block1_yrs, "_Fleet_1"),
          paste0("Block_2_Year_", block2_yrs, "_Fleet_1")
        ),
        fish_fixed_sel_pars_spec = "est_all",
        fish_fixed_sel_pars = fixed_pars_init,
        fish_sel_nonpar_est_bins = list(list(as.list(1:n_ages), as.list(1:n_ages))),
        fish_sel_pen_wts = wts
      )
    }

    input_list_zero <- build_nonpar(list(smooth_bin_curve = 0, smooth_bin_diff = 0, smooth_yr_diff = 0, smooth_dome = 0))
    obj_zero <- RTMB::MakeADFun(
      cmb(SPoRC_rtmb, input_list_zero$data),
      parameters = input_list_zero$par,
      map = input_list_zero$map,
      silent = TRUE
    )
    jnLL_zero <- obj_zero$report(obj_zero$par)$jnLL

    input_list_pen <- build_nonpar(pen_wts)
    obj_pen <- RTMB::MakeADFun(
      cmb(SPoRC_rtmb, input_list_pen$data),
      parameters = input_list_pen$par,
      map = input_list_pen$map,
      silent = TRUE
    )
    rep_pen <- obj_pen$report(obj_pen$par)
    jnLL_pen <- rep_pen$jnLL

    n_yrs <- length(input_list_pen$data$years)
    sel_vals <- array(rep_pen$fish_sel[1, 1, 1:n_yrs, 1, , 1, 1], dim = c(1, n_yrs, n_ages, 1, 1))
    expected_penalty <- Get_Selex_Smoothness_Penalty(
      sel_vals,
      wt_bin_curve = pen_wts$smooth_bin_curve,
      wt_bin_diff = pen_wts$smooth_bin_diff,
      wt_yr_diff = pen_wts$smooth_yr_diff,
      wt_dome = pen_wts$smooth_dome
    )

    expect_equal(jnLL_pen - jnLL_zero, -expected_penalty, tolerance = 1e-6)

    # yr_diff is zero within each block (selectivity kept constant) and only bites at the one
    # actual block transition -- exactly mirroring an ADMB "change year" stability penalty
    expect_equal(as.vector(rep_pen$fish_sel[1, 1, 1, 1, , 1, 1]), as.vector(rep_pen$fish_sel[1, 1, half, 1, , 1, 1]))
    expect_equal(as.vector(rep_pen$fish_sel[1, 1, half + 1, 1, , 1, 1]), as.vector(rep_pen$fish_sel[1, 1, sim_obj$n_years, 1, , 1, 1]))
    expect_false(isTRUE(all.equal(as.vector(rep_pen$fish_sel[1, 1, half, 1, , 1, 1]), as.vector(rep_pen$fish_sel[1, 1, half + 1, 1, , 1, 1]))))
  })

  test_that("a nonpar/blocked fleet with no fish_sel_pen_wts set is unaffected by the penalty-section generalization", {
    n_ages <- sim_obj$n_ages
    half <- floor(sim_obj$n_years / 2)
    block1_yrs <- paste0(1, "-", half)
    block2_yrs <- paste0(half + 1, "-terminal")

    set.seed(301)
    fixed_pars_init <- array(0, dim = c(1, n_ages, 2, 1, 1))
    fixed_pars_init[1, , 1, 1, 1] <- rnorm(n_ages, sd = 0.5)
    fixed_pars_init[1, , 2, 1, 1] <- rnorm(n_ages, sd = 0.5)

    build_nonpar <- function(...) {
      build_input_list(
        fish_sel_model = c("nonpar_Block_1_Fleet_1", "nonpar_Block_2_Fleet_1"),
        fish_sel_blocks = c(
          paste0("Block_1_Year_", block1_yrs, "_Fleet_1"),
          paste0("Block_2_Year_", block2_yrs, "_Fleet_1")
        ),
        fish_fixed_sel_pars_spec = "est_all",
        fish_fixed_sel_pars = fixed_pars_init,
        fish_sel_nonpar_est_bins = list(list(as.list(1:n_ages), as.list(1:n_ages))),
        ...
      )
    }

    input_list_default <- build_nonpar() # fish_sel_pen_wts left at its NULL default
    obj_default <- RTMB::MakeADFun(
      cmb(SPoRC_rtmb, input_list_default$data),
      parameters = input_list_default$par,
      map = input_list_default$map,
      silent = TRUE
    )

    input_list_explicit_zero <- build_nonpar(fish_sel_pen_wts = list(
      smooth_bin_curve = 0,
      smooth_bin_diff = 0,
      smooth_yr_diff = 0,
      smooth_yr_curve = 0,
      smooth_dome = 0,
      smooth_mean_center = 0
    ))
    obj_explicit_zero <- RTMB::MakeADFun(
      cmb(SPoRC_rtmb, input_list_explicit_zero$data),
      parameters = input_list_explicit_zero$par,
      map = input_list_explicit_zero$map,
      silent = TRUE
    )

    expect_equal(obj_default$report(obj_default$par)$jnLL, obj_explicit_zero$report(obj_explicit_zero$par)$jnLL, tolerance = 1e-10)
  })

  # ── Retention selectivity: same bicubic wiring, mirrored for ret_sel_model ───

  test_that("retention bicubic model builds, evaluates, and differentiates without error with a single block", {
    n_bin_nodes <- 4; n_yr_nodes <- 3
    node_par <- matrix(
      c(-1.1, 0.2, 0.8, -0.4,
                         0.5, -0.6, 0.3, 1.0,
                         -0.2, 0.7, -0.5, 0.4),
      nrow = n_yr_nodes,
      ncol = n_bin_nodes,
      byrow = TRUE
    )
    fixed_pars_init <- array(0, dim = c(1, n_bin_nodes * n_yr_nodes, 1, 1, 1))
    fixed_pars_init[1, , 1, 1, 1] <- as.vector(node_par)

    input_list <- build_input_list(
      fish_sel_model = "logist2_Fleet_1",
      ret_sel_model = paste0("bicubic_Bin_", n_bin_nodes, "_Yr_", n_yr_nodes, "_Fleet_1"),
      ret_fixed_sel_pars_spec = "est_all",
      ret_fixed_sel_pars = fixed_pars_init,
      use_fixed_ret_sel = 0
    )

    obj <- RTMB::MakeADFun(
      cmb(SPoRC_rtmb, input_list$data),
      parameters = input_list$par,
      map = input_list$map,
      silent = TRUE
    )

    expect_no_error(obj$fn(obj$par))
    expect_no_error(obj$gr(obj$par))

    rep <- obj$report(obj$par)

    # manual computation using the exact weight matrices Setup_Mod_Fishsel_and_Q constructed
    Wbin <- input_list$data$ret_sel_bicubic_Wbin[1, , 1:n_bin_nodes, 1, 1]
    Wyr  <- input_list$data$ret_sel_bicubic_Wyr[1, , 1:n_yr_nodes, 1, 1]
    expected_surface <- exp(Wyr %*% (node_par %*% t(Wbin))) # n_yrs_total x n_ages

    n_yrs <- length(input_list$data$years)
    for (y in 1:n_yrs) {
      expect_equal(as.vector(rep$ret_sel[1, 1, y, 1, , 1, 1]), expected_surface[y, ],
                  tolerance = 1e-8, label = sprintf("ret_sel year %d", y))
    }
  })

  test_that("retention bicubic n_yr_nodes == 1 gives a time-invariant selectivity surface end-to-end", {
    n_bin_nodes <- 5
    node_par <- c(-0.7, 0.1, 0.9, -0.3, 0.5)
    fixed_pars_init <- array(0, dim = c(1, n_bin_nodes, 1, 1, 1))
    fixed_pars_init[1, , 1, 1, 1] <- node_par

    input_list <- build_input_list(
      fish_sel_model = "logist2_Fleet_1",
      ret_sel_model = paste0("bicubic_Bin_", n_bin_nodes, "_Yr_1_Fleet_1"),
      ret_fixed_sel_pars_spec = "est_all",
      ret_fixed_sel_pars = fixed_pars_init,
      use_fixed_ret_sel = 0
    )

    obj <- RTMB::MakeADFun(
      cmb(SPoRC_rtmb, input_list$data),
      parameters = input_list$par,
      map = input_list$map,
      silent = TRUE
    )
    rep <- obj$report(obj$par)

    n_yrs <- length(input_list$data$years)
    sel_y1 <- as.vector(rep$ret_sel[1, 1, 1, 1, , 1, 1])
    for (y in 2:n_yrs) {
      expect_equal(as.vector(rep$ret_sel[1, 1, y, 1, , 1, 1]), sel_y1, tolerance = 1e-8)
    }

    Wbin <- input_list$data$ret_sel_bicubic_Wbin[1, , 1:n_bin_nodes, 1, 1]
    expect_equal(sel_y1, exp(as.vector(Wbin %*% node_par)), tolerance = 1e-8)
  })

  test_that("retention SelStyr holds pre-fit years constant and fits the spline only from SelStyr onward", {
    n_bin_nodes <- 4; n_yr_nodes <- 3
    selstyr_year <- 5 # years 1:4 should be edge-kept; the spline is fit only over years 5:n_years
    set.seed(303)
    node_par <- matrix(rnorm(n_bin_nodes * n_yr_nodes, sd = 0.5), nrow = n_yr_nodes, ncol = n_bin_nodes)
    fixed_pars_init <- array(0, dim = c(1, n_bin_nodes * n_yr_nodes, 1, 1, 1))
    fixed_pars_init[1, , 1, 1, 1] <- as.vector(node_par)

    input_list <- build_input_list(
      fish_sel_model = "logist2_Fleet_1",
      ret_sel_model = paste0("bicubic_Bin_", n_bin_nodes, "_Yr_", n_yr_nodes, "_SelStyr_", selstyr_year, "_Fleet_1"),
      ret_fixed_sel_pars_spec = "est_all",
      ret_fixed_sel_pars = fixed_pars_init,
      use_fixed_ret_sel = 0
    )

    expect_equal(input_list$data$ret_sel_bicubic_selstyr[1, 1, 1], selstyr_year)

    obj <- RTMB::MakeADFun(
      cmb(SPoRC_rtmb, input_list$data),
      parameters = input_list$par,
      map = input_list$map,
      silent = TRUE
    )
    rep <- obj$report(obj$par)

    n_yrs <- length(input_list$data$years)
    fit_years <- selstyr_year:n_yrs

    Wbin <- input_list$data$ret_sel_bicubic_Wbin[1, , 1:n_bin_nodes, 1, 1]
    Wyr_fit <- Get_Natural_Cubic_Spline_Weights(seq(0, 1, length.out = n_yr_nodes), seq(0, 1, length.out = length(fit_years)))
    expected_fit_surface <- exp(Wyr_fit %*% (node_par %*% t(Wbin)))

    # years before SelStyr are all identical to the SelStyr year's fitted curve
    sel_at_selstyr <- as.vector(rep$ret_sel[1, 1, selstyr_year, 1, , 1, 1])
    for (y in 1:(selstyr_year - 1)) {
      expect_equal(as.vector(rep$ret_sel[1, 1, y, 1, , 1, 1]), sel_at_selstyr, tolerance = 1e-8, label = sprintf("pre-SelStyr year %d", y))
    }
    expect_equal(sel_at_selstyr, expected_fit_surface[1, ], tolerance = 1e-8)

    # years from SelStyr onward follow the spline fit over the sub-range only
    for (y in fit_years) {
      expect_equal(as.vector(rep$ret_sel[1, 1, y, 1, , 1, 1]), expected_fit_surface[y - selstyr_year + 1, ],
                  tolerance = 1e-8, label = sprintf("fit year %d", y))
    }
  })

  test_that("retention NSelBins holds post-fit bins constant (plateau) and fits the spline only over the first NSelBins bins", {
    n_bin_nodes <- 4; n_yr_nodes <- 1
    n_fit_bins <- sim_obj$n_ages - 2 # ages 1:(n_ages-2) are actually fit; the last 2 ages plateau
    set.seed(404)
    node_par <- rnorm(n_bin_nodes, sd = 0.5)
    fixed_pars_init <- array(0, dim = c(1, n_bin_nodes, 1, 1, 1))
    fixed_pars_init[1, , 1, 1, 1] <- node_par

    input_list <- build_input_list(
      fish_sel_model = "logist2_Fleet_1",
      ret_sel_model = paste0("bicubic_Bin_", n_bin_nodes, "_Yr_", n_yr_nodes, "_NSelBins_", n_fit_bins, "_Fleet_1"),
      ret_fixed_sel_pars_spec = "est_all",
      ret_fixed_sel_pars = fixed_pars_init,
      use_fixed_ret_sel = 0
    )

    expect_equal(input_list$data$ret_sel_bicubic_nselbins[1, 1, 1], n_fit_bins)

    obj <- RTMB::MakeADFun(
      cmb(SPoRC_rtmb, input_list$data),
      parameters = input_list$par,
      map = input_list$map,
      silent = TRUE
    )
    rep <- obj$report(obj$par)

    n_ages <- sim_obj$n_ages
    sel_y1 <- as.vector(rep$ret_sel[1, 1, 1, 1, , 1, 1])

    # bins beyond n_fit_bins plateau at the last fitted bin's value
    for (b in (n_fit_bins + 1):n_ages) {
      expect_equal(sel_y1[b], sel_y1[n_fit_bins], tolerance = 1e-8, label = sprintf("plateau bin %d", b))
    }

    # bins 1:n_fit_bins follow the spline fit over the sub-range only
    Wbin_fit <- Get_Natural_Cubic_Spline_Weights(seq(0, 1, length.out = n_bin_nodes), seq(0, 1, length.out = n_fit_bins))
    expected_fit <- exp(as.vector(Wbin_fit %*% node_par))
    expect_equal(sel_y1[1:n_fit_bins], expected_fit, tolerance = 1e-8)
  })

  test_that("retention, survey, and fishery fleets can each independently use bicubic selectivity at once", {
    n_bin_nodes <- 4
    fish_node_par <- c(-0.5, 0.5, 0.9, -0.3)
    srv_node_par  <- c(0.6, -0.8, 0.2, 1.0)
    ret_node_par  <- c(-0.9, 0.4, -0.2, 0.7)

    fish_fixed_pars_init <- array(0, dim = c(1, n_bin_nodes, 1, 1, 1))
    fish_fixed_pars_init[1, , 1, 1, 1] <- fish_node_par
    srv_fixed_pars_init <- array(0, dim = c(1, n_bin_nodes, 1, 1, 1))
    srv_fixed_pars_init[1, , 1, 1, 1] <- srv_node_par
    ret_fixed_pars_init <- array(0, dim = c(1, n_bin_nodes, 1, 1, 1))
    ret_fixed_pars_init[1, , 1, 1, 1] <- ret_node_par

    input_list <- build_input_list(
      fish_sel_model = paste0("bicubic_Bin_", n_bin_nodes, "_Yr_1_Fleet_1"),
      fish_fixed_sel_pars_spec = "est_all",
      fish_fixed_sel_pars = fish_fixed_pars_init,
      srv_sel_model = paste0("bicubic_Bin_", n_bin_nodes, "_Yr_1_Fleet_1"),
      srv_fixed_sel_pars_spec = "est_all",
      srv_extra_args = list(srv_fixed_sel_pars = srv_fixed_pars_init),
      ret_sel_model = paste0("bicubic_Bin_", n_bin_nodes, "_Yr_1_Fleet_1"),
      ret_fixed_sel_pars_spec = "est_all",
      ret_fixed_sel_pars = ret_fixed_pars_init,
      use_fixed_ret_sel = 0
    )

    obj <- RTMB::MakeADFun(
      cmb(SPoRC_rtmb, input_list$data),
      parameters = input_list$par,
      map = input_list$map,
      silent = TRUE
    )
    expect_no_error(obj$fn(obj$par))
    rep <- obj$report(obj$par)

    Wbin_fish <- input_list$data$fish_sel_bicubic_Wbin[1, , 1:n_bin_nodes, 1, 1]
    Wbin_srv  <- input_list$data$srv_sel_bicubic_Wbin[1, , 1:n_bin_nodes, 1, 1]
    Wbin_ret  <- input_list$data$ret_sel_bicubic_Wbin[1, , 1:n_bin_nodes, 1, 1]

    expect_equal(as.vector(rep$fish_sel[1, 1, 1, 1, , 1, 1]), exp(as.vector(Wbin_fish %*% fish_node_par)), tolerance = 1e-8)
    expect_equal(as.vector(rep$srv_sel[1, 1, 1, 1, , 1, 1]), exp(as.vector(Wbin_srv %*% srv_node_par)), tolerance = 1e-8)
    expect_equal(as.vector(rep$ret_sel[1, 1, 1, 1, , 1, 1]), exp(as.vector(Wbin_ret %*% ret_node_par)), tolerance = 1e-8)
  })

  test_that("setting retention bicubic penalty weights matches manual Get_Selex_Smoothness_Penalty computation", {
    n_bin_nodes <- 4; n_yr_nodes <- 3
    set.seed(505)
    node_par <- matrix(rnorm(n_bin_nodes * n_yr_nodes, sd = 0.5), nrow = n_yr_nodes, ncol = n_bin_nodes)
    fixed_pars_init <- array(0, dim = c(1, n_bin_nodes * n_yr_nodes, 1, 1, 1))
    fixed_pars_init[1, , 1, 1, 1] <- as.vector(node_par)

    pen_wts <- list(
      smooth_bin_curve = 10,
      smooth_yr_diff = 1,
      smooth_yr_curve = 300,
      smooth_dome = 30,
      smooth_mean_center = 10000
    )

    build_with_pen <- function(wts) {
      build_input_list(
        fish_sel_model = "logist2_Fleet_1",
        ret_sel_model = paste0("bicubic_Bin_", n_bin_nodes, "_Yr_", n_yr_nodes, "_Fleet_1"),
        ret_fixed_sel_pars_spec = "est_all",
        ret_fixed_sel_pars = fixed_pars_init,
        use_fixed_ret_sel = 0,
        ret_sel_pen_wts = wts
      )
    }

    input_list_zero <- build_with_pen(list(
      smooth_bin_curve = 0,
      smooth_yr_diff = 0,
      smooth_yr_curve = 0,
      smooth_dome = 0,
      smooth_mean_center = 0
    ))
    obj_zero <- RTMB::MakeADFun(
      cmb(SPoRC_rtmb, input_list_zero$data),
      parameters = input_list_zero$par,
      map = input_list_zero$map,
      silent = TRUE
    )
    jnLL_zero <- obj_zero$report(obj_zero$par)$jnLL

    input_list_pen <- build_with_pen(pen_wts)
    obj_pen <- RTMB::MakeADFun(
      cmb(SPoRC_rtmb, input_list_pen$data),
      parameters = input_list_pen$par,
      map = input_list_pen$map,
      silent = TRUE
    )
    rep_pen <- obj_pen$report(obj_pen$par)

    n_yrs <- length(input_list_pen$data$years)
    Wbin <- input_list_pen$data$ret_sel_bicubic_Wbin[1, , 1:n_bin_nodes, 1, 1]
    Wyr  <- input_list_pen$data$ret_sel_bicubic_Wyr[1, , 1:n_yr_nodes, 1, 1]
    surface <- exp(Wyr %*% (node_par %*% t(Wbin)))
    sel_vals <- array(surface, dim = c(1, n_yrs, ncol(surface), 1, 1))

    manual_penalty <- Get_Selex_Smoothness_Penalty(sel_vals,
                                                   wt_bin_curve = pen_wts$smooth_bin_curve,
                                                   wt_yr_diff = pen_wts$smooth_yr_diff,
                                                   wt_yr_curve = pen_wts$smooth_yr_curve,
                                                   wt_dome = pen_wts$smooth_dome,
                                                   wt_mean_center = pen_wts$smooth_mean_center,
                                                   normalize = TRUE)

    expect_equal(rep_pen$jnLL - jnLL_zero, -manual_penalty, tolerance = 1e-6)
  })

})
