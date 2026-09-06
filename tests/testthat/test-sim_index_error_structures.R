library(SPoRC)
library(testthat)

# The simulator draws each index data source under the fleet's estimation-model error
# structure. These tests pin the three layers separately: the factor routines
# that turns a fixed covariance into per-cell draw parameters, the Setup_Sim_*
# validation and storage, and the model-side population-specific blocks that
# share the fleet's LikeType.

test_that("cov_to_factor preserves marginal scale and approximates a one-factor correlation", {

  d <- c(2, 3, 1.5, 2.5, 4)
  lambda <- rep(0.8, 5)
  R <- outer(lambda, lambda); diag(R) <- 1
  S <- outer(d, d) * R
  fac <- SPoRC:::cov_to_factor(S)

  test_that("the marginal scale is the covariance diagonal exactly", {
    expect_equal(fac$d, d, tolerance = 1e-12)
  })

  test_that("loadings stay inside (-1, 1) with a non-negative mean", {
    expect_true(all(abs(fac$lambda) <= 0.99))
    expect_true(mean(fac$lambda) >= 0)
  })

  test_that("the implied correlations approximate the true one-factor correlation", {
    R_hat <- outer(fac$lambda, fac$lambda); diag(R_hat) <- 1
    expect_lt(max(abs(R_hat - R)[upper.tri(R)]), 0.1)
  })

})

test_that("build_idx_factor positions covariance rows the way the model collects observations", {

  # gappy use flags across two regions: the model scans use == 1 in array order
  # (region fastest, then year, then season), so the covariance rows must too
  n_regions <- 2; n_yrs <- 5; n_seas <- 1; n_fleets <- 2
  use <- array(0, dim = c(n_regions, n_yrs, n_seas, n_fleets))
  use[1, c(1, 3, 4), 1, 2] <- 1
  use[2, c(3, 5), 1, 2] <- 1
  S <- diag(seq(1, 2, length.out = 5)^2)

  fac <- SPoRC:::build_idx_factor(list(NULL, S), c(0, 2), use, n_fleets, "SrvIdx_Cov")

  test_that("a non-mvn fleet gets no decomposition", {
    expect_null(fac[[1]])
  })

  test_that("rows follow which(use == 1): region varies fastest, then year", {
    expect_equal(fac[[2]]$row[1, 1, 1], 1)
    expect_equal(fac[[2]]$row[1, 3, 1], 2)
    expect_equal(fac[[2]]$row[2, 3, 1], 3)
    expect_equal(fac[[2]]$row[1, 4, 1], 4)
    expect_equal(fac[[2]]$row[2, 5, 1], 5)
    expect_true(is.na(fac[[2]]$row[2, 1, 1]))
  })

  test_that("the marginal scale is read off the diagonal in row order", {
    expect_equal(fac[[2]]$d, seq(1, 2, length.out = 5), tolerance = 1e-12)
  })

  test_that("a wrong-size covariance is rejected against the use flags", {
    expect_error(SPoRC:::build_idx_factor(list(NULL, diag(3)), c(0, 2), use, n_fleets, "SrvIdx_Cov"),
                 "observations")
  })

})

test_that("resolve_idx_factor falls back to the mean parameters outside the covariance", {

  use <- array(0, dim = c(1, 4, 1, 1))
  use[1, c(1, 2, 4), 1, 1] <- 1
  S <- diag(c(1, 4, 9))
  fac <- SPoRC:::build_idx_factor(list(S), 2, use, 1, "SrvIdx_Cov")[[1]]

  test_that("a cell inside the covariance gets its own scale and loading", {
    got <- SPoRC:::resolve_idx_factor(fac, 1, 2, 1)
    expect_equal(got$d, fac$d[2])
    expect_equal(got$lambda, fac$lambda[2])
  })

  test_that("an unused cell inside the data years gets the means", {
    got <- SPoRC:::resolve_idx_factor(fac, 1, 3, 1)
    expect_equal(got$d, mean(fac$d))
    expect_equal(got$lambda, mean(fac$lambda))
  })

  test_that("a projection year past the covariance gets the means", {
    got <- SPoRC:::resolve_idx_factor(fac, 1, 9, 1)
    expect_equal(got$d, mean(fac$d))
    expect_equal(got$lambda, mean(fac$lambda))
  })

})

test_that("draw_index_obs draws each error structure from the same seed formulae", {

  true <- c(100, 120); se <- c(0.2, 0.3)

  test_that("lognormal multiplies by exp of a log-scale normal", {
    set.seed(5); got <- SPoRC:::draw_index_obs(true, se, 0)
    set.seed(5); expect_equal(got, true * exp(stats::rnorm(2, 0, se)))
  })

  test_that("normal adds an arithmetic-scale normal", {
    set.seed(5); got <- SPoRC:::draw_index_obs(true, se, 1)
    set.seed(5); expect_equal(got, true + stats::rnorm(2, 0, se))
  })

  test_that("mvn combines the shared factor and an independent residual", {
    d <- c(10, 12); lambda <- c(0.6, 0.7); u <- 1.5
    set.seed(5); got <- SPoRC:::draw_index_obs(true, NA, 2, d = d, lambda = lambda, u = u)
    set.seed(5); e <- stats::rnorm(2)
    expect_equal(got, true + d * (lambda * u + sqrt(1 - lambda^2) * e))
  })

  test_that("mvn refuses to draw without the factor parameters", {
    expect_error(SPoRC:::draw_index_obs(true, se, 2), "multivariate normal")
  })

})

test_that("Setup_Sim_Fishing and Setup_Sim_Survey validate and store the mvn routines", {

  n_yrs <- 4; n_ages <- 3; n_sims <- 2
  sim_list <- Setup_Sim_Dim(
    n_sims = n_sims,
    n_yrs = n_yrs,
    n_regions = 1,
    n_ages = n_ages,
    n_lens = NULL,
    n_sexes = 1,
    n_fish_fleets = 1,
    n_srv_fleets = 1,
    n_pop = 1
  )
  sim_list <- Setup_Sim_Containers(sim_list)
  sel <- replicate(n_sims, array(1, dim = c(1, 1, n_yrs, 1, n_ages, 1, 1)))

  use <- array(0, dim = c(1, n_yrs, 1, 1)); use[1, 1:3, 1, 1] <- 1
  S <- diag(3) * 4

  test_that("the lognormal default stores zero codes and no factor routines", {
    sl <- Setup_Sim_Fishing(sim_list = sim_list, fish_sel_input = sel)
    expect_equal(sl$FishIdx_LikeType, 0)
    expect_null(sl$fish_idx_mvn)
    expect_null(sl$fish_idx_u)
  })

  test_that("character codes convert like the estimation model's", {
    sl <- Setup_Sim_Fishing(sim_list = sim_list, fish_sel_input = sel, FishIdx_LikeType = "normal")
    expect_equal(sl$FishIdx_LikeType, 1)
  })

  test_that("an mvn fleet without use flags or with a wrong-size covariance is rejected", {
    expect_error(Setup_Sim_Fishing(
      sim_list = sim_list,
      fish_sel_input = sel,
      FishIdx_LikeType = "mvn",
      FishIdx_Cov = list(S)
    ),
                 "UseFishIdx")
    expect_error(Setup_Sim_Fishing(
      sim_list = sim_list,
      fish_sel_input = sel,
      FishIdx_LikeType = "mvn",
      FishIdx_Cov = list(diag(2)),
      UseFishIdx = use
    ),
                 "observations")
  })

  test_that("an mvn fleet stores factor parameters and the shared-draw container", {
    sl <- Setup_Sim_Fishing(
      sim_list = sim_list,
      fish_sel_input = sel,
      FishIdx_LikeType = "mvn",
      FishIdx_Cov = list(S),
      UseFishIdx = use
    )
    expect_equal(sl$FishIdx_LikeType, 2)
    expect_equal(sl$fish_idx_mvn[[1]]$d, rep(2, 3))
    expect_equal(dim(sl$fish_idx_u), c(1, n_sims))
    expect_true(all(is.na(sl$fish_idx_u)))
  })

  test_that("the survey side mirrors the fishery side", {
    sl <- Setup_Sim_Survey(
      sim_list = sim_list,
      srv_sel_input = sel,
      SrvIdx_LikeType = "mvn",
      SrvIdx_Cov = list(S),
      UseSrvIdx = use
    )
    expect_equal(sl$SrvIdx_LikeType, 2)
    expect_equal(sl$srv_idx_mvn[[1]]$d, rep(2, 3))
    expect_equal(dim(sl$srv_idx_u), c(1, n_sims))
    expect_error(Setup_Sim_Survey(
      sim_list = sim_list,
      srv_sel_input = sel,
      SrvIdx_LikeType = "mvn",
      SrvIdx_Cov = list(S)
    ),
                 "UseSrvIdx")
  })

})

`%||%` <- function(a, b) if(is.null(a)) b else a

# One small operating model shared by the integration tests below: 1 population,
# 1 region, 15 years, 6 ages, 2 replicates. Simulated once per survey error
# structure and cached.
index_error_om <- local({
  cached <- list()

  function(
    SrvIdx_LikeType = NULL,
    SrvIdx_Cov = NULL,
    UseSrvIdx = NULL,
    drop_new_fields = FALSE,
    seed = 321
  ) {

    key <- paste(SrvIdx_LikeType %||% "default", drop_new_fields, seed, sep = "_")
    if(!is.null(cached[[key]])) return(cached[[key]])

    n_yrs <- 15; n_ages <- 6; n_sims <- 2
    sim_list <- Setup_Sim_Dim(
      n_sims = n_sims,
      n_yrs = n_yrs,
      n_regions = 1,
      n_ages = n_ages,
      n_lens = NULL,
      n_sexes = 1,
      n_fish_fleets = 1,
      n_srv_fleets = 1,
      n_pop = 1
    )
    sim_list <- Setup_Sim_Containers(sim_list)

    curve7 <- function(slope, infl, scale = 1) {
      array(rep(scale / (1 + exp(-slope * ((1:n_ages) - infl))), each = n_yrs),
            dim = c(1, 1, n_yrs, 1, n_ages, 1, 1))
    }

    sim_list <- Setup_Sim_Fishing(
      sim_list = sim_list,
      fish_sel_input = replicate(n_sims, curve7(3, 2)),
      ret_sel_input = replicate(n_sims, curve7(3, 2, scale = 1)),
      dmr_input = array(0, dim = c(1, n_yrs, 1, 1, n_sims))
    )

    srv_args <- list(sim_list = sim_list, srv_sel_input = replicate(n_sims, curve7(1, 3)))
    if(!is.null(SrvIdx_LikeType)) srv_args$SrvIdx_LikeType <- SrvIdx_LikeType
    if(!is.null(SrvIdx_Cov)) srv_args$SrvIdx_Cov <- SrvIdx_Cov
    if(!is.null(UseSrvIdx)) srv_args$UseSrvIdx <- UseSrvIdx
    sim_list <- do.call(Setup_Sim_Survey, srv_args)

    biol6 <- function(val) array(rep(val, each = n_yrs), dim = c(1, 1, n_yrs, 1, n_ages, 1))
    waa <- 5 / (1 + exp(-3 * ((1:n_ages) - 3)))
    suppressWarnings(sim_list <- Setup_Sim_Biologicals(
      sim_list = sim_list,
      natmort_input = replicate(n_sims, array(0.3, dim = c(1, 1, n_yrs, n_ages, 1))),
      WAA_input = replicate(n_sims, biol6(waa)),
      WAA_fish_input = replicate(n_sims, array(rep(waa, each = n_yrs), dim = c(1, 1, n_yrs, 1, n_ages, 1, 1))),
      WAA_srv_input = replicate(n_sims, array(rep(waa, each = n_yrs), dim = c(1, 1, n_yrs, 1, n_ages, 1, 1))),
      MatAA_input = replicate(n_sims, biol6(1 / (1 + exp(-3 * ((1:n_ages) - 3)))))
    ))

    sim_list <- Setup_Sim_Tagging(sim_list = sim_list, use_conv_fish_tagging = 0)
    sim_list$Movement <- array(1, dim = c(1, 1, 1, n_yrs, 1, n_ages, 1, n_sims))
    sim_list <- Setup_Sim_Rec(
      sim_list = sim_list,
      R0_input = replicate(n_sims, array(5, dim = c(1, 1, n_yrs))),
      ln_sigmaR = array(log(0.5), dim = c(2, 1, 1)),
      recruitment_opt = "mean_rec",
      init_age_strc = 1
    )

    if(drop_new_fields) {
      sim_list$SrvIdx_LikeType <- NULL
      sim_list$FishIdx_LikeType <- NULL
      sim_list$srv_idx_mvn <- NULL
      sim_list$srv_idx_u <- NULL
      sim_list$fish_idx_mvn <- NULL
      sim_list$fish_idx_u <- NULL
    }

    set.seed(seed)
    out <- Simulate_Pop_Static(sim_list = sim_list, output_path = NULL)
    cached[[key]] <<- out
    out
  }
})

test_that("a simulation list without the new fields reproduces the lognormal draws exactly", {

  om_base <- index_error_om()
  om_old <- index_error_om(drop_new_fields = TRUE)

  # byte-identical draws mean the RNG data source is untouched for existing workflows
  expect_identical(om_base$ObsSrvIdx, om_old$ObsSrvIdx)
  expect_identical(om_base$ObsFishIdx, om_old$ObsFishIdx)
  expect_identical(om_base$ObsSrvIdx_pop, om_old$ObsSrvIdx_pop)
  expect_identical(om_base$ObsFishIdx_pop, om_old$ObsFishIdx_pop)

})

test_that("an mvn survey fleet draws from the covariance with a shared factor per replicate", {

  n_yrs <- 15
  use <- array(0, dim = c(1, n_yrs, 1, 1)); use[1, , 1, 1] <- 1
  lambda <- rep(0.9, n_yrs)
  d <- seq(20, 40, length.out = n_yrs)
  R <- outer(lambda, lambda); diag(R) <- 1
  S <- outer(d, d) * R

  om <- index_error_om(SrvIdx_LikeType = "mvn", SrvIdx_Cov = list(S), UseSrvIdx = use)

  test_that("draws are finite and seed-deterministic", {
    expect_true(all(is.finite(om$ObsSrvIdx)))
    om2 <- index_error_om(SrvIdx_LikeType = "mvn", SrvIdx_Cov = list(S), UseSrvIdx = use, seed = 321)
    expect_identical(om$ObsSrvIdx, om2$ObsSrvIdx)
  })

  test_that("the shared factor compresses within-replicate spread of standardized deviations", {
    # z = lambda * u + sqrt(1 - lambda^2) * e, so var(z) within a replicate is
    # about 1 - 0.9^2 = 0.19 rather than the 1 independent draws would give
    z <- (om$ObsSrvIdx[1, , 1, 1, ] - om$TrueSrvIdx[1, , 1, 1, ]) / d
    expect_true(all(apply(z, 2, stats::var) < 0.6))
  })

  test_that("the population-specific data source keeps lognormal error", {
    expect_true(all(om$ObsSrvIdx_pop > 0))
  })

})

test_that("the population-specific index blocks honor the fleet's LikeType", {

  om <- index_error_om()
  sim_data <- simulation_data_to_SPoRC(sim_env = om, y = om$n_years, sim = 1)
  n_yrs <- om$n_years; n_ages <- om$n_ages

  build_input <- function(SrvIdx_LikeType = "lognormal", SrvIdx_Cov = NULL) {
    input_list <- Setup_Mod_Dim(
      years = 1:n_yrs,
      ages = 1:n_ages,
      lens = om$n_lens,
      n_regions = 1,
      n_sexes = 1,
      n_fish_fleets = 1,
      n_srv_fleets = 1,
      n_pop = 1,
      natal_region = om$natal_region,
      verbose = FALSE
    )
    input_list <- Setup_Mod_Rec(
      input_list = input_list,
      do_rec_bias_ramp = 0,
      sigmaR_switch = 1,
      ln_sigmaR = array(log(0.5), c(2, 1, 1)),
      rec_model = "mean_rec",
      sigmaR_spec = "fix",
      init_age_strc = 1,
      equil_init_age_strc = 2,
      ln_global_R0 = log(5)
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
      Fixed_natmort = array(0.3, dim = c(1, 1, n_yrs, n_ages, 1))
    )
    input_list <- Setup_Mod_Tagging(input_list = input_list, use_conv_fish_tagging = 0)
    input_list <- Setup_Mod_Movement(
      input_list = input_list,
      use_fixed_movement = 1,
      Fixed_Movement = NA,
      do_recruits_move = 0
    )
    suppressWarnings(input_list <- Setup_Mod_Catch_and_F(
      input_list = input_list,
      ObsCatch = sim_data$ObsCatch,
      UseCatch = sim_data$UseCatch,
      Use_F_pen = 1,
      sigmaC_spec = "fix",
      ln_sigmaC = sim_data$ln_sigmaC,
      ln_sigmaF = array(log(1), dim = c(1, 1, 1)),
      ObsDiscard = sim_data$ObsDiscard,
      UseDiscard = sim_data$UseDiscard,
      sigma_dmr_spec = "fix",
      dmr_mean_spec = "est_all",
      ln_sigmaD = sim_data$ln_sigmaD
    ))
    input_list <- Setup_Mod_FishIdx_and_Comps(
      input_list = input_list,
      ObsFishIdx = sim_data$ObsFishIdx,
      ObsFishIdx_SE = sim_data$ObsFishIdx_SE,
      UseFishIdx = array(0, dim = dim(sim_data$UseFishIdx)),
      ObsFishAgeComps = sim_data$ObsFishAgeComps,
      ObsFishLenComps = NULL,
      UseFishAgeComps = sim_data$UseFishAgeComps,
      UseFishLenComps = array(0, dim = dim(sim_data$UseFishAgeComps)),
      ISS_FishAgeComps = sim_data$ISS_FishAgeComps,
      ISS_FishLenComps = NULL,
      fish_idx_type = "biom",
      FishAgeComps_LikeType = "Multinomial",
      FishLenComps_LikeType = "none",
      FishAgeComps_Type = "agg_Year_1-terminal_Fleet_1",
      FishLenComps_Type = "none_Year_1-terminal_Fleet_1"
    )
    input_list <- Setup_Mod_SrvIdx_and_Comps(
      input_list = input_list,
      ObsSrvIdx = sim_data$ObsSrvIdx,
      ObsSrvIdx_SE = sim_data$ObsSrvIdx_SE,
      UseSrvIdx = sim_data$UseSrvIdx,
      ObsSrvIdx_pop = sim_data$ObsSrvIdx_pop,
      ObsSrvIdx_pop_SE = sim_data$ObsSrvIdx_pop_SE,
      UseSrvIdx_pop = sim_data$UseSrvIdx_pop,
      SrvIdx_LikeType = SrvIdx_LikeType,
      SrvIdx_Cov = SrvIdx_Cov,
      ObsSrvAgeComps = sim_data$ObsSrvAgeComps,
      ObsSrvLenComps = NULL,
      UseSrvAgeComps = sim_data$UseSrvAgeComps,
      UseSrvLenComps = array(0, dim = dim(sim_data$UseSrvAgeComps)),
      ISS_SrvAgeComps = sim_data$ISS_SrvAgeComps,
      ISS_SrvLenComps = NULL,
      srv_idx_type = "biom",
      SrvAgeComps_LikeType = "Multinomial",
      SrvLenComps_LikeType = "none",
      SrvAgeComps_Type = "agg_Year_1-terminal_Fleet_1",
      SrvLenComps_Type = "none_Year_1-terminal_Fleet_1"
    )
    input_list <- Setup_Mod_Fishsel_and_Q(
      input_list = input_list,
      fish_sel_model = "logist1_Fleet_1",
      fish_fixed_sel_pars_spec = "est_all",
      fish_q_spec = "est_all",
      use_fixed_ret_sel = 1
    )
    input_list <- Setup_Mod_Srvsel_and_Q(
      input_list = input_list,
      srv_sel_model = "logist1_Fleet_1",
      srv_fixed_sel_pars_spec = "est_all",
      srv_q_spec = "est_all"
    )
    input_list <- Setup_Mod_Weighting(
      input_list = input_list,
      Wt_Catch = 1,
      Wt_FishIdx = 1,
      Wt_SrvIdx = 1,
      Wt_SrvIdx_pop = 1,
      Wt_Rec = 1,
      Wt_F = 1,
      Wt_FishAgeComps = array(1, dim = c(1, n_yrs, 1, 1, 1)),
      Wt_SrvAgeComps = array(1, dim = c(1, n_yrs, 1, 1, 1))
    )
    input_list
  }

  report_of <- function(input_list) {
    obj <- fit_model(
      input_list$data,
      input_list$par,
      input_list$map,
      random = NULL,
      silent = TRUE,
      do_optim = FALSE
    )
    list(rep = obj$report(), data = input_list$data)
  }

  obs_pop <- sim_data$ObsSrvIdx_pop[1, 1, , 1, 1]
  se_pop <- sim_data$ObsSrvIdx_pop_SE[1, 1, , 1, 1]
  obs_reg <- sim_data$ObsSrvIdx[1, , 1, 1]
  se_reg <- sim_data$ObsSrvIdx_SE[1, , 1, 1]

  test_that("lognormal fleets keep the previous population-specific likelihood", {
    m <- report_of(build_input("lognormal"))
    addto <- m$data$addtosrvidx
    pred <- m$rep$PredSrvIdx[1, 1, , 1, 1]
    expect_equal(m$rep$SrvIdx_pop_nLL[1, 1, , 1, 1],
                 -stats::dnorm(log(obs_pop + addto), log(pred + addto), se_pop, log = TRUE),
                 tolerance = 1e-12)
  })

  test_that("a normal fleet evaluates its population data source on the arithmetic scale", {
    m <- report_of(build_input("normal"))
    pred <- m$rep$PredSrvIdx[1, 1, , 1, 1]
    expect_equal(m$rep$SrvIdx_pop_nLL[1, 1, , 1, 1],
                 -stats::dnorm(obs_pop, pred, se_pop, log = TRUE),
                 tolerance = 1e-12)
    expect_equal(m$rep$SrvIdx_nLL[1, , 1, 1],
                 -stats::dnorm(obs_reg, pred, se_reg, log = TRUE),
                 tolerance = 1e-12)
  })

  test_that("an mvn fleet keeps its population data source lognormal and its regional series joint", {
    n_obs <- sum(sim_data$UseSrvIdx == 1)
    lambda <- rep(0.8, n_obs); d <- 0.2 * obs_reg
    R <- outer(lambda, lambda); diag(R) <- 1
    S <- outer(d, d) * R

    m <- report_of(build_input("mvn", SrvIdx_Cov = list(S)))
    addto <- m$data$addtosrvidx
    pred <- m$rep$PredSrvIdx[1, 1, , 1, 1]

    expect_equal(m$rep$SrvIdx_pop_nLL[1, 1, , 1, 1],
                 -stats::dnorm(log(obs_pop + addto), log(pred + addto), se_pop, log = TRUE),
                 tolerance = 1e-12)

    resid <- obs_reg - pred
    full <- -0.5 * (n_obs * log(2 * pi) + determinant(S, logarithm = TRUE)$modulus[1] +
                      sum(resid * solve(S, resid)))
    expect_equal(sum(m$rep$SrvIdx_nLL[1, , 1, 1]), -full, tolerance = 1e-8)
  })

})
