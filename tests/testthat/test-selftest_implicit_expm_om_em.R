library(SPoRC)
library(testthat)
library(Matrix)

# Operating model on the exact matrix exponential, estimation model on the implicit solve.
#
# move_expm_nsub replaces expm(A) with (I - A/n)^-n wherever SPoRC exponentiates a CTMC
# generator. That is a first-order discretization, not a cheaper route to the same
# numbers, so the question this file answers is: if the world runs on the exponential and
# the assessment does not, what does the assessment get wrong, and how fast does that go
# away as substeps are added?
#
# The operating model is the one from test-integration_move_timing_diffusion.R: everything
# the estimation model needs is known exactly (deterministic recruitment, constant
# region-specific F, closed-form selectivity, q = 1, M fixed), the generator is handed to
# the simulator directly as Mrate, and log_move_diffusion_pars is the single free
# parameter. Isolating movement that way is what makes the bias attributable to the
# exponential rather than to a 200-parameter fit wandering off.
#
# The operating model always uses expm_nsub = 0. Only the estimation model varies.

N_YRS <- 25
N_REGIONS <- 3
N_AGES <- 10
TRUE_LOG_THETA <- log(0.30)   # generating CTMC diffusion parameter
START_LOG_THETA <- log(0.10)  # deliberately off the truth, so recovery is not a no-op
OM_R0 <- rep(10, N_REGIONS)
OM_F <- c(0.05, 0.15, 0.30)   # region-varying F is what makes diffusion identifiable
FISH_B50 <- 5; FISH_K <- 3
SRV_B50 <- 3;  SRV_K <- 1
ADJ <- matrix(1L, N_REGIONS, N_REGIONS); diag(ADJ) <- 0L

# Lay a vector of age-specific values into a (pop, region, year, seas, age, ...) array.
# The age dimension is the 5th, so the values have to repeat over the product of every
# earlier dimension. The single-region idiom `rep(v, each = n_yrs)` happens to be correct
# only when n_pop = n_regions = n_seas = 1; with more than one region it silently permutes
# the ages, which shows up as an operating model whose selectivity does not match the
# logistic curve the estimation model estimates.
age_array <- function(v, dims) array(rep(v, each = prod(dims[1:4])), dim = dims)

tag_release_indicator <- function() {
  expand.grid(regions = seq_len(N_REGIONS), tag_years = seq_len(N_YRS), tag_seas = 1)
}

tag_release_platform <- function() {
  matrix(
    c("survey", "1"),
    nrow = nrow(tag_release_indicator()),
    ncol = 2,
    byrow = TRUE,
    dimnames = list(NULL, c("platform", "fleet"))
  )
}

# ---------------------------------------------------------------------------
# Operating model
# ---------------------------------------------------------------------------

build_om <- function(move_timing, seed = 1234) {
  ages <- seq_len(N_AGES)

  sim_list <- Setup_Sim_Dim(
    n_sims = 1,
    n_yrs = N_YRS,
    n_regions = N_REGIONS,
    n_ages = N_AGES,
    n_lens = NULL,
    n_sexes = 1,
    n_fish_fleets = 1,
    n_srv_fleets = 1,
    n_pop = 1
  )
  sim_list <- Setup_Sim_Containers(sim_list)

  sim_list <- Setup_Sim_Fishing(
    sim_list = sim_list,
    Fmort_input = array(rep(OM_F, times = N_YRS), dim = c(N_REGIONS, N_YRS, 1, 1, 1)),
    fish_sel_input = replicate(1, age_array(1 / (1 + exp(-FISH_K * (ages - FISH_B50))),
                                            c(1, N_REGIONS, N_YRS, 1, N_AGES, 1, 1)))
  )

  sim_list <- Setup_Sim_Survey(
    sim_list = sim_list,
    srv_sel_input = replicate(1, age_array(1 / (1 + exp(-SRV_K * (ages - SRV_B50))),
                                           c(1, N_REGIONS, N_YRS, 1, N_AGES, 1, 1)))
  )

  waa <- 5 / (1 + exp(-3 * (ages - 3)))
  mat <- 1 / (1 + exp(-3 * (ages - 3)))
  sim_list <- suppressWarnings(Setup_Sim_Biologicals(
    sim_list = sim_list,
    natmort_input = replicate(1, array(0.3, dim = c(1, N_REGIONS, N_YRS, N_AGES, 1))),
    WAA_input      = replicate(1, age_array(waa, c(1, N_REGIONS, N_YRS, 1, N_AGES, 1))),
    WAA_fish_input = replicate(1, age_array(waa, c(1, N_REGIONS, N_YRS, 1, N_AGES, 1, 1))),
    WAA_srv_input  = replicate(1, age_array(waa, c(1, N_REGIONS, N_YRS, 1, N_AGES, 1, 1))),
    MatAA_input    = replicate(1, age_array(mat, c(1, N_REGIONS, N_YRS, 1, N_AGES, 1)))
  ))

  # Tag data is the observation type that normally identifies movement, so the operating
  # model releases tags in every region and year.
  sim_list <- Setup_Sim_Tagging(
    sim_list = sim_list,
    use_conv_fish_tagging = 1,
    n_tags = 2000,
    conv_tag_max_liberty = N_AGES / 2,
    conv_tag_release_indicator = tag_release_indicator(),
    conv_tag_release_platform = tag_release_platform(),
    conv_tag_t_tagging = 1
  )

  # CTMC movement at a known diffusion parameter, built through Get_Movement so the
  # operating model uses exactly the parameterization the estimation model inverts.
  mv <- Get_Movement(
    move_type = 1, do_recruits_move = 0,
    n_pop = 1, n_regions = N_REGIONS, n_yrs = N_YRS, n_proj_yrs_devs = 0,
    n_ages = N_AGES, n_sexes = 1, n_seas = 1,
    move_pars = array(0, c(1, N_REGIONS, N_REGIONS - 1, N_YRS, 1, N_AGES, 1)),
    move_devs = array(0, c(1, N_REGIONS, N_REGIONS - 1, N_YRS, 1, N_AGES, 1)),
    use_fixed_movement = 0, Fixed_Movement = NULL,
    ctmc_move_dat = expand.grid(
      pop = 1,
      regions = seq_len(N_REGIONS),
      years = seq_len(N_YRS),
      seas = 1,
      ages = seq_len(N_AGES),
      sexes = 1
    ),
    preference_formula = ~ 0, diffusion_formula = ~ 1,
    log_move_diffusion_pars = TRUE_LOG_THETA, move_preference_pars = 0,
    area_r = rep(1, N_REGIONS), adjacency_mat = ADJ, ctmc_diffusion_bounds = 0,
    seasdur = 1, ctmc_scale_by_seasdur = 1, expm_nsub = 0  # OM always uses the exact exponential
  )

  # Recruits do not move, so Get_Movement leaves age 1 at zero. The dynamics substitute an
  # identity transition there, but the simulation arrays still have to be well formed.
  Movement <- mv$Movement
  Mrate <- mv$Mrate
  for (y in seq_len(N_YRS)) {
    Movement[1, , , y, 1, 1, 1] <- diag(N_REGIONS)
    Mrate[1, , , y, 1, 1, 1] <- 0
  }

  sim_list <- Setup_Sim_Rec(
    sim_list = sim_list,
    R0_input = replicate(1, array(rep(OM_R0, times = N_YRS), dim = c(1, N_REGIONS, N_YRS))),
    use_rinit = 0,
    ln_sigmaR = array(log(1e-4), dim = c(2, 1, N_REGIONS)),  # deterministic recruitment
    recruitment_opt = "mean_rec",
    init_age_strc = 2   # iterative equilibrium: the general path, valid with movement
  )

  # Setup_Sim_Rec resets the movement fields, so assert them afterwards
  sim_list$Movement <- replicate(1, Movement)
  sim_list$Mrate <- replicate(1, Mrate)
  sim_list$move_timing <- move_timing
  sim_list$expm_nsub <- 0  # and so do the operating model's own seasonal operators

  # Near-noiseless observations, so any recovery error is attributable to the movement
  # likelihood rather than to observation error
  sim_list$ISS_FishAgeComps[] <- 1e5
  sim_list$ISS_SrvAgeComps[] <- 1e5
  sim_list$ObsFishIdx_SE[] <- 0.05
  sim_list$ObsSrvIdx_SE[] <- 0.05

  set.seed(seed)
  Simulate_Pop_Static(sim_list = sim_list, output_path = NULL)
}

# ---------------------------------------------------------------------------
# Estimation model, matching the operating model's structure
# ---------------------------------------------------------------------------

build_em <- function(om, move_timing, em_expm_nsub = 0) {
  sim_data <- simulation_data_to_SPoRC(sim_env = om, y = om$n_years, sim = 1)

  input_list <- Setup_Mod_Dim(
    years = seq_len(om$n_years),
    ages = seq_len(om$n_ages),
    lens = om$n_lens,
    n_regions = om$n_regions,
    n_sexes = om$n_sexes,
    n_fish_fleets = om$n_fish_fleets,
    n_srv_fleets = om$n_srv_fleets,
    n_pop = om$n_pop,
    natal_region = om$natal_region,
    verbose = FALSE
  )

  input_list <- Setup_Mod_Rec(
    input_list = input_list,
    do_rec_bias_ramp = 0,
    sigmaR_switch = 1,
    ln_sigmaR = array(log(1), c(2, 1, N_REGIONS)),
    rec_model = "mean_rec",
    use_rinit = 0,
    sigmaR_spec = "fix",
    init_age_strc = 2,
    equil_init_age_strc = 2,
    ln_global_R0 = log(sum(OM_R0))
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
    Fixed_natmort = array(0.3, dim = c(1, N_REGIONS, N_YRS, N_AGES, 1))
  )

  input_list <- suppressWarnings(Setup_Mod_Tagging(
    input_list = input_list,
    use_conv_fish_tagging = 1,
    conv_tag_release_indicator = sim_data$conv_tag_release_indicator,
    conv_tag_max_liberty = N_AGES / 2,
    conv_tagged_fish = sim_data$conv_tagged_fish,
    obs_conv_tag_fish_recap = sim_data$obs_conv_tag_fish_recap,
    conv_fish_tag_like = "Poisson",
    conv_tag_t_tagging = 1,
    conv_tagrep_spec = "fix",
    init_conv_tag_mort_spec = "fix",
    conv_tag_shed_spec = "fix",
    conv_fish_tag_attr = "p_a_s",
    conv_tag_release_platform = tag_release_platform()
  ))

  input_list <- Setup_Mod_Movement(
    input_list = input_list,
    move_type = 1,
    do_recruits_move = 0,
    use_fixed_movement = 0,
    ctmc_move_dat = expand.grid(
      pop = 1,
      regions = seq_len(N_REGIONS),
      years = seq_len(N_YRS),
      seas = 1,
      ages = seq_len(N_AGES),
      sexes = 1
    ),
    adjacency_mat = ADJ,
    area_r = rep(1, N_REGIONS),
    diffusion_formula = ~ 1,
    preference_formula = ~ 0,
    move_timing = move_timing,
    log_move_diffusion_pars = START_LOG_THETA,
    move_expm_nsub = em_expm_nsub
  )

  input_list <- suppressWarnings(Setup_Mod_Catch_and_F(
    input_list = input_list,
    ObsCatch = sim_data$ObsCatch,
    UseCatch = sim_data$UseCatch,
    Use_F_pen = 0,
    sigmaC_spec = "fix",
    ln_sigmaC = sim_data$ln_sigmaC,
    ln_sigmaF = array(log(1), dim = c(N_REGIONS, 1, 1))
  ))

  # Region-split compositions. Aggregated ("agg") comps take the observed composition from
  # the first region only while comparing it to an expectation summed over all regions, so
  # they are not the right data spec for a multi-region operating model.
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
    fish_idx_type = "biom",
    FishAgeComps_LikeType = "Multinomial",
    FishLenComps_LikeType = "none",
    FishAgeComps_Type = "spltRjntS_Year_1-terminal_Fleet_1",
    FishLenComps_Type = "none_Year_1-terminal_Fleet_1"
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
    srv_idx_type = "biom",
    SrvAgeComps_LikeType = "Multinomial",
    SrvLenComps_LikeType = "none",
    SrvAgeComps_Type = "spltRjntS_Year_1-terminal_Fleet_1",
    SrvLenComps_Type = "none_Year_1-terminal_Fleet_1"
  )

  input_list <- Setup_Mod_Fishsel_and_Q(
    input_list = input_list,
    fish_sel_model = "logist1_Fleet_1",
    fish_fixed_sel_pars_spec = "est_shared_r",
    fish_q_spec = "est_shared_r"
  )

  input_list <- Setup_Mod_Srvsel_and_Q(
    input_list = input_list,
    srv_sel_model = "logist1_Fleet_1",
    srv_fixed_sel_pars_spec = "est_shared_r",
    srv_q_spec = "est_shared_r"
  )

  Setup_Mod_Weighting(
    input_list = input_list,
    Wt_Catch = 1,
    Wt_FishIdx = 1,
    Wt_SrvIdx = 1,
    Wt_Rec = 1,
    Wt_F = 1,
    Wt_Tagging = 1,
    Wt_FishAgeComps = array(1, dim = c(N_REGIONS, N_YRS, 1, 1, 1)),
    Wt_FishLenComps = array(0, dim = c(N_REGIONS, N_YRS, 1, 1, 1)),
    Wt_SrvAgeComps  = array(1, dim = c(N_REGIONS, N_YRS, 1, 1, 1)),
    Wt_SrvLenComps  = array(0, dim = c(N_REGIONS, N_YRS, 1, 1, 1))
  )
}

# Pin every parameter at the value the operating model generated with, and free only the
# CTMC diffusion parameter. logist1 selectivity is b50 = exp(par1), k = exp(par2); fishing
# mortality is exp(ln_F_mean + ln_F_devs); the operating model used q = 1.
pin_at_truth <- function(input_list, log_theta) {
  par <- input_list$par
  par$ln_global_R0[] <- log(sum(OM_R0))
  if (length(par$rec_region_prop_pars)) par$rec_region_prop_pars[] <- log(OM_R0[-1] / OM_R0[1])
  par$ln_RecDevs[] <- 0
  par$ln_InitDevs[] <- 0
  par$ln_F_mean[] <- rep(log(OM_F), times = prod(dim(par$ln_F_mean)) / N_REGIONS)
  par$ln_F_devs[] <- 0
  par$fish_fixed_sel_pars[, 1, , , ] <- log(FISH_B50)
  par$fish_fixed_sel_pars[, 2, , , ] <- log(FISH_K)
  par$srv_fixed_sel_pars[, 1, , , ] <- log(SRV_B50)
  par$srv_fixed_sel_pars[, 2, , , ] <- log(SRV_K)
  par$ln_fish_q[] <- 0
  par$ln_srv_q[] <- 0
  par$log_move_diffusion_pars[] <- log_theta

  map <- input_list$map
  for (quant_name in names(par)) {
    if (quant_name == "log_move_diffusion_pars") next
    map[[quant_name]] <- factor(rep(NA, length(as.vector(unlist(par[[quant_name]])))))
  }
  map$log_move_diffusion_pars <- factor(seq_along(par$log_move_diffusion_pars))

  list(data = input_list$data, par = par, map = map)
}


# ---------------------------------------------------------------------------
# Test setups. The operating models are built once; only the estimation model varies.
# ---------------------------------------------------------------------------

oms <- list("2" = build_om(2), "0" = build_om(0))

# Report SSB at the generating parameters, with no estimation involved. This isolates the
# forward distortion: how far the implicit operator moves the population itself.
ssb_rel_err <- function(om, move_timing, nsub) {
  em <- build_em(om, move_timing, em_expm_nsub = nsub)
  pinned <- pin_at_truth(em, TRUE_LOG_THETA)
  obj <- fit_model(
    pinned$data,
    pinned$par,
    pinned$map,
    random = NULL,
    silent = TRUE,
    do_optim = FALSE
  )
  truth <- as.vector(om$SSB[, , , 1])
  max(abs(as.vector(obj$report(obj$par)$SSB) - truth) / truth)
}

# ---------------------------------------------------------------------------
# 1. Forward distortion of the population, at the generating parameters
# ---------------------------------------------------------------------------

test_that("the implicit solve leaves SSB alone at nsub = 0 and converges back to it", {
  # nsub = 0 has to be the exact exponential, so this is the same agreement the
  # move_timing test asserts. nsub = 1 is plain solve(I - A) and is expected to be badly
  # off -- pinned here so that a future change cannot quietly make it look harmless.
  err <- vapply(c(0, 1, 8, 512), function(n) ssb_rel_err(oms[["2"]], 2, n), numeric(1))

  expect_lt(err[1], 1e-3)                       # exact
  expect_gt(err[2], 0.05)                       # one backward Euler step is a different model
  expect_true(all(diff(err[-1]) < 0))           # more substeps, less distortion
  expect_lt(err[4], 5e-3)                       # and it converges back
})

test_that("under move_timing 0 only the movement fractions are approximated", {
  # Timings 0 and 1 never exponentiate the generator net of mortality: survival stays
  # elementwise exp(-Z) and only Get_Movement's fractions go through mat_exp. The forward
  # error should therefore be far smaller than the same nsub under continuous movement,
  # where 1/(1 + Z) replaces exp(-Z) and compounds over seasons and ages.
  err_mt0 <- ssb_rel_err(oms[["0"]], 0, 1)
  err_mt2 <- ssb_rel_err(oms[["2"]], 2, 1)

  expect_lt(err_mt0, 0.05)
  expect_gt(err_mt2, 5 * err_mt0)
})

# ---------------------------------------------------------------------------
# 2. Bias in the estimated movement rate
# ---------------------------------------------------------------------------

test_that("an implicit estimation model inflates diffusion, and substeps remove the bias", {
  # The estimation model cannot represent the operating model's movement, so it trades:
  # backward Euler smears abundance further per step than the exponential does at the same
  # rate, and the likelihood answers by pulling the rate up. That direction is the
  # informative part -- a movement rate biased high is not a conservative error.
  fit_theta <- function(nsub) {
    em <- build_em(oms[["2"]], 2, em_expm_nsub = nsub)
    pinned <- pin_at_truth(em, START_LOG_THETA)
    fit <- fit_model(pinned$data, pinned$par, pinned$map, random = NULL, silent = TRUE)
    best <- fit$env$last.par.best
    expect_lt(max(abs(fit$gr(best))), 1e-3)
    exp(as.numeric(best[1]))
  }

  truth <- exp(TRUE_LOG_THETA)
  theta_exact <- fit_theta(0)
  theta_be1 <- fit_theta(1)
  theta_be64 <- fit_theta(64)

  # started at log(0.10), so recovery at nsub = 0 is a real recovery
  expect_equal(theta_exact, truth, tolerance = 0.02)

  # one backward Euler step costs about 12% on the movement rate, biased high
  expect_gt(theta_be1 / truth, 1.05)

  # 64 substeps put the estimate back within the exact model's own recovery error
  expect_lt(abs(theta_be64 - truth), 3 * abs(theta_exact - truth) + 0.005)
  expect_lt(abs(theta_be64 - truth), abs(theta_be1 - truth) / 5)
})
