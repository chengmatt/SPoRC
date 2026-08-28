library(SPoRC)
library(testthat)
library(Matrix)

# Simulate-refit recovery of the CTMC diffusion parameter under every move_timing.
#
# test-transition-operators.R and test-move-timing-selftest.R verify the movement
# operators internally, and test-move-timing-om-em-parity.R verifies that the simulator
# and the estimation model implement the same dynamics. None of them asks the estimation
# question: does the likelihood point back at the diffusion parameter that generated the
# data? That is what this file guards.
#
# The operating model is built so that every quantity the estimation model needs is known
# exactly -- deterministic recruitment, constant region-specific F, logistic selectivity in
# closed form, q = 1, M fixed. The estimation model is then pinned at those generating
# values and log_move_diffusion_pars is the single free parameter, started away from the
# truth. That isolates movement: a failure here is a movement-likelihood defect and not a
# 200-parameter optimization that wandered off.
#

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
  matrix(c("survey", "1"), nrow = nrow(tag_release_indicator()), ncol = 2, byrow = TRUE,
         dimnames = list(NULL, c("platform", "fleet")))
}

# ---------------------------------------------------------------------------
# Operating model
# ---------------------------------------------------------------------------

build_om <- function(move_timing, seed = 1234) {
  ages <- seq_len(N_AGES)

  sim_list <- Setup_Sim_Dim(n_sims = 1, n_yrs = N_YRS, n_regions = N_REGIONS, n_ages = N_AGES,
                           n_lens = NULL, n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = 1,
                           n_pop = 1)
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
    sim_list = sim_list, use_conv_fish_tagging = 1, n_tags = 2000,
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
    ctmc_move_dat = expand.grid(pop = 1, regions = seq_len(N_REGIONS), years = seq_len(N_YRS),
                                seas = 1, ages = seq_len(N_AGES), sexes = 1),
    preference_formula = ~ 0, diffusion_formula = ~ 1,
    log_move_diffusion_pars = TRUE_LOG_THETA, move_preference_pars = 0,
    area_r = rep(1, N_REGIONS), adjacency_mat = ADJ, ctmc_diffusion_bounds = 0,
    seasdur = 1, ctmc_scale_by_seasdur = 1
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

build_em <- function(om, move_timing) {
  sim_data <- simulation_data_to_SPoRC(sim_env = om, y = om$n_years, sim = 1)

  input_list <- Setup_Mod_Dim(
    years = seq_len(om$n_years), ages = seq_len(om$n_ages), lens = om$n_lens,
    n_regions = om$n_regions, n_sexes = om$n_sexes, n_fish_fleets = om$n_fish_fleets,
    n_srv_fleets = om$n_srv_fleets, n_pop = om$n_pop, natal_region = om$natal_region,
    verbose = FALSE
  )

  input_list <- Setup_Mod_Rec(
    input_list = input_list, do_rec_bias_ramp = 0, sigmaR_switch = 1,
    ln_sigmaR = array(log(1), c(2, 1, N_REGIONS)), rec_model = "mean_rec",
    use_rinit = 0, sigmaR_spec = "fix",
    init_age_strc = 2, equil_init_age_strc = 2,
    ln_global_R0 = log(sum(OM_R0))
  )

  input_list <- Setup_Mod_Biologicals(
    input_list = input_list,
    WAA = sim_data$WAA, MatAA = sim_data$MatAA,
    WAA_fish = sim_data$WAA_fish, WAA_srv = sim_data$WAA_srv,
    fit_lengths = 0, AgeingError = sim_data$AgeingError, M_spec = "fix",
    Fixed_natmort = array(0.3, dim = c(1, N_REGIONS, N_YRS, N_AGES, 1))
  )

  input_list <- suppressWarnings(Setup_Mod_Tagging(
    input_list = input_list, use_conv_fish_tagging = 1,
    conv_tag_release_indicator = sim_data$conv_tag_release_indicator,
    conv_tag_max_liberty = N_AGES / 2,
    conv_tagged_fish = sim_data$conv_tagged_fish,
    obs_conv_tag_fish_recap = sim_data$obs_conv_tag_fish_recap,
    conv_fish_tag_like = "Poisson", conv_tag_t_tagging = 1,
    conv_tagrep_spec = "fix", init_conv_tag_mort_spec = "fix", conv_tag_shed_spec = "fix",
    conv_fish_tag_attr = "p_a_s", conv_tag_release_platform = tag_release_platform()
  ))

  input_list <- Setup_Mod_Movement(
    input_list = input_list, move_type = 1, do_recruits_move = 0, use_fixed_movement = 0,
    ctmc_move_dat = expand.grid(pop = 1, regions = seq_len(N_REGIONS),
                                years = seq_len(N_YRS), seas = 1,
                                ages = seq_len(N_AGES), sexes = 1),
    adjacency_mat = ADJ, area_r = rep(1, N_REGIONS),
    diffusion_formula = ~ 1, preference_formula = ~ 0,
    move_timing = move_timing, log_move_diffusion_pars = START_LOG_THETA
  )

  input_list <- suppressWarnings(Setup_Mod_Catch_and_F(
    input_list = input_list, ObsCatch = sim_data$ObsCatch, UseCatch = sim_data$UseCatch,
    Use_F_pen = 0, sigmaC_spec = "fix", ln_sigmaC = sim_data$ln_sigmaC,
    ln_sigmaF = array(log(1), dim = c(N_REGIONS, 1, 1))
  ))

  # Region-split compositions. Aggregated ("agg") comps take the observed composition from
  # the first region only while comparing it to an expectation summed over all regions, so
  # they are not the right data spec for a multi-region operating model.
  input_list <- Setup_Mod_FishIdx_and_Comps(
    input_list = input_list,
    ObsFishIdx = sim_data$ObsFishIdx, ObsFishIdx_SE = sim_data$ObsFishIdx_SE,
    UseFishIdx = sim_data$UseFishIdx,
    ObsFishAgeComps = sim_data$ObsFishAgeComps, ObsFishLenComps = sim_data$ObsFishLenComps,
    UseFishAgeComps = sim_data$UseFishAgeComps, UseFishLenComps = sim_data$UseFishLenComps,
    ISS_FishAgeComps = sim_data$ISS_FishAgeComps, ISS_FishLenComps = sim_data$ISS_FishLenComps,
    fish_idx_type = "biom",
    FishAgeComps_LikeType = "Multinomial", FishLenComps_LikeType = "none",
    FishAgeComps_Type = "spltRjntS_Year_1-terminal_Fleet_1",
    FishLenComps_Type = "none_Year_1-terminal_Fleet_1"
  )

  input_list <- Setup_Mod_SrvIdx_and_Comps(
    input_list = input_list,
    ObsSrvIdx = sim_data$ObsSrvIdx, ObsSrvIdx_SE = sim_data$ObsSrvIdx_SE,
    UseSrvIdx = sim_data$UseSrvIdx,
    ObsSrvAgeComps = sim_data$ObsSrvAgeComps, ObsSrvLenComps = sim_data$ObsSrvLenComps,
    UseSrvAgeComps = sim_data$UseSrvAgeComps, UseSrvLenComps = sim_data$UseSrvLenComps,
    ISS_SrvAgeComps = sim_data$ISS_SrvAgeComps, ISS_SrvLenComps = sim_data$ISS_SrvLenComps,
    srv_idx_type = "biom",
    SrvAgeComps_LikeType = "Multinomial", SrvLenComps_LikeType = "none",
    SrvAgeComps_Type = "spltRjntS_Year_1-terminal_Fleet_1",
    SrvLenComps_Type = "none_Year_1-terminal_Fleet_1"
  )

  input_list <- Setup_Mod_Fishsel_and_Q(
    input_list = input_list, fish_sel_model = "logist1_Fleet_1",
    fish_fixed_sel_pars_spec = "est_shared_r", fish_q_spec = "est_shared_r"
  )

  input_list <- Setup_Mod_Srvsel_and_Q(
    input_list = input_list, srv_sel_model = "logist1_Fleet_1",
    srv_fixed_sel_pars_spec = "est_shared_r", srv_q_spec = "est_shared_r"
  )

  Setup_Mod_Weighting(
    input_list = input_list,
    Wt_Catch = 1, Wt_FishIdx = 1, Wt_SrvIdx = 1, Wt_Rec = 1, Wt_F = 1, Wt_Tagging = 1,
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
  for (nm in names(par)) {
    if (nm == "log_move_diffusion_pars") next
    map[[nm]] <- factor(rep(NA, length(as.vector(unlist(par[[nm]])))))
  }
  map$log_move_diffusion_pars <- factor(seq_along(par$log_move_diffusion_pars))

  list(data = input_list$data, par = par, map = map)
}

# One operating model per timing, reused by all three tests below
fixtures <- lapply(c(0, 1, 2), function(mt) {
  om <- build_om(mt)
  list(move_timing = mt, om = om, em = build_em(om, mt))
})
names(fixtures) <- paste0("move_timing_", c(0, 1, 2))

# ---------------------------------------------------------------------------
# 1. Precondition: the estimation model reproduces the operating model at the truth
# ---------------------------------------------------------------------------

test_that("estimation model reproduces the operating model at the generating parameters", {
  # Without this, recovering the diffusion parameter would be meaningless: the likelihood
  # would be reconciling two different population models rather than estimating movement.
  # This is also the check that catches a mis-laid-out age array in the operating model --
  # scrambling selectivity across ages shows up here as a ~15% SSB discrepancy.
  for (fx in fixtures) {
    pinned <- pin_at_truth(fx$em, TRUE_LOG_THETA)
    obj <- fit_model(pinned$data, pinned$par, pinned$map, random = NULL,
                     silent = TRUE, do_optim = FALSE)
    ssb_om <- as.vector(fx$om$SSB[, , , 1])
    ssb_em <- as.vector(obj$report(obj$par)$SSB)

    expect_equal(ssb_em, ssb_om, tolerance = 1e-3,
                 info = paste("move_timing =", fx$move_timing))
  }
})

# ---------------------------------------------------------------------------
# 2. Recovery of the diffusion parameter
# ---------------------------------------------------------------------------

test_that("CTMC diffusion parameter is recovered at every movement timing", {
  for (fx in fixtures) {
    pinned <- pin_at_truth(fx$em, START_LOG_THETA)
    fit <- fit_model(pinned$data, pinned$par, pinned$map, random = NULL, silent = TRUE)
    best <- fit$env$last.par.best

    # started at log(0.10), so this is a real recovery rather than a fixed point
    expect_equal(as.numeric(best[1]), TRUE_LOG_THETA, tolerance = 0.02,
                 info = paste("move_timing =", fx$move_timing))
    expect_lt(max(abs(fit$gr(best))), 1e-3)
    expect_jnLL_decomposes(fit, label = paste("move_timing =", fx$move_timing))
  }
})

# ---------------------------------------------------------------------------
# 3. The likelihood surface itself is minimized at the truth
# ---------------------------------------------------------------------------

test_that("likelihood over the diffusion parameter is minimized at the generating value", {
  # Guards against the optimizer happening to land on the right answer while the surface
  # is flat or its minimum sits somewhere else.
  grid <- log(c(0.15, 0.20, 0.25, 0.30, 0.35, 0.45, 0.60))
  for (fx in fixtures) {
    pinned <- pin_at_truth(fx$em, TRUE_LOG_THETA)
    obj <- fit_model(pinned$data, pinned$par, pinned$map, random = NULL,
                     silent = TRUE, do_optim = FALSE)
    nll <- vapply(grid, function(lt) obj$fn(lt), numeric(1))

    expect_equal(grid[which.min(nll)], TRUE_LOG_THETA, tolerance = 1e-8,
                 info = paste("move_timing =", fx$move_timing))
  }
})
