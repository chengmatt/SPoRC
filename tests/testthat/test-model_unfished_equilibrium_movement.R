# An unfished three-region run, regions differing in weight and maturity, stays at equilibrium in every year: SSB
# and Dynamic_SSB0 equal S0 from the per recruit table and recruitment equals R0 (1e-10), and the self test's OM agrees.

library(SPoRC)
library(testthat)

quiet <- function(expr) suppressWarnings(suppressMessages(expr))

# Setup ----------------------------------------------------------------------

# Beverton-Holt by region with no catch, started at the equilibrium with movement, and continuous time movement so
# that move_timing 2 has a generator
unfished_input <- function(move_timing, move_expm_nsub = 0, do_recruits_move = 0, rec_lag = 1, n_seas = 1) {

  n_regions <- 3
  n_yrs <- 13
  n_ages <- 7
  biol_dim <- c(1, n_regions, n_yrs, n_seas, n_ages, 2)

  # region 2 at 60 pct of the weight and maturing a year later, region 3 at 130 pct and a year earlier
  wt_scale <- c(1, 0.6, 1.3)
  mat_shift <- c(0, 1, -1)
  WAA <- MatAA <- array(0, dim = biol_dim)
  for(r in 1:n_regions) {
    for(a in 1:n_ages) {
      WAA[,r,,,a,] <- wt_scale[r] * 5 / (1 + exp(-(a - 3))) # weight at age (kg)
      MatAA[,r,,,a,] <- 1 / (1 + exp(-(a - 3 - mat_shift[r]))) # maturity at age
    } # end a loop
  } # end r loop
  if(rec_lag == 0) MatAA[,,,,1,] <- 0 # age 0 recruits cannot be mature

  il <- sweep_input(
    dims = list(n_regions = n_regions, n_seas = n_seas, n_fish_fleets = 1),
    rec = list(
      rec_model = "bh_rec",
      rec_dd = "local", # one curve per region, so every region's S0 is checked
      h_spec = "fix", # steepness kept at its start
      init_age_strc = 2, # matrix geometric series, the equilibrium with movement at every age
      rec_lag = rec_lag,
      spawn_seas = n_seas, # spawning in the last season
      t_spawn = 0.4 # part way into the season, so the reads inside a season are used
    ),
    biol = list(
      WAA = WAA,
      MatAA = MatAA,
      WAA_fish = array(WAA, dim = c(biol_dim, 1)),
      WAA_srv = array(WAA, dim = c(biol_dim, 1)),
      Fixed_natmort = array(0.2, dim = biol_dim)
    ),
    move = list(
      move_type = 1, # continuous time movement
      use_fixed_movement = 0,
      Fixed_Movement = NA,
      do_recruits_move = do_recruits_move,
      ctmc_move_dat = expand.grid(pop = 1, regions = 1:n_regions, years = 1:n_yrs, seas = 1:n_seas, ages = 1:n_ages, sexes = 1:2),
      adjacency_mat = 1 - diag(n_regions), # every region neighbors every other
      area_r = rep(1, n_regions),
      diffusion_formula = ~ 1,
      preference_formula = ~ 0,
      move_timing = move_timing,
      move_expm_nsub = move_expm_nsub,
      log_move_diffusion_pars = log(0.4)
    ),
    catch = list(UseCatch = array(0, dim = c(n_regions, n_yrs, n_seas, 1))), # no catch fit, so F is zero
    srvidx = list(srv_idx_type = "biom")
  )

  # recruits split about 0.59, 0.29 and 0.12 across regions, since even shares under symmetric movement never move
  il$par$rec_region_prop_pars[] <- c(-0.7, -1.6)
  il
}

# unfished spawning biomass by region, from the per recruit table built off what the model reports
unfished_S0 <- function(rep, data) {

  n_regions <- data$n_regions
  n_seas <- data$n_seas
  n_ages <- length(data$ages)
  SR <- data$SR_ref_yr

  sbpr <- SPoRC:::Get_SBPR_Table(
    rec_dd = data$rec_dd,
    rec_region_prop = rep$rec_region_prop,
    rec_seas_prop = rep$rec_seas_prop,
    n_pop = 1,
    n_regions = n_regions,
    n_ages = n_ages,
    n_fish_fleets = 1,
    WAA = array(data$WAA[,,SR,,,1], dim = c(1, n_regions, n_seas, n_ages)),
    MatAA = array(data$MatAA[,,SR,,,1], dim = c(1, n_regions, n_seas, n_ages)),
    natmort = array(rep$natmort[,,SR,,,1], dim = c(1, n_regions, n_seas, n_ages)),
    Movement = array(rep$Movement[,,,SR,,,1], dim = c(1, n_regions, n_regions, n_seas, n_ages)),
    sgl_seas_spawning_movement = array(rep$sgl_seas_spawning_movement[,,,SR,,1], dim = c(1, n_regions, n_regions, n_ages)),
    do_recruits_move = data$do_recruits_move,
    t_spawn = data$t_spawn,
    init_F = array(0, dim = c(n_regions, n_seas, 1)), # unfished
    dmr = array(rep$dmr[,SR,,], dim = c(n_regions, n_seas, 1)),
    fish_sel = array(rep$fish_sel[,,SR,,,1,], dim = c(1, n_regions, n_seas, n_ages, 1)),
    ret_sel = array(rep$ret_sel[,,SR,,,1,], dim = c(1, n_regions, n_seas, n_ages, 1)),
    n_seas = n_seas,
    spawn_seas = data$spawn_seas,
    seasdur = data$seasdur,
    sexratio_f = array(rep$sexratio[,,SR,1], dim = c(1, n_regions)),
    Mrate = if(is.null(rep$Mrate)) NULL else array(rep$Mrate[,,,SR,,,1], dim = c(1, n_regions, n_regions, n_seas, n_ages)),
    move_timing = data$move_timing,
    expm_nsub = data$move_expm_nsub # the matrix exponential the dynamics use
  )

  sbpr$phi0[1,] * rep$R0[1]
}

# largest relative gap from equilibrium over regions and years
unfished_gaps <- function(il) {

  rep <- quiet(fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE))$rep
  yrs <- seq_along(il$data$years)
  S0 <- unfished_S0(rep, il$data)
  R0_r <- rep$R0[1] * rep$rec_region_prop[1,] # unfished recruitment by region

  # refuse a run that would make the comparisons below empty or trivial
  if(!identical(dim(rep$Dynamic_SSB0), c(1L, 3L, length(yrs)))) stop("Dynamic_SSB0 is not population by region by year")
  if(max(rep$Fmort) != 0) stop("the run is fished")
  if(diff(range(S0)) / mean(S0) < 0.1) stop("S0 hardly differs among regions, so nothing tests the regions apart")

  c(
    "SSB / S0" = max(abs(rep$SSB[1,,] / S0 - 1)),
    "Dynamic_SSB0 / S0" = max(abs(rep$Dynamic_SSB0[1,,] / S0 - 1)),
    "Rec / R0" = max(abs(rep$Rec[1,,] / R0_r - 1)),
    "NAA0 against NAA" = max(abs(rep$NAA0[,,yrs,,,] - rep$NAA[,,yrs,,,])) / max(rep$NAA[,,yrs,,,])
  )
}

# Estimation Model -----------------------------------------------------------

test_that("an unfished run under move_timing 0 stays at S0 in every region and year", {

  # unfished numbers move at the season start with the fished ones, so Dynamic_SSB0 is the unfished trajectory
  configs <- expand.grid(do_recruits_move = c(0, 1), rec_lag = c(1, 0), n_seas = 1)
  configs <- rbind(configs, data.frame(do_recruits_move = 0, rec_lag = 1, n_seas = 2))

  for(i in seq_len(nrow(configs))) {
    cfg <- configs[i,]
    gaps <- unfished_gaps(quiet(unfished_input(move_timing = 0, do_recruits_move = cfg$do_recruits_move,
                                               rec_lag = cfg$rec_lag, n_seas = cfg$n_seas)))
    for(gap_name in names(gaps)) {
      expect_lt(gaps[[gap_name]], 1e-10,
                label = sprintf("%s, recruits move %d, rec_lag %d, %d seasons", gap_name, cfg$do_recruits_move, cfg$rec_lag, cfg$n_seas))
    } # end gap_name loop
  } # end i loop
})

test_that("an unfished run under move_timing 2 with move_expm_nsub = 2 stays at S0 in every region and year", {

  # the per recruit table behind S0 takes the same implicit matrix exponential as the dynamics
  configs <- expand.grid(do_recruits_move = c(0, 1), rec_lag = c(1, 0))

  for(i in seq_len(nrow(configs))) {
    cfg <- configs[i,]
    gaps <- unfished_gaps(quiet(unfished_input(move_timing = 2, move_expm_nsub = 2, do_recruits_move = cfg$do_recruits_move,
                                               rec_lag = cfg$rec_lag)))
    for(gap_name in names(gaps)) {
      expect_lt(gaps[[gap_name]], 1e-10,
                label = sprintf("%s, recruits move %d, rec_lag %d", gap_name, cfg$do_recruits_move, cfg$rec_lag))
    } # end gap_name loop
  } # end i loop
})

# Operating Model ------------------------------------------------------------

test_that("the self test's operating model reproduces the unfished run, unfished numbers included", {

  # move_timing 0 checks the unfished numbers move on both sides; move_timing 2 checks both start from the same
  # initial age structure under the implicit matrix exponential
  configs <- data.frame(move_timing = c(0, 2), move_expm_nsub = c(0, 2))

  for(i in seq_len(nrow(configs))) {
    cfg <- configs[i,]
    run <- lockstep_om(quiet(unfished_input(move_timing = cfg$move_timing, move_expm_nsub = cfg$move_expm_nsub)))
    yrs <- seq_along(run$data$years)
    expect_equal(run$sim_list$expm_nsub, cfg$move_expm_nsub) # the OM runs the fit's matrix exponential

    gaps <- c(
      NAA = lockstep_rel(run$om$NAA[,,yrs,,,,1], run$rep$NAA[,,yrs,,,]),
      NAA0 = lockstep_rel(run$om$NAA0[,,yrs,,,,1], run$rep$NAA0[,,yrs,,,]),
      SSB = lockstep_rel(run$om$SSB[,,yrs,1], run$rep$SSB[,,yrs]),
      Dynamic_SSB0 = lockstep_rel(run$om$Dynamic_SSB0[,,yrs,1], run$rep$Dynamic_SSB0[,,yrs])
    )
    for(gap_name in names(gaps)) {
      expect_lt(gaps[[gap_name]], 1e-10, label = sprintf("%s, move_timing %d", gap_name, cfg$move_timing))
    } # end gap_name loop
  } # end i loop
})
