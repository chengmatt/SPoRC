library(SPoRC)
library(testthat)

test_that("Multi-region, population, and seasonal local BH MSY reference points converges to equilibrium", {

  set.seed(555)

  ### Setup Operating Model ---------------------------------------------------
  sim_list <- Setup_Sim_Dim(n_sims = 1, # number of simulations
                            n_yrs = 30, # number of years
                            n_regions = 2,  # number of regions
                            n_ages = 10, # number of ages
                            n_lens = NULL, # number of lengths
                            n_sexes = 1, # number of sexes
                            n_fish_fleets = 1, # number of fishery fleets
                            n_srv_fleets = 1,  # number of survey fleets
                            n_seas = 2,  # number of seasons
                            n_pop = 3, # number of populaitons
                            natal_region = c(1, 1, 2) # natal regions
  )

  ### Setup Simulation Containers ---------------------------------------------
  sim_list <- Setup_Sim_Containers(sim_list)

  ### Setup Fishing Processes -------------------------------------------------
  sim_list <- Setup_Sim_Fishing(
    sim_list = sim_list,

    # logistic selectivity
    fish_sel_input = replicate(
      n = sim_list$n_sims,
      {
        arr <- array(NA, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas,
                                 sim_list$n_ages, sim_list$n_sexes, sim_list$n_fish_fleets))
        for (r in 1:sim_list$n_regions)
          for (y in 1:sim_list$n_yrs)
            for (s in 1:sim_list$n_sexes) {
              for(p in 1:sim_list$n_pop) {
                for(seas in 1:sim_list$n_seas) {
                  arr[p,r, y,seas, , s, 1] <-  1 / (1 + exp(-1.5 * (1:sim_list$n_ages - 4)))
                }
              }
            }
        arr
      }
    ),

    Fmort_input = {
      n = sim_list$n_yrs * sim_list$n_seas * sim_list$n_sims * sim_list$n_fish_fleets
      t = seq(0, 2*pi, length.out = n)
      arr <- array(NA, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas,
                               sim_list$n_fish_fleets, sim_list$n_sims))
      arr[1,,,,] <- 0.15 * exp(sin(t) + rnorm(n, 0, 0.1))   # region 1 higher F, peaks early
      arr[2,,,,] <- 0.05 * exp(-sin(t) + rnorm(n, 0, 0.1))  # region 2 lower F, peaks late
      arr
    },

    # Fishery Age ISS
    ISS_FishAgeComps = array(500, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes,
                                          sim_list$n_fish_fleets, sim_list$n_sims)),
    ISS_FishAgeComps_pop = array(round(500 / sim_list$n_pop), dim = c(sim_list$n_pop, sim_list$n_regions,
                                                                      sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes,
                                                                      sim_list$n_fish_fleets, sim_list$n_sims)),


    # Sigma for catch
    ln_sigmaC = array(log(0.01), dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets)),
    ln_sigmaC_pop = array(log(0.01), dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets))

  )

  ### Setup Survey Process ----------------------------------------------------
  sim_list <- Setup_Sim_Survey(
    sim_list = sim_list,

    # Logistic selectivity
    srv_sel_input = replicate(
      n = sim_list$n_sims,
      {
        arr <- array(NA, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas,
                                 sim_list$n_ages, sim_list$n_sexes, sim_list$n_srv_fleets))
        for (r in 1:sim_list$n_regions)
          for (y in 1:sim_list$n_yrs)
            for (s in 1:sim_list$n_sexes) {
              for(p in 1:sim_list$n_pop) {
                for(seas in 1:sim_list$n_seas) {
                  arr[p,r, y,seas, , s, 1] <-  1 / (1 + exp(-1 * (1:sim_list$n_ages - 2.5)))

                }
              }
            }
        arr
      }
    ),

    # Survey Age ISS
    ISS_SrvAgeComps = array(500, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes,
                                         sim_list$n_srv_fleets, sim_list$n_sims)),
    ISS_SrvAgeComps_pop = array(round(500 / sim_list$n_pop), dim = c(sim_list$n_pop, sim_list$n_regions,
                                                                     sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes,
                                                                     sim_list$n_srv_fleets, sim_list$n_sims)),

    # Sigma for Survey Index
    ObsSrvIdx_SE = array(0.15, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_srv_fleets)),
    ObsSrvIdx_pop_SE = array(0.15, dim = c(sim_list$n_pop, sim_list$n_regions,
                                           sim_list$n_yrs, sim_list$n_seas, sim_list$n_srv_fleets))

  )

  ### Setup Biological Dynamics -----------------------------------------------
  suppressWarnings(
    sim_list <- Setup_Sim_Biologicals(
      sim_list = sim_list,

      # Natural Mortality
      natmort_input = array(0.3, dim = c(sim_list$n_pop, sim_list$n_regions,
                                         sim_list$n_yrs, sim_list$n_ages,
                                         sim_list$n_sexes, sim_list$n_sims)),

      # Weight at age - Same for all pops
      WAA_input = replicate(
        n = sim_list$n_sims,
        array(
          rep(5 / (1 + exp(-3 * ((1:sim_list$n_ages) - 3))),
              each = sim_list$n_pop * sim_list$n_regions * sim_list$n_yrs * sim_list$n_seas,
              times = sim_list$n_sexes),
          dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs,
                  sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes)
        )
      ),

      # Fishery weight at age - same as WAA_input
      WAA_fish_input = replicate(
        n = sim_list$n_sims,
        array(
          rep(5 / (1 + exp(-3 * ((1:sim_list$n_ages) - 3))),
              each = sim_list$n_pop * sim_list$n_regions * sim_list$n_yrs * sim_list$n_seas,
              times = sim_list$n_sexes * sim_list$n_fish_fleets),
          dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes, sim_list$n_fish_fleets)
        )
      ),

      # Survey weight at age - same as WAA_input
      WAA_srv_input = replicate(
        n = sim_list$n_sims,
        array(
          rep(5 / (1 + exp(-3 * ((1:sim_list$n_ages) - 3))),
              each = sim_list$n_pop * sim_list$n_regions * sim_list$n_yrs * sim_list$n_seas,
              times = sim_list$n_sexes * sim_list$n_srv_fleets),
          dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes, sim_list$n_srv_fleets)
        )
      ),

      # Maturity at age
      MatAA_input = replicate(
        n = sim_list$n_sims,
        array(
          rep(1 / (1 + exp(-3 * ((1:sim_list$n_ages) - 3))),
              each = sim_list$n_pop * sim_list$n_regions * sim_list$n_yrs * sim_list$n_seas,
              times = sim_list$n_sexes),
          dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes)
        )
      )
    )
  )

  ### Setup Tagging and Movement -----------------------------------------------------------
  sim_list <- Setup_Sim_Tagging(
    sim_list = sim_list,
    use_conv_fish_tagging = 1,
    n_tags = 1e3,
    conv_tag_max_liberty = 10,
    conv_fish_tag_like = "Poisson"
  )

  # Movement
  sim_list$Movement <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_regions,
                                        sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages,
                                        sim_list$n_sexes, sim_list$n_sims))

  # Fill in movement matrix
  stay_prob <- c(0.7, 0.3, 0.7)  # probability of staying in current region during dispersal season
  disperse_prob <- (1 - stay_prob) / (sim_list$n_regions - 1)  # spread remainder equally
  non_natal_rate <- 0.15 # non natal homing rate

  for (p in seq_len(sim_list$n_pop)) {
    nr <- sim_list$natal_region[p]
    for (r_from in seq_len(sim_list$n_regions)) {

      # Season 1: diffusive dispersal — mostly stay, some movement out
      for (r_to in seq_len(sim_list$n_regions)) {
        prob <- if (r_to == r_from) stay_prob[p] else disperse_prob[p]
        sim_list$Movement[p, r_from, r_to, , 1, , , ] <- prob
      }

      # Season 2: natal return with straying
      for (r_to in seq_len(sim_list$n_regions)) {
        prob <- if (r_to == nr) 1 - non_natal_rate * (sim_list$n_regions - 1) else non_natal_rate
        sim_list$Movement[p, r_from, r_to, , 2, , , ] <- prob
      }
    }
  }

  # Single season movement rates
  sim_list$sgl_seas_spawning_movement <- array(
    NA, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_regions,
                sim_list$n_yrs, sim_list$n_ages, sim_list$n_sexes, sim_list$n_sims)
  )

  ### Setup Recruitment Processes ---------------------------------------------
  sim_list <- Setup_Sim_Rec(
    sim_list = sim_list,
    R0_input = {
      R0_arry <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_sims))
      R0_arry[1, 1, , ] <- 7   # pop 1 recruits to region 1
      R0_arry[2, 1, , ] <- 7   # pop 2 recruits to region 1
      R0_arry[3, 2, , ] <- 7  # pop 3 recruits to region 2
      R0_arry
    },
    ln_sigmaR = array(log(0.5), dim = c(2, sim_list$n_pop, sim_list$n_regions)), # sigma R
    init_age_strc = "matrix", # matrix geometic series to initialize population
    recruitment_opt = 'bh_rec', # BH recruitment
    rec_dd = 'local', # localized DD
    spawn_seas = 2, # spawning season
    t_spawn = 0.5, # spawn timing
    h_input = {
      h_arry <- array(NA, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_sims))
      h_arry[1, 1, , ] <- 0.85 # pop 1 in region 1
      h_arry[2, 1, , ] <- 0.95 # pop 2 in region 1
      h_arry[3, 2, , ] <- 0.9  # pop 3 in region 2
      h_arry
    },
    stray_rate_input = array(0.5, dim = c(sim_list$n_pop, sim_list$n_yrs, sim_list$n_sims))
  )

  ## Simulate Data -----------------------------------------------------------
  sim_pop_obj <- Simulate_Pop_Static(sim_list = sim_list, output_path = NULL) # get simulated datasets

  ## Pull OM Values ----------------------------------------------------------
  n_yrs     <- sim_pop_obj$n_years
  n_regions <- sim_pop_obj$n_regions
  n_ages    <- sim_pop_obj$n_ages
  n_sexes   <- sim_pop_obj$n_sexes
  n_pop     <- sim_pop_obj$n_pop
  n_seas    <- sim_pop_obj$n_seas
  n_fish_fleets <- sim_pop_obj$n_fish_fleets
  seasdur   <- sim_pop_obj$seasdur
  spawn_seas <- sim_pop_obj$spawn_seas
  t_spawn   <- sim_pop_obj$t_spawn
  natal_region <- sim_pop_obj$natal_region
  do_recruits_move <- sim_pop_obj$do_recruits_move
  n_proj_yrs <- 1e3 # projection years

  # Terminal OM state
  terminal_NAA  <- array(sim_pop_obj$NAA[,,n_yrs,,,,1], dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes))
  terminal_NAA0 <- array(sim_pop_obj$NAA0[,,n_yrs,,,,1], dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes))

  # OM bio arrays
  natmort_slice <- sim_pop_obj$natmort[,,n_yrs,1,,,1] # season 1, M is constant within the year here
  natmort <- array(rep(c(natmort_slice), times = n_proj_yrs), dim = c(n_pop, n_regions, n_ages, n_sexes, n_proj_yrs))
  natmort <- aperm(natmort, c(1, 2, 5, 3, 4))
  # packaged report predates seasonal M, hold it across seasons
  natmort <- SPoRC:::expand_natmort_seasons(natmort, n_seas)
  WAA      <- array(rep(sim_pop_obj$WAA[,,n_yrs,,,,1], each = n_proj_yrs),
                    dim = c(n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes))
  WAA_fish <- array(rep(sim_pop_obj$WAA_fish[,,n_yrs,,,,,1], each = n_proj_yrs),
                    dim = c(n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes, n_fish_fleets))
  MatAA    <- array(rep(sim_pop_obj$MatAA[,,n_yrs,,,,1], each = n_proj_yrs),
                    dim = c(n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes))
  fish_sel <- array(rep(sim_pop_obj$fish_sel[,,n_yrs,,,,,1], each = n_proj_yrs),
                    dim = c(n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes, n_fish_fleets))
  ret_sel <- array(rep(sim_pop_obj$ret_sel[,,n_yrs,,,,,1], each = n_proj_yrs),
                    dim = c(n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes, n_fish_fleets))
  terminal_F <- array(sim_pop_obj$Fmort[,n_yrs,,,1], dim = c(n_regions, n_seas, n_fish_fleets))
  terminal_dmr <- array(sim_pop_obj$dmr[,n_yrs,,,1], dim = c(n_regions, n_seas, n_fish_fleets))
  sexratio <- array(1, dim = c(n_pop, n_regions, n_proj_yrs, n_sexes))
  rec_seas_prop <- sim_pop_obj$rec_seas_prop[,,1]
  recruitment <- array(sim_pop_obj$Rec[,,1:(n_yrs - 1), 1], dim = c(n_pop, n_regions, n_yrs - 1))
  sgl_seas_spawning_movement <- array(
    rep(sim_pop_obj$sgl_seas_spawning_movement[,,,n_yrs,,,1], each = n_proj_yrs),
    dim = c(n_pop, n_regions, n_regions, n_proj_yrs, n_ages, n_sexes)
  )
  stray_rate <- array(0.5, dim = c(n_pop, n_proj_yrs))


  # Movement
  moveslice <- sim_pop_obj$Movement[,,,n_yrs,,,,1]
  tmp <- array(rep(moveslice, n_proj_yrs),
               dim = c(n_pop, n_regions, n_regions, n_seas, n_ages, n_sexes, n_proj_yrs))
  Movement <- aperm(tmp, perm = c(1, 2, 3, 7, 4, 5, 6))

  # Get true reference points
  ref_pts <- get_closed_loop_reference_points(
    use_true_values = TRUE,
    sim_env = sim_pop_obj,
    y = sim_pop_obj$n_years,
    sim = 1,
    n_proj_yrs = 1,
    reference_points_opt = list(
      n_avg_yrs = 1,
      type = 'multi_region',
      what = "local_BH_MSY"
    )
  )

  # get BH recruitment options
  srr_opt <- list(
    rec_dd     = 0,
    rec_lag    = 1,
    R0         = apply(sim_pop_obj$R0[,, n_yrs, 1, drop = FALSE], 1, sum),
    h          = array(sim_pop_obj$h[,, n_yrs, 1], dim = c(sim_pop_obj$n_pop, sim_pop_obj$n_regions)),
    rec_region_prop = {
      R0_slice <- sim_pop_obj$R0[,, n_yrs, 1, drop = FALSE]
      row_sums <- rowSums(R0_slice)
      R0_prop <- array(R0_slice / row_sums, dim = c(sim_pop_obj$n_pop, sim_pop_obj$n_regions))
    },
    WAA        = array(sim_pop_obj$WAA[,,1,,,,1], dim = c(n_pop, n_regions, n_seas, n_ages)),
    MatAA      = array(sim_pop_obj$MatAA[,,1,,,,1], dim = c(n_pop, n_regions, n_seas, n_ages)),
    SSB        = sim_pop_obj$SSB[,,,1],
    Movement   = array(Movement[,,,1,,,], dim = c(n_pop, n_regions, n_regions, n_seas, n_ages)),
    do_recruits_move = do_recruits_move,
    sex_ratio_f = array(0.5, dim = c(n_pop, n_regions)),
    sgl_seas_spawning_movement = sgl_seas_spawning_movement[,,,1,,],
    stray_rate = array(0.5, dim = c(n_pop)),
    # M now has seasons, which Get_Det_Recruitment reads
    natmort    = array(natmort[,,1,,,], dim = c(n_pop, n_regions, n_seas, n_ages)),
    fish_sel = array(fish_sel[,,1,,,1,], dim = c(n_pop, n_regions, n_seas, n_ages, n_fish_fleets)),
    ret_sel = array(ret_sel[,,1,,,1,], dim = c(n_pop, n_regions, n_seas, n_ages, n_fish_fleets)),
    init_F = array(0, dim = c(n_regions, n_seas, n_fish_fleets)),
    dmr = array(0, dim = c(n_regions, n_seas, n_fish_fleets))
  )


  # Define HCR to use
  HCR_function <- function(x, frp, brp, alpha = 0.05) {
    stock_status <- x / brp # define stock status
    # If stock status is > 1
    if(stock_status >= 1) f <- frp
    # If stock status is between brp and alpha
    if(stock_status > alpha && stock_status < 1) f <- frp * (stock_status - alpha) / (1 - alpha)
    # If stock status is less than alpha
    if(stock_status < alpha) f <- 0
    return(f)
  }


  # do population projection
  out <- Do_Population_Projection(
    n_proj_yrs = n_proj_yrs,
    n_regions = n_regions,
    n_ages = n_ages,
    n_sexes = n_sexes,
    n_pop = n_pop,
    sexratio = sexratio,
    n_fish_fleets = n_fish_fleets,
    do_recruits_move = do_recruits_move,
    recruitment = recruitment,
    terminal_NAA = terminal_NAA,
    terminal_NAA0 = terminal_NAA0,
    terminal_F = terminal_F,
    dmr = terminal_dmr,
    natmort = natmort,
    rec_seas_prop = rec_seas_prop,
    natal_region = natal_region,
    WAA = WAA,
    WAA_fish = WAA_fish,
    MatAA = MatAA,
    fish_sel = fish_sel,
    ret_sel = ret_sel,
    Movement = Movement,
    f_ref_pt = array(ref_pts$f_ref_pt, dim = c(n_regions, n_proj_yrs)),
    b_ref_pt = NULL,
    HCR_function = HCR_function,
    recruitment_opt = "bh_rec",
    fmort_opt = "Input",
    t_spawn = t_spawn,
    spawn_seas = spawn_seas,
    n_seas = n_seas,
    sgl_seas_spawning_movement = sgl_seas_spawning_movement,
    stray_rate = stray_rate,
    seasdur = seasdur,
    srr_opt = srr_opt
  )

  # Check if F equilibriates back at F40%
  expect_equal(as.numeric(out$proj_F[,n_proj_yrs]), as.numeric(ref_pts$f_ref_pt), tolerance = 0)

  # Check if SSB equilibriates at Bx%
  expect_equal(round(as.numeric(out$proj_SSB[,,n_proj_yrs]), 10),
               round(as.numeric(ref_pts$b_ref_pt), 10), tolerance = 0)

  # Check to see if SSB equilibriates
  expect_equal(round(as.numeric(out$proj_SSB[,,n_proj_yrs]), 10),
               round(as.numeric(out$proj_SSB[,,n_proj_yrs-1]), 10), tolerance = 0)


})

