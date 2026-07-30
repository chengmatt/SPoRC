library(SPoRC)
library(testthat)
data("sgl_rg_sable_rep")
data("sgl_rg_sable_data")

test_that("Single Region BH MSY (mock) Reference Points Sablefish Model Converges to Equilibrium", {

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

  # get reference points
  sgl_ref_pt <- Get_Reference_Points(data = sgl_rg_sable_data,
                                     rep = sgl_rg_sable_rep,
                                     SPR_x = 0.4,
                                     type = 'single_region',
                                     what = 'BH_MSY'
  )

  # set up quantities to use in projection function
  n_proj_yrs <- 500 # number of projection years
  n_regions <- sgl_rg_sable_data$n_regions # number of regions
  n_ages <- length(sgl_rg_sable_data$ages) # number of ages
  n_sexes <- sgl_rg_sable_data$n_sexes # number of sexes
  n_fish_fleets <- sgl_rg_sable_data$n_fish_fleets # number of fishery fleets
  n_seas <- sgl_rg_sable_data$n_seas
  n_pop <- sgl_rg_sable_data$n_pop
  t_spawn <- 0 # spawn timing
  do_recruits_move <- 0 # recruits don't move

  # Numbers
  terminal_NAA <- array(sgl_rg_sable_rep$NAA[,,length(sgl_rg_sable_data$years),,,], dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes)) # terminal numbers at age
  terminal_NAA0 <- array(sgl_rg_sable_rep$NAA0[,,length(sgl_rg_sable_data$years),,,], dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes)) # terminal numbers at age

  # WAA
  WAA <- array(0, dim = c(n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes))
  for(y in 1:n_proj_yrs) WAA[,,y,,,] <- sgl_rg_sable_data$WAA[,,length(sgl_rg_sable_data$years),,,]
  WAA_fish <- array(0, dim = c(n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes, n_fish_fleets))
  for(y in 1:n_proj_yrs) WAA_fish[,,y,,,,] <- sgl_rg_sable_data$WAA[,,length(sgl_rg_sable_data$years),,,]

  # Maturity
  MatAA <- array(0, dim = c(n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes))
  for(y in 1:n_proj_yrs) MatAA[,,y,,,] <- sgl_rg_sable_data$MatAA[,,length(sgl_rg_sable_data$years),,,]

  # Fishery Selectivity
  fish_sel <- array(0, dim = c(n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes, n_fish_fleets))
  for(y in 1:n_proj_yrs) fish_sel[,,y,,,,] <- sgl_rg_sable_rep$fish_sel[,,length(sgl_rg_sable_data$years),,,,]

  # retained Selectivity
  ret_sel <- array(0, dim = c(n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes, n_fish_fleets))
  for(y in 1:n_proj_yrs) ret_sel[,,y,,,,] <- sgl_rg_sable_rep$ret_sel[,,length(sgl_rg_sable_data$years),,,,]

  # Movement
  Movement <- array(1, dim = c(n_pop, n_regions, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes)) # movement

  # Fishing Mortality
  terminal_F <- array(sgl_rg_sable_rep$Fmort[,length(sgl_rg_sable_data$years),,], dim = c(n_regions, n_seas, n_fish_fleets)) # terminal F

  # dmr
  terminal_dmr <- array(sgl_rg_sable_rep$dmr[,length(sgl_rg_sable_data$years),,], dim = c(n_regions, n_seas, n_fish_fleets)) # terminal F

  # Natural Mortality
  natmort_slice <- sgl_rg_sable_rep$natmort[,, length(sgl_rg_sable_data$years), , ]  # [n_pop, n_regions, n_ages, n_sexes]
  natmort <- array(rep(natmort_slice, each = n_proj_yrs), dim = c(n_pop, n_regions, n_proj_yrs, n_ages, n_sexes))

  # Recruitment
  recruitment <- array(sgl_rg_sable_rep$Rec[,,20:(length(sgl_rg_sable_data$years) - 2)],
                       dim = c(n_pop, n_regions, length(20:(length(sgl_rg_sable_data$years) - 2)))) # recruitment values to use for mean recruitment calculations or inverse gaussian parameterization

  # Sex ratio
  sexratio <- array(0.5, dim = c(n_pop, n_regions, n_proj_yrs, n_sexes)) # recruitment sex ratio

  # storage
  sgl_f_proj <- array(0, dim = c(n_regions, n_proj_yrs))
  sgl_ssb_proj <- array(0, dim = c(n_regions, n_proj_yrs))
  sgl_catch_proj <- array(0, dim = c(n_regions, n_proj_yrs, n_fish_fleets))

  # BH recruitment options
  bh_rec_opt <- list(
    rec_dd = 1,
    rec_lag = 2,
    do_recruits_move = do_recruits_move,
    R0 = sgl_rg_sable_rep$R0,
    h = array(sgl_rg_sable_rep$h_trans, dim = c(sgl_rg_sable_data$n_pop, sgl_rg_sable_data$n_regions)),
    rec_region_prop = sgl_rg_sable_rep$rec_region_prop,
    WAA = array(sgl_rg_sable_data$WAA[,,1,,,1], dim = c(sgl_rg_sable_data$n_pop,sgl_rg_sable_data$n_regions,sgl_rg_sable_data$n_seas,length(sgl_rg_sable_data$ages))),
    MatAA = array(sgl_rg_sable_data$MatAA[,,1,,,1], dim = c(sgl_rg_sable_data$n_pop,sgl_rg_sable_data$n_regions,sgl_rg_sable_data$n_seas,length(sgl_rg_sable_data$ages)) ),
    SSB = sgl_rg_sable_rep$SSB,
    Movement = array(Movement[,,,1,,,1], dim = c(sgl_rg_sable_data$n_pop,sgl_rg_sable_data$n_regions,sgl_rg_sable_data$n_regions, sgl_rg_sable_data$n_seas,length(sgl_rg_sable_data$ages)) ),
    sex_ratio_f = array(0.5, dim = c(sgl_rg_sable_data$n_pop, n_regions)),
    sgl_seas_spawning_movement = NULL,
    stray_rate = array(0, dim = c(sgl_rg_sable_data$n_pop)),
    natmort = array(natmort[,,1,,1], dim = c(sgl_rg_sable_data$n_pop, sgl_rg_sable_data$n_regions, length(sgl_rg_sable_data$ages) )),
    fish_sel = array(fish_sel[,,1,,,1,], dim = c(sgl_rg_sable_data$n_pop, sgl_rg_sable_data$n_regions, sgl_rg_sable_data$n_seas, length(sgl_rg_sable_data$ages), sgl_rg_sable_data$n_fish_fleets)),
    ret_sel = array(ret_sel[,,1,,,1,], dim = c(sgl_rg_sable_data$n_pop, sgl_rg_sable_data$n_regions, sgl_rg_sable_data$n_seas, length(sgl_rg_sable_data$ages), sgl_rg_sable_data$n_fish_fleets)),
    init_F = array(0, dim = c(sgl_rg_sable_data$n_regions, sgl_rg_sable_data$n_seas, sgl_rg_sable_data$n_fish_fleets)),
    dmr = array(0, dim = c(sgl_rg_sable_data$n_regions, sgl_rg_sable_data$n_seas, sgl_rg_sable_data$n_fish_fleets))
  )


  # do projection
  out <- Do_Population_Projection(n_proj_yrs = n_proj_yrs,
                                  n_regions = n_regions,
                                  n_ages = n_ages,
                                  n_sexes = n_sexes,
                                  sexratio = sexratio,
                                  n_fish_fleets = n_fish_fleets,
                                  do_recruits_move = do_recruits_move,
                                  recruitment = recruitment,
                                  terminal_NAA = terminal_NAA,
                                  terminal_NAA0 = terminal_NAA0,
                                  terminal_F = terminal_F,
                                  dmr = terminal_dmr,
                                  natmort = natmort,
                                  WAA = WAA,
                                  n_pop = n_pop,
                                  WAA_fish = WAA_fish,
                                  MatAA = MatAA,
                                  fish_sel = fish_sel,
                                  ret_sel = ret_sel,
                                  Movement = Movement,
                                  f_ref_pt = array(sgl_ref_pt$f_ref_pt, dim = c(n_regions, n_proj_yrs)),
                                  b_ref_pt = array(sgl_ref_pt$b_ref_pt, dim = c(n_pop, n_regions, n_proj_yrs)),
                                  HCR_function = HCR_function,
                                  recruitment_opt = "bh_rec",
                                  fmort_opt = "HCR",
                                  t_spawn = t_spawn,
                                  bh_rec_opt = bh_rec_opt
  )

  # Check if F equilibriates back at F40%
  expect_equal(out$proj_F[,n_proj_yrs], sgl_ref_pt$f_ref_pt, tolerance = 0)

  # Check if SSB equilibriates at Bx%
  expect_equal(round(as.numeric(out$proj_SSB[1,1,n_proj_yrs]), 10),
               round(as.numeric(sgl_ref_pt$b_ref_pt), 10), tolerance = 0)

  # Check to see if SSB equilibriates
  expect_equal(round(as.numeric(out$proj_SSB[1,1,n_proj_yrs]), 10),
               round(as.numeric(out$proj_SSB[1,1,n_proj_yrs-1]), 10), tolerance = 0)


})
