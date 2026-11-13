# The growth module against the growth branch of WHAM (timjmiller/wham, branch growth). The
# reference values come from WHAM's own C++, built by dev/wham_growth_bridge/make_wham_growth_reference.R.

library(SPoRC)
library(testthat)

wham_ref <- readRDS(test_path("wham_growth_reference.rds"))
wham_settings <- wham_ref$settings
n_ages <- wham_settings$n_ages
ages <- 1:n_ages
wham_to_sporc_par <- c(3, 2, 1) # WHAM orders K, Linf, L1; this model L1, L2, K

# WHAM's growth in this model's terms: L1 at age one, L2 the asymptote itself, the two spreads read
# as standard deviations interpolated on length, the plus group the decay-weighted mixture
wham_growth_call <- function(case, n_yrs = 1, t_srv = wham_settings$fracyr_index, ...) {

  ln_pars <- log(c(wham_settings$L1, wham_settings$Linf, wham_settings$K, wham_settings$SD1, wham_settings$SDA))
  if(case$growth_model == 2) ln_pars <- c(ln_pars, log(wham_settings$gamma)) # Richards shape last
  args <- list(
    ln_growth_pars = array(ln_pars, dim = c(1, 1, 1, length(ln_pars))),
    growth_A1 = 1, # WHAM's age_L1
    growth_A2 = n_ages, # under the asymptote form the accumulator age bounds the spread
    growth_L2_asymptote = 1,
    growth_L0 = wham_settings$lengths[1], # WHAM's linear phase starts at the first lower edge
    growth_len_lower = wham_settings$lengths,
    growth_cv_type = 0, # interpolated on length
    growth_sd_type = 1, # SD1 and SDA are standard deviations
    growth_dist = 0,
    growth_plus_group = 1,
    derive_waa = 1,
    wt_len_pars = array(wham_settings$lw, dim = c(1, 1, 1, 2)),
    ages = ages,
    seasdur = 1,
    spawn_seas = 1,
    t_spawn = wham_settings$fracyr_SSB,
    n_pop = 1,
    n_regions = 1,
    n_yrs = n_yrs,
    n_seas = 1,
    n_sexes = 1,
    n_fish_fleets = 1, # read at mid season, WHAM's fracyr_catch
    n_srv_fleets = 1,
    t_srv = array(t_srv, dim = c(1, 1, 1))
  )
  do.call(Get_Growth, utils::modifyList(args, list(...)))
}

# cohort growth year by year under one case's year deviations
wham_cohort_call <- function(case) {

  n_yrs <- wham_settings$n_yrs
  n_growth_pars <- if(case$growth_model == 2) 6 else 5
  par_idx <- wham_to_sporc_par[case$dev_par]
  devs <- array(0, dim = c(1, 1, n_yrs, n_growth_pars, 1))
  devs[1, 1, , par_idx, 1] <- case$dev # log link, as WHAM's exp(growth_a + growth_re)
  tv_model <- rep(0, n_growth_pars)
  tv_model[par_idx] <- 1

  growth <- wham_growth_call(case, n_yrs = n_yrs, t_srv = 0, ln_growth_devs = devs, growth_tv_model = tv_model,
                             growth_tv_type = 1, growth_cohort_styr = 1)
  NAA_y <- array(exp(-0.3 * ages), dim = c(1, 1, n_ages, 1)) # read only by the plus group blend
  for(y in 1:n_yrs) {
    growth <- Get_Growth_Year(
      growth = growth,
      y = y,
      NAA_y = NAA_y,
      ln_growth_pars = array(log(c(wham_settings$L1, wham_settings$Linf, wham_settings$K, wham_settings$SD1, wham_settings$SDA,
                                   if(case$growth_model == 2) wham_settings$gamma)), dim = c(1, 1, 1, n_growth_pars)),
      ln_growth_devs = devs,
      growth_tv_model = tv_model,
      growth_tv_link = 0,
      growth_par_bounds = matrix(0, n_growth_pars, 2),
      growth_A1 = 1,
      growth_A2 = n_ages,
      growth_L2_asymptote = 1,
      growth_L0 = wham_settings$lengths[1],
      growth_len_lower = wham_settings$lengths,
      growth_cv_type = 0,
      growth_sd_type = 1,
      growth_dist = 0,
      growth_plus_group = 1,
      derive_waa = 1,
      wt_len_pars = array(wham_settings$lw, dim = c(1, 1, 1, 2)),
      ages = ages,
      seasdur = 1,
      spawn_seas = 1,
      t_spawn = wham_settings$fracyr_SSB,
      n_pop = 1,
      n_regions = 1,
      n_seas = 1,
      n_sexes = 1,
      t_srv = array(0, dim = c(1, 1, 1))
    )
  } # end y loop
  growth
}

# Constant Growth ------------------------------------------------------------

test_that("constant growth matches WHAM at the start of the year, plus group included", {
  for(case_name in c("vb_const", "richards_const")) {
    case <- wham_ref$cases[[case_name]]
    g0 <- wham_growth_call(case, t_srv = 0) # a survey at the start of the year reads WHAM's January key
    expect_equal(g0$mean_LAA_srv[1, 1, 1, 1, , 1, 1], case$LAA[1, ], tolerance = 1e-10)
    expect_equal(g0$sd_LAA_srv[1, 1, 1, 1, , 1, 1], case$SDAA[1, ], tolerance = 1e-10)
    expect_equal(g0$SizeAgeTrans_srv[1, 1, 1, 1, , , 1, 1], case$phi$jan1, tolerance = 1e-10)
    expect_equal(g0$WAA_srv[1, 1, 1, 1, , 1, 1], case$waa$jan1, tolerance = 1e-10)
  } # end case loop
})

test_that("within the year every age below the plus group matches WHAM's key and weight", {
  below_plus <- 1:(n_ages - 1)
  for(case_name in c("vb_const", "richards_const")) {
    case <- wham_ref$cases[[case_name]]
    g <- wham_growth_call(case)

    # spawning, WHAM's fracyr_SSB
    expect_equal(g$SizeAgeTrans_spawn[1, 1, 1, , below_plus, 1], case$phi$ssb[, below_plus], tolerance = 1e-10)
    expect_equal(g$WAA[1, 1, 1, 1, below_plus, 1], case$waa$ssb[below_plus], tolerance = 1e-10)
    # the fishery at mid season, WHAM's fracyr_catch
    expect_equal(g$SizeAgeTrans_fish[1, 1, 1, 1, , below_plus, 1, 1], case$phi$catch[, below_plus], tolerance = 1e-10)
    expect_equal(g$WAA_fish[1, 1, 1, 1, below_plus, 1, 1], case$waa$catch[below_plus], tolerance = 1e-10)
    # the survey, WHAM's fracyr_indices
    expect_equal(g$SizeAgeTrans_srv[1, 1, 1, 1, , below_plus, 1, 1], case$phi$index[, below_plus], tolerance = 1e-10)
    expect_equal(g$WAA_srv[1, 1, 1, 1, below_plus, 1, 1], case$waa$index[below_plus], tolerance = 1e-10)
  } # end case loop
})

test_that("within the year the plus group grows from its adjusted size, where WHAM keeps its January size", {
  for(case_name in c("vb_const", "richards_const")) {
    case <- wham_ref$cases[[case_name]]
    rho <- if(case$growth_model == 2) wham_settings$gamma else 1
    g <- wham_growth_call(case)

    # WHAM: the plus group's mean and spread at the survey are those of January 1, so its key is too
    expect_equal(case$phi$index[, n_ages], case$phi$jan1[, n_ages], tolerance = 1e-12)
    # this model: the January size grown by the curve's increment over the time elapsed
    L_plus_jan1 <- case$LAA[1, n_ages]
    L_plus_srv <- grow_increment(L_plus_jan1, wham_settings$fracyr_index, wham_settings$K, wham_settings$Linf, rho)
    expect_equal(g$mean_LAA_srv[1, 1, 1, 1, n_ages, 1, 1], L_plus_srv, tolerance = 1e-10)
    expect_gt(L_plus_srv, L_plus_jan1)
  } # end case loop
})

# Year Deviations Under Cohort Growth -----------------------------------------

test_that("cohort growth under year deviations on K or Linf matches WHAM below the plus group", {
  below_plus <- 1:(n_ages - 1)
  for(case_name in c("vb_K_y", "vb_Linf_y", "richards_K_y")) {
    case <- wham_ref$cases[[case_name]]
    g <- wham_cohort_call(case)
    expect_equal(g$L_beg[1, 1, , below_plus, 1], case$LAA[, below_plus], tolerance = 1e-10)
  } # end case loop
})

test_that("under year deviations on L1 every age past the first matches WHAM", {
  # age one is left out: at the start of the year this model gives the age at A1 the L1 of its
  # birth year, where WHAM gives it the current year's L1
  case <- wham_ref$cases$vb_L1_y
  g <- wham_cohort_call(case)
  expect_equal(g$L_beg[1, 1, , 2:(n_ages - 1), 1], case$LAA[, 2:(n_ages - 1)], tolerance = 1e-10)
})

test_that("WHAM rebuilds the plus group's size and the spread at age every year, where this model advances both", {
  case <- wham_ref$cases$vb_K_y
  n_yrs <- wham_settings$n_yrs
  K_y <- wham_settings$K * exp(case$dev)
  Linf <- wham_settings$Linf

  # WHAM: the cohort entering the plus group, grown by last year's K, then the decay-weighted mixture
  for(y in 2:n_yrs) {
    entering <- grow_increment(case$LAA[y - 1, n_ages - 1], 1, K_y[y - 1], Linf)
    expect_equal(case$LAA[y, n_ages], plus_group_size(entering, Linf, n_ages), tolerance = 1e-10)
  } # end y loop
  # WHAM: the spread interpolated on the current year's mean length between L1 and Linf
  mid <- 2:(n_ages - 1)
  sd_interp <- wham_settings$SD1 + (case$LAA[, mid] - wham_settings$L1) * (wham_settings$SDA - wham_settings$SD1) / (Linf - wham_settings$L1)
  expect_equal(case$SDAA[, mid], sd_interp, tolerance = 1e-10)

  # this model blends the plus group by numbers and keeps the first year's spread at age, so both
  # differ from WHAM once a deviation has moved the means
  g <- wham_cohort_call(case)
  for(y in 2:n_yrs) expect_equal(g$sd_LAA_srv[1, 1, y, 1, , 1, 1], g$sd_LAA_srv[1, 1, 1, 1, , 1, 1], tolerance = 1e-12)
  expect_false(isTRUE(all.equal(g$sd_LAA_srv[1, 1, n_yrs, 1, mid, 1, 1], case$SDAA[n_yrs, mid])))
  expect_false(isTRUE(all.equal(g$L_beg[1, 1, n_yrs, n_ages, 1], case$LAA[n_yrs, n_ages])))
})
