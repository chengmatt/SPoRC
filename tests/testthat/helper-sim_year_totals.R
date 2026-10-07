# Seasonal and population-specific models the self test, year total and lockstep tests share: totals reported
# once a year in the season the fit holds them in, and at-age data by population and as year totals.

year_total_dims <- list(n_pop = 2, n_regions = 2, n_seas = 2, n_yrs = 8, n_ages = 5, n_sexes = 1)

# use flags on in one season, everywhere else off
use_in_season <- function(dims, seas, pop = FALSE) {
  use <- array(0, dim = dims)
  if(pop) use[,,,seas,] <- 1 else use[,,seas,] <- 1
  use
}

# Catch and discard fractions kept in season two, population discards every season, and three
# compositions reported once a year: by population in season two for the fishery and the survey,
# and the discards' in season one
year_total_input <- function() {

  d <- year_total_dims
  agg_dim <- c(d$n_regions, d$n_yrs, d$n_seas, 1)
  pop_dim <- c(d$n_pop, agg_dim)
  comp_pop_dim <- c(d$n_pop, d$n_regions, d$n_yrs, d$n_seas, d$n_ages, d$n_sexes, 1)
  iss_pop_dim <- c(d$n_pop, d$n_regions, d$n_yrs, d$n_seas, d$n_sexes, 1)

  suppressWarnings(suppressMessages(sweep_input(
    dims = list(n_pop = d$n_pop, n_regions = d$n_regions, n_yrs = d$n_yrs, n_ages = d$n_ages, n_seas = d$n_seas,
                n_sexes = d$n_sexes, n_fish_fleets = 1, n_srv_fleets = 1, natal_region = seq_len(d$n_pop)),
    rec = list(rec_dd = "local", ln_global_R0 = rep(log(1e6), d$n_pop)),
    catch = list(
      ObsCatch = array(1e4, dim = agg_dim),
      UseCatch = use_in_season(agg_dim, 2),
      Catch_seas_Type = "aggSeas",
      ObsDiscard = array(0.2, dim = agg_dim),
      UseDiscard = use_in_season(agg_dim, 2),
      Discard_seas_Type = "aggSeas",
      ObsDiscard_pop = array(0.2, dim = pop_dim),
      UseDiscard_pop = array(1, dim = pop_dim),
      discard_units = "abd_frac",
      sigmaD_spec = "fix",
      sigmaD_pop_spec = "fix"
    ),
    fishidx = list(
      ObsFishAgeComps_pop = array(1 / d$n_ages, dim = comp_pop_dim),
      UseFishAgeComps_pop = use_in_season(pop_dim, 2, pop = TRUE),
      ISS_FishAgeComps_pop = array(100, dim = iss_pop_dim),
      FishAgeComps_pop_LikeType = "Multinomial",
      FishAgeComps_pop_Type = "spltRspltS_Year_1-terminal_Fleet_1",
      FishAgeComps_pop_seas_Type = "aggSeas",
      ObsFishAgeComps_discard = array(1 / d$n_ages, dim = c(d$n_regions, d$n_yrs, d$n_seas, d$n_ages, d$n_sexes, 1)),
      UseFishAgeComps_discard = use_in_season(agg_dim, 1),
      ISS_FishAgeComps_discard = array(80, dim = c(d$n_regions, d$n_yrs, d$n_seas, d$n_sexes, 1)),
      FishAgeComps_discard_LikeType = "Multinomial",
      FishAgeComps_discard_Type = "spltRspltS_Year_1-terminal_Fleet_1",
      FishAgeComps_discard_seas_Type = "aggSeas"
    ),
    srvidx = list(
      ObsSrvAgeComps_pop = array(1 / d$n_ages, dim = comp_pop_dim),
      UseSrvAgeComps_pop = use_in_season(pop_dim, 2, pop = TRUE),
      ISS_SrvAgeComps_pop = array(60, dim = iss_pop_dim),
      SrvAgeComps_pop_LikeType = "Multinomial",
      SrvAgeComps_pop_Type = "spltRspltS_Year_1-terminal_Fleet_1",
      SrvAgeComps_pop_seas_Type = "aggSeas"
    ),
    fishsel = list(use_fixed_ret_sel = 0, ret_sel_model = "logist1_Fleet_1", ret_fixed_sel_pars_spec = "est_all") # some fish discarded
  )))
}

pop_seas_dims <- list(n_pop = 2, n_regions = 2, n_seas = 2, n_yrs = 8, n_ages = 5, n_sexes = 1, n_fish = 2, n_srv = 1)

# Fleet 1 reports catch at age by population every season, fleet 2 a year total of catch at age in
# season one, and the survey a year total of its index at age by population in season two
pop_seas_at_age_input <- function() {

  d <- pop_seas_dims
  aa_fish <- c(d$n_regions, d$n_yrs, d$n_seas, d$n_ages, d$n_sexes, d$n_fish)
  aa_srv <- c(d$n_regions, d$n_yrs, d$n_seas, d$n_ages, d$n_sexes, d$n_srv)

  use_caa <- array(0, dim = aa_fish)
  use_caa[,,1,,1,2] <- 1
  use_caa_pop <- array(0, dim = c(d$n_pop, aa_fish))
  use_caa_pop[,,,,,1,1] <- 1
  use_saa_pop <- array(0, dim = c(d$n_pop, aa_srv))
  use_saa_pop[,,,2,,1,1] <- 1

  il <- sweep_input(
    dims = list(n_pop = d$n_pop, n_regions = d$n_regions, n_yrs = d$n_yrs, n_ages = d$n_ages, n_seas = d$n_seas,
                n_sexes = d$n_sexes, n_fish_fleets = d$n_fish, n_srv_fleets = d$n_srv, natal_region = seq_len(d$n_pop)),
    rec = list(rec_dd = "local", ln_global_R0 = rep(log(1e6), d$n_pop)),
    catch = list(
      UseCatch = array(0, dim = c(d$n_regions, d$n_yrs, d$n_seas, d$n_fish)), # a fleet fits catch at age or catch, not both
      ObsCatchAA = array(1e3, dim = aa_fish),
      UseCatchAA = use_caa,
      CatchAA_seas_Type = c("spltSeas", "aggSeas"),
      sigmaCAA_spec = "fix",
      ObsCatchAA_pop = array(1e3, dim = c(d$n_pop, aa_fish)),
      UseCatchAA_pop = use_caa_pop,
      sigmaCAA_pop_spec = "fix"
    ),
    srvidx = list(
      UseSrvIdx = array(0, dim = c(d$n_regions, d$n_yrs, d$n_seas, d$n_srv)),
      ObsSrvIdxAA_pop = array(1e3, dim = c(d$n_pop, aa_srv)),
      UseSrvIdxAA_pop = use_saa_pop,
      SrvIdxAA_pop_seas_Type = "aggSeas",
      sigmaSrvIdxAA_pop_spec = "fix"
    )
  )

  list(input = il, use_caa = use_caa, use_caa_pop = use_caa_pop, use_saa_pop = use_saa_pop)
}
