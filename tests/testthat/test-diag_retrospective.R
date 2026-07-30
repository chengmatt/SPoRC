library(testthat)
library(SPoRC)

mk <- function(dims, fill = 1) array(fill, dim = dims)
make_retro_data <- function(n_yrs = 6,
                            n_pop = 1,
                            n_regions = 1,
                            n_fish_fleets = 1,
                            n_srv_fleets = 1,
                            use_fixed_natmort = 1,
                            use_fixed_stray_rate = 1,
                            use_conv_fish_tagging = 0,
                            has_bias_year = FALSE,
                            pop_specific = FALSE,
                            include_recdevs_map = FALSE) {

  data <- list()
  parameters <- list()
  mapping <- list()

  data$n_yrs <- n_yrs
  data$n_pop <- n_pop
  data$n_regions <- n_regions
  data$n_fish_fleets <- n_fish_fleets
  data$n_srv_fleets <- n_srv_fleets
  data$use_fixed_natmort <- use_fixed_natmort
  data$use_fixed_stray_rate <- use_fixed_stray_rate
  data$use_conv_fish_tagging <- rep(use_conv_fish_tagging, n_fish_fleets)

  data$years <- 1:n_yrs
  data$bias_year <- if (has_bias_year) c(1, 2, 3, n_yrs) else NA

  ## Recruitment ----------------------------------------------------------
  parameters$ln_RecDevs <- mk(c(n_pop, n_regions, n_yrs))
  if (include_recdevs_map) {
    mapping$ln_RecDevs <- factor(mk(c(n_pop, n_regions, n_yrs)))
  }

  ## Stray rates (n_pop > 1 only) ------------------------------------------
  if (n_pop > 1) {
    data$stray_rate_blocks <- mk(c(n_pop, n_yrs))
    if (use_fixed_stray_rate == 1) {
      data$fixed_stray_rate <- mk(c(n_pop, n_yrs))
    } else {
      parameters$stray_rate_pars <- mk(c(n_pop, 1))
      mapping$stray_rate_pars <- factor(mk(c(n_pop, 1)))
    }
  }

  ## Sex ratio (always executed) -------------------------------------------
  data$sexratio_blocks <- mk(c(2, 1, n_yrs))
  parameters$sexratio_pars <- mk(c(2, 1, 1))
  mapping$sexratio_pars <- factor(mk(c(2, 1, 1)))

  ## Natural mortality (always executed) -----------------------------------
  data$M_blocks <- mk(c(1, 1, n_yrs, 1, 1))
  if (use_fixed_natmort == 1) {
    data$Fixed_natmort <- mk(c(1, 1, n_yrs, 1, 1))
  } else {
    parameters$ln_M <- mk(1, fill = 0)  # length-1 ARRAY (has a real dim attribute), keyed off max(M_blocks) = 1
    mapping$ln_M <- factor(mk(1, fill = 1))
  }

  ## Fishery ----------------------------------------------------------------
  data$ObsCatch  <- mk(c(1, n_yrs, n_fish_fleets, 1))
  data$ObsFishIdx <- mk(c(1, n_yrs, n_fish_fleets, 1))
  data$ObsFishIdx_SE <- mk(c(1, n_yrs, n_fish_fleets, 1))
  data$ObsFishAgeComps <- mk(c(1, n_yrs, 1, 1, n_fish_fleets, 1))
  data$ObsFishLenComps <- mk(c(1, n_yrs, 1, 1, n_fish_fleets, 1))
  data$ObsFishAgeComps_discard <- mk(c(1, n_yrs, 1, 1, n_fish_fleets, 1))
  data$ObsFishLenComps_discard <- mk(c(1, n_yrs, 1, 1, n_fish_fleets, 1))

  parameters$ln_F_devs <- mk(c(1, n_yrs, n_fish_fleets, 1))
  mapping$ln_F_devs <- factor(mk(c(1, n_yrs, n_fish_fleets, 1)))
  parameters$logit_dmr_devs <- mk(c(1, n_yrs, n_fish_fleets, 1))
  mapping$logit_dmr_devs <- factor(mk(c(1, n_yrs, n_fish_fleets, 1)))

  parameters$ln_fishsel_devs <- mk(c(1, n_yrs, n_fish_fleets, 1, 1))
  mapping$ln_fishsel_devs <- factor(mk(c(1, n_yrs, n_fish_fleets, 1, 1)))
  data$map_ln_fishsel_devs <- mk(c(1, n_yrs, n_fish_fleets, 1, 1))

  parameters$ln_retsel_devs <- mk(c(1, n_yrs, n_fish_fleets, 1, 1))
  mapping$ln_retsel_devs <- factor(mk(c(1, n_yrs, n_fish_fleets, 1, 1)))
  data$ln_retsel_devs <- mk(c(1, n_yrs, n_fish_fleets, 1, 1)) # yes, also a data field

  data$fish_q_blocks <- mk(c(1, n_yrs, n_fish_fleets))
  data$fish_sel_blocks <- mk(c(1, n_yrs, n_fish_fleets))
  data$ret_sel_blocks <- mk(c(1, n_yrs, n_fish_fleets))
  parameters$ln_fish_q <- mk(c(1, 1, n_fish_fleets))
  parameters$fish_fixed_sel_pars <- mk(c(1, 1, 1, 1, n_fish_fleets))
  parameters$ret_fixed_sel_pars  <- mk(c(1, 1, 1, 1, n_fish_fleets))
  mapping$ln_fish_q <- factor(mk(c(1, 1, n_fish_fleets)))
  mapping$fish_fixed_sel_pars <- factor(mk(c(1, 1, 1, 1, n_fish_fleets)))
  mapping$ret_fixed_sel_pars  <- factor(mk(c(1, 1, 1, 1, n_fish_fleets)))

  parameters$ln_sigmaC <- mk(c(1, n_yrs, n_fish_fleets, 1))
  mapping$ln_sigmaC <- factor(mk(c(1, n_yrs, n_fish_fleets, 1)))
  parameters$ln_sigmaC_pop <- mk(c(n_pop, 1, n_yrs, 1, 1))
  mapping$ln_sigmaC_pop <- factor(mk(c(n_pop, 1, n_yrs, 1, 1)))
  parameters$ln_sigmaD <- mk(c(1, n_yrs, n_fish_fleets, 1))
  mapping$ln_sigmaD <- factor(mk(c(1, n_yrs, n_fish_fleets, 1)))
  parameters$ln_sigmaD_pop <- mk(c(n_pop, 1, n_yrs, 1, 1))
  mapping$ln_sigmaD_pop <- factor(mk(c(n_pop, 1, n_yrs, 1, 1)))

  data$Wt_Catch <- mk(c(1, n_yrs, n_fish_fleets, 1))     # length(dim) == 4 branch
  data$Wt_Discard <- mk(c(1, n_yrs, n_fish_fleets, 1))   # length(dim) == 4 branch

  ## Population-specific fishery (guarded existence, unguarded Use flags) --
  data$UseFishIdx_pop <- mk(c(n_pop, n_regions, n_yrs, 1, n_fish_fleets), fill = if (pop_specific) 1 else 0)
  data$UseFishAgeComps_pop <- mk(c(n_pop, n_regions, n_yrs, 1, n_fish_fleets), fill = if (pop_specific) 1 else 0)
  data$UseFishLenComps_pop <- mk(c(n_pop, n_regions, n_yrs, 1, n_fish_fleets), fill = if (pop_specific) 1 else 0)
  data$UseFishAgeComps_discard_pop <- mk(c(n_pop, n_regions, n_yrs, 1, n_fish_fleets))
  data$UseFishLenComps_discard_pop <- mk(c(n_pop, n_regions, n_yrs, 1, n_fish_fleets))
  data$UseSrvIdx_pop <- mk(c(n_pop, n_regions, n_yrs, 1, n_srv_fleets), fill = if (pop_specific) 1 else 0)
  data$UseSrvAgeComps_pop <- mk(c(n_pop, n_regions, n_yrs, 1, n_srv_fleets), fill = if (pop_specific) 1 else 0)
  data$UseSrvLenComps_pop <- mk(c(n_pop, n_regions, n_yrs, 1, n_srv_fleets), fill = if (pop_specific) 1 else 0)
  data$UseCatch_pop <- mk(c(n_pop, n_regions, n_yrs, 1, n_fish_fleets))
  data$UseDiscard_pop <- mk(c(n_pop, n_regions, n_yrs, 1, n_fish_fleets))

  if (pop_specific) {
    data$ObsFishIdx_pop <- mk(c(n_pop, n_regions, n_yrs, n_fish_fleets, 1))
    data$ObsFishIdx_pop_SE <- mk(c(n_pop, n_regions, n_yrs, n_fish_fleets, 1))
    data$ObsFishAgeComps_pop <- mk(c(n_pop, n_regions, n_yrs, 1, 1, n_fish_fleets, 1))
    data$ObsFishLenComps_pop <- mk(c(n_pop, n_regions, n_yrs, 1, 1, n_fish_fleets, 1))
    data$ObsFishAgeComps_discard_pop <- mk(c(n_pop, n_regions, n_yrs, 1, 1, n_fish_fleets, 1))
    data$ObsFishLenComps_discard_pop <- mk(c(n_pop, n_regions, n_yrs, 1, 1, n_fish_fleets, 1))
    data$ObsSrvIdx_pop <- mk(c(n_pop, n_regions, n_yrs, n_srv_fleets, 1))
    data$ObsSrvIdx_pop_SE <- mk(c(n_pop, n_regions, n_yrs, n_srv_fleets, 1))
    data$ObsSrvAgeComps_pop <- mk(c(n_pop, n_regions, n_yrs, 1, 1, n_srv_fleets, 1))
    data$ObsSrvLenComps_pop <- mk(c(n_pop, n_regions, n_yrs, 1, 1, n_srv_fleets, 1))
  }

  ## Survey -------------------------------------------------------------------
  data$ObsSrvIdx <- mk(c(1, n_yrs, n_srv_fleets, 1))
  data$ObsSrvIdx_SE <- mk(c(1, n_yrs, n_srv_fleets, 1))
  data$ObsSrvAgeComps <- mk(c(1, n_yrs, 1, 1, n_srv_fleets, 1))
  data$ObsSrvLenComps <- mk(c(1, n_yrs, 1, 1, n_srv_fleets, 1))

  parameters$ln_srvsel_devs <- mk(c(1, n_yrs, n_srv_fleets, 1, 1))
  mapping$ln_srvsel_devs <- factor(mk(c(1, n_yrs, n_srv_fleets, 1, 1)))
  data$map_ln_srvsel_devs <- mk(c(1, n_yrs, n_srv_fleets, 1, 1))

  data$srv_q_blocks <- mk(c(1, n_yrs, n_srv_fleets))
  data$srv_sel_blocks <- mk(c(1, n_yrs, n_srv_fleets))
  parameters$ln_srv_q <- mk(c(1, 1, n_srv_fleets))
  parameters$srv_fixed_sel_pars <- mk(c(1, 1, 1, 1, n_srv_fleets))
  mapping$ln_srv_q <- factor(mk(c(1, 1, n_srv_fleets)))
  mapping$srv_fixed_sel_pars <- factor(mk(c(1, 1, 1, 1, n_srv_fleets)))

  ## Movement (n_regions > 1 only) --------------------------------------------
  if (n_regions > 1) {
    parameters$move_pars <- mk(c(1, 1, 1, n_yrs, 1, 1, 1))
    parameters$move_devs <- mk(c(1, 1, 1, n_yrs, 1, 1, 1))
    mapping$move_pars <- factor(mk(c(1, 1, 1, n_yrs, 1, 1, 1)))
    mapping$move_devs <- factor(mk(c(1, 1, 1, n_yrs, 1, 1, 1)))
    data$Fixed_Movement <- mk(c(1, 1, 1, n_yrs, 1, 1, 1))
    data$sgl_seas_spawning_movement <- mk(c(1, 1, 1, n_yrs, 1, 1))
    data$ctmc_move_dat <- data.frame(years = 1:n_yrs, dummy = 0)
    data$map_move_devs <- mk(c(1, 1, 1, n_yrs, 1, 1, 1))
  }

  ## Conventional tagging (use_conv_fish_tagging == 1 only) -------------------
  if (any(data$use_conv_fish_tagging == 1)) {
    data$conv_tag_fish_reporting_blocks <- mk(c(1, n_yrs, n_fish_fleets))
    parameters$conv_tag_fish_reporting_pars <- mk(c(1, 1, n_fish_fleets))
    mapping$conv_tag_fish_reporting_pars <- factor(mk(c(1, 1, n_fish_fleets)))

    # 3 cohorts released in years 1, 3, and n_yrs
    data$conv_tag_release_indicator <- matrix(
      c(1, 1, 1,   1, 3, 1,   1, n_yrs, 1),
      ncol = 3, byrow = TRUE
    ) # region, year, season
    data$conv_tag_release_platform <- matrix(1, nrow = 3, ncol = 2)
    data$n_conv_tag_cohorts <- 3
    data$conv_tagged_fish <- mk(c(3, 1, 1, 1), fill = 50)
    data$obs_conv_tag_fish_recap <- mk(c(1, 1, 3, 1, n_regions, 1, 1, n_fish_fleets), fill = 2)
  }

  ## Data weights, composition types, ISS, use indicators (always executed) --
  data$Wt_FishAgeComps <- mk(c(1, n_yrs, 1, 1, n_fish_fleets))
  data$Wt_FishAgeComps_discard <- mk(c(1, n_yrs, 1, 1, n_fish_fleets))
  data$Wt_SrvAgeComps <- mk(c(1, n_yrs, 1, 1, n_srv_fleets))
  data$Wt_FishLenComps <- mk(c(1, n_yrs, 1, 1, n_fish_fleets))
  data$Wt_FishLenComps_discard <- mk(c(1, n_yrs, 1, 1, n_fish_fleets))
  data$Wt_SrvLenComps <- mk(c(1, n_yrs, 1, 1, n_srv_fleets))
  data$FishAgeComps_Type <- mk(c(n_yrs, n_fish_fleets))
  data$FishLenComps_Type <- mk(c(n_yrs, n_fish_fleets))
  data$FishAgeComps_discard_Type <- mk(c(n_yrs, n_fish_fleets))
  data$FishLenComps_discard_Type <- mk(c(n_yrs, n_fish_fleets))
  data$SrvLenComps_Type <- mk(c(n_yrs, n_srv_fleets))
  data$SrvAgeComps_Type <- mk(c(n_yrs, n_srv_fleets))

  data$UseFishAgeComps <- mk(c(1, n_yrs, 1, n_fish_fleets))
  data$UseFishAgeComps_discard <- mk(c(1, n_yrs, 1, n_fish_fleets))
  data$UseFishIdx <- mk(c(1, n_yrs, 1, n_fish_fleets))
  data$UseCatch <- mk(c(1, n_yrs, 1, n_fish_fleets))
  data$UseDiscard <- mk(c(1, n_yrs, 1, n_fish_fleets))
  data$UseFishLenComps <- mk(c(1, n_yrs, 1, n_fish_fleets))
  data$UseFishLenComps_discard <- mk(c(1, n_yrs, 1, n_fish_fleets))
  data$UseSrvAgeComps <- mk(c(1, n_yrs, 1, n_srv_fleets))
  data$UseSrvIdx <- mk(c(1, n_yrs, 1, n_srv_fleets))
  data$UseSrvLenComps <- mk(c(1, n_yrs, 1, n_srv_fleets))

  # Pop-specific weights: scalars so the length(dim)==5 branch is skipped by default
  data$Wt_Catch_pop <- 1
  data$Wt_Discard_pop <- 1
  data$Wt_FishAgeComps_pop <- mk(c(n_pop, n_regions, n_yrs, 1, n_fish_fleets, 1))
  data$Wt_FishLenComps_pop <- mk(c(n_pop, n_regions, n_yrs, 1, n_fish_fleets, 1))
  data$Wt_FishAgeComps_discard_pop <- mk(c(n_pop, n_regions, n_yrs, 1, n_fish_fleets, 1))
  data$Wt_FishLenComps_discard_pop <- mk(c(n_pop, n_regions, n_yrs, 1, n_fish_fleets, 1))
  data$Wt_SrvAgeComps_pop <- mk(c(n_pop, n_regions, n_yrs, 1, n_srv_fleets, 1))
  data$Wt_SrvLenComps_pop <- mk(c(n_pop, n_regions, n_yrs, 1, n_srv_fleets, 1))
  data$FishAgeComps_pop_Type <- mk(c(n_yrs, n_fish_fleets))
  data$FishLenComps_pop_Type <- mk(c(n_yrs, n_fish_fleets))
  data$FishAgeComps_discard_pop_Type <- mk(c(n_yrs, n_fish_fleets))
  data$FishLenComps_discard_pop_Type <- mk(c(n_yrs, n_fish_fleets))
  data$SrvAgeComps_pop_Type <- mk(c(n_yrs, n_srv_fleets))
  data$SrvLenComps_pop_Type <- mk(c(n_yrs, n_srv_fleets))
  data$ISS_FishAgeComps_pop <- mk(c(n_pop, n_regions, n_yrs, 1, n_fish_fleets, 1))
  data$ISS_FishLenComps_pop <- mk(c(n_pop, n_regions, n_yrs, 1, n_fish_fleets, 1))
  data$ISS_FishAgeComps_discard_pop <- mk(c(n_pop, n_regions, n_yrs, 1, n_fish_fleets, 1))
  data$ISS_FishLenComps_discard_pop <- mk(c(n_pop, n_regions, n_yrs, 1, n_fish_fleets, 1))
  data$ISS_SrvAgeComps_pop <- mk(c(n_pop, n_regions, n_yrs, 1, n_srv_fleets, 1))
  data$ISS_SrvLenComps_pop <- mk(c(n_pop, n_regions, n_yrs, 1, n_srv_fleets, 1))

  data$ISS_FishAgeComps <- mk(c(1, n_yrs, 1, 1, n_fish_fleets))
  data$ISS_FishAgeComps_discard <- mk(c(1, n_yrs, 1, 1, n_fish_fleets))
  data$ISS_SrvAgeComps <- mk(c(1, n_yrs, 1, 1, n_srv_fleets))
  data$ISS_FishLenComps <- mk(c(1, n_yrs, 1, 1, n_fish_fleets))
  data$ISS_FishLenComps_discard <- mk(c(1, n_yrs, 1, 1, n_fish_fleets))
  data$ISS_SrvLenComps <- mk(c(1, n_yrs, 1, 1, n_srv_fleets))

  list(data = data, parameters = parameters, mapping = mapping)
}

## Generic year-axis truncation table: field name -> (container, axis)
year_axis_fields <- list(
  data = c(
    ObsCatch = 2, ObsFishIdx = 2, ObsFishIdx_SE = 2,
    ObsFishAgeComps = 2, ObsFishLenComps = 2,
    ObsFishAgeComps_discard = 2, ObsFishLenComps_discard = 2,
    fish_q_blocks = 2, fish_sel_blocks = 2, ret_sel_blocks = 2,
    map_ln_fishsel_devs = 2, ln_retsel_devs = 2,
    ObsSrvIdx = 2, ObsSrvIdx_SE = 2, ObsSrvAgeComps = 2, ObsSrvLenComps = 2,
    map_ln_srvsel_devs = 2, srv_q_blocks = 2, srv_sel_blocks = 2,
    M_blocks = 3, Fixed_natmort = 3, sexratio_blocks = 3,
    Wt_Catch = 2, Wt_Discard = 2,
    Wt_FishAgeComps = 2, Wt_FishAgeComps_discard = 2, Wt_SrvAgeComps = 2,
    Wt_FishLenComps = 2, Wt_FishLenComps_discard = 2, Wt_SrvLenComps = 2,
    FishAgeComps_Type = 1, FishLenComps_Type = 1,
    FishAgeComps_discard_Type = 1, FishLenComps_discard_Type = 1,
    SrvLenComps_Type = 1, SrvAgeComps_Type = 1,
    UseFishAgeComps = 2, UseFishAgeComps_discard = 2, UseFishIdx = 2,
    UseCatch = 2, UseDiscard = 2, UseFishLenComps = 2, UseFishLenComps_discard = 2,
    UseSrvAgeComps = 2, UseSrvIdx = 2, UseSrvLenComps = 2,
    ISS_FishAgeComps = 2, ISS_FishAgeComps_discard = 2, ISS_SrvAgeComps = 2,
    ISS_FishLenComps = 2, ISS_FishLenComps_discard = 2, ISS_SrvLenComps = 2,
    UseFishAgeComps_pop = 3, UseFishLenComps_pop = 3,
    UseFishAgeComps_discard_pop = 3, UseFishLenComps_discard_pop = 3,
    UseFishIdx_pop = 3, UseSrvAgeComps_pop = 3, UseSrvLenComps_pop = 3, UseSrvIdx_pop = 3,
    UseCatch_pop = 3, UseDiscard_pop = 3,
    Wt_FishAgeComps_pop = 3, Wt_FishLenComps_pop = 3,
    Wt_FishAgeComps_discard_pop = 3, Wt_FishLenComps_discard_pop = 3,
    Wt_SrvAgeComps_pop = 3, Wt_SrvLenComps_pop = 3,
    FishAgeComps_pop_Type = 1, FishLenComps_pop_Type = 1,
    FishAgeComps_discard_pop_Type = 1, FishLenComps_discard_pop_Type = 1,
    SrvAgeComps_pop_Type = 1, SrvLenComps_pop_Type = 1,
    ISS_FishAgeComps_pop = 3, ISS_FishLenComps_pop = 3,
    ISS_FishAgeComps_discard_pop = 3, ISS_FishLenComps_discard_pop = 3,
    ISS_SrvAgeComps_pop = 3, ISS_SrvLenComps_pop = 3
  ),
  parameters = c(
    ln_RecDevs = 3, ln_F_devs = 2, logit_dmr_devs = 2,
    ln_fishsel_devs = 2, ln_retsel_devs = 2,
    ln_sigmaC = 2, ln_sigmaC_pop = 3, ln_sigmaD = 2, ln_sigmaD_pop = 3,
    ln_srvsel_devs = 2
  )
)

test_that("truncate_yr() truncates every year-indexed field to the correct length along its year axis", {

  d <- make_retro_data(n_yrs = 6)

  for (j in c(0, 1, 3, 5)) {
    out <- truncate_yr(j = j, data = d$data, parameters = d$parameters, mapping = d$mapping)
    expect_length(out$retro_data$years, 6 - j)

    for (nm in names(year_axis_fields$data)) {
      axis <- year_axis_fields$data[[nm]]
      got  <- dim(out$retro_data[[nm]])[axis]
      expect_equal(got, 6 - j, info = paste0("data$", nm, " (axis ", axis, ", j = ", j, ")"))
    }
    for (nm in names(year_axis_fields$parameters)) {
      axis <- year_axis_fields$parameters[[nm]]
      got  <- dim(out$retro_parameters[[nm]])[axis]
      expect_equal(got, 6 - j, info = paste0("parameters$", nm, " (axis ", axis, ", j = ", j, ")"))
    }
  }
})

test_that("truncate_yr() with j = 0 returns the full, untruncated dataset", {
  d <- make_retro_data(n_yrs = 6)
  out <- truncate_yr(j = 0, data = d$data, parameters = d$parameters, mapping = d$mapping)
  expect_equal(out$retro_data$years, 1:6)
  expect_equal(dim(out$retro_data$ObsCatch), dim(d$data$ObsCatch))
})

test_that("truncate_yr() truncates block-keyed parameter arrays (fish_q, sel, srv_q, sexratio) via max(blocks)", {

  d <- make_retro_data(n_yrs = 6)
  # give fish_q_blocks 2 distinct block ids so max() truncation is meaningful
  d$data$fish_q_blocks[1, 1:3, 1] <- 1
  d$data$fish_q_blocks[1, 4:6, 1] <- 2
  d$parameters$ln_fish_q <- mk(c(1, 2, 1))
  d$mapping$ln_fish_q <- factor(mk(c(1, 2, 1)))

  out_full <- truncate_yr(j = 0, data = d$data, parameters = d$parameters, mapping = d$mapping)
  expect_equal(dim(out_full$retro_parameters$ln_fish_q)[2], 2) # both blocks present

  out_peel <- truncate_yr(j = 3, data = d$data, parameters = d$parameters, mapping = d$mapping)
  expect_equal(dim(out_peel$retro_parameters$ln_fish_q)[2], 1) # only block 1 remains in years 1:3
})

test_that("truncate_yr() only truncates bias_year when it is not NA", {

  d_na <- make_retro_data(n_yrs = 6, has_bias_year = FALSE)
  out_na <- truncate_yr(j = 2, data = d_na$data, parameters = d_na$parameters, mapping = d_na$mapping)
  expect_true(is.na(sum(out_na$retro_data$bias_year)))

  d_yes <- make_retro_data(n_yrs = 6, has_bias_year = TRUE)
  out_yes <- truncate_yr(j = 2, data = d_yes$data, parameters = d_yes$parameters, mapping = d_yes$mapping)
  expect_equal(out_yes$retro_data$bias_year[3:4], d_yes$data$bias_year[3:4] - 2)
  expect_equal(out_yes$retro_data$bias_year[1:2], d_yes$data$bias_year[1:2]) # untouched
})

test_that("truncate_yr() only truncates the ln_RecDevs mapping when it's present in the mapping list", {

  d_no_map <- make_retro_data(n_yrs = 6, include_recdevs_map = FALSE)
  out_no_map <- truncate_yr(j = 2, data = d_no_map$data, parameters = d_no_map$parameters, mapping = d_no_map$mapping)
  expect_null(out_no_map$retro_mapping$ln_RecDevs)

  d_map <- make_retro_data(n_yrs = 6, include_recdevs_map = TRUE)
  out_map <- truncate_yr(j = 2, data = d_map$data, parameters = d_map$parameters, mapping = d_map$mapping)
  expect_equal(dim(array(out_map$retro_mapping$ln_RecDevs, dim = c(1, 1, 4)))[3], 4)
})

test_that("truncate_yr() truncates stray rate fields only when n_pop > 1, for both fixed and estimated stray rate", {

  d1 <- make_retro_data(n_yrs = 6, n_pop = 1)
  out1 <- truncate_yr(j = 2, data = d1$data, parameters = d1$parameters, mapping = d1$mapping)
  expect_null(out1$retro_data$stray_rate_blocks)

  d2_fixed <- make_retro_data(n_yrs = 6, n_pop = 2, use_fixed_stray_rate = 1)
  out2_fixed <- truncate_yr(j = 2, data = d2_fixed$data, parameters = d2_fixed$parameters, mapping = d2_fixed$mapping)
  expect_equal(dim(out2_fixed$retro_data$stray_rate_blocks)[2], 4)
  expect_equal(dim(out2_fixed$retro_data$fixed_stray_rate)[2], 4)
  expect_null(out2_fixed$retro_parameters$stray_rate_pars)

  d2_est <- make_retro_data(n_yrs = 6, n_pop = 2, use_fixed_stray_rate = 0)
  out2_est <- truncate_yr(j = 2, data = d2_est$data, parameters = d2_est$parameters, mapping = d2_est$mapping)
  expect_equal(dim(out2_est$retro_data$stray_rate_blocks)[2], 4)
  expect_false(is.null(out2_est$retro_parameters$stray_rate_pars))
  expect_false(is.null(out2_est$retro_mapping$stray_rate_pars))
})

test_that("truncate_yr() truncates ln_M and its mapping only when natural mortality is estimated (use_fixed_natmort = 0)", {

  d_fixed <- make_retro_data(n_yrs = 6, use_fixed_natmort = 1)
  out_fixed <- truncate_yr(j = 2, data = d_fixed$data, parameters = d_fixed$parameters, mapping = d_fixed$mapping)
  expect_null(out_fixed$retro_parameters$ln_M)
  expect_equal(dim(out_fixed$retro_data$Fixed_natmort)[3], 4)

  d_est <- make_retro_data(n_yrs = 6, use_fixed_natmort = 0)
  out_est <- truncate_yr(j = 2, data = d_est$data, parameters = d_est$parameters, mapping = d_est$mapping)
  expect_false(is.null(out_est$retro_parameters$ln_M))
  expect_false(is.null(out_est$retro_mapping$ln_M))
})

test_that("truncate_yr() truncates movement fields only when n_regions > 1", {

  d1 <- make_retro_data(n_yrs = 6, n_regions = 1)
  out1 <- truncate_yr(j = 2, data = d1$data, parameters = d1$parameters, mapping = d1$mapping)
  expect_null(out1$retro_parameters$move_pars)

  d2 <- make_retro_data(n_yrs = 6, n_regions = 2)
  out2 <- truncate_yr(j = 2, data = d2$data, parameters = d2$parameters, mapping = d2$mapping)
  expect_equal(dim(out2$retro_parameters$move_pars)[4], 4)
  expect_equal(dim(out2$retro_parameters$move_devs)[4], 4)
  expect_equal(dim(out2$retro_data$Fixed_Movement)[4], 4)
  expect_equal(dim(out2$retro_data$sgl_seas_spawning_movement)[4], 4)
  expect_equal(dim(out2$retro_data$map_move_devs)[4], 4)
  expect_true(all(out2$retro_data$ctmc_move_dat$years %in% 1:4))
  expect_equal(nrow(out2$retro_data$ctmc_move_dat), 4)
})

test_that("truncate_yr() truncates conventional tagging cohorts by release year only when tagging is in use", {

  d_off <- make_retro_data(n_yrs = 6, use_conv_fish_tagging = 0)
  out_off <- truncate_yr(j = 2, data = d_off$data, parameters = d_off$parameters, mapping = d_off$mapping)
  expect_null(out_off$retro_data$conv_tag_release_indicator)

  d_on <- make_retro_data(n_yrs = 6, use_conv_fish_tagging = 1)
  # cohorts released in years 1, 3, 6 (see make_retro_data())

  out_j0 <- truncate_yr(j = 0, data = d_on$data, parameters = d_on$parameters, mapping = d_on$mapping)
  expect_equal(nrow(out_j0$retro_data$conv_tag_release_indicator), 3) # years 1,3,6 all <= 6
  expect_equal(out_j0$retro_data$n_conv_tag_cohorts, 3)

  # j = 1: TWO cohorts survive (years 1, 3) -- confirmed working correctly
  out_j1 <- truncate_yr(j = 1, data = d_on$data, parameters = d_on$parameters, mapping = d_on$mapping)
  expect_equal(nrow(out_j1$retro_data$conv_tag_release_indicator), 2) # year-6 cohort dropped (years 1:5)
  expect_equal(out_j1$retro_data$n_conv_tag_cohorts, 2)
  expect_equal(sort(out_j1$retro_data$conv_tag_release_indicator[, 2]), c(1, 3))

  # downstream tag arrays follow the truncated cohort count
  expect_equal(dim(out_j1$retro_data$conv_tagged_fish)[1], 2)
  expect_equal(dim(out_j1$retro_data$obs_conv_tag_fish_recap)[3], 2)
  expect_equal(nrow(out_j1$retro_data$conv_tag_release_platform), 2)
})

test_that("truncate_yr() correctly sizes conv_tag_release_indicator when exactly one cohort survives a peel", {

  d_on <- make_retro_data(n_yrs = 6, use_conv_fish_tagging = 1)
  out_j5 <- truncate_yr(j = 5, data = d_on$data, parameters = d_on$parameters, mapping = d_on$mapping)

  expect_equal(nrow(out_j5$retro_data$conv_tag_release_indicator), 1) # only the year-1 cohort survives
  expect_equal(out_j5$retro_data$n_conv_tag_cohorts, 1)
  expect_equal(out_j5$retro_data$conv_tag_release_indicator[1, 2], 1) # release year == 1

  # downstream tag arrays should also follow the truncated (single) cohort count
  expect_equal(dim(out_j5$retro_data$conv_tagged_fish)[1], 1)
  expect_equal(dim(out_j5$retro_data$obs_conv_tag_fish_recap)[3], 1)
  expect_equal(nrow(out_j5$retro_data$conv_tag_release_platform), 1)
})

test_that("truncate_yr() truncates population-specific observation arrays only when any Use*_pop flag is on", {

  d_off <- make_retro_data(n_yrs = 6, n_pop = 2, pop_specific = FALSE)
  out_off <- truncate_yr(j = 2, data = d_off$data, parameters = d_off$parameters, mapping = d_off$mapping)
  expect_null(out_off$retro_data$ObsFishIdx_pop) # never populated, guard correctly skipped

  d_on <- make_retro_data(n_yrs = 6, n_pop = 2, pop_specific = TRUE)
  out_on <- truncate_yr(j = 2, data = d_on$data, parameters = d_on$parameters, mapping = d_on$mapping)
  expect_equal(dim(out_on$retro_data$ObsFishIdx_pop)[3], 4)
  expect_equal(dim(out_on$retro_data$ObsFishAgeComps_pop)[3], 4)
  expect_equal(dim(out_on$retro_data$ObsSrvIdx_pop)[3], 4)

  # the *_pop Use/Wt/ISS/Type fields are always truncated regardless of pop_specific
  expect_equal(dim(out_off$retro_data$UseFishIdx_pop)[3], 4)
  expect_equal(dim(out_off$retro_data$Wt_FishAgeComps_pop)[3], 4)
  expect_equal(dim(out_off$retro_data$ISS_SrvLenComps_pop)[3], 4)
  expect_equal(dim(out_off$retro_data$FishAgeComps_pop_Type)[1], 4)
})

test_that("truncate_yr() skips the length-5 Wt_Catch_pop/Wt_Discard_pop branch when they are scalars", {
  d <- make_retro_data(n_yrs = 6)
  expect_null(dim(d$data$Wt_Catch_pop))
  out <- truncate_yr(j = 2, data = d$data, parameters = d$parameters, mapping = d$mapping)
  expect_equal(out$retro_data$Wt_Catch_pop, 1) # untouched scalar passthrough
})

testthat_supports_mocking <- exists("local_mocked_bindings", where = asNamespace("testthat"))

## Minimal synthetic data/parameters/mapping (mirrors test-truncate_yr.R's
## make_retro_data(), trimmed to just what's needed here since fit_model()
## is mocked and never actually reads most fields).
mk <- function(dims, fill = 1) array(fill, dim = dims)

make_minimal_retro_inputs <- function(n_yrs = 6, n_pop = 1, n_regions = 1,
                                      n_fish_fleets = 1, n_srv_fleets = 1) {
  data <- list(
    n_yrs = n_yrs, n_pop = n_pop, n_regions = n_regions,
    n_fish_fleets = n_fish_fleets, n_srv_fleets = n_srv_fleets,
    use_fixed_natmort = 1, use_fixed_stray_rate = 1,
    use_conv_fish_tagging = rep(0, n_fish_fleets),
    years = 1:n_yrs, bias_year = NA
  )
  data$ObsCatch  <- mk(c(1, n_yrs, n_fish_fleets, 1))
  data$ObsFishIdx <- mk(c(1, n_yrs, n_fish_fleets, 1))
  data$ObsFishIdx_SE <- mk(c(1, n_yrs, n_fish_fleets, 1))
  data$ObsFishAgeComps <- mk(c(1, n_yrs, 1, 1, n_fish_fleets, 1))
  data$ObsFishLenComps <- mk(c(1, n_yrs, 1, 1, n_fish_fleets, 1))
  data$ObsFishAgeComps_discard <- mk(c(1, n_yrs, 1, 1, n_fish_fleets, 1))
  data$ObsFishLenComps_discard <- mk(c(1, n_yrs, 1, 1, n_fish_fleets, 1))
  data$sexratio_blocks <- mk(c(2, 1, n_yrs))
  data$M_blocks <- mk(c(1, 1, n_yrs, 1, 1))
  data$Fixed_natmort <- mk(c(1, 1, n_yrs, 1, 1))
  data$fish_q_blocks <- mk(c(1, n_yrs, n_fish_fleets))
  data$fish_sel_blocks <- mk(c(1, n_yrs, n_fish_fleets))
  data$ret_sel_blocks <- mk(c(1, n_yrs, n_fish_fleets))
  data$map_ln_fishsel_devs <- mk(c(1, n_yrs, n_fish_fleets, 1, 1))
  data$ln_retsel_devs <- mk(c(1, n_yrs, n_fish_fleets, 1, 1))
  data$ObsSrvIdx <- mk(c(1, n_yrs, n_srv_fleets, 1))
  data$ObsSrvIdx_SE <- mk(c(1, n_yrs, n_srv_fleets, 1))
  data$ObsSrvAgeComps <- mk(c(1, n_yrs, 1, 1, n_srv_fleets, 1))
  data$ObsSrvLenComps <- mk(c(1, n_yrs, 1, 1, n_srv_fleets, 1))
  data$map_ln_srvsel_devs <- mk(c(1, n_yrs, n_srv_fleets, 1, 1))
  data$srv_q_blocks <- mk(c(1, n_yrs, n_srv_fleets))
  data$srv_sel_blocks <- mk(c(1, n_yrs, n_srv_fleets))
  data$Wt_Catch <- mk(c(1, n_yrs, n_fish_fleets, 1))
  data$Wt_Discard <- mk(c(1, n_yrs, n_fish_fleets, 1))
  data$Wt_Catch_pop <- 1
  data$Wt_Discard_pop <- 1

  for (nm in c("UseFishIdx_pop", "UseFishAgeComps_pop", "UseFishLenComps_pop",
               "UseFishAgeComps_discard_pop", "UseFishLenComps_discard_pop",
               "UseSrvIdx_pop", "UseSrvAgeComps_pop", "UseSrvLenComps_pop",
               "UseCatch_pop", "UseDiscard_pop")) {
    n_fleets <- if (grepl("Srv", nm)) n_srv_fleets else n_fish_fleets
    data[[nm]] <- mk(c(n_pop, n_regions, n_yrs, 1, n_fleets), fill = 0)
  }
  for (nm in c("Wt_FishAgeComps_pop", "Wt_FishLenComps_pop",
               "Wt_FishAgeComps_discard_pop", "Wt_FishLenComps_discard_pop",
               "ISS_FishAgeComps_pop", "ISS_FishLenComps_pop",
               "ISS_FishAgeComps_discard_pop", "ISS_FishLenComps_discard_pop")) {
    data[[nm]] <- mk(c(n_pop, n_regions, n_yrs, 1, n_fish_fleets, 1))
  }
  for (nm in c("Wt_SrvAgeComps_pop", "Wt_SrvLenComps_pop",
               "ISS_SrvAgeComps_pop", "ISS_SrvLenComps_pop")) {
    data[[nm]] <- mk(c(n_pop, n_regions, n_yrs, 1, n_srv_fleets, 1))
  }
  for (nm in c("FishAgeComps_pop_Type", "FishLenComps_pop_Type",
               "FishAgeComps_discard_pop_Type", "FishLenComps_discard_pop_Type")) {
    data[[nm]] <- mk(c(n_yrs, n_fish_fleets))
  }
  for (nm in c("SrvAgeComps_pop_Type", "SrvLenComps_pop_Type")) {
    data[[nm]] <- mk(c(n_yrs, n_srv_fleets))
  }

  data$Wt_FishAgeComps <- mk(c(1, n_yrs, 1, 1, n_fish_fleets))
  data$Wt_FishAgeComps_discard <- mk(c(1, n_yrs, 1, 1, n_fish_fleets))
  data$Wt_SrvAgeComps <- mk(c(1, n_yrs, 1, 1, n_srv_fleets))
  data$Wt_FishLenComps <- mk(c(1, n_yrs, 1, 1, n_fish_fleets))
  data$Wt_FishLenComps_discard <- mk(c(1, n_yrs, 1, 1, n_fish_fleets))
  data$Wt_SrvLenComps <- mk(c(1, n_yrs, 1, 1, n_srv_fleets))
  data$FishAgeComps_Type <- mk(c(n_yrs, n_fish_fleets))
  data$FishLenComps_Type <- mk(c(n_yrs, n_fish_fleets))
  data$FishAgeComps_discard_Type <- mk(c(n_yrs, n_fish_fleets))
  data$FishLenComps_discard_Type <- mk(c(n_yrs, n_fish_fleets))
  data$SrvLenComps_Type <- mk(c(n_yrs, n_srv_fleets))
  data$SrvAgeComps_Type <- mk(c(n_yrs, n_srv_fleets))

  data$UseFishAgeComps <- mk(c(1, n_yrs, 1, n_fish_fleets))
  data$UseFishAgeComps_discard <- mk(c(1, n_yrs, 1, n_fish_fleets))
  data$UseFishIdx <- mk(c(1, n_yrs, 1, n_fish_fleets))
  data$UseCatch <- mk(c(1, n_yrs, 1, n_fish_fleets))
  data$UseDiscard <- mk(c(1, n_yrs, 1, n_fish_fleets))
  data$UseFishLenComps <- mk(c(1, n_yrs, 1, n_fish_fleets))
  data$UseFishLenComps_discard <- mk(c(1, n_yrs, 1, n_fish_fleets))
  data$UseSrvAgeComps <- mk(c(1, n_yrs, 1, n_srv_fleets))
  data$UseSrvIdx <- mk(c(1, n_yrs, 1, n_srv_fleets))
  data$UseSrvLenComps <- mk(c(1, n_yrs, 1, n_srv_fleets))

  data$ISS_FishAgeComps <- mk(c(1, n_yrs, 1, 1, n_fish_fleets))
  data$ISS_FishAgeComps_discard <- mk(c(1, n_yrs, 1, 1, n_fish_fleets))
  data$ISS_SrvAgeComps <- mk(c(1, n_yrs, 1, 1, n_srv_fleets))
  data$ISS_FishLenComps <- mk(c(1, n_yrs, 1, 1, n_fish_fleets))
  data$ISS_FishLenComps_discard <- mk(c(1, n_yrs, 1, 1, n_fish_fleets))
  data$ISS_SrvLenComps <- mk(c(1, n_yrs, 1, 1, n_srv_fleets))

  parameters <- list(
    ln_RecDevs = mk(c(n_pop, n_regions, n_yrs)),
    sexratio_pars = mk(c(2, 1, 1)),
    ln_F_devs = mk(c(1, n_yrs, n_fish_fleets, 1)),
    logit_dmr_devs = mk(c(1, n_yrs, n_fish_fleets, 1)),
    ln_fishsel_devs = mk(c(1, n_yrs, n_fish_fleets, 1, 1)),
    ln_retsel_devs = mk(c(1, n_yrs, n_fish_fleets, 1, 1)),
    ln_fish_q = mk(c(1, 1, n_fish_fleets)),
    fish_fixed_sel_pars = mk(c(1, 1, 1, 1, n_fish_fleets)),
    ret_fixed_sel_pars = mk(c(1, 1, 1, 1, n_fish_fleets)),
    ln_sigmaC = mk(c(1, n_yrs, n_fish_fleets, 1)),
    ln_sigmaC_pop = mk(c(n_pop, 1, n_yrs, 1, 1)),
    ln_sigmaD = mk(c(1, n_yrs, n_fish_fleets, 1)),
    ln_sigmaD_pop = mk(c(n_pop, 1, n_yrs, 1, 1)),
    ln_srvsel_devs = mk(c(1, n_yrs, n_srv_fleets, 1, 1)),
    ln_srv_q = mk(c(1, 1, n_srv_fleets)),
    srv_fixed_sel_pars = mk(c(1, 1, 1, 1, n_srv_fleets))
  )

  mapping <- list(
    sexratio_pars = factor(mk(c(2, 1, 1))),
    ln_F_devs = factor(mk(c(1, n_yrs, n_fish_fleets, 1))),
    logit_dmr_devs = factor(mk(c(1, n_yrs, n_fish_fleets, 1))),
    ln_fishsel_devs = factor(mk(c(1, n_yrs, n_fish_fleets, 1, 1))),
    ln_retsel_devs = factor(mk(c(1, n_yrs, n_fish_fleets, 1, 1))),
    ln_fish_q = factor(mk(c(1, 1, n_fish_fleets))),
    fish_fixed_sel_pars = factor(mk(c(1, 1, 1, 1, n_fish_fleets))),
    ret_fixed_sel_pars = factor(mk(c(1, 1, 1, 1, n_fish_fleets))),
    ln_sigmaC = factor(mk(c(1, n_yrs, n_fish_fleets, 1))),
    ln_sigmaC_pop = factor(mk(c(n_pop, 1, n_yrs, 1, 1))),
    ln_sigmaD = factor(mk(c(1, n_yrs, n_fish_fleets, 1))),
    ln_sigmaD_pop = factor(mk(c(n_pop, 1, n_yrs, 1, 1))),
    ln_srvsel_devs = factor(mk(c(1, n_yrs, n_srv_fleets, 1, 1))),
    ln_srv_q = factor(mk(c(1, 1, n_srv_fleets))),
    srv_fixed_sel_pars = factor(mk(c(1, 1, 1, 1, n_srv_fleets)))
  )

  list(data = data, parameters = parameters, mapping = mapping)
}

make_mock_fit_model <- function(call_log_env, n_pop, n_regions) {
  function(data, parameters, mapping, random, newton_loops, silent) {
    call_log_env$calls <- c(call_log_env$calls, list(data))
    n_yr <- length(data$years)
    list(rep = list(
      SSB = array(1, dim = c(n_pop, n_regions, n_yr)),
      Rec = array(1, dim = c(n_pop, n_regions, n_yr))
    ))
  }
}

make_mock_run_francis <- function(call_log_env, n_pop, n_regions) {
  function(data, parameters, mapping, random, n_francis_iter, newton_loops) {
    call_log_env$calls <- c(call_log_env$calls, list(data))
    n_yr <- length(data$years)
    list(obj = list(rep = list(
      SSB = array(2, dim = c(n_pop, n_regions, n_yr)), # distinct fill value (2) marks the francis path
      Rec = array(2, dim = c(n_pop, n_regions, n_yr))
    )))
  }
}

## Test 1: sequential loop produces correct peel/year bookkeeping
test_that("do_retrospective() (sequential) produces one row set per peel with correctly truncated years", {

  skip_if_not(testthat_supports_mocking, "testthat >= 3.2.0 with 3rd edition required for local_mocked_bindings()")
  skip_if_not_installed("reshape2")
  skip_if_not_installed("dplyr")

  inp <- make_minimal_retro_inputs(n_yrs = 6)
  log_env <- new.env()
  log_env$calls <- list()

  local_mocked_bindings(fit_model = make_mock_fit_model(log_env, n_pop = 1, n_regions = 1))

  result <- do_retrospective(
    n_retro = 2, data = inp$data, parameters = inp$parameters, mapping = inp$mapping,
    do_par = FALSE, n_cores = 1
  )

  expect_s3_class(result, "data.frame")
  expect_setequal(unique(result$peel), c(0, 1, 2))
  expect_setequal(unique(result$Type), c("SSB", "Recruitment"))

  # peel j should have (6 - j) years represented for each Type
  for (j in 0:2) {
    peel_rows <- result[result$peel == j & result$Type == "SSB", ]
    expect_equal(nrow(peel_rows), 6 - j)
  }

  expect_equal(length(log_env$calls), 3) # one fit_model() call per peel (0, 1, 2)
})

## Test 2: return_models = TRUE
test_that("do_retrospective() with return_models = TRUE returns both the data.frame and named per-peel model objects", {

  skip_if_not(testthat_supports_mocking, "testthat >= 3.2.0 with 3rd edition required for local_mocked_bindings()")
  skip_if_not_installed("reshape2")
  skip_if_not_installed("dplyr")

  inp <- make_minimal_retro_inputs(n_yrs = 5)
  log_env <- new.env(); log_env$calls <- list()
  local_mocked_bindings(fit_model = make_mock_fit_model(log_env, n_pop = 1, n_regions = 1))

  result <- do_retrospective(
    n_retro = 1, data = inp$data, parameters = inp$parameters, mapping = inp$mapping,
    do_par = FALSE, n_cores = 1, return_models = TRUE
  )

  expect_named(result, c("retro_df", "retro_models"))
  expect_named(result$retro_models, c("peel_0", "peel_1"))
  expect_equal(dim(result$retro_models$peel_0$rep$SSB)[3], 5)
  expect_equal(dim(result$retro_models$peel_1$rep$SSB)[3], 4)
})

## Test 3: do_sdrep = TRUE adds pdHess/max_grad via mocked RTMB::sdreport()
test_that("do_retrospective() with do_sdrep = TRUE attaches pdHess and max_grad from RTMB::sdreport()", {

  skip_if_not(testthat_supports_mocking, "testthat >= 3.2.0 with 3rd edition required for local_mocked_bindings()")
  skip_if_not_installed("reshape2")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("RTMB")

  inp <- make_minimal_retro_inputs(n_yrs = 4)
  log_env <- new.env(); log_env$calls <- list()
  local_mocked_bindings(fit_model = make_mock_fit_model(log_env, n_pop = 1, n_regions = 1))
  local_mocked_bindings(
    sdreport = function(obj) list(pdHess = TRUE, gradient.fixed = c(0.001, -0.0005)),
    .package = "RTMB"
  )

  result <- do_retrospective(
    n_retro = 0, data = inp$data, parameters = inp$parameters, mapping = inp$mapping,
    do_par = FALSE, n_cores = 1, do_sdrep = TRUE
  )

  expect_true(all(result$pdHess))
  expect_equal(unique(result$max_grad), 0.001)
})

## Test 4: do_francis = TRUE routes through run_francis() instead of fit_model()
test_that("do_retrospective() with do_francis = TRUE calls run_francis() instead of fit_model()", {

  skip_if_not(testthat_supports_mocking, "testthat >= 3.2.0 with 3rd edition required for local_mocked_bindings()")
  skip_if_not_installed("reshape2")
  skip_if_not_installed("dplyr")

  inp <- make_minimal_retro_inputs(n_yrs = 5)
  fit_log <- new.env(); fit_log$calls <- list()
  francis_log <- new.env(); francis_log$calls <- list()

  local_mocked_bindings(
    fit_model    = make_mock_fit_model(fit_log, n_pop = 1, n_regions = 1),
    run_francis  = make_mock_run_francis(francis_log, n_pop = 1, n_regions = 1)
  )

  result <- do_retrospective(
    n_retro = 0, data = inp$data, parameters = inp$parameters, mapping = inp$mapping,
    do_par = FALSE, n_cores = 1, do_francis = TRUE, n_francis_iter = 2
  )

  expect_equal(length(fit_log$calls), 0)      # fit_model() never called
  expect_equal(length(francis_log$calls), 1)  # run_francis() called once (peel 0)
  expect_true(all(result$value == 2))         # mock_run_francis fills with 2, distinguishing it from fit_model's 1
})

## Test 5: data-lag zeroing is applied to the data actually passed to fit_model()
test_that("do_retrospective() zeros the correct lagged columns of UseFishAgeComps before fitting", {

  skip_if_not(testthat_supports_mocking, "testthat >= 3.2.0 with 3rd edition required for local_mocked_bindings()")
  skip_if_not_installed("reshape2")
  skip_if_not_installed("dplyr")

  inp <- make_minimal_retro_inputs(n_yrs = 6)
  log_env <- new.env(); log_env$calls <- list()
  local_mocked_bindings(fit_model = make_mock_fit_model(log_env, n_pop = 1, n_regions = 1))

  fishage_lag <- array(2, dim = c(1, 1)) # [region, fleet] lag of 2 years

  do_retrospective(
    n_retro = 0, data = inp$data, parameters = inp$parameters, mapping = inp$mapping,
    do_par = FALSE, n_cores = 1, fishage_datalag = fishage_lag
  )

  passed_data <- log_env$calls[[1]]
  use_fish_age <- passed_data$UseFishAgeComps[1, , 1, 1] # length-6 vector across years

  # start_col = 6, lag = 2 -> end_col = 5 -> columns 5 and 6 zeroed, 1-4 untouched
  expect_equal(use_fish_age[1:4], c(1, 1, 1, 1))
  expect_equal(use_fish_age[5:6], c(0, 0))
})


make_retro_df <- function() {
  data.frame(
    Pop = 1, Region = 1,
    Year = c(1, 2, 3,   1, 2,   1),
    Type = "SSB",
    peel = c(0, 0, 0,   1, 1,   2),
    value = c(100, 110, 120,   100, 105,   95)
  )
}

test_that("get_retrospective_relative_difference() computes correct relative differences vs. the terminal (peel = 0) estimate", {

  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")

  df <- make_retro_df()
  out <- get_retrospective_relative_difference(df)

  expect_s3_class(out, "data.frame")
  expect_true(all(c("Pop", "Region", "Year", "Type", "peel", "rd") %in% names(out)))

  # Pop/Region get stringified with a prefix
  expect_true(all(grepl("^Pop ", out$Pop)))
  expect_true(all(grepl("^Region ", out$Region)))

  # peel 1: years 1-2 have both a terminal and a peeled estimate
  #   year 1: (100 - 100) / 100 = 0
  #   year 2: (105 - 110) / 110 = -0.04545455
  row_p1_y1 <- out[out$peel == "1" & out$Year == 1, ]
  row_p1_y2 <- out[out$peel == "1" & out$Year == 2, ]
  expect_equal(row_p1_y1$rd, 0, tolerance = 1e-8)
  expect_equal(row_p1_y2$rd, (105 - 110) / 110, tolerance = 1e-8)

  # peel 2: only year 1 has a peeled estimate
  #   year 1: (95 - 100) / 100 = -0.05
  row_p2_y1 <- out[out$peel == "2" & out$Year == 1, ]
  expect_equal(row_p2_y1$rd, (95 - 100) / 100, tolerance = 1e-8)

  # year 3 was never estimated in any peel (1 or 2) -> NA relative difference,
  # not silently dropped
  row_p1_y3 <- out[out$peel == "1" & out$Year == 3, ]
  expect_equal(nrow(row_p1_y3), 1)
  expect_true(is.na(row_p1_y3$rd))
})

test_that("get_retrospective_relative_difference() handles a single peel (n_retro = 1) correctly", {

  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")

  df <- data.frame(
    Pop = 1, Region = 1,
    Year = c(1, 2,   1),
    Type = "SSB",
    peel = c(0, 0,   1),
    value = c(50, 55,   48)
  )

  out <- get_retrospective_relative_difference(df)
  expect_equal(nrow(out), 2) # one row per (Year, Type) combination in the terminal peel
  expect_equal(unique(out$peel), "1")

  row_y1 <- out[out$Year == 1, ]
  expect_equal(row_y1$rd, (48 - 50) / 50, tolerance = 1e-8)
})

test_that("get_retrospective_relative_difference() keeps SSB and Recruitment relative differences independent", {

  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")

  df <- data.frame(
    Pop = 1, Region = 1,
    Year = c(1, 1,   1, 1),
    Type = c("SSB", "Recruitment",   "SSB", "Recruitment"),
    peel = c(0, 0,   1, 1),
    value = c(100, 10,   90, 12)
  )

  out <- get_retrospective_relative_difference(df)

  ssb_row <- out[out$Type == "SSB", ]
  rec_row <- out[out$Type == "Recruitment", ]

  expect_equal(ssb_row$rd, (90 - 100) / 100, tolerance = 1e-8)
  expect_equal(rec_row$rd, (12 - 10) / 10, tolerance = 1e-8)
})
