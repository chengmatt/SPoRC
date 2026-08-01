library(SPoRC)
library(testthat)

# do_*_mapping helpers call collect_message(), which appends to a
# `messages_list` object via `<<-`. In normal use this is initialized by the
# enclosing Setup_Mod_* wrapper before any mapping helper runs; since these
# tests call the mapping helpers directly, it must be pre-created here.
assign("messages_list", character(0), envir = .GlobalEnv)

# Helper to build a minimal input_list for do_comp_theta_mapping /
# do_comp_corr_pars_mapping (comp_prefix = "FishAge", non-discard).
# n_regions x n_fish_fleets Type matrix, per-fleet LikeType, and a
# [region, year, season, fleet] Use array (or [pop, region, year, season,
# fleet] when has_pop = TRUE).
make_comp_input_list <- function(comp_type_mat, like_type_vec, use_arr,
                                 n_sexes = 2, has_pop = FALSE, n_pop = 1) {

  n_regions <- nrow(comp_type_mat)
  n_fish_fleets <- ncol(comp_type_mat)

  input_list <- list(data = list(), par = list(), map = list())
  input_list$data$n_regions <- n_regions
  input_list$data$n_fish_fleets <- n_fish_fleets
  input_list$data$n_sexes <- n_sexes
  input_list$data$n_pop <- n_pop

  suffix <- if (has_pop) "_pop" else ""
  input_list$data[[paste0("FishAgeComps", suffix, "_Type")]] <- comp_type_mat
  input_list$data[[paste0("FishAgeComps", suffix, "_LikeType")]] <- like_type_vec
  input_list$data[[paste0("UseFishAgeComps", suffix)]] <- use_arr

  theta_dim <- if (has_pop) c(n_pop, n_regions, n_sexes, n_fish_fleets) else c(n_regions, n_sexes, n_fish_fleets)
  theta_agg_dim <- if (has_pop) c(n_pop, n_fish_fleets) else n_fish_fleets
  corr_dim <- if (has_pop) c(n_pop, n_regions, n_sexes, n_fish_fleets, 2) else c(n_regions, n_sexes, n_fish_fleets, 2)
  corr_agg_dim <- theta_agg_dim

  par_prefix <- paste0("FishAge", suffix)
  input_list$par[[paste0("ln_", par_prefix, "_theta")]] <- array(0, dim = theta_dim)
  input_list$par[[paste0("ln_", par_prefix, "_theta_agg")]] <- array(0, dim = theta_agg_dim)
  input_list$par[[paste0(par_prefix, "_corr_pars")]] <- array(0.01, dim = corr_dim)
  input_list$par[[paste0(par_prefix, "_corr_pars_agg")]] <- array(0.01, dim = corr_agg_dim)

  input_list
}

test_that("do_comp_theta_mapping respects the region_has_data guard", {

  # 1 fleet, 2 regions, split-by-sex-and-region (type 1), logistic-normal (LikeType 3).
  # Region 1 has active comps, region 2 does not.
  comp_type_mat <- matrix(1, nrow = 2, ncol = 1)
  like_type_vec <- 3
  use_arr <- array(0, dim = c(2, 1, 1, 1)) # [region, year, season, fleet]
  use_arr[1, , , ] <- 1 # region 1 has data; region 2 does not

  il <- make_comp_input_list(comp_type_mat, like_type_vec, use_arr, n_sexes = 2)
  out <- SPoRC:::do_comp_theta_mapping(il, comp_prefix = "FishAge")

  map_theta <- array(as.integer(out$map$ln_FishAge_theta), dim = c(2, 2, 1))
  expect_false(any(is.na(map_theta[1, , ]))) # region 1 (has data): both sexes estimated
  expect_true(all(is.na(map_theta[2, , ])))  # region 2 (no data): fixed at NA
  expect_true(all(is.na(out$map$ln_FishAge_theta_agg))) # type 0 (agg) never requested
})

test_that("do_comp_theta_mapping activates the aggregated parameter for type 0", {

  comp_type_mat <- matrix(0, nrow = 1, ncol = 1)
  like_type_vec <- 3
  use_arr <- array(1, dim = c(1, 1, 1, 1))

  il <- make_comp_input_list(comp_type_mat, like_type_vec, use_arr, n_sexes = 1)
  out <- SPoRC:::do_comp_theta_mapping(il, comp_prefix = "FishAge")

  expect_false(is.na(out$map$ln_FishAge_theta_agg))
  expect_true(all(is.na(out$map$ln_FishAge_theta))) # type 0 only, not 1 or 2
})

test_that("do_comp_theta_mapping maps everything to NA for a multinomial (LikeType 0) fleet", {

  comp_type_mat <- matrix(1, nrow = 2, ncol = 1)
  like_type_vec <- 0 # multinomial
  use_arr <- array(1, dim = c(2, 1, 1, 1))

  il <- make_comp_input_list(comp_type_mat, like_type_vec, use_arr, n_sexes = 2)
  out <- SPoRC:::do_comp_theta_mapping(il, comp_prefix = "FishAge")

  expect_true(all(is.na(out$map$ln_FishAge_theta)))
  expect_true(all(is.na(out$map$ln_FishAge_theta_agg)))
})

test_that("do_comp_theta_mapping: has_pop = TRUE applies the guard per (pop, region)", {

  comp_type_mat <- matrix(1, nrow = 2, ncol = 1)
  like_type_vec <- 3
  # [pop, region, year, season, fleet]; pop 1 has data in region 1 only,
  # pop 2 has data in region 2 only.
  use_arr <- array(0, dim = c(2, 2, 1, 1, 1))
  use_arr[1, 1, , , ] <- 1
  use_arr[2, 2, , , ] <- 1

  il <- make_comp_input_list(comp_type_mat, like_type_vec, use_arr, n_sexes = 1, has_pop = TRUE, n_pop = 2)
  out <- SPoRC:::do_comp_theta_mapping(il, comp_prefix = "FishAge", has_pop = TRUE)

  map_theta <- array(as.integer(out$map$ln_FishAge_pop_theta), dim = c(2, 2, 1, 1))
  expect_false(is.na(map_theta[1, 1, , ])) # pop 1, region 1: has data
  expect_true(is.na(map_theta[1, 2, , ]))  # pop 1, region 2: no data
  expect_true(is.na(map_theta[2, 1, , ]))  # pop 2, region 1: no data
  expect_false(is.na(map_theta[2, 2, , ])) # pop 2, region 2: has data
})

test_that("do_comp_corr_pars_mapping activates one element for LikeType 3 and two for LikeType 4", {

  comp_type_mat <- matrix(2, nrow = 1, ncol = 1) # joint by sex, split by region
  use_arr <- array(1, dim = c(1, 1, 1, 1))

  il_1d <- make_comp_input_list(comp_type_mat, like_type_vec = 3, use_arr, n_sexes = 2)
  out_1d <- SPoRC:::do_comp_corr_pars_mapping(il_1d, comp_prefix = "FishAge")
  corr_1d <- array(as.integer(out_1d$map$FishAge_corr_pars), dim = c(1, 2, 1, 2))
  expect_false(is.na(corr_1d[1, 1, 1, 1])) # AR1 coefficient active
  expect_true(is.na(corr_1d[1, 1, 1, 2]))  # sex-correlation slot unused for LikeType 3

  il_2d <- make_comp_input_list(comp_type_mat, like_type_vec = 4, use_arr, n_sexes = 2)
  out_2d <- SPoRC:::do_comp_corr_pars_mapping(il_2d, comp_prefix = "FishAge")
  corr_2d <- array(as.integer(out_2d$map$FishAge_corr_pars), dim = c(1, 2, 1, 2))
  expect_false(is.na(corr_2d[1, 1, 1, 1])) # AR1 coefficient
  expect_false(is.na(corr_2d[1, 1, 1, 2])) # sex correlation, both active for LikeType 4
})

test_that("do_comp_corr_pars_mapping maps everything to NA when there is no data", {

  comp_type_mat <- matrix(1, nrow = 1, ncol = 1)
  like_type_vec <- 3
  use_arr <- array(0, dim = c(1, 1, 1, 1)) # no active comps anywhere

  il <- make_comp_input_list(comp_type_mat, like_type_vec, use_arr, n_sexes = 1)
  out <- SPoRC:::do_comp_corr_pars_mapping(il, comp_prefix = "FishAge")

  expect_true(all(is.na(out$map$FishAge_corr_pars)))
  expect_true(all(is.na(out$map$FishAge_corr_pars_agg)))
})

test_that("do_sigmaC_pop_mapping / do_sigmaD_pop_mapping correctly share across dimensions (regression test)", {

  # Regression test for a bug where the hand-rolled pop-variant mapping used
  # `grepl("p", spec) && p > 1` to decide whether to broadcast a shared
  # dimension -- since this is always FALSE at the exact index (1) where the
  # broadcast write happens, sharing over ANY dimension silently left every
  # cell past the first index as NA (fixed) instead of tied to a shared
  # parameter. Both functions were rewritten to delegate to
  # build_shared_spec_map(), which does not have this bug.

  n_pop <- 2; n_regions <- 2; n_years <- 1; n_seas <- 1; n_fish_fleets <- 1
  il <- list(
    data = list(n_pop = n_pop, n_regions = n_regions, years = 1:n_years,
               n_seas = n_seas, n_fish_fleets = n_fish_fleets),
    par = list(
      ln_sigmaC_pop = array(0, dim = c(n_pop, n_regions, n_years, n_seas, n_fish_fleets)),
      ln_sigmaD_pop = array(0, dim = c(n_pop, n_regions, n_years, n_seas, n_fish_fleets))
    ),
    map = list()
  )

  out_C <- SPoRC:::do_sigmaC_pop_mapping(il, sigmaC_pop_spec = "est_shared_pop")
  map_C <- array(as.integer(out_C$map$ln_sigmaC_pop), dim = c(n_pop, n_regions, n_years, n_seas, n_fish_fleets))
  expect_false(anyNA(map_C)) # sharing must not leave cells fixed at NA
  expect_equal(map_C[1, 1, 1, 1, 1], map_C[2, 1, 1, 1, 1]) # shared across pop...
  expect_false(map_C[1, 1, 1, 1, 1] == map_C[1, 2, 1, 1, 1]) # ...but unique per region

  out_D <- SPoRC:::do_sigmaD_pop_mapping(il, sigmaD_pop_spec = "est_shared_pop_r")
  map_D <- array(as.integer(out_D$map$ln_sigmaD_pop), dim = c(n_pop, n_regions, n_years, n_seas, n_fish_fleets))
  expect_false(anyNA(map_D))
  expect_equal(length(unique(as.vector(map_D))), 1) # shared across both pop and region -> single parameter
})

test_that("do_comp_theta_mapping / do_comp_corr_pars_mapping also serve survey comps via fleet_field", {

  # Same helpers now back setup_survey_comps.R's SrvAge/SrvLen mapping too, selected
  # via comp_prefix + fleet_field = "n_srv_fleets" instead of the fishery
  # default "n_fish_fleets". Region 1 has data, region 2 does not -- this also
  # regression-tests that the region_has_data guard (previously present only
  # for the non-pop SrvAge variant) is now applied consistently to SrvLen and
  # the pop variants too.
  n_regions <- 2; n_srv_fleets <- 1; n_sexes <- 2
  comp_type_mat <- matrix(1, nrow = n_regions, ncol = n_srv_fleets)
  use_arr <- array(0, dim = c(n_regions, 1, 1, n_srv_fleets)) # [region, year, season, fleet]
  use_arr[1, , , ] <- 1

  il <- list(data = list(), par = list(), map = list())
  il$data$n_regions <- n_regions
  il$data$n_srv_fleets <- n_srv_fleets
  il$data$n_sexes <- n_sexes
  il$data$SrvLenComps_Type <- comp_type_mat
  il$data$SrvLenComps_LikeType <- 3
  il$data$UseSrvLenComps <- use_arr
  il$par$ln_SrvLen_theta <- array(0, dim = c(n_regions, n_sexes, n_srv_fleets))
  il$par$ln_SrvLen_theta_agg <- array(0, dim = n_srv_fleets)

  out <- SPoRC:::do_comp_theta_mapping(il, comp_prefix = "SrvLen", fleet_field = "n_srv_fleets")
  map_theta <- array(as.integer(out$map$ln_SrvLen_theta), dim = c(n_regions, n_sexes, n_srv_fleets))
  expect_false(any(is.na(map_theta[1, , ]))) # region 1 (has data): estimated
  expect_true(all(is.na(map_theta[2, , ])))  # region 2 (no data): fixed at NA
})

test_that("do_q_mapping estimates q when only population-specific index data is used (regression test)", {

  # Regression test for a bug where `sum(UseSrvIdx_pop[,r,,,f] == 0)` (misplaced
  # parenthesis, counting zero-cells) was used instead of `sum(UseSrvIdx_pop[,r,,,f]) == 0`
  # (checking whether the sum itself is zero). With any mix of used/unused
  # populations, the buggy version was truthy far too often, incorrectly fixing
  # q at NA even when population-specific data existed. Same class of bug as the
  # sigmaC_pop/sigmaD_pop one above, found while merging do_fish_q_mapping and
  # do_srv_q_mapping into do_q_mapping.

  n_regions <- 1; n_years <- 1; n_seas <- 1; n_fish_fleets <- 1; n_pop <- 2

  il <- list(
    data = list(
      n_regions = n_regions, n_fish_fleets = n_fish_fleets,
      UseFishIdx = array(0, dim = c(n_regions, n_years, n_seas, n_fish_fleets)), # aggregate index NOT used
      UseFishIdx_pop = array(c(1, 0), dim = c(n_pop, n_regions, n_years, n_seas, n_fish_fleets)), # pop 1 uses it, pop 2 doesn't
      fish_q_blocks = array("none_Fleet_1", dim = c(n_regions, n_years, n_fish_fleets))
    ),
    par = list(ln_fish_q = array(0, dim = c(n_regions, 1, n_fish_fleets))),
    map = list()
  )

  out <- SPoRC:::do_q_mapping(il, q_spec = "est_all", prefix = "fish", fleet_field = "n_fish_fleets", fleet_label = "fishery fleet")
  expect_false(is.na(out$map$ln_fish_q[1])) # must be estimated, not fixed, since pop 1 uses this index
})

test_that("do_q_mapping fixes q when no index data (aggregate or pop-specific) is used at all", {

  n_regions <- 1; n_years <- 1; n_seas <- 1; n_fish_fleets <- 1; n_pop <- 2

  il <- list(
    data = list(
      n_regions = n_regions, n_fish_fleets = n_fish_fleets,
      UseFishIdx = array(0, dim = c(n_regions, n_years, n_seas, n_fish_fleets)),
      UseFishIdx_pop = array(0, dim = c(n_pop, n_regions, n_years, n_seas, n_fish_fleets)),
      fish_q_blocks = array("none_Fleet_1", dim = c(n_regions, n_years, n_fish_fleets))
    ),
    par = list(ln_fish_q = array(0, dim = c(n_regions, 1, n_fish_fleets))),
    map = list()
  )

  out <- SPoRC:::do_q_mapping(il, q_spec = "est_all", prefix = "fish", fleet_field = "n_fish_fleets", fleet_label = "fishery fleet")
  expect_true(is.na(out$map$ln_fish_q[1]))
})

test_that("do_q_mapping also serves survey catchability via prefix = 'srv'", {

  n_regions <- 1; n_years <- 1; n_seas <- 1; n_srv_fleets <- 1; n_pop <- 1

  il <- list(
    data = list(
      n_regions = n_regions, n_srv_fleets = n_srv_fleets,
      UseSrvIdx = array(1, dim = c(n_regions, n_years, n_seas, n_srv_fleets)),
      UseSrvIdx_pop = array(0, dim = c(n_pop, n_regions, n_years, n_seas, n_srv_fleets)),
      srv_q_blocks = array("none_Fleet_1", dim = c(n_regions, n_years, n_srv_fleets))
    ),
    par = list(ln_srv_q = array(0, dim = c(n_regions, 1, n_srv_fleets))),
    map = list()
  )

  out <- SPoRC:::do_q_mapping(il, q_spec = "est_all", prefix = "srv", fleet_field = "n_srv_fleets", fleet_label = "survey fleet")
  expect_false(is.na(out$map$ln_srv_q[1]))
})

test_that("do_sigma_dmr_mapping / do_dmr_mean_mapping / do_sigmaD_mapping match build_shared_spec_map", {

  n_regions <- 2; n_years <- 2; n_seas <- 2; n_fish_fleets <- 2
  il <- list(
    data = list(n_regions = n_regions, years = 1:n_years, n_seas = n_seas, n_fish_fleets = n_fish_fleets),
    par = list(
      ln_sigma_dmr = array(0, dim = c(n_regions, n_seas, n_fish_fleets)),
      logit_dmr_mean = array(0, dim = c(n_regions, n_seas, n_fish_fleets)),
      ln_sigmaD = array(0, dim = c(n_regions, n_years, n_seas, n_fish_fleets))
    ),
    map = list()
  )

  dims_3d <- c(region = n_regions, season = n_seas, fleet = n_fish_fleets)
  abbrev_3d <- c(r = "region", seas = "season", f = "fleet")
  expected_3d <- SPoRC:::build_shared_spec_map(dims_3d, "est_shared_r_f", abbrev_3d)

  out_dmr <- SPoRC:::do_sigma_dmr_mapping(il, sigma_dmr_spec = "est_shared_r_f")
  expect_equal(as.integer(out_dmr$map$ln_sigma_dmr), as.integer(expected_3d))

  out_mean <- SPoRC:::do_dmr_mean_mapping(il, dmr_mean_spec = "est_shared_r_f")
  expect_equal(as.integer(out_mean$map$logit_dmr_mean), as.integer(expected_3d))

  dims_4d <- c(region = n_regions, year = n_years, season = n_seas, fleet = n_fish_fleets)
  abbrev_4d <- c(r = "region", y = "year", seas = "season", f = "fleet")
  expected_4d <- SPoRC:::build_shared_spec_map(dims_4d, "est_shared_r_y_f", abbrev_4d)

  out_D <- SPoRC:::do_sigmaD_mapping(il, sigmaD_spec = "est_shared_r_y_f")
  expect_equal(as.integer(out_D$map$ln_sigmaD), as.integer(expected_4d))
})
