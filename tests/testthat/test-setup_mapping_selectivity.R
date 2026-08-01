library(SPoRC)
library(testthat)

# do_*_mapping helpers call collect_message(), which appends to a
# `messages_list` object via `<<-`. See test-setup_fishery_mapping.R for why
# this must be pre-created when calling the helpers directly (bypassing the
# enclosing Setup_Mod_* wrapper that normally initializes it).
assign("messages_list", character(0), envir = .GlobalEnv)

# Minimal single-region, single-fleet, single-sex, single-block input_list
# with a logistic (model 0, 2 parameters) selectivity form, used across all
# three do_*_mapping smoke tests below. `prefix` selects fish/ret/srv field
# naming; `use_field` mirrors what setup_fishery_selectivity.R/setup_survey_selectivity.R pass in
# (fish/ret -> "Catch", srv -> "SrvIdx").
make_sel_input_list <- function(prefix, use_field, n_bins = 2, n_years = 1) {

  n_regions <- 1; n_fleets <- 1; n_sexes <- 1

  il <- list(data = list(), par = list(), map = list())
  il$data$n_regions <- n_regions
  il$data$n_sexes <- n_sexes
  il$data$years <- seq_len(n_years)
  il$data$ages <- seq_len(n_bins)
  il$data$n_proj_yrs_devs <- 0
  il$data[[paste0("n_", ifelse(prefix == "srv", "srv", "fish"), "_fleets")]] <- n_fleets
  il$data[[paste0("use_fixed_", prefix, "_sel")]] <- 0
  il$data[[paste0(prefix, "_sel_blocks")]] <- array(1, dim = c(n_regions, n_years, n_fleets))
  il$data[[paste0(prefix, "_sel_model")]] <- array(0, dim = c(n_regions, n_years, n_fleets)) # logistic -> 2 pars
  il$data[[paste0("cont_tv_", prefix, "_sel")]] <- array(1, dim = c(n_regions, n_fleets)) # iid time-variation
  il$data[[paste0(prefix, "_selex_type")]] <- 0
  il$data[[paste0("Use", use_field)]] <- array(1, dim = c(n_regions, n_years, 1, n_fleets))
  il$data[[paste0("Use", use_field, "_pop")]] <- array(0, dim = c(1, n_regions, n_years, 1, n_fleets))

  il$par[[paste0(prefix, "_fixed_sel_pars")]] <- array(0, dim = c(n_regions, n_bins, 1, n_sexes, n_fleets))
  il$par[[paste0(prefix, "sel_pe_pars")]] <- array(0, dim = c(n_regions, 4, n_sexes, n_fleets))
  il$par[[paste0("ln_", prefix, "sel_devs")]] <- array(0, dim = c(n_regions, n_years, n_bins, n_sexes, n_fleets))

  il
}

test_that("do_fixed_sel_pars_mapping resolves field names correctly for fish/ret/srv", {

  for (cfg in list(c(prefix = "fish", use_field = "Catch", fleet_field = "n_fish_fleets"),
                   c(prefix = "ret",  use_field = "Catch", fleet_field = "n_fish_fleets"),
                   c(prefix = "srv",  use_field = "SrvIdx", fleet_field = "n_srv_fleets"))) {

    il <- make_sel_input_list(cfg["prefix"], cfg["use_field"])
    out <- SPoRC:::do_fixed_sel_pars_mapping(
      il, sel_pars_spec = "est_all", bins = 2, sel_nonpar_est_bins = NULL,
      prefix = cfg["prefix"], fleet_field = cfg["fleet_field"], use_field = cfg["use_field"],
      fleet_label = "test fleet"
    )

    par_nm <- paste0(cfg["prefix"], "_fixed_sel_pars")
    map_vals <- as.integer(out$map[[par_nm]])
    expect_equal(sort(map_vals), 1:2) # logistic (model 0) has 2 parameters, both estimated (est_all)
  }
})

test_that("do_sel_pe_pars_mapping resolves field names correctly for fish/ret/srv", {

  for (cfg in list(c(prefix = "fish", use_field = "Catch", fleet_field = "n_fish_fleets"),
                   c(prefix = "ret",  use_field = "Catch", fleet_field = "n_fish_fleets"),
                   c(prefix = "srv",  use_field = "SrvIdx", fleet_field = "n_srv_fleets"))) {

    il <- make_sel_input_list(cfg["prefix"], cfg["use_field"])
    out <- SPoRC:::do_sel_pe_pars_mapping(
      il, pe_pars_spec = "est_all", corr_opt_semipar = NULL, bins = 2,
      prefix = cfg["prefix"], fleet_field = cfg["fleet_field"], use_field = cfg["use_field"],
      fleet_label = "test fleet"
    )

    par_nm <- paste0(cfg["prefix"], "sel_pe_pars")
    map_vals <- as.integer(out$map[[par_nm]])
    expect_equal(sort(map_vals[!is.na(map_vals)]), 1:2) # iid time-variation, logistic -> 2 pars estimated
  }
})

test_that("do_sel_pe_pars_mapping fixes all parameters when there is no time-variation", {

  il <- make_sel_input_list("fish", "Catch")
  il$data$cont_tv_fish_sel[1, 1] <- 0 # no time-variation

  out <- SPoRC:::do_sel_pe_pars_mapping(
    il, pe_pars_spec = "est_all", corr_opt_semipar = NULL, bins = 2,
    prefix = "fish", fleet_field = "n_fish_fleets", use_field = "Catch", fleet_label = "fishery fleet"
  )

  expect_true(all(is.na(out$map$fishsel_pe_pars)))
})

test_that("do_sel_devs_mapping resolves field names correctly for fish/ret/srv", {

  for (cfg in list(c(prefix = "fish", use_field = "Catch", fleet_field = "n_fish_fleets"),
                   c(prefix = "ret",  use_field = "Catch", fleet_field = "n_fish_fleets"),
                   c(prefix = "srv",  use_field = "SrvIdx", fleet_field = "n_srv_fleets"))) {

    il <- make_sel_input_list(cfg["prefix"], cfg["use_field"])
    out <- SPoRC:::do_sel_devs_mapping(
      il, sel_devs_spec = "est_all", sel_devs_shared_bins = NULL, bins = 2,
      prefix = cfg["prefix"], fleet_field = cfg["fleet_field"], use_field = cfg["use_field"],
      fleet_label = "test fleet"
    )

    par_nm <- paste0("ln_", cfg["prefix"], "sel_devs")
    map_vals <- as.integer(out$map[[par_nm]])
    expect_equal(sort(map_vals[!is.na(map_vals)]), 1:2) # 1 year, logistic -> 2 pars estimated

    # the equivalent integer array should also be attached to $data
    data_nm <- paste0("map_ln_", cfg["prefix"], "sel_devs")
    expect_true(!is.null(out$data[[data_nm]]))
  }
})

test_that("do_fixed_sel_pars_mapping fleet-sharing copies the reference fleet's map", {

  n_regions <- 1; n_fleets <- 2; n_sexes <- 1

  il <- list(data = list(), par = list(), map = list())
  il$data$n_regions <- n_regions
  il$data$n_sexes <- n_sexes
  il$data$n_fish_fleets <- n_fleets
  il$data$use_fixed_fish_sel <- c(0, 0)
  il$data$fish_sel_blocks <- array(1, dim = c(n_regions, 1, n_fleets))
  il$data$fish_sel_model <- array(0, dim = c(n_regions, 1, n_fleets)) # logistic -> 2 pars
  il$data$UseCatch <- array(1, dim = c(n_regions, 1, 1, n_fleets))
  il$data$UseCatch_pop <- array(0, dim = c(1, n_regions, 1, 1, n_fleets))
  il$par$fish_fixed_sel_pars <- array(0, dim = c(n_regions, 2, 1, n_sexes, n_fleets))

  out <- SPoRC:::do_fixed_sel_pars_mapping(
    il, sel_pars_spec = c("est_all", "est_shared_f_1"), bins = 2, sel_nonpar_est_bins = NULL,
    prefix = "fish", fleet_field = "n_fish_fleets", use_field = "Catch", fleet_label = "fishery fleet"
  )

  map_arr <- array(as.integer(out$map$fish_fixed_sel_pars), dim = c(n_regions, 2, 1, n_sexes, n_fleets))
  expect_equal(map_arr[,,,,1], map_arr[,,,,2]) # fleet 2 copies fleet 1's mapping
})
