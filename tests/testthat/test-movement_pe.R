library(SPoRC)
library(testthat)

# ── shared helpers ───────────────────────────────────────────────────────────

# Minimal input_list stub with just enough of $data/$par for
# do_cont_vary_move_mapping() and Get_move_PE_loglik() to run, without going
# through the full Setup_Mod_* pipeline (no fishery/survey/tagging data
# required).
make_move_input_list <- function(n_pop = 2, n_regions = 3, n_yrs = 4, n_proj = 0,
                                  n_seas = 2, n_ages = 5, n_sexes = 2,
                                  do_recruits_move = 0, move_type = 0,
                                  cont_vary_movement_val,
                                  adjacency_collapsed = matrix(1, n_regions, n_regions - 1)) {

  n_yrs_devs <- n_yrs + n_proj

  list(
    data = list(
      n_pop = n_pop, n_regions = n_regions,
      years = 1:n_yrs, n_proj_yrs_devs = n_proj,
      n_seas = n_seas, ages = 1:n_ages, n_sexes = n_sexes,
      do_recruits_move = do_recruits_move,
      use_fixed_movement = 0,
      move_type = move_type,
      cont_vary_movement = cont_vary_movement_val,
      adjacency_collapsed = adjacency_collapsed
    ),
    par = list(
      move_devs = array(0, dim = c(n_pop, n_regions, n_regions - 1, n_yrs_devs, n_seas, n_ages, n_sexes)),
      move_pe_pars = array(0, dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes))
    ),
    map = list()
  )
}

cont_move_map <- data.frame(
  type = c("none", "iid_y", "iid_a", "iid_y_a", "iid_y_a_s", "iid_y_seas_a_s",
           "iid_p_y", "iid_p_a", "iid_p_y_a", "iid_p_y_a_s", "iid_p_y_seas_a_s"),
  num = 0:10
)

# ── do_cont_vary_move_mapping: sharing structure ────────────────────────────

test_that("do_cont_vary_move_mapping builds the expected sharing structure", {

  n_pop <- 2; n_regions <- 3; n_yrs <- 4; n_seas <- 2; n_ages <- 5; n_sexes <- 2
  n_movable_ages <- n_ages - 1 # do_recruits_move = 0 -> age 1 excluded
  n_pairs <- n_regions * (n_regions - 1) # (from, to) combinations in the collapsed array

  # expected number of distinct estimated groups PER (from, to) region pair
  expected_groups_per_pair <- c(
    iid_y            = n_yrs,
    iid_a            = n_movable_ages,
    iid_y_a          = n_yrs * n_movable_ages,
    iid_y_a_s        = n_yrs * n_movable_ages * n_sexes,
    iid_y_seas_a_s   = n_yrs * n_seas * n_movable_ages * n_sexes,
    iid_p_y          = n_pop * n_yrs,
    iid_p_a          = n_pop * n_movable_ages,
    iid_p_y_a        = n_pop * n_yrs * n_movable_ages,
    iid_p_y_a_s      = n_pop * n_yrs * n_movable_ages * n_sexes,
    iid_p_y_seas_a_s = n_pop * n_yrs * n_seas * n_movable_ages * n_sexes
  )

  for (spec_name in names(expected_groups_per_pair)) {
    val <- cont_move_map$num[cont_move_map$type == spec_name]
    il <- make_move_input_list(n_pop = n_pop, n_regions = n_regions, n_yrs = n_yrs,
                               n_seas = n_seas, n_ages = n_ages, n_sexes = n_sexes,
                               cont_vary_movement_val = val)
    il <- SPoRC:::do_cont_vary_move_mapping(il, spec_name, "fix")

    map_arr <- il$data$map_move_devs
    non_na_ids <- unique(as.vector(map_arr))
    non_na_ids <- non_na_ids[!is.na(non_na_ids)]

    expect_equal(length(non_na_ids), n_pairs * unname(expected_groups_per_pair[spec_name]),
                info = paste("spec:", spec_name))

    # Age is only excluded for recruits (age 1) when age is itself a key
    # dimension of the spec (i.e. "a" appears in the name); when age is a
    # broadcast dim (e.g. "iid_y"), age 1 shares the same tied parameter as
    # every other age rather than being masked to NA.
    has_age_dim <- grepl("(^|_)a(_|$)", spec_name)
    if (has_age_dim) {
      expect_true(all(is.na(map_arr[, , , , , 1, ])), info = paste("recruit age NA, spec:", spec_name))
    } else {
      expect_true(all(!is.na(map_arr[, , , , , 1, ])), info = paste("recruit age shared, spec:", spec_name))
    }
  }
})

test_that("do_cont_vary_move_mapping returns all-NA map when cont_vary_movement is 'none'", {
  il <- make_move_input_list(cont_vary_movement_val = 0)
  il <- SPoRC:::do_cont_vary_move_mapping(il, "none", "fix")
  expect_true(all(is.na(il$data$map_move_devs)))
})

test_that("do_cont_vary_move_mapping respects CTMC adjacency masking", {
  n_regions <- 3
  adjacency_collapsed <- matrix(1, n_regions, n_regions - 1)
  adjacency_collapsed[1, 1] <- 0 # region 1 -> its first "other" region is not adjacent

  il <- make_move_input_list(n_regions = n_regions, move_type = 1,
                             cont_vary_movement_val = 1,
                             adjacency_collapsed = adjacency_collapsed)
  il <- SPoRC:::do_cont_vary_move_mapping(il, "iid_y", "fix")

  map_arr <- il$data$map_move_devs
  expect_true(all(is.na(map_arr[, 1, 1, , , , ])))
  expect_false(all(is.na(map_arr[, 1, 2, , , , ])))
})

# ── Get_move_PE_loglik: likelihood values against a hand-computed baseline ──

test_that("Get_move_PE_loglik matches a hand-computed dnorm sum for each PE_model", {

  n_pop <- 2; n_regions <- 2; n_yrs <- 3; n_seas <- 2; n_ages <- 4; n_sexes <- 2
  adjacency_collapsed <- matrix(1, n_regions, n_regions - 1)

  dims <- c(n_pop, n_regions, n_regions - 1, n_yrs, n_seas, n_ages, n_sexes)
  set.seed(42)
  move_devs <- array(rnorm(prod(dims), sd = 0.3), dim = dims)
  PE_pars <- array(log(seq(0.2, 0.6, length.out = n_pop * n_regions * n_seas * n_ages * n_sexes)),
                   dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes))

  hand_ll <- function(PE_model, do_recruits_move) {
    age_start <- if (do_recruits_move == 0 && n_ages >= 2) 2 else 1
    ll <- 0
    for (rr in 1:(n_regions - 1)) {
      for (r in 1:n_regions) {
        if (PE_model == 1) for (y in 1:n_yrs) ll <- ll + dnorm(move_devs[1, r, rr, y, 1, 1, 1], 0, exp(PE_pars[1, r, 1, 1, 1]), TRUE)
        if (PE_model == 2) for (a in age_start:n_ages) ll <- ll + dnorm(move_devs[1, r, rr, 1, 1, a, 1], 0, exp(PE_pars[1, r, 1, a, 1]), TRUE)
        if (PE_model == 3) for (y in 1:n_yrs) for (a in age_start:n_ages) ll <- ll + dnorm(move_devs[1, r, rr, y, 1, a, 1], 0, exp(PE_pars[1, r, 1, a, 1]), TRUE)
        if (PE_model == 4) for (y in 1:n_yrs) for (a in age_start:n_ages) for (s in 1:n_sexes) ll <- ll + dnorm(move_devs[1, r, rr, y, 1, a, s], 0, exp(PE_pars[1, r, 1, a, s]), TRUE)
        if (PE_model == 5) for (y in 1:n_yrs) for (seas in 1:n_seas) for (a in age_start:n_ages) for (s in 1:n_sexes) ll <- ll + dnorm(move_devs[1, r, rr, y, seas, a, s], 0, exp(PE_pars[1, r, seas, a, s]), TRUE)
        if (PE_model == 6) for (p in 1:n_pop) for (y in 1:n_yrs) ll <- ll + dnorm(move_devs[p, r, rr, y, 1, 1, 1], 0, exp(PE_pars[p, r, 1, 1, 1]), TRUE)
        if (PE_model == 7) for (p in 1:n_pop) for (a in age_start:n_ages) ll <- ll + dnorm(move_devs[p, r, rr, 1, 1, a, 1], 0, exp(PE_pars[p, r, 1, a, 1]), TRUE)
        if (PE_model == 8) for (p in 1:n_pop) for (y in 1:n_yrs) for (a in age_start:n_ages) ll <- ll + dnorm(move_devs[p, r, rr, y, 1, a, 1], 0, exp(PE_pars[p, r, 1, a, 1]), TRUE)
        if (PE_model == 9) for (p in 1:n_pop) for (y in 1:n_yrs) for (a in age_start:n_ages) for (s in 1:n_sexes) ll <- ll + dnorm(move_devs[p, r, rr, y, 1, a, s], 0, exp(PE_pars[p, r, 1, a, s]), TRUE)
        if (PE_model == 10) for (p in 1:n_pop) for (y in 1:n_yrs) for (seas in 1:n_seas) for (a in age_start:n_ages) for (s in 1:n_sexes)
          ll <- ll + dnorm(move_devs[p, r, rr, y, seas, a, s], 0, exp(PE_pars[p, r, seas, a, s]), TRUE)
      }
    }
    ll
  }

  for (PE_model in 1:10) {
    do_recruits_move <- 0
    got <- SPoRC:::Get_move_PE_loglik(
      PE_model = PE_model, PE_pars = PE_pars, move_devs = move_devs,
      map_move_devs = array(0, dim = dims), # unused by the current implementation's math, only dims are read
      do_recruits_move = do_recruits_move,
      adjacency_collapsed = adjacency_collapsed, move_type = 0
    )
    expect_equal(got, hand_ll(PE_model, do_recruits_move), tolerance = 1e-10,
                info = paste("PE_model:", PE_model))
  }
})
