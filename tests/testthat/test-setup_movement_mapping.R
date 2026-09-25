library(SPoRC)
library(testthat)

# ── shared helpers ───────────────────────────────────────────────────────────

# Minimal input_list stub with just enough of $data/$par for
# do_cont_vary_move_mapping() and Get_move_PE_loglik() to run, without going
# through the full Setup_Mod_* pipeline (no fishery/survey/tagging data
# required).
make_move_input_list <- function(
  n_pop = 2,
  n_regions = 3,
  n_yrs = 4,
  n_proj = 0,
  n_seas = 2,
  n_ages = 5,
  n_sexes = 2,
  do_recruits_move = 0,
  move_type = 0,
  cont_vary_movement,
  adjacency_collapsed = matrix(1, n_regions, n_regions - 1),
  adjacency_mat = matrix(1, n_regions, n_regions) - diag(n_regions)
) {

  n_yrs_devs <- n_yrs + n_proj
  n_dev_to <- ifelse(move_type == 1, 1, n_regions - 1) # CTMC deviations sit on a region's preference

  list(
    data = list(
      n_pop = n_pop,
      n_regions = n_regions,
      years = 1:n_yrs,
      n_proj_yrs_devs = n_proj,
      n_seas = n_seas,
      ages = 1:n_ages,
      n_sexes = n_sexes,
      do_recruits_move = do_recruits_move,
      use_fixed_movement = 0,
      move_type = move_type,
      cont_vary_movement = cont_vary_movement,
      adjacency_collapsed = adjacency_collapsed,
      adjacency_mat = adjacency_mat
    ),
    par = list(
      move_devs = array(0, dim = c(n_pop, n_regions, n_dev_to, n_yrs_devs, n_seas, n_ages, n_sexes)),
      move_pe_pars = array(0, dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes))
    ),
    map = list()
  )
}

# ── do_cont_vary_move_mapping: sharing structure ────────────────────────────

test_that("do_cont_vary_move_mapping builds the expected sharing structure", {

  n_pop <- 2
  n_regions <- 3
  n_yrs <- 4
  n_seas <- 2
  n_ages <- 5
  n_sexes <- 2
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
    il <- make_move_input_list(
      n_pop = n_pop,
      n_regions = n_regions,
      n_yrs = n_yrs,
      n_seas = n_seas,
      n_ages = n_ages,
      n_sexes = n_sexes,
      cont_vary_movement = spec_name
    )
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
  il <- make_move_input_list(cont_vary_movement = "none")
  il <- SPoRC:::do_cont_vary_move_mapping(il, "none", "fix")
  expect_true(all(is.na(il$data$map_move_devs)))
})

test_that("a CTMC deviation belongs to a region, and only an isolated region loses it", {
  # region 3 is cut off from both others, so shifting its preference moves nothing. regions 1 and 2
  # still trade, and each keeps its own deviation
  n_regions <- 3
  adjacency_mat <- matrix(0, n_regions, n_regions)
  adjacency_mat[1, 2] <- adjacency_mat[2, 1] <- 1

  il <- make_move_input_list(
    n_regions = n_regions,
    move_type = 1,
    cont_vary_movement = "iid_y",
    adjacency_mat = adjacency_mat
  )
  il <- SPoRC:::do_cont_vary_move_mapping(il, "iid_y", "fix")

  map_arr <- il$data$map_move_devs
  expect_equal(dim(map_arr)[3], 1) # one deviation per region, not per pair
  expect_true(all(is.na(map_arr[, 3, , , , , ])))
  expect_false(any(is.na(map_arr[, 1:2, , , , , ])))
})

test_that("a CTMC deviation reaches every edge of its own region", {
  # the payoff of holding the deviation on preference: raising region 1's preference pulls fish in
  # from both neighbours at once, where a diffusion deviation would have moved one edge
  A <- matrix(1, 3, 3) - diag(3)
  move_ctmc <- function(devs) {
    mv <- SPoRC:::Get_Movement(
      move_type = 1, do_recruits_move = 1, n_pop = 1, n_regions = 3, n_yrs = 1, n_proj_yrs_devs = 0,
      n_ages = 1, n_sexes = 1, n_seas = 1, move_pars = NULL,
      move_devs = array(devs, dim = c(1, 3, 1, 1, 1, 1, 1)), use_fixed_movement = 0,
      ctmc_move_dat = expand.grid(pop = 1, regions = 1:3, years = 1, seas = 1, ages = 1, sexes = 1),
      preference_formula = ~0, diffusion_formula = ~1,
      log_move_diffusion_pars = array(log(0.3), dim = c(1, 1)),
      move_preference_pars = array(0, dim = c(1, 1)), area_r = rep(1, 3), adjacency_mat = A,
      ctmc_diffusion_bounds = 0
    )
    list(M = matrix(mv$Movement[1, , , 1, 1, 1, 1], 3, 3), # [origin, destination]
         Q = matrix(mv$Mrate[1, , , 1, 1, 1, 1], 3, 3))
  }

  flat <- move_ctmc(rep(0, 3))
  lifted <- move_ctmc(c(0.02, 0, 0)) # well under the 0.09 diffusion rate, so every rate stays positive

  expect_equal(rowSums(lifted$M), rep(1, 3), tolerance = 1e-10)
  expect_true(lifted$Q[2, 1] > flat$Q[2, 1]) # region 2 sends faster to region 1
  expect_true(lifted$Q[3, 1] > flat$Q[3, 1]) # and so does region 3, off the one deviation
  expect_true(lifted$Q[1, 2] < flat$Q[1, 2]) # while region 1 lets go of fish more slowly
  expect_equal(lifted$Q[2, 3], flat$Q[2, 3], tolerance = 1e-12) # the edge that misses region 1 keeps its rate
  expect_true(lifted$M[1, 1] > flat$M[1, 1]) # region 1 holds a larger share of itself

  # a constant added to every region is a level shift of the surface, which the gradient drops
  expect_equal(move_ctmc(rep(0.4, 3))$M, flat$M, tolerance = 1e-12)
})

# ── Get_move_PE_loglik: likelihood values against a hand-computed baseline ──

test_that("Get_move_PE_loglik matches a hand-computed dnorm sum for every form", {

  n_pop <- 2
  n_regions <- 2
  n_yrs <- 3
  n_seas <- 2
  n_ages <- 4
  n_sexes <- 2
  adjacency_collapsed <- matrix(1, n_regions, n_regions - 1)

  dims <- c(n_pop, n_regions, n_regions - 1, n_yrs, n_seas, n_ages, n_sexes)
  set.seed(42)
  move_devs <- array(rnorm(prod(dims), sd = 0.3), dim = dims)
  PE_pars <- array(log(seq(0.2, 0.6, length.out = n_pop * n_regions * n_seas * n_ages * n_sexes)),
                   dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes))

  hand_ll <- function(form, do_recruits_move) {
    age_start <- if (do_recruits_move == 0 && n_ages >= 2) 2 else 1
    ll <- 0
    for (rr in 1:(n_regions - 1)) {
      for (r in 1:n_regions) {
        if (form == "iid_y") for (y in 1:n_yrs) ll <- ll + dnorm(move_devs[1, r, rr, y, 1, 1, 1], 0, exp(PE_pars[1, r, 1, 1, 1]), TRUE)
        if (form == "iid_a") for (a in age_start:n_ages) ll <- ll + dnorm(move_devs[1, r, rr, 1, 1, a, 1], 0, exp(PE_pars[1, r, 1, a, 1]), TRUE)
        if (form == "iid_y_a") for (y in 1:n_yrs) for (a in age_start:n_ages) ll <- ll + dnorm(move_devs[1, r, rr, y, 1, a, 1], 0, exp(PE_pars[1, r, 1, a, 1]), TRUE)
        if (form == "iid_y_a_s") for (y in 1:n_yrs) for (a in age_start:n_ages) for (s in 1:n_sexes) ll <- ll + dnorm(move_devs[1, r, rr, y, 1, a, s], 0, exp(PE_pars[1, r, 1, a, s]), TRUE)
        if (form == "iid_y_seas_a_s") for (y in 1:n_yrs) for (seas in 1:n_seas) for (a in age_start:n_ages) for (s in 1:n_sexes) ll <- ll + dnorm(move_devs[1, r, rr, y, seas, a, s], 0, exp(PE_pars[1, r, seas, a, s]), TRUE)
        if (form == "iid_p_y") for (p in 1:n_pop) for (y in 1:n_yrs) ll <- ll + dnorm(move_devs[p, r, rr, y, 1, 1, 1], 0, exp(PE_pars[p, r, 1, 1, 1]), TRUE)
        if (form == "iid_p_a") for (p in 1:n_pop) for (a in age_start:n_ages) ll <- ll + dnorm(move_devs[p, r, rr, 1, 1, a, 1], 0, exp(PE_pars[p, r, 1, a, 1]), TRUE)
        if (form == "iid_p_y_a") for (p in 1:n_pop) for (y in 1:n_yrs) for (a in age_start:n_ages) ll <- ll + dnorm(move_devs[p, r, rr, y, 1, a, 1], 0, exp(PE_pars[p, r, 1, a, 1]), TRUE)
        if (form == "iid_p_y_a_s") for (p in 1:n_pop) for (y in 1:n_yrs) for (a in age_start:n_ages) for (s in 1:n_sexes) ll <- ll + dnorm(move_devs[p, r, rr, y, 1, a, s], 0, exp(PE_pars[p, r, 1, a, s]), TRUE)
        if (form == "iid_p_y_seas_a_s") for (p in 1:n_pop) for (y in 1:n_yrs) for (seas in 1:n_seas) for (a in age_start:n_ages) for (s in 1:n_sexes)
          ll <- ll + dnorm(move_devs[p, r, rr, y, seas, a, s], 0, exp(PE_pars[p, r, seas, a, s]), TRUE)
      }
    }
    ll
  }

  forms <- c("iid_y", "iid_a", "iid_y_a", "iid_y_a_s", "iid_y_seas_a_s",
             "iid_p_y", "iid_p_a", "iid_p_y_a", "iid_p_y_a_s", "iid_p_y_seas_a_s")

  for (form in forms) {
    do_recruits_move <- 0
    got <- SPoRC:::Get_move_PE_loglik(
      cont_vary_movement = form,
      PE_pars = PE_pars,
      move_devs = move_devs,
      map_move_devs = array(0, dim = dims), # unused by the current implementation's math, only dims are read
      do_recruits_move = do_recruits_move,
      adjacency_collapsed = adjacency_collapsed,
      move_type = 0
    )
    expect_equal(got, hand_ll(form, do_recruits_move), tolerance = 1e-10,
                info = paste("form:", form))
  }
})
