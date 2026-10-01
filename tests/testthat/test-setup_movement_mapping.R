# Movement deviations: which cells share one, that a region's deviation reaches every edge
# touching it, and the Kronecker process error density against a normal density worked by hand.

library(SPoRC)
library(testthat)

# Helpers --------------------------------------------------------------------

# a minimal input list holding just enough for the mapping and the movement penalty to run,
# without going through the setup functions or needing any fishery or survey data
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
  adjacency_collapsed = matrix(1, n_regions, n_regions - 1),
  adjacency_mat = matrix(1, n_regions, n_regions) - diag(n_regions),
  move_re_pops = 1:n_pop,
  move_re_years = 1:n_yrs,
  move_re_seas = 1:n_seas,
  move_re_ages = if(do_recruits_move == 0) 2:n_ages else 1:n_ages,
  move_re_sexes = 1:n_sexes
) {

  n_dev_to <- ifelse(move_type == 1, 1, n_regions - 1) # CTMC deviations sit on a region's preference
  n_corr <- function(n) max(1, n * (n - 1) / 2)

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
      adjacency_collapsed = adjacency_collapsed,
      adjacency_mat = adjacency_mat,
      move_re_pops = move_re_pops,
      move_re_years = move_re_years,
      move_re_seas = move_re_seas,
      move_re_ages = move_re_ages,
      move_re_sexes = move_re_sexes
    ),
    par = list(
      move_devs = array(0, dim = c(n_pop, n_regions, n_dev_to, n_yrs + n_proj, n_seas, n_ages, n_sexes)),
      move_pe_pars = array(0, dim = c(n_regions, n_dev_to, 3)),
      move_pop_corr_pars = rep(0, n_corr(n_pop)),
      move_seas_corr_pars = rep(0, n_corr(n_seas)),
      move_sex_corr_pars = rep(0, n_corr(n_sexes))
    ),
    map = list()
  )
}

n_levels <- function(m) length(unique(stats::na.omit(as.vector(m))))

# every switch off unless named
map_re <- function(il, year = "none", age = "none", pop = "none", seas = "none", sex = "none", pe_spec = "est_all") {
  SPoRC:::do_move_re_mapping(il, year, age, pop, seas, sex, pe_spec)
}

# log density of a mean zero multivariate normal, by its Cholesky factor
mvn_logpdf <- function(x, S) {
  L <- chol(S)
  z <- backsolve(L, x, transpose = TRUE)
  -0.5 * (length(x) * log(2 * pi) + 2 * sum(log(diag(L))) + sum(z^2))
}
ar1_corr <- function(n, rho) rho^abs(outer(seq_len(n), seq_len(n), "-"))
us_corr <- function(pars, n) { L <- SPoRC:::build_us_chol(pars, n); L %*% t(L) }

# Which Movement Deviations Are Shared ---------------------------------------

test_that("each switch says what varies, and a block list shares within its blocks", {

  il <- make_move_input_list()
  n_pairs <- 3 * 2 # origin and destination in the collapsed array
  n_yrs <- 4
  n_movable <- 5 - 1 # recruits do not move

  # year deviations, one series per pair, shared across everything else
  m <- map_re(il, year = "iid")
  expect_equal(nrow(m$data$move_pairs), n_pairs)
  expect_equal(n_levels(m$data$map_move_devs), n_pairs * n_yrs)
  expect_true(all(is.na(m$data$map_move_devs[,,,,,1,]))) # age one holds no deviation when recruits do not move
  expect_equal(m$data$map_move_devs[1,1,1,,1,2,1], m$data$map_move_devs[2,1,1,,2,4,2]) # the same deviation in every shared slot
  expect_equal(n_levels(m$map$move_pe_pars), n_pairs) # one log sd per pair, no correlation
  expect_true(all(is.na(array(as.integer(m$map$move_pe_pars), dim = c(3, 2, 3))[,,2:3])))
  expect_true(all(is.na(m$map$move_sex_corr_pars)))

  # each switch multiplies the levels by its active count
  expect_equal(n_levels(map_re(il, age = "iid")$data$map_move_devs), n_pairs * n_movable)
  expect_equal(n_levels(map_re(il, year = "iid", age = "iid")$data$map_move_devs), n_pairs * n_yrs * n_movable)
  expect_equal(n_levels(map_re(il, year = "iid", pop = "iid")$data$map_move_devs), n_pairs * n_yrs * 2)
  expect_equal(n_levels(map_re(il, year = "iid", seas = "iid", sex = "iid")$data$map_move_devs), n_pairs * n_yrs * 2 * 2)
  expect_equal(n_levels(map_re(il, year = "iid", age = "iid", pop = "us", seas = "us", sex = "us")$data$map_move_devs), n_pairs * n_yrs * n_movable * 8)

  # the process error can be shared across pairs in blocks, the deviations stay per pair
  m <- map_re(il, year = "iid", pe_spec = list(1:3, 4:6))
  expect_equal(n_levels(m$data$map_move_devs), n_pairs * n_yrs)
  expect_equal(n_levels(m$map$move_pe_pars), 2)

  # blocks within a dim: the sexes in a block share, and a block per level is iid
  expect_equal(map_re(il, year = "iid", sex = list(1:2))$data$map_move_devs, map_re(il, year = "iid", sex = "none")$data$map_move_devs)
  expect_equal(map_re(il, year = "iid", pop = list(1, 2))$data$map_move_devs, map_re(il, year = "iid", pop = "iid")$data$map_move_devs)
  three_seas <- make_move_input_list(n_seas = 3)
  m <- map_re(three_seas, year = "iid", seas = list(1:2, 3))
  expect_equal(n_levels(m$data$map_move_devs), n_pairs * n_yrs * 2)
  expect_equal(m$data$map_move_devs[1,1,1,,1,2,1], m$data$map_move_devs[1,1,1,,2,2,1])
  expect_false(any(m$data$map_move_devs[1,1,1,,1,2,1] == m$data$map_move_devs[1,1,1,,3,2,1]))

  # a correlation parameter exists only where its dim asks for one, and the unstructured ones are shared by every pair
  pe <- array(as.integer(map_re(il, year = "ar1")$map$move_pe_pars), dim = c(3, 2, 3))
  expect_true(all(is.na(pe[,,2])) && all(!is.na(pe[,,3])))
  pe <- array(as.integer(map_re(il, age = "ar1")$map$move_pe_pars), dim = c(3, 2, 3))
  expect_true(all(!is.na(pe[,,2])) && all(is.na(pe[,,3])))
  m <- map_re(il, year = "ar1", age = "ar1", pop = "us", seas = "us", sex = "us")
  expect_equal(n_levels(m$map$move_pe_pars), 3 * n_pairs)
  expect_equal(n_levels(m$map$move_pop_corr_pars), 1)
  expect_equal(n_levels(m$map$move_seas_corr_pars), 1)
  expect_equal(n_levels(m$map$move_sex_corr_pars), 1)

  # one sd and one pair of correlations for every pair
  m <- map_re(il, year = "ar1", age = "ar1", pe_spec = "est_shared")
  expect_equal(n_levels(m$map$move_pe_pars), 3)
  expect_equal(n_levels(m$data$map_move_devs), n_pairs * n_yrs * n_movable) # the deviations themselves stay per pair

  # every switch off is no deviations at all
  m <- map_re(il)
  expect_true(all(is.na(m$data$map_move_devs)))
  expect_equal(nrow(m$data$move_pairs), 0)
  expect_true(all(is.na(m$map$move_pe_pars)))

  # a dsem splits every dim with more than one level, since it links one series each, and reads no process error
  m <- map_re(il, year = "dsem")
  expect_equal(n_levels(m$data$map_move_devs), n_pairs * n_yrs * n_movable * 2 * 2 * 2)
  expect_true(all(is.na(m$map$move_pe_pars)))

  # the active sets: only those levels hold deviations, and a shared dim reads its first active slot
  sub <- make_move_input_list(move_re_years = 2:3, move_re_ages = 3:5)
  m <- map_re(sub, year = "iid", age = "iid")
  expect_equal(n_levels(m$data$map_move_devs), n_pairs * 2 * 3)
  expect_true(all(is.na(m$data$map_move_devs[,,,c(1, 4),,,])))
  expect_true(all(is.na(m$data$map_move_devs[,,,,,1:2,])))
  m <- map_re(sub, year = "iid")
  expect_true(all(is.na(m$data$map_move_devs[,,,,,1:2,])) && all(!is.na(m$data$map_move_devs[,,,2:3,,3:5,])))
  one_sex <- make_move_input_list(move_re_sexes = 2, move_re_pops = 1)
  m <- map_re(one_sex, year = "iid", sex = "iid", pop = "iid")
  expect_true(all(is.na(m$data$map_move_devs[,,,,,,1])) && all(is.na(m$data$map_move_devs[2,,,,,,])))
  expect_equal(n_levels(m$data$map_move_devs), n_pairs * n_yrs)

  # a pair with no edge has no deviations
  cut <- make_move_input_list(adjacency_collapsed = matrix(c(1, 1, 1, 0, 1, 1), 3, 2))
  m <- map_re(cut, year = "iid")
  expect_equal(nrow(m$data$move_pairs), n_pairs - 1)
  expect_true(all(is.na(m$data$map_move_devs[,1,2,,,,])))

})

test_that("the setup refuses structure it cannot hold", {

  move_on <- list(use_fixed_movement = 0, Fixed_Movement = NA)
  expect_error(suppressMessages(sweep_input(move = c(move_on, list(move_year_re = "iid", move_sex_re = list(1))))), "cover each active level")
  expect_error(suppressMessages(sweep_input(move = c(move_on, list(move_year_re = "dsem", move_sex_re = list(1, 2))))), "block list")
  expect_error(suppressMessages(sweep_input(move = c(move_on, list(move_year_re = list(1:6, 7:13))))), "Movement_yearblk_spec")
  expect_error(suppressMessages(sweep_input(move = c(move_on, list(move_year_re = "iid", move_sex_re = "us", move_re_sexes = 1)))), "two active sexes")
  ok <- suppressMessages(sweep_input(move = c(move_on, list(move_year_re = "ar1", move_sex_re = list(1:2), move_pe_spec = list(1:3, 4:6)))))
  expect_equal(length(levels(ok$map$move_pe_pars)), 2 * 2) # a log sd and a year correlation per pair block

})

test_that("a CTMC deviation belongs to a region, and only an isolated region loses it", {
  # region 3 is cut off from both others, so shifting its preference moves nothing. regions 1 and 2
  # still trade, and each keeps its own deviation
  n_regions <- 3
  adjacency_mat <- matrix(0, n_regions, n_regions)
  adjacency_mat[1, 2] <- adjacency_mat[2, 1] <- 1

  il <- make_move_input_list(n_regions = n_regions, move_type = 1, adjacency_mat = adjacency_mat)
  m <- map_re(il, year = "iid")

  map_arr <- m$data$map_move_devs
  expect_equal(unname(dim(map_arr)[3]), 1) # one deviation per region, not per pair
  expect_true(all(is.na(map_arr[, 3, , , , , ])))
  expect_false(any(is.na(map_arr[, 1:2, , , , 2:5, ])))
  expect_equal(nrow(m$data$move_pairs), 2)
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

# The Movement Penalty against a Hand Calculation ----------------------------

test_that("Get_move_PE_loglik matches the density by hand for every form", {

  il <- make_move_input_list(n_pop = 2, n_regions = 2, n_yrs = 3, n_seas = 2, n_ages = 4, n_sexes = 2)
  dims <- dim(il$par$move_devs)
  set.seed(42)
  move_devs <- array(rnorm(prod(dims), sd = 0.3), dim = dims)
  sigma <- 0.25 # the conditional sd the parameter holds
  rho_a <- SPoRC:::rho_trans(0.4)
  rho_y <- SPoRC:::rho_trans(-0.6)
  s_y <- sigma / sqrt(1 - rho_y^2) # marginal sd under an ar1 over years
  s_a <- sigma / sqrt(1 - rho_a^2) # and over ages
  s_ya <- sigma / sqrt(1 - rho_y^2) / sqrt(1 - rho_a^2) # and over both
  pe <- il$par$move_pe_pars
  pe[,,1] <- log(sigma); pe[,,2] <- 0.4; pe[,,3] <- -0.6
  pop_corr <- 0.3; seas_corr <- -0.5; sex_corr <- 0.8 # one unconstrained parameter each for two levels
  ages <- 2:4 # recruits do not move
  code <- c(none = 0, iid = 1, ar1 = 2, us = 2)

  # the penalty as the objective calls it, with the blocks the mapping wrote
  pen <- function(year = "none", age = "none", pop = "none", seas = "none", sex = "none", pe_spec = "est_all", dsem = 0) {
    d <- map_re(il, year, age, pop, seas, sex, pe_spec)$data
    SPoRC:::Get_move_PE_loglik(move_year_re = code[[year]], move_age_re = code[[age]], move_pop_re = code[[pop]], move_seas_re = code[[seas]], move_sex_re = code[[sex]],
                               PE_pars = pe, move_pop_corr_pars = pop_corr, move_seas_corr_pars = seas_corr, move_sex_corr_pars = sex_corr,
                               move_devs = move_devs, map_move_devs = d$map_move_devs, move_pairs = d$move_pairs,
                               move_pe_block = d$move_pe_block, move_pop_block = d$move_pop_block, move_year_block = d$move_year_block,
                               move_seas_block = d$move_seas_block, move_age_block = d$move_age_block, move_sex_block = d$move_sex_block, move_dsem = dsem)
  }
  # a shared dim is read at its first active slot: population 1, season 1, age 2, sex 1
  over_pairs <- function(f) { ll <- 0; for(r in 1:2) for(rr in 1) ll <- ll + f(r, rr); ll }

  # independent everywhere: a plain normal sum over the penalized cells
  expect_equal(pen(year = "iid"), over_pairs(function(r, rr) sum(dnorm(move_devs[1,r,rr,,1,2,1], 0, sigma, TRUE))), tolerance = 1e-10)
  expect_equal(pen(age = "iid"), over_pairs(function(r, rr) sum(dnorm(move_devs[1,r,rr,1,1,ages,1], 0, sigma, TRUE))), tolerance = 1e-10)
  expect_equal(pen(year = "iid", age = "iid", sex = "iid"), over_pairs(function(r, rr) sum(dnorm(move_devs[1,r,rr,,1,ages,], 0, sigma, TRUE))), tolerance = 1e-10)

  # a correlation on one dim: an ar1 multivariate normal with marginal sd sigma
  expect_equal(pen(year = "ar1"), over_pairs(function(r, rr) mvn_logpdf(move_devs[1,r,rr,,1,2,1], s_y^2 * ar1_corr(3, rho_y))), tolerance = 1e-10)
  expect_equal(pen(age = "ar1"), over_pairs(function(r, rr) mvn_logpdf(move_devs[1,r,rr,1,1,ages,1], s_a^2 * ar1_corr(3, rho_a))), tolerance = 1e-10)
  expect_equal(pen(sex = "us"), over_pairs(function(r, rr) mvn_logpdf(move_devs[1,r,rr,1,1,2,], sigma^2 * us_corr(sex_corr, 2))), tolerance = 1e-10)

  # years and ages: the separable covariance, years running fastest in the vector
  S_ya <- s_ya^2 * kronecker(ar1_corr(3, rho_a), ar1_corr(3, rho_y))
  expect_equal(pen(year = "ar1", age = "ar1"), over_pairs(function(r, rr) mvn_logpdf(as.vector(move_devs[1,r,rr,,1,ages,1]), S_ya)), tolerance = 1e-10)
  S_y <- s_y^2 * kronecker(diag(3), ar1_corr(3, rho_y))
  expect_equal(pen(year = "ar1", age = "iid"), over_pairs(function(r, rr) mvn_logpdf(as.vector(move_devs[1,r,rr,,1,ages,1]), S_y)), tolerance = 1e-10)

  # every dim at once, with one shared sd: sex, age, season, year and population from slowest to fastest, pairs independent
  S_all <- s_ya^2 * kronecker(us_corr(sex_corr, 2), kronecker(ar1_corr(3, rho_a), kronecker(us_corr(seas_corr, 2), kronecker(ar1_corr(3, rho_y), us_corr(pop_corr, 2)))))
  expect_equal(pen("ar1", "ar1", "us", "us", "us", pe_spec = "est_shared"), over_pairs(function(r, rr) mvn_logpdf(as.vector(move_devs[,r,rr,,,ages,]), S_all)), tolerance = 1e-10)

  # nothing to penalize
  expect_equal(pen(), 0)
  expect_equal(pen("ar1", "ar1", dsem = 1), 0)

})
