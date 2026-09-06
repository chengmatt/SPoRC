library(SPoRC)
library(testthat)

# Helper to build minimal NAA, fish_sel, srv_sel arrays for tests.
# Dimensions: NAA[n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes]
#             fish_sel[n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_fish_fleets]
#             srv_sel [n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_srv_fleets]
make_arrays <- function(
  n_pop = 2,
  n_ages = 3,
  n_sexes = 2,
  n_regions = 2,
  n_yrs = 5,
  n_seas = 1,
  n_fish_fleets = 1,
  n_srv_fleets = 1,
  naa_values = NULL
) {
  naa_dim  <- c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes)
  fsel_dim <- c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_fish_fleets)
  ssel_dim <- c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_srv_fleets)

  NAA      <- if (is.null(naa_values)) array(1, dim = naa_dim) else
    array(naa_values, dim = naa_dim)
  fish_sel <- array(1, dim = fsel_dim)
  srv_sel  <- array(1, dim = ssel_dim)
  list(
    NAA = NAA,
    fish_sel = fish_sel,
    srv_sel = srv_sel,
    n_pop = n_pop,
    n_ages = n_ages,
    n_sexes = n_sexes
  )
}

test_that("release_conv_tag_attr: p_a_s returns tagged_fish unchanged", {
  arrs <- make_arrays()
  tagged_fish <- array(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
                       dim = c(arrs$n_pop, arrs$n_ages, arrs$n_sexes))

  result <- release_conv_tag_attr(
    tagged_fish          = tagged_fish,
    tag_attr             = "p_a_s",
    tag_release_platform = c("population", NA),
    srv_sel              = arrs$srv_sel,
    fish_sel             = arrs$fish_sel,
    NAA                  = arrs$NAA,
    ty = 1,
    tseas = 1,
    tr = 1,
    n_pop = arrs$n_pop,
    n_ages = arrs$n_ages,
    n_sexes = arrs$n_sexes
  )

  expect_equal(result, array(tagged_fish, dim = c(arrs$n_pop, arrs$n_ages, arrs$n_sexes)))
})

test_that("release_conv_tag_attr: output sums match input totals for all tag_attr combos", {
  # With uniform weights, total fish should be preserved regardless of attr
  arrs <- make_arrays(n_pop = 2, n_ages = 3, n_sexes = 2)
  total_tags <- 120

  attrs_and_dims <- list(
    list(attr = "p",     dim = c(2, 1, 1)),
    list(attr = "a",     dim = c(1, 3, 1)),
    list(attr = "s",     dim = c(1, 1, 2)),
    list(attr = "p_a",   dim = c(2, 3, 1)),
    list(attr = "p_s",   dim = c(2, 1, 2)),
    list(attr = "a_s",   dim = c(1, 3, 2))
  )

  for (cfg in attrs_and_dims) {
    tagged_fish <- array(total_tags / prod(cfg$dim), dim = cfg$dim)
    result <- release_conv_tag_attr(
      tagged_fish          = tagged_fish,
      tag_attr             = cfg$attr,
      tag_release_platform = c("population", NA),
      srv_sel              = arrs$srv_sel,
      fish_sel             = arrs$fish_sel,
      NAA                  = arrs$NAA,
      ty = 1,
      tseas = 1,
      tr = 1,
      n_pop = arrs$n_pop,
      n_ages = arrs$n_ages,
      n_sexes = arrs$n_sexes
    )
    expect_equal(sum(result), total_tags, tolerance = 1e-10,
                 label = paste("total preserved for tag_attr =", cfg$attr))
  }
})

test_that("release_conv_tag_attr: marginals along attended dims are preserved", {
  # With uniform NAA weights, marginal totals along each attended dimension
  # in the output must exactly match the input values for that dimension.
  arrs <- make_arrays(n_pop = 2, n_ages = 3, n_sexes = 2)

  # tag_attr = "p": only population is attended; each pop has a distinct total
  tagged_p <- array(0, dim = c(2, 1, 1))
  tagged_p[1, 1, 1] <- 30
  tagged_p[2, 1, 1] <- 70
  result_p <- release_conv_tag_attr(
    tagged_fish = tagged_p,
    tag_attr = "p",
    tag_release_platform = c("population", NA),
    srv_sel = arrs$srv_sel,
    fish_sel = arrs$fish_sel,
    NAA = arrs$NAA,
    ty = 1,
    tseas = 1,
    tr = 1,
    n_pop = arrs$n_pop,
    n_ages = arrs$n_ages,
    n_sexes = arrs$n_sexes
  )
  # Sum across age and sex for each pop should equal input
  expect_equal(sum(result_p[1, , ]), 30, tolerance = 1e-10)
  expect_equal(sum(result_p[2, , ]), 70, tolerance = 1e-10)

  # tag_attr = "a": only age is attended; each age bin has a distinct total
  tagged_a <- array(0, dim = c(1, 3, 1))
  tagged_a[1, 1, 1] <- 10
  tagged_a[1, 2, 1] <- 40
  tagged_a[1, 3, 1] <- 50
  result_a <- release_conv_tag_attr(
    tagged_fish = tagged_a,
    tag_attr = "a",
    tag_release_platform = c("population", NA),
    srv_sel = arrs$srv_sel,
    fish_sel = arrs$fish_sel,
    NAA = arrs$NAA,
    ty = 1,
    tseas = 1,
    tr = 1,
    n_pop = arrs$n_pop,
    n_ages = arrs$n_ages,
    n_sexes = arrs$n_sexes
  )
  expect_equal(sum(result_a[, 1, ]), 10, tolerance = 1e-10)
  expect_equal(sum(result_a[, 2, ]), 40, tolerance = 1e-10)
  expect_equal(sum(result_a[, 3, ]), 50, tolerance = 1e-10)

  # tag_attr = "s": only sex is attended
  tagged_s <- array(0, dim = c(1, 1, 2))
  tagged_s[1, 1, 1] <- 60
  tagged_s[1, 1, 2] <- 40
  result_s <- release_conv_tag_attr(
    tagged_fish = tagged_s,
    tag_attr = "s",
    tag_release_platform = c("population", NA),
    srv_sel = arrs$srv_sel,
    fish_sel = arrs$fish_sel,
    NAA = arrs$NAA,
    ty = 1,
    tseas = 1,
    tr = 1,
    n_pop = arrs$n_pop,
    n_ages = arrs$n_ages,
    n_sexes = arrs$n_sexes
  )
  expect_equal(sum(result_s[, , 1]), 60, tolerance = 1e-10)
  expect_equal(sum(result_s[, , 2]), 40, tolerance = 1e-10)

  # tag_attr = "p_a": both pop and age attended; check each [p, a] marginal
  tagged_pa <- array(1:6, dim = c(2, 3, 1))
  result_pa <- release_conv_tag_attr(
    tagged_fish = tagged_pa,
    tag_attr = "p_a",
    tag_release_platform = c("population", NA),
    srv_sel = arrs$srv_sel,
    fish_sel = arrs$fish_sel,
    NAA = arrs$NAA,
    ty = 1,
    tseas = 1,
    tr = 1,
    n_pop = arrs$n_pop,
    n_ages = arrs$n_ages,
    n_sexes = arrs$n_sexes
  )
  for (p in 1:2) for (a in 1:3) {
    expect_equal(sum(result_pa[p, a, ]), tagged_pa[p, a, 1], tolerance = 1e-10,
                 label = paste("p_a marginal preserved at p =", p, "a =", a))
  }
})

test_that("release_conv_tag_attr: uniform weights distribute evenly", {
  # Uniform NAA => each [p, a, s] cell should receive equal share
  arrs <- make_arrays(n_pop = 2, n_ages = 3, n_sexes = 2)
  total_tags <- 120
  # No attended dims means all weight spread equally across all 12 cells
  tagged_none <- array(total_tags, dim = c(1, 1, 1))
  result <- release_conv_tag_attr(
    tagged_fish          = tagged_none,
    tag_attr             = "", # no attended dims
    tag_release_platform = c("population", NA),
    srv_sel = arrs$srv_sel,
    fish_sel = arrs$fish_sel,
    NAA = arrs$NAA,
    ty = 1,
    tseas = 1,
    tr = 1,
    n_pop = arrs$n_pop,
    n_ages = arrs$n_ages,
    n_sexes = arrs$n_sexes
  )
  n_cells <- arrs$n_pop * arrs$n_ages * arrs$n_sexes   # 12
  expected_per_cell <- total_tags / n_cells
  expect_true(all(abs(as.vector(result) - expected_per_cell) < 1e-10))
})

test_that("release_conv_tag_attr: output has correct dimensions", {
  arrs <- make_arrays(n_pop = 3, n_ages = 4, n_sexes = 2)
  tagged_fish <- array(10, dim = c(1, 1, 1))
  result <- release_conv_tag_attr(
    tagged_fish          = tagged_fish,
    tag_attr             = "",
    tag_release_platform = c("population", NA),
    srv_sel = arrs$srv_sel,
    fish_sel = arrs$fish_sel,
    NAA = arrs$NAA,
    ty = 1,
    tseas = 1,
    tr = 1,
    n_pop = arrs$n_pop,
    n_ages = arrs$n_ages,
    n_sexes = arrs$n_sexes
  )
  expect_equal(dim(result), c(3, 4, 2))
})

test_that("release_conv_tag_attr: platform = 'fishery' uses fish_sel weights", {
  # Make fish_sel zero for pop 1 and nonzero for pop 2; tags should go entirely
  # to pop 2 regardless of NAA being uniform.
  arrs <- make_arrays(n_pop = 2, n_ages = 3, n_sexes = 2)
  arrs$fish_sel[1, , , , , , 1] <- 0   # zero out pop 1 fishery selectivity
  arrs$fish_sel[2, , , , , , 1] <- 1

  tagged_fish <- array(100, dim = c(1, 1, 1))
  result <- release_conv_tag_attr(
    tagged_fish          = tagged_fish,
    tag_attr             = "",
    tag_release_platform = c("fishery", "1"),
    srv_sel = arrs$srv_sel,
    fish_sel = arrs$fish_sel,
    NAA = arrs$NAA,
    ty = 1,
    tseas = 1,
    tr = 1,
    n_pop = arrs$n_pop,
    n_ages = arrs$n_ages,
    n_sexes = arrs$n_sexes
  )
  # All tags should be in pop 2
  expect_equal(sum(result[1, , ]), 0,   tolerance = 1e-10)
  expect_equal(sum(result[2, , ]), 100, tolerance = 1e-10)
})

test_that("release_conv_tag_attr: platform = 'survey' uses srv_sel weights", {
  # Mirror of the fishery test but using survey selectivity
  arrs <- make_arrays(n_pop = 2, n_ages = 3, n_sexes = 2)
  arrs$srv_sel[1, , , , , , 1] <- 0
  arrs$srv_sel[2, , , , , , 1] <- 1

  tagged_fish <- array(100, dim = c(1, 1, 1))
  result <- release_conv_tag_attr(
    tagged_fish          = tagged_fish,
    tag_attr             = "",
    tag_release_platform = c("survey", "1"),
    srv_sel = arrs$srv_sel,
    fish_sel = arrs$fish_sel,
    NAA = arrs$NAA,
    ty = 1,
    tseas = 1,
    tr = 1,
    n_pop = arrs$n_pop,
    n_ages = arrs$n_ages,
    n_sexes = arrs$n_sexes
  )
  expect_equal(sum(result[1, , ]), 0,   tolerance = 1e-10)
  expect_equal(sum(result[2, , ]), 100, tolerance = 1e-10)
})

test_that("release_conv_tag_attr: non-uniform NAA allocates proportionally", {
  # n_pop = 1 so the only axis with variation is age
  arrs <- make_arrays(n_pop = 1, n_ages = 3, n_sexes = 1)
  # NAA has values 1, 2, 3 across ages at [pop=1, region=1, yr=1, seas=1]
  arrs$NAA[1, 1, 1, 1, 1, 1] <- 1
  arrs$NAA[1, 1, 1, 1, 2, 1] <- 2
  arrs$NAA[1, 1, 1, 1, 3, 1] <- 3

  total_tags <- 60
  tagged_fish <- array(total_tags, dim = c(1, 1, 1))
  result <- release_conv_tag_attr(
    tagged_fish          = tagged_fish,
    tag_attr             = "",
    tag_release_platform = c("population", NA),
    srv_sel = arrs$srv_sel,
    fish_sel = arrs$fish_sel,
    NAA = arrs$NAA,
    ty = 1,
    tseas = 1,
    tr = 1,
    n_pop = 1,
    n_ages = 3,
    n_sexes = 1
  )
  # Expected shares: 1/6, 2/6, 3/6 of total
  expect_equal(result[1, 1, 1], 60 * 1/6, tolerance = 1e-10)
  expect_equal(result[1, 2, 1], 60 * 2/6, tolerance = 1e-10)
  expect_equal(result[1, 3, 1], 60 * 3/6, tolerance = 1e-10)
})

test_that("release_conv_tag_attr: n_sexes = 1 does not error", {
  arrs <- make_arrays(n_pop = 2, n_ages = 3, n_sexes = 1)
  tagged_fish <- array(50, dim = c(1, 1, 1))
  result <- release_conv_tag_attr(
    tagged_fish          = tagged_fish,
    tag_attr             = "",
    tag_release_platform = c("population", NA),
    srv_sel = arrs$srv_sel,
    fish_sel = arrs$fish_sel,
    NAA = arrs$NAA,
    ty = 1,
    tseas = 1,
    tr = 1,
    n_pop = 2,
    n_ages = 3,
    n_sexes = 1
  )
  expect_equal(dim(result), c(2, 3, 1))
  expect_equal(sum(result), 50, tolerance = 1e-10)
})

test_that("release_conv_tag_attr: NAA = 0 for some cells yields 0 in those cells", {
  # If NAA is zero for pop 1 entirely, those cells should get 0 regardless of
  # what tag_attr is used (unattended pop dim)
  arrs <- make_arrays(n_pop = 2, n_ages = 3, n_sexes = 2)
  arrs$NAA[1, , , , , ] <- 0

  tagged_fish <- array(100, dim = c(1, 1, 1))
  result <- release_conv_tag_attr(
    tagged_fish          = tagged_fish,
    tag_attr             = "",
    tag_release_platform = c("population", NA),
    srv_sel = arrs$srv_sel,
    fish_sel = arrs$fish_sel,
    NAA = arrs$NAA,
    ty = 1,
    tseas = 1,
    tr = 1,
    n_pop = 2,
    n_ages = 3,
    n_sexes = 2
  )
  expect_true(all(result[1, , ] == 0))
  expect_equal(sum(result[2, , ]), 100, tolerance = 1e-10)
})

test_that("release_conv_tag_attr: all output values are non-negative", {
  arrs <- make_arrays(n_pop = 2, n_ages = 5, n_sexes = 2)
  tagged_fish <- array(200, dim = c(1, 1, 1))
  result <- release_conv_tag_attr(
    tagged_fish          = tagged_fish,
    tag_attr             = "",
    tag_release_platform = c("population", NA),
    srv_sel = arrs$srv_sel,
    fish_sel = arrs$fish_sel,
    NAA = arrs$NAA,
    ty = 1,
    tseas = 1,
    tr = 1,
    n_pop = 2,
    n_ages = 5,
    n_sexes = 2
  )
  expect_true(all(result >= 0))
})

test_that("release_conv_tag_attr: result is array with correct type", {
  arrs <- make_arrays()
  tagged_fish <- array(10, dim = c(1, 1, 1))
  result <- release_conv_tag_attr(
    tagged_fish          = tagged_fish,
    tag_attr             = "",
    tag_release_platform = c("population", NA),
    srv_sel = arrs$srv_sel,
    fish_sel = arrs$fish_sel,
    NAA = arrs$NAA,
    ty = 1,
    tseas = 1,
    tr = 1,
    n_pop = arrs$n_pop,
    n_ages = arrs$n_ages,
    n_sexes = arrs$n_sexes
  )
  expect_true(is.array(result))
  expect_true(is.numeric(result))
})

