# Checks peel_dsem_years rebuilds the linked cells rather than truncating them. Recruitment's year dim
# is last so its cells do not move; growth's is third so they do, and truncating would read the wrong ones.

peel_setup <- function(par_name, dims, yr_dim, idx, n_grid_yrs) {
  stride <- cumprod(c(1, dims[-length(dims)]))
  cell <- vapply(seq_len(dims[yr_dim]), function(y) { i <- idx; i[yr_dim] <- y; as.integer(1 + sum((i - 1) * stride)) }, 1L)
  vars <- c("x", "linked")
  model <- read_dsem_arrows(c("x -> linked, 0, b", "x <-> x, 0, sx", "linked <-> linked, 0, sl"), vars)
  data <- list(dsem_model = model, dsem_var_names = vars, dsem_n_grid_yrs = n_grid_yrs,
               dsem_cells = get_dsem_cells(model, n_grid_yrs),
               dsem_cov_obs = matrix(0, n_grid_yrs, 1),
               dsem_link_par = par_name, dsem_link_col = 2L, dsem_link_row = list(seq_len(dims[yr_dim])),
               dsem_link_cell = list(cell), dsem_link_idx = list(replace(idx, yr_dim, 0L)), dsem_link_yr_dim = yr_dim)
  pars <- list(dsem_x = matrix(0, n_grid_yrs, 2))
  pars[[par_name]] <- array(0, dim = dims)
  list(data = data, pars = pars, map = list(dsem_x = factor(seq_len(n_grid_yrs * 2))))
}

test_that("a peel to the same length changes nothing", {

  s <- peel_setup("ln_RecDevs", c(1, 3, 12), 3L, c(1L, 2L, 1L), 12)
  out <- peel_dsem_years(s$data, s$pars, s$map, 12)

  expect_equal(out$data$dsem_link_cell, s$data$dsem_link_cell)
  expect_equal(out$data$dsem_n_grid_yrs, 12)
  expect_equal(nrow(out$parameters$dsem_x), 12)

})

test_that("recruitment cells stay put under a peel, and the grid shrinks", {

  s <- peel_setup("ln_RecDevs", c(1, 3, 12), 3L, c(1L, 2L, 1L), 12)
  s$pars$ln_RecDevs <- array(0, dim = c(1, 3, 9)) # as the retrospective peels it
  out <- peel_dsem_years(s$data, s$pars, s$map, 9)

  # the year dim is last, so the first nine cells are the same numbers
  expect_equal(out$data$dsem_link_cell[[1]], s$data$dsem_link_cell[[1]][1:9])
  expect_equal(out$data$dsem_n_grid_yrs, 9)
  expect_equal(nrow(out$parameters$dsem_x), 9)
  expect_equal(out$data$dsem_cells$n_grid_yrs, 9)
  expect_equal(nrow(out$data$dsem_cov_obs), 9)

})

test_that("growth cells are rebuilt, and truncating them would be wrong", {

  dims <- c(1, 3, 12, 2, 2) # pop, region, year, par, sex
  idx <- c(1L, 2L, 1L, 2L, 2L)
  s <- peel_setup("ln_growth_devs", dims, 3L, idx, 12)
  s$pars$ln_growth_devs <- array(0, dim = c(1, 3, 9, 2, 2))
  out <- peel_dsem_years(s$data, s$pars, s$map, 9)

  # what a fresh build on the peeled array gives
  new_dims <- c(1, 3, 9, 2, 2)
  stride <- cumprod(c(1, new_dims[-length(new_dims)]))
  want <- vapply(1:9, function(y) { i <- idx; i[3] <- y; as.integer(1 + sum((i - 1) * stride)) }, 1L)

  expect_equal(out$data$dsem_link_cell[[1]], want)
  expect_false(identical(out$data$dsem_link_cell[[1]], s$data$dsem_link_cell[[1]][1:9])) # truncation differs
  expect_true(all(out$data$dsem_link_cell[[1]] <= prod(new_dims)))                        # and stays in bounds

})

test_that("a model with no dsem passes through untouched", {

  data <- list(years = 1:10)
  out <- peel_dsem_years(data, list(a = 1), list(), 8)
  expect_equal(out$data, data)
  expect_equal(out$parameters, list(a = 1))

})
