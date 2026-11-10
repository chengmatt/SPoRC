# A series with an sd of zero is solved out with sparse matrices throughout: the two flags setup works out, the
# precision against dense blocks by hand (1e-10), and stored entries that grow with years, not years squared.

# the pieces of one model at arrow values given by name
sparse_rr_pieces <- function(sem, vars, n_t, values, flags = TRUE) {
  m <- read_dsem_arrows(sem, vars)
  cells <- get_dsem_cells(m, n_t)
  if(!flags) cells$solve_uu_identity <- cells$solve_has_ou <- NULL # a cell list built before the flags existed
  x <- matrix(stats::rnorm(n_t * length(vars), 0, 0.5), n_t, length(vars))
  mu <- matrix(0, n_t, length(vars))
  beta <- if(length(m$beta_names)) unname(values[m$beta_names]) else numeric(0)
  ln_sd <- if(length(m$ln_sd_names)) log(unname(values[m$ln_sd_names])) else numeric(0)
  grid <- get_dsem_grid(beta, ln_sd, x, mu, m, cells)
  Q <- get_dsem_Q_oo(grid$parts$IminusB, grid$parts, grid$solve_mat, cells)
  list(model = m, cells = cells, grid = grid, Q = Q, x = x, mu = mu,
       nll = as.numeric(get_dsem_nLL(beta, ln_sd, x, mu, m, cells)))
}

# the same precision from dense blocks
dense_rr_Q <- function(o) {
  ImB <- as.matrix(o$grid$parts$IminusB)
  u <- o$cells$unobs_idx
  ob <- o$cells$obs_idx
  V <- if(is.null(o$grid$parts$sd_cell)) as.matrix(o$grid$parts$V) else diag(o$grid$parts$sd_cell^2)
  S <- ImB[ob,ob,drop = FALSE] - ImB[ob,u,drop = FALSE] %*% solve(ImB[u,u,drop = FALSE], ImB[u,ob,drop = FALSE])
  t(S) %*% solve(V[ob,ob,drop = FALSE]) %*% S
}

rr_values <- c(sd_x = 0.9, rho = 0.4, b_xy = 0.7, b_yz = -0.5, sd_z = 0.6, c_xz = 0.3, sd_w = 1.1)

# y has no innovation: set by x alone, then feeding z, then with a path of its own from the year before
sem_leaf <- c("x <-> x, 0, sd_x", "x -> x, 1, rho", "x -> y, 0, b_xy", "y <-> y, 0, NA, 0", "z <-> z, 0, sd_z")
sem_feeds <- c(sem_leaf, "y -> z, 0, b_yz")
sem_own_lag <- c(sem_feeds, "y -> y, 1, NA, 0.5")
rr_vars <- c("x", "y", "z")

test_that("setup flags when the solved cells need no inversion and when they feed nothing", {
  leaf <- get_dsem_cells(read_dsem_arrows(sem_leaf, rr_vars), 6)
  feeds <- get_dsem_cells(read_dsem_arrows(sem_feeds, rr_vars), 6)
  own_lag <- get_dsem_cells(read_dsem_arrows(sem_own_lag, rr_vars), 6)

  expect_true(leaf$solve_uu_identity)
  expect_false(leaf$solve_has_ou)
  expect_true(feeds$solve_uu_identity)
  expect_true(feeds$solve_has_ou)
  expect_false(own_lag$solve_uu_identity)
  expect_true(own_lag$solve_has_ou)

  # nothing solved out, so nothing to flag either way
  none <- get_dsem_cells(read_dsem_arrows(c("x <-> x, 0, sd_x", "z <-> z, 0, sd_z"), c("x", "z")), 6)
  expect_length(none$unobs_idx, 0)
})

test_that("the sparse precision is the dense one, with and without a covariance line", {
  set.seed(3)
  n_t <- 7

  for(sem in list(sem_leaf, sem_feeds, sem_own_lag)) {
    for(with_cov in c(FALSE, TRUE)) {

      arrows <- if(with_cov) c(sem, "x <-> z, 0, c_xz") else sem
      o <- sparse_rr_pieces(arrows, rr_vars, n_t, rr_values)
      by_hand <- dense_rr_Q(o)

      expect_s4_class(o$Q, "sparseMatrix")
      expect_equal(unname(as.matrix(o$Q)), unname(by_hand), tolerance = 1e-10)

      # the density is the normal one at that precision
      off <- as.vector(o$x - o$mu)[o$cells$obs_idx]
      nll_by_hand <- -0.5 * as.numeric(determinant(by_hand)$modulus) + 0.5 * sum(off * (by_hand %*% off)) + 0.5 * length(off) * log(2 * pi)
      expect_equal(o$nll, nll_by_hand, tolerance = 1e-9)

      # a cell list without the flags takes the longer way to the same matrix
      set.seed(3)
      no_flags <- sparse_rr_pieces(arrows, rr_vars, n_t, rr_values, flags = FALSE)
      expect_equal(unname(as.matrix(no_flags$Q)), unname(as.matrix(sparse_rr_pieces(arrows, rr_vars, n_t, rr_values)$Q)), tolerance = 1e-12)

    } # end with_cov loop
  } # end sem loop
})

test_that("stored entries grow with the number of years, also under a covariance line", {
  # a fourth series w tied to nothing, so a covariance line elsewhere must not fill its rows
  sem <- c(sem_feeds, "x <-> z, 0, c_xz", "w <-> w, 0, sd_w")
  vars <- c(rr_vars, "w")
  stored_20 <- Matrix::nnzero(sparse_rr_pieces(sem, vars, 20, rr_values)$Q)
  stored_40 <- Matrix::nnzero(sparse_rr_pieces(sem, vars, 40, rr_values)$Q)

  expect_lt(stored_40, 2.2 * stored_20) # twice the years, about twice the entries
  expect_lt(stored_40, 0.1 * (3 * 40)^2) # far from every pair of the 120 cells that keep an innovation
})

test_that("the sparse route differentiates", {
  m <- read_dsem_arrows(c(sem_feeds, "x <-> z, 0, c_xz"), rr_vars)
  cells <- get_dsem_cells(m, 6)
  set.seed(5)
  x <- matrix(stats::rnorm(18, 0, 0.5), 6, 3)
  mu <- matrix(0, 6, 3)
  start <- c(rr_values[m$beta_names], log(rr_values[m$ln_sd_names]))
  n_beta <- length(m$beta_names)

  nll <- function(p) get_dsem_nLL(p[seq_len(n_beta)], p[-seq_len(n_beta)], x, mu, m, cells)
  tape <- RTMB::MakeTape(nll, start)
  by_differences <- numeric(length(start))
  for(i in seq_along(start)) {
    up <- dn <- start
    up[i] <- start[i] + 1e-6
    dn[i] <- start[i] - 1e-6
    by_differences[i] <- (as.numeric(nll(up)) - as.numeric(nll(dn))) / 2e-6
  } # end i loop
  expect_equal(as.numeric(tape$jacobian(start)), by_differences, tolerance = 1e-5)
})
