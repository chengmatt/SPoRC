# Checks the marginal variance form: without covariance lines it is the diagonal form to the bit, with them
# every cell lands on its sd line squared while the diagonal form only lands near it, and the tape differentiates.

# the grid's covariance under the innovations the builder settled on, from the same parts the density reads
grid_covariance <- function(parts) {
  IminusB <- as.matrix(Matrix::sparseMatrix(i = parts$IminusB@i + 1L, p = parts$IminusB@p, x = as.numeric(parts$IminusB@x), dims = dim(parts$IminusB)))
  V <- solve(as.matrix(Matrix::sparseMatrix(i = parts$Vinv@i + 1L, p = parts$Vinv@p, x = as.numeric(parts$Vinv@x), dims = dim(parts$Vinv))))
  A <- solve(IminusB)
  list(Sigma = A %*% V %*% t(A), V = V)
}

test_that("without covariance lines the marginal form is the diagonal form", {

  n_yrs <- 12
  arrows <- c("env -> rec, 0, b0", "env -> rec, 1, b1", "env -> env, 1, rho", "env <-> env, 0, sd_env", "rec <-> rec, 0, sd_rec")
  m_diag <- read_dsem_arrows(arrows, c("env", "rec"), variance = "diagonal")
  m_marg <- read_dsem_arrows(arrows, c("env", "rec"), variance = "marginal")
  cells <- get_dsem_cells(m_diag, n_yrs)
  beta <- c(0.4, 0.25, 0.55)
  ln_sd <- log(c(0.9, 0.65))
  set.seed(41)
  x <- matrix(rnorm(2 * n_yrs), n_yrs, 2)
  mu <- matrix(0, n_yrs, 2)

  expect_identical(get_dsem_nLL(beta, ln_sd, x, mu, m_marg, cells), get_dsem_nLL(beta, ln_sd, x, mu, m_diag, cells))
  expect_identical(as.numeric(get_dsem_matrices(beta, ln_sd, m_marg, cells)$sd_cell), as.numeric(get_dsem_matrices(beta, ln_sd, m_diag, cells)$sd_cell))

})

test_that("with a covariance line the marginal form lands every cell on its sd line and keeps the innovation correlation", {

  n_yrs <- 10
  case <- dsem_test_cases$covariance # a -> b at lag 1, a an AR1, and a <-> b
  m_cond <- read_dsem_arrows(case$arrows, case$variables)
  m_diag <- read_dsem_arrows(case$arrows, case$variables, variance = "diagonal")
  m_marg <- read_dsem_arrows(case$arrows, case$variables, variance = "marginal")
  cells <- get_dsem_cells(m_cond, n_yrs)
  pars <- dsem_case_pars(m_cond, case$values)
  target <- rep(c(case$values[["sd_a"]], case$values[["sd_b"]])^2, each = n_yrs)

  cond <- grid_covariance(get_dsem_matrices(pars$dsem_beta, pars$ln_dsem_sd, m_cond, cells))
  diag_form <- grid_covariance(get_dsem_matrices(pars$dsem_beta, pars$ln_dsem_sd, m_diag, cells))
  marg <- grid_covariance(get_dsem_matrices(pars$dsem_beta, pars$ln_dsem_sd, m_marg, cells))

  # the conditional form drifts off the sd lines, the diagonal form gets close, the marginal form lands on them
  expect_gt(max(abs(diag(cond$Sigma) - target) / target), 0.1)
  miss_diag <- max(abs(diag(diag_form$Sigma) - target) / target)
  expect_gt(miss_diag, 1e-6)
  expect_lt(miss_diag, 0.2)
  expect_equal(diag(marg$Sigma), target, tolerance = 1e-8, ignore_attr = TRUE)

  # the innovation correlation the lines wrote is the one the marginal form keeps; the diagonal form bends it
  corr_of <- function(V) stats::cov2cor(V)
  expect_equal(corr_of(marg$V), corr_of(cond$V), tolerance = 1e-8, ignore_attr = TRUE)
  expect_gt(max(abs(corr_of(diag_form$V) - corr_of(cond$V))), 1e-3)

  # the density is the dense MVN with that covariance, and the tape agrees with finite differences
  set.seed(42)
  x <- matrix(rnorm(2 * n_yrs), n_yrs, 2)
  mu <- matrix(0.1, n_yrs, 2)
  L <- t(chol(marg$Sigma))
  z <- forwardsolve(L, as.vector(x - mu))
  by_hand <- sum(log(diag(L))) + 0.5 * sum(z^2) + 0.5 * length(x) * log(2 * pi)
  plain <- get_dsem_nLL(pars$dsem_beta, pars$ln_dsem_sd, x, mu, m_marg, cells)
  expect_equal(plain, by_hand, tolerance = 1e-8)
  f <- function(p) get_dsem_nLL(p$dsem_beta, p$ln_dsem_sd, p$x, mu, m_marg, cells)
  obj <- RTMB::MakeADFun(f, list(dsem_beta = pars$dsem_beta, ln_dsem_sd = pars$ln_dsem_sd, x = x), silent = TRUE)
  expect_equal(obj$fn(obj$par), plain, tolerance = 1e-10)
  expect_equal(as.vector(obj$gr(obj$par)), numDeriv::grad(obj$fn, obj$par), tolerance = 1e-6)

  # and the correction reads the same settled variances
  mv <- get_dsem_margvar(pars$dsem_beta, pars$ln_dsem_sd, matrix(0, n_yrs, 2), m_marg, cells, rep(FALSE, 2 * n_yrs))
  expect_equal(as.numeric(mv), target, tolerance = 1e-8)

})

test_that("Setup_Mod_DSEM takes the marginal form and refuses what it cannot settle", {

  plain <- suppressMessages(sweep_input(rec = list(RecDevs_model = "dsem"), dims = list(n_regions = 1)))
  n_yrs <- length(plain$data$years)
  latent <- data.frame(year = plain$data$years, F1 = NA_real_)
  arrows <- c("F1 -> rec, 0, NA, 0.5", "F1 -> F1, 1, NA, 0.6", "F1 <-> F1, 0, NA, 1", "rec <-> rec, 0, NA, 0.7", "F1 <-> rec, 0, NA, 0.2")

  marg <- suppressMessages(Setup_Mod_DSEM(plain, arrows, latent, dsem_mu_spec = "fix", dsem_variance = "marginal"))
  expect_equal(marg$data$dsem_model$variance, "marginal")
  obj <- fit_model(marg$data, marg$par, marg$map, random = NULL, do_optim = FALSE, silent = TRUE)
  expect_true(is.finite(obj$fn(obj$par)))
  expect_equal(as.numeric(obj$rep$dsem_margvar_grid[,2]), rep(0.7^2, n_yrs), tolerance = 1e-8) # on the sd line, covariance line and all

  # a random walk is refused under either form
  expect_error(suppressMessages(Setup_Mod_DSEM(plain, c("rec -> rec, 1, NA, 1", "rec <-> rec, 0, sd_rec, 0.7"), NULL, dsem_variance = "marginal")), "no innovation variance")

})
