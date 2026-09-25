# Checks the diagonal variance form: an sd line is the marginal sd, so the innovation sd is the AR1 closed
# form, every cell's variance comes back at its sd line squared, the density is the dense MVN, and the refusals.

# I - B as a dense numeric matrix from the cell table, and the grid's covariance A diag(v) t(A) from per-cell sds
dense_IminusB <- function(dsem_model, dsem_cells) {
  arrow_value <- get_dsem_arrow_values(numeric(0), numeric(0), dsem_model)
  m <- matrix(0, dsem_cells$n_cells, dsem_cells$n_cells)
  m[cbind(dsem_cells$IminusB$entry_row, dsem_cells$IminusB$entry_col)] <- c(1, -arrow_value)[dsem_cells$IminusB$entry_arrow + 1]
  m
}
mvn_nLL <- function(x, mu, Sigma) {
  L <- t(chol(Sigma))
  z <- forwardsolve(L, as.vector(x - mu))
  sum(log(diag(L))) + 0.5 * sum(z^2) + 0.5 * length(x) * log(2 * pi)
}

test_that("an AR1 under the diagonal form has the closed form innovation sd and a flat marginal variance", {

  n_yrs <- 10
  rho <- 0.6
  s <- 0.8
  held <- hold_arrows(c("x -> x, 1, rho", "x <-> x, 0, s"), c(rho = rho, s = s))
  diag_model <- read_dsem_arrows(held, "x", variance = "diagonal")
  cond_model <- read_dsem_arrows(held, "x")
  expect_equal(cond_model$variance, "conditional")
  cells <- get_dsem_cells(diag_model, n_yrs)

  # year one takes the whole sd, every later year what the self path leaves
  parts <- get_dsem_matrices(numeric(0), numeric(0), diag_model, cells)
  expect_equal(as.numeric(parts$sd_cell), c(s, rep(s * sqrt(1 - rho^2), n_yrs - 1)), tolerance = 1e-12)
  expect_equal(as.numeric(get_dsem_matrices(numeric(0), numeric(0), cond_model, cells)$sd_cell), rep(s, n_yrs)) # the conditional form is untouched

  # with nothing known, every cell's marginal variance is the sd line squared; under the conditional form it settles at s^2 / (1 - rho^2)
  x <- matrix(0, n_yrs, 1)
  expect_equal(as.numeric(get_dsem_margvar(numeric(0), numeric(0), x, diag_model, cells, rep(FALSE, n_yrs))), rep(s^2, n_yrs), tolerance = 1e-12)
  cond_var <- as.numeric(get_dsem_margvar(numeric(0), numeric(0), x, cond_model, cells, rep(FALSE, n_yrs)))
  expect_equal(cond_var[1], s^2)
  expect_equal(cond_var[n_yrs], s^2 / (1 - rho^2), tolerance = 1e-4)

  # the density is the dense MVN built from those innovation sds
  set.seed(31)
  x <- matrix(rnorm(n_yrs), n_yrs, 1)
  mu <- matrix(0.2, n_yrs, 1)
  A <- solve(dense_IminusB(diag_model, cells))
  Sigma <- A %*% diag(as.numeric(parts$sd_cell)^2) %*% t(A)
  expect_equal(diag(Sigma), rep(s^2, n_yrs), tolerance = 1e-12)
  expect_equal(get_dsem_nLL(numeric(0), numeric(0), x, mu, diag_model, cells), mvn_nLL(x, mu, Sigma), tolerance = 1e-10)

})

test_that("with paths between series every cell still comes back at its own sd line, and the tape differentiates", {

  n_yrs <- 12
  arrows <- c("env -> rec, 0, b0", "env -> rec, 1, b1", "env -> env, 1, rho", "env <-> env, 0, sd_env", "rec <-> rec, 0, sd_rec")
  diag_model <- read_dsem_arrows(arrows, c("env", "rec"), variance = "diagonal")
  cells <- get_dsem_cells(diag_model, n_yrs)
  beta <- c(0.4, 0.25, 0.55)
  ln_sd <- log(c(0.9, 0.65))

  # every cell's marginal variance is its sd line squared, whatever feeds it
  x <- matrix(0, n_yrs, 2)
  mv <- get_dsem_margvar(beta, ln_sd, x, diag_model, cells, rep(FALSE, 2 * n_yrs))
  expect_equal(as.numeric(mv[,1]), rep(0.9^2, n_yrs), tolerance = 1e-12)
  expect_equal(as.numeric(mv[,2]), rep(0.65^2, n_yrs), tolerance = 1e-12)

  # the recruitment innovation shrinks as the paths carry more of the spread
  parts <- get_dsem_matrices(beta, ln_sd, diag_model, cells)
  sd_rec <- as.numeric(parts$sd_cell)[n_yrs + 1:n_yrs]
  expect_equal(sd_rec[1], sqrt(0.65^2 - 0.4^2 * 0.9^2), tolerance = 1e-12) # year one: only the lag-zero path into it
  expect_lt(sd_rec[n_yrs], sd_rec[1])

  # the density matches the dense MVN from the solved innovation variances
  set.seed(32)
  x <- matrix(rnorm(2 * n_yrs), n_yrs, 2)
  mu <- matrix(c(0.1, 0), n_yrs, 2, byrow = TRUE)
  m <- matrix(0, 2 * n_yrs, 2 * n_yrs)
  m[cbind(cells$IminusB$entry_row, cells$IminusB$entry_col)] <- c(1, -get_dsem_arrow_values(beta, ln_sd, diag_model))[cells$IminusB$entry_arrow + 1]
  A <- solve(m)
  Sigma <- A %*% diag(as.numeric(parts$sd_cell)^2) %*% t(A)
  expect_equal(diag(Sigma), rep(c(0.9^2, 0.65^2), each = n_yrs), tolerance = 1e-12)
  plain <- get_dsem_nLL(beta, ln_sd, x, mu, diag_model, cells)
  expect_equal(plain, mvn_nLL(x, mu, Sigma), tolerance = 1e-10)

  # on the tape, states included, with the gradient against finite differences
  f <- function(p) get_dsem_nLL(p$beta, p$ln_sd, p$x, mu, diag_model, cells)
  obj <- RTMB::MakeADFun(f, list(beta = beta, ln_sd = ln_sd, x = x), silent = TRUE)
  expect_equal(obj$fn(obj$par), plain, tolerance = 1e-12)
  expect_equal(as.vector(obj$gr(obj$par)), numDeriv::grad(obj$fn, obj$par), tolerance = 1e-6)

  # and the precision the simulator draws from has that covariance
  Q <- get_dsem_precision(beta, ln_sd, diag_model, cells)
  expect_equal(as.matrix(solve(Q)), Sigma, tolerance = 1e-10, ignore_attr = TRUE)

})

test_that("the default is bit-identical to before, and the diagonal form refuses what it cannot solve", {

  case <- dsem_test_cases$two_series
  n_yrs <- 12
  m_default <- read_dsem_arrows(case$arrows, case$variables)
  m_cond <- read_dsem_arrows(case$arrows, case$variables, variance = "conditional")
  cells <- get_dsem_cells(m_default, n_yrs)
  pars <- dsem_case_pars(m_default, case$values)
  set.seed(33)
  x <- matrix(rnorm(n_yrs * 2), n_yrs, 2)
  mu <- matrix(0, n_yrs, 2)
  expect_identical(get_dsem_nLL(pars$dsem_beta, pars$ln_dsem_sd, x, mu, m_default, cells), get_dsem_nLL(pars$dsem_beta, pars$ln_dsem_sd, x, mu, m_cond, cells))

  expect_error(read_dsem_arrows("x <-> x, 0, s", "x", variance = "full"), "conditional")
  expect_error(read_dsem_arrows(c("x -> y, 0, NA, 2", "x -> x, 1, rx", "x <-> x, 0, sx", "y <-> y, 0, NA, 0"), c("x", "y"), variance = "diagonal"), "sd of zero")
  expect_error(read_dsem_arrows(c("r -> r, 1, rho_r", "r <-> r, 0, sd_r", "y <-> y, 0, r"), c("r", "y"), variance = "diagonal"), "moderated sd")

  # a random walk carries every cell past its sd line from year two on, which setup catches at the start
  plain <- suppressMessages(sweep_input(rec = list(RecDevs_model = "dsem"), dims = list(n_regions = 1)))
  expect_error(suppressMessages(Setup_Mod_DSEM(plain, c("rec -> rec, 1, NA, 1", "rec <-> rec, 0, sd_rec, 0.7"), NULL, dsem_variance = "diagonal")), "no innovation variance")

})

test_that("in a fit the linked recruitment correction is half the sd line squared under the diagonal form", {

  # a latent factor on an AR1 driving recruitment. under the conditional form the recruitment cell's
  # variance is b^2 var(F1) + sd_rec^2; under the diagonal form it is sd_rec^2 by construction
  plain <- suppressMessages(sweep_input(rec = list(RecDevs_model = "dsem"), dims = list(n_regions = 1)))
  n_yrs <- length(plain$data$years)
  latent <- data.frame(year = plain$data$years, F1 = NA_real_)
  arrows <- c("F1 -> rec, 0, NA, 0.5", "F1 -> F1, 1, NA, 0.6", "F1 <-> F1, 0, NA, 1", "rec <-> rec, 0, NA, 0.7")
  value_of <- function(il) {
    obj <- fit_model(il$data, il$par, il$map, random = NULL, do_optim = FALSE, silent = TRUE)
    list(fn = as.numeric(obj$fn(obj$par)), margvar = obj$rep$dsem_margvar_grid)
  }

  cond <- suppressMessages(Setup_Mod_DSEM(plain, arrows, latent, dsem_mu_spec = "fix"))
  diag <- suppressMessages(Setup_Mod_DSEM(plain, arrows, latent, dsem_mu_spec = "fix", dsem_variance = "diagonal"))
  expect_equal(diag$data$dsem_model$variance, "diagonal")
  out_cond <- value_of(cond)
  out_diag <- value_of(diag)
  expect_true(is.finite(out_diag$fn))
  expect_false(isTRUE(all.equal(out_diag$fn, out_cond$fn)))

  expect_equal(as.numeric(out_diag$margvar[,2]), rep(0.7^2, n_yrs), tolerance = 1e-10)
  expect_equal(as.numeric(out_cond$margvar[1,2]), 0.5^2 * 1 + 0.7^2, tolerance = 1e-10) # year one: the factor's own sd
  expect_equal(as.numeric(out_cond$margvar[n_yrs,2]), 0.5^2 / (1 - 0.6^2) + 0.7^2, tolerance = 1e-3) # settled

})
