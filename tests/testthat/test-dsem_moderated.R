# Checks arrows named after a series, so a path coefficient or an sd changes year to year: against a
# dense covariance (1e-9), against dsem's C++ (1e-10), on the tape, and the refusals.

test_that("a moderated arrow matches the dense covariance and tapes", {

  set.seed(20)
  for(case_name in names(dsem_moderated_cases)) {

    case <- dsem_moderated_cases[[case_name]]
    n_v <- length(case$variables)
    dsem_model <- read_dsem_arrows(case$arrows, case$variables, mod_var_logscale = case$mod_var_logscale)
    x_grid <- matrix(rnorm(15 * n_v), 15, n_v)
    mu_grid <- matrix(0.2, 15, n_v)
    dsem_cells <- get_dsem_cells(dsem_model, 15)
    case_pars <- dsem_case_pars(dsem_model, case$values)

    plain <- get_dsem_nLL(case_pars$dsem_beta, case_pars$ln_dsem_sd, x_grid, mu_grid, dsem_model, dsem_cells)
    dense <- dense_dsem_nLL(case_pars$dsem_beta, case_pars$ln_dsem_sd, x_grid, mu_grid, dsem_model)
    expect_equal(plain, dense, tolerance = 1e-9, label = case_name)

    # the grid is a parameter here, so the gradient runs through the moderated entries as well
    f <- function(p) get_dsem_nLL(p$dsem_beta, p$ln_dsem_sd, p$x_grid, mu_grid, dsem_model, dsem_cells)
    obj <- RTMB::MakeADFun(f, list(dsem_beta = case_pars$dsem_beta, ln_dsem_sd = case_pars$ln_dsem_sd, x_grid = x_grid), silent = TRUE)
    expect_equal(obj$fn(obj$par), plain, tolerance = 1e-12, label = case_name)
    expect_equal(as.vector(obj$gr(obj$par)), numDeriv::grad(obj$fn, obj$par), tolerance = 1e-6, label = case_name)

  } # end case_name loop

})

test_that("a moderated arrow matches dsem's C++ template", {

  skip_if_not_installed("dsem")

  set.seed(21)
  for(case_name in names(dsem_moderated_cases)) {

    case <- dsem_moderated_cases[[case_name]]
    n_v <- length(case$variables)
    dsem_model <- read_dsem_arrows(case$arrows, case$variables, mod_var_logscale = case$mod_var_logscale)
    x_grid <- matrix(rnorm(25 * n_v), 25, n_v, dimnames = list(NULL, case$variables))
    case_pars <- dsem_case_pars(dsem_model, case$values)

    y_na <- stats::ts(x_grid)
    y_na[] <- NA
    family <- stats::setNames(rep(list(dsem::fixed()), n_v), case$variables)
    fit <- dsem::dsem(sem = case$arrows, tsdata = y_na, family = family,
                      control = dsem::dsem_control(run_model = FALSE, quiet = TRUE,
                                                   gmrf_parameterization = "full",
                                                   logscale_moderating_variance = case$mod_var_logscale))
    cpp_par <- fit$obj$env$last.par
    cpp_par[names(cpp_par) == "beta_z"] <- case$values[case$cpp_pars]
    cpp_par[names(cpp_par) == "x_tj"] <- as.vector(x_grid)

    ours <- get_dsem_nLL(case_pars$dsem_beta, case_pars$ln_dsem_sd, x_grid, 0 * x_grid, dsem_model, get_dsem_cells(dsem_model, 25))
    expect_equal(ours, fit$obj$env$f(cpp_par, order = 0), tolerance = 1e-10, label = case_name)

  } # end case_name loop

})

test_that("a moderated coefficient is the moderating series, year by year", {

  case <- dsem_moderated_cases$moderated_path
  dsem_model <- read_dsem_arrows(case$arrows, case$variables, mod_var_logscale = FALSE)
  dsem_cells <- get_dsem_cells(dsem_model, 8)
  case_pars <- dsem_case_pars(dsem_model, case$values)

  set.seed(22)
  x_grid <- matrix(rnorm(8 * 3), 8, 3)
  Q <- get_dsem_precision(case_pars$dsem_beta, case_pars$ln_dsem_sd, dsem_model, dsem_cells, x_grid)
  Sigma <- solve(as.matrix(Q))

  # recruitment given the environment is normal with mean b_t times env_t, so the conditional sd is sd_rec
  # and the covariance of rec with env in the same year is b_t times the variance of env
  cond_sd <- 1 / sqrt(diag(Q)[16 + 1:8]) # the rec block is third in the cell numbering
  expect_equal(unname(cond_sd), rep(0.7, 8), tolerance = 1e-10)
  expect_equal(unname(Sigma[cbind(16 + 1:8, 1:8)]), unname(x_grid[,2] * Sigma[cbind(1:8, 1:8)]), tolerance = 1e-9)

  # the same model with a fixed coefficient is what a constant moderating series gives
  flat <- x_grid
  flat[,2] <- 0.45
  Q_flat <- get_dsem_precision(case_pars$dsem_beta, case_pars$ln_dsem_sd, dsem_model, dsem_cells, flat)
  fixed_model <- read_dsem_arrows(sub("env -> rec, 0, b", "env -> rec, 0, NA, 0.45", case$arrows, fixed = TRUE), case$variables)
  Q_fixed <- get_dsem_precision(dsem_case_pars(fixed_model, case$values)$dsem_beta,
                                dsem_case_pars(fixed_model, case$values)$ln_dsem_sd, fixed_model, get_dsem_cells(fixed_model, 8))
  expect_equal(as.matrix(Q_flat), as.matrix(Q_fixed), tolerance = 1e-12)

})

test_that("the reader refuses a moderator the arrow sets itself", {

  v <- c("a", "b", "c")
  sds <- c("a <-> a, 0, sa", "b <-> b, 0, sb", "c <-> c, 0, sc")

  # c moderates the arrow into c, so its value would set itself in the same year
  expect_error(read_dsem_arrows(c("a -> c, 0, c", sds), v), "sets in the same year")

  # b moderates the arrow into c, and c affects b in the same year
  expect_error(read_dsem_arrows(c("a -> c, 0, b", "c -> b, 0, p", sds), v), "sets in the same year")

  # a loop the same-year paths make on their own is read, since dgmrf works out the determinant
  loop_sds <- c("a <-> a, 0, sa", "b <-> b, 0, sb", "m <-> m, 0, sm")
  loop <- read_dsem_arrows(c("a -> b, 0, m", "b -> a, 0, p", loop_sds), c("a", "b", "m"))
  expect_null(loop$series_order)
  expect_false(get_dsem_cells(loop, 6)$det_is_one)
  expect_true(get_dsem_cells(loop, 6)$needs_dense_logdet)

  # the same moderator with the path running the other way is read
  m <- read_dsem_arrows(c("a -> c, 0, b", "b -> c, 0, p", sds), v)
  expect_equal(m$arrows$mod_idx, c(2L, 0L, 0L, 0L, 0L))
  expect_equal(m$beta_names, "p")

  # a lagged moderated arrow still reads its moderator in the year it points to
  m_lag <- read_dsem_arrows(c("a -> c, 2, b", sds), v)
  expect_equal(m_lag$arrows$mod_idx[1], 2L)
  expect_true(get_dsem_cells(m_lag, 6)$has_mod)

})

test_that("a series with an sd of zero is projected off what points into it", {

  variables <- c("ones", "logX", "logX1", "logX2", "X")
  arrows <- c("ones -> ones, 1, NA, 1", "ones -> logX, 0, alpha", "logX -> logX, 1, NA, 1
             logX -> logX1, 0, NA, 1", "logX1 -> logX2, 0, logX
             ones -> X, 0, NA, 1", "logX1 -> X, 0, NA, 1", "logX2 -> X, 0, NA, 0.5
             ones <-> ones, 0, NA, 0.001", "logX <-> logX, 0, sd_logX
             logX1 <-> logX1, 0, NA, 0", "logX2 <-> logX2, 0, NA, 0", "X <-> X, 0, NA, 0")
  dsem_model <- read_dsem_arrows(arrows, variables)
  expect_equal(variables[dsem_model$project_k], c("logX1", "logX2", "X"))

  n_t <- 10
  x_grid <- matrix(0, n_t, length(variables))
  x_grid[,1] <- 1
  x_grid[,2] <- seq(-0.6, 0.6, length.out = n_t)
  mu_grid <- matrix(0, n_t, length(variables))
  dsem_cells <- get_dsem_cells(dsem_model, n_t)
  filled <- get_dsem_grid(0.05, log(0.3), x_grid, mu_grid, dsem_model, dsem_cells)$x_grid

  # the chain is the quadratic approximation to exp, so it agrees to the term left out
  expect_equal(filled[,5], 1 + x_grid[,2] + 0.5 * x_grid[,2]^2, tolerance = 1e-12)
  expect_lt(max(abs(filled[,5] - exp(x_grid[,2]))), 0.05)

  # the projected cells owe nothing to the density, which is the two series that keep an innovation
  ours <- get_dsem_nLL(0.05, log(0.3), x_grid, mu_grid, dsem_model, dsem_cells)
  by_hand <- -dnorm(x_grid[1,1], 0, 0.001, TRUE) - sum(dnorm(diff(x_grid[,1]), 0, 0.001, TRUE)) -
    dnorm(x_grid[1,2] - 0.05 * x_grid[1,1], 0, 0.3, TRUE) -
    sum(dnorm(diff(x_grid[,2]) - 0.05 * x_grid[-1,1], 0, 0.3, TRUE))
  expect_equal(ours, by_hand, tolerance = 1e-8)

  expect_error(read_dsem_arrows(sub("logX -> logX1, 0, NA, 1", "", arrows, fixed = TRUE), variables), "nothing sets their value")

})
