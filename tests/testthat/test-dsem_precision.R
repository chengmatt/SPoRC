# Checks the dgmrf density against a dense covariance built from the arrows (1e-9), on the tape as in
# plain R, gradients against finite differences (1e-6), and the arrow reader's refusals.

test_that("dsem density matches the dense covariance for every arrow type", {

  set.seed(10)
  for(case_name in names(dsem_test_cases)) {

    case <- dsem_test_cases[[case_name]]
    dsem_model <- read_dsem_arrows(case$arrows, case$variables)
    x_grid <- matrix(rnorm(12 * length(case$variables)), 12, length(case$variables))
    mu_grid <- matrix(0.3, 12, length(case$variables))
    dsem_cells <- get_dsem_cells(dsem_model, 12)
    case_pars <- dsem_case_pars(dsem_model, case$values)

    plain <- get_dsem_nLL(case_pars$dsem_beta, case_pars$ln_dsem_sd, x_grid, mu_grid, dsem_model, dsem_cells)
    dense <- dense_dsem_nLL(case_pars$dsem_beta, case_pars$ln_dsem_sd, x_grid, mu_grid, dsem_model)
    expect_equal(plain, dense, tolerance = 1e-9, label = case_name)

    # on the tape, with the states as parameters too
    f <- function(p) get_dsem_nLL(p$dsem_beta, p$ln_dsem_sd, p$x_grid, mu_grid, dsem_model, dsem_cells)
    obj <- RTMB::MakeADFun(f, list(dsem_beta = case_pars$dsem_beta, ln_dsem_sd = case_pars$ln_dsem_sd, x_grid = x_grid), silent = TRUE)
    expect_equal(obj$fn(obj$par), plain, tolerance = 1e-12, label = case_name)
    expect_equal(as.vector(obj$gr(obj$par)), numDeriv::grad(obj$fn, obj$par), tolerance = 1e-6, label = case_name)

  } # end case_name loop

})

test_that("dsem density matches dsem's C++ template", {

  skip_if_not_installed("dsem")
  case <- dsem_test_cases$two_series
  dsem_model <- read_dsem_arrows(case$arrows, case$variables)
  set.seed(11)
  x_grid <- matrix(rnorm(30 * 2), 30, 2, dimnames = list(NULL, case$variables))
  case_pars <- dsem_case_pars(dsem_model, case$values)

  y_na <- stats::ts(x_grid)
  y_na[] <- NA
  fit <- dsem::dsem(sem = case$arrows, tsdata = y_na, family = list(env = dsem::fixed(), rec = dsem::fixed()),
                    control = dsem::dsem_control(run_model = FALSE, quiet = TRUE))
  cpp_par <- fit$obj$env$last.par
  cpp_par[names(cpp_par) == "beta_z"] <- case$values[c("b0", "b1", "rho", "sd_env", "sd_rec")]
  cpp_par[names(cpp_par) == "x_tj"] <- as.vector(x_grid)

  ours <- get_dsem_nLL(case_pars$dsem_beta, case_pars$ln_dsem_sd, x_grid, 0 * x_grid, dsem_model, get_dsem_cells(dsem_model, 30))
  expect_equal(ours, fit$obj$env$f(cpp_par, order = 0), tolerance = 1e-10)

})

test_that("arrow reader refuses what it cannot read", {

  v <- c("a", "b")
  expect_error(read_dsem_arrows(c("a -> b, 0, p", "a <-> a, 1, sa", "b <-> b, 0, sb"), v), "cannot be lagged")
  expect_equal(read_dsem_arrows(c("a -> b, 0, p", "a <-> a, 0, sa"), v)$ln_sd_names, c("sa", "V[b]")) # the missing sd line is added
  expect_error(read_dsem_arrows(c("a -> a, 0, p", "a <-> a, 0, sa", "b <-> b, 0, sb"), v), "same year")
  expect_error(read_dsem_arrows(c("a -> b, 1, p", "a -> b, 1, q", "a <-> a, 0, sa", "b <-> b, 0, sb"), v), "given twice")
  expect_error(read_dsem_arrows(c("a -> b, 1, NA", "a <-> a, 0, sa", "b <-> b, 0, sb"), v), "needs a value")
  expect_error(read_dsem_arrows(c("a -> b, 1, s", "a <-> a, 0, s", "b <-> b, 0, sb"), v), "both an sd")
  expect_error(read_dsem_arrows(c("a -> z, 1, p", "a <-> a, 0, sa", "b <-> b, 0, sb"), v), "not defined")

  # comments, blank lines and a shared name are read
  m <- read_dsem_arrows(c("# header", "a -> b, 1, p  # lagged", "b -> b, 1, p", "a <-> a, 0, sd", "b <-> b, 0, sd"), v)
  expect_equal(m$beta_names, "p")
  expect_equal(m$ln_sd_names, "sd")
  expect_equal(m$arrows$par, c(1L, 1L, 1L, 1L))

})

test_that("leads, covs groups and automatic sd lines read as dsem reads them", {

  skip_if_not_installed("dsem")

  variables <- c("a", "b")
  n_t <- 20
  set.seed(12)
  x_grid <- matrix(rnorm(n_t * 2), n_t, 2, dimnames = list(NULL, variables))

  cases <- list(
    lead_only = list(arrows = c("a -> b, -1, p", "a -> a, 1, r", "a <-> a, 0, sa", "b <-> b, 0, sb"),
                     values = c(p = 0.5, r = 0.4, sa = 0.9, sb = 0.6)),
    lag_and_lead = list(arrows = c("a -> a, 1, r", "a -> a, -1, r2", "b -> b, 1, rb", "a <-> a, 0, sa", "b <-> b, 0, sb"),
                        values = c(r = 0.3, r2 = 0.2, rb = 0.5, sa = 0.9, sb = 0.6)),
    covs_group = list(arrows = c("a -> b, 1, p", "a -> a, 1, r"), covs = "a, b",
                      values = c(p = 0.5, r = 0.4, `C[a,b]` = 0.3, `V[a]` = 0.9, `V[b]` = 0.6)))

  for(case_name in names(cases)) {

    case <- cases[[case_name]]
    dsem_model <- read_dsem_arrows(case$arrows, variables, covs = case$covs)
    ours <- get_dsem_nLL(unname(case$values[dsem_model$beta_names]), log(unname(case$values[dsem_model$ln_sd_names])),
                         x_grid, 0 * x_grid, dsem_model, get_dsem_cells(dsem_model, n_t))

    y_na <- stats::ts(x_grid)
    y_na[] <- NA
    fit <- dsem::dsem(sem = case$arrows, tsdata = y_na, covs = case$covs,
                      family = list(a = dsem::fixed(), b = dsem::fixed()),
                      control = dsem::dsem_control(run_model = FALSE, quiet = TRUE, gmrf_parameterization = "full"))
    ram <- dsem::make_dsem_ram(case$arrows, times = 1:n_t, variables = variables, covs = case$covs, quiet = TRUE)
    cpp_par <- fit$obj$env$last.par
    cpp_par[names(cpp_par) == "beta_z"] <- case$values[unique(stats::na.omit(ifelse(ram$model$parameter > 0, ram$model$name, NA)))]
    cpp_par[names(cpp_par) == "x_tj"] <- as.vector(x_grid)

    expect_equal(ours, fit$obj$env$f(cpp_par, order = 0), tolerance = 1e-10, label = case_name)

  } # end case_name loop

  # a lead alone still orders the cells, a lead against a lag does not, and the innovation route
  # needs that ordering, so it is refused rather than giving the wrong number
  lead <- read_dsem_arrows(c("a -> b, -1, p", "a -> a, 1, r", "a <-> a, 0, sa", "b <-> b, 0, sb"), variables)
  both <- read_dsem_arrows(c("a -> a, 1, r", "a -> a, -1, r2", "a <-> a, 0, sa", "b <-> b, 0, sb"), variables)
  expect_true(get_dsem_cells(lead, 8)$det_is_one)
  expect_false(get_dsem_cells(both, 8)$det_is_one)
  expect_true(get_dsem_cells(lead, 8)$has_lead)

  # a series with no sd line of its own gets one, named the way dsem names it
  filled <- read_dsem_arrows("a -> b, 1, p", variables)
  expect_equal(filled$arrows$name[filled$arrows$type == "sd"], c("V[a]", "V[b]"))

})

test_that("a first year offset is propagated through the paths as dsem does", {

  skip_if_not_installed("dsem")

  variables <- c("a", "b")
  n_t <- 15
  set.seed(14)
  x_grid <- matrix(rnorm(n_t * 2), n_t, 2, dimnames = list(NULL, variables))
  arrows <- c("a -> b, 0, p", "a -> a, 1, r", "b -> b, 1, rb", "a <-> a, 0, sa", "b <-> b, 0, sb")
  values <- c(p = 0.5, r = 0.6, rb = 0.3, sa = 0.9, sb = 0.7)
  delta0 <- c(0.8, -0.4)

  dsem_model <- read_dsem_arrows(arrows, variables)
  dsem_cells <- get_dsem_cells(dsem_model, n_t)
  beta <- unname(values[dsem_model$beta_names])
  ln_sd <- log(unname(values[dsem_model$ln_sd_names]))
  ours <- get_dsem_nLL(beta, ln_sd, x_grid, 0 * x_grid, dsem_model, dsem_cells, delta0 = delta0)

  y_na <- stats::ts(x_grid)
  y_na[] <- NA
  fit <- dsem::dsem(sem = arrows, tsdata = y_na, family = list(a = dsem::fixed(), b = dsem::fixed()),
                    estimate_delta0 = TRUE,
                    control = dsem::dsem_control(run_model = FALSE, quiet = TRUE, gmrf_parameterization = "full"))
  cpp_par <- fit$obj$env$last.par
  cpp_par[names(cpp_par) == "beta_z"] <- values
  cpp_par[names(cpp_par) == "x_tj"] <- as.vector(x_grid)
  cpp_par[names(cpp_par) == "delta0_j"] <- delta0
  expect_equal(ours, fit$obj$env$f(cpp_par, order = 0), tolerance = 1e-10)

  # offsets of zero leave the density exactly where it was, and the offset decays through the ar1
  expect_equal(get_dsem_nLL(beta, ln_sd, x_grid, 0 * x_grid, dsem_model, dsem_cells, delta0 = c(0, 0)),
               get_dsem_nLL(beta, ln_sd, x_grid, 0 * x_grid, dsem_model, dsem_cells), tolerance = 1e-12)

  parts <- get_dsem_matrices(beta, ln_sd, dsem_model, dsem_cells, need_Vinv = FALSE)
  delta_cell <- rep(0, dsem_cells$n_cells)
  delta_cell[1] <- 1
  propagated <- as.vector(solve(parts$IminusB, delta_cell))[1:n_t]
  expect_equal(propagated, 0.6^(0:(n_t - 1)), tolerance = 1e-10) # an impulse on a, decaying at its own rho

})

test_that("the variance given the known cells matches a dense Schur complement and the AR1 sum", {

  case <- dsem_test_cases$two_series # an AR1 covariate with a same-year and a lagged path into rec
  dsem_model <- read_dsem_arrows(case$arrows, case$variables)
  dsem_cells <- get_dsem_cells(dsem_model, 12)
  x_grid <- matrix(0, 12, 2)
  case_pars <- dsem_case_pars(dsem_model, case$values)

  # nothing known: the covariate's variance is the AR1 sum from an unconditioned first year
  none <- matrix(FALSE, 12, 2)
  v_none <- get_dsem_margvar(case_pars$dsem_beta, case_pars$ln_dsem_sd, x_grid, dsem_model, dsem_cells, as.vector(none))
  expect_equal(as.numeric(v_none[,1]), 0.9^2 * cumsum(0.55^(2 * (0:11))), tolerance = 1e-10)
  expect_equal(as.numeric(v_none), as.numeric(dense_dsem_margvar(case_pars$dsem_beta, case_pars$ln_dsem_sd, x_grid, dsem_model, none)), tolerance = 1e-10)

  # the covariate known in the first eight years: recruitment's variance there is its own sd squared, and
  # once the covariate is unknown again both paths feed it
  known <- none
  known[1:8, 1] <- TRUE
  v_known <- get_dsem_margvar(case_pars$dsem_beta, case_pars$ln_dsem_sd, x_grid, dsem_model, dsem_cells, as.vector(known))
  expect_equal(as.numeric(v_known[1:8, 1]), rep(0, 8))
  expect_equal(as.numeric(v_known[1:8, 2]), rep(0.65^2, 8), tolerance = 1e-10)
  expect_true(all(v_known[9:12, 2] > 0.65^2))
  expect_equal(as.numeric(v_known), as.numeric(dense_dsem_margvar(case_pars$dsem_beta, case_pars$ln_dsem_sd, x_grid, dsem_model, known)), tolerance = 1e-10)

  # and on the tape it differentiates
  f <- function(p) sum(get_dsem_margvar(p$dsem_beta, p$ln_dsem_sd, x_grid, dsem_model, dsem_cells, as.vector(known)))
  obj <- RTMB::MakeADFun(f, list(dsem_beta = case_pars$dsem_beta, ln_dsem_sd = case_pars$ln_dsem_sd), silent = TRUE)
  expect_equal(obj$fn(obj$par), sum(v_known), tolerance = 1e-10)
  expect_equal(as.vector(obj$gr(obj$par)), numDeriv::grad(obj$fn, obj$par), tolerance = 1e-6)

})
