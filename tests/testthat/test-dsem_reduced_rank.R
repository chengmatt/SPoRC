# Checks a reduced rank dsem, where a series with an sd of zero is projected off the others, against the dsem
# package's gmrf_project: the precision Q_oo entry by entry, the density, and the cells themselves.

# dsem's own reduced rank fit at a state we set. A latent series is a column of NA; every other series takes a
# family so its cells stay free, and the state goes in as those free cells
dsem_project <- function(sem, vars, latent, n_t, state_o, delta0 = NULL) {
  cols <- stats::setNames(lapply(vars, function(v) if(v %in% latent) rep(NA_real_, n_t) else stats::rnorm(n_t)), vars)
  fam <- stats::setNames(lapply(vars, function(v) if(v %in% latent) dsem::fixed() else stats::gaussian()), vars)
  d <- dsem::dsem(tsdata = stats::ts(as.data.frame(cols)), sem = sem, family = fam, estimate_delta0 = !is.null(delta0),
                  control = dsem::dsem_control(run_model = FALSE, quiet = TRUE, gmrf_parameterization = "gmrf_project"))
  p <- d$obj$env$last.par
  expect_equal(sum(names(p) == "x_tj"), length(state_o)) # the cells that keep an innovation, and no others
  p[names(p) == "x_tj"] <- state_o
  if(!is.null(delta0)) p[names(p) == "delta0_j"] <- delta0
  r <- d$obj$report(p)
  expect_true(!is.null(r$jnll_gmrf) && is.finite(r$jnll_gmrf))
  list(nll = as.numeric(r$jnll_gmrf), z = matrix(as.vector(r$z_tj), n_t), Q_oo = as.matrix(r$Q_oo))
}

# the same model through us, arrow values read out of the model text by the name each one was given
ours <- function(sem, vars, n_t, values, state_o, delta0 = NULL) {
  m <- read_dsem_arrows(sem, vars)
  cells <- get_dsem_cells(m, n_t)
  x <- matrix(0, n_t, length(vars)); x[, !m$project_k] <- state_o
  mu <- matrix(0, n_t, length(vars))
  beta <- if(length(m$beta_names)) unname(values[m$beta_names]) else numeric(0)
  ln_sd <- if(length(m$ln_sd_names)) log(unname(values[m$ln_sd_names])) else numeric(0)
  grid <- get_dsem_grid(beta, ln_sd, x, mu, m, cells, delta0 = delta0)
  list(nll = as.numeric(get_dsem_nLL(beta, ln_sd, x, mu, m, cells, delta0 = delta0)),
       Q_oo = as.matrix(get_dsem_Q_oo(grid$parts$IminusB, grid$parts, grid$solve_mat, cells)),
       z = grid$x_grid, projected = vars[m$project_k], model = m, cells = cells, grid = grid)
}

# and the precision worked out the long way, from dense blocks, as a check on both
by_hand <- function(o, n_t) {
  ImB <- as.matrix(o$grid$parts$IminusB)
  u <- o$cells$unobs_idx; ob <- o$cells$obs_idx
  V <- if(is.null(o$grid$parts$sd_cell)) as.matrix(o$grid$parts$V) else diag(o$grid$parts$sd_cell^2)
  S <- ImB[ob,ob,drop = FALSE] - ImB[ob,u,drop = FALSE] %*% solve(ImB[u,u,drop = FALSE], ImB[u,ob,drop = FALSE])
  t(S) %*% solve(V[ob,ob,drop = FALSE]) %*% S
}

set.seed(11)
n_t <- 8
vals <- c(sd_x = 0.9, sd_e = 0.6, b_xy = 0.7, b_yx = 0.4, b_fx = 0.8, b_fy = 0.5, c_xe = 0.2)
state <- stats::rnorm(n_t, 0, 0.8)

test_that("a one factor dynamic factor analysis matches dsem's gmrf_project", {
  skip_if_not_installed("dsem")

  # the factor holds all the process variance, so both manifest series are their loadings on it and nothing else
  sem <- c("F1 <-> F1, 0, NA, 1", "F1 -> F1, 1, NA, 1", "F1 -> x, 0, b_fx", "F1 -> y, 0, b_fy",
           "x <-> x, 0, NA, 0", "y <-> y, 0, NA, 0")
  o <- ours(sem, c("F1", "x", "y"), n_t, vals, state)
  d <- dsem_project("F1 <-> F1, 0, NA, 1
                     F1 -> F1, 1, NA, 1
                     F1 -> x, 0, b_fx, 0.8
                     F1 -> y, 0, b_fy, 0.5
                     x <-> x, 0, NA, 0
                     y <-> y, 0, NA, 0", c("F1", "x", "y"), "F1", n_t, state)

  expect_equal(o$projected, c("x", "y"))
  expect_equal(o$cells$unobs_idx, n_t + 1:(2 * n_t)) # the two manifest series, the factor left in the density
  expect_equal(o$Q_oo, d$Q_oo, tolerance = 1e-10)
  expect_equal(o$nll, d$nll, tolerance = 1e-10)
  expect_equal(unname(o$z), unname(d$z), tolerance = 1e-10)
  expect_equal(o$z[,2], 0.8 * state, tolerance = 1e-12) # the loading times the factor, exactly
  expect_equal(o$Q_oo, by_hand(o, n_t), tolerance = 1e-10)

})

test_that("a same-year loop through a projected series matches dsem, determinant and all", {
  skip_if_not_installed("dsem")

  sem <- c("x <-> x, 0, sd_x", "y <-> y, 0, NA, 0", "x -> y, 0, b_xy", "y -> x, 0, b_yx")
  o <- ours(sem, c("x", "y"), n_t, vals, state)
  d <- dsem_project("x <-> x, 0, sd_x, 0.9
                     y <-> y, 0, NA, 0
                     x -> y, 0, b_xy, 0.7
                     y -> x, 0, b_yx, 0.4", c("x", "y"), character(0), n_t, state)

  expect_equal(o$Q_oo, d$Q_oo, tolerance = 1e-10)
  expect_equal(o$nll, d$nll, tolerance = 1e-10)

  # the loop puts det(S) at (1 - b b')^T, which the precision accounts for and the innovations alone do not
  innovation <- as.vector(as.matrix(o$grid$parts$IminusB) %*% as.vector(o$z))[o$cells$obs_idx]
  expect_equal(o$nll, -sum(stats::dnorm(innovation, 0, 0.9, log = TRUE)) - n_t * log(1 - 0.7 * 0.4), tolerance = 1e-10)
  expect_equal(det(by_hand(o, n_t))^0.5, (1 - 0.7 * 0.4)^n_t / 0.9^n_t, tolerance = 1e-8)

})

test_that("a lead into a projected series matches dsem", {
  skip_if_not_installed("dsem")
  sem <- c("x <-> x, 0, sd_x", "y <-> y, 0, NA, 0", "x -> y, -1, b_xy")
  o <- ours(sem, c("x", "y"), n_t, vals, state)
  d <- dsem_project("x <-> x, 0, sd_x, 0.9
                     y <-> y, 0, NA, 0
                     x -> y, -1, b_xy, 0.7", c("x", "y"), character(0), n_t, state)
  expect_equal(o$Q_oo, d$Q_oo, tolerance = 1e-10)
  expect_equal(o$nll, d$nll, tolerance = 1e-10)
})

test_that("a covariance arrow alongside a projected series matches dsem", {
  skip_if_not_installed("dsem")
  state2 <- cbind(stats::rnorm(n_t, 0, 0.8), stats::rnorm(n_t, 0, 0.5))
  sem <- c("x <-> x, 0, sd_x", "e <-> e, 0, sd_e", "x <-> e, 0, c_xe", "y <-> y, 0, NA, 0", "x -> y, 0, b_fx")
  o <- ours(sem, c("x", "e", "y"), n_t, vals, state2)
  d <- dsem_project("x <-> x, 0, sd_x, 0.9
                     e <-> e, 0, sd_e, 0.6
                     x <-> e, 0, c_xe, 0.2
                     y <-> y, 0, NA, 0
                     x -> y, 0, b_fx, 0.8", c("x", "e", "y"), character(0), n_t, as.vector(state2))
  expect_equal(o$Q_oo, d$Q_oo, tolerance = 1e-10)
  expect_equal(o$nll, d$nll, tolerance = 1e-10)

  # a covariance line onto the projected series itself has no innovation to covary with, so it is refused
  expect_error(read_dsem_arrows(c(sem, "y <-> e, 0, c_ye"), c("x", "e", "y")), "no innovation to covary")
})

test_that("a first year offset works alongside a projected series", {
  skip_if_not_installed("dsem")
  sem <- c("F1 <-> F1, 0, NA, 1", "F1 -> F1, 1, NA, 0.7", "F1 -> x, 0, b_fx", "x <-> x, 0, NA, 0")
  o <- ours(sem, c("F1", "x"), n_t, vals, state, delta0 = c(0.5, 0))
  d <- dsem_project("F1 <-> F1, 0, NA, 1
                     F1 -> F1, 1, NA, 0.7
                     F1 -> x, 0, b_fx, 0.8
                     x <-> x, 0, NA, 0", c("F1", "x"), "F1", n_t, state, delta0 = c(0.5, 0))
  expect_equal(o$nll, d$nll, tolerance = 1e-10)
  expect_equal(unname(o$z), unname(d$z), tolerance = 1e-10) # the offset reaches the projected cells too
})

test_that("a chain of projected series through a moderated arrow keeps its pattern", {

  # the quadratic approximation to exp: three series with an sd of zero, one reading the next, and the arrow
  # between them moderated by a series that does keep an innovation
  variables <- c("ones", "logX", "logX1", "logX2", "X")
  sem <- c("ones -> ones, 1, NA, 1", "ones -> logX, 0, alpha", "logX -> logX, 1, NA, 1",
           "logX -> logX1, 0, NA, 1", "logX1 -> logX2, 0, logX", "ones -> X, 0, NA, 1",
           "logX1 -> X, 0, NA, 1", "logX2 -> X, 0, NA, 0.5", "ones <-> ones, 0, NA, 0.001",
           "logX <-> logX, 0, sd_logX", "logX1 <-> logX1, 0, NA, 0", "logX2 <-> logX2, 0, NA, 0",
           "X <-> X, 0, NA, 0")
  m <- read_dsem_arrows(sem, variables)
  expect_equal(variables[m$project_k], c("logX1", "logX2", "X"))

  n_grid <- 10
  cells <- get_dsem_cells(m, n_grid)
  x <- matrix(0, n_grid, 5); x[,1] <- 1; x[,2] <- seq(-0.6, 0.6, length.out = n_grid)
  mu <- matrix(0, n_grid, 5)
  grid <- get_dsem_grid(0.05, log(0.3), x, mu, m, cells)

  # the chain is that approximation, so it agrees with exp to the term left out
  expect_equal(grid$x_grid[,5], 1 + x[,2] + 0.5 * x[,2]^2, tolerance = 1e-12)
  expect_lt(max(abs(grid$x_grid[,5] - exp(x[,2]))), 0.05)

  # the projected cells owe nothing to the density, which is the two series that keep an innovation
  ours_nll <- as.numeric(get_dsem_nLL(0.05, log(0.3), x, mu, m, cells))
  by_hand_nll <- -stats::dnorm(x[1,1], 0, 0.001, TRUE) - sum(stats::dnorm(diff(x[,1]), 0, 0.001, TRUE)) -
    stats::dnorm(x[1,2] - 0.05 * x[1,1], 0, 0.3, TRUE) - sum(stats::dnorm(diff(x[,2]) - 0.05 * x[-1,2] * 0 - 0.05 * x[-1,1], 0, 0.3, TRUE))
  expect_equal(ours_nll, by_hand_nll, tolerance = 1e-8)

  # the pattern has to cover every entry the dense product puts in, or dgmrf would see a different matrix
  Q <- as.matrix(get_dsem_Q_oo(grid$parts$IminusB, grid$parts, grid$solve_mat, cells))
  expect_equal(Q, by_hand(list(grid = grid, cells = cells), n_grid), tolerance = 1e-10)

  expect_error(read_dsem_arrows(sub("logX -> logX1, 0, NA, 1", "", sem, fixed = TRUE), variables), "nothing sets their value")

})

test_that("the projection differentiates", {
  sem <- c("x <-> x, 0, sd_x", "y <-> y, 0, NA, 0", "x -> y, 0, b_xy", "y -> x, 0, b_yx")
  m <- read_dsem_arrows(sem, c("x", "y")); cells <- get_dsem_cells(m, n_t)
  f <- function(p) {
    xg <- matrix(0, n_t, 2); xg[,1] <- p[3:(n_t + 2)]
    get_dsem_nLL(p[1:2], log(p[n_t + 3]), xg, matrix(0, n_t, 2), m, cells)
  }
  p0 <- c(0.7, 0.4, state, 0.9)
  tape <- RTMB::MakeTape(f, p0)
  fd <- sapply(seq_along(p0), function(i) { e <- rep(0, length(p0)); e[i] <- 1e-6; (f(p0 + e) - f(p0 - e)) / 2e-6 })
  expect_equal(as.vector(tape$jacobian(p0)), fd, tolerance = 1e-6)
  expect_true(all(is.finite(tape$jacfun()$jacobian(p0)))) # the precision's log determinant has a second derivative too
})

test_that("a projected series cannot set the coefficient on an arrow", {

  # m1 has an sd of zero and moderates x's own lagged path, so I - B would need m1's cells before it can be
  # built and m1's cells come out of I - B. dsem refuses the same pairing
  sem <- c("z <-> z, 0, sd_z", "x <-> x, 0, sd_x", "m1 <-> m1, 0, NA, 0", "z -> m1, 0, b_zm", "x -> x, 1, m1")
  expect_error(read_dsem_arrows(sem, c("x", "z", "m1")), "also set the coefficient on an arrow")

  # the same arrows with an sd on m1 are fine, since its cells are then estimated
  ok <- read_dsem_arrows(sub("m1 <-> m1, 0, NA, 0", "m1 <-> m1, 0, sd_m", sem, fixed = TRUE), c("x", "z", "m1"))
  expect_false(any(ok$project_k))

})

test_that("a covariate dynamic factor analysis reproduces dsem's fit", {
  skip_if_not_installed("dsem")

  # two series with no process error of their own, both loading on one latent factor and both observed with
  # error. The first loading is fixed at one, which sets the factor's scale and sign so the fits are comparable
  set.seed(21)
  n_t <- 30
  factor_true <- cumsum(stats::rnorm(n_t, 0, 1))
  obs_mat <- cbind(env1 = 1.0 * factor_true + stats::rnorm(n_t, 0, 0.4) + 2,
                   env2 = 0.6 * factor_true + stats::rnorm(n_t, 0, 0.7) - 1)
  vars <- c("F1", "env1", "env2")
  arrows <- c("F1 <-> F1, 0, NA, 1", "F1 -> F1, 1, NA, 1", "F1 -> env1, 0, NA, 1", "F1 -> env2, 0, load2",
              "env1 <-> env1, 0, NA, 0", "env2 <-> env2, 0, NA, 0")

  theirs <- suppressWarnings(dsem::dsem(
    sem = "F1 <-> F1, 0, NA, 1
           F1 -> F1, 1, NA, 1
           F1 -> env1, 0, NA, 1
           F1 -> env2, 0, load2, 0.5
           env1 <-> env1, 0, NA, 0
           env2 <-> env2, 0, NA, 0",
    tsdata = stats::ts(data.frame(F1 = rep(NA_real_, n_t), env1 = obs_mat[,1], env2 = obs_mat[,2])),
    family = list(F1 = dsem::fixed(), env1 = stats::gaussian(), env2 = stats::gaussian()),
    estimate_mu = c("env1", "env2"),
    control = dsem::dsem_control(quiet = TRUE, getsd = FALSE, use_REML = FALSE,
                                 gmrf_parameterization = "gmrf_project")))

  # ours: the process density through the projection, plus the covariate observations read off the projected
  # cells, which is the wiring the objective uses. The factor is integrated out on both sides
  m <- read_dsem_arrows(arrows, vars)
  cells <- get_dsem_cells(m, n_t)

  ours_nll <- function(p) {
    RTMB::getAll(p)
    "c" <- RTMB::ADoverload("c") # testthat's env has the namespace as parent, where c is base and drops the AD class
    x_grid <- matrix(0, n_t, 3)
    x_grid[cells$obs_idx] <- factor_cells
    mu_grid <- matrix(rep(c(0, mu), each = n_t), n_t, 3)
    grid <- get_dsem_grid(load2, numeric(0), x_grid, mu_grid, m, cells)
    nll <- get_dsem_nLL(load2, numeric(0), x_grid, mu_grid, m, cells)
    for(k in 1:2) {
      nll <- nll + get_dsem_obs_nLL(y = obs_mat[,k], x = grid$x_grid[,k + 1], family = 1,
                                    link = dsem_default_link(1), obs_sd = exp(ln_obs_sd[k]), tweedie_p = 1.5)
    }
    nll
  }

  obj <- RTMB::MakeADFun(ours_nll, list(load2 = 0.5, ln_obs_sd = c(0, 0), mu = c(0, 0),
                                        factor_cells = rep(0, n_t)), random = "factor_cells", silent = TRUE)
  fit <- stats::nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 1000, eval.max = 2000))

  # dsem holds the loading, the two log observation sds and the two means in that order, on the same scales
  expect_equal(length(obj$env$random), length(theirs$obj$env$random))
  expect_equal(as.numeric(fit$objective), as.numeric(theirs$opt$objective), tolerance = 1e-8) # dsem's carries a logarithm attribute
  expect_equal(unname(fit$par), unname(theirs$opt$par), tolerance = 1e-5)

  # the factor is what a dynamic factor analysis is fitted for, so it has to agree too
  our_factor <- obj$env$last.par.best[names(obj$env$last.par.best) == "factor_cells"]
  their_factor <- as.matrix(theirs$obj$report(theirs$obj$env$last.par.best)$z_tj)[,1]
  expect_equal(unname(our_factor), unname(their_factor), tolerance = 1e-5)
  expect_gt(stats::cor(our_factor, factor_true), 0.9) # and it recovers the series it was simulated from

})
