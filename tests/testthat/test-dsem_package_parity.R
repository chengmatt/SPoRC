# Checks our density against the dsem package's C++ at the same values (1e-10) and our fit against its fit (1e-4),
# every series observed as data so both sides evaluate the same joint density. Skipped without dsem.

dsem_joint <- function(sem, tsdata, dsem_par, covs = NULL, mod_var_logscale = FALSE) {
  vars <- colnames(tsdata)
  fam <- stats::setNames(lapply(vars, function(v) dsem::fixed()), vars)
  ctl <- dsem::dsem_control(quiet = TRUE, nlminb_loops = 0, newton_loops = 0, getsd = FALSE, gmrf_parameterization = "full",
                            use_REML = FALSE, logscale_moderating_variance = mod_var_logscale)
  d <- suppressWarnings(dsem::dsem(sem = sem, tsdata = tsdata, family = fam, estimate_mu = vars, covs = if(is.null(covs)) vars else covs, control = ctl))
  d$obj$fn(dsem_par)
}

our_joint <- function(sem, tsdata, beta, ln_sd, mu, covs = NULL, mod_var_logscale = FALSE) {
  vars <- colnames(tsdata); n_t <- nrow(tsdata)
  m <- read_dsem_arrows(sem, vars, covs = covs, mod_var_logscale = mod_var_logscale)
  as.numeric(get_dsem_nLL(beta, ln_sd, x_grid = unclass(tsdata)[,], mu_grid = matrix(mu, n_t, length(vars), byrow = TRUE),
                          dsem_model = m, dsem_cells = get_dsem_cells(m, n_t)))
}

set.seed(1)
n_t <- 20
x <- as.numeric(stats::arima.sim(list(ar = 0.5), n_t)); y <- 0.4 * x + rnorm(n_t, 0, 0.7); z <- rnorm(n_t)
xy <- stats::ts(data.frame(x = x, y = y), start = 1)
xyz <- stats::ts(data.frame(x = x, y = y, z = z), start = 1)

test_that("lagged and same-year paths match dsem's C++", {
  skip_if_not_installed("dsem")
  sem <- c("x -> y, 0, b", "x -> x, 1, rho", "x <-> x, 0, sx", "y <-> y, 0, sy")
  expect_equal(our_joint(sem, xy, c(0.4, 0.5), log(c(0.9, 0.7)), c(0.1, -0.2)), dsem_joint(sem, xy, c(0.4, 0.5, 0.9, 0.7, 0.1, -0.2)), tolerance = 1e-10)
})

test_that("a covariance arrow matches dsem's C++", {
  skip_if_not_installed("dsem")
  sem <- c("x -> x, 1, rho", "x <-> x, 0, sx", "y <-> y, 0, sy", "x <-> y, 0, cxy")
  expect_equal(our_joint(sem, xy, c(0.5, 0.3), log(c(0.9, 0.7)), c(0, 0)), dsem_joint(sem, xy, c(0.5, 0.9, 0.7, 0.3, 0, 0)), tolerance = 1e-10)
})

test_that("a lead matches dsem's C++", {
  skip_if_not_installed("dsem")
  sem <- c("x -> y, -1, b", "x -> x, 1, rho", "x <-> x, 0, sx", "y <-> y, 0, sy")
  expect_equal(our_joint(sem, xy, c(0.4, 0.5), log(c(0.9, 0.7)), c(0, 0)), dsem_joint(sem, xy, c(0.4, 0.5, 0.9, 0.7, 0, 0)), tolerance = 1e-10)
})

test_that("a moderated path and a moderated sd match dsem's C++", {
  skip_if_not_installed("dsem")
  sem_p <- c("x -> y, 0, z", "x -> x, 1, rho", "z -> z, 1, rz", "x <-> x, 0, sx", "y <-> y, 0, sy", "z <-> z, 0, sz")
  expect_equal(our_joint(sem_p, xyz, c(0.5, 0.3), log(c(0.9, 0.7, 0.5)), c(0, 0, 0)), dsem_joint(sem_p, xyz, c(0.5, 0.3, 0.9, 0.7, 0.5, 0, 0, 0)), tolerance = 1e-10)
  sem_v <- c("x -> y, 0, b", "x -> x, 1, rho", "z -> z, 1, rz", "x <-> x, 0, sx", "y <-> y, 0, z", "z <-> z, 0, sz")
  expect_equal(our_joint(sem_v, xyz, c(0.4, 0.5, 0.3), log(c(0.9, 0.5)), c(0, 0, 0), mod_var_logscale = TRUE),
               dsem_joint(sem_v, xyz, c(0.4, 0.5, 0.3, 0.9, 0.5, 0, 0, 0), mod_var_logscale = TRUE), tolerance = 1e-10)
})

test_that("fitting the same dsem gives the dsem package's estimates", {

  skip_if_not_installed("dsem")
  sem <- c("x -> y, 0, b", "x -> x, 1, rho", "x <-> x, 0, sx", "y <-> y, 0, sy")
  vars <- colnames(xy)

  # dsem's fit, every series observed without error and the means estimated
  fam <- stats::setNames(lapply(vars, function(v) dsem::fixed()), vars)
  ctl <- dsem::dsem_control(quiet = TRUE, getsd = FALSE, gmrf_parameterization = "full", use_REML = FALSE)
  theirs <- suppressWarnings(dsem::dsem(sem = sem, tsdata = xy, family = fam, estimate_mu = vars, control = ctl))

  # ours: the same density taped by RTMB and minimized the same way. dsem holds the sds themselves in
  # beta_z, so ours are exponentiated for the comparison
  m <- read_dsem_arrows(sem, vars)
  cells <- get_dsem_cells(m, n_t)
  x_grid <- unclass(xy)[,]
  ours_nll <- function(p) {
    RTMB::getAll(p)
    get_dsem_nLL(beta, ln_sd, x_grid = x_grid, mu_grid = matrix(mu, n_t, length(vars), byrow = TRUE), dsem_model = m, dsem_cells = cells)
  }
  obj <- RTMB::MakeADFun(ours_nll, list(beta = c(0, 0), ln_sd = c(0, 0), mu = c(0, 0)), silent = TRUE)
  ours <- stats::nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 1000, eval.max = 2000))
  ours_as_theirs <- c(ours$par[1:2], exp(ours$par[3:4]), ours$par[5:6])

  expect_equal(unname(ours_as_theirs), unname(theirs$opt$par), tolerance = 1e-4)
  expect_equal(ours$objective, theirs$opt$objective, tolerance = 1e-8)
  expect_equal(obj$fn(c(theirs$opt$par[1:2], log(theirs$opt$par[3:4]), theirs$opt$par[5:6])), theirs$opt$objective, tolerance = 1e-8) # ours at their optimum

})
