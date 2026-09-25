# Checks the diagonal and marginal sd forms against the dsem package's C++ at the same values (1e-10),
# and that conditional, diagonal and marginal are one model under iid arrows, in the density and in a fit.
# The marginal form parts ways with dsem's only when a covariance line couples the innovations.

# dsem's joint density and precision at one set of arrow values, every series observed without error.
# dsem holds one entry per named arrow in declaration order, and the means after them
dsem_variance_joint <- function(sem, tsdata, values, mu, constant_variance) {

  vars <- colnames(tsdata)
  fam <- stats::setNames(lapply(vars, function(v) dsem::fixed()), vars)
  ctl <- dsem::dsem_control(quiet = TRUE, nlminb_loops = 0, newton_loops = 0, getsd = FALSE, gmrf_parameterization = "full",
                            use_REML = FALSE, constant_variance = constant_variance)
  d <- suppressWarnings(dsem::dsem(sem = sem, tsdata = tsdata, family = fam, estimate_mu = vars, covs = vars, control = ctl))

  first <- !duplicated(d$sem_full$parameter) & d$sem_full$parameter > 0 # one number per estimated arrow
  nLL <- as.numeric(d$obj$fn(c(unname(values[d$sem_full$name[first]]), mu)))

  return(list(nLL = nLL, Q = d$obj$report()$Q_kk))

}

# the same density and precision from our reader, at the same arrow values
our_variance_joint <- function(sem, tsdata, values, mu, variance) {

  vars <- colnames(tsdata)
  n_t <- nrow(tsdata)
  m <- read_dsem_arrows(sem, vars, variance = variance)
  cells <- get_dsem_cells(m, n_t)
  pars <- dsem_case_pars(m, values)

  nLL <- as.numeric(get_dsem_nLL(pars$dsem_beta,
                                 pars$ln_dsem_sd,
                                 x_grid = unclass(tsdata)[,],
                                 mu_grid = matrix(mu, n_t, length(vars), byrow = TRUE),
                                 dsem_model = m,
                                 dsem_cells = cells))

  return(list(nLL = nLL, Q = get_dsem_precision(pars$dsem_beta, pars$ln_dsem_sd, m, cells)))

}

cell_variance <- function(Q) diag(solve(as.matrix(Q))) # each cell's variance over the whole grid

set.seed(51)
n_t <- 20
x <- as.numeric(stats::arima.sim(list(ar = 0.5), n_t))
y <- 0.4 * x + rnorm(n_t, 0, 0.7)
xy <- stats::ts(data.frame(x = x, y = y), start = 1)
mu <- c(0.1, -0.2)

sem_plain <- c("x -> y, 0, b", "x -> x, 1, rho", "x <-> x, 0, sx", "y <-> y, 0, sy") # paths, no covariance line
sem_cov <- c(sem_plain, "x <-> y, 0, cxy") # the same with the innovations correlated
sem_iid <- c("x <-> x, 0, sx", "y <-> y, 0, sy") # no path at all
values <- c(b = 0.4, rho = 0.5, sx = 0.9, sy = 0.7, cxy = 0.3)

test_that("the diagonal form matches dsem's C++, covariance line or not", {

  skip_if_not_installed("dsem")

  # both solve each cell's innovation variance from the sd lines as if the innovations were independent,
  # so the density and every cell's variance agree to the bit
  for(sem in list(sem_plain, sem_cov)) {
    ours <- our_variance_joint(sem, xy, values, mu, "diagonal")
    theirs <- dsem_variance_joint(sem, xy, values, mu, "diagonal")
    expect_equal(ours$nLL, theirs$nLL, tolerance = 1e-10)
    expect_equal(cell_variance(ours$Q), cell_variance(theirs$Q), tolerance = 1e-8, ignore_attr = TRUE)
  } # end sem loop

})

test_that("without a covariance line the marginal form matches dsem's C++ and lands on the sd lines", {

  skip_if_not_installed("dsem")
  ours <- our_variance_joint(sem_plain, xy, values, mu, "marginal")
  theirs <- dsem_variance_joint(sem_plain, xy, values, mu, "marginal")
  expect_equal(ours$nLL, theirs$nLL, tolerance = 1e-10)

  # with independent innovations there is nothing between the two forms, on either side
  expect_equal(ours$nLL, our_variance_joint(sem_plain, xy, values, mu, "diagonal")$nLL, tolerance = 1e-12)
  expect_equal(theirs$nLL, dsem_variance_joint(sem_plain, xy, values, mu, "diagonal")$nLL, tolerance = 1e-10)

  # and both leave every cell at its sd line squared
  target <- rep(c(values[["sx"]], values[["sy"]])^2, each = n_t)
  expect_equal(cell_variance(ours$Q), target, tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(cell_variance(theirs$Q), target, tolerance = 1e-8, ignore_attr = TRUE)

  # the conditional form is the one that differs: a self path of 0.5 widens x to sx^2 / (1 - rho^2)
  cond <- our_variance_joint(sem_plain, xy, values, mu, "conditional")
  expect_equal(cond$nLL, dsem_variance_joint(sem_plain, xy, values, mu, "conditional")$nLL, tolerance = 1e-10)
  expect_equal(cell_variance(cond$Q)[n_t], values[["sx"]]^2 / (1 - values[["rho"]]^2), tolerance = 1e-4, ignore_attr = TRUE)

})

test_that("with a covariance line the marginal form settles on the sd lines where dsem's rescale does not", {

  skip_if_not_installed("dsem")
  ours <- our_variance_joint(sem_cov, xy, values, mu, "marginal")
  theirs <- dsem_variance_joint(sem_cov, xy, values, mu, "marginal")
  target <- rep(c(values[["sx"]], values[["sy"]])^2, each = n_t)

  # ours scales the innovations until every cell's variance, cross terms included, is its sd line squared
  expect_equal(cell_variance(ours$Q), target, tolerance = 1e-8, ignore_attr = TRUE)

  # dsem scales the rows of I - B and of Gamma once rather than solving, so the cross terms are left
  # over and x settles about 20% above its sd line
  expect_gt(max(abs(cell_variance(theirs$Q) - target) / target), 0.05)
  expect_false(isTRUE(all.equal(ours$nLL, theirs$nLL)))

  # the conditional and diagonal forms still agree with dsem under the same covariance line
  expect_equal(our_variance_joint(sem_cov, xy, values, mu, "conditional")$nLL,
               dsem_variance_joint(sem_cov, xy, values, mu, "conditional")$nLL, tolerance = 1e-10)
  expect_equal(our_variance_joint(sem_cov, xy, values, mu, "diagonal")$nLL,
               dsem_variance_joint(sem_cov, xy, values, mu, "diagonal")$nLL, tolerance = 1e-10)

})

test_that("under iid arrows the three forms are one model, ours and dsem's", {

  skip_if_not_installed("dsem")

  # with no path, I - B is the identity: the innovation sd and the marginal sd are the same number,
  # so there is nothing for the diagonal or marginal solve to change
  forms <- c("conditional", "diagonal", "marginal")
  ours <- theirs <- stats::setNames(rep(0, 3), forms)

  for(v in forms) {
    ours[v] <- our_variance_joint(sem_iid, xy, values, mu, v)$nLL
    theirs[v] <- dsem_variance_joint(sem_iid, xy, values, mu, v)$nLL
  } # end v loop

  expect_equal(ours[["diagonal"]], ours[["conditional"]], tolerance = 1e-12)
  expect_equal(ours[["marginal"]], ours[["conditional"]], tolerance = 1e-12)
  expect_equal(theirs[["diagonal"]], theirs[["conditional"]], tolerance = 1e-12)
  expect_equal(theirs[["marginal"]], theirs[["conditional"]], tolerance = 1e-12)
  expect_equal(unname(ours), unname(theirs), tolerance = 1e-10)

  # and that one number is the plain sum of normal densities, one per cell at its own sd line
  by_hand <- -sum(stats::dnorm(x, mu[1], values[["sx"]], TRUE)) - sum(stats::dnorm(y, mu[2], values[["sy"]], TRUE))
  expect_equal(ours[["conditional"]], by_hand, tolerance = 1e-10)

})

test_that("in a fit the three forms read iid recruitment arrows the same way", {

  # the same check through Setup_Mod_DSEM and fit_model, where the sd line reaches the recruitment
  # correction as well as the density
  plain <- suppressMessages(sweep_input(rec = list(RecDevs_model = "dsem"), dims = list(n_regions = 1)))
  arrows <- "rec <-> rec, 0, sd_rec, 0.7"
  fits <- list()

  for(v in c("conditional", "diagonal", "marginal")) {
    il <- suppressMessages(Setup_Mod_DSEM(plain, arrows, NULL, dsem_variance = v))
    obj <- fit_model(il$data, il$par, il$map, random = NULL, do_optim = FALSE, silent = TRUE)
    fits[[v]] <- list(fn = as.numeric(obj$fn(obj$par)), gr = as.numeric(obj$gr(obj$par)), margvar = obj$rep$dsem_margvar_grid)
  } # end v loop

  expect_equal(fits$diagonal$fn, fits$conditional$fn, tolerance = 1e-12)
  expect_equal(fits$marginal$fn, fits$conditional$fn, tolerance = 1e-12)
  expect_equal(fits$diagonal$gr, fits$conditional$gr, tolerance = 1e-10)
  expect_equal(fits$marginal$gr, fits$conditional$gr, tolerance = 1e-10)

  # the correction reads 0.7^2 under every form, since no path adds to it
  for(v in names(fits)) expect_equal(as.numeric(fits[[v]]$margvar), rep(0.7^2, length(fits[[v]]$margvar)), tolerance = 1e-10)

})
