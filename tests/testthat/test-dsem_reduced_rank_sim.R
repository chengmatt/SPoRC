# Checks the operating model draw for a reduced rank grid: a projected cell is set from the drawn
# cells rather than drawn (1e-12), those cells carry Q_oo's covariance (5 SE), conditioned years are
# kept, and a full rank grid draws as the whole grid precision did (1e-14).

# a minimal operating model environment, the way test-dsem_simulation.R builds its Monte Carlo one
rr_sim_env <- function(dsem_model,
                       dsem_beta,
                       ln_dsem_sd,
                       n_yrs,
                       n_sims,
                       n_cond = 0,
                       x_fit = NULL) {

  n_vars <- length(dsem_model$variables)
  e <- new.env()

  e$dsem_model <- dsem_model
  e$dsem_beta <- dsem_beta
  e$ln_dsem_sd <- ln_dsem_sd
  e$n_yrs <- n_yrs
  e$n_sims <- n_sims

  # no covariate series and no linked process, so every mean is zero and nothing is written to an array
  e$dsem_cov_var_idx <- integer(0)
  e$dsem_mu <- rep(0, n_vars)
  e$dsem_link_par <- character(0)
  e$dsem_link_sim_par <- character(0)
  e$dsem_link_col <- integer(0)
  e$rec_bias_correct <- 0
  e$dsem_rec_corr_on <- FALSE

  e$dsem_n_cond_yrs <- n_cond
  e$dsem_x_fit <- if(is.null(x_fit)) matrix(0, n_yrs, n_vars) else x_fit
  e$dsem_x_known <- matrix(0, n_yrs, n_vars)

  return(e)

} # end function

# one factor on a unit random walk, two series loaded off it with no innovation of their own
rr_arrows <- "
F1 <-> F1, 0, NA, 1
F1 -> F1, 1, NA, 1
F1 -> env_1, 0, load_1
F1 -> env_2, 0, load_2
env_1 <-> env_1, 0, NA, 0
env_2 <-> env_2, 0, NA, 0
"
rr_vars <- c("F1", "env_1", "env_2")
rr_loads <- c(load_1 = 0.8, load_2 = -0.5)
rr_model <- read_dsem_arrows(rr_arrows, rr_vars)
rr_beta <- unname(rr_loads[rr_model$beta_names])

# the precision the density reads, built the way the objective builds it
rr_Q_oo <- function(n_yrs) {
  cells <- get_dsem_cells(rr_model, n_yrs)
  parts <- get_dsem_matrices(rr_beta, numeric(0), rr_model, cells, matrix(0, n_yrs, length(rr_vars)),
                             need_Vinv = FALSE, need_V = FALSE)
  get_dsem_Q_oo(parts$IminusB, parts, get_dsem_solve_mat(parts$IminusB, cells), cells)
} # end function


test_that("both loaded series are projected off the factor rather than drawn", {

  n_yrs <- 8
  n_sims <- 200

  expect_equal(rr_vars[rr_model$project_k], c("env_1", "env_2")) # the two with an sd of zero

  om <- rr_sim_env(rr_model, rr_beta, numeric(0), n_yrs, n_sims)
  set.seed(101)
  draw_dsem_sim(om)
  x <- om$dsem_x_sim

  expect_equal(dim(x), c(n_yrs, length(rr_vars), n_sims))

  # the whole grid precision put an infinity on a zero sd, so every cell used to come back NaN
  expect_true(all(is.finite(x)))

  # each loaded series is its loading times the factor, since its mean is zero and only F1 points in
  expect_equal(x[,2,], rr_loads[["load_1"]] * x[,1,], tolerance = 1e-12)
  expect_equal(x[,3,], rr_loads[["load_2"]] * x[,1,], tolerance = 1e-12)

  # and the factor itself is not flat, so the identity above is not holding at zero
  expect_gt(stats::sd(as.vector(x[,1,])), 0.5)
})


test_that("the cells that keep an innovation are drawn with Q_oo's covariance", {

  n_yrs <- 6
  n_sims <- 20000

  om <- rr_sim_env(rr_model, rr_beta, numeric(0), n_yrs, n_sims)
  set.seed(102)
  draw_dsem_sim(om)

  Sigma <- solve(as.matrix(rr_Q_oo(n_yrs)))
  drawn <- t(om$dsem_x_sim[,1,]) # the factor's years, one replicate per row

  # mean and covariance within 5 SE of the field the density describes
  se_mean <- sqrt(diag(Sigma) / n_sims)
  expect_lt(max(abs(colMeans(drawn) / se_mean)), 5)

  se_cov <- sqrt((outer(diag(Sigma), diag(Sigma)) + Sigma^2) / n_sims)
  expect_lt(max(abs((stats::cov(drawn) - Sigma) / se_cov)), 5)

  # a random walk with unit innovations, so the last year is the most variable
  expect_equal(unname(diag(Sigma)), as.numeric(1:n_yrs), tolerance = 1e-8)
})


test_that("conditioned years are kept and the projected cells still read off the factor", {

  n_yrs <- 8
  n_cond <- 3
  n_sims <- 50

  set.seed(100)
  x_fit <- matrix(0, n_yrs, length(rr_vars))
  x_fit[,1] <- cumsum(stats::rnorm(n_yrs))
  x_fit[,2] <- rr_loads[["load_1"]] * x_fit[,1]
  x_fit[,3] <- rr_loads[["load_2"]] * x_fit[,1]

  om <- rr_sim_env(rr_model, rr_beta, numeric(0), n_yrs, n_sims, n_cond = n_cond, x_fit = x_fit)
  set.seed(103)
  draw_dsem_sim(om)
  x <- om$dsem_x_sim

  # the factor's fitted years come back at the values given, and the later years do not
  for(sim in 1:n_sims) expect_equal(x[1:n_cond,1,sim], x_fit[1:n_cond,1], tolerance = 1e-12)
  expect_gt(stats::sd(x[n_yrs,1,]), 0.5)

  # every year of a loaded series reads off the factor, the conditioned ones included
  expect_equal(x[,2,], rr_loads[["load_1"]] * x[,1,], tolerance = 1e-12)
  expect_equal(x[,3,], rr_loads[["load_2"]] * x[,1,], tolerance = 1e-12)
})


test_that("a full rank grid draws as the whole grid precision did", {

  fr_arrows <- c("env <-> env, 0, sd_env", "env -> env, 1, rho",
                 "rec <-> rec, 0, sd_rec", "env -> rec, 0, b_env")
  fr_vars <- c("env", "rec")
  fr_model <- read_dsem_arrows(paste(fr_arrows, collapse = "\n"), fr_vars)
  fr_beta <- unname(c(rho = 0.6, b_env = 0.4)[fr_model$beta_names])
  fr_ln_sd <- log(unname(c(sd_env = 0.7, sd_rec = 0.5)[fr_model$ln_sd_names]))

  n_yrs <- 10
  n_sims <- 50
  n_cells <- n_yrs * length(fr_vars)

  expect_false(any(fr_model$project_k)) # nothing projected, so Q_oo is the whole grid precision

  om <- rr_sim_env(fr_model, fr_beta, fr_ln_sd, n_yrs, n_sims)
  set.seed(104)
  draw_dsem_sim(om)

  # the same draw taken from the whole grid precision, which is the route the simulator used before
  set.seed(104)
  Q <- get_dsem_precision(fr_beta, fr_ln_sd, fr_model, get_dsem_cells(fr_model, n_yrs))
  cond <- get_dsem_conditional(Q, rep(0, n_cells), integer(0), numeric(0))
  flat <- array(0, dim = c(n_cells, n_sims))
  flat[cond$unknown_cell,] <- draw_dsem_conditional(cond, n_sims)

  expect_equal(om$dsem_x_sim, array(flat, dim = c(n_yrs, length(fr_vars), n_sims)), tolerance = 1e-14)
})
