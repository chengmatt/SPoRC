# The numbers at age state penalty whitens its population, region, season and sex correlations by hand.
# Checks that against one separable density with multivariate normal factors, and the normal by hand (1e-8).

library(SPoRC)
library(testthat)

test_that("the state penalty is the Kronecker normal, written as dseparable with dmvnorm factors", {

  set.seed(8)
  np <- 2; nr <- 2; ny <- 4; nk <- 2; na <- 3; ns <- 2
  dims <- c(np, nr, ny, nk, na, ns) # population, region, year, season, age, sex
  ln_NAA <- array(rnorm(prod(dims), 0, 0.3), dim = dims)
  pred <- array(1, dim = dims) # so the innovation is the state itself
  sd <- 0.25 # conditional sd
  rho_a <- SPoRC:::rho_trans(0.4); rho_y <- SPoRC:::rho_trans(-0.5)
  pe <- array(0, dim = c(np, nr, 3, ns)); pe[,,1,] <- 0.4; pe[,,2,] <- -0.5 # the same ar1 on every population, region and sex
  pop_c <- 0.3; sex_c <- -0.6; reg_c <- 0.5; seas_c <- 0.2 # one unconstrained parameter each for two levels
  reg_pars <- array(reg_c, dim = c(np, 1, ns)); seas_pars <- array(seas_c, dim = c(np, 1, ns))

  nll <- SPoRC:::Get_NAA_state_penalty(ln_NAA, pred, array(sd, dim = dims), naa_re_ages = 1:na, naa_re_yrs = 1:ny, naa_re_seas = 1:nk,
                                       NAA_re = 4, NAA_pe_pars = pe, bias_correct = 0,
                                       NAA_re_region = 1, NAA_region_corr_pars = reg_pars, NAA_re_pop = 1, NAA_pop_corr_pars = pop_c,
                                       NAA_re_sex = 1, NAA_sex_corr_pars = sex_c, NAA_re_season = 1, NAA_season_corr_pars = seas_pars)

  # the normal by hand: sex, age, season, year, region and population from slowest to fastest, at the marginal sd
  C <- function(p, k) SPoRC:::build_us_corr(p, k)
  ar1 <- function(k, r) r^abs(outer(1:k, 1:k, "-"))
  marginal2 <- sd^2 / (1 - rho_a^2) / (1 - rho_y^2)
  S <- marginal2 * kronecker(C(sex_c, ns), kronecker(ar1(na, rho_a), kronecker(C(seas_c, nk), kronecker(ar1(ny, rho_y), kronecker(C(reg_c, nr), C(pop_c, np))))))
  mvn <- function(v, S) { L <- chol(S); z <- backsolve(L, v, transpose = TRUE); -0.5 * (length(v) * log(2 * pi) + 2 * sum(log(diag(L))) + sum(z^2)) }
  expect_equal(nll, -mvn(as.vector(ln_NAA), S), tolerance = 1e-8)

  # the same density as one dseparable call, a dmvnorm factor on each unstructured dim
  f_us <- function(Cm) { force(Cm); function(v) RTMB::dmvnorm(v, Sigma = Cm, log = TRUE) }
  f_ar <- function(r) { force(r); function(v) RTMB::dautoreg(v, phi = r, log = TRUE) }
  sep <- RTMB::dseparable(f_us(C(pop_c, np)), f_us(C(reg_c, nr)), f_ar(rho_y), f_us(C(seas_c, nk)), f_ar(rho_a), f_us(C(sex_c, ns)))
  expect_equal(nll, -sep(ln_NAA, scale = sqrt(marginal2)), tolerance = 1e-8)

})
