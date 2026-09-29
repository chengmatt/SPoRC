# The conditional factorizations behind the one-step-ahead composition residuals. Each peeled bin
# has to be the exact conditional of the likelihood being fit at every candidate value
# oneStepPredict sweeps over, not only at the observed one, which is why these run on a tape: off
# the tape the frozen observed counts and the candidate are the same numbers and nothing is tested.

library(SPoRC)
library(testthat)

test_that("the one-step-ahead conditional densities", {

  p     <- c(0.30, 0.25, 0.18, 0.12, 0.09, 0.06)
  x     <- c(28, 22, 19, 14, 11, 6)
  N     <- sum(x)
  A     <- length(p)
  alpha <- p * 0.7 * N

  # an observation with one bin peeled, as oneStepPredict hands it to the density
  peel <- function(values, active, ord = NULL) {
    keep <- matrix(0, length(values), 1)
    keep[active, 1] <- 1
    if (!is.null(ord)) attr(keep, "ord") <- ord
    methods::new("osa", x = values, keep = keep)
  }

  # the density of the peeled bin as a function of the candidate, with the rest of the
  # composition frozen at the observed counts the way the tape freezes them
  conditional <- function(dens, active) {
    RTMB::MakeTape(function(candidate) {
      "[<-" <- RTMB::ADoverload("[<-")
      values <- RTMB::advector(x)
      values[active] <- candidate
      dens(peel(values, active))
    }, x[active])
  }

  # trials left for a bin once the earlier ones are spent
  trials <- function(a) N - sum(x[seq_len(a - 1)])

  # the two conditionals, written out rather than called from the package
  ref_binom <- function(k, a) stats::dbinom(k, trials(a), p[a] / sum(p[a:A]), log = TRUE)
  ref_betabinom <- function(k, a) {
    n <- trials(a)
    s1 <- alpha[a]
    s2 <- sum(alpha[seq.int(a + 1, A)])
    lchoose(n, k) + lbeta(k + s1, n - k + s2) - lbeta(s1, s2)
  }

  test_that("with nothing peeled they are the ordinary densities", {
    expect_equal(dmultinom_osa(x, p), stats::dmultinom(x, prob = p, log = TRUE))
    expect_equal(ddirmult_osa(x, alpha),
                 lgamma(N + 1) - sum(lgamma(x + 1)) +
                   lgamma(sum(alpha)) - lgamma(N + sum(alpha)) +
                   sum(lgamma(x + alpha) - lgamma(alpha)))
  })

  test_that("every candidate on the support gives the exact conditional", {
    for (a in seq_len(A - 1)) {
      support <- 0:trials(a)
      mn <- conditional(function(o) dmultinom_osa(o, p), a)
      dm <- conditional(function(o) ddirmult_osa(o, alpha), a)
      expect_equal(vapply(support, mn, numeric(1)), ref_binom(support, a),
                   label = sprintf("multinomial bin %d", a))
      expect_equal(vapply(support, dm, numeric(1)), ref_betabinom(support, a),
                   label = sprintf("Dirichlet-multinomial bin %d", a))
    }
  })

  test_that("a candidate past the last trial carries no mass", {
    for (a in seq_len(A - 1)) {
      beyond <- trials(a) + c(1, 5, 40)
      mn <- conditional(function(o) dmultinom_osa(o, p), a)
      dm <- conditional(function(o) ddirmult_osa(o, alpha), a)
      expect_equal(vapply(beyond, mn, numeric(1)), rep(-Inf, 3),
                   label = sprintf("multinomial bin %d", a))
      expect_equal(vapply(beyond, dm, numeric(1)), rep(-Inf, 3),
                   label = sprintf("Dirichlet-multinomial bin %d", a))
    }
  })

  test_that("whole-numbered concentrations do not turn that into NaN", {
    # the pole in -lgamma(rest + 1) is what refuses the candidate, and a whole-numbered
    # concentration for the remaining bins is where a second pole could cancel it
    whole <- c(20, 16, 12, 8, 6, 4)
    dm <- RTMB::MakeTape(function(candidate) {
      "[<-" <- RTMB::ADoverload("[<-")
      values <- RTMB::advector(x)
      values[2] <- candidate
      ddirmult_osa(peel(values, 2), whole)
    }, x[2])
    expect_equal(vapply(trials(2) + c(1, 5, 40), dm, numeric(1)), rep(-Inf, 3))
  })

  test_that("the residuals are the ones RTMB's own multinomial method gives", {
    residual <- function(dens) {
      dat <- list(obs = x)
      f <- function(par) {
        RTMB::getAll(dat, par)
        obs <- RTMB::OBS(obs)
        -dens(obs)
      }
      obj <- RTMB::MakeADFun(f, list(dummy = 0), silent = TRUE)
      RTMB::oneStepPredict(obj, observation.name = "obs", method = "oneStepGeneric",
                           discrete = TRUE, discreteSupport = 0:N, subset = seq_len(A - 1),
                           seed = 1, trace = FALSE)$residual
    }
    expect_equal(residual(function(o) dmultinom_osa(o, p, log = TRUE)),
                 residual(function(o) RTMB::dmultinom(o, prob = p, log = TRUE)),
                 tolerance = 1e-8)
  })

  test_that("the bins are taken in the order they are peeled in", {
    expect_equal(osa_order(x, A), seq_len(A))
    expect_equal(osa_order(peel(x, 2, ord = seq_len(A)), A), seq_len(A))
    expect_equal(osa_order(peel(x, 2, ord = rev(seq_len(A))), A), rev(seq_len(A)))

    # peeling bin 2 under a reversed order is peeling bin A - 1 of the reversed composition
    expect_equal(dmultinom_osa(peel(x, 2, ord = rev(seq_len(A))), p),
                 dmultinom_osa(peel(rev(x), A - 1), rev(p)))
    expect_equal(ddirmult_osa(peel(x, 2, ord = rev(seq_len(A))), alpha),
                 ddirmult_osa(peel(rev(x), A - 1), rev(alpha)))
  })

})
