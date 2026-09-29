# Stage 2 of 3: objective function
#
# The quantile and randomization routines that turn a fitted composition or tag likelihood into a
# residual that is standard normal under a correct model. Only meaningful inside an OSA evaluation.

#' Squeeze a probability onto the open unit interval
#'
#' Maps values in the closed interval \eqn{[0, 1]} onto the open interval to
#' avoid boundary evaluations (e.g. \code{log(0)} or division by zero) inside
#' the conditional composition likelihoods. Identical to TMB's
#' \code{convenience.hpp} \code{squeeze()} (with \code{eps} equal to machine
#' epsilon), so this implementation agrees with WHAM's OSA composition
#' likelihood to the last digit.
#'
#' @param u Numeric or AD scalar/vector of probabilities in \eqn{[0, 1]}.
#'
#' @return The input mapped into the open interval \eqn{(eps, 1 - eps)}.
#' @keywords internal
osa_squeeze <- function(u) {
  eps <- .Machine$double.eps
  (1 - eps) * (u - 0.5) + 0.5
}

#' Extract frozen numeric values from an OSA observation
#'
#' Returns the numeric \emph{values} of an observation slice, detached from the
#' AD tape. This is the R analog of WHAM's \code{asDouble()}: it is used to
#' build the running composition remainder so that, when
#' \code{\link[RTMB]{oneStepPredict}} peels and perturbs a single bin, the
#' "all other bins" count does not move with the perturbation. Without this the
#' running remainder can be driven negative on candidate values, producing
#' \code{NaN} density evaluations (most visibly with random effects switched
#' off, where the generic integrator sweeps the full support).
#'
#' @param xobs Either an object of class \code{"osa"} (with slots \code{@x} and
#'   \code{@keep}) supplied by \code{oneStepPredict}, or a plain numeric vector
#'   used during ordinary fitting.
#'
#' @return A plain numeric vector of observed counts.
#' @keywords internal
osa_extract_values <- function(xobs) {
  x <- if (methods::is(xobs, "osa")) xobs@x else xobs
  v <- try(RTMB:::getValues(x), silent = TRUE)
  if (inherits(v, "try-error")) v <- as.numeric(x)
  as.numeric(v)
}

#' Extract the keep indicator from an OSA observation
#' @param xobs An \code{"osa"} object or plain numeric vector.
#' @param n Length for the default all-ones indicator.
#' @return Numeric/AD keep vector.
#' @keywords internal
osa_extract_keep <- function(xobs, n) {
  if (methods::is(xobs, "osa")) {
    K <- xobs@keep
    if (is.matrix(K)) K[, 1] else K
  } else rep(1, n)
}


#' Extract the lower/upper CDF indicators from an OSA observation
#'
#' Returns the \code{cdf_lower} / \code{cdf_upper} data indicators that the
#' \code{"cdf"} method of \code{\link[RTMB]{oneStepPredict}} toggles. When
#' \code{xobs} is a plain numeric vector (ordinary fitting) these are all zero,
#' so the CDF terms contribute nothing and the density reduces to the ordinary
#' likelihood.
#'
#' @param xobs An \code{"osa"} object or plain numeric vector.
#' @param n Length for the default all-zero indicators.
#' @return A list with numeric/AD vectors \code{lower} and \code{upper}.
#' @keywords internal
osa_extract_cdf <- function(xobs, n) {
  if (methods::is(xobs, "osa")) {
    K <- xobs@keep
    if (is.matrix(K) && ncol(K) >= 3) {
      return(list(lower = K[, 2], upper = K[, 3]))
    }
  }
  list(lower = rep(0, n), upper = rep(0, n))
}

#' Extract the (AD) values slot from an OSA observation
#' @param xobs An \code{"osa"} object or plain numeric vector.
#' @return AD (or numeric) values.
#' @keywords internal
osa_extract_x <- function(xobs) {
  if (methods::is(xobs, "osa")) xobs@x else xobs
}

#' Put an OSA observation slice into the order it is peeled in
#'
#' \code{\link[RTMB]{oneStepPredict}} predicts observations in the order its
#' \code{ord} attribute gives rather than in storage order, and the conditional
#' factorizations below are the right ones only when they follow that order.
#' RTMB's own \code{dmultinom} OSA method and WHAM's \code{age_comp_osa.hpp}
#' apply the same permutation. Under the ascending subset
#' \code{\link{osa_keep_subset}} builds the two orders already agree, so this is
#' the identity on every call the package makes; it is what keeps a
#' caller-supplied \code{subset} correct.
#'
#' @param xobs An \code{"osa"} object or plain numeric vector.
#' @param n Length of the slice.
#' @return Integer permutation, the identity when there is nothing to reorder.
#' @keywords internal
osa_order <- function(xobs, n) {
  if (!methods::is(xobs, "osa")) return(seq_len(n))
  ord <- attr(xobs@keep, "ord")
  if (is.null(ord)) return(seq_len(n))
  order(ord)
}

#' Binomial CDF, P(X <= x), via the regularized incomplete beta
#'
#' \eqn{P(X \le x) = 1 - I_p(x + 1, n - x)}, matching Trijoulet et al. (2023)
#' \code{dists::pbinom}. AD-safe through \code{RTMB::pbeta}.
#'
#' \code{oneStepPredict} sweeps \code{x} past \code{n}, where the CDF is one.
#' \code{pbeta} needs a positive second shape to say so, so the trials left are
#' held at zero and nudged off it.
#'
#' @param x Count.
#' @param n Number of trials.
#' @param prob Success probability.
#' @return Lower-tail CDF value.
#' @keywords internal
osa_pbinom <- function(x, n, prob) {
  rest <- n - x
  rest <- 0.5 * (rest + abs(rest))
  1 - RTMB::pbeta(prob, x + 1, rest + .Machine$double.eps)
}

#' Keep-aware multinomial log-density for OSA residuals (cdf-capable)
#'
#' Conditional-binomial decomposition of an \eqn{A}-bin multinomial for
#' \code{\link[RTMB]{oneStepPredict}}, following Trijoulet et al. (2023). Each
#' of the first \eqn{A-1} bins is a binomial conditional on the running
#' remainder, selected by its \code{keep} element; the final bin is fixed by the
#' sum-to-\eqn{N} constraint. In addition to the density term, the analytic
#' conditional binomial CDF is accumulated through the \code{cdf_lower} /
#' \code{cdf_upper} indicators, so this density supports \strong{both}
#' \code{method = "cdf"} and \code{method = "oneStepGeneric"}.
#'
#' The conditional trial count is the total left after the earlier bins, frozen
#' at the observed counts. Given what has already been peeled it is a constant,
#' so it must not move with the candidate value \code{oneStepPredict} sweeps
#' over: tying it to the candidate turns the binomial into a size-weighted
#' negative binomial and mis-calibrates the residuals, most visibly at a small
#' sample size across many bins. RTMB's own \code{dmultinom} OSA method and
#' WHAM's \code{age_comp_osa.hpp} hold it fixed the same way. Nothing after the
#' peeled bin moves either, so no remaining count can be driven negative.
#'
#' @param xobs An \code{"osa"} object from \code{oneStepPredict}, or a plain
#'   numeric count vector (length \eqn{A}) during fitting.
#' @param p Predicted proportions (length \eqn{A}); normalized internally.
#' @param log Logical; return the log-density (default) or the density.
#' @return Scalar (log-)density contribution.
#' @keywords internal
dmultinom_osa <- function(xobs, p, log = TRUE) {
  "[<-" <- RTMB::ADoverload("[<-")
  "c"   <- RTMB::ADoverload("c")

  perm <- osa_order(xobs, length(p))
  xobs <- xobs[perm]
  p    <- p[perm]

  x    <- osa_extract_x(xobs)
  kk   <- osa_extract_keep(xobs, length(p))
  cdfi <- osa_extract_cdf(xobs, length(p))
  l    <- cdfi$lower
  h    <- cdfi$upper
  xval <- osa_extract_values(xobs)          # frozen values

  A    <- length(xval)
  p_x  <- p / (sum(p) + 1e-300)
  Ntot <- sum(xval)

  logres <- 0
  pUsed  <- 0
  for (i in 1:A) {
    if (i != A) {
      q  <- osa_squeeze(p_x[i]) / osa_squeeze(1 - pUsed)
      q  <- osa_squeeze(q)
      # trials left after the earlier bins: data, so the candidate value cannot move it
      n_i <- Ntot - sum(xval[seq_len(i - 1)])
      logres <- logres + kk[i] * RTMB::dbinom(x[i], size = n_i, prob = q, log = TRUE)
      # analytic conditional CDF hook (Trijoulet dists::pbinom)
      cdf <- osa_squeeze(osa_pbinom(x[i], n_i, q))
      logres <- logres + l[i] * log(cdf) + h[i] * log(1 - cdf)
      pUsed  <- osa_squeeze(pUsed + p_x[i])
    } else {
      logres <- logres + kk[i] * 0            # last bin: density fixed by sum-to-N
      cdf <- osa_squeeze(1)                    # CDF of a determined bin is 1
      logres <- logres + l[i] * log(cdf) + h[i] * log(1 - cdf)
    }
  }
  if (log) logres else exp(logres)
}

#' Two-category Dirichlet-multinomial (beta-binomial) log-density
#'
#' Taken on the number of trials rather than on a pair of counts, so that the
#' trials stay fixed while \code{\link[RTMB]{oneStepPredict}} sweeps \code{x}
#' over its support. A candidate above \code{size} carries no mass, and
#' \code{-lgamma(rest + 1)} is what says so: it is the one term with a pole
#' there, so the shape term is held at zero rather than allowed to answer with a
#' second pole that would cancel it and leave \code{NaN} behind.
#'
#' @param x Count in the first category.
#' @param size Number of trials.
#' @param shape1,shape2 Concentration of the first category and of the rest.
#' @return Scalar log-density.
#' @keywords internal
osa_dbetabinom <- function(x, size, shape1, shape2) {
  rest     <- size - x
  rest_pos <- 0.5 * (rest + abs(rest))
  lgamma(size + 1) - lgamma(x + 1) - lgamma(rest + 1) +
    lgamma(shape1 + shape2) - lgamma(size + shape1 + shape2) +
    lgamma(x + shape1) - lgamma(shape1) +
    lgamma(rest_pos + shape2) - lgamma(shape2)
}

#' Beta-binomial CDF, P(X <= x), by summation
#'
#' Sums two-category Dirichlet-multinomial (beta-binomial) probabilities over
#' \eqn{0, \ldots, x}, matching Trijoulet et al. (2023) \code{pbetabinom}.
#' Clamped at 1 to guard against floating-point overshoot.
#'
#' @param x Count (upper summation limit; taken from frozen observed values).
#' @param N Number of trials.
#' @param alpha First beta-binomial shape.
#' @param beta Second beta-binomial shape.
#' @return Lower-tail CDF value.
#' @keywords internal
osa_pbetabinom <- function(x, N, alpha, beta) {
  x_int <- as.integer(round(osa_extract_values(x)))   # summation limit is data
  Fx <- 0
  for (i in 0:x_int) {
    Fx <- Fx + exp(osa_dbetabinom(i, N, alpha, beta))
  }
  # clamp to 1 without an AD comparison
  Fx <- Fx - osa_squeeze(0) * 0   # no-op to keep AD class
  pmin_ad <- function(a, b) 0.5 * (a + b - abs(a - b))
  pmin_ad(Fx, 1)
}

#' Keep-aware Dirichlet-multinomial log-density for OSA residuals (cdf-capable)
#'
#' Conditional beta-binomial decomposition of an \eqn{A}-bin
#' Dirichlet-multinomial for \code{\link[RTMB]{oneStepPredict}}, following
#' Trijoulet et al. (2023). Each of the first \eqn{A-1} bins is a beta-binomial
#' conditional on the running remainder, selected by its \code{keep} element; the
#' analytic conditional beta-binomial CDF is accumulated through
#' \code{cdf_lower} / \code{cdf_upper}, so this density supports \strong{both}
#' \code{method = "cdf"} and \code{method = "oneStepGeneric"}.
#'
#' The conditional trial count is held fixed exactly as it is under
#' \code{\link{dmultinom_osa}}, and for the same reason.
#'
#' @param xobs An \code{"osa"} object from \code{oneStepPredict}, or a plain
#'   numeric count vector (length \eqn{A}) during fitting.
#' @param alpha Concentration parameters (length \eqn{A}). Typically
#'   \eqn{\alpha = \hat{p}\,\exp(\ln\theta)\,N_{total}}; must match the fitting
#'   parameterization.
#' @param log Logical; return the log-density (default) or the density.
#' @return Scalar (log-)density contribution.
#' @keywords internal
ddirmult_osa <- function(xobs, alpha, log = TRUE) {
  "[<-" <- RTMB::ADoverload("[<-")
  "c"   <- RTMB::ADoverload("c")

  perm  <- osa_order(xobs, length(alpha))
  xobs  <- xobs[perm]
  alpha <- alpha[perm]

  obs  <- osa_extract_x(xobs)
  kk   <- osa_extract_keep(xobs, length(alpha))
  cdfi <- osa_extract_cdf(xobs, length(alpha))
  l    <- cdfi$lower
  h    <- cdfi$upper
  oval <- osa_extract_values(xobs)          # frozen values

  A       <- length(oval)
  Ntot    <- sum(oval)
  alp_rem <- sum(alpha)

  loglik <- 0
  for (a in 1:A) {
    if (a != A) {
      # trials left after the earlier bins: data, so the candidate value cannot move it
      n_a     <- Ntot - sum(oval[seq_len(a - 1)])
      alp_rem <- alp_rem - alpha[a]
      loglik  <- loglik + kk[a] * osa_dbetabinom(obs[a], n_a, alpha[a], alp_rem)
      # analytic conditional beta-binomial CDF hook
      cdf <- osa_squeeze(osa_pbetabinom(obs[a], n_a, alpha[a], alp_rem))
      loglik <- loglik + l[a] * log(cdf) + h[a] * log(1 - cdf)
    } else {
      loglik <- loglik + kk[a] * 0
      cdf <- osa_squeeze(1)
      loglik <- loglik + l[a] * log(cdf) + h[a] * log(1 - cdf)
    }
  }
  if (log) loglik else exp(loglik)
}
