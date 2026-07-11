#' Evaluate a symmetric beta log-density
#'
#' Computes a log-density that penalises a parameter value toward the midpoint
#' of \code{[p_lb, p_ub]} using a symmetric beta-like kernel. The penalty
#' strengthens as \code{p_prsd} increases and diffuses as \code{p_prsd}
#' decreases. Used in SPoRC as a prior for tag reporting rates to discourage
#' values near the boundaries of the unit interval.
#'
#' @param p_val Numeric. Parameter value to evaluate; must lie in
#'   \code{(p_lb, p_ub)}.
#' @param p_ub Numeric. Upper bound of the parameter support.
#' @param p_lb Numeric. Lower bound of the parameter support.
#' @param p_prsd Numeric. Pseudo-standard-deviation controlling prior
#'   concentration. Larger values produce a stronger penalty toward the
#'   midpoint \eqn{(p_{ub} + p_{lb}) / 2}; smaller values produce a more
#'   diffuse prior.
#' @param log Logical. If \code{TRUE} (default), returns the log-density;
#'   otherwise returns the density on the probability scale.
#'
#' @return Numeric. Log-density (or density if \code{log = FALSE}).
#'
#'
#' @keywords internal
dbeta_symmetric <- function(p_val, p_ub, p_lb, p_prsd, log = TRUE) {
  # Calculate mu term
  mu <- p_prsd * log((p_ub + p_lb)/2 - p_lb) - p_prsd * log(0.5)
  # Calculate the prior likelihood components
  term1 <- mu
  term2 <- p_prsd * log(p_val - p_lb + 1e-4)
  term3 <- p_prsd * log(1 - (p_val - p_lb - 1e-4)/(p_ub - p_lb))
  # Combine terms to get final prior likelihood
  nLL <- term1 + term2 + term3
  if(log == TRUE) return(nLL) else return(exp(nLL))
}


#' Evaluate a Dirichlet log-density
#'
#' Computes the Dirichlet log-density
#' \eqn{\log \Gamma(\sum \alpha) - \sum \log \Gamma(\alpha_k) +
#' \sum (\alpha_k - 1) \log x_k} for a compositional vector \code{x} and
#' concentration parameter vector \code{alpha}. Used in SPoRC as a prior for
#' movement rates and recruitment regional apportionment.
#'
#' @param x Numeric vector of compositional values summing to 1; all elements
#'   must be strictly positive.
#' @param alpha Numeric vector of Dirichlet concentration parameters of the
#'   same length as \code{x}. Larger values of \eqn{\sum \alpha} concentrate
#'   mass near \eqn{\alpha / \sum \alpha}.
#' @param log Logical. If \code{TRUE} (default), returns the log-density;
#'   otherwise returns the density.
#'
#' @return Numeric. Log-density (or density if \code{log = FALSE}).
#'
#' @keywords internal
ddirichlet <- function(x, alpha, log = TRUE) {
  logres = lgamma(sum(alpha)) - sum(lgamma(alpha)) + sum((alpha - 1) * log(x))
  if(log == TRUE) return(logres) else return(exp(logres))
} # end function

#' Evaluate a Dirichlet-multinomial log-likelihood
#'
#' Computes the Dirichlet-multinomial log-likelihood following the
#' parameterisation of Thorson et al. (CCSRA). The concentration parameters
#' are \eqn{\alpha_k = \exp(\ln\theta) \times N \times \hat{p}_k}, so
#' \eqn{\exp(\ln\theta)} is the per-observation overdispersion scalar: values
#' near zero approach the multinomial and larger values increase variance.
#' Non-integer observed counts are supported via \code{lgamma}.
#'
#' @param obs Numeric vector of observed proportions of length \eqn{K}
#'   (need not sum exactly to 1 after \code{addtotag} offsets).
#' @param pred Numeric vector of predicted proportions of length \eqn{K};
#'   must sum to 1.
#' @param Ntotal Numeric. Total count (input sample size \eqn{N}).
#' @param ln_theta Numeric. Log overdispersion parameter. The Dirichlet
#'   concentration is \eqn{\exp(\ln\theta) \times N \times \hat{p}_k}.
#' @param give_log Logical. If \code{TRUE} (default), returns the
#'   log-likelihood; otherwise returns the likelihood.
#'
#' @return Numeric. Log-likelihood (or likelihood if \code{give_log = FALSE}).
#'
#'
#' @keywords internal
ddirmult = function(obs, pred, Ntotal, ln_theta, give_log = TRUE) {
  # Set up function variables
  n_c = length(obs) # number of categories
  p_exp = pred # expected values container
  p_obs = obs # observed values container
  dirichlet_Parm = exp(ln_theta) * Ntotal # Dirichlet alpha parameters

  # set up pdf
  logres = lgamma(Ntotal + 1)
  for(c in 1:n_c) logres = logres - lgamma(Ntotal*p_obs[c]+1) # integration constant
  logres = logres + lgamma(dirichlet_Parm) - lgamma(Ntotal+dirichlet_Parm) # 2nd term in formula

  # Summation in 3rd term in formula
  for(c in 1:n_c) {
    logres = logres + lgamma(Ntotal*p_obs[c] + dirichlet_Parm*p_exp[c])
    logres = logres - lgamma(dirichlet_Parm * p_exp[c])
  } # end c

  if(give_log == TRUE) return(logres)
  else return(exp(logres))
} # end function

#' Evaluate a logistic-normal log-likelihood
#'
#' Applies the additive log-ratio (ALR) transformation to observed and
#' predicted compositions (removing the last reference bin) and evaluates a multivariate normal log-density
#' (\code{RTMB::dmvnorm}, using a covariance matrix). The ALR mean vector
#' is \eqn{\mu_k = \log(\hat{p}_k / \hat{p}_K)}, \eqn{k = 1,\ldots,K-1}.
#'
#' @param obs Numeric vector of observed composition values of length
#'   \eqn{K}. Need not sum to 1; the last element is the ALR reference.
#' @param pred Numeric vector of predicted proportions of length \eqn{K};
#'   the last element is the ALR reference.
#' @param Sigma Covariance matrix \eqn{\Sigma}
#' @param give_log Logical. If \code{TRUE} (default), returns the
#'   log-likelihood; otherwise returns the likelihood.
#'
#' @return Numeric. Log-likelihood (or likelihood if \code{give_log = FALSE}).
#'
#' @import RTMB
#'
#' @keywords internal
dlogistnormal = function(obs, pred, Sigma, give_log = TRUE) {
  # do logistic transformation on observed values
  tmp_Obs = log(obs[-length(obs)])
  tmp_Obs = tmp_Obs - log(obs[length(obs)])
  # do logistic transformation on expected values
  mu = log(pred[-length(pred)]) # remove last bin since it's known
  mu = mu - log(pred[length(pred)]) # calculate log ratio
  res = RTMB::dmvnorm(x = as.vector(tmp_Obs), mu = as.vector(mu), Sigma = Sigma, log = give_log)
  return(res)
}

#' Evaluate a robust negative binomial log-likelihood
#'
#' Computes the negative binomial log-likelihood using a
#' \eqn{(\mu, \sigma^2 - \mu)} reparameterisation that remains valid for
#' non-integer observations via \code{lgamma}. The overdispersion parameter
#' is recovered as \eqn{k = \mu^2 / (\sigma^2 - \mu)}.
#'
#' @param x Numeric. Observed count (may be non-integer).
#' @param log_mu Numeric. Log of the mean parameter \eqn{\mu}.
#' @param log_var_minus_mu Numeric. Log of the excess variance
#'   \eqn{\sigma^2 - \mu}. Must satisfy \eqn{\sigma^2 > \mu} (i.e.,
#'   overdispersion); the implied size parameter is
#'   \eqn{k = \mu^2 / (\sigma^2 - \mu)}.
#' @param give_log Logical. If \code{TRUE} (default), returns the
#'   log-likelihood; otherwise returns the likelihood.
#'
#' @return Numeric. Log-likelihood (or likelihood if \code{give_log = FALSE}).
#'
#' @keywords internal
dnbinom_robust_noint <- function(x, log_mu, log_var_minus_mu, give_log = TRUE) {
  mu = exp(log_mu)
  var_minus_mu = exp(log_var_minus_mu)
  k = mu^2 / var_minus_mu # get overdispersion
  logres = lgamma(k+x)-lgamma(k)-lgamma(x+1)+k*log(k)-k*log(mu+k)+x*log(mu)-x*log(mu+k)
  if(give_log) return(logres) else return(exp(logres))
}


#' Evaluate a Poisson log-likelihood for non-integer counts
#'
#' Computes \eqn{-\lambda + x \log \lambda - \log \Gamma(x+1)}, which
#' reduces to the standard Poisson log-likelihood for integer \code{x} and
#' extends it continuously to non-integer values via \code{lgamma}.
#'
#' @param x Numeric. Observed count (may be non-integer).
#' @param pred Numeric. Predicted mean \eqn{\lambda > 0}.
#' @param give_log Logical. If \code{TRUE} (default), returns the
#'   log-likelihood; otherwise returns the likelihood.
#'
#' @return Numeric. Log-likelihood (or likelihood if \code{give_log = FALSE}).
#'
#' @keywords internal
dpois_noint <- function(x, pred, give_log = TRUE) {
  logres <- -pred + x*log(pred) - lgamma(x+1)
  if(give_log == TRUE) return(logres) else return(exp(logres))
}

#' Compute shape parameters for a scaled beta distribution
#'
#' Converts a mean and standard deviation expressed on the \code{[low, high]}
#' scale to the \eqn{\alpha} and \eqn{\beta} shape parameters of a beta
#' distribution defined on \code{[0, 1]} after location-scale transformation.
#' Used in SPoRC to construct beta priors for steepness (\eqn{h}) bounded to
#' \eqn{[0.2, 1]}.
#'
#' @param low Numeric. Lower bound of the parameter support.
#' @param high Numeric. Upper bound of the parameter support.
#' @param mu Numeric. Prior mean on the \code{[low, high]} scale.
#' @param sigma Numeric. Prior standard deviation on the \code{[low, high]}
#'   scale. Must satisfy \eqn{\sigma^2 < \bar{\mu}(1-\bar{\mu})} where
#'   \eqn{\bar{\mu} = (\mu - \text{low}) / (\text{high} - \text{low})}.
#'
#' @return Numeric vector \code{c(alpha, beta, low, scale)} where
#'   \code{scale = high - low}.
#'
#'
#' @keywords internal
get_beta_scaled_pars <- function(low, high, mu, sigma) {
  # convert mean and sd to alpha and beta
  scale = high - low
  mean = (mu - low) / scale
  var = (sigma / scale) ^ 2
  # stopifnot(var < mean * (1 - mean))
  var = mean * (1 - mean) / var - 1
  # alpha and beta conversion
  a = mean * var
  b = (1 - mean) * var
  return(c(a,b,low,scale))
}

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
#' AD tape. This is the R analogue of WHAM's \code{asDouble()}: it is used to
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
  x <- if (is(xobs, "osa")) xobs@x else xobs
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
  if (is(xobs, "osa")) {
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
  if (is(xobs, "osa")) {
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
  if (is(xobs, "osa")) xobs@x else xobs
}

#' Binomial CDF, P(X <= x), via the regularized incomplete beta
#'
#' \eqn{P(X \le x) = 1 - I_p(x + 1, n - x)}, matching Trijoulet et al. (2023)
#' \code{dists::pbinom}. AD-safe through \code{RTMB::pbeta}.
#'
#' @param x Count.
#' @param n Number of trials.
#' @param prob Success probability.
#' @return Lower-tail CDF value.
#' @keywords internal
osa_pbinom <- function(x, n, prob) {
  1 - RTMB::pbeta(prob, x + 1, n - x)
}

#' Keep-aware multinomial log-density for OSA residuals (cdf-capable)
#'
#' Conditional-binomial decomposition of an \eqn{A}-bin multinomial for
#' \code{\link[RTMB]{oneStepPredict}}, following Trijoulet et al. (2023). Each
#' of the first \eqn{A-1} bins is a binomial conditional on the running
#' remainder, gated by its \code{keep} element; the final bin is fixed by the
#' sum-to-\eqn{N} constraint. In addition to the density term, the analytic
#' conditional binomial CDF is accumulated through the \code{cdf_lower} /
#' \code{cdf_upper} indicators, so this density supports \strong{both}
#' \code{method = "cdf"} and \code{method = "oneStepGeneric"}.
#'
#' The running remainder is frozen to the observed total so peeling a late bin
#' cannot drive the remaining count negative.
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
      rem_fixed <- Ntot - sum(xval[seq_len(i)])
      q  <- osa_squeeze(p_x[i]) / osa_squeeze(1 - pUsed)
      q  <- osa_squeeze(q)
      # density: binomial( x_i ; nUnused , q ), nUnused frozen to remaining count
      n_i <- x[i] + rem_fixed                 # = frozen remaining total at bin i
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
#' @param obs2 Length-2 count vector \code{c(count_a, count_remaining)}.
#' @param alpha2 Length-2 concentration \code{c(alpha_a, alpha_remaining)}.
#' @return Scalar log-density.
#' @keywords internal
ddirmult2 <- function(obs2, alpha2) {
  "c" <- RTMB::ADoverload("c")
  N  <- sum(obs2)
  A0 <- sum(alpha2)
  lgamma(N + 1) - sum(lgamma(obs2 + 1)) +
    lgamma(A0) - lgamma(N + A0) +
    sum(lgamma(obs2 + alpha2) - lgamma(alpha2))
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
  "c" <- RTMB::ADoverload("c")
  x_int <- as.integer(round(osa_extract_values(x)))   # summation limit is data
  Fx <- 0
  for (i in 0:x_int) {
    Fx <- Fx + exp(ddirmult2(c(i, N - i), c(alpha, beta)))
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
#' conditional on the running remainder, gated by its \code{keep} element; the
#' analytic conditional beta-binomial CDF is accumulated through
#' \code{cdf_lower} / \code{cdf_upper}, so this density supports \strong{both}
#' \code{method = "cdf"} and \code{method = "oneStepGeneric"}.
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

  obs  <- osa_extract_x(xobs)
  kk   <- osa_extract_keep(xobs, length(alpha))
  cdfi <- osa_extract_cdf(xobs, length(alpha))
  l    <- cdfi$lower
  h    <- cdfi$upper
  oval <- osa_extract_values(xobs)          # frozen values

  A       <- length(oval)
  Ntot    <- sum(oval)
  alp_rem <- sum(alpha)

  ll <- 0
  for (a in 1:A) {
    if (a != A) {
      obs_rem_fixed <- Ntot - sum(oval[seq_len(a)])
      alp_rem       <- alp_rem - alpha[a]
      obs2   <- c(obs[a],   obs_rem_fixed)
      alpha2 <- c(alpha[a], alp_rem)
      ll <- ll + kk[a] * ddirmult2(obs2, alpha2)
      # analytic conditional beta-binomial CDF hook
      cdf <- osa_squeeze(osa_pbetabinom(obs[a], obs[a] + obs_rem_fixed, alpha[a], alp_rem))
      ll <- ll + l[a] * log(cdf) + h[a] * log(1 - cdf)
    } else {
      ll <- ll + kk[a] * 0
      cdf <- osa_squeeze(1)
      ll <- ll + l[a] * log(cdf) + h[a] * log(1 - cdf)
    }
  }
  if (log) ll else exp(ll)
}
