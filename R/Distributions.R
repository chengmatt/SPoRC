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
