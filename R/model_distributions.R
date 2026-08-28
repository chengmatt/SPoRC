# Stage 2 of 3: objective function
#
# Density functions RTMB needs that base R does not provide, or does not provide
# in a form that differentiates cleanly: symmetric beta, Dirichlet, Dirichlet
# multinomial, logistic normal, and the count distributions used for tag
# recaptures.

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

#' Evaluate an age-disaggregated observation likelihood
#'
#' Every at-age observation stream in the model routes through here: retained
#' catch at age, discard catch at age, and the fishery and survey indices at age,
#' each in an aggregated and a population-specific form. They differ only in
#' which array supplies the prediction and which parameter supplies the standard
#' deviation, so the likelihood itself is written once.
#'
#' Ages within one cell may be independent, correlated as an AR(1) across ages,
#' or correlated through an unstructured matrix, the three structures ICES
#' age-structured assessments allow. A correlated cell contributes its whole
#' density to the first age present, leaving the remaining ages at zero, so the
#' returned vector still sums to the cell's contribution.
#'
#' An AR(1) is a statement about age distance, not about position in the
#' observed vector. A fleet that observes ages 2, 3, 5 and 6 has a gap, and lag
#' one across that gap is not lag one, so the covariance is built from the ages
#' themselves whenever they are not consecutive. Consecutive ages take the
#' autoregressive recursion instead, which is the same density at lower cost.
#'
#' @param obs_t Numeric vector of observations for the ages present in one cell,
#'   on the scale the fleet's likelihood is written on.
#' @param pred_t Vector of predictions, matching \code{obs_t}.
#' @param sigma Vector of standard deviations, matching \code{obs_t}.
#' @param corr_type Integer. \code{0} is \code{"iid"}, \code{1} is
#'   \code{"1dar1"}, \code{2} is \code{"us"}.
#' @param rho Correlation for \code{corr_type = 1}, on the natural scale.
#' @param ages Integer vector of the ages present, used to space the AR(1). The
#'   position in the vector is assumed when this is \code{NULL}.
#' @param corr_mat Correlation matrix for \code{corr_type = 2}, already subset to
#'   the ages present.
#'
#' @return Numeric vector the length of \code{obs_t}, the negative log
#'   likelihood contribution of each age.
#'
#' @keywords internal
get_at_age_nLL = function(obs_t, pred_t, sigma, corr_type = 0, rho = 0,
                          ages = NULL, corr_mat = NULL) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  n = length(obs_t)
  nLL = rep(0, n)

  # a single observation carries no correlation to describe
  if(corr_type == 0 || n == 1) {
    nLL = -1 * RTMB::dnorm(obs_t, pred_t, sigma, TRUE)
    return(nLL)
  }

  resid = obs_t - pred_t

  if(corr_type == 1) {
    if(is.null(ages)) ages = seq_len(n)
    if(all(diff(ages) == 1)) { # consecutive ages are the recursion itself
      nLL[1] = -1 * RTMB::dautoreg(resid, mu = 0, phi = rho, log = TRUE, scale = sigma)
      return(nLL)
    }
    corr_mat = matrix(0, n, n) # gapped ages need the distance stated explicitly
    for(i in 1:n) {
      for(j in 1:n) corr_mat[i,j] = rho^abs(ages[i] - ages[j])
    } # end j loop, end i loop
  }

  # standardising keeps the correlation matrix free of the standard deviations,
  # which is what lets an unstructured matrix be built once per fleet and reused
  z = resid / sigma
  nLL[1] = -1 * (RTMB::dmvnorm(z, 0, corr_mat, log = TRUE) - sum(log(sigma)))

  return(nLL)
}

#' Evaluate an at-age observation block correlated over both age and year
#'
#' A separable first-order autoregression treats the residual surface as an
#' AR(1) across ages and an AR(1) across years, with the covariance the Kronecker
#' product of the two. It is defined over a complete grid, so the caller supplies
#' a rectangular block and the whole block's density is returned as one number.
#'
#' The residuals are standardised before the separable density is applied, since
#' the standard deviations vary by age while the correlation does not. The
#' determinant of that scaling is added back.
#'
#' The correlations arrive untransformed and are constrained here, before the two
#' closures are defined. \code{\link[RTMB]{dseparable}} evaluates those closures
#' in a context of its own, and a constraint left as an unevaluated argument is
#' forced inside that context, where the tape reports an invalid \code{advector}
#' rather than a wrong number.
#'
#' @param resid Matrix of residuals, years by ages.
#' @param sigma Matrix of standard deviations, shaped like \code{resid}.
#' @param trans_rho_age Unconstrained correlation between adjacent ages.
#' @param trans_rho_year Unconstrained correlation between adjacent years.
#'
#' @return The negative log likelihood of the whole block, a scalar.
#'
#' @keywords internal
get_at_age_2dar1_nLL = function(resid, sigma, trans_rho_age, trans_rho_year) {

  z = resid / sigma

  rho_age = rho_trans(trans_rho_age)   # forced before the closures capture them
  rho_year = rho_trans(trans_rho_year)

  f_year = function(x) RTMB::dautoreg(x, mu = 0, phi = rho_year, log = TRUE)
  f_age = function(x) RTMB::dautoreg(x, mu = 0, phi = rho_age, log = TRUE)

  ll = RTMB::dseparable(f_year, f_age)(z) - sum(log(sigma))

  return(-1 * ll)
}

#' Combine reported index standard errors with an estimated component
#'
#' An index observation carries a standard error from its own survey design, and
#' an assessment may additionally estimate a component covering everything that
#' design does not. This returns the total standard deviation the index likelihood should use.
#'
#' Additive represents a parameter added to an existing estiamte of the standard error. Quadrature treats the two as independent
#' variances. Replacement discards the reported errors and estimates a rate.
#'
#' @param se Reported standard errors, any shape.
#' @param extra Estimated component on the natural scale, conformable with
#'   \code{se}.
#' @param form Integer. \code{0} returns \code{se} unchanged, \code{1} additive,
#'   \code{2} in quadrature, \code{3} replacement.
#'
#' @return The total standard deviation, shaped like \code{se}.
#'
#' @keywords internal
combine_idx_sd = function(se, extra, form) {
  if(form == 1) return(se + extra)
  if(form == 2) return(sqrt(se^2 + extra^2))
  if(form == 3) return(0 * se + extra) # keeps the shape of se while dropping its values
  se
}

#' Build the total index standard deviation from a per-fleet estimated component
#'
#' Applies \code{\link{combine_idx_sd}} fleet by fleet, where the fleet is the
#' last dimension of \code{se}. Handles both the aggregated arrays, indexed
#' region by year by season by fleet, and the population-specific arrays, which
#' carry a leading population dimension.
#'
#' @param se Reported standard errors, with fleet as the last dimension.
#' @param ln_sigma Log-scale estimated component, one value per fleet.
#' @param form Integer form code, see \code{\link{combine_idx_sd}}.
#'
#' @return An array shaped like \code{se}.
#'
#' @keywords internal
build_idx_sd = function(se, ln_sigma, form) {

  "[<-" <- RTMB::ADoverload("[<-")

  if(form == 0 || is.null(se)) return(se)

  d = dim(se)
  n_fleets = d[length(d)]
  out = se

  for(f in 1:n_fleets) {
    extra = exp(ln_sigma[f])
    if(length(d) == 4) out[,,,f] = combine_idx_sd(se[,,,f], extra, form)
    if(length(d) == 5) out[,,,,f] = combine_idx_sd(se[,,,,f], extra, form)
  } # end f loop

  return(out)
}

#' Evaluate an index negative log-likelihood under a chosen error structure
#'
#' Calls an abundance or biomass index to use one of three error structures.
#' \code{like_type = 0} is lognormal, where \code{sigma} is a log-scale
#' standard deviation and \code{const} is added inside both logarithms.
#' \code{like_type = 1} is normal on the arithmetic scale, where \code{sigma}
#' is an arithmetic standard deviation. \code{like_type = 2} is multivariate
#' normal on the arithmetic scale with a fixed covariance supplied through
#' \code{Sigma}, which places the whole series in a single density and is what
#' a survey-provided covariance across years calls for, since a diagonal
#' likelihood would treat correlated residuals as independent information.
#'
#'
#' @param obs Numeric vector of observed index values.
#' @param pred Numeric vector of predicted index values.
#' @param sigma Numeric vector of standard deviations. Ignored when
#'   \code{like_type = 2}.
#' @param like_type Integer. \code{0} lognormal, \code{1} normal, \code{2}
#'   multivariate normal.
#' @param Sigma Covariance matrix. Required when \code{like_type = 2}. Must be
#'   symmetric and positive definite; \code{\link[RTMB]{dmvnorm}} reads only the
#'   lower triangle and returns \code{NaN} without warning otherwise, so this is
#'   checked during setup rather than here.
#' @param const Numeric constant added inside the logarithms of the lognormal
#'   form. Default \code{0}.
#'
#' @return Numeric vector of negative log-likelihood contributions the same
#'   length as \code{obs}. The multivariate normal cannot be split across
#'   observations, so its total is returned in the first element with zeros
#'   elsewhere.
#'
#' @import RTMB
#'
#' @keywords internal
get_index_nLL = function(obs, pred, sigma, like_type, Sigma = NULL, const = 0) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  nLL = rep(0, length(obs))

  if(like_type == 0) nLL = -1 * RTMB::dnorm(log(obs + const), log(pred + const), sigma, TRUE) # lognormal
  if(like_type == 1) nLL = -1 * RTMB::dnorm(obs, pred, sigma, TRUE) # normal

  if(like_type == 2) {
    if(is.null(Sigma)) stop("A multivariate normal index likelihood requires a covariance matrix.") # multivariate normal
    nLL[1] = -1 * RTMB::dmvnorm(as.vector(obs), as.vector(pred), Sigma, log = TRUE)
  }

  return(nLL)
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
