# Operating model
#
# Random variate generators the operating model needs that base R does not
# supply: logistic normal, Dirichlet multinomial and inverse Gaussian
# recruitment.

#' Simulate from a logistic-normal distribution
#'
#' Draws a single composition vector of length \eqn{K} from a logistic-normal
#' distribution using the additive log-ratio (ALR) parameterization. The last
#' category is treated as the reference bin: the mean vector of the underlying
#' multivariate normal is \eqn{\mu_k = \log(p_k / p_K)}, \eqn{k = 1, \ldots,
#' K-1}. The draw is back-transformed via the additive softmax so the returned
#' vector sums to 1. Four covariance structures are supported, controlled by
#' \code{comp_like}.
#'
#' @param exp Numeric vector of length \eqn{K}. Expected (predicted) composition
#'   proportions; must be positive and sum to 1. The last element is used as
#'   the ALR reference bin and is not directly simulated.
#' @param pars Numeric vector of parameters governing the covariance structure.
#'   Required elements depend on \code{comp_like}:
#'   \describe{
#'     \item{\code{comp_like = 2} (iid)}{1 parameter: \code{pars[1]} = marginal
#'       standard deviation \eqn{\sigma}.}
#'     \item{\code{comp_like = 3} (AR1 by bin)}{2 parameters: \code{pars[1]} =
#'       \eqn{\sigma}, \code{pars[2]} = AR1 correlation \eqn{\rho} across
#'       bins. The covariance matrix is \eqn{\sigma^2 / (1 - \rho^2)} times the
#'       AR1 correlation matrix, with the last row/column removed.}
#'     \item{\code{comp_like = 4} (AR1 × constant sex)}{3 parameters:
#'       \code{pars[1]} = \eqn{\sigma}, \code{pars[2]} = AR1 correlation across
#'       bins, \code{pars[3]} = constant correlation across sexes. The
#'       covariance is a Kronecker product
#'       \eqn{C_{\text{sex}} \otimes C_{\text{AR1}}} scaled by
#'       \eqn{\sigma^2 / ((1 - \rho_{\text{bin}}^2)(1 - \rho_{\text{sex}}^2))},
#'       with the last row/column removed.}
#'   }
#' @param comp_like Integer. Covariance structure: \code{2} = iid,
#'   \code{3} = AR1 by bin, \code{4} = AR1 by bin with constant sex
#'   correlation.
#' @param n_sexes Integer. Number of sexes. Used only when
#'   \code{comp_like = 4} to construct the Kronecker product covariance.
#'
#' @return Numeric vector of length \eqn{K} summing to 1, representing a
#'   single draw from the logistic-normal distribution with the specified mean
#'   and covariance structure.
#'
#' @importFrom MASS mvrnorm
#' @keywords internal
rlogistnormal <- function(exp,
                          pars,
                          comp_like,
                          n_sexes
                          ) {
  # set up expected value vector
  mu <- log(exp[-length(exp)]) # remove last bin since it's known
  mu <- mu - log(exp[length(exp)]) # calculate log ratio

  # if iid logistic normal
  if(comp_like == 2) {
    Sigma <- diag(length(exp)-1) # set up sigma
    diag(Sigma) <- pars[1]^2 # input parameter
  } # end if iid logistic normal

  # if logistic normal, AR1 by bin
  if(comp_like == 3) {
    Sigma <- get_AR1_CorrMat(length(exp), pars[2]) * (pars[1]^2 / (1 - pars[2]^2))
    Sigma <- Sigma[-nrow(Sigma), -ncol(Sigma)] # remove last row and column
  } # end if iid logistic normal

  # if logistic normal, AR1 by bin, constant correlation by sex
  if(comp_like == 4) {
    Sigma <- kronecker(get_Constant_CorrMat(n_sexes, pars[3]), get_AR1_CorrMat(length(exp) / n_sexes, pars[2])) * (pars[1]^2 / (1 - pars[2]^2) / (1 - pars[3]^2))
    Sigma <- Sigma[-nrow(Sigma), -ncol(Sigma)] # remove last row and column
  }

  x <- MASS::mvrnorm(1, mu, Sigma) # simulate from mvnorm (does not sum to 1) and length k
  p <- exp(x)/(1 + sum(exp(x))) # do additive transformation length k and does not sum to 1
  p <- c(p, 1 - sum(p)) # output now so it sums to 1

  return(p)
}


#' Draw samples from a Dirichlet-multinomial distribution
#'
#' Generates \code{n} independent draws from a Dirichlet-multinomial
#' distribution with total count \code{N} and concentration parameter vector
#' \code{alpha}. Each draw is produced by first sampling a Dirichlet-
#' distributed probability vector from \code{Gamma(alpha)} variates and then
#' drawing multinomial counts conditioned on that probability vector.
#'
#' The Dirichlet-multinomial arises naturally as an overdispersed alternative
#' to the multinomial: the marginal variance of each category count is
#' \eqn{N \bar{p}_k (1 - \bar{p}_k) (N + \theta) / (1 + \theta)}, where
#' \eqn{\theta = \sum \alpha_k} controls overdispersion. In SPoRC,
#' \code{alpha} is typically supplied as
#' \eqn{\exp(\ln\theta) \times N \times \hat{p}} so that
#' \eqn{\exp(\ln\theta)} is the per-observation overdispersion scalar.
#'
#' @param n Integer. Number of independent draws to generate.
#' @param N Integer. Total count per draw (multinomial sample size).
#' @param alpha Numeric vector of length \eqn{K}. Dirichlet concentration
#'   parameters. All elements must be positive. Larger values relative to
#'   \eqn{N} produce draws closer to the expected proportions
#'   \eqn{\alpha / \sum \alpha}.
#'
#' @return Integer matrix of dimensions \eqn{K \times n}. Each column is one
#'   draw: a vector of category counts summing to \code{N}.
#'
#' @importFrom stats rgamma rmultinom
#' @keywords internal
rdirM <- function(n, N, alpha) {

  # Get dirichlet draws
  rdirichlet <- function(alpha) {
    x <- stats::rgamma(length(alpha), shape = alpha, scale = 1)
    return(x / sum(x))
  }

  # Generate DM samples
  result <- replicate(n, {
    p <- rdirichlet(alpha)
    counts <- stats::rmultinom(1, size = N, prob = p)
    as.vector(counts)
  })

  return(result)
}

#' Sample recruitment from an inverse-Gaussian distribution
#'
#' Generates \code{sims} random recruitment values from an inverse-Gaussian
#' distribution whose parameters are estimated from a historical recruitment
#' vector using the method of moments. The arithmetic mean \eqn{\bar{R}_a}
#' and harmonic mean \eqn{\bar{R}_h} of \code{recruitment} are used to derive
#' the shape parameter \eqn{\delta = 1 / (\bar{R}_a / \bar{R}_h - 1)} and
#' scale parameter \eqn{\beta = \bar{R}_a}. Random draws are generated via
#' the Michael-Schucany-Haas (1976) acceptance-mixture algorithm: a squared
#' standard normal variate \eqn{\psi} yields two candidate roots \eqn{\omega}
#' and \eqn{\zeta}, and a uniform draw selects between them with probability
#' \eqn{\beta / (\beta + \omega)}.
#'
#' The inverse-Gaussian is appropriate for recruitment time series that are
#' right-skewed and strictly positive, and is used in SPoRC's resampling
#' recruitment option (\code{recruitment_opt = 999}) during projection years.
#'
#' @param sims Integer. Number of random recruitment values to generate.
#' @param recruitment Numeric vector of historical recruitment values (must be
#'   strictly positive). Used to estimate \eqn{\bar{R}_a} and \eqn{\bar{R}_h}.
#'
#' @return Numeric vector of length \code{sims} of positive random variates
#'   drawn from the fitted inverse-Gaussian distribution.
#'
#' @importFrom stats rnorm runif
#' @keywords internal
rinvgauss_rec <- function(sims,
                          recruitment
                          ) {

  a_meanRec <- mean(recruitment) # get arithmetic mean recruitment
  h_meanRec <- 1 / mean(1 / recruitment) # get harmonic mean recruitment

  # define parameters for inverse gaussian
  gamma <- a_meanRec / h_meanRec
  gi_beta <- a_meanRec
  delta <- 1 / (gamma - 1)
  cvrec <- sqrt(1 / delta)

  # Generate random variables with transformation
  psi <- stats::rnorm(sims,0,1)^2 # generate squared random normal
  omega <- gi_beta * (1 + (psi - sqrt(4 * delta * psi + psi^2)) / (2 * delta))
  zeta <- gi_beta * (1 + (psi + sqrt(4 * delta * psi + psi^2)) / (2 * delta))
  gtheta <- gi_beta / (gi_beta + omega)

  unifs <- stats::runif(sims) # generate uniform
  rv <- ifelse(unifs <= gtheta, omega, zeta) # generate random draws based on mixture

  return(rv)
}

