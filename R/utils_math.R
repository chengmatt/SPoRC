# Shared helpers
#
# Numerical helpers used inside the objective function and the operating model:
# correlation matrices, natural cubic spline weights, the correlation transform,
# and the logistic normal covariance.

#' Construct an AR(1) correlation matrix
#'
#' Builds an \eqn{n \times n} correlation matrix whose \eqn{(i,j)} element
#' equals \eqn{\rho^{|i-j|}}, corresponding to a stationary AR(1) process
#' with autocorrelation parameter \eqn{\rho}.
#'
#' @param n Integer. Matrix dimension (number of bins, ages, or lengths).
#' @param rho Numeric. AR(1) autocorrelation parameter in \eqn{(-1, 1)}.
#'   Values close to 1 produce strong positive correlation between adjacent
#'   bins; values close to 0 approach the identity matrix.
#'
#' @return Numeric \eqn{n \times n} correlation matrix.
#'
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' get_AR1_CorrMat(10, 0.5)
#' }
get_AR1_CorrMat <- function(n, rho) {
  corrMatrix <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      # Calculate the correlation based on the lag distance
      corrMatrix[i, j] <- rho^(abs(i - j))
    } # end i
  } # end j
  return(corrMatrix)
}

#' Construct a constant (exchangeable) correlation matrix
#'
#' Builds an \eqn{n \times n} correlation matrix with 1 on the diagonal and
#' \eqn{\rho} on all off-diagonal elements, corresponding to a compound
#' symmetry (exchangeable) covariance structure. Used in SPoRC to model
#' constant correlation across sexes in the 2D logistic-normal composition
#' likelihood (\code{comp_like = 4}).
#'
#' @param n Integer. Matrix dimension (typically \code{n_sexes}).
#' @param rho Numeric. Off-diagonal correlation in \eqn{(-1, 1)}. A value
#'   of 0 produces the identity matrix; a value approaching 1 produces
#'   near-perfect correlation across sexes.
#'
#' @return Numeric \eqn{n \times n} correlation matrix.
#'
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' get_Constant_CorrMat(2, 0.5)
#' }
get_Constant_CorrMat <- function(n, rho) {
  corrMatrix <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      if(i != j) corrMatrix[i, j] <- rho
      else corrMatrix[i, j] <- 1
    } # end i
  } # end j
  return(corrMatrix)
}

#' Construct a natural cubic spline interpolation weight matrix
#'
#' Precomputes a linear operator \eqn{W} such that \code{W \%*\% y_nodes}
#' reproduces natural cubic spline interpolation (zero second derivative at
#' the endpoints) of the node values \code{y_nodes} (defined at
#' \code{x_nodes}) onto the query points \code{x_out}. Because natural cubic
#' spline interpolation is linear in the node values for fixed node/query
#' positions, \eqn{W} depends only on \code{x_nodes} and \code{x_out} (both
#' treated as fixed data), never on the node values themselves. This lets
#' bicubic/cubic selectivity splines (see \code{\link{Get_Selex}},
#' \code{Selex_Model == 8}) be evaluated as a pair of matrix multiplications
#' against AD parameter vectors, rather than re-solving a spline system on
#' every function evaluation.
#'
#' @param x_nodes Numeric vector of strictly increasing node (knot)
#'   positions, length \eqn{n \ge 1}. Typically bin or year positions
#'   rescaled to \eqn{[0,1]}.
#' @param x_out Numeric vector of query positions at which the spline is to
#'   be evaluated, length \eqn{m}. Values are clamped to the innermost
#'   spline segment if they fall outside \code{range(x_nodes)} (no
#'   extrapolation).
#'
#' @return Numeric \eqn{m \times n} weight matrix \eqn{W}. When \eqn{n == 1}
#'   (a single node), \eqn{W} is a column of ones (the interpolated curve is
#'   constant, equal to the single node value). When \eqn{n == 2}, natural
#'   boundary conditions force the spline to reduce to linear interpolation.
#'
#' @details
#' Natural cubic spline second derivatives \eqn{M} at the nodes solve a
#' tridiagonal linear system \eqn{A M = R y} where \eqn{A} and \eqn{R} depend
#' only on the node spacing \code{diff(x_nodes)} (data), so
#' \eqn{M = A^{-1} R y = \text{Mmat} \, y} is itself linear in \eqn{y}.
#' Substituting into the standard piecewise-cubic evaluation formula for each
#' query point yields one row of \eqn{W} per query point, each a fixed linear
#' combination of the node basis vectors and rows of \code{Mmat}.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' x_nodes <- seq(0, 1, length.out = 4)
#' x_out <- seq(0, 1, length.out = 20)
#' W <- Get_Natural_Cubic_Spline_Weights(x_nodes, x_out)
#' y_nodes <- c(0.1, 0.8, 0.6, 0.3)
#' W %*% y_nodes # interpolated curve at x_out
#' }
Get_Natural_Cubic_Spline_Weights <- function(x_nodes, x_out) {

  n <- length(x_nodes)

  if(n == 1) return(matrix(1, nrow = length(x_out), ncol = 1)) # constant curve

  if(is.unsorted(x_nodes, strictly = TRUE)) stop("x_nodes must be strictly increasing")

  h <- diff(x_nodes) # node spacing, length n - 1

  # Tridiagonal system A %*% M = Rm %*% y for natural (M_1 = M_n = 0) second derivatives
  A <- matrix(0, n, n)
  Rm <- matrix(0, n, n)
  A[1,1] <- 1
  A[n,n] <- 1
  if(n > 2) {
    for(i in 2:(n-1)) {
      A[i, i-1] <- h[i-1]
      A[i, i]   <- 2 * (h[i-1] + h[i])
      A[i, i+1] <- h[i]

      Rm[i, i-1] <- 6 / h[i-1]
      Rm[i, i]   <- -6 * (1 / h[i-1] + 1 / h[i])
      Rm[i, i+1] <- 6 / h[i]
    } # end i loop
  } # end if

  Mmat <- solve(A, Rm) # n x n operator mapping node values y -> second derivatives M

  m <- length(x_out)
  W <- matrix(0, m, n)
  for(k in 1:m) {
    x0 <- min(max(x_out[k], x_nodes[1]), x_nodes[n]) # clamp query to node range (no extrapolation)
    i <- findInterval(x0, x_nodes, all.inside = TRUE) # clamp segment index to innermost segment
    hi <- h[i]
    t <- (x0 - x_nodes[i]) / hi
    a <- 1 - t
    b <- t

    row <- (a^3 - a) * hi^2 / 6 * Mmat[i, ] + (b^3 - b) * hi^2 / 6 * Mmat[i+1, ]
    row[i]   <- row[i]   + a
    row[i+1] <- row[i+1] + b
    W[k, ] <- row
  } # end k loop

  return(W)
}

#' Transform a real-valued parameter to the interval (-1, 1)
#'
#' Applies the scaled logistic transformation
#' \eqn{2 / (1 + e^{-2x}) - 1} to map an unconstrained real value to
#' \eqn{(-1, 1)}, suitable for parameterizing correlation coefficients.
#' Used in SPoRC to back-transform raw correlation parameters before
#' constructing AR(1) and constant covariance matrices for logistic-normal
#' composition likelihoods.
#'
#' @param x Numeric. Unconstrained real-valued parameter.
#'
#' @return Numeric. Transformed value in \eqn{(-1, 1)}.
#'
#'
#' @export rho_trans
#' @family Utility
rho_trans <- function(x){
  2/(1+ exp(-2 * x)) - 1 # constraint between -1 and 1
}

#' Build an unstructured correlation matrix from unconstrained parameters
#'
#' An unstructured correlation across ages is also somtimes desirable as a check. It
#' places no shape on how ages covary, at the cost of \eqn{n(n-1)/2} parameters.
#'
#' The parameters fill the strict lower triangle of a matrix whose diagonal is
#' one. Normalizing each row to unit length makes that matrix a Cholesky factor,
#' so the product with its transpose is a correlation matrix with ones on the
#' diagonal and is positive definite for any parameter values.
#'
#' @param pars Numeric vector of length \eqn{n(n-1)/2}, unconstrained.
#' @param n Number of ages.
#'
#' @return An \eqn{n \times n} correlation matrix.
#'
#' @keywords internal
build_us_corr <- function(pars, n) {
  L <- build_us_chol(pars, n)
  return(L %*% t(L))
}

#' Cholesky factor of an unstructured correlation matrix
#'
#' The factor behind \code{\link{build_us_corr}}, returned in its own right for
#' callers that need to whiten rather than to form the correlation. Whitening a
#' margin, \eqn{z = L^{-1} x}, is what lets a correlation over one margin compose
#' with an arbitrary structure over the others: the transformed margin is
#' independent, so each slice can then be scored by whatever density the
#' remaining margins call for, separable or not.
#'
#' Rows are normalized to unit length, which makes the lower triangular matrix a
#' Cholesky factor of a correlation matrix for any parameter values, so no
#' positive-definiteness constraint is needed on the parameters themselves.
#'
#' @param pars Numeric vector of length \eqn{n(n-1)/2}, unconstrained.
#' @param n Dimension of the margin.
#'
#' @return An \eqn{n 	imes n} lower triangular matrix with unit-length rows.
#'
#' @keywords internal
build_us_chol <- function(pars, n) {

  "[<-" <- RTMB::ADoverload("[<-")

  if(n == 1) return(matrix(1, 1, 1))

  L = matrix(0, n, n)
  L[1,1] = 1
  k = 1
  for(i in 2:n) {
    L[i,i] = 1
    for(j in 1:(i - 1)) {
      L[i,j] = pars[k]
      k = k + 1
    } # end j loop
  } # end i loop

  for(i in 1:n) L[i,] = L[i,] / sqrt(sum(L[i,]^2)) # unit rows make L a Cholesky factor

  return(L)
}

#' Construct a logistic-normal covariance matrix
#'
#' Builds the covariance matrix \eqn{\Sigma} used in logistic-normal
#' composition likelihoods for a given correlation structure. Three structures
#' are supported, matching the \code{comp_like} codes used throughout SPoRC:
#' iid (\code{2}), AR(1) across bins (\code{3}), and AR(1) across bins with
#' constant correlation across sexes via a Kronecker product (\code{4}).
#'
#' @param comp_like Integer. Covariance structure: \code{2} = iid (diagonal
#'   \eqn{\theta^2 I}), \code{3} = AR(1) across bins
#'   (\eqn{\theta^2 C_{\text{AR1}}}), \code{4} = Kronecker product of
#'   constant sex correlation and AR(1) bin correlation
#'   (\eqn{\theta^2 (C_{\text{sex}} \otimes C_{\text{AR1}})}).
#' @param n_bins Integer. Number of composition categories (ages or lengths).
#'   The resulting matrix has dimension \code{n_bins} for \code{comp_like}
#'   \code{2} and \code{3}, or \code{n_bins × n_sexes} for \code{comp_like
#'   = 4}.
#' @param n_sexes Integer. Number of sexes. Required for \code{comp_like = 4};
#'   ignored otherwise.
#' @param theta Numeric. Marginal standard deviation \eqn{\theta > 0}
#'   controlling the overall scale of \eqn{\Sigma}.
#' @param corr_b Numeric. AR(1) correlation across bins in \eqn{(-1, 1)}.
#'   Required for \code{comp_like} \code{3} and \code{4}; ignored for
#'   \code{comp_like = 2}.
#' @param corr_s Numeric. Constant (exchangeable) correlation across sexes in
#'   \eqn{(-1, 1)}. Required for \code{comp_like = 4}; ignored otherwise.
#'
#' @return Numeric covariance matrix \eqn{\Sigma} of dimension
#'   \code{n_bins × n_bins} (\code{comp_like} \code{2}, \code{3}) or
#'   \code{(n_bins × n_sexes) × (n_bins × n_sexes)} (\code{comp_like = 4}).
#'
#'
#' @export get_logistN_Sigma
#' @family Utility
#'
#' @examples
#' \dontrun{
#' # iid
#' get_logistN_Sigma(comp_like = 2, n_bins = 5, n_sexes = NULL, theta = 0.5)
#'
#' # AR(1) across bins
#' get_logistN_Sigma(comp_like = 3, n_bins = 5, n_sexes = NULL,
#'                   theta = 0.5, corr_b = 0.3)
#'
#' # AR(1) across bins x constant across sexes
#' get_logistN_Sigma(comp_like = 4, n_bins = 5, n_sexes = 2,
#'                   theta = 0.5, corr_b = 0.3, corr_s = 0.2)
#' }
get_logistN_Sigma <- function(comp_like,
                              n_bins,
                              n_sexes,
                              theta,
                              corr_b = NULL,
                              corr_s = NULL
                              ) {

  # iid
  if(comp_like == 2) Sigma <- diag(rep(theta^2, n_bins))

  # 1dar1 across
  if(comp_like == 3) {
    # Construct Sigma matrix
    LN_corr_b <- corr_b # correlation by age / length
    Sigma <- get_AR1_CorrMat(n_bins, LN_corr_b)
    Sigma <- Sigma * theta^2
  }

  # 2dar1 across
  if(comp_like == 4) {
    # Construct Sigma matrix
    LN_corr_b <- corr_b
    LN_corr_s <- corr_s
    Sigma <- kronecker(get_Constant_CorrMat(n_sexes, LN_corr_s), get_AR1_CorrMat(n_bins, LN_corr_b))
    Sigma <- Sigma * theta^2
  }

  return(Sigma)
}

#' Matrix exponential of a movement generator, exactly or by implicit solve
#'
#' Turns an instantaneous rate matrix into transition fractions. This is the single
#' point where \code{SPoRC} decides how a matrix exponential is evaluated, so that
#' \code{\link{Get_Movement}} and every \code{move_timing = 2} operator in
#' \code{model_transition.R} share one convention.
#'
#' @param A Square matrix, dense or sparse, numeric or \code{advector}. The
#'   generator whose exponential is wanted, already scaled by whatever time step
#'   the caller intends (i.e. this returns \eqn{e^{A}}, not \eqn{e^{A\Delta}}).
#' @param expm_nsub Integer. \code{0} (default) evaluates \eqn{e^{A}} with
#'   \code{Matrix::expm}. A power of two \eqn{n \ge 1} uses the implicit (backward
#'   Euler) scheme with \eqn{n} substeps (faster but loses acurracy).
#'
#' @return A plain dense matrix of the same dimension as \code{A}.
#'
#' @keywords internal
#' @import RTMB
mat_exp <- function(A, expm_nsub = 0) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  if(expm_nsub == 0) {
    if(is.matrix(A) && is.numeric(A)) A <- methods::as(A, "sparseMatrix")
    return(as.matrix(Matrix::expm(A)))
  }

  # backwards euler
  A <- as.matrix(A)
  out <- solve(diag(1, nrow(A)) - A / expm_nsub)

  reps <- expm_nsub
  while(reps > 1) {
    out <- out %*% out
    reps <- reps / 2 # keep raising to powerr till done
  } # end reps loop

  out
}

