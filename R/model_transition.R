# Stage 2 of 3: objective function
#
# The lowest level of the population dynamics: what happens to a vector of
# abundance over one season. Every module that moves fish through time goes
# through here, which is what keeps movement timing consistent between the
# estimation model, the reference points and the operating model.
#
# build_seas_operator and advance_seas step abundance forward; spawn_state stops
# partway through a season for spawning; integrate_seas_abundance and
# catch_at_age give the within season averages the catch equation needs;
# survey_state gives abundance at the survey timing.

#' Build a seasonal transition operator combining movement and survival
#'
#' Returns the linear operator that advances a numbers-at-age vector across one
#' season, folding together movement and total mortality according to
#' \code{move_timing}. The operator is returned in the package's row-vector
#' convention, i.e. \code{N_new = t(N) \%*\% T} where \code{T[from, to]}, matching
#' the storage convention of \code{Movement}.
#'
#' @param Move Square \code{[n_regions x n_regions]} movement matrix in row
#'   convention (\code{Move[from, to]}), as stored in \code{Movement}. Required
#'   for \code{move_timing} 0 and 1; ignored for \code{move_timing = 2}.
#' @param Z Numeric vector of length \code{n_regions} giving total mortality for
#'   this season, already scaled by season duration (i.e. the same quantity
#'   stored in \code{ZAA}). Pass zeros for a movement-only step.
#' @param Q Square \code{[n_regions x n_regions]} instantaneous rate matrix
#'   (generator) in row convention, as stored in \code{Mrate}. Required when
#'   \code{move_timing = 2}; ignored otherwise.
#' @param dur Season duration used to scale \code{Q}. Should be
#'   \code{seasdur[seas]}. Only used when \code{move_timing = 2}.
#' @param move_timing Integer flag for the movement/mortality ordering:
#'   \code{0} = movement then mortality (default, historical SPoRC behaviour),
#'   \code{1} = mortality then movement, \code{2} = continuous (simultaneous)
#'   movement and mortality.
#' @param expm_nsub Integer controlling how the matrix exponential is evaluated
#'   under \code{move_timing = 2}: \code{0} (default) uses \code{Matrix::expm},
#'   a value \eqn{n \ge 1} uses the implicit backward Euler scheme
#'   \eqn{(I - A/n)^{-n}}. See \code{\link{mat_exp}}.
#'
#' @return A square \code{[n_regions x n_regions]} matrix in row convention.
#'
#' @details
#' Writing \eqn{s = \exp(-Z)} and letting \eqn{M} be the row-convention movement
#' matrix, the three operators are
#' \deqn{T_0 = M \, \mathrm{diag}(s), \qquad T_1 = \mathrm{diag}(s) \, M, \qquad
#'       T_2 = \left[\exp\left(Q^\top \Delta - \mathrm{diag}(Z)\right)\right]^\top}
#' where \eqn{\Delta} is \code{dur}. \eqn{Q^\top} converts the stored row-convention
#' generator back to the column convention the exponential is taken in, which is
#' \code{Matrix::expm} or the implicit solve of \code{\link{mat_exp}} according to
#' \code{expm_nsub}.
#'
#' The three agree exactly when \code{Z} is constant across regions, because a
#' scalar multiple of the identity commutes with the generator. They also agree
#' when \code{Z} is all zero (movement only) and when movement is absent
#' (\code{Move} the identity, \code{Q} all zero).
#'
#' @keywords internal
#' @import RTMB
build_seas_operator <- function(Move, Z, Q = NULL, dur = 1, move_timing = 0, expm_nsub = 0) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  n_regions <- length(Z)

  if(move_timing == 0) {
    # movement then mortality
    Move %*% diag(exp(-Z), n_regions)
  } else if(move_timing == 1) {
    # mortality then movement
    diag(exp(-Z), n_regions) %*% Move
  } else if(move_timing == 2) {
    # continuous: movement and mortality act simultaneously
    if(is.null(Q)) stop("Q (instantaneous rate matrix) is required when move_timing = 2.")
    A_ss <- t(Q) * dur - diag(Z, n_regions) # column convention generator, net of mortality
    t(mat_exp(A_ss, expm_nsub))
  } else stop("move_timing must be 0 (movement then mortality), 1 (mortality then movement), or 2 (continuous).")
}

#' Advance a numbers-at-region vector across one season
#'
#' Thin convenience wrapper around \code{build_seas_operator} that applies the
#' seasonal transition to a single numbers-at-region vector. For
#' \code{move_timing} 0 and 1 the operator is never formed explicitly, which
#' keeps the AD tape smaller than the equivalent matrix product.
#'
#' @inheritParams build_seas_operator
#' @param N Numeric vector of length \code{n_regions}.
#'
#' @return Numeric vector of length \code{n_regions}.
#'
#' @keywords internal
#' @import RTMB
advance_seas <- function(N, Move, Z, Q = NULL, dur = 1, move_timing = 0, expm_nsub = 0) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  if(move_timing == 0) {
    as.vector(t(N) %*% Move) * exp(-Z)
  } else if(move_timing == 1) {
    as.vector(t(N * exp(-Z)) %*% Move)
  } else {
    as.vector(t(N) %*% build_seas_operator(Move, Z, Q, dur, move_timing, expm_nsub))
  }
}

#' Numbers at the spawning point within a season
#'
#' Propagates a numbers-at-region vector from the start of a season to the
#' spawning point \code{t_spawn} within it, consistently with \code{move_timing}.
#' Each timing has its own natural spawning state, so no additional convention is
#' imposed:
#' \describe{
#'   \item{\code{move_timing = 0}}{Fish move at the start of the season and then
#'     experience \code{t_spawn} worth of mortality, so spawning happens at the
#'     post-movement location. This reproduces the historical SPoRC calculation.}
#'   \item{\code{move_timing = 1}}{Movement occurs at the end of the season, so
#'     spawning happens at the pre-movement location after \code{t_spawn} worth of
#'     mortality.}
#'   \item{\code{move_timing = 2}}{The population is propagated a fraction
#'     \code{t_spawn} of the way through the season under the combined
#'     movement-mortality generator, so spawners are partially redistributed.}
#' }
#'
#' @inheritParams build_seas_operator
#' @param N Numeric vector of length \code{n_regions} at the start of the season.
#' @param t_spawn Fraction of the season elapsed before spawning, in \code{[0, 1]}.
#'
#' @return Numeric vector of length \code{n_regions}.
#'
#' @keywords internal
#' @import RTMB
spawn_state <- function(N, Move, Z, Q = NULL, dur = 1, t_spawn = 0, move_timing = 0, expm_nsub = 0) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  if(move_timing == 0) {
    as.vector(t(N) %*% Move) * exp(-t_spawn * Z)
  } else if(move_timing == 1) {
    N * exp(-t_spawn * Z)
  } else {
    n_regions <- length(Z)
    A_ss <- (t(Q) * dur - diag(Z, n_regions)) * t_spawn
    as.vector(mat_exp(A_ss, expm_nsub) %*% N)
  }
}

#' Integrate abundance over a season under continuous movement
#'
#' Computes \eqn{\int_0^1 e^{A\tau} d\tau \, N(0)}, the season-integrated
#' abundance required for the spatial Baranov catch equation when movement and
#' mortality act simultaneously. Under \code{move_timing = 2} the usual
#' \eqn{F/Z \, (1 - e^{-Z}) N} form is not valid, because fish redistribute among
#' regions while they are dying.
#'
#' The integral runs over the unit interval, not over \eqn{[0, \Delta]}, because
#' \eqn{A = Q^\top \Delta - \mathrm{diag}(Z)} already carries the season duration in
#' both of its terms. The integration variable is elapsed \emph{fraction} of the
#' season, so a full season is \eqn{\tau = 1}.
#'
#' @inheritParams build_seas_operator
#' @param N Numeric vector of length \code{n_regions} at the start of the season.
#'
#' @return Numeric vector of length \code{n_regions} giving season-integrated
#'   abundance by region. Multiplying elementwise by a fishing mortality vector
#'   gives catch.
#'
#' @details
#' Van Loan's block identity
#' \deqn{\exp\left(\begin{bmatrix} A & I \\ 0 & 0\end{bmatrix}\right) =
#'       \begin{bmatrix} e^{A} & \int_0^1 e^{A\tau}d\tau \\ 0 & I\end{bmatrix}}
#' which recovers the integral from a single matrix exponential of twice the
#' dimension. This is preferred over \eqn{A^{-1}(e^{A} - I)} because it avoids an
#' explicit inverse and stays well conditioned as \eqn{A} approaches singularity.
#' That is not an edge case here: the generator has a zero eigenvalue (its rows sum
#' to zero), so \eqn{A} becomes exactly singular as total mortality approaches zero,
#' precisely where the integral itself remains perfectly well behaved.
#'
#' @keywords internal
#' @import RTMB
integrate_seas_abundance <- function(N, Z, Q, dur = 1, expm_nsub = 0) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  as.vector(seas_operator_and_integral(Z, Q, dur, expm_nsub)$Integral %*% N)
}

#' Catch-at-age over one season, consistent with the movement timing
#'
#' Applies the catch equation to a numbers-at-region vector in the way each
#' \code{move_timing} requires. Catch is taken where the fish actually are during the
#' season: at their post-movement locations under \code{move_timing = 0} (movement
#' happens at the start of the season), at their pre-movement locations under
#' \code{move_timing = 1} (movement happens at the end). Under \code{move_timing = 2}
#' fish redistribute among regions while they are being caught, so the region-local
#' \eqn{F/Z\,(1 - e^{-Z})\,N} form is invalid and the season-integrated (spatial
#' Baranov) abundance is used instead.
#'
#' @inheritParams build_seas_operator
#' @param N Numeric vector of length \code{n_regions} at the start of the season.
#' @param F_landed Numeric vector of length \code{n_regions} giving the seasonal
#'   fishing mortality to apply (retained, discarded or landed, per the caller).
#'
#' @return Numeric vector of length \code{n_regions} giving catch by region.
#'
#' @keywords internal
#' @import RTMB
catch_at_age <- function(N, Move, Z, Q = NULL, dur = 1, F_landed, move_timing = 0, expm_nsub = 0) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  if(move_timing == 2) {
    F_landed * integrate_seas_abundance(N, Z, Q, dur, expm_nsub)
  } else {
    # under timing 0 movement precedes the fishery; under timing 1 it follows it
    N_at_risk <- if(move_timing == 0) as.vector(t(N) %*% Move) else N
    N_at_risk * (F_landed / Z) * (1 - exp(-Z))
  }
}

#' Seasonal transition operator and abundance integral from one matrix exponential
#'
#' Returns both the seasonal transition operator and the season-integrated abundance
#' operator for \code{move_timing = 2}, at the cost of a single matrix exponential.
#'
#' The Van Loan block used to recover the integral already contains the transition
#' operator in its top-left corner,
#' \deqn{\exp\left(\begin{bmatrix} A & I \\ 0 & 0\end{bmatrix}\right) =
#'       \begin{bmatrix} e^{A} & \int_0^1 e^{A\tau}d\tau \\ 0 & I\end{bmatrix},}
#' so callers that need the population step \emph{and} the catch integral for the same
#' stratum, which is every fished stratum, since the dynamics advance the numbers and
#' the Baranov equation integrates them over the identical \eqn{A}, should take both
#' from here rather than exponentiating \eqn{A} once and the block again. Under
#' reverse-mode AD the adjoint of a matrix exponential is far more expensive than its
#' forward evaluation, so halving the number of exponentials on the tape is worth more
#' than the flop count alone suggests.
#'
#' @inheritParams build_seas_operator
#'
#' @return A list with
#' \describe{
#'   \item{\code{T}}{The seasonal transition operator in row convention, identical to
#'     \code{build_seas_operator(..., move_timing = 2)}.}
#'   \item{\code{Integral}}{\eqn{\int_0^1 e^{A\tau}d\tau} in column convention;
#'     \code{Integral \%*\% N} is the season-integrated abundance.}
#' }
#'
#' @keywords internal
#' @import RTMB
seas_operator_and_integral <- function(Z, Q, dur = 1, expm_nsub = 0) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  n_regions <- length(Z)
  A_ss <- t(Q) * dur - diag(Z, n_regions) # column convention generator, net of mortality

  # Van Loan block: [[A, I], [0, 0]], exponentiated to recover both blocks at once
  blk <- matrix(0, 2 * n_regions, 2 * n_regions)
  blk[1:n_regions, 1:n_regions] <- A_ss
  blk[1:n_regions, (n_regions + 1):(2 * n_regions)] <- diag(1, n_regions)
  E <- mat_exp(blk, expm_nsub)

  list(T = t(E[1:n_regions, 1:n_regions]),                          # row convention
       Integral = E[1:n_regions, (n_regions + 1):(2 * n_regions)])  # top-right block
}

#' Numbers at the survey point within a season
#'
#' Propagates a numbers-at-region vector from the start of a season to the survey
#' timing \code{t_srv} within it, consistently with \code{move_timing}. A survey
#' index is a snapshot at a point inside the season, not an accumulation over it,
#' so this uses the same partial propagation as \code{\link{spawn_state}} rather
#' than the season integral used for catch.
#'
#' @inheritParams build_seas_operator
#' @param N Numeric vector of length \code{n_regions} at the start of the season.
#' @param t_srv Numeric vector of length \code{n_regions} giving the fraction of
#'   the season elapsed before the survey, per region. A scalar is recycled.
#'
#' @return Numeric vector of length \code{n_regions}.
#'
#' @details
#' Under \code{move_timing} 0 and 1 the fish occupy one region for the whole
#' season (post-movement and pre-movement respectively), so the historical
#' \eqn{N \exp(-t_{srv} Z)} form is exact and is used unchanged.
#'
#' Under \code{move_timing = 2} the population is propagated by
#' \eqn{\exp(A \, t_{srv})} with
#' \eqn{A = Q^\top \Delta - \mathrm{diag}(Z)}.
#'
#' \code{t_srv} may differ by region, which a single propagation operator cannot
#' represent: fish observed in region \eqn{r} arrived from regions whose elapsed
#' times differ. The convention adopted here is that the survey in region \eqn{r}
#' observes the population propagated to \emph{that region's} survey time, i.e.
#' element \eqn{r} of \eqn{\exp(A \, t_{srv,r}) N}. When \code{t_srv} is constant
#' across regions (the usual case) this reduces to a single propagation and one
#' matrix exponential.
#'
#' @keywords internal
#' @import RTMB
survey_state <- function(N, Move, Z, Q = NULL, dur = 1, t_srv = 0, move_timing = 0, expm_nsub = 0) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  n_regions <- length(Z)
  if(length(t_srv) == 1) t_srv <- rep(t_srv, n_regions)

  # Movement has already been resolved for this season under the discrete timings,
  # so the fish are stationary and the elementwise discount is exact
  if(move_timing != 2) return(N * exp(-t_srv * Z))

  A_ss <- t(Q) * dur - diag(Z, n_regions) # column convention generator, net of mortality

  # Region-specific survey timing: propagate to each region's own survey time and
  # read off that region's entry
  out <- rep(0, n_regions)
  for(r in 1:n_regions) {
    out[r] <- (mat_exp(A_ss * t_srv[r], expm_nsub) %*% N)[r]
  } # end r loop
  out
}
