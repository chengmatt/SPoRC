#' Calculate Selectivity
#'
#' Computes selectivity-at-bin using a suite of parametric, semi-parametric,
#' and non-parametric formulations. Supports constant, time-varying, and
#' fully flexible selectivity structures.
#'
#' @param Selex_Model Integer specifying the selectivity model:
#'   \describe{
#'     \item{0}{Logistic (b50, slope): \eqn{1 / (1 + \exp(-k(\text{bin} - b_{50})))}}
#'     \item{1}{Gamma-shaped dome (bin-at-peak \eqn{b_{\max}}, curvature \eqn{\delta}).}
#'     \item{2}{Power function (monotonic decreasing): \eqn{1 / \text{bin}^{\text{power}}}.}
#'     \item{3}{Logistic (b50, b95 parameterization).}
#'     \item{4}{Double-normal dome with plateau and flexible tails (6 parameters).}
#'     \item{5}{Non-parametric selectivity: bin-level logit parameters mapped via \code{plogis},
#'              optionally modified by time-varying deviations.}
#'     \item{6}{Logistic selectivity with asymptote:
#'              \eqn{\alpha / (1 + \exp(-k(\text{bin} - b_{50})))}.
#'              Allows maximum selectivity \eqn{\alpha \in (0,1)}.}
#'     \item{7}{Logistic selectivity with asymptote (b50, b95 parameterization):
#'              \eqn{\alpha / (1 + 19^{(b_{50} - \text{bin})/b_{95}})}.
#'              Equivalent to Model 3 scaled by asymptote \eqn{\alpha}.}
#'     \item{8}{Bicubic spline over a bin-node x year-node grid (see
#'              \code{Wbin_bicubic}, \code{Wyr_bicubic}). One generalized form
#'              covers three cases depending on how the caller constructs the
#'              node grid and interpolation weights: a single smooth bin x
#'              year surface (bicubic), a time-invariant bin-only spline
#'              (\code{n_yr_nodes == 1}), or a bin-only spline re-fit
#'              independently within each of several year blocks
#'              (\code{n_yr_nodes == 1} within each of SPoRC's existing
#'              selectivity blocks). No \code{TimeVary_Model} deviation
#'              layering applies to this model.}
#'   }
#'
#' @param TimeVary_Model Integer specifying temporal structure:
#'   \describe{
#'     \item{0}{No time variation.}
#'     \item{1}{IID deviations applied multiplicatively to model parameters.}
#'     \item{2}{Random walk deviations applied multiplicatively to model parameters.}
#'     \item{3}{3D GMRF (marginal variance): deviations applied multiplicatively at bin level.}
#'     \item{4}{3D GMRF (conditional variance): deviations applied multiplicatively at bin level.}
#'     \item{5}{Separable 2D AR(1): deviations applied multiplicatively at bin level.}
#'   }
#'
#' @param pars Numeric vector of log-scale selectivity parameters.
#' Parameters are exponentiated or transformed depending on model specification:
#'   \describe{
#'     \item{Model 0}{\code{c(ln_b50, ln_slope)}}
#'     \item{Model 1}{\code{c(ln_bmax, ln_delta)}}
#'     \item{Model 2}{\code{c(ln_power)}}
#'     \item{Model 3}{\code{c(ln_b50, ln_b95)}}
#'     \item{Model 4}{\code{c(p1, p2, p3, p4, p5, p6)}}
#'     \item{Model 5}{\code{c(logit_sel_1, ..., logit_sel_nbins)}}
#'     \item{Model 6}{\code{c(logit_alpha, ln_b50, ln_k)}}
#'     \item{Model 7}{\code{c(logit_alpha, ln_b50, ln_b95)}}
#'     \item{Model 8}{Flattened bin-node x year-node log-selectivity grid,
#'       length \code{ncol(Wyr_bicubic) * ncol(Wbin_bicubic)}, filled
#'       column-major into a \code{ncol(Wyr_bicubic)} (rows, year nodes) by
#'       \code{ncol(Wbin_bicubic)} (columns, bin nodes) matrix.}
#'   }
#'
#' @param ln_seldevs Array of log-scale selectivity deviations with dimension
#'   \code{[n_regions, n_years, n_parameters_or_bins, n_sexes, 1]}.
#'
#'   For \code{TimeVary_Model = 1–2}: deviations apply to selectivity parameters
#'   on the natural scale after exponentiation.
#'
#'   For \code{TimeVary_Model = 3–5}: deviations apply multiplicatively at the
#'   bin level to the constructed selectivity curve.
#'
#'   For \code{Selex_Model = 5}: deviations act directly on bin-level logit
#'   selectivity parameters prior to logistic transformation.
#'
#' @param Region Integer region index.
#' @param Year Integer year index (absolute, i.e. a row index into
#'   \code{Wyr_bicubic}). Only used directly by \code{Selex_Model == 8};
#'   otherwise only used to index \code{ln_seldevs}.
#' @param Bin Numeric vector of bins (ages or lengths).
#' @param Sex Integer sex index.
#' @param Wbin_bicubic Numeric \code{length(Bin) x n_bin_nodes} natural cubic
#'   spline weight matrix (see \code{\link{Get_Natural_Cubic_Spline_Weights}}),
#'   mapping bin-node log-selectivity values onto \code{Bin}. Only used when
#'   \code{Selex_Model == 8}; ignored (may be \code{NULL}) otherwise. Zero
#'   padding in unused columns (e.g. when a shared parameter array is padded
#'   to a common width across fleets/blocks) contributes nothing, since it is
#'   multiplied through to zero.
#' @param Wyr_bicubic Numeric \code{n_yrs_total x n_yr_nodes} interpolation
#'   weight matrix mapping year-node log-selectivity values onto every
#'   absolute model year (rows beyond the fitted block are typically
#'   constructed to hold the boundary node constant). Row \code{Year} is used
#'   for this call. Only used when \code{Selex_Model == 8}; ignored (may be
#'   \code{NULL}) otherwise. A single-column matrix of all-1s (\code{n_yr_nodes
#'   == 1}) yields a time-invariant bin-only spline, since every year maps
#'   onto the same single node.
#' @param n_bin_nodes_bicubic,n_yr_nodes_bicubic Integer. This fleet/block's
#'   own true number of bin nodes / year nodes. Only used when
#'   \code{Selex_Model == 8}. \strong{Must} be supplied whenever
#'   \code{Wbin_bicubic}/\code{Wyr_bicubic} may have been zero-padded wider
#'   than this specific block's own grid (e.g. because some \emph{other}
#'   fleet/block shares the same padded storage array but uses a larger
#'   bicubic grid) -- \code{ncol(Wbin_bicubic)}/\code{ncol(Wyr_bicubic)} give
#'   the padded (shared) width, not this block's true node counts, and using
#'   the padded width to reshape \code{pars} would misassign which flattened
#'   parameter values land in which (bin-node, year-node) cell. Default
#'   \code{NULL} falls back to \code{ncol(Wbin_bicubic)}/\code{ncol(Wyr_bicubic)}
#'   for backward compatibility when no padding-width mismatch is possible
#'   (e.g. a single bicubic block/fleet, or direct unit testing).
#'
#' @return Numeric vector of selectivity values corresponding to \code{Bin}.
#' Values are on the natural scale and are not normalized unless specified
#' in downstream components.
#'
#' @details
#'
#' For \code{TimeVary_Model = 0}, only the base parametric form is evaluated.
#'
#' For \code{TimeVary_Model = 1–2}, deviations modify model parameters
#' multiplicatively on the natural scale after transformation.
#'
#' For \code{TimeVary_Model = 3–5}, deviations act directly on the constructed
#' selectivity curve as multiplicative log-normal perturbations:
#' \deqn{\text{selex} = \text{selex} \cdot \exp(\delta_{r,y,b,s})}.
#'
#' For \code{Selex_Model = 5}, selectivity is fully non-parametric:
#' bin-specific logit parameters are optionally adjusted by time-varying
#' deviations and transformed via:
#' \deqn{\text{selex}_b = \text{logit}^{-1}(\eta_b)}.
#'
#' Models 6 and 7 extend logistic selectivity by introducing an asymptote
#' parameter \eqn{\alpha \in (0,1)} that allows selectivity to saturate below
#' full vulnerability.
#'
#'
#' @keywords internal
Get_Selex = function(Selex_Model,
                     TimeVary_Model,
                     pars,
                     ln_seldevs,
                     Region,
                     Year,
                     Bin,
                     Sex,
                     Wbin_bicubic = NULL,
                     Wyr_bicubic = NULL,
                     n_bin_nodes_bicubic = NULL,
                     n_yr_nodes_bicubic = NULL) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  selex = rep(0, length(Bin)) # Temporary container vector

  if(Selex_Model == 0) { # logistic selectivity (b50 and slope)
    # Extract out and exponentiate the parameters here
    b50 = exp(pars[1]); # b50
    k = exp(pars[2]); # slope

    if(TimeVary_Model %in% c(1:2)) {
      b50 = b50 * exp(ln_seldevs[Region, Year, 1, Sex, 1]) # b50 parameter varying
      k = k * exp(ln_seldevs[Region, Year, 2, Sex, 1]) # slope parameter varying
    } # end if iid or random walk

    selex = 1 / (1 + exp(-k * (Bin - b50))) # return parmetric form
  }

  if(Selex_Model == 1) { # gamma dome-shaped selectivity
    # Extract out and exponentiate the parameters here
    bmax = exp(pars[1]) # Bin at max selex
    delta = exp(pars[2]) # slope parameter

    if(TimeVary_Model %in% c(1:2)) {
      bmax = bmax * exp(ln_seldevs[Region, Year, 1, Sex, 1]) # bmax parameter varying
      delta = delta * exp(ln_seldevs[Region, Year, 2, Sex, 1]) # delta parameter varying
    } # end if iid or random walk

    # Now, calculate/derive power parameter + selex values
    p = 0.5 * (sqrt( bmax^2 + (4 * delta^2)) - bmax)
    selex = (Bin / bmax)^(bmax/p) * exp( (bmax - Bin) / p ) # return parametric form
  }

  if(Selex_Model == 2) { # power function selectivity
    # Extract out and exponentiate the parameters here
    power = exp(pars[1]); # power parameter

    if(TimeVary_Model %in% c(1:2)) {
      power = power * exp(ln_seldevs[Region, Year, 1, Sex, 1]) # power parameter varying
    } # end if iid or random walk

    selex = 1 / Bin^power # return parametric form
  }

  if(Selex_Model == 3) { # logistic selectivity (b50 and b95)

    # Extract out and exponentiate the parameters here
    b50 = exp(pars[1]); # b50
    b95 = exp(pars[2]); # b95

    if(TimeVary_Model %in% c(1:2)) {
      b50 = b50 * exp(ln_seldevs[Region, Year, 1, Sex, 1]) # b50 parameter varying
      b95 = b95 * exp(ln_seldevs[Region, Year, 2, Sex, 1]) # b95 parameter varying
    } # end if iid or random walk

    selex = 1 / (1+19^((b50-Bin)/b95)) # 19 b/c 0.95 / (1 - 0.95) return parametric form
  }

  if(Selex_Model == 4) {

    # define bin ranges for double normal here
    midbin <- Bin

    # Extract and transform parameters here
    p1trans <- min(Bin) + (max(Bin) - min(Bin)) * RTMB::plogis(pars[1]) # peak bin at plateau
    p2trans <- p1trans + 1 + (0.99 + max(Bin) - p1trans - 1)/(1 + exp(-1.0 * pars[2])) # width of plateau
    p3trans <- exp(pars[3]) # ascending width
    p4trans <- exp(pars[4]) # descending width
    p5trans <- 1/(1 + exp(-1.0 * pars[5])) # selectivity at first bin
    p6trans <- 1/(1 + exp(-1.0 * pars[6])) # selectivity at last bin

    if(TimeVary_Model %in% c(1:2)) {
      p1trans = p1trans * exp(ln_seldevs[Region, Year, 1, Sex, 1]) # p1 parameter varying
      p2trans = p2trans * exp(ln_seldevs[Region, Year, 2, Sex, 1]) # p2 parameter varying
      p3trans = p3trans * exp(ln_seldevs[Region, Year, 3, Sex, 1]) # p3 parameter varying
      p4trans = p4trans * exp(ln_seldevs[Region, Year, 4, Sex, 1]) # p4 parameter varying
      p5trans = p5trans * exp(ln_seldevs[Region, Year, 5, Sex, 1]) # p5 parameter varying
      p6trans = p6trans * exp(ln_seldevs[Region, Year, 6, Sex, 1]) # p6 parameter varying
    } # end if iid or random walk

    # construct selectivity function
    asc <- exp(-((midbin - p1trans)^2/p3trans))
    asc.scaled <- (p5trans + (1 - p5trans) * (asc - 0)/(1 - 0))
    desc <- exp(-((midbin - p2trans)^2/p4trans))
    stj <- exp(-((40 - p2trans)^2/p4trans))
    des.scaled <- (1 + (p6trans - 1) * (desc - 1) /(stj - 1))
    join1 <- 1/(1 + exp(-(20 * (midbin - p1trans)/(1 + abs(midbin - p1trans))))) # joiner functions
    join2 <- 1/(1 + exp(-(20 * (midbin - p2trans)/(1 + abs(midbin - p2trans))))) # joiner functions
    selex <- asc.scaled * (1 - join1) + join1 * ((1 - join2) + des.scaled * join2) # return parameteric form
    selex[1] <- p5trans # return parameteric form
  }

  if(Selex_Model == 5) {
    if(TimeVary_Model %in% c(1:2)) {
      pars = pars + ln_seldevs[Region, Year, , Sex, 1] # non-parametric parameters varying
    } # end if iid or random walk
    selex = RTMB::plogis(pars) # non parametric parameters
  } # end if non-parametric selectivity

  if(Selex_Model == 6) {

    alpha = RTMB::plogis(pars[1]) # asymptote parameter
    b50 = exp(pars[2]) # b50 parameter
    k = exp(pars[3]) # slope parameter

    if(TimeVary_Model %in% c(1:2)) {
      alpha = alpha * exp(ln_seldevs[Region, Year, 1, Sex, 1]) # alpha parameter varying
      b50 = b50 * exp(ln_seldevs[Region, Year, 2, Sex, 1]) # b50 parameter varying
      k = k * exp(ln_seldevs[Region, Year, 3, Sex, 1]) # slope parameter varying
    } # end if iid or random walk

    selex = alpha / (1 + exp(-k * (Bin - b50))) # return parmetric form

  } # end if logistic selex with an asymptotic parameter, w/ b50 k parameterization

  if(Selex_Model == 7) {

    alpha = RTMB::plogis(pars[1]) # asymptote parameter
    b50 = exp(pars[2]) # b50 parameter
    b95 = exp(pars[3]) # b95 parameter

    if(TimeVary_Model %in% c(1:2)) {
      alpha = alpha * exp(ln_seldevs[Region, Year, 1, Sex, 1]) # alpha parameter varying
      b50 = b50 * exp(ln_seldevs[Region, Year, 2, Sex, 1]) # b50 parameter varying
      b95 = b95 * exp(ln_seldevs[Region, Year, 3, Sex, 1]) # slope parameter varying
    } # end if iid or random walk

    selex = alpha / (1+19^((b50-Bin)/b95)) # 19 b/c 0.95 / (1 - 0.95) return parametric form

  } # end if logistic selex with an asymptotic parameter, w/ b50 b95 parameterization

  if(Selex_Model == 8) { # bicubic spline (bin node x year node grid); generalizes bicubic, bin-only, and blocked-bin-only cubic spline forms

    n_bin_nodes = if(is.null(n_bin_nodes_bicubic)) ncol(Wbin_bicubic) else n_bin_nodes_bicubic
    n_yr_nodes = if(is.null(n_yr_nodes_bicubic)) ncol(Wyr_bicubic) else n_yr_nodes_bicubic

    # flattened log-selectivity grid -> (n_yr_nodes rows) x (n_bin_nodes cols) matrix
    node_par = matrix(pars[1:(n_yr_nodes * n_bin_nodes)], nrow = n_yr_nodes, ncol = n_bin_nodes)

    # trim Wbin_bicubic/Wyr_bicubic down to this block's true node-column width before
    # multiplying, in case they were zero-padded wider than n_bin_nodes/n_yr_nodes above
    bin_interp = node_par %*% t(Wbin_bicubic[, 1:n_bin_nodes, drop = FALSE]) # spline across bin nodes, evaluated at Bin, for every year node row: n_yr_nodes x length(Bin)
    yr_row = Wyr_bicubic[Year, 1:n_yr_nodes, drop = FALSE] # spline (or constant/blocked) weights across year nodes for this specific Year: 1 x n_yr_nodes
    log_sel = yr_row %*% bin_interp # 1 x length(Bin)

    selex = exp(as.vector(log_sel)) # return spline

  } # end if bicubic spline selectivity

  # 3dgmrf model or 2dar1 (sel devs dimensioned as region, year, bin, sex)
  if(TimeVary_Model %in% c(3:5)) selex = selex * exp(ln_seldevs[Region,Year,,Sex, 1]) # varies semi-parametriclly

  return(selex)
} # end function

