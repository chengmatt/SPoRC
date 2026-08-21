# Stage 2 of 3: objective function
#
# Evaluates selectivity at age or length for every parametric and non parametric
# form the package supports, including the time varying ones.

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
#'     \item{4}{Double-normal dome with plateau and flexible tails (6 parameters).
#'              \eqn{p_{1}} is the bin at which the plateau starts, carried on the bin
#'              scale rather than transformed, \eqn{p_{5}} is the selectivity at the
#'              first bin and \eqn{p_{6}} the selectivity at the last bin.}
#'     \item{5}{Non-parametric selectivity: bin-level logit parameters mapped via \code{plogis},
#'              optionally modified by time-varying deviations.}
#'     \item{6}{Logistic selectivity with asymptote:
#'              \eqn{\alpha / (1 + \exp(-k(\text{bin} - b_{50})))}.
#'              Allows maximum selectivity \eqn{\alpha \in (0,1)}.}
#'     \item{7}{Logistic selectivity with asymptote (b50, b95 parameterization):
#'              \eqn{\alpha / (1 + 19^{(b_{50} - \text{bin})/b_{95}})}.
#'              Equivalent to Model 3 scaled by asymptote \eqn{\alpha}.}
#'     \item{9}{Non-parametric selectivity on the log scale, standardized so
#'              each year's selectivity averages to one across bins. Differs
#'              from model 5 in both respects: 5 bounds every raw value below
#'              one through \code{plogis} and standardizes over years and bins
#'              jointly, whereas 9 leaves the scale free and centers within the
#'              year.}
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
#'     \item{Model 9}{\code{c(ln_sel_1, ..., ln_sel_nbins)}}
#'     \item{Model 8}{Flattened bin-node x year-node log-selectivity grid,
#'       length \code{ncol(Wyr_bicubic) * ncol(Wbin_bicubic)}, filled
#'       column-major into a \code{ncol(Wyr_bicubic)} (rows, year nodes) by
#'       \code{ncol(Wbin_bicubic)} (columns, bin nodes) matrix.}
#'   }
#'
#' @param ln_seldevs Array of log-scale selectivity deviations with dimension
#'   \code{[n_regions, n_years, n_parameters_or_bins, n_sexes, 1]}.
#'
#'   For \code{TimeVary_Model = 1-2}: deviations apply to selectivity parameters
#'   on the natural scale after exponentiation.
#'
#'   For \code{TimeVary_Model = 3-5}: deviations apply multiplicatively at the
#'   bin level to the constructed selectivity curve.
#'
#'   For \code{Selex_Model = 5}: deviations act directly on bin-level logit
#'   selectivity parameters prior to logistic transformation.
#'
#' @param bin_devs Array of log-scale bin-override deviations with dimension
#'   \code{[n_regions, n_years, n_bins, n_sexes, 1]}, or \code{NULL}. Supplies
#'   the value for every bin named in \code{bin_dev_bins}.
#' @param sel_norm_bins Integer vector of bins the mean-one standardization
#'   averages over (\code{Selex_Model = 9}), or \code{NULL} to use every bin.
#'   A gear whose catchability is defined against only part of the age range
#'   standardizes over that part, which shifts the scale absorbed by q.
#' @param bin_dev_bins Integer vector of bins whose selectivity is replaced by
#'   \code{exp(bin_devs[...])} rather than taken from the functional form, or
#'   \code{NULL} for none. The override is applied after everything else,
#'   including any standardization the form performs internally, so the named
#'   bins are governed entirely by their own deviations while the rest of the
#'   curve keeps its parametric shape.
#' @param apical Numeric. For the double normal (\code{Selex_Model == 4}), the
#'   height the ascending and descending limbs are built up to and the plateau
#'   sits at. \code{1} (the default) is the ordinary curve. A sex carrying an
#'   apical offset takes \code{exp(ln_*sel_sex_scale)} here, which moves the
#'   middle of its curve and leaves the selectivity at the first and last bins
#'   where their own parameters put them. Ignored by every other form.
#' @param n_sel_bins Integer or \code{NULL}/\code{0} for none. Bins beyond this
#'   one are held at its computed value rather than evaluated through the
#'   functional form, the \code{NSelBins} plateau convention several existing
#'   assessments apply (e.g. \code{nselages}). Applied after the form and its
#'   parameter deviations, but before the bin-level semi-parametric deviations and
#'   the bin overrides. 
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
#'   bicubic grid), \code{ncol(Wbin_bicubic)}/\code{ncol(Wyr_bicubic)} give
#'   the padded (shared) width, not this block's true node counts, and using
#'   the padded width to reshape \code{pars} would misassign which flattened
#'   parameter values land in which (bin-node, year-node) cell. Default
#'   \code{NULL} falls back to \code{ncol(Wbin_bicubic)}/\code{ncol(Wyr_bicubic)}
#'   which is correct whenever no padding mismatch is possible (a single bicubic
#'   block or fleet, or direct unit testing).
#'
#' @return Numeric vector of selectivity values corresponding to \code{Bin}.
#' Values are on the natural scale and are not normalized unless specified
#' in downstream components.
#'
#' @details
#'
#' For \code{TimeVary_Model = 0}, only the base parametric form is evaluated.
#'
#' For \code{TimeVary_Model = 1-2}, deviations modify model parameters
#' multiplicatively on the natural scale after transformation.
#'
#' For \code{TimeVary_Model = 3-5}, deviations act directly on the constructed
#' selectivity curve as multiplicative log-normal perturbations:
#' \deqn{\text{selex} = \text{selex} \cdot \exp(\delta_{r,y,b,s})}.
#'
#' For \code{Selex_Model = 5}, selectivity is fully non-parametric:
#' bin-specific logit parameters are optionally adjusted by time-varying
#' deviations and transformed via:
#' \deqn{\text{selex}_b = \text{logit}^{-1}(\eta_b)}.
#'
#' For \code{Selex_Model = 9}, selectivity is likewise fully non-parametric but
#' held on the log scale and standardized within each year:
#' \deqn{\text{selex}_b = \exp(\eta_b) / \overline{\exp(\eta)}}.
#' Only the differences among \eqn{\eta} within a year are identified, so the
#' level of \eqn{\eta} is free and is absorbed by catchability or fishing
#' mortality.
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
                     n_yr_nodes_bicubic = NULL,
                     bin_devs = NULL,
                     bin_dev_bins = NULL,
                     sel_norm_bins = NULL,
                     n_sel_bins = NULL,
                     apical = 1) {

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

    # double normal with plateau. The ascending limb is rescaled by its own
    # value at the first bin so that p5 is the selectivity there, and the
    # descending limb by its value at the last bin so that p6 is the
    # selectivity there, with logistic joiners carrying the curve between the
    # limbs and the plateau.
    binwidth <- if(length(Bin) > 1) Bin[2] - Bin[1] else 1

    p1trans <- pars[1] # bin at the start of the plateau, on the bin scale
    p2trans <- p1trans + binwidth + (0.99 * max(Bin) - p1trans - binwidth)/(1 + exp(-1.0 * pars[2])) # bin at the end of the plateau
    p3trans <- exp(pars[3]) # ascending width
    p4trans <- exp(pars[4]) # descending width
    p5trans <- 1/(1 + exp(-1.0 * pars[5])) # selectivity at the first bin
    p6trans <- 1/(1 + exp(-1.0 * pars[6])) # selectivity at the last bin

    if(TimeVary_Model %in% c(1:2)) {
      p1trans = p1trans * exp(ln_seldevs[Region, Year, 1, Sex, 1]) # p1 parameter varying
      p2trans = p2trans * exp(ln_seldevs[Region, Year, 2, Sex, 1]) # p2 parameter varying
      p3trans = p3trans * exp(ln_seldevs[Region, Year, 3, Sex, 1]) # p3 parameter varying
      p4trans = p4trans * exp(ln_seldevs[Region, Year, 4, Sex, 1]) # p4 parameter varying
      p5trans = p5trans * exp(ln_seldevs[Region, Year, 5, Sex, 1]) # p5 parameter varying
      p6trans = p6trans * exp(ln_seldevs[Region, Year, 6, Sex, 1]) # p6 parameter varying
    } # end if iid or random walk

    # construct selectivity function. apical is the height the two limbs are
    # built up to and the plateau sits at, one for the reference sex; the
    # endpoints are anchors that do not move with it
    asc_min <- exp(-((min(Bin) - p1trans)^2/p3trans)) # ascending limb evaluated at the first bin
    dsc_min <- exp(-((max(Bin) - p2trans)^2/p4trans)) # descending limb evaluated at the last bin
    asc <- p5trans + (apical - p5trans) * (exp(-((Bin - p1trans)^2/p3trans)) - asc_min)/(1 - asc_min)
    dsc <- apical + (p6trans - apical) * (exp(-((Bin - p2trans)^2/p4trans)) - 1)/(dsc_min - 1)
    join1 <- 1/(1 + exp(-(20 * (Bin - p1trans)/(1 + abs(Bin - p1trans))))) # joiner between the ascending limb and the plateau
    join2 <- 1/(1 + exp(-(20 * (Bin - p2trans)/(1 + abs(Bin - p2trans))))) # joiner between the plateau and the descending limb
    selex <- asc * (1 - join1) + join1 * (apical * (1 - join2) + dsc * join2) # return parametric form
  }

  if(Selex_Model == 5) {

    if(TimeVary_Model %in% c(1:2)) {
      pars = pars + ln_seldevs[Region, Year, , Sex, 1] # non-parametric parameters varying
    } # end if iid or random walk
    selex = RTMB::plogis(pars) # non parametric parameters

  } # end if non-parametric selectivity, logit space, standardized across all years

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

  if(Selex_Model == 9) {

    # Non-parametric on the log scale, standardized so each year's selectivity
    # averages to one across bins. The logit-scale form (Selex_Model 5) bounds
    # every raw value below one before standardizing and centers over years and
    # bins jointly; this one leaves the scale free and centers within the year. 
    if(TimeVary_Model %in% c(1:2)) pars = pars + ln_seldevs[Region, Year, , Sex, 1]
    selex = exp(pars)
    # Figure out which bins to normalize across
    norm_bins = if(is.null(sel_norm_bins) || length(sel_norm_bins) == 0) seq_along(selex) else sel_norm_bins
    selex = selex / mean(selex[norm_bins])

  } # end if non-parametric on the log scale, standardized within year

  # Hold bins beyond n_sel_bins at the last computed bin's value.
  if(!is.null(n_sel_bins) && n_sel_bins > 0 && n_sel_bins < length(Bin) && Selex_Model != 8) {
    selex[(n_sel_bins + 1):length(Bin)] = selex[n_sel_bins]
  }

  # 3dgmrf model or 2dar1 (sel devs dimensioned as region, year, bin, sex)
  if(TimeVary_Model %in% c(3:5)) selex = selex * exp(ln_seldevs[Region,Year,,Sex, 1]) # varies semi-parametriclly

  # Named bins take a free annual value instead of whatever the form produced.
  # Applied last, so it overrides any within-form standardization 
  if(!is.null(bin_dev_bins) && length(bin_dev_bins) > 0) {
    selex[bin_dev_bins] = exp(bin_devs[Region, Year, bin_dev_bins, Sex, 1])
  }

  return(selex)
} # end function

#' Build the full selectivity array for one selectivity type (retention, fishery, survey)
#'
#' Evaluates \code{\link{Get_Selex}} across regions, years, fleets, and sexes
#' for a single selectivity type, and mean standardizes the result where the
#' time-varying or non-parametric form calls for it. Total fishery, retention,
#' and survey selectivity all share this code path and differ only in which
#' data arrays are handed in.
#'
#' @param selex_type Integer switch (0 = age-based, 1 = length-based).
#' @param bins Numeric vector of bins; \code{ages} when \code{selex_type == 0}
#'   and \code{lens} when \code{selex_type == 1}.
#' @param sel_blocks Integer array \code{n_regions x n_yrs x n_fleets} of
#'   selectivity block indices.
#' @param sel_model Integer array \code{n_regions x n_yrs x n_fleets} of
#'   selectivity functional forms (see \code{\link{Get_Selex}}).
#' @param fixed_sel_pars Array \code{n_regions x n_pars x n_blocks x n_sexes x
#'   n_fleets} of fixed-effect selectivity parameters.
#' @param cont_tv_sel Integer matrix \code{n_regions x n_fleets} of continuous
#'   time-varying forms.
#' @param ln_seldevs Array \code{n_regions x n_yrs_total x n_bins x n_sexes x
#'   n_fleets} of selectivity deviations; sliced per fleet internally.
#' @param use_fixed_sel Integer vector of length \code{n_fleets} (0 = estimate,
#'   1 = use \code{sel_input}).
#' @param sel_input Array \code{n_pop x n_regions x n_yrs x n_seas x n_bins x
#'   n_sexes x n_fleets} of fixed selectivity values.
#' @param bicubic_Wbin,bicubic_Wyr Bicubic spline bin-node and year-node weight
#'   arrays (\code{Selex_Model == 8} only).
#' @param bicubic_binnodes,bicubic_yrnodes Integer arrays \code{n_regions x
#'   n_yrs x n_fleets} giving each block's true node counts
#'   (\code{Selex_Model == 8} only).
#' @param n_pop,n_regions,n_yrs,n_proj_yrs_devs,n_seas,n_ages,n_lens,n_sexes,n_fleets
#'   Model dimensions.
#' @param sex_par_offset Integer vector of length \code{n_fleets} (0/1), or
#'   \code{NULL} for all zero. For a fleet flagged 1, the stored fixed-effect
#'   selectivity parameters of every sex beyond the first are additive offsets
#'   on the first sex's stored parameters, evaluated as
#'   \code{pars[sex 1] + pars[sex s]} on the transformed (usually log) scale.
#'   Because most forms exponentiate their parameters, an offset of \eqn{\delta}
#'   makes the sex-\eqn{s} natural-scale parameter the first sex's value times
#'   \eqn{e^{\delta}}, which is the male-offset convention several existing
#'   assessments use. An offset held at zero reproduces sex-shared parameters.
#' @param sex_scale_offset Integer vector of length \code{n_fleets} (0/1), or
#'   \code{NULL} for all zero. For a fleet flagged 1, the realized selectivity
#'   curve of every sex beyond the first is multiplied by
#'   \code{exp(sex_scale[region, block, sex, fleet])}, a constant log-scale
#'   offset on the whole curve. The scaled curve may exceed one, matching the
#'   convention of a male selectivity offset applied in log space. Not
#'   meaningful for forms whose scale is standardized away downstream
#'   (non-parametric forms and the semi-parametric time-varying structures);
#'   setup refuses those combinations.
#' @param sex_apical_offset Integer vector of length \code{n_fleets} (0/1), or
#'   \code{NULL}. Where it is 1, every sex beyond the first builds its double
#'   normal up to \code{exp(sex_scale[region, block, sex, fleet])} rather than
#'   to one, leaving the first and last bins where their own parameters put
#'   them. Mutually exclusive with a scale offset on the same fleet.
#' @param sex_scale Array \code{n_regions x n_blocks x n_sexes x n_fleets} of
#'   log-scale curve offsets read when \code{sex_scale_offset} or
#'   \code{sex_apical_offset} flags a fleet, or \code{NULL} when no fleet is
#'   flagged.
#' @param nselbins Integer array \code{n_regions x n_yrs x n_fleets} naming
#'   each cell's plateau bin (\code{0} = none), or \code{NULL}. Bins beyond it
#'   are held at its value for the parametric forms (see \code{n_sel_bins} in
#'   \code{\link{Get_Selex}}); the bicubic form's own setup-built plateau is
#'   unaffected.
#'
#' @return List with \code{sel}, the age-based array (\code{n_pop x n_regions x
#'   n_yrs_total x n_seas x n_ages x n_sexes x n_fleets}), and \code{sel_l}, the
#'   length-based array (\code{n_regions x n_yrs_total x n_lens x n_sexes x
#'   n_fleets}). Only the array matching \code{selex_type} is populated; when
#'   \code{selex_type == 1}, \code{sel} is filled downstream by multiplying
#'   \code{sel_l} through the size-age transition matrix.
#'
#' @keywords internal
Get_Selex_Array = function(selex_type,
                           bins,
                           sel_blocks,
                           sel_model,
                           fixed_sel_pars,
                           cont_tv_sel,
                           ln_seldevs,
                           use_fixed_sel,
                           sel_input,
                           bicubic_Wbin,
                           bicubic_Wyr,
                           bicubic_binnodes,
                           bicubic_yrnodes,
                           n_pop,
                           n_regions,
                           n_yrs,
                           n_proj_yrs_devs,
                           n_seas,
                           n_ages,
                           n_lens,
                           n_sexes,
                           n_fleets,
                           bin_devs = NULL,
                           bin_dev_bins = NULL,
                           sel_norm_bins = NULL,
                           sex_par_offset = NULL,
                           sex_scale_offset = NULL,
                           sex_apical_offset = NULL,
                           sex_scale = NULL,
                           nselbins = NULL) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  sel = array(data = 0, dim = c(n_pop, n_regions, n_yrs + n_proj_yrs_devs, n_seas, n_ages, n_sexes, n_fleets)) # age-based selectivity
  sel_l = array(data = 0, dim = c(n_regions, n_yrs + n_proj_yrs_devs, n_lens, n_sexes, n_fleets)) # length-based selectivity

  for(r in 1:n_regions) {
    for(y in 1:(n_yrs + n_proj_yrs_devs)) {
      for(f in 1:n_fleets) {

        y_idx = min(y, n_yrs) # projection years reuse the terminal year's block structure

        # estimating selectivity
        if(use_fixed_sel[f] == 0) {
          for(s in 1:n_sexes) {

            # Extract variables
            sel_blk_idx = sel_blocks[r,y_idx,f] # selectivity block index
            tmp_sel_model = sel_model[r,y_idx,f] # selectivity model
            tmp_n_bin_nodes = bicubic_binnodes[r,y_idx,f] # true bin-node count for this block (Selex_Model == 8 only)
            tmp_n_yr_nodes = bicubic_yrnodes[r,y_idx,f] # true year-node count for this block (Selex_Model == 8 only)

            # Extract out fixed-effect selectivity parameters for a given block. For a sex parameter offset, sexes beyond the first store additive offsets on the first sex's transformed-scale parameters
            tmp_sel_vec = fixed_sel_pars[r,,sel_blk_idx,s,f]
            if(!is.null(sex_par_offset) && s > 1 && sex_par_offset[f] == 1) tmp_sel_vec = fixed_sel_pars[r,,sel_blk_idx,1,f] + tmp_sel_vec

            # Extract bicubic spline interpolation weight matrices for this block
            tmp_Wbin = array(bicubic_Wbin[r,,,sel_blk_idx,f], dim = dim(bicubic_Wbin)[c(2,3)])
            tmp_Wyr = array(bicubic_Wyr[r,,,sel_blk_idx,f], dim = dim(bicubic_Wyr)[c(2,3)])

            # Compute selectivity functional form
            tmp_sel = Get_Selex(Selex_Model = tmp_sel_model, # selectivity model
                                TimeVary_Model = cont_tv_sel[r,f], # time varying model
                                pars = tmp_sel_vec, # fixed effect selectivity parameters
                                ln_seldevs = ln_seldevs[,,,,f, drop = FALSE], # Selectivity deviations
                                Region = r, # region index
                                Year = y, # year index
                                Bin = bins, # bin vector
                                Sex = s, # sex index
                                Wbin_bicubic = tmp_Wbin, # bicubic spline bin-node weight matrix (Selex_Model == 8 only)
                                Wyr_bicubic = tmp_Wyr, # bicubic spline year-node weight matrix (Selex_Model == 8 only)
                                n_bin_nodes_bicubic = tmp_n_bin_nodes, # true bin-node count for this block (Selex_Model == 8 only)
                                n_yr_nodes_bicubic = tmp_n_yr_nodes, # true year-node count for this block (Selex_Model == 8 only)
                                bin_devs = bin_devs, # bin-override deviations
                                bin_dev_bins = if(is.null(bin_dev_bins)) NULL else which(bin_dev_bins[,f] == 1), # bins this fleet overrides
                                sel_norm_bins = if(is.null(sel_norm_bins)) NULL else which(sel_norm_bins[,f] == 1), # bins this fleet standardizes over
                                n_sel_bins = if(is.null(nselbins)) NULL else nselbins[r,y_idx,f], # plateau bin for parametric forms (0 = none)
                                # a sex carrying an apical offset builds its limbs up to its own
                                # height rather than to one, which leaves the first and last bins
                                # where their own parameters put them
                                apical = if(!is.null(sex_apical_offset) && s > 1 && sex_apical_offset[f] == 1) exp(sex_scale[r,sel_blk_idx,s,f]) else 1
            )

            # For a sex scale offset, sexes beyond the first carry a constant log-scale offset on the whole realized curve
            if(!is.null(sex_scale_offset) && s > 1 && sex_scale_offset[f] == 1) tmp_sel = tmp_sel * exp(sex_scale[r,sel_blk_idx,s,f])

            # Compute selectivity
            if(selex_type == 0) {
              for(p in 1:n_pop) {
                for(seas in 1:n_seas) {
                  sel[p,r,y,seas,,s,f] = tmp_sel # age-based selectivity
                } # end seas loop
              } # end p loop
            }
            if(selex_type == 1) sel_l[r,y,,s,f] = tmp_sel # input into length-based selectivity

          } # end s loop
        } else {

          # Input fixed selectivity; sel_input only spans n_yrs, so projection years reuse the terminal year
          if(selex_type == 0) {
            for(p in 1:n_pop) {
              for(seas in 1:n_seas) {
                sel[p,r,y,seas,,,f] = sel_input[p,r,y_idx,seas,,,f] # age-based selectivity
              } # end seas loop
            } # end p loop
          }
          # length-based selectivity carries no population or season dimension, so read the first of each
          if(selex_type == 1) sel_l[r,y,,,f] = sel_input[1,r,y_idx,1,,,f]

        } # end if else for whether or not using fixed selex inputs

      } # end f loop
    } # end y loop
  } # end r loop

  # Mean standardizing to help with interpretability
  for(r in 1:n_regions) {
    for(f in 1:n_fleets) {
      # Selex_Model 9 already standardizes within each year, so it is excluded
      # here rather than being centered a second time across years and bins.
      if((cont_tv_sel[r,f] %in% 3:5 || any(sel_model[r,,f] == 5)) && !any(sel_model[r,,f] == 9)) {
        for(s in 1:n_sexes) {

          if(selex_type == 0) {
            tmp_mean = log(mean(sel[1,r,,1,,s,f])) # indexing 1 for pop and season because those dims not estimated (depends on growth transition)
            for(p in 1:n_pop) for(seas in 1:n_seas)
              sel[p,r,,seas,,s,f] = exp(log(sel[p,r,,seas,,s,f]) - tmp_mean)
          }
          if(selex_type == 1) sel_l[r,,,s,f] = exp(log(sel_l[r,,,s,f]) - log(mean(sel_l[r,,,s,f]))) # length-based selectivity

        } # end s loop
      } # end if mean standardizing
    } # end f loop
  } # end r loop

  return(list(sel = sel, sel_l = sel_l))
} # end function

