# Stage 2 of 3: objective function
#
# Every prior and penalty in the objective, one small function each. These return positive log
# likelihoods and the objective negates them.

# Selectivity Smoothness ----------------------------------------------------

#' Compute a model-agnostic selectivity smoothness / dome-shape penalty (Positive Scale)
#'
#' Regularization penalty operating directly on a realized selectivity-at-bin-at-year
#' surface, rather than on any particular parameterization's deviations. Because it
#' only ever looks at the resulting selectivity values, it applies uniformly to any
#' selectivity functional form and any fleet, called once per fleet from the
#' "Selectivity Smoothness Penalty" section of \code{SPoRC_rtmb.R}.
#'
#' @param sel_vals Array of selectivity values dimensioned \code{[1, year, bin, sex, 1]}.
#'   Evaluated on the log scale internally.
#' @param wt_bin_curve Non-negative weight on the age/bin curvature penalty: the
#'   sum of squared second differences of log-selectivity across bins, within each
#'   year, normalized by the number of bins. Penalizes jagged (non-smooth)
#'   selectivity-at-age curves. \code{0} (default) disables this term. Requires at
#'   least 3 bins to have any effect.
#' @param wt_bin_diff Non-negative weight on the unconditional bin first-difference
#'   penalty: the sum of squared first differences of log-selectivity across bins,
#'   within each year, normalized by the number of bins. Unlike \code{wt_dome}
#'   (which only penalizes decreases), both increases and decreases contribute. Requires at least 2 bins to have any
#'   effect.
#' @param wt_yr_diff Non-negative weight on the inter-annual first-difference
#'   penalty: the sum of squared first differences of log-selectivity across years,
#'   within each bin, normalized by the number of years. Penalizes abrupt year-to-year jumps in
#'   selectivity-at-bin. \code{0} (default) disables this term. Requires at least 2
#'   years to have any effect.
#' @param wt_yr_curve Non-negative weight on the inter-annual second-difference
#'   (smoothness) penalty: the sum of squared second differences of log-selectivity
#'   across years, within each bin, normalized by the number of years. Penalizes jagged (non-smooth) year-to-year
#'   selectivity trajectories. \code{0} (default) disables this term. Requires at
#'   least 3 years to have any effect.
#' @param wt_dome Non-negative weight on the dome-shape (non-monotonicity) penalty:
#'   for each year, penalizes any decrease in log-selectivity moving from one bin to
#'   the next (i.e. discourages, but does not forbid, dome shaped dynamics. \code{0} (default) disables this term.
#' @param wt_mean_center Non-negative weight on a per-year mean-centering
#'   (sum-to-zero) regularization: for each year, penalizes the squared mean of
#'   log-selectivity across bins. \code{0} (default) disables this term;
#'   set to \code{10000}.
#' @param normalize Logical. If \code{TRUE} (default), \code{wt_bin_curve} is
#'   divided by the number of bins the penalties act over and
#'   \code{wt_yr_diff}/\code{wt_yr_curve} are divided by the number of years.
#'   \code{SPoRC_rtmb.R} always calls this with \code{normalize = TRUE}.
#' @param bin_range Length-two vector giving the first and last bin the
#'   penalties act over, or \code{NULL} (default) for every bin. Restricting the
#'   range is how a shape penalty is confined to the older ages where a curve is
#'   expected to flatten, without constraining the ascending limb.
#'
#' @details Every \code{wt_} argument accepts either a single number applied to
#'   all years, or a vector with one value per year. A per-year vector lets a
#'   penalty act only in the years where selectivity is allowed to change, or
#'   act with a different strength in each year, which is how a random walk with
#'   a year-specific standard deviation is expressed: set the year's weight to
#'   \code{1 / (2 * sigma^2)} and pass \code{normalize = FALSE}. Years whose
#'   weight is zero are skipped entirely.
#'
#' @return Numeric scalar: the positive log-likelihood contribution from the
#'   requested penalty terms. Negated externally to form the negative log-likelihood.
#'
#' @keywords internal
#' @import RTMB
Get_Selex_Smoothness_Penalty <- function(
  sel_vals,
  wt_bin_curve = 0,
  wt_bin_diff = 0,
  wt_yr_diff = 0,
  wt_yr_curve = 0,
  wt_dome = 0,
  wt_mean_center = 0,
  normalize = TRUE,
  bin_range = NULL,
  yr_diff_ref = NULL
) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  loglik = 0 # initialize likelihood (positive scale, negated by the caller)

  n_yrs = dim(sel_vals)[2]
  n_bins = dim(sel_vals)[3]
  n_sexes = dim(sel_vals)[4]

  # A weight is either one number for the whole series or one per year, so that a
  # penalty can act only in some years or act with a different strength in each.
  expand_wt = function(w) if(length(w) == 1) rep(w, n_yrs) else w
  wt_bin_curve = expand_wt(wt_bin_curve)
  wt_bin_diff = expand_wt(wt_bin_diff)
  wt_yr_diff = expand_wt(wt_yr_diff)
  wt_yr_curve = expand_wt(wt_yr_curve)
  wt_dome = expand_wt(wt_dome)
  wt_mean_center = expand_wt(wt_mean_center)

  # bin_range and normalize are either one setting shared by every term or a named list giving each
  # its own, so a shape penalty on the older ages and a walk over every age can be set together
  get_bins = function(term) {
    br = if(is.list(bin_range)) bin_range[[term]] else bin_range
    if(is.null(br)) return(c(1, n_bins))
    c(max(1, br[1]), min(n_bins, br[2]))
  }
  get_norm = function(term) {
    nz = if(is.list(normalize)) normalize[[term]] else normalize
    if(is.null(nz)) TRUE else nz
  }

  bins = get_bins("smooth_bin_curve"); b_lo = bins[1]; b_hi = bins[2]
  bin_norm = if(get_norm("smooth_bin_curve")) b_hi - b_lo + 1 else 1
  if(any(wt_bin_curve != 0) && (b_hi - b_lo + 1) >= 3) { # age/bin curvature (second difference across bins)
    for(s in 1:n_sexes) {
      for(y in 1:n_yrs) {
        if(wt_bin_curve[y] == 0) next
        for(b in (b_lo + 1):(b_hi - 1)) {
          bin_penalty = log(sel_vals[1,y,b+1,s,1]) - 2 * log(sel_vals[1,y,b,s,1]) + log(sel_vals[1,y,b-1,s,1])
          loglik = loglik - wt_bin_curve[y] / bin_norm * bin_penalty^2
        } # end b loop
      } # end y loop
    } # end s loop
  }

  bins = get_bins("smooth_bin_diff"); b_lo = bins[1]; b_hi = bins[2]
  bin_norm = if(get_norm("smooth_bin_diff")) b_hi - b_lo + 1 else 1
  if(any(wt_bin_diff != 0) && (b_hi - b_lo + 1) >= 2) { # unconditional bin first-difference (both directions penalized, unlike wt_dome)
    for(s in 1:n_sexes) {
      for(y in 1:n_yrs) {
        if(wt_bin_diff[y] == 0) next
        for(b in b_lo:(b_hi - 1)) {
          bin_diff_penalty = log(sel_vals[1,y,b,s,1]) - log(sel_vals[1,y,b+1,s,1])
          loglik = loglik - wt_bin_diff[y] / bin_norm * bin_diff_penalty^2
        } # end b loop
      } # end y loop
    } # end s loop
  }

  bins = get_bins("smooth_yr_diff"); b_lo = bins[1]; b_hi = bins[2]
  yr_norm = if(get_norm("smooth_yr_diff")) n_yrs else 1

  # the walk has no previous value in its first year, so that year is normally unpenalized. a
  # reference holds the first penalized year toward yr_diff_ref on the log scale instead
  yr_ref_first = if(is.null(yr_diff_ref)) 0 else which(wt_yr_diff != 0)[1]
  yr_ref = if(is.null(yr_diff_ref)) NULL else rep(yr_diff_ref, length.out = n_bins)
  if(any(wt_yr_diff != 0) && n_yrs >= 2) { # inter-annual first difference
    for(s in 1:n_sexes) {
      for(b in b_lo:b_hi) {
        for(y in 1:n_yrs) {
          if(wt_yr_diff[y] == 0) next
          if(y == 1 && y != yr_ref_first) next
          prev_sel = if(y == yr_ref_first) yr_ref[b] else log(sel_vals[1,y-1,b,s,1])
          yr_diff_penalty = log(sel_vals[1,y,b,s,1]) - prev_sel
          loglik = loglik - wt_yr_diff[y] / yr_norm * yr_diff_penalty^2
        } # end y loop
      } # end b loop
    } # end s loop
  }

  bins = get_bins("smooth_yr_curve"); b_lo = bins[1]; b_hi = bins[2]
  yr_norm = if(get_norm("smooth_yr_curve")) n_yrs else 1
  if(any(wt_yr_curve != 0) && n_yrs >= 3) { # inter-annual second difference / smoothness
    for(s in 1:n_sexes) {
      for(b in b_lo:b_hi) {
        for(y in 2:(n_yrs - 1)) {
          if(wt_yr_curve[y] == 0) next
          year_penalty = log(sel_vals[1,y+1,b,s,1]) - 2 * log(sel_vals[1,y,b,s,1]) + log(sel_vals[1,y-1,b,s,1])
          loglik = loglik - wt_yr_curve[y] / yr_norm * year_penalty^2
        } # end y loop
      } # end b loop
    } # end s loop
  }

  bins = get_bins("smooth_dome"); b_lo = bins[1]; b_hi = bins[2]
  if(any(wt_dome != 0) && (b_hi - b_lo + 1) >= 2) { # dome-shape / non-monotonicity, within each year
    for(s in 1:n_sexes) {
      for(y in 1:n_yrs) {
        if(wt_dome[y] == 0) next
        for(b in b_lo:(b_hi - 1)) {
          decrease = max(log(sel_vals[1,y,b,s,1]) - log(sel_vals[1,y,b+1,s,1]), 0) # only decreases contribute
          loglik = loglik - wt_dome[y] * decrease^2
        } # end b loop
      } # end y loop
    } # end s loop
  }

  if(any(wt_mean_center != 0)) { # per-year mean-centering / sum-to-zero regularization
    bins = get_bins("smooth_mean_center")
    for(s in 1:n_sexes) {
      for(y in 1:n_yrs) {
        if(wt_mean_center[y] == 0) next
        z = mean(log(sel_vals[1,y,bins[1]:bins[2],s,1]))
        loglik = loglik - wt_mean_center[y] * z^2
      } # end y loop
    } # end s loop
  }

  return(loglik)
} # return log likelihood

# Process Error Penalties ---------------------------------------------------

#' Compute Process Error Log-Likelihood for a Deviation Surface (Positive Scale)
#'
#' Calculates the positive log-likelihood contribution for a surface of
#' deviations indexed by year and by some second dimension, under a variety of
#' temporal and spatiotemporal structures. Selectivity deviations use it over
#' years and bins, growth's semi-parametric deviations over years and ages, and
#' a time-varying growth parameter over years alone (a surface one column wide).
#' The argument names still read \code{bin} for that second dimension.
#'
#' The function supports:
#' \itemize{
#'   \item IID process error
#'   \item Random walk process error
#'   \item 3D Gaussian Markov Random Field (GMRF) models (marginal or conditional variance)
#'   \item Separable 2D AR(1) models
#' }
#'
#'
#' \strong{Note:} The returned value is on the \emph{positive} log-likelihood scale.
#' It must be negated to obtain a negative log-likelihood contribution, which is
#' handled outside this function.
#'
#' @param PE_model Integer specifying the process error structure:
#' \itemize{
#'   \item 1 = IID: deviations drawn independently as \eqn{N(0, \sigma^2)}.
#'   \item 2 = Random walk: deviations follow a first-order random walk initialized
#'     with a diffuse prior (\eqn{\sigma = 5}) at \code{y = 1}.
#'   \item 3 = 3D GMRF with marginal variance parameterization.
#'   \item 4 = 3D GMRF with conditional variance parameterization.
#'   \item 5 = Separable 2D AR(1) across bins and years.
#' }
#'
#' @param PE_pars Array of process error parameters dimensioned
#'   \code{[1, par_index, sex, 1]}. The \code{par_index} slot meaning
#'   depends on \code{PE_model}:
#' \itemize{
#'   \item Models 1-2: \code{[1,1,s,1]} = log standard deviation (\eqn{\log \sigma})
#'     for sex \code{s}, indexed by bin/age.
#'   \item Models 3-4: \code{[1,1,s,1]} = unconstrained partial correlation by age/bin;
#'     \code{[1,2,s,1]} = unconstrained partial correlation by year;
#'     \code{[1,3,s,1]} = unconstrained partial correlation by cohort;
#'     \code{[1,4,s,1]} = log variance.
#'   \item Model 5: \code{[1,1,s,1]} = unconstrained bin correlation (transformed via
#'     \eqn{2/(1+e^{-2x})-1}); \code{[1,2,s,1]} = unconstrained year correlation;
#'     \code{[1,4,s,1]} = log standard deviation.
#' }
#'
#' @param ln_devs Array of log-scale selectivity deviations dimensioned
#'   \code{[1, year, bin, sex, 1]}.
#'
#' @param map_sel_devs Integer array dimensioned \code{[fleet, year, bin, sex]}
#'   mapping deviations to unique estimated parameters. Shared deviations
#'   hold the same integer value; \code{NA} entries are treated as fixed
#'   and excluded from likelihood evaluation.
#'
#' @param map_sel_devs_full The same map across every unit this penalty is
#'   evaluated over, that unit being the first dim: regions for selectivity,
#'   populations by region for growth. A deviation shared over those units is
#'   one parameter appearing in each of their slices, and this function runs one
#'   unit at a time, so its contribution is divided by the number holding it.
#'   Without the split, a series shared over \code{n} units is penalized
#'   \code{n} times, an implicit \eqn{\sigma / \sqrt{n}}. A deviation that is
#'   not shared appears once, divides by one, and is unaffected.
#' @param rw_init_sigma Standard deviation given to the first year of a random
#'   walk. A number (5 by default) leaves that year effectively unconstrained;
#'   \code{NA} starts the walk at zero under its own sigma instead.
#' @param min_sel_devs_shared_bins Integer vector. Indices of the reference (minimum) bin
#'   within each shared deviation group, used to subset the bin dimension when
#'   evaluating GMRF or 2D AR(1) likelihoods (PE models 3-5). When no bin sharing
#'   is specified, defaults to \code{1:n_bins} (i.e., all bins are included).
#'
#' @return Numeric scalar: the positive log-likelihood contribution from selectivity
#'   process error. Negated externally to form the negative log-likelihood.
#'
#' @keywords internal
#' @import RTMB
Get_PE_loglik <- function(PE_model,
                              PE_pars,
                              ln_devs,
                              map_sel_devs,
                              map_sel_devs_full,
                              min_sel_devs_shared_bins,
                              rw_init_sigma = 5
                              ) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  # Note that the likelihood calculations are positive within the function,
  # because it gets converted to negative outside the wrapper function

  loglik = 0 # initialize likelihood

  # find unique selectivity deviations to penalize (sort drops NAs)
  unique_sel_devs = sort(unique(as.vector(map_sel_devs)))

  # Exit out fxn if this region x fleet slice has no estimated deviations at all.
  if(length(unique_sel_devs) == 0) return(loglik)
  n_yrs = dim(map_sel_devs)[2] # get years for indexing
  n_bins = dim(map_sel_devs)[3] # get bins / pars for indexing
  n_sexes = dim(map_sel_devs)[4] # get sexes for indexing

  # a shared deviation is one parameter appearing in every sharing unit's slice, and this runs one
  # unit at a time, so its penalty is split over the units holding it. dim 1 is the unit
  unit_of = slice.index(map_sel_devs_full, 1) # unit index of every cell
  is_est = !is.na(map_sel_devs_full)
  n_units_sharing = tapply(unit_of[is_est], map_sel_devs_full[is_est],
                           function(x) length(unique(x)))

  if(PE_model %in% c(1, 2)) {

    for(dev_idx in 1:length(unique_sel_devs)) {

      # figure out where unique sel devs first occur
      idx = which(map_sel_devs == unique_sel_devs[dev_idx], arr.ind = TRUE)[1,]
      y = idx[2] # get unique year index
      i = idx[3] # get unique age or parmeter index
      s = idx[4] # get unique sex index
      share = as.numeric(n_units_sharing[as.character(unique_sel_devs[dev_idx])]) # units holding it

      if(PE_model == 1) {
        if(y >= 1) loglik = loglik + RTMB::dnorm(ln_devs[1,y,i,s,1], 0, exp(PE_pars[1,i,s,1]), TRUE) / share
      } # iid process error

      if(PE_model == 2) {
        # the walk needs a distribution for its first year. a wide sigma leaves it unconstrained;
        # NA starts the walk at zero under its own sigma
        if(y == 1) {
          init_sd = if(is.na(rw_init_sigma)) exp(PE_pars[1,i,s,1]) else rw_init_sigma
          loglik = loglik + RTMB::dnorm(ln_devs[1,y,i,s,1], 0, init_sd, TRUE) / share
        }
        else loglik = loglik + RTMB::dnorm(ln_devs[1,y,i,s,1], ln_devs[1,y-1,i,s,1], exp(PE_pars[1,i,s,1]), TRUE) / share
      } # end random walk process error

    } # end dev_idx loop

  } # end iid or random walk process error

  if(PE_model %in% c(3,4,5)) {

    if(PE_model == 3) Var_Type = 0 # marginal variance
    if(PE_model == 4) Var_Type = 1 # conditional variance

    # get first unique combination
    unique_comb = which(map_sel_devs == unique_sel_devs[1], arr.ind = TRUE)[1,]
    # cbind to get all unique combinations here (cbinding first one, so loop starts at 2)
    for(i in 2:length(unique_sel_devs)) unique_comb = cbind(unique_comb, which(map_sel_devs == unique_sel_devs[i], arr.ind = TRUE)[1,])

    # Next, get unique sex deviations
    unique_s = unique(unique_comb[4,])

    for(idx in 1:length(unique_s)) {

      s = unique_s[idx] # get sex index

      # the whole density is evaluated once per unit holding these deviations, so split it as the
      # iid and walk forms do. sharing is set per fleet, so the levels here all share alike
      s_levels = unique_sel_devs[unique_comb[4,] == s] # levels this sex evaluates
      share = max(as.numeric(n_units_sharing[as.character(s_levels)])) # units holding them

      # Construct precision matrix for 3d gmrf
      if(PE_model %in% c(3,4)) {

        # the precision matrix spans the bins the deviations are actually evaluated
        # over, which is one per shared group, not the full bin dimension
        Q = Get_3d_precision(n_ages = length(min_sel_devs_shared_bins), # number of ages
                             n_yrs = n_yrs,  # number of years
                             pcorr_age = PE_pars[1,1,s,1], # unconstrained partial correlation by age
                             pcorr_year = PE_pars[1,2,s,1], # unconstrained partial correlation by year
                             pcorr_cohort = PE_pars[1,3,s,1], # unconstrained partial correlation by cohort
                             ln_var_value = PE_pars[1,4,s,1], # log variance
                             Var_Type = Var_Type
                             ) # variance type, == 0 (marginal), == 1 (conditional)

        # apply gmrf likelihood
        eps_ay = as.vector(t(ln_devs[1,,min_sel_devs_shared_bins,s,1])) # convert to vector
        loglik = loglik + RTMB::dgmrf(x = eps_ay, mu = 0, Q = Q, log = TRUE) / share
      } # end if

      # 2dar1 model
      if(PE_model == 5) {
        # Function to constrain values between -1 and 1
        rho_trans = function(x) 2/(1+ exp(-2 * x)) - 1

        # Extract out varaibles and transform into appropriate space
        eps_ya = ln_devs[1,,min_sel_devs_shared_bins,s,1] # needs to be in matrix format for dseparable
        rho_b = rho_trans(PE_pars[1,1,s,1]) # correlation across bins
        rho_y = rho_trans(PE_pars[1,2,s,1]) # correlation across years
        sigma2 = exp(PE_pars[1,4,s,1])^2 # get sigma

        # Define 2d scale
        scale = sqrt(sigma2) / sqrt(1 - rho_y^2) / sqrt(1 - rho_b^2)

        # Define ar1 separable functions
        f1 = function(x) RTMB::dautoreg(x, mu = 0, phi = rho_y, log = TRUE)
        f2 = function(x) RTMB::dautoreg(x, mu = 0, phi = rho_b, log = TRUE)
        loglik = loglik + RTMB::dseparable(f1, f2)(eps_ya, scale = scale) / share
      } # end if
    } # end idx loop

  } # end 3dgrmf or 2dar1 process error

  return(loglik)
} # return log likelihood

#' Compute Movement Process Error Log-Likelihood (Positive Scale)
#'
#' Calculates the positive log-likelihood contribution for movement
#' process error deviations under multiple IID structural assumptions.
#' Deviations are penalized as \eqn{N(0, \sigma^2)} where \eqn{\sigma}
#' is drawn from \code{PE_pars} according to the selected model structure.
#' Only origin-destination pairs that are adjacent (non-zero in
#' \code{adjacency_collapsed}) contribute to the likelihood.
#'
#' \strong{Note:} The returned value is on the \emph{positive} log-likelihood
#' scale. It must be negated externally to form the negative log-likelihood.
#'
#' @param PE_model Integer specifying the movement process error structure.
#'   All models are IID; they differ in which dimensions share a common
#'   standard deviation. Models 1-5 are single-population (fix \code{pop = 1});
#'   models 6-10 estimate separate parameters per population:
#' \itemize{
#'   \item \strong{1}: IID across years (single \eqn{\sigma} per origin region)
#'   \item \strong{2}: IID across ages (single \eqn{\sigma} per origin region and age)
#'   \item \strong{3}: IID across years and ages
#'   \item \strong{4}: IID across years, ages, and sexes
#'   \item \strong{5}: IID across years, seasons, ages, and sexes
#'   \item \strong{6}: IID across populations and years
#'   \item \strong{7}: IID across populations and ages
#'   \item \strong{8}: IID across populations, years, and ages
#'   \item \strong{9}: IID across populations, years, ages, and sexes
#'   \item \strong{10}: IID across populations, years, seasons, ages, and sexes
#' }
#'
#' @param PE_pars Array of movement process error parameters (log standard
#'   deviations) dimensioned \code{[pop, from_region, seas, age, sex]}.
#'   Exponentiated internally to obtain \eqn{\sigma}. Which dimensions
#'   are active depends on \code{PE_model}; unused dimensions should be
#'   fixed at a constant (e.g., index 1) via the parameter map.
#'
#' @param move_devs Movement deviation array dimensioned
#'   \code{[pop, from_region, to_region, year, seas, age, sex]}.
#'
#' @param map_move_devs Integer array dimensioned
#'   \code{[pop, from_region, to_region, year, seas, age, sex]}
#'   mapping deviations to unique estimated parameters. Shared deviations
#'   hold the same integer value; dimensions are extracted from this
#'   array to determine loop bounds.
#'
#' @param do_recruits_move Integer (0/1). If \code{0} and \code{n_ages >= 2},
#'   age-1 recruits are excluded from the likelihood (loop starts at age 2).
#'   If \code{1}, all ages including recruits are penalized.
#'
#' @param adjacency_collapsed Square \code{[n_regions x n_regions]} matrix
#'   of allowable movement connections among regions, excluding self-retention
#'   (diagonal entries should be 0). Origin-destination pairs with a value
#'   of 0 are skipped and contribute nothing to the likelihood.
#'
#' @param move_type Integer specifying the movement formulation:
#' \itemize{
#'   \item \strong{0} = Unstructured multinomial logit movement
#'   \item \strong{1} = CTMC-based movement
#' }
#'   Currently used for dispatch context; likelihood computation is
#'   identical across movement types within this function.
#'
#' @return Numeric scalar: the positive log-likelihood contribution from
#'   movement process error deviations. Negated externally to form the
#'   negative log-likelihood.
#'
#' @keywords internal
#' @import RTMB
Get_move_PE_loglik <- function(PE_model,
                               PE_pars,
                               move_devs,
                               map_move_devs,
                               do_recruits_move,
                               adjacency_collapsed,
                               move_type
                               ) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  # Note that the likelihood calculations are positive within the function,
  # because it gets converted to negative outside the wrapper function

  loglik = 0 # initialize likelihood

  # Get dimensions for penalty
  n_pop = dim(map_move_devs)[1]
  n_regions_from = dim(map_move_devs)[2]
  n_regions_to = dim(map_move_devs)[3]
  n_yrs = dim(map_move_devs)[4]
  n_seas = dim(map_move_devs)[5]
  n_ages = dim(map_move_devs)[6]
  n_sexes = dim(map_move_devs)[7]

  # whether recruits move
  age_start = ifelse(do_recruits_move == 0 && n_ages >= 2, 2, 1)

  # dims named in each PE_model's spec (see cont_move_map in Setup_Movement.R for the type to
  # integer correspondence); dims absent from a model's name are shared across, kept at index 1
  key_dims_by_model <- list(
    `1` = "year", `2` = "age", `3` = c("year","age"), `4` = c("year","age","sex"),
    `5` = c("year","season","age","sex"),
    `6` = c("pop","year"), `7` = c("pop","age"), `8` = c("pop","year","age"),
    `9` = c("pop","year","age","sex"), `10` = c("pop","year","season","age","sex")
  )
  key_dims <- key_dims_by_model[[as.character(PE_model)]]

  pop_idx  = if("pop"    %in% key_dims) 1:n_pop   else 1
  yr_idx   = if("year"   %in% key_dims) 1:n_yrs   else 1
  seas_idx = if("season" %in% key_dims) 1:n_seas  else 1
  age_idx  = if("age"    %in% key_dims) age_start:n_ages else 1
  sex_idx  = if("sex"    %in% key_dims) 1:n_sexes else 1

  # Penalize Deviations
  for(rr in 1:n_regions_to) {
    for(r in 1:n_regions_from) {

      if(adjacency_collapsed[r,rr] == 0) next # skip

      for(p in pop_idx) {
        for(y in yr_idx) {
          for(seas in seas_idx) {
            for(a in age_idx) {
              for(s in sex_idx) {
                loglik = loglik + RTMB::dnorm(move_devs[p,r,rr,y,seas,a,s], 0, exp(PE_pars[p,r,seas,a,s]), TRUE)
              } # end s loop
            } # end a loop
          } # end seas loop
        } # end y loop
      } # end p loop

    } # end r loop
  } # end rr loop

  return(loglik)
}

#' Compute Fishing Mortality Deviation Process Error Log-Likelihood (Negative Scale)
#'
#' Calculates the negative log-likelihood contribution for fishing mortality
#' deviations (\code{ln_F_devs}) under an iid, random walk, or AR1 process
#' error structure.
#'
#' \strong{Note:} Unlike \code{\link{Get_PE_loglik}} and
#' \code{\link{Get_move_PE_loglik}}, which return a single positive
#' log-likelihood scalar to be negated by the caller, this function returns
#' an already-negated array with the same dimensions as \code{ln_F_devs}
#' (one value per region/year/season/fleet cell, \code{0} where catch is not
#' used), matching the existing \code{Fmort_nLL} reporting convention.
#'
#' Random walk and AR1 do not require catch-active years to be contiguous.
#' Instead, the transition between two active years is taken over the
#' elapsed gap \eqn{d} between them, exactly the marginal transition you
#' would get from estimating deviations for the closed years in between and
#' integrating them out, without actually estimating them:
#' \describe{
#'   \item{Random walk}{\eqn{\delta_t \mid \delta_s \sim N(\delta_s, d\sigma^2)}}
#'   \item{AR1}{\eqn{\delta_t \mid \delta_s \sim N(\rho^d \delta_s, \sigma^2
#'     \sum_{i=0}^{d-1} \rho^{2i})}, where the sum has closed form
#'     \eqn{(1 - \rho^{2d}) / (1 - \rho^2)}}
#' }
#' Both reduce exactly to the standard single-step transition when \eqn{d = 1}.
#'
#' @param PE_model Integer specifying the process error structure: \code{1} =
#'   IID (deviations drawn independently as \eqn{N(0, \sigma^2)}); \code{2} =
#'   random walk (first active year initialized with a diffuse \eqn{N(0, 5)}
#'   prior); \code{3} = AR1 (first active year drawn from its stationary
#'   marginal distribution \eqn{N(0, \sigma^2 / (1 - \rho^2))}).
#' @param ln_sigmaF Array \code{[n_regions x n_seas x n_fish_fleets]} of
#'   log-scale process error SD.
#' @param Fdev_rho Array \code{[n_regions x n_seas x n_fish_fleets]} of
#'   unconstrained AR1 partial correlation (only used when \code{PE_model ==
#'   3}); transformed to \eqn{(-1, 1)} via \eqn{2 / (1 + e^{-2x}) - 1}.
#' @param ln_F_devs Array \code{[n_regions x n_years x n_seas x
#'   n_fish_fleets]} of log-scale fishing mortality deviations.
#' @param map_ln_F_devs Array \code{[n_regions x n_years x n_seas x
#'   n_fish_fleets]} mirroring \code{$map$ln_F_devs}: an estimation index
#'   where a deviation is estimated, \code{NA} where it is fixed. Only
#'   estimated deviations are penalized, and they alone form the active
#'   sequence, so a fixed cell is skipped and widens the gap \eqn{d} between
#'   the deviations either side of it. \code{\link{do_Fmort_mapping}} builds
#'   this from the catch-usage indicators, estimating a deviation wherever
#'   aggregated or any population-specific catch is used, or where the
#'   aggregate catch observation (\code{ObsCatch}) is missing (\code{NA})
#'   rather than a true recorded zero (fishing presumably continued, we
#'   simply lack a value to fit), while a true recorded zero is a real
#'   closure and is excluded.
#'
#' @return Array with the same dimensions as \code{ln_F_devs}: the negative
#'   log-likelihood contribution per cell.
#'
#' @keywords internal
#' @import RTMB
Get_Fdev_PE_loglik <- function(PE_model, ln_sigmaF, Fdev_rho, ln_F_devs, map_ln_F_devs, Fdev_pen_center = 0) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  rho_trans <- function(x) 2 / (1 + exp(-2 * x)) - 1 # constrain to (-1, 1)

  n_regions <- dim(ln_F_devs)[1]
  n_yrs <- dim(ln_F_devs)[2]
  n_seas <- dim(ln_F_devs)[3]
  n_fish_fleets <- dim(ln_F_devs)[4]

  # only estimated deviations are penalized, so mapping one off by hand removes
  # its penalty as well and drops it from the active sequence
  is_estimated <- !is.na(map_ln_F_devs)

  Fmort_nLL <- array(0, dim = dim(ln_F_devs))

  for(f in 1:n_fish_fleets) {
    for(r in 1:n_regions) {
      for(seas in 1:n_seas) {

        sigma <- exp(ln_sigmaF[r,seas,f])
        if(PE_model == 3) rho <- rho_trans(Fdev_rho[r,seas,f])

        # centering on the deviations' own mean penalizes only their spread. with a mean-plus-
        # deviations parameterization the mean already holds the level, so it is not penalized twice
        dev_mu <- 0
        if(Fdev_pen_center == 1) {
          act <- which(is_estimated[r,,seas,f])
          if(length(act) > 1) dev_mu <- sum(ln_F_devs[r,act,seas,f]) / length(act)
        }

        last_active_y <- NA # calendar year of the previous active year (NA until the first one)
        for(y in 1:n_yrs) {

          if(!is_estimated[r,y,seas,f]) next # skip cells whose deviation is fixed

          if(PE_model == 1) { # iid
            Fmort_nLL[r,y,seas,f] <- -RTMB::dnorm(ln_F_devs[r,y,seas,f], dev_mu, sigma, TRUE)
          }

          if(PE_model == 2) { # random walk
            if(is.na(last_active_y)) Fmort_nLL[r,y,seas,f] <- -RTMB::dnorm(ln_F_devs[r,y,seas,f], 0, 5, TRUE) # diffuse init
            else {
              d <- y - last_active_y # elapsed years since the previous active year
              Fmort_nLL[r,y,seas,f] <- -RTMB::dnorm(ln_F_devs[r,y,seas,f], ln_F_devs[r,last_active_y,seas,f], sigma * sqrt(d), TRUE)
            }
          }

          if(PE_model == 3) { # ar1
            if(is.na(last_active_y)) Fmort_nLL[r,y,seas,f] <- -RTMB::dnorm(ln_F_devs[r,y,seas,f], 0, sigma / sqrt(1 - rho^2), TRUE) # stationary marginal sd
            else {
              d <- y - last_active_y # elapsed years since the previous active year
              trans_sd <- sigma * sqrt((1 - rho^(2*d)) / (1 - rho^2))
              Fmort_nLL[r,y,seas,f] <- -RTMB::dnorm(ln_F_devs[r,y,seas,f], rho^d * ln_F_devs[r,last_active_y,seas,f], trans_sd, TRUE)
            }
          }

          last_active_y <- y

        } # end y loop
      } # end seas loop
    } # end r loop
  } # end f loop

  return(Fmort_nLL)
}

# Discard and Selectivity Penalties -----------------------------------------

#' Discard mortality rate deviation penalty
#'
#' IID penalty on every estimated discard mortality rate deviation, called once
#' from the "Discard Mortality Rate (Penalty)" section of \code{SPoRC_rtmb.R}.
#'
#' The penalized set is read off \code{map_logit_dmr_devs} rather than
#' recomputed, so it always matches what is estimated, including deviations
#' mapped off by hand after setup. By default
#' \code{\link{do_dmr_dev_mapping}} estimates a deviation in every cell the
#' objective does not treat as a true closure. Discard observations are not the
#' boundary: \code{dmr} is identified through total mortality (\code{ZAA})
#' wherever a cell is fished and retention is less than one, and it cancels out
#' of \code{PredDiscard} itself.
#'
#' @param logit_dmr_devs Array \code{[region, year, season, fish_fleet]} of
#'   discard mortality rate deviations on the logit scale.
#' @param ln_sigma_dmr Array \code{[region, season, fish_fleet]} of log-sigma
#'   for the deviation penalty.
#' @param map_logit_dmr_devs Array \code{[region, year, season, fish_fleet]}
#'   mirroring \code{$map$logit_dmr_devs}: an estimation index where a deviation
#'   is estimated, \code{NA} where it is fixed.
#' @param n_fish_fleets,n_yrs,n_regions,n_seas Dimension sizes.
#'
#' @return Array \code{[region, year, season, fish_fleet]} of negative
#'   log-likelihood penalties (0 where the deviation is not estimated).
#'
#' @keywords internal
#' @import RTMB
get_dmr_penalty <- function(logit_dmr_devs, ln_sigma_dmr, map_logit_dmr_devs,
                             n_fish_fleets, n_yrs, n_regions, n_seas) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  dmr_nLL <- array(0, dim = dim(logit_dmr_devs))

  # only estimated deviations are penalized, so mapping one off by hand removes
  # its penalty as well and dmr falls back on logit_dmr_mean
  is_estimated <- !is.na(map_logit_dmr_devs)

  for(f in 1:n_fish_fleets) {
    for(y in 1:n_yrs) {
      for(r in 1:n_regions) {
        for(seas in 1:n_seas) {

          if(is_estimated[r,y,seas,f]) {
            dmr_nLL[r,y,seas,f] <- -RTMB::dnorm(logit_dmr_devs[r,y,seas,f], 0, exp(ln_sigma_dmr[r,seas,f]), TRUE)
          } # end if fished

        } # end seas loop
      } # end r loop
    } # y loop
  } # f loop

  return(dmr_nLL)
}

#' Prior on selectivity, on the parameters or on realized values
#'
#' Shared across the total fishery, retained fishery, and survey "Selectivity
#' (Prior)" blocks in \code{SPoRC_rtmb.R} since all three prior tables and
#' their corresponding parameter arrays share the same
#' \code{[region, par, block, sex, fleet]} layout. Each row of the table is one
#' prior, and its optional \code{type} column selects what the row constrains:
#' \describe{
#'   \item{\code{"par"} (the default when the column is absent)}{A lognormal
#'     prior on one fixed selectivity parameter,
#'     \code{dnorm(pars[region,par,block,sex,fleet], log(mu), sd)}, with
#'     \code{mu} on the natural scale and \code{sd} on the log scale.}
#'   \item{\code{"value"}}{A normal prior on the realized selectivity value at
#'     one bin, \code{dnorm(sel[bin], mu, sd)}, with both hyperparameters on
#'     the natural scale. \code{par} instead names the bin, on the grid the
#'     data source's selectivity is parameterized on (ages or lengths per its
#'     selectivity type), and the value is read at the first model year of
#'     \code{block} (blocked and time-invariant selectivity are constant within
#'     a block). This is a constraint on a derived quantity rather than on the
#'     parameters (the ADMB rockfish convention of pinning survey selectivity
#'     at a reference age near one is its motivating case), so it can express
#'     statements no set of independent parameter priors can, e.g. the rank-one
#'     ridge in (a50, slope) space implied by constraining a logistic curve's
#'     value at one age.}
#' }
#'
#' @param selex_prior Data frame with columns \code{region}, \code{par},
#'   \code{block}, \code{sex}, \code{fleet}, \code{mu}, \code{sd}, and
#'   optionally \code{type} (\code{"par"}/\code{"value"}), one row per prior.
#' @param fixed_sel_pars Array \code{[region, par, block, sex, fleet]} of fixed
#'   selectivity parameters on the log scale, read by \code{"par"} rows.
#' @param sel Array \code{[pop, region, year, seas, age, sex, fleet]} of
#'   realized age-based selectivity, read by \code{"value"} rows at pop 1 and
#'   season 1, matching the smoothness penalties.
#' @param sel_l Array \code{[region, year, len, sex, fleet]} of realized
#'   length-based selectivity, read by \code{"value"} rows instead of
#'   \code{sel} when the data source is length-based.
#' @param selex_type Integer. \code{0} reads \code{sel}, \code{1} reads
#'   \code{sel_l}.
#' @param sel_blocks Integer array \code{[region, year, fleet]} mapping model
#'   years to selectivity blocks, used to resolve a \code{"value"} row's
#'   \code{block} to the first year in it.
#'
#' @return Numeric scalar negative log-likelihood contribution, summed across
#'   all rows of \code{selex_prior}.
#'
#' @keywords internal
#' @import RTMB
get_selex_prior <- function(selex_prior, fixed_sel_pars, sel, sel_l, selex_type, sel_blocks) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  row_type <- if(is.null(selex_prior$type)) rep("par", nrow(selex_prior)) else selex_prior$type

  nLL <- 0
  for(i in 1:nrow(selex_prior)) {
    r <- selex_prior$region[i]
    p <- selex_prior$par[i]
    b <- selex_prior$block[i]
    s <- selex_prior$sex[i]
    f <- selex_prior$fleet[i]
    if(row_type[i] == "value") {
      y <- min(which(sel_blocks[r,,f] == b)) # first model year in the block
      sel_val <- if(selex_type == 0) sel[1,r,y,1,p,s,f] else sel_l[r,y,p,s,f]
      nLL <- nLL - RTMB::dnorm(sel_val, selex_prior$mu[i], selex_prior$sd[i], TRUE)
    } else {
      nLL <- nLL - RTMB::dnorm(fixed_sel_pars[r,p,b,s,f], log(selex_prior$mu[i]), selex_prior$sd[i], TRUE)
    }
  } # end i loop

  return(nLL)
}

#' Centering penalty on a set of selectivity fixed-effect parameters
#'
#' Penalizes the squared log of the mean exponentiated value of a named set of
#' selectivity parameters, \eqn{w \, (\log \overline{e^{\theta}})^2}. A
#' non-parametric selectivity curve is only identified up to a scalar once
#' catchability or fishing mortality is free to absorb it, and this pins that
#' scalar by pushing the set's average selectivity toward one, which is a softer
#' constraint than fixing a bin outright.
#'
#' Each row of the table names one set, so the penalty applies to the group
#' jointly rather than to each parameter separately. Because the expression
#' averages on the natural scale, it is meant for parameter sets kept on the log
#' scale; a set stored on the logit scale (the non-parametric form, or the
#' asymptote of the asymptotic logistic forms) would not average to anything
#' interpretable as selectivity.
#'
#' @param selex_penalty Data frame with columns \code{region}, \code{fleet},
#'   \code{block}, \code{sex}, \code{par}, and \code{wt}, one row per penalized
#'   set. \code{par} is a list column of integer vectors naming the parameter
#'   indices in the set; \code{wt} is the weight.
#' @param fixed_sel_pars Array \code{[region, par, block, sex, fleet]} of
#'   selectivity fixed effects.
#'
#' @return Numeric scalar negative log-likelihood contribution, summed across
#'   all rows of \code{selex_penalty}.
#'
#' @keywords internal
#' @import RTMB
get_selex_fixed_penalty <- function(selex_penalty, fixed_sel_pars) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  nLL <- 0
  for(i in 1:nrow(selex_penalty)) {
    r <- selex_penalty$region[i]
    b <- selex_penalty$block[i]
    s <- selex_penalty$sex[i]
    f <- selex_penalty$fleet[i]
    p <- selex_penalty$par[[i]]
    avgsel <- log(mean(exp(fixed_sel_pars[r,p,b,s,f])))
    nLL <- nLL + selex_penalty$wt[i] * avgsel^2
  } # end i loop

  return(nLL)
}

# Recruitment Penalties -----------------------------------------------------

#' Process error density for one series of recruitment deviations
#'
#' The recruitment deviations of one population and region are either
#' independent draws about a supplied mean, a random walk, or an AR1 process.
#' The walk and the AR1 always step from the previous year, whether or not that
#' year's deviation is estimated: every year has a recruitment, so a deviation
#' fixed at a value is a year the walk passes through rather than a gap in the
#' series. A year whose own deviation is fixed contributes no density.
#'
#' A walk has no stationary distribution to start from, so year one is given a
#' diffuse normal. An AR1 starts from its stationary marginal standard deviation
#' \eqn{\sigma / \sqrt{1 - \rho^2}}. Fixing year one instead leaves the series
#' with no penalty on its level at all, which is what SAM's flat prior on the
#' first year amounts to.
#'
#' @param devs Numeric vector of deviations for one population and region, by year.
#' @param is_est Numeric vector the same length, \code{1} where the deviation is
#'   estimated and \code{0} where it is fixed.
#' @param sigma Numeric vector the same length of standard deviations, so the
#'   early and late regimes can differ. A step reads the standard deviation of
#'   the year it lands on.
#' @param dev_mu Numeric vector the same length of prior means. Only read when
#'   \code{PE_model = 1}, since a walk's mean is the previous deviation.
#' @param PE_model Integer. \code{1} independent, \code{2} random walk,
#'   \code{3} AR1.
#' @param rho AR1 correlation, already on the natural scale. Only read when
#'   \code{PE_model = 3}.
#' @param init_sd Standard deviation given to year one of a random walk. Default
#'   \code{5}, which leaves the level effectively free. \code{NA} instead starts
#'   the walk at zero under its own sigma.
#'
#' @return Numeric vector of negative log-likelihood contributions, zero where
#'   the deviation is fixed.
#'
#' @keywords internal
#' @import RTMB
get_recdev_pe_nLL <- function(devs, is_est, sigma, dev_mu, PE_model, rho = 0, init_sd = 5) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  n_yrs <- length(devs)
  nLL <- rep(0, n_yrs) * devs[1] # keep the AD type of the deviations

  for(y in 1:n_yrs) {

    if(is_est[y] == 0) next # skip cells whose deviation is fixed

    if(PE_model == 1) { # independent
      nLL[y] <- -RTMB::dnorm(devs[y], dev_mu[y], sigma[y], TRUE)
    }

    if(PE_model == 2) { # random walk
      if(y == 1) {
        # the walk needs a distribution for year one. a wide sigma leaves the level of the
        # series effectively free; NA starts the walk at zero under its own sigma
        first_sd <- if(is.na(init_sd)) sigma[y] else init_sd
        nLL[y] <- -RTMB::dnorm(devs[y], 0, first_sd, TRUE)
      }
      else nLL[y] <- -RTMB::dnorm(devs[y], devs[y-1], sigma[y], TRUE)
    }

    if(PE_model == 3) { # ar1
      if(y == 1) nLL[y] <- -RTMB::dnorm(devs[y], 0, sigma[y] / sqrt(1 - rho^2), TRUE) # stationary marginal sd
      else nLL[y] <- -RTMB::dnorm(devs[y], rho * devs[y-1], sigma[y], TRUE)
    }

  } # end y loop

  return(nLL)
}


#' Share of one penalty owed by each cell of a mapped deviation array
#'
#' Cells sharing a map level hold one parameter between them, so penalizing
#' every cell would count that parameter once per cell. Each cell takes one over
#' the number of cells at its level, and a cell mapped off takes zero.
#'
#' @param map Numeric array of map levels, \code{NA} where the cell is fixed, or
#'   \code{NULL} when every cell holds its own parameter.
#' @param dims Dimensions the weights are returned on.
#'
#' @return Numeric array of weights over \code{dims}, all ones when \code{map} is
#'   \code{NULL}.
#' @keywords internal
dev_share_weights <- function(map, dims) {
  if(is.null(map)) return(array(1, dim = dims)) # every cell holds its own parameter
  levels_by_cell <- as.vector(map) # one map level per cell, NA where fixed
  cells_at_level <- table(levels_by_cell[!is.na(levels_by_cell)]) # how many cells share each level
  wt <- ifelse(is.na(levels_by_cell), 0, 1 / as.numeric(cells_at_level[as.character(levels_by_cell)]))
  return(array(wt, dim = dims))
}

#' A deviation series' own weighted mean
#'
#' The center a penalty takes when the level of the series is left to the rest of
#' the model rather than fixed at the bias-corrected mean. Penalizing about it
#' constrains only the spread, which is what a sum of squares about the series'
#' own mean amounts to. Fewer than two penalized cells leaves no spread to
#' measure, so the center falls back to zero.
#'
#' @param devs Vector of deviations, on the log scale.
#' @param wt Numeric vector the same length, zero where the cell is not
#'   penalized.
#'
#' @return The weighted mean, or \code{0} when fewer than two cells are
#'   penalized.
#' @keywords internal
dev_own_mean <- function(devs, wt) {
  if(sum(wt) < 2) return(0) # one cell has no spread of its own to measure
  return(sum(devs * wt) / sum(wt))
}

#' Initial age deviation penalties
#'
#' Population and region specific penalties on the initial age deviations
#' (\code{ln_InitDevs}), plus the tie holding each later sex's curve near the
#' first sex's. Called from \code{\link{get_recruitment_penalty}}.
#'
#' Each penalized deviation is Gaussian on the log scale with the early
#' recruitment sigma, \eqn{-\log \phi(d_{a,s} \mid \mu, \sigma_{R,1})}, where
#'
#' \itemize{
#'   \item \eqn{d_{a,s}} is the deviation at age \eqn{a} and sex \eqn{s}, log scale, estimated.
#'   \item \eqn{\sigma_{R,1} = \exp(\code{ln_sigmaR[1,p,r]})} is the early recruitment sigma, log scale.
#'   \item \eqn{\mu} is the center, either the bias-corrected mean
#'     \eqn{-\sigma_{R,1}^2 b_a / 2} with \eqn{b_a} the bias ramp read at the year age
#'     \eqn{a} was born, or the deviations' own weighted mean pooled over ages and sexes.
#' }
#'
#' @param n_pop,n_regions,n_ages Dimension sizes.
#' @param rec_region_prop_spec Integer switch; when \code{1}, populations and
#'   regions with a fixed zero recruitment proportion are skipped.
#' @param rec_region_prop Array \code{[pop, region]} of recruitment regional
#'   apportionment.
#' @param equil_init_age_strc Integer switch naming which initial age deviations
#'   are penalized (\code{0}: none, \code{1}: all but the plus group, \code{2}:
#'   all, \code{3}: the shared subset).
#' @param ln_InitDevs Array \code{[pop, region, age, sex]} of initial age
#'   deviations. A 3-D \code{[pop, region, age]} array, the layout before the sex
#'   dim existed, is accepted and treated as one shared curve.
#' @param init_age_devs_shared Integer vector of shared initial age deviation
#'   indices, read when \code{equil_init_age_strc == 3}.
#' @param ln_sigmaR Array \code{[early/late, pop, region]} of log sigma. Initial
#'   ages read the early one.
#' @param bias_ramp Numeric vector \code{[year]} of bias ramp adjustment factors.
#' @param InitDevs_pen_center Integer. \code{1} centers on the deviations' own
#'   weighted mean, \code{0} on the bias-corrected mean.
#' @param init_devs_pen_use Array of 0/1 matching \code{ln_InitDevs}, naming
#'   which cells are penalized. Sexes sharing one parameter keep only the first
#'   sex's copy flagged so the shared parameter is not penalized twice;
#'   sex-specific deviations flag every sex. \code{NULL} penalizes only the first
#'   sex's slice, which is the pre-sex-dim behavior.
#' @param Use_init_sex_pen Integer (0/1). Whether each later sex's deviations are
#'   tied to the first sex's through a Gaussian on their difference at every
#'   penalized age. Only meaningful when the sexes have their own curves.
#' @param ln_sigma_init_sex Log standard deviation of that tie.
#' @param init_bias_ramp Numeric vector of length \code{n_ages - 1}, the bias
#'   ramp read at the year each initial age was born (deviation index
#'   \code{1 - age}). \code{NULL} reads the first model year's ramp value at
#'   every age.
#' @param map_ln_InitDevs Numeric array matching \code{ln_InitDevs} of map levels
#'   (\code{NA} where fixed). Cells sharing a level hold one parameter and split
#'   one penalty between them. \code{NULL} penalizes every cell in full.
#'
#' @return List with \code{Init_Rec_nLL} (array \code{[pop, region, age, sex]})
#'   and \code{Init_Sex_nLL} (the same layout, the between-sex tie, zero for the
#'   first sex and whenever the tie is off), each holding negative log-likelihood
#'   penalties and zero where nothing is penalized.
#'
#' @keywords internal
#' @import RTMB
get_init_devs_penalty <- function(
  n_pop,
  n_regions,
  n_ages,
  rec_region_prop_spec,
  rec_region_prop,
  equil_init_age_strc,
  ln_InitDevs,
  init_age_devs_shared,
  ln_sigmaR,
  bias_ramp,
  InitDevs_pen_center = 0,
  init_devs_pen_use = NULL,
  Use_init_sex_pen = 0,
  ln_sigma_init_sex = 0,
  init_bias_ramp = NULL,
  map_ln_InitDevs = NULL
) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  # the sex dim was added after this array, so a 3-D one is a single shared curve
  if(length(dim(ln_InitDevs)) == 3) ln_InitDevs <- array(ln_InitDevs, dim = c(dim(ln_InitDevs), 1))
  n_init_ages <- dim(ln_InitDevs)[3] # ages an initial deviation is estimated at
  n_init_sexes <- dim(ln_InitDevs)[4] # sexes with their own initial age curve

  # with no flags supplied only the first sex is penalized
  if(is.null(init_devs_pen_use)) {
    init_devs_pen_use <- array(0, dim = dim(ln_InitDevs))
    init_devs_pen_use[,,,1] <- 1
  }

  Init_Rec_nLL <- array(0, dim = dim(ln_InitDevs)) # deviations against their center
  Init_Sex_nLL <- array(0, dim = dim(ln_InitDevs)) # later sexes against the first sex

  # equil_init_age_strc 0 and 4 penalizes nothing, so the containers go back empty
  if(!equil_init_age_strc %in% c(1,2,3)) return(list(Init_Rec_nLL = Init_Rec_nLL, Init_Sex_nLL = Init_Sex_nLL))

  # ages the penalty covers: all but the plus group, all of them, or the shared subset
  if(equil_init_age_strc == 1) init_idx <- 1:(n_ages - 2)
  else if(equil_init_age_strc == 2) init_idx <- 1:n_init_ages
  else init_idx <- unique(init_age_devs_shared[!is.na(init_age_devs_shared)])

  # a cell outside the penalty owns no share of it, so it drops out of the map before the split
  init_map_active <- map_ln_InitDevs
  if(!is.null(init_map_active)) {
    init_map_active[init_devs_pen_use == 0] <- NA # sexes that share the first sex's parameter
    unpenalized_ages <- setdiff(seq_len(n_init_ages), init_idx) # ages this setting leaves out
    if(length(unpenalized_ages) > 0) init_map_active[,,unpenalized_ages,] <- NA
  }
  init_wt <- dev_share_weights(init_map_active, dim(ln_InitDevs)) # a shared deviation splits one penalty

  # the ramp read at the year each initial age was born, or the first model year's value at every age
  ramp_init <- if(is.null(init_bias_ramp)) rep(bias_ramp[1], length(init_idx)) else init_bias_ramp[init_idx]

  for(p in 1:n_pop) {
    for(r in 1:n_regions) {

      if(rec_region_prop_spec == 1 && as.numeric(rec_region_prop[p,r]) == 0) next # no recruits here, no penalty

      sigma_init <- exp(ln_sigmaR[1,p,r]) # initial ages read the early recruitment sigma

      # the center is the deviations' own mean, pooled over ages and sexes, or the bias-corrected mean
      if(InitDevs_pen_center == 1) init_mu <- dev_own_mean(ln_InitDevs[p,r,init_idx,], init_devs_pen_use[p,r,init_idx,])
      else init_mu <- -sigma_init^2 / 2 * ramp_init

      # each sex's curve against that center, at its share of one penalty
      for(s_init in 1:n_init_sexes) {
        init_dnorm <- -RTMB::dnorm(ln_InitDevs[p,r,init_idx,s_init], init_mu, sigma_init, TRUE)
        init_share <- init_devs_pen_use[p,r,init_idx,s_init] * init_wt[p,r,init_idx,s_init] # penalized, and its share
        Init_Rec_nLL[p,r,init_idx,s_init] <- init_dnorm * init_share
      } # end s_init loop

      # and the tie pulling each later sex's curve toward the first sex's at the same ages
      if(Use_init_sex_pen == 1 && n_init_sexes > 1) {
        for(s_init in 2:n_init_sexes) {
          sex_diff <- ln_InitDevs[p,r,init_idx,s_init] - ln_InitDevs[p,r,init_idx,1] # departure from the first sex
          sex_dnorm <- -RTMB::dnorm(sex_diff, 0, exp(ln_sigma_init_sex), TRUE)
          Init_Sex_nLL[p,r,init_idx,s_init] <- sex_dnorm * init_devs_pen_use[p,r,init_idx,s_init]
        } # end s_init loop
      } # end between-sex tie

    } # end r loop
  } # end p loop

  return(list(Init_Rec_nLL = Init_Rec_nLL, Init_Sex_nLL = Init_Sex_nLL))
}

#' Recruitment deviation penalties
#'
#' Population and region specific penalties on the recruitment deviations
#' (\code{ln_RecDevs}), independent or as a process over time. Called from
#' \code{\link{get_recruitment_penalty}}.
#'
#' Under \code{RecDevs_model = 1} the deviations are independent and split into
#' an early and a late sigma regime at \code{sigmaR_switch}. Both regimes are the
#' same penalty read at a different sigma over different years,
#' \eqn{-\log \phi(\varepsilon_y \mid \mu_y, \sigma_{R,k})} with an extra
#' \eqn{-(1 - b_y / 2)\log \sigma_{R,k}} when the bias ramp is on, where
#'
#' \itemize{
#'   \item \eqn{\varepsilon_y} is the deviation in year \eqn{y}, log scale, estimated.
#'   \item \eqn{k} is 1 for years before \code{sigmaR_switch} and 2 from it on.
#'   \item \eqn{\sigma_{R,k} = \exp(\code{ln_sigmaR[k,p,r]})} is that regime's sigma, log scale.
#'   \item \eqn{b_y} is the Methot and Taylor bias ramp in year \eqn{y}, between 0 and 1, data.
#'   \item \eqn{\mu_y = -\sigma_{R,k}^2 b_y / 2} is the bias-corrected mean, or the deviations'
#'     own weighted mean over the regime's years.
#' }
#'
#' The \eqn{\log \sigma} term is what makes the sigma estimable: without it a
#' larger sigma always reduces the penalty. Under \code{RecDevs_model = 2} or
#' \code{3} the deviations are a walk or an AR1 process instead, each year
#' centered on the one before it, so neither the bias ramp nor the own-mean
#' center applies.
#'
#' @param n_pop,n_regions,n_est_rec_devs Dimension sizes.
#' @param rec_region_prop_spec Integer switch; when \code{1}, populations and
#'   regions with a fixed zero recruitment proportion are skipped.
#' @param rec_region_prop Array \code{[pop, region]} of recruitment regional
#'   apportionment.
#' @param ln_sigmaR Array \code{[early/late, pop, region]} of log sigma.
#' @param bias_ramp Numeric vector \code{[year]} of bias ramp adjustment factors.
#' @param sigmaR_switch Integer year index at which the deviations switch from
#'   the early to the late sigma regime.
#' @param ln_RecDevs Array \code{[pop, region, year]} of recruitment deviations.
#' @param sigmaR2_early,sigmaR2_late Arrays \code{[pop, region]} of squared sigma
#'   used for the bias-corrected mean.
#' @param do_rec_bias_ramp Integer switch enabling the bias ramp log sigma term.
#' @param map_ln_RecDevs Array \code{[pop, region, year]} mirroring
#'   \code{map$ln_RecDevs}. Cells that are \code{NA} are fixed rather than
#'   estimated and go unpenalized; cells sharing a level split one penalty.
#'   \code{NULL} penalizes every cell in full.
#' @param RecDevs_model Integer process error structure: \code{1} independent,
#'   \code{2} random walk, \code{3} AR1.
#' @param RecDevs_rho Array \code{[pop, region]} of unconstrained AR1
#'   correlations, transformed to \eqn{(-1, 1)} here. Read when
#'   \code{RecDevs_model = 3}.
#' @param RecDevs_rw_init_sigma Standard deviation given to year one of a random
#'   walk. Default \code{5}, which leaves the level of the series effectively
#'   free. \code{NA} starts the walk at zero under its own sigma. Read when
#'   \code{RecDevs_model = 2}.
#' @param RecDevs_pen_center Integer. \code{1} centers on the deviations' own
#'   weighted mean, \code{0} on the bias-corrected mean. Read under
#'   \code{RecDevs_model = 1} only.
#'
#' @return Array \code{[pop, region, year]} of negative log-likelihood penalties,
#'   zero where nothing is penalized.
#'
#' @keywords internal
#' @import RTMB
get_rec_devs_penalty <- function(
  n_pop,
  n_regions,
  n_est_rec_devs,
  rec_region_prop_spec,
  rec_region_prop,
  ln_sigmaR,
  bias_ramp,
  sigmaR_switch,
  ln_RecDevs,
  sigmaR2_early,
  sigmaR2_late,
  do_rec_bias_ramp,
  map_ln_RecDevs = NULL,
  RecDevs_model = 1,
  RecDevs_rho = NULL,
  RecDevs_rw_init_sigma = 5,
  RecDevs_pen_center = 0
) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  Rec_nLL <- array(0, dim = dim(ln_RecDevs))

  # a deviation mapped off by hand is fixed rather than estimated, so it loses its penalty too
  is_est <- if(is.null(map_ln_RecDevs)) array(1, dim = dim(ln_RecDevs)) else array(as.numeric(!is.na(map_ln_RecDevs)), dim = dim(ln_RecDevs))
  rec_wt <- dev_share_weights(map_ln_RecDevs, dim(ln_RecDevs)) # a shared deviation splits one penalty

  early_idx <- seq_len(sigmaR_switch - 1) # years under the early sigma, empty when the switch is year one
  late_idx <- sigmaR_switch:n_est_rec_devs # years under the late sigma
  pen_idx <- seq_len(n_est_rec_devs) # every year, for a walk or an ar1

  # the ramp is a statement about the whole series, so an all-zero ramp turns the term off everywhere
  ramp_is_on <- do_rec_bias_ramp == 1 && any(bias_ramp != 0)

  for(p in 1:n_pop) {
    for(r in 1:n_regions) {

      if(rec_region_prop_spec == 1 && as.numeric(rec_region_prop[p,r]) == 0) next # no recruits here, no penalty

      if(RecDevs_model == 1) {

        # the early and late regimes are the same penalty read at a different sigma over different years
        for(regime in 1:2) {

          idx <- if(regime == 1) early_idx else late_idx
          if(length(idx) == 0) next # no early years when the switch is year one

          ln_sigma <- ln_sigmaR[regime,p,r] # this regime's log sigma
          sigmaR2 <- if(regime == 1) sigmaR2_early[p,r] else sigmaR2_late[p,r] # and its squared sigma

          # the center is the deviations' own mean over these years, or the bias-corrected mean
          if(RecDevs_pen_center == 1) dev_mu <- dev_own_mean(ln_RecDevs[p,r,idx], is_est[p,r,idx])
          else dev_mu <- -sigmaR2 / 2 * bias_ramp[idx]

          Rec_nLL[p,r,idx] <- -RTMB::dnorm(ln_RecDevs[p,r,idx], dev_mu, exp(ln_sigma), TRUE)

          # the log sigma term a bias-corrected year owes, which is what makes sigmaR estimable
          if(ramp_is_on) Rec_nLL[p,r,idx] <- Rec_nLL[p,r,idx] - (1 - 0.5 * bias_ramp[idx]) * ln_sigma

        } # end regime loop

      } # end independent recruitment deviations

      else {

        # a walk still reads the early and late sigma, so a regime switch stays available.
        # a step takes the sigma of the year it lands on
        sigma_yr <- exp(ln_sigmaR[2,p,r]) * rep(1, n_est_rec_devs)
        if(length(early_idx) > 0) sigma_yr[early_idx] <- exp(ln_sigmaR[1,p,r])

        rho <- if(RecDevs_model == 3) 2 / (1 + exp(-2 * RecDevs_rho[p,r])) - 1 else 0 # constrain to (-1, 1)

        Rec_nLL[p,r,pen_idx] <- get_recdev_pe_nLL(
          devs = ln_RecDevs[p,r,pen_idx],   # the deviations themselves
          is_est = is_est[p,r,pen_idx],     # only estimated years are stepped through
          sigma = sigma_yr,                 # early or late sigma, by year
          dev_mu = rep(0, n_est_rec_devs),  # unused by a walk, whose mean is the previous deviation
          PE_model = RecDevs_model,         # 2 random walk, 3 ar1
          rho = rho,                        # ar1 correlation, natural scale
          init_sd = RecDevs_rw_init_sigma   # sd on the first estimated year of a walk
        )

      } # end random walk or ar1 recruitment deviations

      # drop the penalty on deviations that are fixed rather than estimated, and
      # give a deviation shared across cells its share of one penalty
      Rec_nLL[p,r,] <- Rec_nLL[p,r,] * is_est[p,r,] * rec_wt[p,r,]

    } # end r loop
  } # end p loop

  return(Rec_nLL)
}

#' Recruitment and initial age deviation penalties
#'
#' The two deviation penalties the recruitment section owes, gathered so the
#' objective reads them in one place: the initial age deviations from
#' \code{\link{get_init_devs_penalty}} and the recruitment deviations from
#' \code{\link{get_rec_devs_penalty}}. Called once from the "Recruitment
#' (Penalty)" section of \code{SPoRC_rtmb.R}.
#'
#' @inheritParams get_init_devs_penalty
#' @inheritParams get_rec_devs_penalty
#'
#' @return List with \code{Init_Rec_nLL} and \code{Init_Sex_nLL} (arrays
#'   \code{[pop, region, age, sex]}) and \code{Rec_nLL} (array \code{[pop,
#'   region, year]}), each holding negative log-likelihood penalties and zero
#'   where nothing is penalized.
#'
#' @keywords internal
#' @import RTMB
get_recruitment_penalty <- function(
  n_pop,
  n_regions,
  n_ages,
  n_est_rec_devs,
  rec_region_prop_spec,
  rec_region_prop,
  equil_init_age_strc,
  ln_InitDevs,
  init_age_devs_shared,
  ln_sigmaR,
  bias_ramp,
  sigmaR_switch,
  ln_RecDevs,
  sigmaR2_early,
  sigmaR2_late,
  do_rec_bias_ramp,
  map_ln_RecDevs = NULL,
  RecDevs_model = 1,
  RecDevs_rho = NULL,
  RecDevs_rw_init_sigma = 5,
  RecDevs_pen_center = 0,
  InitDevs_pen_center = 0,
  init_devs_pen_use = NULL,
  Use_init_sex_pen = 0,
  ln_sigma_init_sex = 0,
  init_bias_ramp = NULL,
  map_ln_InitDevs = NULL
) {

  init_pen <- get_init_devs_penalty(
    n_pop = n_pop,
    n_regions = n_regions,
    n_ages = n_ages,
    rec_region_prop_spec = rec_region_prop_spec, # regions with no recruits are skipped
    rec_region_prop = rec_region_prop, # recruitment regional apportionment
    equil_init_age_strc = equil_init_age_strc, # which initial ages are penalized
    ln_InitDevs = ln_InitDevs, # initial age deviations
    init_age_devs_shared = init_age_devs_shared, # shared subset, read under setting 3
    ln_sigmaR = ln_sigmaR, # initial ages read the early sigma
    bias_ramp = bias_ramp, # bias ramp by year
    InitDevs_pen_center = InitDevs_pen_center, # own mean or bias-corrected mean
    init_devs_pen_use = init_devs_pen_use, # cells that are penalized
    Use_init_sex_pen = Use_init_sex_pen, # whether later sexes are tied to the first
    ln_sigma_init_sex = ln_sigma_init_sex, # log sd of that tie
    init_bias_ramp = init_bias_ramp, # ramp read at the year each age was born
    map_ln_InitDevs = map_ln_InitDevs # map levels, for the shared-penalty split
  )

  Rec_nLL <- get_rec_devs_penalty(
    n_pop = n_pop,
    n_regions = n_regions,
    n_est_rec_devs = n_est_rec_devs,
    rec_region_prop_spec = rec_region_prop_spec, # regions with no recruits are skipped
    rec_region_prop = rec_region_prop, # recruitment regional apportionment
    ln_sigmaR = ln_sigmaR, # early and late log sigma
    bias_ramp = bias_ramp, # bias ramp by year
    sigmaR_switch = sigmaR_switch, # year the early regime ends
    ln_RecDevs = ln_RecDevs, # recruitment deviations
    sigmaR2_early = sigmaR2_early, # squared early sigma, for the bias-corrected mean
    sigmaR2_late = sigmaR2_late, # squared late sigma, for the bias-corrected mean
    do_rec_bias_ramp = do_rec_bias_ramp, # whether the ramp's log sigma term is on
    map_ln_RecDevs = map_ln_RecDevs, # map levels, for fixed cells and the shared split
    RecDevs_model = RecDevs_model, # 1 independent, 2 random walk, 3 ar1
    RecDevs_rho = RecDevs_rho, # ar1 correlation, unconstrained scale
    RecDevs_rw_init_sigma = RecDevs_rw_init_sigma, # sd on year one of a walk
    RecDevs_pen_center = RecDevs_pen_center # own mean or bias-corrected mean
  )

  return(list(Init_Rec_nLL = init_pen$Init_Rec_nLL, Init_Sex_nLL = init_pen$Init_Sex_nLL, Rec_nLL = Rec_nLL))
}

#' Penalty on the level of the recruitment series itself
#'
#' Penalizes the log recruitment time series directly, independently of any
#' stock-recruit residual penalty. Under a stock-recruit relationship the
#' deviations are residuals about the predicted curve, so a model that also
#' wants to keep the recruitment series itself from wandering has nowhere to say
#' so; this is that second statement. Centering on the series' own mean penalizes
#' only its variability and leaves its level to the rest of the model.
#'
#' @param Rec Array \code{[pop, region, year]} of recruitment.
#' @param sigma Numeric standard deviation of the penalty.
#' @param center Integer. \code{1} centers on the mean of the log series,
#'   \code{0} centers on zero.
#' @param yrs Integer vector of years the penalty applies over, or \code{NULL}
#'   for every year.
#'
#' @return Array \code{[pop, region, year]} of negative log-likelihood
#'   contributions, zero outside \code{yrs}.
#'
#' @keywords internal
#' @import RTMB
get_rec_level_penalty <- function(Rec, sigma, center = 1, yrs = NULL) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  nLL <- array(0, dim = dim(Rec))
  use_y <- if(is.null(yrs)) 1:dim(Rec)[3] else yrs

  for(p in 1:dim(Rec)[1]) {
    for(r in 1:dim(Rec)[2]) {
      ln_rec <- log(Rec[p,r,use_y])
      mu <- if(center == 1) sum(ln_rec) / length(ln_rec) else 0
      nLL[p,r,use_y] <- -RTMB::dnorm(ln_rec, mu, sigma, TRUE)
    } # end r loop
  } # end p loop

  return(nLL)
}

#' Stock-recruit residual penalty under mean recruitment
#'
#' Compares the recruitment series against a stock-recruit curve without letting
#' the curve generate it. Under \code{rec_model = "mean_rec"} recruitment is
#' \eqn{R_y = \exp(\mu + \varepsilon_y)} and the curve enters only here, as a
#' Gaussian on the log residual \eqn{\log R_y - \log \widehat{R}_y}. That is a
#' different statement from the deviation penalty applied when the curve
#' generates recruitment: there the residual is the parameter, here it is a
#' derived quantity and the deviations remain free.
#'
#' Several AFSC models are written this way to reflect that a weakly determined SR relationship
# should inform the recruitment series rather than completly dictate it
#'
#'
#' @param Rec Array \code{[pop, region, year]} of realized recruitment.
#' @param SR_pred Array \code{[pop, region, year]} of the curve's prediction,
#'   computed alongside the population projection.
#' @param sigma Numeric standard deviation of the residual.
#' @param yrs Integer vector of years the penalty applies over, or \code{NULL}
#'   for every year. Years outside it contribute zero and stay free.
#'
#' @return Array \code{[pop, region, year]} of negative log-likelihood
#'   contributions, zero outside \code{yrs}.
#'
#' @keywords internal
#' @import RTMB
get_sr_penalty <- function(Rec, SR_pred, sigma, yrs = NULL) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  nLL <- array(0, dim = dim(Rec))
  use_y <- if(is.null(yrs)) 1:dim(Rec)[3] else yrs

  for(p in 1:dim(Rec)[1]) {
    for(r in 1:dim(Rec)[2]) {
      resid <- log(Rec[p,r,use_y]) - log(SR_pred[p,r,use_y])
      nLL[p,r,use_y] <- -RTMB::dnorm(resid, 0, sigma, TRUE)
    } # end r loop
  } # end p loop

  return(nLL)
}

# Priors --------------------------------------------------------------------

#' Normal prior on log catchability
#'
#' Shared across the fishery and survey catchability prior blocks in
#' \code{SPoRC_rtmb.R}, since both prior tables and their corresponding
#' catchability arrays share the same \code{[region, block, fleet]} layout.
#'
#' @param q_prior Data frame with columns \code{region}, \code{block},
#'   \code{fleet}, \code{mu} (prior mean, natural scale), \code{sd} (prior SD,
#'   log scale), one row per penalized parameter.
#' @param ln_q Array \code{[region, block, fleet]} of log catchability.
#'
#' @return Numeric scalar negative log-likelihood contribution, summed across
#'   all rows of \code{q_prior}.
#'
#' @keywords internal
#' @import RTMB
get_q_prior <- function(q_prior, ln_q) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  nLL <- 0
  for(i in 1:nrow(q_prior)) {
    r <- q_prior$region[i]
    b <- q_prior$block[i]
    f <- q_prior$fleet[i]
    nLL <- nLL - sum(RTMB::dnorm(ln_q[r,b,f], log(q_prior$mu[i]), q_prior$sd[i], TRUE))
  } # end i loop

  return(nLL)
}

#' Normal prior on natural mortality
#'
#' Called once from the "Natural Mortality (Prior)" section of
#' \code{SPoRC_rtmb.R}.
#'
#' @param M_prior Data frame with columns \code{popblk}, \code{regionblk},
#'   \code{yearblk}, \code{ageblk}, \code{sexblk} (block indices into
#'   \code{M_blocks}), \code{mu} (prior mean, natural scale), \code{sd} (prior
#'   SD, log scale), one row per penalized parameter. An optional
#'   \code{seasblk} column names the season block the prior applies to; without
#'   it the prior reads the block covering the first season, which is every
#'   season when mortality does not vary within the year.
#' @param ln_M Vector of estimated log natural mortality values, indexed by
#'   \code{M_blocks}.
#' @param M_blocks Array \code{[pop, region, year, season, age, sex]} mapping each
#'   population/region/year/season/age/sex cell to an index into \code{ln_M}.
#'
#' @return Numeric scalar negative log-likelihood contribution, summed across
#'   all rows of \code{M_prior}.
#'
#' @keywords internal
#' @import RTMB
get_natmort_prior <- function(M_prior, ln_M, M_blocks) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  nLL <- 0
  for(i in 1:nrow(M_prior)) {
    p <- M_prior$popblk[i]
    r <- M_prior$regionblk[i]
    b <- M_prior$yearblk[i]
    a <- M_prior$ageblk[i]
    s <- M_prior$sexblk[i]
    # priors written before seasonal M have no seasblk, so use the first
    seas <- if(is.null(M_prior$seasblk)) 1 else M_prior$seasblk[i]
    idx <- M_blocks[p,r,b,seas,a,s]
    nLL <- nLL + -RTMB::dnorm(ln_M[idx], log(M_prior$mu[i]), M_prior$sd[i], TRUE) # TMB likelihood
  } # end i loop

  return(nLL)
}

#' Scaled beta prior on steepness
#'
#' Called once from the "Steepness (Prior)" section of \code{SPoRC_rtmb.R}.
#'
#' @param h_prior Data frame with optional columns \code{lb} and \code{ub}
#'   giving the beta's support (defaulting to \code{0.2} and \code{1}), and
#'   columns \code{pop}, \code{region}, \code{mu}
#'   (prior mean steepness, natural scale), \code{sd} (prior SD, natural scale)
#'   with one row per penalized parameter.
#' @param h_trans Array \code{[pop, region]} of steepness on its transformed
#'   (0.2, 1) scale.
#'
#' @return Numeric scalar negative log-likelihood contribution, summed across
#'   all rows of \code{h_prior}.
#'
#' @keywords internal
#' @import RTMB
get_steepness_prior <- function(h_prior, h_trans) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  nLL <- 0
  for(i in 1:nrow(h_prior)) {
    p <- h_prior$pop[i]
    r <- h_prior$region[i]
    # bounding for steepness, but also can be user specified
    lb <- if(is.null(h_prior$lb)) 0.2 else h_prior$lb[i]
    ub <- if(is.null(h_prior$ub)) 1 else h_prior$ub[i]
    beta_pars <- get_beta_scaled_pars(low = lb, high = ub, mu = h_prior$mu[i], sigma = h_prior$sd[i]) # get alpha and beta parameters
    h_trans_i <- (h_trans[p,r] - lb) / (ub - lb) # transform random variable
    nLL <- nLL - RTMB::dbeta(x = h_trans_i, shape1 = beta_pars[1], shape2 = beta_pars[2], log = TRUE) # penalize
  } # end i loop

  return(nLL)
}

#' Dirichlet prior on movement rates
#'
#' Called once from the "Movement Rates (Prior)" section of \code{SPoRC_rtmb.R}.
#'
#' For CTMC movement (\code{move_type = 1}) the prior is placed on the
#' \emph{annual} movement fractions \eqn{\exp(\dot{\mathbf{Q}})}, not on the
#' seasonal fractions \eqn{\exp(\dot{\mathbf{Q}}\,\Delta t)} stored in
#' \code{Movement}. Those differ once \code{ctmc_scale_by_seasdur = 1}, and the
#' difference is not benign: a season's movement matrix approaches the identity as
#' the season shortens, so a fixed \code{alpha} silently becomes a much stronger
#' constraint as \code{n_seas} grows. On a three-region test setup the same
#' \code{alpha = 3} prior cost 1.04 nLL units at \code{n_seas = 1} but 9.91 at
#' \code{n_seas = 12}. Evaluating on the annual matrix makes \code{alpha} mean the
#' same thing regardless of seasonal structure, and matches how such priors are
#' elicited ("what fraction of fish move per year"). Under
#' \code{ctmc_scale_by_seasdur = 0} the two coincide, so nothing changes there.
#'
#' Unstructured movement (\code{move_type = 0}) has no generator; its
#' \code{Movement} entries are transition fractions in their own right and are used
#' directly, unchanged.
#'
#' @param Movement_prior Data frame with columns \code{pop}, \code{region_from},
#'   \code{year}, \code{seas}, \code{age}, \code{sex}, and \code{alpha} (list
#'   column of Dirichlet concentration vectors), one row per penalized
#'   movement-from vector. For CTMC movement, \code{seas} selects which season's
#'   generator supplies the annual fractions.
#' @param Movement Array \code{[pop, region_from, region_to, year, season, age,
#'   sex]} of movement rates.
#' @param Mrate Array of instantaneous movement rates (generator \eqn{\dot{Q}}),
#'   dimensioned like \code{Movement} and stored unscaled by season duration.
#'   \code{NULL} for unstructured movement, in which case \code{Movement} is used
#'   as-is.
#'
#' @return Numeric scalar negative log-likelihood contribution, summed across
#'   all rows of \code{Movement_prior}.
#'
#' @keywords internal
#' @import RTMB
get_movement_dirichlet_prior <- function(Movement_prior, Movement, Mrate = NULL) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  nLL <- 0
  for(i in 1:nrow(Movement_prior)) {
    p <- Movement_prior$pop[i] # population
    region_from <- Movement_prior$region_from[i] # region from
    y <- Movement_prior$year[i] # year
    seas <- Movement_prior$seas[i] # seas
    a <- Movement_prior$age[i] # age
    s <- Movement_prior$sex[i] # sex
    alpha <- Movement_prior$alpha[[i]] # get prior values

    if(is.null(Mrate)) {
      frac <- Movement[p, region_from,,y,seas,a,s] # unstructured: already fractions
    } else {
      # Mrate is stored row-convention (Mrate = t(Q_ss)), and t(expm(t(A))) = expm(A),
      # so exponentiating it directly gives the annual fractions in row convention.
      Q_ss <- methods::as(Mrate[p,,,y,seas,a,s], "sparseMatrix")
      frac <- as.matrix(Matrix::expm(Q_ss))[region_from, ]
    }

    nLL <- nLL - ddirichlet(x = frac, alpha = alpha, log = TRUE) # dirichlet prior
  } # end i loop

  return(nLL)
}

#' Normal prior on global R0
#'
#' Called once from the "Recruitment R0 (Prior)" section of \code{SPoRC_rtmb.R}.
#' Returns a scalar contribution added directly (unweighted by \code{Wt_Rec})
#' into the joint negative log-likelihood, alongside the other scalar priors
#' (\code{M_nLL}, \code{h_nLL}, etc.); it is not part of the \code{Rec_nLL}
#' recruitment-deviation array, since it penalizes a single global parameter
#' rather than a per-year deviation.
#'
#' @param r0_prior Data frame with columns \code{pop}, \code{mu} (prior mean R0,
#'   natural scale), \code{sd} (prior SD, log scale), one row per penalized
#'   population.
#' @param ln_global_R0 Vector \code{[pop]} of log mean recruitment (R0).
#'
#' @return Numeric scalar negative log-likelihood contribution, summed across
#'   all rows of \code{r0_prior}.
#'
#' @keywords internal
#' @import RTMB
get_r0_prior <- function(r0_prior, ln_global_R0) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  nLL <- 0
  for(i in 1:nrow(r0_prior)) {
    p <- r0_prior$pop[i] # population
    nLL <- nLL - RTMB::dnorm(ln_global_R0[p], log(r0_prior$mu[i]), r0_prior$sd[i], TRUE) # normal prior
  } # end i loop

  return(nLL)
}

#' Dirichlet/beta priors on recruitment apportionment
#'
#' Combines the recruitment regional apportionment prior, seasonal
#' apportionment prior, and stray rate prior, since all three feed the single
#' \code{rec_prop_nLL} accumulator in \code{SPoRC_rtmb.R}. Called once from the
#' "Recruitment Proportions (Prior)" / "Stray Rates (Prior)" sections.
#'
#' @param use_rec_region_prop_prior Integer (0/1) switch for the regional
#'   apportionment prior.
#' @param rec_region_prop_prior Data frame with columns \code{pop} and
#'   \code{alpha} (list column of Dirichlet concentration vectors).
#' @param rec_region_prop Array \code{[pop, region]} of recruitment regional
#'   apportionment.
#' @param use_rec_seas_prop_prior,use_fixed_rec_seas_prop Integer (0/1)
#'   switches; the seasonal apportionment prior is skipped when seasonal
#'   apportionment is fixed rather than estimated.
#' @param rec_seas_prop_prior Data frame with columns \code{pop} and
#'   \code{alpha} (list column of Dirichlet concentration vectors).
#' @param rec_seas_prop Array \code{[pop, season]} of recruitment seasonal
#'   apportionment.
#' @param rec_lag,spawn_seas,n_seas Integers controlling which seasons are
#'   structurally zero (age-0 recruits before the spawning event) and so
#'   excluded from the seasonal Dirichlet prior.
#' @param use_stray_rate_prior Integer (0/1) switch for the stray rate prior.
#' @param stray_rate_prior Data frame with columns \code{pop}, \code{block},
#'   \code{mu} (prior mean, natural scale), \code{sd} (prior SD, natural scale).
#' @param stray_rate_pars Array \code{[pop, block]} of stray rate parameters on
#'   the logit scale.
#'
#' @return Numeric scalar negative log-likelihood contribution, summed across
#'   all three prior sources.
#'
#' @keywords internal
#' @import RTMB
get_recruitment_proportion_priors <- function(use_rec_region_prop_prior, rec_region_prop_prior, rec_region_prop,
                                               use_rec_seas_prop_prior, use_fixed_rec_seas_prop, rec_seas_prop_prior,
                                               rec_seas_prop, rec_lag, spawn_seas, n_seas,
                                               use_stray_rate_prior, stray_rate_prior, stray_rate_pars) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  nLL <- 0

  if(use_rec_region_prop_prior == 1) { # recruitment regional apportionment
    for(i in 1:nrow(rec_region_prop_prior)) {
      p <- rec_region_prop_prior$pop[i] # population
      alpha <- rec_region_prop_prior$alpha[[i]] # get concentration values
      nLL <- nLL - ddirichlet(x = rec_region_prop[p,], alpha = alpha, log = TRUE) # dirichlet prior
    }
  }

  if(use_rec_seas_prop_prior == 1 && use_fixed_rec_seas_prop == 0) { # recruitment seasonal apportionment
    for(i in 1:nrow(rec_seas_prop_prior)) { # recruitment seasonal apportionment
      p <- rec_seas_prop_prior$pop[i] # population
      alpha <- rec_seas_prop_prior$alpha[[i]] # get concentration values
      if(rec_lag == 0 && spawn_seas > 1) {
        # seasons before spawn_seas are structurally zero (age-0 recruits
        # can't predate the spawning event that produced them); not evaluating then ...
        nLL <- nLL - ddirichlet(x = rec_seas_prop[p, spawn_seas:n_seas], alpha = alpha, log = TRUE) # dirichlet prior
      } else {
        nLL <- nLL - ddirichlet(x = rec_seas_prop[p,], alpha = alpha, log = TRUE) # dirichlet prior
      }
    }
  }

  if(use_stray_rate_prior == 1) {
    for(i in 1:nrow(stray_rate_prior)) {
      # extract indices
      p <- stray_rate_prior$pop[i]
      b <- stray_rate_prior$block[i]
      # extract beta pars
      mu <- stray_rate_prior$mu[i]
      sd <- stray_rate_prior$sd[i]
      # derive beta pars
      concentration <- mu * (1 - mu) / sd^2 - 1
      alpha <- mu * concentration
      beta <- (1 - mu) * concentration
      # extract values
      stray_rate_val <- 1e-4 + (1 - 2*1e-4) * RTMB::plogis(stray_rate_pars[p,b])
      nLL <- nLL - RTMB::dbeta(x = stray_rate_val, shape1 = alpha, shape2 = beta, log = TRUE) # penalize
    }
  }

  return(nLL)
}

#' Beta prior on tag reporting rate
#'
#' Called once from the "Tag Reporting Rate (Prior)" section of
#' \code{SPoRC_rtmb.R}. Supports both a symmetric-beta parameterization
#' (\code{type == 0}) and a mean/sd beta parameterization (\code{type == 1}).
#'
#' @param conv_tag_fishrep_prior Data frame with columns \code{region},
#'   \code{block}, \code{fleet}, \code{type} (0 = symmetric beta, 1 = mean/sd
#'   beta), \code{mu}, \code{sd}, one row per penalized parameter.
#' @param conv_tag_fish_reporting_pars Array \code{[region, block, fish_fleet]}
#'   of tag reporting rate parameters on the logit scale.
#'
#' @return Numeric scalar negative log-likelihood contribution, summed across
#'   all rows of \code{conv_tag_fishrep_prior}.
#'
#' @keywords internal
#' @import RTMB
get_tagrep_prior <- function(conv_tag_fishrep_prior, conv_tag_fish_reporting_pars) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  nLL <- 0
  for(i in 1:nrow(conv_tag_fishrep_prior)) {

    # Extract indices
    r <- conv_tag_fishrep_prior$region[i]
    b <- conv_tag_fishrep_prior$block[i]
    f <- conv_tag_fishrep_prior$fleet[i]

    conv_tag_fishrep_val <- RTMB::plogis(conv_tag_fish_reporting_pars[r,b,f]) # extract tag reporting rate value
    if(conv_tag_fishrep_prior$type[i] == 0) {
      nLL <- nLL - dbeta_symmetric(
        p_val = conv_tag_fishrep_val,
        p_ub = 1,
        p_lb = 0,
        p_prsd = conv_tag_fishrep_prior$sd[i],
        log = TRUE
      ) # penalize
    } # end if symmetric beta

    if(conv_tag_fishrep_prior$type[i] == 1) {
      # extract pars
      mu <- conv_tag_fishrep_prior$mu[i]
      sd <- conv_tag_fishrep_prior$sd[i]
      # derive beta pars
      concentration <- mu * (1 - mu) / sd^2 - 1
      alpha <- mu * concentration
      beta <- (1 - mu) * concentration
      nLL <- nLL - RTMB::dbeta(x = conv_tag_fishrep_val, shape1 = alpha, shape2 = beta, log = TRUE) # penalize
    } # end if for full beta

  } # end i loop

  return(nLL)
}

# State-Space Numbers at Age ------------------------------------------------

#' State-space numbers at age
#'
#' Penalizes the realized innovation of the centered numbers-at-age state,
#' \eqn{\eta = \log N - \log \hat{N}}, where \eqn{\hat{N}} is the deterministic
#' mortality and ageing prediction the dynamics computed into \code{NAA_pred}
#' before the state overwrote \code{NAA}.
#'
#' The state is a level rather than a deviation, so the prediction is subtracted
#' here instead of multiplying a deviation onto a value the dynamics still
#' computes, which makes the random-effects Hessian block-tridiagonal in
#' year: every term is supported on a two-year block, because \eqn{\eta_y}
#' depends on the states in years \eqn{y} and \eqn{y-1} and on nothing earlier.
#'
#' @param ln_NAA Array \code{[pop, region, year, season, age, sex]} of log numbers
#'   at the start of each season.
#' @param NAA_pred Array of the same shape holding the deterministic prediction.
#' @param sigmaNAA Array of the same shape holding the process error standard
#'   deviation for each cell, already expanded from its blocking structure.
#' @param naa_re_ages Integer vector of age indices the state is active over.
#' @param naa_re_yrs Integer vector of year indices the state is active over.
#' @param naa_re_seas Integer vector of season indices the state is active over.
#'   Season one alone is the annual state, with the numbers deterministic between
#'   seasons.
#' @param NAA_re Integer code for the structure over the age-year grid.
#'   \code{1} independent, \code{2} AR(1)
#'   over ages, \code{3} AR(1) over years with ages independent, \code{4}
#'   separable AR(1) over ages and years, \code{5} and \code{6} the
#'   three-dimensional Gaussian Markov random field on the conditional and the
#'   marginal variance respectively.
#' @param NAA_pe_pars Array \code{[pop, region, 3, sex]} of correlation
#'   parameters on the unconstrained scale, read as age, year and cohort. Unused
#'   under \code{NAA_re = 1}.
#' @param NAA_re_region Integer code for the structure across regions. \code{0}
#'   independent, \code{1} unstructured.
#' @param NAA_region_corr_pars Array \code{[pop, n_regions(n_regions-1)/2, sex]}
#'   of unconstrained parameters for the region correlation.
#' @param NAA_re_pop,NAA_re_sex Integer codes for the structure across populations
#'   and across sexes. \code{0} independent, \code{1} unstructured.
#' @param NAA_pop_corr_pars,NAA_sex_corr_pars Numeric vectors of unconstrained
#'   parameters for those correlations, one per pair. Both are global to the
#'   model rather than varying over the other dims, which is what keeps a
#'   two-level dim at exactly one parameter.
#' @param NAA_re_season Integer code for the structure across the active seasons.
#'   \code{0} independent, \code{1} unstructured.
#' @param NAA_season_corr_pars Array \code{[pop, n_k(n_k-1)/2, sex]} of
#'   unconstrained parameters for the season correlation, over the \eqn{n_k}
#'   active seasons.
#'
#' @details
#' Only the independent form admits a standard deviation that varies cell by
#' cell. Every other structure is separable or Markov in a dim, and a
#' per-cell variance is neither, so the setup function holds the year and age
#' standard deviation blocks to one apiece whenever a correlated form is chosen
#' and this function reads one standard deviation per population, region, season
#' and sex. The season dim is whitened before the age and year density is
#' reached, so it keeps its own standard deviation under every form, and gives it
#' up only under \code{NAA_re_season = 1}.
#'
#' @return Scalar negative log likelihood.
#'
#' @keywords internal
#' @import RTMB
Get_NAA_state_penalty <- function(
  ln_NAA,
  NAA_pred,
  sigmaNAA,
  naa_re_ages,
  naa_re_yrs,
  naa_re_seas,
  NAA_re = 1,
  NAA_pe_pars = NULL,
  NAA_re_region = 0,
  NAA_region_corr_pars = NULL,
  NAA_re_pop = 0,
  NAA_pop_corr_pars = NULL,
  NAA_re_sex = 0,
  NAA_sex_corr_pars = NULL,
  NAA_re_season = 0,
  NAA_season_corr_pars = NULL,
  naa_re_where = NULL
) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  d <- dim(ln_NAA)
  n_pop <- d[1]; n_regions <- d[2]; n_sexes <- d[6]
  ny <- length(naa_re_yrs); na <- length(naa_re_ages); nk <- length(naa_re_seas)

  # a population that never occupies a region holds no fish there, so it has no state to
  # penalize and the logarithm below would be taken on a structural zero
  if(is.null(naa_re_where)) naa_re_where <- base::matrix(1, n_pop, n_regions)
  keep <- array(rep(as.vector(naa_re_where), times = ny * nk * na * n_sexes), dim = c(n_pop, n_regions, ny, nk, na, n_sexes))

  pred <- NAA_pred[,,naa_re_yrs,naa_re_seas,naa_re_ages,,drop = FALSE]
  pred[keep == 0] <- 1 # never read, and it keeps the logarithm and value finite

  # compute the epsilon
  eta <- ln_NAA[,,naa_re_yrs,naa_re_seas,naa_re_ages,,drop = FALSE] - log(pred)
  eta[keep == 0] <- 0
  sig <- sigmaNAA[,,naa_re_yrs,naa_re_seas,naa_re_ages,,drop = FALSE] # get sigma NAA

  # Compute nLL for independent deviatiosn on every dimension
  if(NAA_re == 1 && NAA_re_region == 0 && NAA_re_pop == 0 && NAA_re_sex == 0 && NAA_re_season == 0)
    return(-sum(RTMB::dnorm(as.vector(eta), 0, as.vector(sig), TRUE) * as.vector(keep)))

  if(any(naa_re_where == 0) && (NAA_re_region > 0 || NAA_re_pop > 0))
    stop("naa_re_where drops a population and region cell while the numbers at age state correlates ",
         "across regions or populations. A dropped cell sits inside that joint density, so it cannot ",
         "be left out of it. Give the state every region, or turn the region and population ",
         "correlations off.")

  # initialize nll for other cases (beyond independent cases)
  nll <- 0

  # get correaltion by population (unstructured)
  if(NAA_re_pop > 0) {
    Lp <- build_us_chol(NAA_pop_corr_pars, n_pop)
    flat <- array(eta, dim = c(n_pop, n_regions * ny * nk * na * n_sexes))
    eta <- array(solve(Lp, flat), dim = c(n_pop, n_regions, ny, nk, na, n_sexes))
    nll <- nll + n_regions * ny * nk * na * n_sexes * sum(log(diag(Lp)))
  }

  # get correaltion by sexes (unstructured)
  if(NAA_re_sex > 0) {
    Ls <- build_us_chol(NAA_sex_corr_pars, n_sexes)
    flat <- array(eta, dim = c(n_pop * n_regions * ny * nk * na, n_sexes))
    eta <- array(t(solve(Ls, t(flat))), dim = c(n_pop, n_regions, ny, nk, na, n_sexes))
    nll <- nll + n_pop * n_regions * ny * nk * na * sum(log(diag(Ls)))
  }

  for(p in 1:n_pop) {
    for(s in 1:n_sexes) {

      # get epsilon
      eps <- array(eta[p,,,,,s], dim = c(n_regions, ny, nk, na))

      # get correlation by region (unstructured)
      if(NAA_re_region > 0) {
        Lc <- build_us_chol(NAA_region_corr_pars[p,,s], n_regions)
        flat <- array(eps, dim = c(n_regions, ny * nk * na))
        eps <- array(solve(Lc, flat), dim = c(n_regions, ny, nk, na))
        nll <- nll + ny * nk * na * sum(log(diag(Lc)))
      }

      # get correlation by season (unstructured). The dim is gathered by explicit slicing
      # rather than aperm, which the AD types do not have.
      if(NAA_re_season > 0) {
        Lk <- build_us_chol(NAA_season_corr_pars[p,,s], nk)
        flat <- array(0, dim = c(nk, n_regions * ny * na))
        for(k in 1:nk) flat[k,] <- as.vector(array(eps[,,k,], dim = c(n_regions, ny, na)))
        flat <- solve(Lk, flat)
        for(k in 1:nk) eps[,,k,] <- array(flat[k,], dim = c(n_regions, ny, na))
        nll <- nll + n_regions * ny * na * sum(log(diag(Lk)))
      }

      for(r in 1:n_regions) {
        if(naa_re_where[p,r] == 0) next # no state in a region this population never occupies
        for(k in 1:nk) {
          sd_prsk <- sig[p,r,1,k,1,s]
          eps_ya <- array(eps[r,,k,], dim = c(ny, na)) # year by age, matching the deviation surfaces
          nll <- nll + penalize_naa_age_year(eps_ya, sd_prsk, NAA_re, if(is.null(NAA_pe_pars)) NULL else NAA_pe_pars[p,r,,s], ny, na)
        } # end k loop
      } # end r loop

    } # end s loop
  } # end p loop

  nll
}

#' Compare one region's age-by-year innovation surface
#'
#' The age and year half of \code{\link{Get_NAA_state_penalty}}, split out so the
#' region correlation can whiten its dim and then reuse this unchanged for
#' every structure, including the three-dimensional field whose cohort term makes
#' it non-separable.
#'
#' @param eps_ya Matrix \code{[year, age]} of innovations, already whitened across
#'   regions when a region correlation is active.
#' @param sd_prs Standard deviation for this population, region and sex.
#' @param NAA_re Integer structure code, as in \code{\link{Get_NAA_state_penalty}}.
#' @param pe Numeric vector of three correlation parameters on the unconstrained
#'   scale, read as age, year and cohort.
#' @param ny,na Number of active years and ages.
#'
#' @return Scalar negative log likelihood.
#'
#' @keywords internal
#' @import RTMB
penalize_naa_age_year <- function(eps_ya, sd_prs, NAA_re, pe, ny, na) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  # 1 = independent over ages and years
  if(NAA_re == 1) return(-sum(RTMB::dnorm(as.vector(eps_ya), 0, sd_prs, TRUE)))

  # 3 = autoregression over ages, 7 = over years, 4 = separable over both. The dim a form
  # leaves alone takes an independent standard normal, a valid mean zero unit variance factor
  if(NAA_re %in% c(2, 3, 4)) {
    rho_a <- if(NAA_re %in% c(2, 4)) rho_trans(pe[1]) else 0 # 1dar1 over ages
    rho_y <- if(NAA_re %in% c(3, 4)) rho_trans(pe[2]) else 0 # 1dar1_y over years, 2dar1 over both
    scale <- sd_prs / sqrt(1 - rho_y^2) / sqrt(1 - rho_a^2) # get unit scale
    iid_f <- function(x) sum(RTMB::dnorm(x, 0, 1, TRUE)) # the dim left independent
    f_yr <- if(NAA_re %in% c(3, 4)) function(x) RTMB::dautoreg(x, mu = 0, phi = rho_y, log = TRUE) else iid_f
    f_ag <- if(NAA_re %in% c(2, 4)) function(x) RTMB::dautoreg(x, mu = 0, phi = rho_a, log = TRUE) else iid_f
    return(-RTMB::dseparable(f_yr, f_ag)(eps_ya, scale = scale))
  }

  # 5 and 6 = three-dimensional field on the conditional and the marginal variance
  if(NAA_re %in% c(5, 6)) {
    Q <- Get_3d_precision(na, ny, rho_trans(pe[1]), rho_trans(pe[2]), rho_trans(pe[3]),
                          2 * log(sd_prs), Var_Type = if(NAA_re == 5) 1 else 0)
    # Get_3d_precision numbers its nodes age fastest, so the vector is laid out the same way
    return(-RTMB::dgmrf(x = as.vector(t(eps_ya)), mu = 0, Q = Q, log = TRUE))
  }

  stop("NAA_re code ", NAA_re, " has no penalty branch.")
}
