# Stage 2 of 3: objective function
#
# Every prior and penalty that contributes to the objective, one small function
# each: the process error likelihoods for selectivity, movement and fishing
# mortality deviations, and priors on natural mortality, steepness, catchability,
# R0, recruitment apportionment, movement, tag reporting rate, discard mortality
# rate and fixed selectivity.
#
# These return positive log likelihoods; the objective negates them.

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
Get_Selex_Smoothness_Penalty <- function(sel_vals, wt_bin_curve = 0, wt_bin_diff = 0,
                                         wt_yr_diff = 0, wt_yr_curve = 0, wt_dome = 0, wt_mean_center = 0,
                                         normalize = TRUE, bin_range = NULL, yr_diff_ref = NULL) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  ll = 0 # initialize likelihood (positive scale, negated by the caller)

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

  # bin_range and normalize are either one setting shared by every term, or a
  # named list giving each term its own. A shape penalty confined to the older
  # ages and a random walk spanning every age are then specified together.
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
          ll = ll - wt_bin_curve[y] / bin_norm * bin_penalty^2
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
          ll = ll - wt_bin_diff[y] / bin_norm * bin_diff_penalty^2
        } # end b loop
      } # end y loop
    } # end s loop
  }

  bins = get_bins("smooth_yr_diff"); b_lo = bins[1]; b_hi = bins[2]
  yr_norm = if(get_norm("smooth_yr_diff")) n_yrs else 1
  # The walk has no previous value to penalize in in its first year, so that year is
  # normally left unpenalized. A reference supplies one: the first penalized
  # year is then held toward yr_diff_ref on the log scale, which anchors an
  # otherwise free series to a known selectivity before the data begin.
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
          ll = ll - wt_yr_diff[y] / yr_norm * yr_diff_penalty^2
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
          ll = ll - wt_yr_curve[y] / yr_norm * year_penalty^2
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
          ll = ll - wt_dome[y] * decrease^2
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
        ll = ll - wt_mean_center[y] * z^2
      } # end y loop
    } # end s loop
  }

  return(ll)
} # return log likelihood

#' Compute Selectivity Process Error Log-Likelihood (Positive Scale)
#'
#' Calculates the positive log-likelihood contribution for selectivity
#' process error deviations under a variety of temporal/spatiotemporal
#' structures.
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
#'   carry the same integer value; \code{NA} entries are treated as fixed
#'   and excluded from likelihood evaluation.
#'
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
Get_sel_PE_loglik <- function(PE_model,
                              PE_pars,
                              ln_devs,
                              map_sel_devs,
                              min_sel_devs_shared_bins,
                              rw_init_sigma = 5
                              ) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  # Note that the likelihood calculations are positive within the function,
  # because it gets converted to negative outside the wrapper function

  ll = 0 # initialize likelihood

  # find unique selectivity deviations to penalize (sort drops NAs)
  unique_sel_devs = sort(unique(as.vector(map_sel_devs)))

  # Exit out fxn if this region x fleet slice has no estimated deviations at all.
  if(length(unique_sel_devs) == 0) return(ll)
  n_yrs = dim(map_sel_devs)[2] # get years for indexing
  n_bins = dim(map_sel_devs)[3] # get bins / pars for indexing
  n_sexes = dim(map_sel_devs)[4] # get sexes for indexing

  if(PE_model %in% c(1, 2)) {

    for(dev_idx in 1:length(unique_sel_devs)) {

      # figure out where unique sel devs first occur
      idx = which(map_sel_devs == unique_sel_devs[dev_idx], arr.ind = TRUE)[1,]
      y = idx[2] # get unique year index
      i = idx[3] # get unique age or parmeter index
      s = idx[4] # get unique sex index

      if(PE_model == 1) {
        if(y >= 1) ll = ll + RTMB::dnorm(ln_devs[1,y,i,s,1], 0, exp(PE_pars[1,i,s,1]), TRUE)
      } # iid process error

      if(PE_model == 2) {
        # The walk needs a distribution for its first year. A wide sigma leaves
        # that year essentially free; NA instead starts the walk at zero under
        # its own sigma, which is what a first difference taken against a
        # selectivity of one amounts to.
        if(y == 1) {
          init_sd = if(is.na(rw_init_sigma)) exp(PE_pars[1,i,s,1]) else rw_init_sigma
          ll = ll + RTMB::dnorm(ln_devs[1,y,i,s,1], 0, init_sd, TRUE)
        }
        else ll = ll + RTMB::dnorm(ln_devs[1,y,i,s,1], ln_devs[1,y-1,i,s,1], exp(PE_pars[1,i,s,1]), TRUE)
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

      # Construct precision matrix for 3d gmrf
      if(PE_model %in% c(3,4)) {
        Q = Get_3d_precision(n_ages = n_bins[min_sel_devs_shared_bins], # number of ages
                             n_yrs = n_yrs,  # number of years
                             pcorr_age = PE_pars[1,1,s,1], # unconstrained partial correlation by age
                             pcorr_year = PE_pars[1,2,s,1], # unconstrained partial correlation by year
                             pcorr_cohort = PE_pars[1,3,s,1], # unconstrained partial correlation by cohort
                             ln_var_value = PE_pars[1,4,s,1], # log variance
                             Var_Type = Var_Type
                             ) # variance type, == 0 (marginal), == 1 (conditional)

        # apply gmrf likelihood
        eps_ay = as.vector(t(ln_devs[1,,min_sel_devs_shared_bins,s,1])) # convert to vector
        ll = ll + RTMB::dgmrf(x = eps_ay, mu = 0, Q = Q, log = TRUE)
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
        ll = ll + RTMB::dseparable(f1, f2)(eps_ya, scale = scale)
      } # end if
    } # end idx loop

  } # end 3dgrmf or 2dar1 process error

  return(ll)
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
#'   carry the same integer value; dimensions are extracted from this
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

  ll = 0 # initialize likelihood

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

  # Dimensions named in each PE_model's spec (see cont_move_map in
  # Setup_Movement.R for the type <-> integer correspondence); dims absent
  # from a model's name are shared/broadcast across (held fixed at index 1
  # for PE_pars, which has no region_to or year dimension of its own).
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
                ll = ll + RTMB::dnorm(move_devs[p,r,rr,y,seas,a,s], 0, exp(PE_pars[p,r,seas,a,s]), TRUE)
              } # end s loop
            } # end a loop
          } # end seas loop
        } # end y loop
      } # end p loop

    } # end r loop
  } # end rr loop

  return(ll)
}

#' Compute Fishing Mortality Deviation Process Error Log-Likelihood (Negative Scale)
#'
#' Calculates the negative log-likelihood contribution for fishing mortality
#' deviations (\code{ln_F_devs}) under an iid, random walk, or AR1 process
#' error structure.
#'
#' \strong{Note:} Unlike \code{\link{Get_sel_PE_loglik}} and
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

        # Centring on the deviations' own mean penalizes only their spread, leaving
        # their level free. With a mean-plus-deviations parameterization the level
        # is already carried by the mean, so this is the form that does not
        # penalize it twice.
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
#'     stream's selectivity is parameterized on (ages or lengths per its
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
#'   \code{sel} when the stream is length-based.
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
#' averages on the natural scale, it is meant for parameter sets held on the log
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

#' Recruitment and initial age deviation penalties
#'
#' Population/region-specific IID penalties on initial age deviations
#' (\code{ln_InitDevs}) and on early/late recruitment deviations
#' (\code{ln_RecDevs}), including the Methot & Taylor bias-ramp adjustment.
#' Called once from the "Recruitment (Penalty)" section of \code{SPoRC_rtmb.R}.
#'
#' @param n_pop,n_regions,n_ages,n_est_rec_devs Dimension sizes.
#' @param rec_dd Integer recruitment density-dependence switch (used only to
#'   pick \code{sigma_idx} when \code{n_pop == 1}).
#' @param natal_region Integer vector \code{[pop]} of natal region indices.
#' @param rec_region_prop_spec Integer switch; when \code{1}, populations/regions
#'   with a fixed zero recruitment proportion are skipped.
#' @param rec_region_prop Array \code{[pop, region]} of recruitment regional
#'   apportionment.
#' @param equil_init_age_strc Integer switch selecting which initial age
#'   deviations are penalized (\code{1}: all but plus group, \code{2}: all,
#'   \code{3}: shared subset).
#' @param ln_InitDevs Array \code{[pop, region, age]} of initial age deviations.
#' @param init_age_devs_shared Integer vector of shared initial age deviation
#'   indices (used when \code{equil_init_age_strc == 3}).
#' @param ln_sigmaR Array \code{[early/late, pop, region]} of log-sigma for
#'   recruitment deviations.
#' @param bias_ramp Numeric vector \code{[year]} of bias-ramp adjustment factors.
#' @param sigmaR_switch Integer year index at which recruitment deviations
#'   switch from the early to the late sigma regime.
#' @param ln_RecDevs Array \code{[pop, region, year]} of recruitment deviations.
#' @param sigmaR2_early,sigmaR2_late Arrays \code{[pop, region]} of squared
#'   sigma used for the bias-ramp mean offset.
#' @param do_rec_bias_ramp Integer switch enabling the bias-ramp log-sigma
#'   adjustment.
#' @param map_ln_RecDevs Array \code{[pop, region, year]} mirroring
#'   \code{map$ln_RecDevs}; cells that are \code{NA} are fixed rather than
#'   estimated and are left unpenalized. \code{NULL} penalizes every cell.
#'
#' @return List with elements \code{Init_Rec_nLL} (array \code{[pop, region,
#'   age]}) and \code{Rec_nLL} (array \code{[pop, region, year]}), each holding
#'   negative log-likelihood penalties (0 where not penalized).
#'
#' @keywords internal
#' @import RTMB
get_recruitment_penalty <- function(n_pop, n_regions, n_ages, n_est_rec_devs, rec_dd, natal_region,
                                     rec_region_prop_spec, rec_region_prop, equil_init_age_strc,
                                     ln_InitDevs, init_age_devs_shared, ln_sigmaR, bias_ramp,
                                     sigmaR_switch, ln_RecDevs, sigmaR2_early, sigmaR2_late,
                                     do_rec_bias_ramp, map_ln_RecDevs = NULL,
                                     RecDevs_pen_center = 0, InitDevs_pen_center = 0) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  Init_Rec_nLL <- array(0, dim = dim(ln_InitDevs))
  Rec_nLL <- array(0, dim = dim(ln_RecDevs))

  # only estimated deviations are penalized, so mapping one off by hand removes its penalty as well.
  is_est <- if(is.null(map_ln_RecDevs)) array(1, dim = dim(ln_RecDevs)) else array(as.numeric(!is.na(map_ln_RecDevs)), dim = dim(ln_RecDevs))

  # The prior mean of a deviation is either asserted (zero, or the bias-corrected
  # -sigma^2/2) or estimated from the deviations themselves. The second form
  # penalizes only the spread and leaves the level unconstrained, which is what a
  # sum of squares about the series' own mean amounts to.
  own_mean <- function(x, w) {
    if(sum(w) < 2) return(0)
    sum(x * w) / sum(w)
  }

  for(p in 1:n_pop) {
    for(r in 1:n_regions) {

      # get sigma index
      sigma_idx <- ifelse(n_pop == 1 && rec_dd == 0, r, natal_region[p])

      # Skip penalty if no dispersal and p = r has no recruits
      if(rec_region_prop_spec == 1 && as.numeric(rec_region_prop[p,r]) == 0) next

      # Initial age deviations (if equil_init_age_strc == 0; don't penalize at all).
      # init_idx is only computed inside the guard: the shared-subset case reads
      # init_age_devs_shared, which a model without that option never carries.
      if(equil_init_age_strc %in% c(1,2,3)) {
        init_idx <- if(equil_init_age_strc == 1) 1:(n_ages - 2) else if(equil_init_age_strc == 2) 1:dim(ln_InitDevs)[3] else unique(init_age_devs_shared[!is.na(init_age_devs_shared)])
        init_mu <- if(InitDevs_pen_center == 1) own_mean(ln_InitDevs[p,r,init_idx], rep(1, length(init_idx))) else -exp(ln_sigmaR[1,p,sigma_idx])^2/2 * bias_ramp[1]
        Init_Rec_nLL[p,r,init_idx] <- -RTMB::dnorm(ln_InitDevs[p,r,init_idx], init_mu, exp(ln_sigmaR[1,p,sigma_idx]), TRUE)
      }

      # Early recruitment deviations
      if(sigmaR_switch > 1) {
        e_idx <- 1:(sigmaR_switch-1)
        e_mu <- if(RecDevs_pen_center == 1) own_mean(ln_RecDevs[p,r,e_idx], is_est[p,r,e_idx]) else -sigmaR2_early[p,sigma_idx]/2 * bias_ramp[e_idx]
        Rec_nLL[p,r,e_idx] <- -RTMB::dnorm(ln_RecDevs[p,r,e_idx], e_mu, exp(ln_sigmaR[1,p,sigma_idx]), TRUE)
        if(do_rec_bias_ramp == 1 && any(bias_ramp != 0)) Rec_nLL[p,r,1:(sigmaR_switch-1)] <- Rec_nLL[p,r,1:(sigmaR_switch-1)] - (1 - 0.5 * bias_ramp[1:(sigmaR_switch-1)]) * ln_sigmaR[1,p,sigma_idx] # adjust w/ bias correction
      }

      # Late recruitment deviations
      l_idx <- sigmaR_switch:n_est_rec_devs
      l_mu <- if(RecDevs_pen_center == 1) own_mean(ln_RecDevs[p,r,l_idx], is_est[p,r,l_idx]) else -sigmaR2_late[p,sigma_idx]/2 * bias_ramp[l_idx]
      Rec_nLL[p,r,l_idx] <- -RTMB::dnorm(ln_RecDevs[p,r,l_idx], l_mu, exp(ln_sigmaR[2,p,sigma_idx]), TRUE)
      if(do_rec_bias_ramp == 1 && any(bias_ramp != 0)) Rec_nLL[p,r,sigmaR_switch:n_est_rec_devs] <- Rec_nLL[p,r,sigmaR_switch:n_est_rec_devs] - (1 - 0.5 * bias_ramp[sigmaR_switch:n_est_rec_devs]) * ln_sigmaR[2,p,sigma_idx] # adjust w/ bias correction

      # drop the penalty on deviations that are fixed rather than estimated
      Rec_nLL[p,r,] <- Rec_nLL[p,r,] * is_est[p,r,]

    } # end r loop
  } # end p loop

  return(list(Init_Rec_nLL = Init_Rec_nLL, Rec_nLL = Rec_nLL))
}

#' Penalty on the level of the recruitment series itself
#'
#' Penalizes the log recruitment time series directly, independently of any
#' stock-recruit residual penalty. Under a stock-recruit relationship the
#' deviations are residuals about the predicted curve, so a model that also
#' wants to keep the recruitment series itself from wandering has nowhere to say
#' so; this is that second statement. Centring on the series' own mean penalizes
#' only its variability and leaves its level to the rest of the model.
#'
#' @param Rec Array \code{[pop, region, year]} of recruitment.
#' @param sigma Numeric standard deviation of the penalty.
#' @param center Integer. \code{1} centres on the mean of the log series,
#'   \code{0} centres on zero.
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
#'   SD, log scale), one row per penalized parameter.
#' @param ln_M Vector of estimated log natural mortality values, indexed by
#'   \code{M_blocks}.
#' @param M_blocks Array \code{[pop, region, year, age, sex]} mapping each
#'   population/region/year/age/sex cell to an index into \code{ln_M}.
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
    idx <- M_blocks[p,r,b,a,s]
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
    # The support defaults to the usual (0.2, 1) for steepness, but is settable,
    # because a beta placed on (0, 1) instead is a different function of h rather
    # than the same one shifted: it carries log(h) where the rescaled form
    # carries log(h - 0.2), so no choice of shape reconciles the two.
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
#' constraint as \code{n_seas} grows. On a three-region fixture the same
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
      nLL <- nLL - dbeta_symmetric(p_val = conv_tag_fishrep_val, p_ub = 1, p_lb = 0, p_prsd = conv_tag_fishrep_prior$sd[i], log = TRUE) # penalize
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
