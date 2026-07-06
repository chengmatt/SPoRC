#' Compute a model-agnostic selectivity smoothness / dome-shape penalty (Positive Scale)
#'
#' Regularization penalty operating directly on a realized selectivity-at-bin-at-year
#' surface, rather than on any particular parameterization's deviations. Because it
#' only ever looks at the resulting selectivity values, it applies uniformly to any
#' selectivity functional form -- semi-parametric process-error models
#' (\code{TimeVary_Model} 3--5, via \code{\link{Get_sel_PE_loglik}}) and the bicubic /
#' cubic spline (\code{Selex_Model == 8}, via \code{\link{Get_Selex}}) alike. This is
#' the "modular" building block behind SPoRC's independently-weighted selectivity
#' penalty terms: each term below is switched off by setting its weight to 0, and
#' terms can be combined freely without one implicit shared on/off flag.
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
#'   the next (i.e. discourages, but does not forbid, dome shapes), matching ADMB's
#'   \code{sel_like} dome penalty (\code{lambda(3)}) for double-logistic / spline
#'   selectivity forms. \code{0} (default) disables this term. Uses \code{max(., 0)}
#'   (an RTMB/CppAD-safe smooth hinge; direct \code{if()} branching on AD types is
#'   unsupported) so only decreases, not increases, are penalized.
#' @param wt_mean_center Non-negative weight on a per-year mean-centering
#'   (sum-to-zero) regularization: for each year, penalizes the squared mean of
#'   log-selectivity across bins. \code{0} (default) disables this term;
#'   set to \code{10000}.
#' @param normalize Logical. If \code{TRUE} (default), \code{wt_bin_curve} is
#'   divided by the number of bins and \code{wt_yr_diff}/\code{wt_yr_curve} are
#'   divided by the number of years. Set to \code{FALSE} to reproduce the
#'   older, unnormalized bin/year curvature penalty used by
#'   \code{\link{Get_sel_PE_loglik}} for the semi-parametric process-error
#'   models (\code{TimeVary_Model} 3--5) -- \code{Get_sel_PE_loglik} always
#'   calls this with \code{normalize = FALSE}.
#'
#' @return Numeric scalar: the positive log-likelihood contribution from the
#'   requested penalty terms. Negated externally to form the negative log-likelihood.
#'
#' @keywords internal
#' @import RTMB
Get_Selex_Smoothness_Penalty <- function(sel_vals, wt_bin_curve = 0, wt_bin_diff = 0, wt_yr_diff = 0, wt_yr_curve = 0, wt_dome = 0, wt_mean_center = 0, normalize = TRUE) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  ll = 0 # initialize likelihood (positive scale, negated by the caller)

  n_yrs = dim(sel_vals)[2]
  n_bins = dim(sel_vals)[3]
  n_sexes = dim(sel_vals)[4]

  bin_norm = if(normalize) n_bins else 1
  yr_norm = if(normalize) n_yrs else 1

  if(wt_bin_curve != 0 && n_bins >= 3) { # age/bin curvature (second difference across bins), normalized by n_bins as in ADMB
    for(s in 1:n_sexes) {
      for(y in 1:n_yrs) {
        for(b in 2:(n_bins - 1)) {
          bin_penalty = log(sel_vals[1,y,b+1,s,1]) - 2 * log(sel_vals[1,y,b,s,1]) + log(sel_vals[1,y,b-1,s,1])
          ll = ll - wt_bin_curve / bin_norm * bin_penalty^2
        } # end b loop
      } # end y loop
    } # end s loop
  }

  if(wt_bin_diff != 0 && n_bins >= 2) { # unconditional bin first-difference (both directions penalized, unlike wt_dome), normalized by n_bins
    for(s in 1:n_sexes) {
      for(y in 1:n_yrs) {
        for(b in 1:(n_bins - 1)) {
          bin_diff_penalty = log(sel_vals[1,y,b,s,1]) - log(sel_vals[1,y,b+1,s,1])
          ll = ll - wt_bin_diff / bin_norm * bin_diff_penalty^2
        } # end b loop
      } # end y loop
    } # end s loop
  }

  if(wt_yr_diff != 0 && n_yrs >= 2) { # inter-annual first difference, normalized by n_yrs as in ADMB
    for(s in 1:n_sexes) {
      for(b in 1:n_bins) {
        for(y in 2:n_yrs) {
          yr_diff_penalty = log(sel_vals[1,y,b,s,1]) - log(sel_vals[1,y-1,b,s,1])
          ll = ll - wt_yr_diff / yr_norm * yr_diff_penalty^2
        } # end y loop
      } # end b loop
    } # end s loop
  }

  if(wt_yr_curve != 0 && n_yrs >= 3) { # inter-annual second difference / smoothness, normalized by n_yrs as in ADMB
    for(s in 1:n_sexes) {
      for(b in 1:n_bins) {
        for(y in 2:(n_yrs - 1)) {
          year_penalty = log(sel_vals[1,y+1,b,s,1]) - 2 * log(sel_vals[1,y,b,s,1]) + log(sel_vals[1,y-1,b,s,1])
          ll = ll - wt_yr_curve / yr_norm * year_penalty^2
        } # end y loop
      } # end b loop
    } # end s loop
  }

  if(wt_dome != 0 && n_bins >= 2) { # dome-shape / non-monotonicity, within each year
    for(s in 1:n_sexes) {
      for(y in 1:n_yrs) {
        for(b in 1:(n_bins - 1)) {
          decrease = max(log(sel_vals[1,y,b,s,1]) - log(sel_vals[1,y,b+1,s,1]), 0) # only decreases contribute
          ll = ll - wt_dome * decrease^2
        } # end b loop
      } # end y loop
    } # end s loop
  }

  if(wt_mean_center != 0) { # per-year mean-centering / sum-to-zero regularization
    for(s in 1:n_sexes) {
      for(y in 1:n_yrs) {
        z = mean(log(sel_vals[1,y,,s,1]))
        ll = ll - wt_mean_center * z^2
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
#' Independently-weighted regularization penalties can also be applied via
#' \code{pen_wts} (see \code{\link{Get_Selex_Smoothness_Penalty}} for the bin/year
#' curvature terms, which this function delegates to for \code{PE_model} 3--5):
#' \itemize{
#'   \item For \code{PE_model} 1--2: \code{pen_wts["yr_devs"]} weights a first-difference
#'     penalty on log-deviations across years.
#'   \item For \code{PE_model} 3--5: \code{pen_wts["bin_curve"]} and \code{pen_wts["yr_curve"]}
#'     independently weight second-difference (smoothness) penalties on log-selectivity
#'     across bins and across years, respectively.
#' }
#' Each weight defaults to \code{0} (off); set any subset of them to apply only the
#' penalty terms desired, rather than one shared on/off flag.
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
#'   \item Models 1--2: \code{[1,1,s,1]} = log standard deviation (\eqn{\log \sigma})
#'     for sex \code{s}, indexed by bin/age.
#'   \item Models 3--4: \code{[1,1,s,1]} = unconstrained partial correlation by age/bin;
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
#' @param sel_vals Array of selectivity values dimensioned \code{[1, year, bin, sex, 1]},
#'   used on the log scale when computing bin and year smoothness penalties
#'   (\code{do_sel_pen = TRUE}).
#'
#' @param pen_wts Named numeric vector with elements \code{"yr_devs"}, \code{"bin_curve"},
#'   \code{"yr_curve"} (any missing name is treated as \code{0}). Independently weights
#'   the additional regularization penalties applied beyond the process error
#'   likelihood: \code{"yr_devs"} weights a first-difference-across-years penalty on
#'   \code{ln_devs} (models 1--2 only); \code{"bin_curve"} and \code{"yr_curve"} weight
#'   second-difference (curvature) penalties on log-selectivity across bins and across
#'   years, respectively (models 3--5 only, via \code{\link{Get_Selex_Smoothness_Penalty}}).
#'   A weight of \code{0} disables that term.
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
                              sel_vals,
                              pen_wts,
                              min_sel_devs_shared_bins
                              ) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  # Note that the likelihood calculations are positive within the function,
  # because it gets converted to negative outside the wrapper function

  # Named weights default to 0 (off) for any term not supplied
  wt_yr_devs = if("yr_devs" %in% names(pen_wts)) pen_wts[["yr_devs"]] else 0
  wt_bin_curve = if("bin_curve" %in% names(pen_wts)) pen_wts[["bin_curve"]] else 0
  wt_yr_curve = if("yr_curve" %in% names(pen_wts)) pen_wts[["yr_curve"]] else 0

  ll = 0 # initialize likelihood

  # find unique selectivity deviations to penalize (sort drops NAs)
  unique_sel_devs = sort(unique(as.vector(map_sel_devs)))
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
        if(y == 1) ll = ll + RTMB::dnorm(ln_devs[1,y,i,s,1], 0, 5, TRUE) # initialize w/ big value
        else ll = ll + RTMB::dnorm(ln_devs[1,y,i,s,1], ln_devs[1,y-1,i,s,1], exp(PE_pars[1,i,s,1]), TRUE)
      } # end random walk process error

    } # end dev_idx loop

    # Temporal Stability Penalty (independently weighted; 0 disables it)
    if(wt_yr_devs != 0) {
      for(y in 2:n_yrs) {
        for(s in 1:n_sexes) {
          for(b in 1:n_bins) {
            year_penalty = ln_devs[1,y,b,s,1] - ln_devs[1,y-1,b,s,1]
            ll = ll - wt_yr_devs * year_penalty^2
          } # end b loop
        } # end s loop
      } # end y loop
    }

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

    # Age/bin curvature and year curvature penalties, independently weighted (0 disables a term);
    # delegated to the model-agnostic Get_Selex_Smoothness_Penalty so the same terms can also be
    # applied directly to non-devs-based selectivity forms (e.g. the bicubic spline). normalize =
    # FALSE preserves the exact (unnormalized) penalty magnitude this always used, predating the
    # ADMB-aligned smoothness penalty package (smooth_*), which normalizes by n_bins/n_yrs (normalize = TRUE).
    ll = ll + Get_Selex_Smoothness_Penalty(sel_vals, wt_bin_curve = wt_bin_curve, wt_yr_curve = wt_yr_curve, normalize = FALSE)

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
#'   standard deviation. Models 1--5 are single-population (fix \code{pop = 1});
#'   models 6--10 estimate separate parameters per population:
#' \itemize{
#'   \item \strong{1} – IID across years (single \eqn{\sigma} per origin region)
#'   \item \strong{2} – IID across ages (single \eqn{\sigma} per origin region and age)
#'   \item \strong{3} – IID across years and ages
#'   \item \strong{4} – IID across years, ages, and sexes
#'   \item \strong{5} – IID across years, seasons, ages, and sexes
#'   \item \strong{6} – IID across populations and years
#'   \item \strong{7} – IID across populations and ages
#'   \item \strong{8} – IID across populations, years, and ages
#'   \item \strong{9} – IID across populations, years, ages, and sexes
#'   \item \strong{10} – IID across populations, years, seasons, ages, and sexes
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
#' elapsed gap \eqn{d} between them -- exactly the marginal transition you
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
#' @param UseCatch,UseCatch_pop Same binary catch-usage indicators as
#'   elsewhere; a cell is penalized if aggregated catch or any
#'   population-specific catch is used.
#' @param missing_catch Logical array \code{[n_regions x n_years x n_seas x
#'   n_fish_fleets]}, \code{TRUE} where the aggregate catch observation
#'   (\code{ObsCatch}) is missing (\code{NA}) rather than a true recorded
#'   zero. A cell with \code{UseCatch == 0} (and no population-specific
#'   catch used) is still penalized as an ordinary active year when
#'   \code{missing_catch} is \code{TRUE} there (fishing presumably
#'   continued, we simply lack a value to fit), as opposed to a true
#'   recorded zero, which is treated as a real closure and excluded from
#'   the gap count.
#'
#' @return Array with the same dimensions as \code{ln_F_devs}: the negative
#'   log-likelihood contribution per cell.
#'
#' @keywords internal
#' @import RTMB
Get_Fdev_PE_loglik <- function(PE_model, ln_sigmaF, Fdev_rho, ln_F_devs, UseCatch, UseCatch_pop, missing_catch) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  rho_trans <- function(x) 2 / (1 + exp(-2 * x)) - 1 # constrain to (-1, 1)

  n_regions <- dim(ln_F_devs)[1]
  n_yrs <- dim(ln_F_devs)[2]
  n_seas <- dim(ln_F_devs)[3]
  n_fish_fleets <- dim(ln_F_devs)[4]

  # a cell is penalized if aggregated OR any population-specific catch is
  # used, or if the aggregate observation is missing (see missing_catch above)
  has_catch <- UseCatch == 1 | apply(UseCatch_pop == 1, c(2,3,4,5), any) | missing_catch

  Fmort_nLL <- array(0, dim = dim(ln_F_devs))

  for(f in 1:n_fish_fleets) {
    for(r in 1:n_regions) {
      for(seas in 1:n_seas) {

        sigma <- exp(ln_sigmaF[r,seas,f])
        if(PE_model == 3) rho <- rho_trans(Fdev_rho[r,seas,f])

        last_active_y <- NA # calendar year of the previous catch-active year (NA until the first one)
        for(y in 1:n_yrs) {

          if(!has_catch[r,y,seas,f]) next # skip cells without catch

          if(PE_model == 1) { # iid
            Fmort_nLL[r,y,seas,f] <- -RTMB::dnorm(ln_F_devs[r,y,seas,f], 0, sigma, TRUE)
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
