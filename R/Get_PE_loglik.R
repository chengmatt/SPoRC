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
#' When \code{do_sel_pen = TRUE}, additional regularization penalties are applied:
#' \itemize{
#'   \item For \code{PE_model} 1--2: first-difference penalty on log-deviations across years.
#'   \item For \code{PE_model} 3--5: second-difference (smoothness) penalty on log-selectivity
#'     across bins and across years.
#' }
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
#' @param do_sel_pen Logical. If \code{TRUE}, applies additional regularization penalties
#'   beyond the process error likelihood. For models 1--2, penalizes first differences
#'   of log-deviations across years. For models 3--5, penalizes second differences
#'   (curvature) of log-selectivity across bins and across years.
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
                              do_sel_pen,
                              min_sel_devs_shared_bins
                              ) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  # Note that the likelihood calculations are positive within the function,
  # because it gets converted to negative outside the wrapper function

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

    # Temporal Stability Penalty
    if(do_sel_pen == TRUE) {
      for(y in 2:n_yrs) {
        for(s in 1:n_sexes) {
          for(b in 1:n_bins) {
            year_penalty = ln_devs[1,y,b,s,1] - ln_devs[1,y-1,b,s,1]
            ll = ll - year_penalty^2
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

    if(do_sel_pen == TRUE) {
      # Regularity on bins
      for(s in 1:n_sexes) {
        for (y in 1:n_yrs) {
          if (n_bins >= 3) {
            for (b in 2:(n_bins - 1)) {
              bin_penalty = log(sel_vals[1,y,b+1,s,1]) - 2 * log(sel_vals[1,y,b,s,1]) + log(sel_vals[1,y,b-1,s,1])
              ll = ll - bin_penalty^2
            } # end b loop
          } # end if
        } # end y loop
      } # end s loop

      # Regularity on years
      for(s in 1:n_sexes) {
        for(b in 1:n_bins) {
          if(n_yrs >= 3) {
            for(y in 2:(n_yrs - 1)) {
              year_penalty = log(sel_vals[1,y+1,b,s,1]) - 2 * log(sel_vals[1,y,b,s,1]) + log(sel_vals[1,y-1,b,s,1])
              ll = ll - year_penalty^2
            } # end y loop
          } # end if n_yrs >= 3
        } # end a loop
      }
    }

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

  # Penalize Deviations
  for(rr in 1:n_regions_to) {
    for(r in 1:n_regions_from) {

      if(adjacency_collapsed[r,rr] == 0) next # skip

      if(PE_model == 1) {
        for(y in 1:n_yrs) {
          ll = ll + RTMB::dnorm(move_devs[1,r,rr,y,1,1,1], 0, exp(PE_pars[1,r,1,1,1]), TRUE)
        } # end y loop
      } # iid_y

      if(PE_model == 2) {
        for(a in age_start:n_ages) {
          ll = ll + RTMB::dnorm(move_devs[1,r,rr,1,1,a,1], 0, exp(PE_pars[1,r,1,a,1]), TRUE)
        } # end a loop
      } # iid_a

      if(PE_model == 3) {
        for(y in 1:n_yrs) {
          for(a in age_start:n_ages) {
            ll = ll + RTMB::dnorm(move_devs[1,r,rr,y,1,a,1], 0, exp(PE_pars[1,r,1,a,1]), TRUE)
          } # end a loop
        } # end y loop
      } # iid_y_a

      if(PE_model == 4) {
        for(y in 1:n_yrs) {
          for(a in age_start:n_ages) {
            for(s in 1:n_sexes) {
              ll = ll + RTMB::dnorm(move_devs[1,r,rr,y,1,a,s], 0, exp(PE_pars[1,r,1,a,s]), TRUE)
            } # end s loop
          } # end a loop
        } # end y loop
      } # iid_y_a_s

      if(PE_model == 5) {
        for(y in 1:n_yrs) {
          for(a in age_start:n_ages) {
            for(s in 1:n_sexes) {
              for(seas in 1:n_seas) {
                ll = ll + RTMB::dnorm(move_devs[1,r,rr,y,seas,a,s], 0, exp(PE_pars[1,r,seas,a,s]), TRUE)
              } # end seas loop
            } # end s loop
          } # end a loop
        } # end y loop
      } # iid_y_seas_a_s

      if(PE_model == 6) {
        for(p in 1:n_pop) {
          for(y in 1:n_yrs) {
            ll = ll + RTMB::dnorm(move_devs[p,r,rr,y,1,1,1], 0, exp(PE_pars[p,r,1,1,1]), TRUE)
          } # end y loop
        } # end p loop
      } # iid_p_y

      if(PE_model == 7) {
        for(p in 1:n_pop) {
          for(a in age_start:n_ages) {
            ll = ll + RTMB::dnorm(move_devs[p,r,rr,1,1,a,1], 0, exp(PE_pars[p,r,1,a,1]), TRUE)
          } # end a loop
        } # end p loop
      } # iid_p_a

      if(PE_model == 8) {
        for(p in 1:n_pop) {
          for(y in 1:n_yrs) {
            for(a in age_start:n_ages) {
              ll = ll + RTMB::dnorm(move_devs[p,r,rr,y,1,a,1], 0, exp(PE_pars[p,r,1,a,1]), TRUE)
            } # end a loop
          } # end y loop
        } # end p loop
      } # iid_p_y_a

      if(PE_model == 9) {
        for(p in 1:n_pop) {
          for(y in 1:n_yrs) {
            for(a in age_start:n_ages) {
              for(s in 1:n_sexes) {
                ll = ll + RTMB::dnorm(move_devs[p,r,rr,y,1,a,s], 0, exp(PE_pars[p,r,1,a,s]), TRUE)
              } # end s loop
            } # end a loop
          } # end y loop
        } # end p loop
      } # iid_p_y_a_s

      if(PE_model == 10) {
        for(p in 1:n_pop) {
          for(y in 1:n_yrs) {
            for(a in age_start:n_ages) {
              for(s in 1:n_sexes) {
                for(seas in 1:n_seas) {
                  ll = ll + RTMB::dnorm(move_devs[p,r,rr,y,seas,a,s], 0, exp(PE_pars[p,r,seas,a,s]), TRUE)
                } # end seas loop
              } # end s loop
            } # end a loop
          } # end y loop
        } # end p loop
      } # iid_p_y_seas_a_s

    } # end r loop
  } # end rr loop

  return(ll)
}
