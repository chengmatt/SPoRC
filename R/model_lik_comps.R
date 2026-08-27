# Stage 2 of 3: objective function
#
# Composition likelihood assembly for every fishery and survey composition
# stream. Get_Comp_Likelihoods evaluates the likelihood itself; pack_comp_osa and
# eval_comp_osa carry the one step ahead residual bookkeeping, which has to be
# built during the objective evaluation and cannot be reconstructed afterwards.

#' Composition Data Likelihood
#'
#' Computes the negative log-likelihood contribution for composition data
#' (age or length) for a single year and fleet. The function supports multiple
#' composition parameterizations and likelihood families commonly used in
#' stock assessment models, including multinomial, Dirichlet-multinomial,
#' and logistic-normal likelihoods.
#'
#' Expected and observed compositions are provided as arrays indexed by
#' region, composition bin (age or length), and sex. The function optionally
#' applies ageing error, aggregates compositions depending on the
#' parameterization type, and evaluates the likelihood for each region
#' and/or sex.
#'
#' @param Exp Expected composition values (e.g., predicted catch-at-age or
#'   survey age compositions), structured as an array indexed by
#'   \eqn{[region \times model\_bins \times sex]}.
#' @param Obs Observed composition counts indexed by
#'   \eqn{[region \times observed\_bins \times sex]}.
#' @param ISS Input sample size for the composition data,
#'   indexed by \eqn{[region \times sex]}.
#' @param Wt_Mltnml Multinomial weighting applied to the effective sample
#'   size, indexed by \eqn{[region \times sex]}.
#' @param ln_theta_agg Log overdispersion parameter used when compositions
#'   are aggregated (\code{Comp_Type = 0}).
#' @param ln_theta Log overdispersion parameters used for
#'   Dirichlet-multinomial or logistic-normal likelihoods, indexed by
#'   \eqn{[region \times sex]}.
#' @param LN_corr_pars Logistic-normal correlation parameters used for
#'   correlated logistic-normal likelihoods, dimensioned by
#'   \eqn{[region \times sex \times parameters]}.
#' @param LN_corr_pars_agg Logistic-normal correlation parameters used
#'   when compositions are aggregated.
#' @param Comp_Type Integer specifying the composition parameterization:
#'   \itemize{
#'     \item \code{0}: Aggregated compositions across sexes and regions.
#'     \item \code{1}: Compositions split by sex and region (no implicit
#'     sex or region ratio information).
#'     \item \code{2}: Joint compositions across sexes but split by region
#'     (implicit sex ratio information).
#'   }
#' @param Likelihood_Type Integer specifying the likelihood family:
#'   \itemize{
#'     \item \code{0}: Multinomial.
#'     \item \code{1}: Dirichlet-multinomial.
#'     \item \code{2}: Logistic-normal with independent bins.
#'     \item \code{3}: Logistic-normal with AR(1) correlation across bins.
#'     \item \code{4}: Logistic-normal with AR(1) correlation across bins
#'     and constant correlation across sexes.
#'   }
#' @param n_regions Number of regions modeled.
#' @param n_model_bins Number of composition bins used internally in the model.
#' @param n_obs_bins Number of observed composition bins.
#' @param n_sexes Number of sexes modeled.
#' @param age_or_len Indicator for composition type:
#'   \itemize{
#'     \item \code{0}: Age compositions.
#'     \item \code{1}: Length compositions.
#'   }
#' @param AgeingError Ageing error matrix used to map model age bins to
#'   observed age bins. For length compositions, either \code{NA} (the model
#'   and observed bins coincide) or a matrix \code{[n_model_bins x n_obs_bins]}
#'   mapping the model's length bins onto the observed ones, the same way, when
#'   the compositions are recorded on coarser bins than the model carries.
#' @param use Integer vector indicating which regions have observations
#'   (\code{1} = use data, \code{0} = ignore).
#' @param comp_bins Integer vector of bins the composition is fitted over, or
#'   \code{NULL} (default) for all of them. Both the observed and expected
#'   compositions are restricted to these bins and renormalized within them, so
#'   bins outside the range are left out of the likelihood rather than being
#'   forced to be explained. Indices refer to observed bins, that is after any
#'   ageing error or length bin map has mapped model bins onto observed ones.
#'   The restriction applies to every composition type: for the sex-joint comps
#'   (\code{Comp_Type = 2}) the named bins are dropped from each sex's block of
#'   the joint stack, so the sex ratio the joint comps carry is the ratio within
#'   the fitted bins. Logistic-normal covariances are built over all observed
#'   bins and then cut down to the fitted ones, so a gap in \code{comp_bins}
#'   still counts towards the AR1 lag between the bins on either side of it.
#' @param addtocomp Small constant added to compositions to avoid numerical
#'   issues when zeros are present.
#' @param comp_const_obs Integer (0 or 1). Whether \code{addtocomp} is added to
#'   the observed proportions that weight the multinomial and enter the
#'   Dirichlet-multinomial. For the Dirichlet-multinomial a constant on the
#'   observed side is not neutral: a bin with no observed and no expected mass
#'   contributes \eqn{\log(\theta/(1+\theta))} regardless of the constant's size,
#'   so compositions with many structurally empty bins (conditional
#'   age-at-length above all) potentially bias \eqn{\theta} upward under \code{1}. Use
#'   \code{0} there. It is also added inside the
#'   logarithms. \code{1} (default) is the unbiased choice: the stationary point
#'   of the likelihood is exactly \code{p = obs}. \code{0} weights by the raw
#'   observed proportions. 
#' @keywords internal
Get_Comp_Likelihoods = function(Exp,
                                Obs,
                                ISS,
                                Wt_Mltnml,
                                ln_theta_agg,
                                ln_theta,
                                LN_corr_pars = 0,
                                LN_corr_pars_agg = 0,
                                Comp_Type,
                                Likelihood_Type,
                                n_regions,
                                n_model_bins,
                                n_obs_bins,
                                n_sexes,
                                age_or_len,
                                AgeingError,
                                use,
                                addtocomp,
                                comp_bins = NULL,
                                comp_const_obs = 1
                                ) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  rho_trans = function(x) 2/(1+ exp(-2 * x)) - 1 # constraint between -1 and 1
  comp_nLL = array(0, dim = c(n_regions, n_sexes)) # initialize nLL here
  const = addtocomp # small constant
  # Filter expectation and observations to regions that have observations
  n_regions_obs_use = sum(use == 1) # get number of regions that have observations

  # Making sure things are correctly formatted (and regions are not dropped)
  Obs = array(Obs, dim = c(n_regions, n_obs_bins, n_sexes))
  dim(Obs) = c(n_regions, n_obs_bins, n_sexes)

  Exp = array(Exp, dim = c(n_regions, n_model_bins, n_sexes)) # using n_obs_bins, because non-square ageing error matrix will collapse to n_obs_bins
  ISS = array(ISS, dim = c(n_regions, n_sexes))
  Wt_Mltnml = array(Wt_Mltnml, dim = c(n_regions, n_sexes))
  ln_theta = array(ln_theta, dim = c(n_regions, n_sexes))
  LN_corr_pars = array(LN_corr_pars, dim = c(n_regions, n_sexes, 3))

  # Filter every per-region input down to the regions that have observations, so the
  # loop index r lines up across Obs/Exp/ISS/Wt_Mltnml/ln_theta/LN_corr_pars
  used = which(use == 1)
  Obs = Obs[used,,,drop = FALSE]
  Exp = Exp[used,,,drop = FALSE]
  ISS = ISS[used,,drop = FALSE]
  Wt_Mltnml = Wt_Mltnml[used,,drop = FALSE]
  ln_theta = ln_theta[used,,drop = FALSE]
  LN_corr_pars = LN_corr_pars[used,,,drop = FALSE]

  # Bin restriction here
  fit_bins = if(is.null(comp_bins)) seq_len(n_obs_bins) else comp_bins
  n_fit_bins = length(fit_bins)
  fit_bins_joint = as.vector(outer(fit_bins, (seq_len(n_sexes) - 1) * n_obs_bins, "+"))
  # Comparing against the full run of bins rather than just counting them, so a
  # reordered comp_bins is honoured rather than passing as unrestricted
  restrict = !identical(as.integer(fit_bins), seq_len(n_obs_bins))
  if(restrict) Obs = Obs[,fit_bins,,drop = FALSE]

  # Aggregated comps by sex and region
  if(Comp_Type == 0) {

    # Nothing in the fitted bins means nothing to fit, and normalizing would give
    # NaN. Same guard the split and joint comps carry.
    if(!any(is.finite(Obs[1,,1])) || sum(Obs[1,,1], na.rm = TRUE) == 0) return(comp_nLL)

    tmp_Exp = matrix(rowSums(matrix(Exp, nrow = n_model_bins)) / (n_sexes * n_regions), nrow = 1) # aggregate
    tmp_Exp = tmp_Exp / sum(tmp_Exp) # normalize

    # Expected age bins get collapsed to observed age bins if ageing error is
    # non-square; length bins likewise when a bin map is handed in
    if(age_or_len == 0 || is.matrix(AgeingError)) {
      tmp_Exp = tmp_Exp %*% AgeingError # apply ageing error, or the length bin map
      tmp_Exp = as.vector((tmp_Exp) / sum(tmp_Exp)) # renormalize
    } else tmp_Exp = as.vector((tmp_Exp) / sum(tmp_Exp)) # renormalize (lengths)

    # Restrict the expectation to the bins being fit and renormalize within them
    if(restrict) tmp_Exp = tmp_Exp[fit_bins] / sum(tmp_Exp[fit_bins])

    # Multinomial likelihood
    if(Likelihood_Type == 0) { # Note that this indexes 1 because it's only a single sex and single region
      tmp_Obs = (Obs[1,,1]) / sum(Obs[1,,1]) # Normalize observed values
      ESS = ISS[1,1] * Wt_Mltnml[1,1] # Effective sample size
      obs_w = if(comp_const_obs == 1) tmp_Obs + const else tmp_Obs # add composition constant to observed or not
      comp_nLL[1,1] = -1 * ESS * sum((obs_w * log(tmp_Exp + const))) # ADMB multinomial likelihood
      comp_nLL[1,1] = comp_nLL[1,1] - -1 * ESS * sum((obs_w * log(tmp_Obs + const))) # Multinomial offset (subtract offset from actual likelihood)
    } # end if multinomial likelihood

    if(Likelihood_Type == 1) {
      tmp_Obs = (Obs[1,,1]) / sum(Obs[1,,1]) # Normalize observed values
      comp_nLL[1,1] = -1 * ddirmult(obs = tmp_Obs, pred = tmp_Exp, Ntotal = ISS[1,1], ln_theta = ln_theta_agg, TRUE) # Dirichlet Multinomial likelihood
    } # end if dirichlet multinomial

    if(Likelihood_Type == 2) {

      # Dealing with zeros
      tmp_Obs = Obs[1,,1] / sum(Obs[1,,1])
      zeros = which(tmp_Obs == 0)

      # Construct Sigma
      Sigma = diag(rep(1/exp(ln_theta_agg)^2, length(tmp_Obs)))

      if(length(zeros) > 0) {
        # Remove zeros and renormalize
        tmp_Obs = tmp_Obs[-zeros]
        tmp_Exp = tmp_Exp[-zeros]
        tmp_Obs = tmp_Obs / sum(tmp_Obs)
        tmp_Exp = tmp_Exp / sum(tmp_Exp)
        # Adjust Sigma
        Sigma = Sigma[-zeros, -zeros]
      }

      Sigma = Sigma[-nrow(Sigma), -ncol(Sigma)] # remove last row and column
      comp_nLL[1,1] = -1 * dlogistnormal(obs = tmp_Obs, pred = tmp_Exp, Sigma = Sigma, TRUE) # Logistic Normal likelihood (iid)
    } # end if logistic normal (iid)

    if(Likelihood_Type == 3) {
      # Dealing with zeros
      tmp_Obs = Obs[1,,1] / sum(Obs[1,,1])
      zeros = which(tmp_Obs == 0)

      # Construct Sigma
      LN_corr_b = rho_trans(LN_corr_pars_agg) # correlation by age / length
      Sigma =  get_AR1_CorrMat(n_obs_bins, LN_corr_b)
      Sigma = Sigma * (exp(ln_theta_agg)^2 / (1 - LN_corr_b^2))
      if(restrict) Sigma = Sigma[fit_bins, fit_bins] # cut to the bins being fit, lags measured over the full range

      if(length(zeros) > 0) {
        # Remove zeros and renormalize
        tmp_Obs = tmp_Obs[-zeros]
        tmp_Exp = tmp_Exp[-zeros]
        tmp_Obs = tmp_Obs / sum(tmp_Obs)
        tmp_Exp = tmp_Exp / sum(tmp_Exp)
        # Adjust Sigma
        Sigma = Sigma[-zeros, -zeros]
      }

      Sigma = Sigma[-nrow(Sigma), -ncol(Sigma)] # remove last row and column
      comp_nLL[1,1] = -1 * dlogistnormal(obs = tmp_Obs, pred = tmp_Exp, Sigma = Sigma, TRUE) # Logistic Normal likelihood (1dar1)
    } # end if logistic normal (1dar1)

  } # end if aggregated comps across sexes and regions

  # 'Split' comps by sex and region (no implicit sex ratio information)
  if(Comp_Type == 1) {
    for(s in 1:n_sexes) {
      for(r in 1:n_regions_obs_use) {

        if(!any(is.finite(Obs[r,,s])) || sum(Obs[r,,s], na.rm = TRUE) == 0) next # skipping zeros since nromalizing would cause Inf
        # Expected Values
        if(age_or_len == 0 || is.matrix(AgeingError)) {
          tmp_Exp = ((Exp[r,,s]) / sum(Exp[r,,s])) %*% AgeingError # Normalize temporary variable (ages), or map the length bins
          tmp_Exp = tmp_Exp / sum(tmp_Exp) # renormalize
        } else tmp_Exp = (Exp[r,,s]) / sum(Exp[r,,s]) # Normalize lengths

        # Restrict the expectation to the bins being fit and renormalize within them
        if(restrict) tmp_Exp = tmp_Exp[fit_bins] / sum(tmp_Exp[fit_bins])

        # Multinomial likelihood
        if(Likelihood_Type == 0) {
          tmp_Obs = (Obs[r,,s]) / sum(Obs[r,,s]) # Normalize observed temporary variable
          ESS = ISS[r,s] * Wt_Mltnml[r,s] # Effective sample size
          obs_w = if(comp_const_obs == 1) tmp_Obs + const else tmp_Obs
          comp_nLL[r,s] = -1 * ESS * sum((obs_w * log(tmp_Exp + const))) # ADMB multinomial likelihood
          comp_nLL[r,s] = comp_nLL[r,s] - -1 * ESS * sum((obs_w * log(tmp_Obs + const))) # Multinomial offset (subtract offset from actual likelihood)
        } # end if multinomial likelihood

        if(Likelihood_Type == 1) {
          tmp_Obs = if(comp_const_obs == 1) (Obs[r,,s] + const) / sum(Obs[r,,s] + const) else Obs[r,,s] / sum(Obs[r,,s]) # Normalize observed temporary variable
          tmp_Exp = (tmp_Exp + const) / sum(tmp_Exp + const) # Normalize expected temporary variable (lgamma can't take 0s)
          comp_nLL[r,s] = -1 * ddirmult(tmp_Obs, tmp_Exp, ISS[r,s], ln_theta[r,s], TRUE) # Dirichlet Multinomial likelihood
        } # end if dirichlet multinomial

        if(Likelihood_Type == 2) {

          # Dealing with zeros
          tmp_Obs = Obs[r,,s] / sum(Obs[r,,s])
          zeros = which(tmp_Obs == 0)

          # Construct Sigma
          Sigma = diag(rep(exp(ln_theta[r,s])^2, length(tmp_Obs)))

          if(length(zeros) > 0) {
            # Remove zeros and renormalize
            tmp_Obs = tmp_Obs[-zeros]
            tmp_Exp = tmp_Exp[-zeros]
            tmp_Obs = tmp_Obs / sum(tmp_Obs)
            tmp_Exp = tmp_Exp / sum(tmp_Exp)
            # Adjust Sigma
            Sigma = Sigma[-zeros, -zeros]
          }

          Sigma = Sigma[-nrow(Sigma), -ncol(Sigma)] # remove last row and column
          comp_nLL[r,s] = -1 * dlogistnormal(obs = tmp_Obs, pred = tmp_Exp, Sigma = Sigma, TRUE) # Logistic Normal likelihood (iid)
        } # end if logistic normal

        if(Likelihood_Type == 3) {
          # Dealing with zeros
          tmp_Obs = Obs[r,,s] / sum(Obs[r,,s])
          zeros = which(tmp_Obs == 0)

          # Construct Sigma
          LN_corr_b = rho_trans(LN_corr_pars[r,s,1]) # correlation by age / length
          Sigma =  get_AR1_CorrMat(n_obs_bins, LN_corr_b)
          Sigma = Sigma * (exp(ln_theta[r,s])^2 / (1 - LN_corr_b^2))
          if(restrict) Sigma = Sigma[fit_bins, fit_bins] # cut to the bins being fit, lags measured over the full range

          if(length(zeros) > 0) {
            # Remove zeros and renormalize
            tmp_Obs = tmp_Obs[-zeros]
            tmp_Exp = tmp_Exp[-zeros]
            tmp_Obs = tmp_Obs / sum(tmp_Obs)
            tmp_Exp = tmp_Exp / sum(tmp_Exp)
            # Adjust Sigma
            Sigma = Sigma[-zeros, -zeros]
          }

          Sigma = Sigma[-nrow(Sigma), -ncol(Sigma)] # remove last row and column
          comp_nLL[r,s] = -1 * dlogistnormal(obs = tmp_Obs, pred = tmp_Exp, Sigma = Sigma, TRUE) # Logistic Normal likelihood (1dar1)
        } # end if logistic normal (1dar1)

      } # end r loop
    } # end s loop
  } # end if 'Split' comps by sex and region

  # Joint by sex, Split by region
  if(Comp_Type == 2) {
    for(r in 1:n_regions_obs_use) {

      # A region whose fitted block holds no fish contributes nothing, and
      # normalizing it would give Inf. Matches the guard the split comps carry,
      # and matters more once a bin restriction can empty a block that the full
      # composition filled.
      if(!any(is.finite(Obs[r,,])) || sum(Obs[r,,], na.rm = TRUE) == 0) next

      # Expected values
      if(age_or_len == 0 || is.matrix(AgeingError)) { # if ages, or lengths with a bin map
        tmp_Exp = t(as.vector((Exp[r,,])/ sum(Exp[r,,]))) %*% kronecker(diag(n_sexes), AgeingError) # apply ageing error, or the length bin map
        tmp_Exp = as.vector((tmp_Exp) / sum(tmp_Exp)) # renormalize to make sure sum to 1
      } else tmp_Exp = as.vector((Exp[r,,]) / sum((Exp[r,,]))) # Normalize temporary variable (lengths)

      # Restrict the expectation to the bins being fit and renormalize within them
      if(restrict) tmp_Exp = tmp_Exp[fit_bins_joint] / sum(tmp_Exp[fit_bins_joint])

      # Multinomial likelihood
      if(Likelihood_Type == 0) { # Indexing by r for a given region since it's 'Split' by region and 1 for sex since it's 'Joint' for sex
        tmp_Obs = as.vector((Obs[r,,]) / sum(Obs[r,,])) # Normalize observed temporary variable
        ESS = ISS[r,1] * Wt_Mltnml[r,1] # Effective sample size
        obs_w = if(comp_const_obs == 1) tmp_Obs + const else tmp_Obs
        comp_nLL[r,1] = -1 * ESS * sum((obs_w * log(tmp_Exp + const))) # ADMB multinomial likelihood
        comp_nLL[r,1] = comp_nLL[r,1] - -1 * ESS * sum((obs_w * log(tmp_Obs + const))) # Multinomial offset (subtract offset from actual likelihood)
      } # end if multinomial likelihood

      if(Likelihood_Type == 1) {
        # see the split by sex case for why the constant is optional on the observed side
        tmp_Obs = if(comp_const_obs == 1) as.vector((Obs[r,,] + const) / sum(Obs[r,,] + const)) else as.vector(Obs[r,,] / sum(Obs[r,,])) # Normalize observed temporary variable
        tmp_Exp = (tmp_Exp + const) / sum(tmp_Exp + const) # normalize temporary expected variable
        comp_nLL[r,1] = -1 * ddirmult(tmp_Obs, tmp_Exp, ISS[r,1], ln_theta[r,1], TRUE) # Dirichlet Multinomial likelihood
      } # end if dirichlet multinomial

      if(Likelihood_Type == 2) {
        # Dealing with zeros
        tmp_Obs = Obs[r,,] / sum(Obs[r,,])
        zeros = which(tmp_Obs == 0)

        # Construct Sigma
        Sigma = diag(rep(exp(ln_theta[r,1])^2, length(tmp_Obs)))

        if(length(zeros) > 0) {
          # Remove zeros and renormalize
          tmp_Obs = tmp_Obs[-zeros]
          tmp_Exp = tmp_Exp[-zeros]
          tmp_Obs = tmp_Obs / sum(tmp_Obs)
          tmp_Exp = tmp_Exp / sum(tmp_Exp)
          # Adjust Sigma
          Sigma = Sigma[-zeros, -zeros]
        }

        Sigma = Sigma[-nrow(Sigma), -ncol(Sigma)] # remove last row and column

        comp_nLL[r,1] = -1 * dlogistnormal(obs = tmp_Obs, pred = tmp_Exp, Sigma = Sigma, TRUE) # Logistic Normal likelihood (iid)
      } # end if logistic normal (iid)

      if(Likelihood_Type == 3) {

        # Dealing with zeros
        tmp_Obs = Obs[r,,] / sum(Obs[r,,])
        zeros = which(tmp_Obs == 0)

        # Construct Sigma
        LN_corr_b = rho_trans(LN_corr_pars[r,1,1]) # correlation by age / length
        Sigma =  get_AR1_CorrMat(n_obs_bins * n_sexes, LN_corr_b)
        Sigma = Sigma * (exp(ln_theta[r,1])^2 / (1 - LN_corr_b^2))
        if(restrict) Sigma = Sigma[fit_bins_joint, fit_bins_joint] # cut to the bins being fit, lags measured over the full range

        if(length(zeros) > 0) {
          # Remove zeros and renormalize
          tmp_Obs = tmp_Obs[-zeros]
          tmp_Exp = tmp_Exp[-zeros]
          tmp_Obs = tmp_Obs / sum(tmp_Obs)
          tmp_Exp = tmp_Exp / sum(tmp_Exp)
          # Adjust Sigma
          Sigma = Sigma[-zeros, -zeros]
        }

        Sigma = Sigma[-nrow(Sigma), -ncol(Sigma)] # remove last row and column

        comp_nLL[r,1] = -1 * dlogistnormal(obs = tmp_Obs, pred = tmp_Exp, Sigma = Sigma, TRUE) # Logistic Normal likelihood (1dar1)
      } # end if logistic normal (1dar1)

      if(Likelihood_Type == 4) {

        # Dealing with zeros
        tmp_Obs = Obs[r,,] / sum(Obs[r,,])
        zeros = which(tmp_Obs == 0)

        # Construct Sigma
        LN_corr_b = rho_trans(LN_corr_pars[r,1,1])
        LN_corr_s = rho_trans(LN_corr_pars[r,1,2])
        mat1 = get_Constant_CorrMat(n_sexes, LN_corr_s)
        mat2 = get_AR1_CorrMat(n_obs_bins, LN_corr_b)
        Sigma =  Matrix::kronecker(mat1, mat2)  * (exp(ln_theta[r,1])^2 / (1 - LN_corr_s^2) / (1 - LN_corr_b^2))
        if(restrict) Sigma = Sigma[fit_bins_joint, fit_bins_joint] # cut to the bins being fit, lags measured over the full range

        if(length(zeros) > 0) {
          # Remove zeros and renormalize
          tmp_Obs = tmp_Obs[-zeros]
          tmp_Exp = tmp_Exp[-zeros]
          tmp_Obs = tmp_Obs / sum(tmp_Obs)
          tmp_Exp = tmp_Exp / sum(tmp_Exp)
          # Adjust Sigma
          Sigma = Sigma[-zeros, -zeros]
        }

        Sigma = Sigma[-nrow(Sigma), -ncol(Sigma)] # remove last bins

        # likelihood
        comp_nLL[r,1] = -1 * dlogistnormal(obs = tmp_Obs, pred = tmp_Exp, Sigma = Sigma, TRUE) # Logistic Normal likelihood (1dar1 by age, constant corr by sex)
      }

    } # end r loop
  } # end if 'Joint' comps by sex, but 'Split' by region

  # Input the split-comp results (Comp_Type 1/2) from the compacted region slots
  # 1:n_regions_obs_use, that the loop above wrote to, back to their true region slots.
  if(Comp_Type %in% c(1, 2) && n_regions_obs_use < n_regions) {
    reordered = RTMB::AD(array(0, dim = c(n_regions, n_sexes)))
    dim(reordered) = c(n_regions, n_sexes)
    reordered[used, ] = comp_nLL[seq_len(n_regions_obs_use), ]
    comp_nLL = reordered
  }

  return(comp_nLL) # return negative log likelihood
  
} # end function

#' Composition Data Likelihood (OSA variant)
#'
#' Computes multinomial (0), Dirichlet-multinomial (1), and logistic-normal
#' (2 iid, 3 AR1, 4 2D‑AR1) composition likelihoods for one‑step‑ahead (OSA)
#' residuals using \code{RTMB::oneStepPredict}. The function evaluates the
#' likelihood for a single flat tracked OBS vector, respecting the reduced
#' logistic‑normal block lengths used during packing.
#'
#' The tracked \code{Obs} vector is **never reshaped**. All expectation‑side
#' quantities (\code{Exp}, \code{ISS}, \code{ln_theta}, \code{LN_corr_pars},
#' ageing error) are reshaped and filtered by \code{use}, exactly as in fitting.
#'
#' **Logistic‑normal note:** Because \code{RTMB::OBS()} cannot be altered after
#' tracking, the additive‑log‑ratio (ALR) transform of the *observation* is
#' performed in the packer. Thus, \code{Obs[idx]} is **already ALR‑transformed**
#' (last bin dropped). Here we only ALR‑transform the expectation, construct the
#' covariance matrix \code{Sigma} (dropping its last row/column), and evaluate
#' the multivariate normal density.
#'
#' Reduced LN block lengths:
#'   * Comp_Type 0: \code{n_fit_bins - 1}
#'   * Comp_Type 1: \code{n_fit_bins - 1} per region/sex
#'   * Comp_Type 2: \code{n_fit_bins * n_sexes - 1} per region (one joint reference)
#'
#' where \code{n_fit_bins} is the number of bins named by \code{comp_bins},
#' equal to \code{n_obs_bins} when the fleet fits every bin. The packer applies
#' the same restriction before transforming, so the ALR reference is the last
#' fitted bin rather than the last observed one.
#'
#' @param Exp Expected proportions [n_regions × n_model_bins × n_sexes].
#' @param Obs Flat tracked observation vector (already ALR‑transformed for LN).
#' @param ISS Input sample size [n_regions × n_sexes].
#' @param ln_theta Log overdispersion [n_regions × n_sexes].
#' @param ln_theta_agg Log overdispersion scalar for aggregated comps.
#' @param LN_corr_pars LN correlation parameters [n_regions × n_sexes × 3].
#' @param LN_corr_pars_agg LN aggregated correlation scalar(s).
#' @inheritParams Get_Comp_Likelihoods
#' @keywords internal
Get_Comp_Likelihoods_OSA = function(Exp,
                                    Obs,
                                    ISS,
                                    ln_theta,
                                    ln_theta_agg,
                                    LN_corr_pars = 0,
                                    LN_corr_pars_agg = 0,
                                    Comp_Type,
                                    Likelihood_Type,
                                    n_regions,
                                    n_model_bins,
                                    n_obs_bins,
                                    n_sexes,
                                    age_or_len,
                                    AgeingError,
                                    use,
                                    addtocomp,
                                    comp_bins = NULL) {

  "c"   <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  rho_trans = function(x) 2 / (1 + exp(-2 * x)) - 1
  alr_mu    = function(p) log(p[-length(p)]) - log(p[length(p)])

  # Bin restriction, mirroring Get_Comp_Likelihoods
  fit_bins = if(is.null(comp_bins)) seq_len(n_obs_bins) else comp_bins
  n_fit_bins = length(fit_bins)
  fit_bins_joint = as.vector(outer(fit_bins, (seq_len(n_sexes) - 1) * n_obs_bins, "+"))
  restrict = !identical(as.integer(fit_bins), seq_len(n_obs_bins))

  const = addtocomp
  used  = which(use == 1)          # observed region rows
  n_ru  = length(used)             # n_regions_obs_use

  # comp_nLL AD from birth so subassignment keeps class
  comp_nLL = array(0, dim = c(n_regions, n_sexes))
  comp_nLL = RTMB::AD(comp_nLL)
  dim(comp_nLL) = c(n_regions, n_sexes)

  # Exp-side quantities: reshape + filter to observed regions (NOT tracked)
  Exp          = array(Exp,          dim = c(n_regions, n_model_bins, n_sexes))[used, , , drop = FALSE]
  ISS          = array(ISS,          dim = c(n_regions, n_sexes))[used, , drop = FALSE]
  ln_theta     = array(ln_theta,     dim = c(n_regions, n_sexes))[used, , drop = FALSE]
  LN_corr_pars = array(LN_corr_pars, dim = c(n_regions, n_sexes, 3))[used, , , drop = FALSE]

  # Comp_Type 0 (aggregated)
  if(Comp_Type == 0) {

    tmp_Exp = rowSums(matrix(Exp, nrow = n_model_bins)) / (n_sexes * n_ru)
    tmp_Exp = tmp_Exp / sum(tmp_Exp)
    if(age_or_len == 0 || is.matrix(AgeingError)) {
      tmp_Exp = as.vector(matrix(tmp_Exp, nrow = 1) %*% AgeingError)
      tmp_Exp = tmp_Exp / sum(tmp_Exp)
    }
    if(restrict) tmp_Exp = tmp_Exp[fit_bins] / sum(tmp_Exp[fit_bins])
    tmp_Exp = tmp_Exp + const
    tmp_Exp = tmp_Exp / sum(tmp_Exp)

    if(Likelihood_Type == 0) { # Multinomial
      comp_nLL[1,1] = -dmultinom_osa(Obs, tmp_Exp, log = TRUE)
    }
    if(Likelihood_Type == 1) { # Dirichlet-multinomial
      comp_nLL[1,1] = -ddirmult_osa(Obs, tmp_Exp * exp(ln_theta_agg) * ISS[1,1], log = TRUE)
    }
    if(Likelihood_Type %in% c(2,3)) { # Logistic-normal (Obs already ALR)
      if(Likelihood_Type == 2) {
        Sigma = diag(rep(exp(ln_theta_agg)^2, n_fit_bins))
      } else {
        LN_corr_b = rho_trans(LN_corr_pars_agg)
        Sigma = get_AR1_CorrMat(n_obs_bins, LN_corr_b) * (exp(ln_theta_agg)^2 / (1 - LN_corr_b^2))
        if(restrict) Sigma = Sigma[fit_bins, fit_bins]
      }
      Sigma = Sigma[-nrow(Sigma), -ncol(Sigma)]
      mu = alr_mu(tmp_Exp)
      comp_nLL[1,1] = -RTMB::dmvnorm(x = Obs, mu = mu, Sigma = Sigma, log = TRUE)
    }

  } # end Comp_Type 0

  # Comp_Type 1 (split sex & region)
  if(Comp_Type == 1) {
    for(s in 1:n_sexes) {
      for(r in 1:n_ru) {

        if(age_or_len == 0 || is.matrix(AgeingError)) {
          tmp_Exp = as.vector((Exp[r,,s] / sum(Exp[r,,s])) %*% AgeingError)
          tmp_Exp = tmp_Exp / sum(tmp_Exp)
        } else {
          tmp_Exp = Exp[r,,s] / sum(Exp[r,,s])
        }
        if(restrict) tmp_Exp = tmp_Exp[fit_bins] / sum(tmp_Exp[fit_bins])
        tmp_Exp = tmp_Exp + const
        tmp_Exp = tmp_Exp / sum(tmp_Exp)

        if(Likelihood_Type %in% c(0,1)) {
          idx = seq(from = r + (s - 1) * n_ru * n_fit_bins, by = n_ru, length.out = n_fit_bins)
        }
        if(Likelihood_Type %in% c(2,3)) {
          idx = seq(from = r + (s - 1) * n_ru * (n_fit_bins - 1), by = n_ru, length.out = n_fit_bins - 1)
        }
        if(Likelihood_Type == 0) { # Multinomial
          comp_nLL[used[r],s] = -dmultinom_osa(Obs[idx], tmp_Exp, log = TRUE)
        }
        if(Likelihood_Type == 1) { # Dirichlet-multinomial
          comp_nLL[used[r],s] = -ddirmult_osa(Obs[idx], tmp_Exp * exp(ln_theta[r,s]) * ISS[r,s], log = TRUE)
        }
        if(Likelihood_Type %in% c(2,3)) { # Logistic-normal (Obs already ALR)
          if(Likelihood_Type == 2) {
            Sigma = diag(rep(exp(ln_theta[r,s])^2, n_fit_bins))
          } else {
            LN_corr_b = rho_trans(LN_corr_pars[r,s,1])
            Sigma = get_AR1_CorrMat(n_obs_bins, LN_corr_b) * (exp(ln_theta[r,s])^2 / (1 - LN_corr_b^2))
            if(restrict) Sigma = Sigma[fit_bins, fit_bins]
          }
          Sigma = Sigma[-nrow(Sigma), -ncol(Sigma)]
          mu = alr_mu(tmp_Exp)
          comp_nLL[used[r],s] = -RTMB::dmvnorm(x = Obs[idx], mu = mu, Sigma = Sigma, log = TRUE)
        }

      } # r
    } # s
  } # end Comp_Type 1

  # Comp_Type 2 (joint sex, split region)
  if(Comp_Type == 2) {
    for(r in 1:n_ru) {

      if(age_or_len == 0 || is.matrix(AgeingError)) {
        tmp_Exp = as.vector(t(as.vector(Exp[r,,] / sum(Exp[r,,]))) %*% kronecker(diag(n_sexes), AgeingError))
        tmp_Exp = tmp_Exp / sum(tmp_Exp)
      } else {
        tmp_Exp = as.vector(Exp[r,,] / sum(Exp[r,,]))
      }
      if(restrict) tmp_Exp = tmp_Exp[fit_bins_joint] / sum(tmp_Exp[fit_bins_joint])
      tmp_Exp = tmp_Exp + const
      tmp_Exp = tmp_Exp / sum(tmp_Exp)

      if(Likelihood_Type %in% c(0,1)) {
        idx = seq(from = r, by = n_ru, length.out = n_fit_bins * n_sexes)
      }
      if(Likelihood_Type %in% c(2,3,4)) {
        # joint drops ONE reference for the whole [bin x sex] stack
        Lred = n_fit_bins * n_sexes - 1
        idx  = seq(from = r, by = n_ru, length.out = Lred)

        if(Likelihood_Type == 2) { # iid
          Sigma = diag(rep(exp(ln_theta[r,1])^2, n_fit_bins * n_sexes))
        }
        if(Likelihood_Type == 3) { # 1D AR1
          LN_corr_b = rho_trans(LN_corr_pars[r,1,1])
          Sigma = get_AR1_CorrMat(n_obs_bins * n_sexes, LN_corr_b)
          Sigma = Sigma * (exp(ln_theta[r,1])^2 / (1 - LN_corr_b^2))
          if(restrict) Sigma = Sigma[fit_bins_joint, fit_bins_joint]
        }
        if(Likelihood_Type == 4) { # 2D AR1 (age/len x sex)
          LN_corr_b = rho_trans(LN_corr_pars[r,1,1])
          LN_corr_s = rho_trans(LN_corr_pars[r,1,2])
          mat1 = get_Constant_CorrMat(n_sexes, LN_corr_s)
          mat2 = get_AR1_CorrMat(n_obs_bins, LN_corr_b)
          Sigma = Matrix::kronecker(mat1, mat2) * (exp(ln_theta[r,1])^2 / (1 - LN_corr_s^2) / (1 - LN_corr_b^2))
          if(restrict) Sigma = Sigma[fit_bins_joint, fit_bins_joint]
        }
        Sigma = Sigma[-nrow(Sigma), -ncol(Sigma)]  # drop last row/col -> Lred x Lred
        mu = alr_mu(tmp_Exp)                        # ALR of expectation -> length Lred
      }

      if(Likelihood_Type == 0) { # Multinomial
        comp_nLL[used[r],1] = -dmultinom_osa(Obs[idx], tmp_Exp, log = TRUE)
      }
      if(Likelihood_Type == 1) { # Dirichlet-multinomial
        comp_nLL[used[r],1] = -ddirmult_osa(Obs[idx], tmp_Exp * exp(ln_theta[r,1]) * ISS[r,1], log = TRUE)
      }
      if(Likelihood_Type %in% c(2,3,4)) { # Logistic-normal (Obs already ALR)
        comp_nLL[used[r],1] = -RTMB::dmvnorm(x = Obs[idx], mu = mu, Sigma = Sigma, log = TRUE)
      }

    } # r
  } # end Comp_Type 2

  return(comp_nLL)
}

#' Pack observed composition data into a single flat OBS vector (OSA)
#'
#' Produces the flat tracked OBS vector required by \code{RTMB::oneStepPredict}.
#' Population is the outermost dimension, so the entire result is one continuous
#' vector with a single pointer.
#'
#' Discrete families (LikeType 0,1):
#' \itemize{
#'   \item Multinomial (0): counts = round(prop x ISS x Wt)
#'   \item Dirichlet-multinomial (1): counts = round(prop x ISS)
#' }
#'
#' Continuous families (LikeType 2,3,4): logistic-normal
#' The ALR transform is performed here, because the tracked OBS vector
#' cannot be modified later. Proportions receive \code{+addtocomp}, are
#' renormalized, then transformed to \code{log(p_k / p_K)} for k = 1..K-1.
#' The last bin is the ALR reference and is dropped:
#' \itemize{
#'   \item Comp_Type 0: length = \code{n_obs_bins - 1}
#'   \item Comp_Type 1: length = \code{n_ru x (n_obs_bins - 1) x n_sexes}
#'   \item Comp_Type 2: joint ALR of the full [bin x sex] stack ->
#'         length = \code{n_obs_bins x n_sexes - 1}
#'         (Joint drops one reference for the whole stack -> length \code{n_obs_bins * n_sexes - 1})
#' }
#'
#' The resulting vector is ordered region-fastest so that the likelihood
#' evaluator can use simple strided indexing.
#'
#' @param ObsArr Observed proportions or counts.
#' @param ISSArr Input sample sizes.
#' @param WtArr Optional weighting for multinomial.
#' @param UseArr Region-use flags.
#' @param TypeMat Composition type matrix (0,1,2).
#' @param LikeTypeVec Likelihood type per fleet.
#' @param n_yrs Number of model years.
#' @param n_seas Number of seasons per year.
#' @param n_fleets Total number of fishing fleets.
#' @param n_sexes Number of biological sexes.
#' @param addtocomp Small constant added to proportions before normalization.
#' @param family Character string specifying the likelihood type, either "discrete" or "continuous".
#' @param pop Logical; if TRUE, the population dimension is treated as the outermost layer.
#' @param n_pop Number of population structures or pools.
#' @param BinsArr Optional \code{[n_obs_bins x n_fleets]} 0/1 array naming the
#'   observed bins each fleet is fitted over, or \code{NULL} (default) for all
#'   bins. Restricted fleets pack a shorter block, and \code{eval_comp_osa} must
#'   be handed the same array so its strides stay in step with the packer.
#' @param return_labels Logical; if TRUE, also builds a per-element label
#'   data.frame identifying the origin (pop, region, year, season, fleet, sex,
#'   bin, comp_type, likelihood_type, family, last_in_group) of every entry in
#'   the tracked vector, in the same order. Intended for post-hoc relabeling of
#'   \code{TMB::oneStepPredict()} residuals (see [get_osa()]); left \code{FALSE}
#'   (default) inside the model itself to avoid the extra bookkeeping cost.
#'
#' @return If \code{return_labels = FALSE} (default): flat OBS vector, or
#'   \code{NULL} if no fleet of this family is present (unchanged behavior).
#'   If \code{return_labels = TRUE}: a list with elements \code{vec} (the flat
#'   OBS vector) and \code{labels} (a data.frame with one row per element of
#'   \code{vec}), or \code{NULL} if no fleet of this family is present.
#' @keywords internal
pack_comp_osa = function(ObsArr, ISSArr, WtArr, UseArr, TypeMat, LikeTypeVec,
                         n_yrs, n_seas, n_fleets, n_sexes, addtocomp,
                         family = "discrete", pop = FALSE, n_pop = 1,
                         return_labels = FALSE, BinsArr = NULL) {

  "c"   <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  n_obs_bins = if(pop) dim(ObsArr)[5] else dim(ObsArr)[4]

  # Observed bins a fleet is actually fitted over. Everything downstream is sized
  # on these, so a restricted fleet packs a shorter block and its ALR reference
  # becomes the last fitted bin rather than the last observed one.
  # A fleet column is only meaningful against the bins the observations carry.
  # A mismatched row count would desynchronize the packer from the evaluator
  # silently, which is the one failure this machinery cannot survive.
  bins_of = function(f) {
    if(is.null(BinsArr)) return(seq_len(n_obs_bins))
    if(nrow(BinsArr) != n_obs_bins) {
      stop("bin selection array has ", nrow(BinsArr), " rows but the observations carry ", n_obs_bins,
           " bins. The *_bins argument must be indexed on the observed bins of the stream it restricts.")
    }
    which(BinsArr[,f] == 1)
  }

  fam_of = function(lt) if(lt %in% c(0,1)) "discrete" else "continuous"

  # additive log-ratio: drop last element as reference
  alr = function(p) log(p[-length(p)]) - log(p[length(p)])

  # build the per-element label rows for one accepted (p,y,f,seas) group,
  # mirroring the exact element order used to build 'g' below
  make_labels = function(used, n_ru, ct, lt, p, y, seas, f, len) {
    fit_bins = bins_of(f)   # true observed bin numbers, so labels stay comparable across fleets
    n_bins   = length(fit_bins)
    if(lt %in% c(0,1)) { # discrete: every fitted bin retained
      if(ct == 0) {
        region = rep(used[1], n_bins); sex = rep(1L, n_bins); bin = fit_bins
        last_in_group = (bin == fit_bins[n_bins])
      } else {
        region = rep(used, times = n_bins * n_sexes)
        bin    = rep(rep(fit_bins, each = n_ru), times = n_sexes)
        sex    = rep(1:n_sexes, each = n_ru * n_bins)
        if(ct == 1) {
          # split by sex: each (region,sex) is its own independent multinomial,
          # so each has its own redundant/determined bin.
          last_in_group = (bin == fit_bins[n_bins])
        } else {
          # ct == 2, joint by sex: the whole [bin x sex] stack per region is
          # one multinomial (matches the joint scaling in the packing step
          # above), so only the single last cell is redundant/determined,
          # not one per sex. 
          last_in_group = (bin == fit_bins[n_bins]) & (sex == n_sexes)
        }
      }
    } else { # continuous (LN): reference bin already ALR-dropped during packing
      if(ct == 0) {
        region = rep(used[1], n_bins - 1); sex = rep(1L, n_bins - 1); bin = fit_bins[-n_bins]
      } else if(ct == 1) {
        region = rep(used, times = (n_bins - 1) * n_sexes)
        bin    = rep(rep(fit_bins[-n_bins], each = n_ru), times = n_sexes)
        sex    = rep(1:n_sexes, each = n_ru * (n_bins - 1))
      } else {
        Lred = n_bins * n_sexes - 1
        region = rep(used, times = Lred)
        pos_in_stack = rep(1:Lred, each = n_ru)
        bin = fit_bins[((pos_in_stack - 1) %% n_bins) + 1]
        sex = ((pos_in_stack - 1) %/% n_bins) + 1
      }
      last_in_group = rep(FALSE, len)
    }
    data.frame(pop = p, region = region, year = y, season = seas, fleet = f,
               sex = sex, bin = bin, comp_type = ct, likelihood_type = lt,
               family = family, last_in_group = last_in_group)
  }

  clean = list()
  label_rows = list()
  for(p in 1:n_pop) {
    for(y in 1:n_yrs) {
      for(f in 1:n_fleets) {
        if(fam_of(LikeTypeVec[f]) != family) next
        for(seas in 1:n_seas) {
          use_vec = if(pop) UseArr[p, , y, seas, f] else UseArr[, y, seas, f]
          if(sum(use_vec) >= 1) {
            ct   = TypeMat[y, f]
            lt   = LikeTypeVec[f]
            used = which(use_vec == 1)
            n_ru = length(used)
            if(pop) {
              obs_slice = ObsArr[p, used, y, seas, , , f, drop = FALSE]
            } else {
              obs_slice = ObsArr[used, y, seas, , , f, drop = FALSE]
            }
            dim(obs_slice) = c(n_ru, n_obs_bins, n_sexes)

            fit_bins = bins_of(f)
            n_bins   = length(fit_bins)
            if(n_bins < n_obs_bins) obs_slice = obs_slice[, fit_bins, , drop = FALSE]

            if(lt %in% c(0,1)) {
              # Discrete: scale to counts
              if(ct == 2) {
                # Joint by sex: the true sample is one joint draw across sexes
                for(rr in 1:n_ru) {
                  r_orig = used[rr]
                  iss = if(pop) ISSArr[p, r_orig, y, seas, 1, f] else ISSArr[r_orig, y, seas, 1, f]
                  wt  = if(pop) WtArr[p, r_orig, y, seas, 1, f]  else WtArr[r_orig, y, seas, 1, f]
                  v   = as.vector(obs_slice[rr, , ])             # bin-fastest-then-sex
                  pr  = (v + addtocomp) / sum(v + addtocomp)
                  if(lt == 0) v = round(pr * iss * wt)  # multinomial
                  if(lt == 1) v = round(pr * iss)       # DM (no Wt)
                  for(s in 1:n_sexes) obs_slice[rr, , s] = v[((s - 1) * n_bins + 1):(s * n_bins)]
                }
              } else {
                # Aggregated (ct 0) / split by region and sex (ct 1): each
                # region/sex cell is its own independent sample.
                for(rr in 1:n_ru) {
                  r_orig = used[rr]
                  for(s in 1:n_sexes) {
                    iss = if(pop) ISSArr[p, r_orig, y, seas, s, f] else ISSArr[r_orig, y, seas, s, f]
                    wt  = if(pop) WtArr[p, r_orig, y, seas, s, f]  else WtArr[r_orig, y, seas, s, f]
                    pr  = (obs_slice[rr, , s] + addtocomp) / sum(obs_slice[rr, , s] + addtocomp)
                    if(lt == 0) obs_slice[rr, , s] = round(pr * iss * wt)  # multinomial
                    if(lt == 1) obs_slice[rr, , s] = round(pr * iss)       # DM (no Wt)
                  }
                }
              }
              g = if(ct == 0) as.vector(obs_slice[1, , 1]) else as.vector(obs_slice)

            } else {
              # Continuous (LN): ALR-transform per block, drop reference bin
              # Layout mirrors the discrete path so the likelihood(just n_bins-1 instead of n_bins):
              #   ct 0: single vector, length n_bins-1
              #   ct 1: [n_ru, n_bins-1, n_sexes] as.vector (region-fastest)
              #   ct 2: per region, ALR of the full [bin x sex] stack (length n_bins*n_sexes-1), region-fastest via seq
              if(ct == 0) {
                pr = (obs_slice[1, , 1] + addtocomp) / sum(obs_slice[1, , 1] + addtocomp)
                g  = alr(as.vector(pr))
              } else if(ct == 1) {
                # ALR each (region,sex) -> [n_ru, n_bins-1, n_sexes], as.vector col-major
                arr = array(0, dim = c(n_ru, n_bins - 1, n_sexes))
                arr = RTMB::AD(arr); dim(arr) = c(n_ru, n_bins - 1, n_sexes)
                for(rr in 1:n_ru) {
                  for(s in 1:n_sexes) {
                    pr = (obs_slice[rr, , s] + addtocomp) / sum(obs_slice[rr, , s] + addtocomp)
                    arr[rr, , s] = alr(as.vector(pr))
                  }
                }
                g = as.vector(arr)          # region-fastest, matches strided idx
              } else {
                # joint: per region ALR of [bin x sex] stack. Store as
                # [n_ru, (n_bins*n_sexes - 1)] and as.vector region-fastest,
                # matching idx = seq(from=r, by=n_ru, length.out=(n_bins-1)*n_sexes)
                Lred = n_bins * n_sexes - 1
                arr = array(0, dim = c(n_ru, Lred))
                arr = RTMB::AD(arr); dim(arr) = c(n_ru, Lred)
                for(rr in 1:n_ru) {
                  v  = as.vector(obs_slice[rr, , ])           # bin-fastest-then-sex
                  pr = (v + addtocomp) / sum(v + addtocomp)
                  arr[rr, ] = alr(pr)
                }
                g = as.vector(arr)          # region-fastest
              }
            }

            clean[[length(clean) + 1]] = g
            if(return_labels) {
              label_rows[[length(label_rows) + 1]] = make_labels(used, n_ru, ct, lt, p, y, seas, f, length(g))
            }
          }
        }
      }
    }
  }
  if(length(clean) == 0) return(NULL)
  vec = do.call(c, clean)
  if(!return_labels) return(vec)
  list(vec = vec, labels = do.call(rbind, label_rows))
}

#' Evaluate OSA composition negative log-likelihood from a flat tracked vector
#'
#' Walks the same group order used by \code{pack_comp_osa()}, keeping the
#' pointer \code{k} synchronized with the packed slice lengths. Evaluates the
#' multinomial, Dirichlet-multinomial, or logistic-normal likelihood for each
#' region/sex/fleet/season block.
#'
#' Slice lengths must match the packer exactly:
#'
#' Discrete (LikeType 0,1):
#' \itemize{
#'   \item Comp_Type 0: \code{n_fit_bins}
#'   \item Comp_Type 1/2: \code{n_ru x n_fit_bins x n_sexes}
#' }
#'
#' Logistic-normal (LikeType 2,3,4):
#' \itemize{
#'   \item Comp_Type 0: \code{n_fit_bins - 1}
#'   \item Comp_Type 1: \code{n_ru x (n_fit_bins - 1) x n_sexes}
#'   \item Comp_Type 2: \code{n_ru x (n_fit_bins x n_sexes - 1)}
#' }
#'
#' \code{n_fit_bins} is the number of bins the fleet is fitted over, taken from
#' \code{BinsArr} and equal to \code{n_obs_bins} when the fleet fits every bin.
#'
#' These reduced lengths reflect that the tracked \code{Obs} vector is already
#' ALR-transformed (the last reference bin is dropped).
#'
#' @param nLL_arr Array receiving negative log-likelihood contributions.
#' @param tracked Flat tracked OBS vector.
#' @param ExpArrFn Function returning expected proportions for (p,y,seas,f).
#' @param UseArr Region-use flags.
#' @param TypeMat Composition type matrix.
#' @param LikeTypeVec Likelihood type per fleet.
#' @param ISSArr Input sample sizes.
#' @param lnThetaArr Log overdispersion.
#' @param lnThetaAggVec Aggregated log overdispersion.
#' @param LNcorrArr LN correlation parameters.
#' @param LNcorrAggVec Aggregated LN correlation parameters.
#' @param n_regions Total number of structural regions.
#' @param n_yrs Number of model years.
#' @param n_seas Number of seasons per year.
#' @param n_fleets Total number of fishing fleets.
#' @param n_sexes Number of biological sexes.
#' @param n_model_bins Number of internal model bins.
#' @param n_obs_bins Number of observational bins.
#' @param age_or_len Flag indicating age-based or length-based composition.
#' @param AgeingErrorFn Function \code{(y, f)} returning the ageing error matrix
#'   for a given year and fleet, or the length bin map, which ignores both. Fleet
#'   specific because a fishery and a survey need not read ages the same way.
#' @param addtocomp Small constant added to proportions before normalization.
#' @param BinsArr Optional \code{[n_obs_bins x n_fleets]} 0/1 array naming the
#'   observed bins each fleet is fitted over, or \code{NULL} (default) for all
#'   bins. Must be the same array handed to \code{\link{pack_comp_osa}}, since
#'   the strides walked here are sized on it.
#' @param family Character string specifying the likelihood type, either "discrete" or "continuous".
#' @param zero_init Logical; whether to zero out the nLL array on entry.
#' @param pop Logical; if TRUE, evaluations account for the population structure layer.
#' @param n_pop Number of population structures or pools.
#'
#' @return Updated \code{nLL_arr} containing the evaluated negative log-likelihood values.
#' @keywords internal
eval_comp_osa = function(nLL_arr, tracked, ExpArrFn,
                         UseArr, TypeMat, LikeTypeVec,
                         ISSArr, lnThetaArr, lnThetaAggVec,
                         LNcorrArr, LNcorrAggVec,
                         n_regions, n_yrs, n_seas, n_fleets, n_sexes,
                         n_model_bins, n_obs_bins, age_or_len,
                         AgeingErrorFn, addtocomp, BinsArr = NULL,
                         family = "discrete", zero_init = TRUE,
                         pop = FALSE, n_pop = 1) {

  "c"   <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  if(is.null(tracked)) return(nLL_arr)

  fam_of = function(lt) if(lt %in% c(0,1)) "discrete" else "continuous"

  # Must match pack_comp_osa exactly, or the pointer k walks out of step
  bins_of = function(f) {
    if(is.null(BinsArr)) return(seq_len(n_obs_bins))
    if(nrow(BinsArr) != n_obs_bins) {
      stop("bin selection array has ", nrow(BinsArr), " rows but the observations carry ", n_obs_bins,
           " bins. The *_bins argument must be indexed on the observed bins of the stream it restricts.")
    }
    which(BinsArr[,f] == 1)
  }

  d = dim(nLL_arr)
  if(zero_init) nLL_arr = RTMB::AD(as.vector(nLL_arr) * 0)
  else          nLL_arr = RTMB::AD(as.vector(nLL_arr))
  dim(nLL_arr) = d

  k = 1
  for(p in 1:n_pop) {
    for(y in 1:n_yrs) {
      for(f in 1:n_fleets) {
        if(fam_of(LikeTypeVec[f]) != family) next
        for(seas in 1:n_seas) {

          use_vec = if(pop) UseArr[p, , y, seas, f] else UseArr[, y, seas, f]
          if(sum(use_vec) >= 1) {

            ct   = TypeMat[y, f]
            lt   = LikeTypeVec[f]
            n_ru = sum(use_vec == 1)

            fit_bins = bins_of(f)
            n_bins   = length(fit_bins)

            # slice length depends on family + comp type (LN drops ALR reference)
            if(lt %in% c(0,1)) {
              slice_length = if(ct == 0) n_bins else n_ru * n_bins * n_sexes
            } else { # LN
              if(ct == 0)      slice_length = (n_bins - 1)
              else if(ct == 1) slice_length = n_ru * (n_bins - 1) * n_sexes
              else             slice_length = n_ru * (n_bins * n_sexes - 1)
            }

            active_obs_slice = tracked[seq.int(from = k, length.out = slice_length)]

            AE = if(is.null(AgeingErrorFn)) NA else AgeingErrorFn(y, f)

            ISS_g          = if(pop) ISSArr[p, , y, seas, , f] else ISSArr[, y, seas, , f]
            ln_theta_g     = if(pop) lnThetaArr[p, , , f]      else lnThetaArr[, , f]
            ln_theta_agg_g = if(pop) lnThetaAggVec[p, f]       else lnThetaAggVec[f]
            LNcorr_g       = if(pop) LNcorrArr[p, , , f, ]     else LNcorrArr[, , f, ]
            LNcorr_agg_g   = if(pop) LNcorrAggVec[p, f]        else LNcorrAggVec[f]

            comp_out = Get_Comp_Likelihoods_OSA(
              Exp = ExpArrFn(p, y, seas, f), Obs = active_obs_slice,
              ISS = ISS_g, ln_theta = ln_theta_g, ln_theta_agg = ln_theta_agg_g,
              LN_corr_pars = LNcorr_g, LN_corr_pars_agg = LNcorr_agg_g,
              Comp_Type = ct, Likelihood_Type = lt,
              n_regions = n_regions, n_model_bins = n_model_bins, n_obs_bins = n_obs_bins,
              n_sexes = n_sexes, age_or_len = age_or_len, AgeingError = AE,
              use = use_vec, addtocomp = addtocomp,
              comp_bins = if(n_bins < n_obs_bins) fit_bins else NULL
            )

            if(pop) nLL_arr[p, , y, seas, , f] = comp_out
            else    nLL_arr[, y, seas, , f]    = comp_out

            k = k + slice_length
          }
        }
      }
    }
  }
  nLL_arr
}
