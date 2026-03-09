#' Composition Data Likelihood
#'
#' Computes the negative log-likelihood contribution for composition data
#' (age or length) for a single year and fleet. The function supports multiple
#' composition parameterizations and likelihood families commonly used in
#' stock assessment models, including multinomial, Dirichlet–multinomial,
#' and logistic–normal likelihoods.
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
#' @param ln_theta Log overdispersion parameters used for Dirichlet–
#'   multinomial or logistic–normal likelihoods, indexed by
#'   \eqn{[region \times sex]}.
#' @param LN_corr_pars Logistic–normal correlation parameters used for
#'   correlated logistic–normal likelihoods, dimensioned by
#'   \eqn{[region \times sex \times parameters]}.
#' @param LN_corr_pars_agg Logistic–normal correlation parameters used
#'   when compositions are aggregated.
#' @param Comp_Type Integer specifying the composition parameterization:
#'   \itemize{
#'     \item \code{0} – Aggregated compositions across sexes and regions.
#'     \item \code{1} – Compositions split by sex and region (no implicit
#'     sex or region ratio information).
#'     \item \code{2} – Joint compositions across sexes but split by region
#'     (implicit sex ratio information).
#'   }
#' @param Likelihood_Type Integer specifying the likelihood family:
#'   \itemize{
#'     \item \code{0} – Multinomial.
#'     \item \code{1} – Dirichlet–multinomial.
#'     \item \code{2} – Logistic–normal with independent bins.
#'     \item \code{3} – Logistic–normal with AR(1) correlation across bins.
#'     \item \code{4} – Logistic–normal with AR(1) correlation across bins
#'     and constant correlation across sexes.
#'   }
#' @param n_regions Number of regions modeled.
#' @param n_model_bins Number of composition bins used internally in the model.
#' @param n_obs_bins Number of observed composition bins.
#' @param n_sexes Number of sexes modeled.
#' @param age_or_len Indicator for composition type:
#'   \itemize{
#'     \item \code{0} – Age compositions.
#'     \item \code{1} – Length compositions.
#'   }
#' @param AgeingError Ageing error matrix used to map model age bins to
#'   observed age bins.
#' @param use Integer vector indicating which regions have observations
#'   (\code{1} = use data, \code{0} = ignore).
#' @param addtocomp Small constant added to compositions to avoid numerical
#'   issues when zeros are present.
#'
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
                                addtocomp
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
  Exp = array(Exp, dim = c(n_regions, n_model_bins, n_sexes)) # using n_obs_bins, because non-square ageing error matrix will collapse to n_obs_bins
  ISS = array(ISS, dim = c(n_regions, n_sexes))
  Wt_Mltnml = array(Wt_Mltnml, dim = c(n_regions, n_sexes))
  ln_theta = array(ln_theta, dim = c(n_regions, n_sexes))
  LN_corr_pars = array(LN_corr_pars, dim = c(n_regions, n_sexes, 3))

  # filter regions that have obs
  Obs = Obs[which(use == 1),,,drop = FALSE]
  Exp = Exp[which(use == 1),,,drop = FALSE]

  # Aggregated comps by sex and region
  if(Comp_Type == 0) {
    tmp_Exp = matrix(rowSums(matrix(Exp, nrow = n_model_bins)) / (n_sexes * n_regions), nrow = 1) # aggregate
    tmp_Exp = tmp_Exp / sum(tmp_Exp) # normalize

    # Expected age bins get collapsed to observed age bins if ageing error is non-square
    if(age_or_len == 0) {
      tmp_Exp = tmp_Exp %*% AgeingError # apply ageing error
      tmp_Exp = as.vector((tmp_Exp) / sum(tmp_Exp)) # renormalize
    }

    if(age_or_len == 1) tmp_Exp = as.vector((tmp_Exp) / sum(tmp_Exp)) # renormalize (lengths)

    # Multinomial likelihood
    if(Likelihood_Type == 0) { # Note that this indexes 1 because it's only a single sex and single region
      tmp_Obs = (Obs[1,,1]) / sum(Obs[1,,1]) # Normalize observed values
      ESS = ISS[1,1] * Wt_Mltnml[1,1] # Effective sample size
      comp_nLL[1,1] = -1 * ESS * sum(((tmp_Obs + const) * log(tmp_Exp + const))) # ADMB multinomial likelihood
      comp_nLL[1,1] = comp_nLL[1,1] - -1 * ESS * sum(((tmp_Obs + const) * log(tmp_Obs + const))) # Multinomial offset (subtract offset from actual likelihood)
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
      comp_nLL[1,1] = -1 * dlogistnormal(obs = tmp_Obs, pred = tmp_Exp, Sigma_or_Q = Sigma, type = 'dmvnorm', TRUE) # Logistic Normal likelihood (iid)
    } # end if logistic normal (iid)

    if(Likelihood_Type == 3) {
      # Dealing with zeros
      tmp_Obs = Obs[1,,1] / sum(Obs[1,,1])
      zeros = which(tmp_Obs == 0)

      # Construct Sigma
      LN_corr_b = rho_trans(LN_corr_pars_agg) # correlation by age / length
      Sigma =  get_AR1_CorrMat(n_obs_bins, LN_corr_b)
      Sigma = Sigma * (exp(ln_theta_agg)^2 / (1 - LN_corr_b^2))

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
      comp_nLL[1,1] = -1 * dlogistnormal(obs = tmp_Obs, pred = tmp_Exp, Sigma_or_Q = Sigma, type = 'dmvnorm', TRUE) # Logistic Normal likelihood (1dar1)
    } # end if logistic normal (1dar1)

  } # end if aggregated comps across sexes and regions

  # 'Split' comps by sex and region (no implicit sex ratio information)
  if(Comp_Type == 1) {
    for(s in 1:n_sexes) {
      for(r in 1:n_regions_obs_use) {
        # Expected Values
        if(age_or_len == 0) {
          tmp_Exp = ((Exp[r,,s]) / sum(Exp[r,,s])) %*% AgeingError # Normalize temporary variable (ages)
          tmp_Exp = tmp_Exp / sum(tmp_Exp) # renormalize
        }
        if(age_or_len == 1) tmp_Exp = (Exp[r,,s]) / sum(Exp[r,,s]) # Normalize lengths

        # Multinomial likelihood
        if(Likelihood_Type == 0) {
          tmp_Obs = (Obs[r,,s]) / sum(Obs[r,,s]) # Normalize observed temporary variable
          ESS = ISS[r,s] * Wt_Mltnml[r,s] # Effective sample size
          comp_nLL[r,s] = -1 * ESS * sum(((tmp_Obs + const) * log(tmp_Exp + const))) # ADMB multinomial likelihood
          comp_nLL[r,s] = comp_nLL[r,s] - -1 * ESS * sum(((tmp_Obs + const) * log(tmp_Obs + const))) # Multinomial offset (subtract offset from actual likelihood)
        } # end if multinomial likelihood

        if(Likelihood_Type == 1) {
          tmp_Obs = (Obs[r,,s] + const) / sum(Obs[r,,s] + const) # Normalize observed temporary variable
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
          comp_nLL[r,s] = -1 * dlogistnormal(obs = tmp_Obs, pred = tmp_Exp, Sigma_or_Q = Sigma, type = 'dmvnorm', TRUE) # Logistic Normal likelihood (iid)
        } # end if logistic normal

        if(Likelihood_Type == 3) {
          # Dealing with zeros
          tmp_Obs = Obs[r,,s] / sum(Obs[r,,s])
          zeros = which(tmp_Obs == 0)

          # Construct Sigma
          LN_corr_b = rho_trans(LN_corr_pars[r,s,1]) # correlation by age / length
          Sigma =  get_AR1_CorrMat(n_obs_bins, LN_corr_b)
          Sigma = Sigma * (exp(ln_theta[r,s])^2 / (1 - LN_corr_b^2))

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
          comp_nLL[r,s] = -1 * dlogistnormal(obs = tmp_Obs, pred = tmp_Exp, Sigma_or_Q = Sigma, type = 'dmvnorm', TRUE) # Logistic Normal likelihood (1dar1)
        } # end if logistic normal (1dar1)

      } # end r loop
    } # end s loop
  } # end if 'Split' comps by sex and region

  # Joint by sex, Split by region
  if(Comp_Type == 2) {
    for(r in 1:n_regions_obs_use) {

      # Expected values
      if(age_or_len == 0) { # if ages
        tmp_Exp = t(as.vector((Exp[r,,])/ sum(Exp[r,,]))) %*% kronecker(diag(n_sexes), AgeingError) # apply ageing error
        tmp_Exp = as.vector((tmp_Exp) / sum(tmp_Exp)) # renormalize to make sure sum to 1
      } # if ages
      if(age_or_len == 1) tmp_Exp = as.vector((Exp[r,,]) / sum((Exp[r,,]))) # Normalize temporary variable (lengths)

      # Multinomial likelihood
      if(Likelihood_Type == 0) { # Indexing by r for a given region since it's 'Split' by region and 1 for sex since it's 'Joint' for sex
        tmp_Obs = as.vector((Obs[r,,]) / sum(Obs[r,,])) # Normalize observed temporary variable
        ESS = ISS[r,1] * Wt_Mltnml[r,1] # Effective sample size
        comp_nLL[r,1] = -1 * ESS * sum(((tmp_Obs + const) * log(tmp_Exp + const))) # ADMB multinomial likelihood
        comp_nLL[r,1] = comp_nLL[r,1] - -1 * ESS * sum(((tmp_Obs + const) * log(tmp_Obs + const))) # Multinomial offset (subtract offset from actual likelihood)
      } # end if multinomial likelihood

      if(Likelihood_Type == 1) {
        tmp_Obs = as.vector((Obs[r,,] + const) / sum(Obs[r,,] + const)) # Normalize observed temporary variable
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

        comp_nLL[r,1] = -1 * dlogistnormal(obs = tmp_Obs, pred = tmp_Exp, Sigma_or_Q = Sigma, type = 'dmvnorm', TRUE) # Logistic Normal likelihood (iid)
      } # end if logistic normal (iid)

      if(Likelihood_Type == 3) {

        # Dealing with zeros
        tmp_Obs = Obs[r,,] / sum(Obs[r,,])
        zeros = which(tmp_Obs == 0)

        # Construct Sigma
        LN_corr_b = rho_trans(LN_corr_pars[r,1,1]) # correlation by age / length
        Sigma =  get_AR1_CorrMat(n_obs_bins * n_sexes, LN_corr_b)
        Sigma = Sigma * (exp(ln_theta[r,1])^2 / (1 - LN_corr_b^2))

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

        comp_nLL[r,1] = -1 * dlogistnormal(obs = tmp_Obs, pred = tmp_Exp, Sigma_or_Q = Sigma, type = 'dmvnorm', TRUE) # Logistic Normal likelihood (1dar1)
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
        comp_nLL[r,1] = -1 * dlogistnormal(obs = tmp_Obs, pred = tmp_Exp, Sigma_or_Q = Sigma, type = 'dmvnorm', TRUE) # Logistic Normal likelihood (1dar1 by age, constant corr by sex)
      }

    } # end r loop
  } # end if 'Joint' comps by sex, but 'Split' by region

  return(comp_nLL) # return negative log likelihood
} # end function
