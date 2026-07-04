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
  dim(Obs) = c(n_regions, n_obs_bins, n_sexes)

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


  return(comp_nLL) # return negative log likelihood
} # end function

#' Composition Data Likelihood (OSA variant)
#'
#' Computes multinomial (0), Dirichlet–multinomial (1), and logistic–normal
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
#'   * Comp_Type 0: \code{n_obs_bins - 1}
#'   * Comp_Type 1: \code{n_obs_bins - 1} per region/sex
#'   * Comp_Type 2: \code{n_obs_bins * n_sexes - 1} per region (one joint reference)
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
                                    addtocomp) {

  "c"   <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  rho_trans = function(x) 2 / (1 + exp(-2 * x)) - 1
  alr_mu    = function(p) log(p[-length(p)]) - log(p[length(p)])

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
    if(age_or_len == 0) {
      tmp_Exp = as.vector(matrix(tmp_Exp, nrow = 1) %*% AgeingError)
      tmp_Exp = tmp_Exp / sum(tmp_Exp)
    }
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
        Sigma = diag(rep(exp(ln_theta_agg)^2, n_obs_bins))
      } else {
        LN_corr_b = rho_trans(LN_corr_pars_agg)
        Sigma = get_AR1_CorrMat(n_obs_bins, LN_corr_b) * (exp(ln_theta_agg)^2 / (1 - LN_corr_b^2))
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

        if(age_or_len == 0) {
          tmp_Exp = as.vector((Exp[r,,s] / sum(Exp[r,,s])) %*% AgeingError)
          tmp_Exp = tmp_Exp / sum(tmp_Exp)
        } else {
          tmp_Exp = Exp[r,,s] / sum(Exp[r,,s])
        }
        tmp_Exp = tmp_Exp + const
        tmp_Exp = tmp_Exp / sum(tmp_Exp)

        if(Likelihood_Type %in% c(0,1)) {
          idx = seq(from = r + (s - 1) * n_ru * n_obs_bins, by = n_ru, length.out = n_obs_bins)
        }
        if(Likelihood_Type %in% c(2,3)) {
          idx = seq(from = r + (s - 1) * n_ru * (n_obs_bins - 1), by = n_ru, length.out = n_obs_bins - 1)
        }
        if(Likelihood_Type == 0) { # Multinomial
          comp_nLL[r,s] = -dmultinom_osa(Obs[idx], tmp_Exp, log = TRUE)
        }
        if(Likelihood_Type == 1) { # Dirichlet-multinomial
          comp_nLL[r,s] = -ddirmult_osa(Obs[idx], tmp_Exp * exp(ln_theta[r,s]) * ISS[r,s], log = TRUE)
        }
        if(Likelihood_Type %in% c(2,3)) { # Logistic-normal (Obs already ALR)
          if(Likelihood_Type == 2) {
            Sigma = diag(rep(exp(ln_theta[r,s])^2, n_obs_bins))
          } else {
            LN_corr_b = rho_trans(LN_corr_pars[r,s,1])
            Sigma = get_AR1_CorrMat(n_obs_bins, LN_corr_b) * (exp(ln_theta[r,s])^2 / (1 - LN_corr_b^2))
          }
          Sigma = Sigma[-nrow(Sigma), -ncol(Sigma)]
          mu = alr_mu(tmp_Exp)
          comp_nLL[r,s] = -RTMB::dmvnorm(x = Obs[idx], mu = mu, Sigma = Sigma, log = TRUE)
        }

      } # r
    } # s
  } # end Comp_Type 1

  # Comp_Type 2 (joint sex, split region)
  if(Comp_Type == 2) {
    for(r in 1:n_ru) {

      if(age_or_len == 0) {
        tmp_Exp = as.vector(t(as.vector(Exp[r,,] / sum(Exp[r,,]))) %*% kronecker(diag(n_sexes), AgeingError))
        tmp_Exp = tmp_Exp / sum(tmp_Exp)
      } else {
        tmp_Exp = as.vector(Exp[r,,] / sum(Exp[r,,]))
      }
      tmp_Exp = tmp_Exp + const
      tmp_Exp = tmp_Exp / sum(tmp_Exp)

      if(Likelihood_Type %in% c(0,1)) {
        idx = seq(from = r, by = n_ru, length.out = n_obs_bins * n_sexes)
      }
      if(Likelihood_Type %in% c(2,3,4)) {
        # joint drops ONE reference for the whole [bin x sex] stack
        Lred = n_obs_bins * n_sexes - 1
        idx  = seq(from = r, by = n_ru, length.out = Lred)

        if(Likelihood_Type == 2) { # iid
          Sigma = diag(rep(exp(ln_theta[r,1])^2, n_obs_bins * n_sexes))
        }
        if(Likelihood_Type == 3) { # 1D AR1
          LN_corr_b = rho_trans(LN_corr_pars[r,1,1])
          Sigma = get_AR1_CorrMat(n_obs_bins * n_sexes, LN_corr_b)
          Sigma = Sigma * (exp(ln_theta[r,1])^2 / (1 - LN_corr_b^2))
        }
        if(Likelihood_Type == 4) { # 2D AR1 (age/len x sex)
          LN_corr_b = rho_trans(LN_corr_pars[r,1,1])
          LN_corr_s = rho_trans(LN_corr_pars[r,1,2])
          mat1 = get_Constant_CorrMat(n_sexes, LN_corr_s)
          mat2 = get_AR1_CorrMat(n_obs_bins, LN_corr_b)
          Sigma = Matrix::kronecker(mat1, mat2) * (exp(ln_theta[r,1])^2 / (1 - LN_corr_s^2) / (1 - LN_corr_b^2))
        }
        Sigma = Sigma[-nrow(Sigma), -ncol(Sigma)]  # drop last row/col -> Lred x Lred
        mu = alr_mu(tmp_Exp)                        # ALR of expectation -> length Lred
      }

      if(Likelihood_Type == 0) { # Multinomial
        comp_nLL[r,1] = -dmultinom_osa(Obs[idx], tmp_Exp, log = TRUE)
      }
      if(Likelihood_Type == 1) { # Dirichlet-multinomial
        comp_nLL[r,1] = -ddirmult_osa(Obs[idx], tmp_Exp * exp(ln_theta[r,1]) * ISS[r,1], log = TRUE)
      }
      if(Likelihood_Type %in% c(2,3,4)) { # Logistic-normal (Obs already ALR)
        comp_nLL[r,1] = -RTMB::dmvnorm(x = Obs[idx], mu = mu, Sigma = Sigma, log = TRUE)
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
#'
#' @return Flat OBS vector or NULL if no fleet of this family is present.
#' @keywords internal
pack_comp_osa = function(ObsArr, ISSArr, WtArr, UseArr, TypeMat, LikeTypeVec,
                         n_yrs, n_seas, n_fleets, n_sexes, addtocomp,
                         family = "discrete", pop = FALSE, n_pop = 1) {

  "c"   <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  n_obs_bins = if(pop) dim(ObsArr)[5] else dim(ObsArr)[4]

  fam_of = function(lt) if(lt %in% c(0,1)) "discrete" else "continuous"

  # additive log-ratio: drop last element as reference
  alr = function(p) log(p[-length(p)]) - log(p[length(p)])

  clean = list()
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

            if(lt %in% c(0,1)) {
              # Discrete: scale to counts, per region/sex
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
              g = if(ct == 0) as.vector(obs_slice[1, , 1]) else as.vector(obs_slice)

            } else {
              # Continuous (LN): ALR-transform per block, drop reference bin
              # Layout mirrors the discrete path so the likelihood's strided idx
              # works unchanged (just n_obs_bins-1 instead of n_obs_bins):
              #   ct 0: single vector, length n_obs_bins-1
              #   ct 1: [n_ru, n_obs_bins-1, n_sexes] as.vector (region-fastest)
              #   ct 2: per region, ALR of the full [bin x sex] stack
              #         (length n_obs_bins*n_sexes-1), region-fastest via seq stride
              if(ct == 0) {
                pr = (obs_slice[1, , 1] + addtocomp) / sum(obs_slice[1, , 1] + addtocomp)
                g  = alr(as.vector(pr))
              } else if(ct == 1) {
                # ALR each (region,sex) -> [n_ru, n_obs_bins-1, n_sexes], as.vector col-major
                arr = array(0, dim = c(n_ru, n_obs_bins - 1, n_sexes))
                arr = RTMB::AD(arr); dim(arr) = c(n_ru, n_obs_bins - 1, n_sexes)
                for(rr in 1:n_ru) {
                  for(s in 1:n_sexes) {
                    pr = (obs_slice[rr, , s] + addtocomp) / sum(obs_slice[rr, , s] + addtocomp)
                    arr[rr, , s] = alr(as.vector(pr))
                  }
                }
                g = as.vector(arr)          # region-fastest, matches strided idx
              } else {
                # joint: per region ALR of [bin x sex] stack. Store as
                # [n_ru, (n_obs_bins*n_sexes - 1)] and as.vector region-fastest,
                # matching idx = seq(from=r, by=n_ru, length.out=(n_obs_bins-1)*n_sexes)
                # NOTE: joint drops ONE ref for the whole stack -> length
                #       n_obs_bins*n_sexes - 1, NOT (n_obs_bins-1)*n_sexes.
                Lred = n_obs_bins * n_sexes - 1
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
          }
        }
      }
    }
  }
  if(length(clean) == 0) return(NULL)
  do.call(c, clean)
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
#'   \item Comp_Type 0: \code{n_obs_bins}
#'   \item Comp_Type 1/2: \code{n_ru x n_obs_bins x n_sexes}
#' }
#'
#' Logistic-normal (LikeType 2,3,4):
#' \itemize{
#'   \item Comp_Type 0: \code{n_obs_bins - 1}
#'   \item Comp_Type 1: \code{n_ru x (n_obs_bins - 1) x n_sexes}
#'   \item Comp_Type 2: \code{n_ru x (n_obs_bins x n_sexes - 1)}
#' }
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
#' @param AgeingErrorFn Function returning the ageing error matrix for a given year.
#' @param addtocomp Small constant added to proportions before normalization.
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
                         AgeingErrorFn, addtocomp,
                         family = "discrete", zero_init = TRUE,
                         pop = FALSE, n_pop = 1) {

  "c"   <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  if(is.null(tracked)) return(nLL_arr)

  fam_of = function(lt) if(lt %in% c(0,1)) "discrete" else "continuous"

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

            # slice length depends on family + comp type (LN drops ALR reference)
            if(lt %in% c(0,1)) {
              slice_length = if(ct == 0) n_obs_bins else n_ru * n_obs_bins * n_sexes
            } else { # LN
              if(ct == 0)      slice_length = (n_obs_bins - 1)
              else if(ct == 1) slice_length = n_ru * (n_obs_bins - 1) * n_sexes
              else             slice_length = n_ru * (n_obs_bins * n_sexes - 1)
            }

            active_obs_slice = tracked[k:(k + slice_length - 1)]

            AE = if(is.null(AgeingErrorFn)) NA else AgeingErrorFn(y)

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
              use = use_vec, addtocomp = addtocomp
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
