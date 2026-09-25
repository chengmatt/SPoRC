# At-Age Observation Draws --------------------------------------------------

#' Draw one at-age data source for one region, year, season and fleet
#'
#' The operating model states an at-age observation the way the estimation model
#' reads it: summed over whichever of regions and sexes the fleet reports
#' together, read through the fleet's ageing error onto the observed ages, with
#' the fleet's own density and its own standard deviation. A data source summed
#' over regions is one number, so it is drawn once, when the region loop reaches
#' region one.
#'
#' @param numbers Array \code{[n_pop, n_regions, n_ages, n_sexes]} of the
#'   quantity at model age for this year, season and fleet.
#' @param weight Array shaped like \code{numbers}, read when \code{use_weight}.
#' @param use Integer array \code{[n_regions, n_obs_ages, n_sexes]} of use flags.
#' @param se Reported standard errors shaped like \code{use}.
#' @param ln_sigma Log-scale observation error, \code{[n_obs_ages, n_sexes]}.
#' @param type_code,like_code,form_code The fleet's aggregation, density and
#'   error-source codes.
#' @param use_weight Logical, \code{TRUE} for an observation in weight.
#' @param r Region the loop is on.
#' @param ageing_error Matrix \code{[n_ages, n_obs_ages]} reading model ages as
#'   observed ages for this year and fleet, or \code{NULL} for the identity.
#'
#' @return A list with \code{true} and \code{obs}, both \code{[n_obs_ages, n_sexes]}
#'   and \code{NA} wherever nothing was drawn.
#'
#' @keywords internal
sim_at_age_cell <- function(numbers, weight, use, se, ln_sigma, type_code, like_code,
                            form_code, use_weight, r, ageing_error = NULL) {

  n_regions <- dim(numbers)[2] # regions
  n_ages <- dim(numbers)[3] # model ages
  n_sexes <- dim(numbers)[4] # sexes
  n_obs_ages <- dim(use)[2] # ages the observations are recorded on

  if(is.null(ageing_error)) ageing_error <- base::diag(n_ages)
  if(!identical(as.integer(dim(ageing_error)), as.integer(c(n_ages, n_obs_ages)))) {
    stop("The ageing error for an at-age draw is ", paste(dim(ageing_error), collapse = " by "),
         " where ", n_ages, " model ages by ", n_obs_ages, " observed ages is expected. Set n_obs_ages ",
         "in Setup_Sim_Dim to the observed ages of AgeingError_input.")
  }

  split <- at_age_split(type_code)
  out_true <- array(NA_real_, dim = c(n_obs_ages, n_sexes))
  out_obs <- out_true

  # summed over regions the observation is one number, drawn once
  if(!split$region && r != 1) return(list(true = out_true, obs = out_obs))
  r_idx <- if(split$region) r else seq_len(n_regions)

  for(s in seq_len(n_sexes)) {

    if(!split$sex && s > 1) next
    s_idx <- if(split$sex) s else seq_len(n_sexes)

    # quantity at each model age, summed over the dims this fleet reports together
    at_model_age <- vapply(seq_len(n_ages), function(a) {
      if(use_weight) sum(numbers[,r_idx,a,s_idx] * weight[,r_idx,a,s_idx]) else sum(numbers[,r_idx,a,s_idx])
    }, numeric(1))

    for(a in seq_len(n_obs_ages)) {
      if(use[r,a,s] != 1) next
      read_as <- which(ageing_error[,a] != 0) # model ages read as this observed age
      val <- sum(at_model_age[read_as] * ageing_error[read_as,a])
      sd <- exp(ln_sigma[a,s])
      if(form_code != 0) sd <- at_age_obs_sd(se[r,a,s], sd, form_code)
      out_true[a,s] <- val
      out_obs[a,s] <- if(like_code == 0) val * exp(stats::rnorm(1, 0, sd))
                      else val + stats::rnorm(1, 0, sd)
    } # end a loop

  } # end s loop

  return(list(true = out_true, obs = out_obs))
}

#' Write a simulated at-age cell into its true and observed containers
#'
#' @param sim_env Environment holding the simulation containers.
#' @param data_source Data source tag, e.g. \code{"CatchAA"}.
#' @param drawn List returned by \code{\link{sim_at_age_cell}}.
#' @param r,y,seas,f,sim Region, year, season, fleet and replicate.
#'
#' @return \code{invisible(NULL)}, called for its side effect.
#'
#' @keywords internal
store_at_age_cell <- function(sim_env, data_source, drawn, r, y, seas, f, sim) {

  drawn_cells <- which(!is.na(drawn$true), arr.ind = TRUE)
  for(k in seq_len(nrow(drawn_cells))) {
    a <- drawn_cells[k,1]
    s <- drawn_cells[k,2]
    sim_env[[paste0("True", data_source)]][r,y,seas,a,s,f,sim] <- drawn$true[a,s]
    sim_env[[paste0("Obs", data_source)]][r,y,seas,a,s,f,sim] <- drawn$obs[a,s]
  } # end k loop

  return(invisible(NULL))
}

# Index Error Structures ----------------------------------------------------

#' Convert an index covariance matrix to common-factor parameters
#'
#' A multivariate normal index likelihood is supplied as a full covariance over
#' the observation vector, which cannot be drawn from one year at a time and does
#' not extend past the years it covers, so it is unusable for closed loop as
#' given. This decomposes it into a marginal scale and a single common factor,
#' \code{obs_t = pred_t + d_t (lambda_t u + sqrt(1 - lambda_t^2) e_t)}, with
#' \code{u} shared across years and \code{e_t} independent. Both are then drawn
#' per year, and a projection year past the end of the covariance simply reuses
#' the mean loading and scale with the same \code{u}.
#'
#'
#' @param S Covariance matrix over a fleet's observation vector.
#' @return List with \code{d} (marginal sd by observation) and \code{lambda}
#'   (factor loading by observation, in (-1, 1)).
#' @keywords internal
cov_to_factor <- function(S) {
  d <- sqrt(diag(S))
  R <- S / outer(d, d)
  e <- eigen(R, symmetric = TRUE)
  lam <- e$vectors[, 1] * sqrt(max(e$values[1], 0))
  if(mean(lam) < 0) lam <- -lam                 # sign of an eigenvector is arbitrary
  lam <- pmin(pmax(lam, -0.99), 0.99)           # keep 1 - lambda^2 strictly positive
  list(d = as.vector(d), lambda = as.vector(lam))
}

#' Precompute common-factor draw parameters for multivariate normal index fleets
#'
#' Validates each multivariate normal fleet's covariance against its use flags
#' with \code{\link{parse_idx_cov}}, factor-decomposes it with
#' \code{\link{cov_to_factor}}, and records where each simulated cell sits in
#' the covariance. Row \code{i} of the covariance is the \code{i}-th cell with
#' a use flag of 1 when scanning in array order (region varies fastest, then
#' year, then season), which is the same order the estimation model collects
#' the observation vector in, so the simulated and fitted series line up.
#'
#' @param cov_list List with one element per fleet, each a covariance matrix or
#'   \code{NULL}.
#' @param like_type_vals Integer vector of index likelihood codes; only fleets
#'   coded \code{2} are decomposed.
#' @param use_arr Array \code{[region, year, season, fleet]} of use flags. Its
#'   year dimension may be shorter than the simulation when projection years
#'   extend past the data.
#' @param n_fleets Integer. Number of fleets.
#' @param what Character. Name used in error messages.
#'
#' @return List with one element per fleet: \code{NULL} for non-mvn fleets, and
#'   otherwise \code{d} and \code{lambda} from \code{\link{cov_to_factor}}, a
#'   \code{row} lookup array \code{[region, year, season]} holding each used
#'   cell's covariance row (\code{NA} elsewhere), and the means \code{d_mean}
#'   and \code{lambda_mean} used for cells outside the covariance.
#' @keywords internal
build_idx_factor <- function(cov_list, like_type_vals, use_arr, n_fleets, what) {
  cov_parsed <- parse_idx_cov(cov_list, like_type_vals, use_arr, n_fleets, what)
  out <- vector("list", n_fleets)
  for(f in seq_len(n_fleets)) {
    if(like_type_vals[f] != 2) next
    fac <- cov_to_factor(cov_parsed[[f]])
    use_f <- array(use_arr[,,,f], dim = dim(use_arr)[1:3])
    row_arr <- array(NA_real_, dim = dim(use_f))
    row_arr[which(use_f == 1)] <- seq_len(sum(use_f == 1))
    out[[f]] <- list(
      d = fac$d,
      lambda = fac$lambda,
      row = row_arr,
      d_mean = mean(fac$d),
      lambda_mean = mean(fac$lambda)
    )
  } # end f loop
  return(out)
}

#' Resolve the factor scale and loading for one simulated index cell
#'
#' Looks a cell up in a fleet's factor decomposition from
#' \code{\link{build_idx_factor}}. A cell inside the covariance gets its own
#' marginal scale and loading; a cell outside it (a projection year, or a cell
#' the fleet never fit) gets the mean scale and loading, so closed loop can
#' keep drawing the series past the end of the data.
#'
#' @param mvn One fleet's element from \code{\link{build_idx_factor}}.
#' @param r,y,seas Integer indices of the simulated cell.
#' @return List with scalars \code{d} and \code{lambda}.
#' @keywords internal
resolve_idx_factor <- function(mvn, r, y, seas) {
  if(y <= dim(mvn$row)[2]) {
    i <- mvn$row[r, y, seas]
    if(!is.na(i)) return(list(d = mvn$d[i], lambda = mvn$lambda[i]))
  }
  list(d = mvn$d_mean, lambda = mvn$lambda_mean)
}

#' Draw an observed index given its likelihood family
#'
#' The estimation model supports lognormal, arithmetic-scale normal, and
#' multivariate normal index likelihoods. The simulator drew lognormal error
#' unconditionally, so a self-test of a normal or MVN index generated data under
#' the wrong error structure and then reported the mismatch as estimation bias.
#'
#' @param true True (error-free) index value(s).
#' @param se Observation standard deviation, on the log scale for lognormal and
#'   the arithmetic scale for normal. Unused for MVN, which takes its scale from
#'   the covariance instead; note the two are not interchangeable, the pollock
#'   trawl covariance has a diagonal about twice the reported SEs.
#' @param like_type 0 lognormal, 1 normal, 2 multivariate normal.
#' @param d,lambda Common-factor scale and loading for this observation, from
#'   \code{\link{cov_to_factor}}. MVN only.
#' @param u Shared factor draw for this fleet and replicate, kept constant across
#'   years. MVN only.
#' @keywords internal
draw_index_obs <- function(true, se, like_type = 0, d = NULL, lambda = NULL, u = NULL) {
  n <- length(true)
  if(like_type == 2) {
    if(is.null(d) || is.null(lambda) || is.null(u)) stop("A multivariate normal index draw needs d, lambda and u from cov_to_factor().")
    return(true + d * (lambda * u + sqrt(1 - lambda^2) * stats::rnorm(n, 0, 1)))
  }
  # guad against non-finite values
  se_n <- rep_len(se, n)
  ok <- is.finite(se_n) & se_n >= 0
  eps <- rep(NA, n)
  if(any(ok)) eps[ok] <- stats::rnorm(sum(ok), 0, se_n[ok])

  if(like_type == 1) return(true + eps)
  if(like_type == 0) return(true * exp(eps))
}

# Operating model
#
# Generates the data an assessment would see from the operating model's true state: catch and
# survey compositions and indices, and tag releases and recaptures, each with its own error.

# Composition and Age-at-Length Draws ---------------------------------------

#' Simulate age or length compositions
#'
#' Draws composition samples for one region, year, fleet, season and replicate
#' under the multinomial, Dirichlet-multinomial or logistic-normal likelihoods,
#' applying ageing error after the draw for age compositions. A
#' \code{comp_type} or \code{comp_like} of \code{999} returns \code{Obs}
#' unchanged.
#'
#' Joint compositions (\code{comp_type = 2}) apply ageing error across the
#' combined age by sex vector through the Kronecker product
#' \code{diag(n_sexes)} and \code{AgeingError}. Aggregated compositions
#' (\code{comp_type = 0}) are drawn only on the final region pass, from expected
#' proportions marginalized over regions and sexes. Under
#' \code{pop_specific = TRUE} each population is drawn separately from its own
#' sample sizes, dispersion and correlations, and that aggregation happens within
#' a population.
#'
#' @param r,y,f,seas,sim Region, year, fleet, season and replicate indices.
#' @param Exp Expected compositions \code{[n_pop × n_regions × n_yrs × n_seas ×
#'   n_cat × n_sexes × n_fleets × n_sims]}.
#' @param ISS Integer sample sizes \code{[n_regions × n_yrs × n_seas × n_sexes ×
#'   n_fleets × n_sims]}, read when \code{pop_specific = FALSE}.
#' @param AgeingError Ageing error matrices \code{[n_yrs × n_obs_ages × n_ages ×
#'   n_sims]}, ignored when \code{age_or_len = 1}.
#' @param comp_like Integer vector \code{[n_fleets]} of the likelihood per fleet:
#'   \code{0} multinomial, \code{1} Dirichlet-multinomial, \code{2}-\code{4} the
#'   logistic-normal forms.
#' @param ln_theta Log overdispersion \code{[n_regions × n_sexes × n_fleets]},
#'   read when \code{pop_specific = FALSE}.
#' @param corr_pars Logistic-normal correlation parameters \code{[n_regions ×
#'   n_sexes × n_fleets × n_corr_pars]}.
#' @param ln_theta_agg,corr_pars_agg Their counterparts for the aggregated
#'   compositions, each of length \code{n_fleets}.
#' @param comp_type Integer matrix \code{[n_yrs × n_fleets]}: \code{0} aggregated
#'   across regions, \code{1} split by sex, \code{2} joint across sexes,
#'   \code{999} no data.
#' @param n_sexes,n_pop,n_regions,n_cat Model dimensions, \code{n_cat} being the
#'   number of ages or lengths.
#' @param Obs Observed composition container, dimensioned like \code{Exp} and
#'   written in place.
#' @param pop_specific Logical. \code{TRUE} simulates each population separately
#'   from population-specific inputs.
#' @param age_or_len Integer. \code{0} for age compositions, which take ageing
#'   error, \code{1} for length compositions, which do not.
#' @param ISS_pop Population-specific sample sizes \code{[n_pop × n_regions ×
#'   n_yrs × n_seas × n_sexes × n_fleets × n_sims]}, read when
#'   \code{pop_specific = TRUE}.
#' @param pop_comp_like,pop_comp_type The likelihood and aggregation structure for
#'   the population-specific compositions, shaped as their aggregate
#'   counterparts.
#' @param ln_pop_theta Log overdispersion \code{[n_pop × n_regions × n_sexes ×
#'   n_fleets]}.
#' @param pop_corr_pars Logistic-normal correlation parameters \code{[n_pop ×
#'   n_regions × n_sexes × n_fleets × n_corr_pars]}.
#' @param ln_pop_theta_agg,pop_corr_pars_agg Their counterparts for the
#'   population-specific aggregated compositions, each \code{[n_pop ×
#'   n_fleets]}.
#'
#' @return \code{Obs} with the draws filled in at \code{[r, y, seas, , , f, sim]},
#'   or \code{[p, r, y, seas, , , f, sim]} under \code{pop_specific = TRUE}. Every
#'   other slice is unchanged.
#'
#' @keywords internal
simulate_comps <- function(r,
                           y,
                           f,
                           seas,
                           sim,
                           Exp,
                           ISS  = NULL,
                           AgeingError,
                           comp_like  = NULL,
                           ln_theta  = NULL,
                           corr_pars  = NULL,
                           ln_theta_agg  = NULL,
                           corr_pars_agg  = NULL,
                           comp_type = NULL,
                           n_sexes,
                           n_pop = NULL,
                           n_regions,
                           n_cat,
                           Obs,
                           pop_specific = FALSE,
                           ISS_pop = NULL,
                           pop_comp_like = NULL,
                           pop_comp_type = NULL,
                           ln_pop_theta = NULL,
                           pop_corr_pars = NULL,
                           ln_pop_theta_agg = NULL,
                           pop_corr_pars_agg = NULL,
                           age_or_len = 0) {

  if(!pop_specific && (comp_type[y,f] == 999 || comp_like[f] == 999)) return(Obs)
  if(pop_specific && (pop_comp_type[y,f] == 999 || pop_comp_like[f] == 999)) return(Obs)

  # helper functions
  get_expected <- function(prob_vec) prob_vec / sum(prob_vec)
  apply_error <- function(mat, age_or_len, AgeingError) {
    if(age_or_len == 0) return(mat %*% AgeingError)
    if(age_or_len == 1) return(mat)
  }

  if(!pop_specific) {
    if(age_or_len == 0) {
      if(comp_type[y,f] %in% c(0,1)) age_error_mat <- AgeingError[y,,,sim]
      if(comp_type[y,f] == 2) age_error_mat <- kronecker(diag(n_sexes), AgeingError[y,,,sim])
    }
  } else {
    if(age_or_len == 0) {
      if(pop_comp_type[y,f] %in% c(0,1)) age_error_mat <- AgeingError[y,,,sim]
      if(pop_comp_type[y,f] == 2) age_error_mat <- kronecker(diag(n_sexes), AgeingError[y,,,sim])
    }
  }

  if(pop_specific == FALSE) {
    # Split by sex
    if(comp_type[y,f] == 1) {
      for(s in 1:n_sexes) {

        tmp_prob <- apply(Exp[,r,y,seas,,s,f,sim, drop = FALSE], 5, sum) # extract compositions

        # multinomial
        if(comp_like[f] == 0) {
          Obs[r,y,seas,,s,f,sim] <- array(
            apply_error(as.vector(
              stats::rmultinom(n = 1, ISS[r,y,seas,s,f,sim], get_expected(tmp_prob))), age_or_len, age_error_mat),
            dim = dim(Obs[r,y,seas,,s,f,sim, drop = FALSE])
          )

          # dirichlet-multinomial
        } else if(comp_like[f] == 1) {
          Obs[r,y,seas,,s,f,sim] <- array(
            apply_error(as.vector(
              rdirM(
                n = 1,
                N = ISS[r,y,seas,s,f,sim],
                alpha = (exp(ln_theta[r,s,f]) * ISS[r,y,seas,s,f,sim]) * get_expected(tmp_prob)
              )
            ), age_or_len, age_error_mat),
            dim = dim(Obs[r,y,seas,,s,f,sim, drop = FALSE])
          )

          # logistic normal
        } else if(comp_like[f] %in% 2:7) {
          Obs[r,y,seas,,s,f,sim] <- array(
            apply_error(as.vector(
              rlogistnormal(
                exp = get_expected(tmp_prob),
                pars = c(exp(ln_theta[r,s,f]), comp_corr_natural(corr_pars[r,s,f,], comp_like[f])),
                comp_like = comp_like[f],
                n_sexes = n_sexes,
                ISS = ISS[r,y,seas,s,f,sim]
              )
            ), age_or_len, age_error_mat),
            dim = dim(Obs[r,y,seas,,s,f,sim, drop = FALSE])
          )
        }

      } # end s loop
    } # end split by sex

    # Joint compositions
    if(comp_type[y,f] == 2) {

      tmp_prob <- apply(Exp[,r,y,seas,,,f,sim, drop = FALSE], c(5,6), sum) # extract compositions

      # multinomial
      if(comp_like[f] == 0) {
        Obs[r,y,seas,,,f,sim] <- array(
          apply_error(as.vector(stats::rmultinom(1, ISS[r,y,seas,1,f,sim], get_expected(tmp_prob))),
                      age_or_len, age_error_mat),
          dim = dim(Obs[r,y,seas,,,f,sim, drop = FALSE])
        )

        # dirichlet-multinomial
      } else if(comp_like[f] == 1) {
        Obs[r,y,seas,,,f,sim] <- array(
          apply_error(as.vector(
            rdirM(
              n = 1,
              N = ISS[r,y,seas,1,f,sim],
              alpha = (exp(ln_theta[r,1,f]) * ISS[r,y,seas,1,f,sim]) * get_expected(tmp_prob)
            )
          ), age_or_len, age_error_mat),
          dim = dim(Obs[r,y,seas,,,f,sim, drop = FALSE])
        )

        # logistic normal
      } else if(comp_like[f] %in% 2:7) {
        Obs[r,y,seas,,,f,sim] <- array(
          apply_error(as.vector(
            rlogistnormal(
              exp = get_expected(tmp_prob),
              pars = c(exp(ln_theta[r,1,f]), comp_corr_natural(corr_pars[r,1,f,], comp_like[f])),
              comp_like = comp_like[f],
              n_sexes = n_sexes,
              ISS = ISS[r,y,seas,1,f,sim]
            )
          ), age_or_len, age_error_mat),
          dim = dim(Obs[r,y,seas,,,f,sim, drop = FALSE])
        )
      }

    } # end joint compositions

    # Aggregated comps across regions
    if(r == n_regions && comp_type[y,f] == 0) {

      # extract compositions
      tmp_prob <- apply(Exp[,,y,seas,,,f,sim, drop = FALSE], 5, sum)
      tmp_prob <- tmp_prob / sum(tmp_prob)

      # multinomial
      if(comp_like[f] == 0) {
        Obs[1,y,seas,,1,f,sim] <- array(
          apply_error(as.vector(stats::rmultinom(1, ISS[1,y,seas,1,f,sim], get_expected(tmp_prob))), age_or_len, age_error_mat),
          dim = dim(Obs[1,y,seas,,1,f,sim, drop = FALSE])
        )

        # dirichlet-multinomial
      } else if(comp_like[f] == 1) {
        Obs[1,y,seas,,1,f,sim] <- array(
          apply_error(as.vector(
            rdirM(
              n = 1,
              N = ISS[1,y,seas,1,f,sim],
              alpha = (exp(ln_theta_agg[f]) * ISS[1,y,seas,1,f,sim]) * get_expected(tmp_prob)
            )
          ), age_or_len, age_error_mat),
          dim = dim(Obs[1,y,seas,,1,f,sim, drop = FALSE])
        )

        # logistic normal
      } else if(comp_like[f] %in% 2:7) {
        Obs[1,y,seas,,1,f,sim] <- array(
          apply_error(as.vector(
            rlogistnormal(
              exp = get_expected(tmp_prob),
              pars = c(exp(ln_theta_agg[f]), comp_corr_natural(corr_pars_agg[f], comp_like[f])),
              comp_like = comp_like[f],
              n_sexes = n_sexes,
              ISS = ISS[1,y,seas,1,f,sim]
            )
          ), age_or_len, age_error_mat),
          dim = dim(Obs[1,y,seas,,1,f,sim, drop = FALSE])
        )
      }
    }
  } # end if not pop-specific

  if(pop_specific == TRUE) {
    for(p in 1:n_pop) {
      # Split by sex
      if(pop_comp_type[y,f] == 1) {
        for(s in 1:n_sexes) {

          tmp_prob <- Exp[p,r,y,seas,,s,f,sim, drop = FALSE]# extract compositions

          # multinomial
          if(pop_comp_like[f] == 0) {
            Obs[p,r,y,seas,,s,f,sim] <- array(
              apply_error(as.vector(
                stats::rmultinom(n = 1, ISS_pop[p,r,y,seas,s,f,sim], get_expected(tmp_prob))), age_or_len, age_error_mat),
              dim = dim(Obs[p,r,y,seas,,s,f,sim, drop = FALSE])
            )

            # dirichlet-multinomial
          } else if(pop_comp_like[f] == 1) {
            Obs[p,r,y,seas,,s,f,sim] <- array(
              apply_error(as.vector(
                rdirM(
                  n = 1,
                  N = ISS_pop[p,r,y,seas,s,f,sim],
                  alpha = (exp(ln_pop_theta[p,r,s,f]) * ISS_pop[p,r,y,seas,s,f,sim]) * get_expected(tmp_prob)
                )
              ), age_or_len, age_error_mat),
              dim = dim(Obs[p,r,y,seas,,s,f,sim, drop = FALSE])
            )

            # logistic normal
          } else if(pop_comp_like[f] %in% 2:7) {
            Obs[p,r,y,seas,,s,f,sim] <- array(
              apply_error(as.vector(
                rlogistnormal(
                  exp = get_expected(tmp_prob),
                  pars = c(exp(ln_pop_theta[p,r,s,f]), comp_corr_natural(pop_corr_pars[p,r,s,f,], pop_comp_like[f])),
                  comp_like = pop_comp_like[f],
                  n_sexes = n_sexes,
                  ISS = ISS_pop[p,r,y,seas,s,f,sim]
                )
              ), age_or_len, age_error_mat),
              dim = dim(Obs[p,r,y,seas,,s,f,sim, drop = FALSE])
            )
          }

        } # end s loop
      } # end split by sex

      # Joint compositions
      if(pop_comp_type[y,f] == 2) {

        tmp_prob <- Exp[p,r,y,seas,,,f,sim, drop = FALSE] # extract compositions

        # multinomial
        if(pop_comp_like[f] == 0) {
          Obs[p,r,y,seas,,,f,sim] <- array(
            apply_error(as.vector(stats::rmultinom(1, ISS_pop[p,r,y,seas,1,f,sim], get_expected(tmp_prob))),
                        age_or_len, age_error_mat),
            dim = dim(Obs[p,r,y,seas,,,f,sim, drop = FALSE])
          )

          # dirichlet-multinomial
        } else if(pop_comp_like[f] == 1) {
          Obs[p,r,y,seas,,,f,sim] <- array(
            apply_error(as.vector(
              rdirM(
                n = 1,
                N = ISS_pop[p,r,y,seas,1,f,sim],
                alpha = (exp(ln_pop_theta[p,r,1,f]) * ISS_pop[p,r,y,seas,1,f,sim]) * get_expected(tmp_prob)
              )
            ), age_or_len, age_error_mat),
            dim = dim(Obs[p,r,y,seas,,,f,sim, drop = FALSE])
          )

          # logistic normal
        } else if(pop_comp_like[f] %in% 2:7) {
          Obs[p,r,y,seas,,,f,sim] <- array(
            apply_error(as.vector(
              rlogistnormal(
                exp = get_expected(tmp_prob),
                pars = c(exp(ln_pop_theta[p,r,1,f]), comp_corr_natural(pop_corr_pars[p,r,1,f,], pop_comp_like[f])),
                comp_like = pop_comp_like[f],
                n_sexes = n_sexes,
                ISS = ISS_pop[p,r,y,seas,1,f,sim]
              )
            ), age_or_len, age_error_mat),
            dim = dim(Obs[p,r,y,seas,,,f,sim, drop = FALSE])
          )
        }

      } # end joint compositions

      # Aggregated comps across regions
      if(r == n_regions && pop_comp_type[y,f] == 0) {

        # extract compositions
        tmp_prob <- apply(Exp[p,,y,seas,,,f,sim, drop = FALSE], 5, sum)
        tmp_prob <- tmp_prob / sum(tmp_prob)

        # multinomial
        if(pop_comp_like[f] == 0) {
          Obs[p,1,y,seas,,1,f,sim] <- array(
            apply_error(as.vector(stats::rmultinom(1, ISS_pop[p,1,y,seas,1,f,sim], get_expected(tmp_prob))), age_or_len, age_error_mat),
            dim = dim(Obs[p,1,y,seas,,1,f,sim, drop = FALSE])
          )

          # dirichlet-multinomial
        } else if(pop_comp_like[f] == 1) {
          Obs[p,1,y,seas,,1,f,sim] <- array(
            apply_error(as.vector(
              rdirM(
                n = 1,
                N = ISS_pop[p,1,y,seas,1,f,sim],
                alpha = (exp(ln_pop_theta_agg[p,f]) * ISS_pop[p,1,y,seas,1,f,sim]) * get_expected(tmp_prob)
              )
            ), age_or_len, age_error_mat),
            dim = dim(Obs[p,1,y,seas,,1,f,sim, drop = FALSE])
          )

          # logistic normal
        } else if(pop_comp_like[f] %in% 2:7) {
          Obs[p,1,y,seas,,1,f,sim] <- array(
            apply_error(as.vector(
              rlogistnormal(
                exp = get_expected(tmp_prob),
                pars = c(exp(ln_pop_theta_agg[p,f]), comp_corr_natural(pop_corr_pars_agg[p,f], pop_comp_like[f])),
                comp_like = pop_comp_like[f],
                n_sexes = n_sexes,
                ISS = ISS_pop[p,1,y,seas,1,f,sim]
              )
            ), age_or_len, age_error_mat),
            dim = dim(Obs[p,1,y,seas,,1,f,sim, drop = FALSE])
          )
        }
      }
    } # end p loop
  } # end if pop-specific

  return(Obs)
}

#' Simulate conditional age-at-length observations
#'
#' Draws one age composition per length bin from the joint distribution of
#' length and age implied by the size-age transition matrix and the true numbers
#' at age (catch at age for the fishery, index at age for the survey). The joint
#' for a region, length bin and sex is \eqn{P(l \mid a) N_a} summed over
#' populations, and the draw for bin \eqn{l} is a multinomial (or
#' Dirichlet-multinomial) of \code{ISS[l]} fish across ages with that row as
#' the probability, which is the conditional \eqn{P(a \mid l)} by construction.
#' Ageing error is applied to the drawn counts the same way
#' \code{\link{simulate_comps}} applies it to marginal age compositions.
#'
#' Composition types follow \code{simulate_comps}: split by region and sex (1)
#' draws each sex separately, joint by sex (2) draws one sample across the age
#' by sex stack, and aggregated (0) pools regions and sexes and is drawn once
#' when the last region is reached. Only the multinomial and Dirichlet
#' multinomial families exist for CAAL.
#'
#' @param r,y,f,seas,sim Region, year, fleet, season and replicate indices.
#' @param SizeAgeTrans Array \code{[pop, region, year, season, len, age, sex,
#'   sim]} of \eqn{P(l \mid a)}.
#' @param AtAge Array \code{[pop, region, year, season, age, sex, fleet, sim]}
#'   of true numbers at age for this fleet type.
#' @param ISS Array \code{[region, year, season, len, sex, fleet, sim]} of fish
#'   aged per length bin. A zero skips the bin.
#' @param AgeingError Array \code{[year, model_age, obs_age, sim]}.
#' @param comp_like Likelihood code per fleet (0 multinomial, 1 DM, 999 none).
#' @param ln_theta Array \code{[region, sex, fleet]} of DM log overdispersion.
#' @param ln_theta_agg Vector of aggregated DM log overdispersion per fleet.
#' @param comp_type Matrix \code{[year, fleet]} of composition type codes.
#' @param n_sexes,n_regions,n_lens Dimension sizes.
#' @param Obs Array \code{[region, year, season, len, obs_age, sex, fleet, sim]}
#'   the draws are written into.
#'
#' @return The updated \code{Obs} array.
#'
#' @keywords internal
simulate_caal <- function(r, y, f, seas, sim, SizeAgeTrans, AtAge, ISS, AgeingError,
                          comp_like, ln_theta, ln_theta_agg, comp_type,
                          n_sexes, n_regions, n_lens, Obs) {

  if(comp_type[y,f] == 999 || comp_like[f] == 999) return(Obs)

  n_pop <- dim(AtAge)[1]
  ae <- AgeingError[y,,,sim] # model age by observed age
  if(is.null(dim(ae))) ae <- matrix(ae, nrow = 1) # a single model age arrives as a vector

  # P(l, a) for one region, length bin and sex, summed over populations
  joint_row <- function(rr, l, s) {
    out <- 0
    for(p in 1:n_pop) out <- out + SizeAgeTrans[p,rr,y,seas,l,,s,sim] * AtAge[p,rr,y,seas,,s,f,sim]
    return(out)
  }

  # one draw of N fish across the cells of prob, under the fleet's family
  draw <- function(N, prob, theta) {
    N <- round(N)
    if(N <= 0 || sum(prob) <= 0) return(rep(0, length(prob)))
    prob <- prob / sum(prob)
    if(comp_like[f] == 0) return(as.vector(stats::rmultinom(1, N, prob)))
    return(as.vector(rdirM(n = 1, N = N, alpha = (exp(theta) * N) * prob)))
  }

  for(l in 1:n_lens) {

    # Split by region and sex: each sex in this bin is its own sample
    if(comp_type[y,f] == 1) {
      for(s in 1:n_sexes) {
        counts <- draw(ISS[r,y,seas,l,s,f,sim], joint_row(r, l, s), ln_theta[r,s,f])
        Obs[r,y,seas,l,,s,f,sim] <- as.vector(counts %*% ae)
      } # end s loop
    }

    # Joint by sex: one sample across the age by sex stack, bin fastest then sex
    if(comp_type[y,f] == 2) {
      prob <- as.vector(sapply(1:n_sexes, function(s) joint_row(r, l, s)))
      counts <- draw(ISS[r,y,seas,l,1,f,sim], prob, ln_theta[r,1,f])
      counts <- as.vector(counts %*% kronecker(diag(n_sexes), ae))
      Obs[r,y,seas,l,,,f,sim] <- array(counts, dim = c(ncol(ae), n_sexes))
    }

    # Aggregated across regions and sexes, drawn once when the last region arrives
    if(r == n_regions && comp_type[y,f] == 0) {
      prob <- 0
      for(rr in 1:n_regions) for(s in 1:n_sexes) prob <- prob + joint_row(rr, l, s)
      counts <- draw(ISS[1,y,seas,l,1,f,sim], prob, ln_theta_agg[f])
      Obs[1,y,seas,l,,1,f,sim] <- as.vector(counts %*% ae)
    }

  } # end l loop

  return(Obs)
}

# Tag Recapture Draws -------------------------------------------------------

#' Simulate conventional tag recaptures for fishery fleets
#'
#' Draws observed recapture counts for one liberty, season and cohort cell from
#' the predicted recaptures, under the Poisson, negative binomial, or the
#' release- and recovery-conditioned multinomial and Dirichlet-multinomial. Dims
#' absent from \code{tag_recaptures_attr} are summed over and the recaptures are
#' written into index 1 of each.
#'
#' The release-conditioned forms express the predictions as proportions of the
#' tags released, appending a not-recaptured bin to complete the probability
#' vector before the draw and removing it afterwards. The recovery-conditioned
#' forms condition on the total predicted recaptures and need no such bin.
#'
#' @param conv_fish_tag_like Integer likelihood: \code{0} Poisson, \code{1}
#'   negative binomial, \code{2} and \code{3} the release- and
#'   recovery-conditioned multinomial, \code{4} and \code{5} the two
#'   Dirichlet-multinomials.
#' @param tag_recaptures_attr Which dims are attended in the recapture
#'   likelihood, built from \code{"p"}, \code{"a"} and \code{"s"} joined by
#'   underscores. Region and fleet are always kept, and the rest are summed over
#'   into index 1.
#' @param conv_tagged_fish Released tagged fish \code{[n_conv_tag_cohorts ×
#'   n_pop × n_ages × n_sexes × n_sims]}, the release sample size for the
#'   release-conditioned forms.
#' @param pred_conv_tag_fish_recap Predicted recaptures
#'   \code{[conv_tag_max_liberty × n_seas × n_conv_tag_cohorts × n_pop ×
#'   n_regions × n_ages × n_sexes × n_fish_fleets × n_sims]}.
#' @param obs_conv_tag_fish_recap Observed recaptures on the same dims, written
#'   in place at the \code{[ry, rseas, tc, ...]} slice.
#' @param ln_conv_fish_tag_theta Log overdispersion: the negative binomial size
#'   is \code{exp(ln_conv_fish_tag_theta)} and the Dirichlet-multinomial
#'   concentration \code{exp(ln_conv_fish_tag_theta) × N × p}. Ignored by the
#'   Poisson and the multinomial.
#' @param ry,rseas,tc,sim Years at liberty, recovery season, tag cohort and
#'   replicate indices.
#' @param n_pop,n_regions,n_ages,n_sexes,n_fish_fleets Dimension sizes.
#'
#' @return \code{obs_conv_tag_fish_recap} with the draws filled in at
#'   \code{[ry, rseas, tc, pop_idx, reg_idx, age_idx, sex_idx, flt_idx, sim]},
#'   the summed dims fixed at index 1.
#'
#' @keywords internal
simulate_conv_tag_fish_recaptures <- function(conv_fish_tag_like,
                                              tag_recaptures_attr,
                                              conv_tagged_fish,
                                              pred_conv_tag_fish_recap,
                                              obs_conv_tag_fish_recap,
                                              ln_conv_fish_tag_theta,
                                              ry,
                                              rseas,
                                              tc,
                                              sim,
                                              n_pop,
                                              n_regions,
                                              n_ages,
                                              n_sexes,
                                              n_fish_fleets
                                              ) {

  # get full dimensions of tag recaptures we simulate
  full_dims  <- c(n_pop, n_regions, n_ages, n_sexes, n_fish_fleets)
  attr_parts <- strsplit(tag_recaptures_attr, "_")[[1]]

  # Which of the 5 free dims (pop=1, region=2, age=3, sex=4, fleet=5) to retain
  # Region and fleet are always kept; pop/age/sex kept only if present in attr string
  keep_dims <- c(
    if("p" %in% attr_parts) 1,
    2,                            # region always kept
    if("a" %in% attr_parts) 3,
    if("s" %in% attr_parts) 4,
    5                             # fleet always kept
  )

  # get obs slice indices: full range if dim is kept, fixed at 1 if marginalized out
  pop_idx <- if("p" %in% attr_parts) seq_len(n_pop)   else 1
  age_idx <- if("a" %in% attr_parts) seq_len(n_ages)  else 1
  sex_idx <- if("s" %in% attr_parts) seq_len(n_sexes) else 1
  reg_idx <- seq_len(n_regions)   # always full
  flt_idx <- seq_len(n_fish_fleets) # always full

  # Function to marginalize tag recaptures. If dims == 5, then keep dimensions, otherwise,
  # marginalize and retain the kept dimensions
  marginalize <- function(vals) {
    tmp <- array(vals, dim = full_dims)
    if(length(keep_dims) < 5) apply(tmp, keep_dims, sum) else tmp
  }

  # Poisson or Neg Bin
  if(conv_fish_tag_like %in% c(0, 1)) {
    lambda <- marginalize(pred_conv_tag_fish_recap[ry, rseas, tc, , , , , , sim]) # get lambda / mu parameter
    # input and simulate
    obs_conv_tag_fish_recap[ry, rseas, tc, pop_idx, reg_idx, age_idx, sex_idx, flt_idx, sim] <-
      if(conv_fish_tag_like == 0) {
        stats::rpois(n = length(lambda), lambda = lambda)
      } else {
        stats::rnbinom(n = length(lambda), mu = lambda, size = exp(ln_conv_fish_tag_theta))
      }
  }

  # Multinomial or Dirichlet Multinomial (Release conditioned)
  if(conv_fish_tag_like %in% c(2, 4)) {
    tmp_n_tags_rel <- round(sum(conv_tagged_fish[tc, , , , sim])) # get sample size
    tmp_recap      <- marginalize(pred_conv_tag_fish_recap[ry, rseas, tc, , , , , , sim] / tmp_n_tags_rel) # marginalize
    tmp_probs      <- c(tmp_recap, 1 - sum(tmp_recap)) # get non-recaptured state

    # simualte
    tmp_sim_recap <-
      if(conv_fish_tag_like == 2) {
        stats::rmultinom(1, tmp_n_tags_rel, tmp_probs)
      } else {
        rdirM(n = 1, N = tmp_n_tags_rel, exp(ln_conv_fish_tag_theta) * tmp_n_tags_rel * tmp_probs)
      }

    # Drop the "not recaptured" bin and restore array shape
    tmp_sim_recap <- array(tmp_sim_recap[-length(tmp_sim_recap)], dim(tmp_recap))
    obs_conv_tag_fish_recap[ry, rseas, tc, pop_idx, reg_idx, age_idx, sex_idx, flt_idx, sim] <- tmp_sim_recap
  }

  # Multinomial or Dirichlet Multinomial (Recovery conditioned)
  if(conv_fish_tag_like %in% c(3, 5)) {
    tmp_n_tags_recap <- round(sum(pred_conv_tag_fish_recap[ry, rseas, tc, , , , , , sim])) # get sample size
    tmp_probs        <- marginalize(pred_conv_tag_fish_recap[ry, rseas, tc, , , , , , sim] / tmp_n_tags_recap)

    # simulate
    tmp_sim_recap <-
      if(conv_fish_tag_like == 3) {
        stats::rmultinom(1, tmp_n_tags_recap, c(tmp_probs))
      } else {
        rdirM(n = 1, N = tmp_n_tags_recap, exp(ln_conv_fish_tag_theta) * tmp_n_tags_recap * c(tmp_probs))
      }

    # input
    obs_conv_tag_fish_recap[ry, rseas, tc, pop_idx, reg_idx, age_idx, sex_idx, flt_idx, sim] <- tmp_sim_recap
  }

  return(obs_conv_tag_fish_recap)

}

#' Marginalize conventional fishery tag arrays across unattended dimensions
#'
#' Collapses population, age, and/or sex dimensions of a tag count array
#' \code{[n_pop × n_ages × n_sexes]} by summing over dimensions absent from
#' \code{tag_recaptures_attr}, placing the result into index 1 of the
#' corresponding dimension and zeroing all other indices. Region and fleet are
#' not handled here (they are managed at the calling level). The function is
#' used to align the release array \code{conv_tagged_fish} with the attended
#' resolution of the recapture likelihood.
#'
#' @param vals Numeric vector or array of tag counts, interpreted as a
#'   \code{[n_pop × n_ages × n_sexes]} array.
#' @param tag_recaptures_attr Character string specifying attended dimensions.
#'   Same format as \code{conv_fish_tag_attr} in
#'   \code{\link{Setup_Sim_Tagging}}: any combination of \code{"p"},
#'   \code{"a"}, \code{"s"} joined by underscores.
#' @param n_pop Integer. Number of populations.
#' @param n_ages Integer. Number of age classes.
#' @param n_sexes Integer. Number of sexes.
#'
#' @return Array \code{[n_pop × n_ages × n_sexes]} with unattended dimensions
#'   summed into index 1 and all other indices set to zero.
#'
#'
#' @keywords internal
marginalize_conv_fish_tags <- function(vals,
                                       tag_recaptures_attr,
                                       n_pop,
                                       n_ages,
                                       n_sexes) {

  full_dims  <- c(n_pop, n_ages, n_sexes)
  attr_parts <- strsplit(tag_recaptures_attr, "_")[[1]]
  tmp <- array(vals, dim = full_dims)  # temporary array

  if (!("p" %in% attr_parts) && n_pop > 1) {
    summed        <- apply(tmp, c(2, 3), sum) # collapse pop
    tmp[]         <- 0 # zero out stuff
    tmp[1, , ]  <- summed # input into 1st pop
  }
  if (!("a" %in% attr_parts) && n_ages > 1) {
    summed        <- apply(tmp, c(1, 3), sum) # collapse ages
    tmp[]         <- 0 # zero out
    tmp[, 1, ]  <- summed # input into 1st age
  }
  if (!("s" %in% attr_parts) && n_sexes > 1) {
    summed        <- apply(tmp, c(1, 2), sum)  # collapse sexes
    tmp[]         <- 0 # zero out
    tmp[, , 1]  <- summed # input into 1st sex
  }

  return(tmp)
}

# Annual Observation Generators ---------------------------------------------

#' Generate fishery catches, compositions, and indices in simulation
#'
#' Takes retained and dead discard catch at age from Baranov's equation for every
#' population, region, season and fleet, converts them to catch at length when a
#' size-age key is available, and draws the observed catch, discard and fishery
#' indices under lognormal error together with the age and length compositions of
#' both retained and discarded catch. The composition draws go through
#' \code{\link{simulate_comps}} and follow the likelihood and aggregation type set
#' in \code{sim_env}.
#'
#' Composition draws are skipped where \code{Fmort = 0}, and the discard ones also
#' where retention selectivity is fully 1, so nothing is discarded. Discard
#' indices come in four units: abundance, biomass, abundance fraction and biomass
#' fraction. Under \code{ISS_FishAgeComps_fill = "F_pattern"} with feedback
#' active, the sample sizes for the current and prior years are rescaled by
#' \code{\link{predict_sim_fish_iss_fmort}} before sampling.
#'
#' @param y Integer. Year index.
#' @param sim Integer. Simulation replicate index.
#' @param sim_env Simulation environment from \code{\link{Setup_sim_env}},
#'   modified in place. It gains the retained and dead discard catch at age
#'   \code{CAA} and \code{DAA}, their at-length counterparts \code{CAL} and
#'   \code{DAL} when a size-age key is present, the true and observed regional
#'   catch, discard and fishery indices with their population-specific
#'   counterparts, the observed retained and discard age and length compositions
#'   with theirs, and the input sample sizes of all eight composition data
#'   sources.
#'
#' @return \code{invisible(NULL)}; everything is modified by reference within
#'   \code{sim_env}.
#'
#' @keywords internal
generate_fishery_catch_comp_idx <- function(y, sim, sim_env) {

  sim_env$y   <- y
  sim_env$sim <- sim

  with(sim_env, {
    for(seas in 1:n_seas) {

      # Season-integrated abundance for the spatial Baranov under continuous movement.
      # Computed once per season across all regions, since the integral couples them.
      if(move_timing == 2) {
        NAA_int <- array(0, dim = c(n_pop, n_regions, n_ages, n_sexes))
        for(p in 1:n_pop) for(a in 1:n_ages) for(s in 1:n_sexes) {
          NAA_int[p,,a,s] <- integrate_seas_abundance(NAA[p,,y,seas,a,s,sim], ZAA[p,,y,seas,a,s,sim],
                                                      Mrate[p,,,y,seas,a,s,sim], seasdur[seas], expm_nsub = expm_nsub)
        }
      }

      for(r in 1:n_regions) {
        for(f in 1:n_fish_fleets) {

          for(p in 1:n_pop) {

            if(move_timing == 2) {
              # Spatial Baranov: under continuous movement fish redistribute among regions
              # while dying, so catch uses the season-integrated abundance
              sim_env$CAA[p,r,y,seas,,,f,sim] <- (Fmort[r,y,seas,f,sim] * fish_sel[p,r,y,seas,,,f,sim] * ret_sel[p,r,y,seas,,,f,sim]) *
                NAA_int[p,r,,]
              sim_env$DAA[p,r,y,seas,,,f,sim] <- (Fmort[r,y,seas,f,sim] * fish_sel[p,r,y,seas,,,f,sim] * (1 - ret_sel[p,r,y,seas,,,f,sim]) * dmr[r,y,seas,f,sim]) *
                NAA_int[p,r,,]
            } else {
              # Baranov's catch equation (retained catch-at-age)
              sim_env$CAA[p,r,y,seas,,,f,sim] <- (Fmort[r,y,seas,f,sim] * fish_sel[p,r,y,seas,,,f,sim] * ret_sel[p,r,y,seas,,,f,sim]) / ZAA[p,r,y,seas,,,sim] *
                NAA[p,r,y,seas,,,sim] * (1 - exp(-ZAA[p,r,y,seas,,,sim]))

              # Baranov's catch equation (dead discard catch-at-age)
              sim_env$DAA[p,r,y,seas,,,f,sim] <- (Fmort[r,y,seas,f,sim] * fish_sel[p,r,y,seas,,,f,sim] * (1 - ret_sel[p,r,y,seas,,,f,sim]) * dmr[r,y,seas,f,sim]) / ZAA[p,r,y,seas,,,sim] *
                NAA[p,r,y,seas,,,sim] * (1 - exp(-ZAA[p,r,y,seas,,,sim]))
            }

            # Catch-at-length
            if((exists("SizeAgeTrans") && !is.null(SizeAgeTrans)) || (exists("SizeAgeTrans_fish") && !is.null(SizeAgeTrans_fish))) for(s in 1:n_sexes) sim_env$CAL[p,r,y,seas,,s,f,sim] <- (if(exists("SizeAgeTrans_fish") && !is.null(SizeAgeTrans_fish)) SizeAgeTrans_fish[p,r,y,seas,,,s,f,sim] else SizeAgeTrans[p,r,y,seas,,,s,sim]) %*% CAA[p,r,y,seas,,s,f,sim] # Retained Catch at length
            if((exists("SizeAgeTrans") && !is.null(SizeAgeTrans)) || (exists("SizeAgeTrans_fish") && !is.null(SizeAgeTrans_fish))) for(s in 1:n_sexes) sim_env$DAL[p,r,y,seas,,s,f,sim] <- (if(exists("SizeAgeTrans_fish") && !is.null(SizeAgeTrans_fish)) SizeAgeTrans_fish[p,r,y,seas,,,s,f,sim] else SizeAgeTrans[p,r,y,seas,,,s,sim]) %*% DAA[p,r,y,seas,,s,f,sim] # Discarded Catch at length

          } # end p loop


          # Regional Retained Catch
          if(catch_units[f] == 0) sim_env$TrueCatch[r,y,seas,f,sim] <- sum(CAA[,r,y,seas,,,f,sim]) # abundance
          if(catch_units[f] == 1) sim_env$TrueCatch[r,y,seas,f,sim] <- sum(CAA[,r,y,seas,,,f,sim] * WAA_fish[,r,y,seas,,,f,sim]) # biomass
          sim_env$ObsCatch[r,y,seas,f,sim] <- TrueCatch[r,y,seas,f,sim] * exp(stats::rnorm(1, 0, exp(ln_sigmaC[r,y,seas,f]))) # Observed Catch w/ lognormal deviations
          aa_dim <- c(n_pop, n_regions, n_ages, n_sexes)
          aa_obs_dim <- c(n_regions, n_obs_ages, n_sexes) # use flags and errors sit on the observed ages
          aa_wt <- array(WAA_fish[,,y,seas,,,f,sim], dim = aa_dim)

          if(exists("use_catch_aa") && use_catch_aa[f] == 1) {
            store_at_age_cell(sim_env, "CatchAA",
                              sim_at_age_cell(array(CAA[,,y,seas,,,f,sim], dim = aa_dim), aa_wt,
                                              array(UseCatchAA[,y,seas,,,f], dim = aa_obs_dim),
                                              array(ObsCatchAA_SE[,y,seas,,,f], dim = aa_obs_dim),
                                              array(ln_sigmaCAA[,,f], dim = c(n_obs_ages, n_sexes)),
                                              CatchAA_Type[y,f], CatchAA_LikeType[f], CatchAA_sigma_form[f],
                                              catch_units[f] == 1, r,
                                              ageing_error = array(AgeingError_fish[y,,,f,sim], dim = dim(AgeingError_fish)[2:3])),
                              r, y, seas, f, sim)
          }

          if(exists("use_discard_aa") && use_discard_aa[f] == 1) {
            daa <- array(DAA[,,y,seas,,,f,sim], dim = aa_dim)
            for(rr in 1:n_regions) daa[,rr,,] <- daa[,rr,,] / dmr[rr,y,seas,f,sim] # dead discards raised to the total
            store_at_age_cell(sim_env, "DiscardAA",
                              sim_at_age_cell(daa, aa_wt,
                                              array(UseDiscardAA[,y,seas,,,f], dim = aa_obs_dim),
                                              array(ObsDiscardAA_SE[,y,seas,,,f], dim = aa_obs_dim),
                                              array(ln_sigmaDAA[,,f], dim = c(n_obs_ages, n_sexes)),
                                              DiscardAA_Type[y,f], DiscardAA_LikeType[f], DiscardAA_sigma_form[f],
                                              discard_units[f] == 1, r,
                                              ageing_error = array(AgeingError_fish[y,,,f,sim], dim = dim(AgeingError_fish)[2:3])),
                              r, y, seas, f, sim)
          }

          # Population Specific Catch
          if(catch_units[f] == 0) sim_env$TrueCatch_pop[,r,y,seas,f,sim] <- apply(CAA[,r,y,seas,,,f,sim, drop = FALSE], 1, sum)  # abundance
          if(catch_units[f] == 1) sim_env$TrueCatch_pop[,r,y,seas,f,sim] <- apply(CAA[,r,y,seas,,,f,sim, drop = FALSE] * WAA_fish[,r,y,seas,,,f,sim, drop = FALSE], 1, sum)  # biomass
          sim_env$ObsCatch_pop[,r,y,seas,f,sim] <- sim_env$TrueCatch_pop[,r,y,seas,f,sim] * exp(stats::rnorm(n_pop, 0, exp(ln_sigmaC_pop[,r,y,seas,f])))

          # Regional Discards
          if(discard_units[f] == 0) sim_env$TrueDiscard[r,y,seas,f,sim] <- sum(DAA[,r,y,seas,,,f,sim]  / dmr[r,y,seas,f,sim]) # abd
          if(discard_units[f] == 1) sim_env$TrueDiscard[r,y,seas,f,sim] <- sum((DAA[,r,y,seas,,,f,sim]  / dmr[r,y,seas,f,sim]) * WAA_fish[,r,y,seas,,,f,sim]) # biom
          if(discard_units[f] == 2) { # abd frac
            total_catch <- CAA[,r,y,seas,,,f,sim, drop = FALSE] + DAA[,r,y,seas,,,f,sim, drop = FALSE] / dmr[r,y,seas,f,sim]
            sim_env$TrueDiscard[r,y,seas,f,sim] <- 1 - sum(CAA[,r,y,seas,,,f,sim]) / sum(total_catch)
          }
          if(discard_units[f] == 3) { # biom frac
            total_catch <- CAA[,r,y,seas,,,f,sim, drop = FALSE] + DAA[,r,y,seas,,,f,sim, drop = FALSE] / dmr[r,y,seas,f,sim]
            sim_env$TrueDiscard[r,y,seas,f,sim] <- 1 - sum(CAA[,r,y,seas,,,f,sim] * WAA_fish[,r,y,seas,,,f,sim]) / sum(total_catch * WAA_fish[,r,y,seas,,,f,sim, drop = FALSE])
          }

          # lognormal
          sim_env$ObsDiscard[r,y,seas,f,sim] <- TrueDiscard[r,y,seas,f,sim] * exp(stats::rnorm(1, 0, exp(ln_sigmaD[r,y,seas,f])))

          # Population Specific Discards
          if(discard_units[f] == 0) sim_env$TrueDiscard_pop[,r,y,seas,f,sim] <- apply(DAA[,r,y,seas,,,f,sim, drop = FALSE]  / dmr[r,y,seas,f,sim], 1, sum) #abd
          if(discard_units[f] == 1) sim_env$TrueDiscard_pop[,r,y,seas,f,sim] <- apply((DAA[,r,y,seas,,,f,sim, drop = FALSE] / dmr[r,y,seas,f,sim]) * WAA_fish[,r,y,seas,,,f,sim, drop = FALSE], 1, sum) # biom
          if(discard_units[f] == 2) { # abd frac
            total_catch_pop <- CAA[,r,y,seas,,,f,sim, drop = FALSE] + DAA[,r,y,seas,,,f,sim, drop = FALSE] / dmr[r,y,seas,f,sim]
            tmp_c <- apply(CAA[,r,y,seas,,,f,sim, drop = FALSE], 1, sum)
            tmp_total <- apply(total_catch_pop, 1, sum)
            sim_env$TrueDiscard_pop[,r,y,seas,f,sim] <- 1 - tmp_c / tmp_total
          }
          if(discard_units[f] == 3) { # biom frac
            total_catch_pop <- CAA[,r,y,seas,,,f,sim, drop = FALSE] + DAA[,r,y,seas,,,f,sim, drop = FALSE] / dmr[r,y,seas,f,sim]
            tmp_c <- apply(CAA[,r,y,seas,,,f,sim, drop = FALSE] * WAA_fish[,r,y,seas,,,f,sim, drop = FALSE], 1, sum)
            tmp_total <- apply(total_catch_pop * WAA_fish[,r,y,seas,,,f,sim, drop = FALSE], 1, sum)
            sim_env$TrueDiscard_pop[,r,y,seas,f,sim] <- 1 - tmp_c / tmp_total
          }

          # Lognormal
          sim_env$ObsDiscard_pop[,r,y,seas,f,sim] <- sim_env$TrueDiscard_pop[,r,y,seas,f,sim] * exp(stats::rnorm(n_pop, 0, exp(ln_sigmaD_pop[,r,y,seas,f])))

          # Fishery Index
          tmp_NAA <- NAA[,r,y,seas,,,sim, drop = FALSE]
          if(any(t_fish[r,seas,f] != 0))
            tmp_NAA <- tmp_NAA * exp(-t_fish[r,seas,f] * ZAA[,r,y,seas,,,sim, drop = FALSE])
          tmp_expl_abd <- sweep(tmp_NAA, c(1,5,6), fish_sel[,r,y,seas,,,f,sim, drop = FALSE] * ret_sel[,r,y,seas,,,f,sim, drop = FALSE], "*")
          tmp_expl_biom <- sweep(tmp_expl_abd, c(1,5,6), WAA_fish[,r,y,seas,,,f,sim, drop = FALSE], "*") # get exploitable abundance
          if(fish_idx_type[f] == 0) sim_env$TrueFishIdx[r,y,seas,f,sim] <- fish_q[r,y,f,sim] * sum(tmp_expl_abd) # True Fishery Index (abundance)
          if(fish_idx_type[f] == 1) sim_env$TrueFishIdx[r,y,seas,f,sim] <- fish_q[r,y,f,sim] * sum(tmp_expl_biom) # True Fishery Index (biomass)

          # observed index. an mvn fleet takes its scale from the covariance's factor decomposition
          # rather than the SE array, with one factor draw shared across its series per replicate
          fidx_like <- if(exists("FishIdx_LikeType")) FishIdx_LikeType[f] else 0
          if(fidx_like == 2) {
            if(is.na(fish_idx_u[f,sim])) sim_env$fish_idx_u[f,sim] <- stats::rnorm(1)
            fidx_fac <- resolve_idx_factor(fish_idx_mvn[[f]], r, y, seas)
            sim_env$ObsFishIdx[r,y,seas,f,sim] <- draw_index_obs(
              TrueFishIdx[r,y,seas,f,sim],
              NA,
              2,
              d = fidx_fac$d,
              lambda = fidx_fac$lambda,
              u = fish_idx_u[f,sim]
            )
          } else {
            sim_env$ObsFishIdx[r,y,seas,f,sim] <- draw_index_obs(TrueFishIdx[r,y,seas,f,sim], ObsFishIdx_SE[r,y,seas,f], fidx_like)
          }

          # population-specific index. the covariance describes the regional series only, so an mvn
          # fleet's population data source keeps lognormal error, mirroring the estimation model
          if(fish_idx_type[f] == 0) sim_env$TrueFishIdx_pop[,r,y,seas,f,sim] <- fish_q[r,y,f,sim] * apply(tmp_expl_abd[,1,1,1,,,1, drop = FALSE], 1, sum)  # abundance
          if(fish_idx_type[f] == 1) sim_env$TrueFishIdx_pop[,r,y,seas,f,sim] <- fish_q[r,y,f,sim] * apply(tmp_expl_biom[,1,1,1,,,1, drop = FALSE], 1, sum)  # biomass
          sim_env$ObsFishIdx_pop[,r,y,seas,f,sim] <- draw_index_obs(sim_env$TrueFishIdx_pop[,r,y,seas,f,sim], ObsFishIdx_pop_SE[,r,y,seas,f], if(fidx_like == 1) 1 else 0)

          # Fishery Compositions
          if(Fmort[r,y,seas,f,sim] > 0) { # only simulate if Fishing Mortality > 0

            # Retained Compositions
            # Age Compositions (Dynamic ISS based on feedback fishing mortality)
            if(exists("ISS_FishAgeComps_fill") && isTRUE(ISS_FishAgeComps_fill == "F_pattern") && isTRUE(run_feedback) && y >= feedback_start_yr + 1 && r == 1 && f == 1) {
              sim_env$ISS_FishAgeComps[,1:y,seas,,,sim] <- predict_sim_fish_iss_fmort(ISS_FishComps = ISS_FishAgeComps, Fmort = Fmort, y = y, sim = sim, seas = seas)
            }
            if(exists("ISS_FishAgeComps_pop_fill") && isTRUE(ISS_FishAgeComps_fill == "F_pattern") && isTRUE(run_feedback) && y >= feedback_start_yr + 1 && r == 1 && f == 1) {
              for(p in 1:n_pop) sim_env$ISS_FishAgeComps_pop[p,,1:y,seas,,,sim] <- predict_sim_fish_iss_fmort(
                ISS_FishComps = array(ISS_FishAgeComps_pop[p,,1:y,,,,sim],
                                                                                                                                    dim = c(n_regions, length(1:y), n_seas, n_sexes, n_fish_fleets, n_sims)),
                Fmort = Fmort,
                y = y,
                sim = sim,
                seas = seas
              )
            }

            # Length Compositions (Dynamic ISS based on feedback fishing mortality)
            if(exists("ISS_FishLenComps_fill") && isTRUE(ISS_FishLenComps_fill == "F_pattern") && isTRUE(run_feedback) && y >= feedback_start_yr + 1 && r == 1 && f == 1) {
              sim_env$ISS_FishLenComps[,1:y,seas,,,sim] <- predict_sim_fish_iss_fmort(ISS_FishComps = ISS_FishLenComps, Fmort = Fmort, y = y, sim = sim, seas = seas)
            }
            if(exists("ISS_FishLenComps_pop_fill") && isTRUE(ISS_FishLenComps_fill == "F_pattern") && isTRUE(run_feedback) && y >= feedback_start_yr + 1 && r == 1 && f == 1) {
              for(p in 1:n_pop) sim_env$ISS_FishLenComps_pop[p,,1:y,seas,,,sim] <- predict_sim_fish_iss_fmort(
                ISS_FishComps = array(ISS_FishLenComps_pop[p,,1:y,,,,sim],
                                                                                                                                    dim = c(n_regions, length(1:y), n_seas, n_sexes, n_fish_fleets, n_sims)),
                Fmort = Fmort,
                y = y,
                sim = sim,
                seas = seas
              )
            }

            # Sample fishery ages (non-population specific, retained compositions)
            sim_env$ObsFishAgeComps <- simulate_comps(r = r,
                                                      y = y,
                                                      seas = seas,
                                                      f = f,
                                                      sim = sim,
                                                      Exp = CAA,
                                                      ISS = ISS_FishAgeComps,
                                                      AgeingError = array(AgeingError_fish[,,,f,], dim = dim(AgeingError)),
                                                      comp_like = comp_fishage_like,
                                                      ln_theta = ln_FishAge_theta,
                                                      ln_theta_agg = ln_FishAge_theta_agg,
                                                      corr_pars = FishAge_corr_pars,
                                                      corr_pars_agg = FishAge_corr_pars_agg,
                                                      comp_type = FishAgeComps_Type,
                                                      n_sexes = n_sexes,
                                                      n_regions = n_regions,
                                                      n_cat = n_ages,
                                                      Obs = ObsFishAgeComps,
                                                      pop_specific = FALSE,
                                                      age_or_len = 0)

            # Sample fishery ages (population specific, retained compositions)
            sim_env$ObsFishAgeComps_pop <- simulate_comps(r = r,
                                                          y = y,
                                                          seas = seas,
                                                          f = f,
                                                          sim = sim,
                                                          Exp = CAA,
                                                          ISS_pop = ISS_FishAgeComps_pop,
                                                          AgeingError = array(AgeingError_fish[,,,f,], dim = dim(AgeingError)),
                                                          pop_comp_like = comp_fishage_pop_like,
                                                          ln_pop_theta = ln_FishAge_pop_theta,
                                                          ln_pop_theta_agg = ln_FishAge_pop_theta_agg,
                                                          pop_corr_pars = FishAge_pop_corr_pars,
                                                          pop_corr_pars_agg = FishAge_pop_corr_pars_agg,
                                                          pop_comp_type = FishAgeComps_pop_Type,
                                                          n_sexes = n_sexes,
                                                          n_regions = n_regions,
                                                          n_pop = n_pop,
                                                          n_cat = n_ages,
                                                          Obs = ObsFishAgeComps_pop,
                                                          pop_specific = TRUE,
                                                          age_or_len = 0)

            # Sample fishery lengths (retained compositions)
            if((exists("SizeAgeTrans") && !is.null(SizeAgeTrans)) || (exists("SizeAgeTrans_fish") && !is.null(SizeAgeTrans_fish))) {
              sim_env$ObsFishLenComps <- simulate_comps(r = r,
                                                        y = y,
                                                        seas = seas,
                                                        f = f,
                                                        sim = sim,
                                                        Exp = CAL,
                                                        ISS = ISS_FishLenComps,
                                                        AgeingError = NULL,
                                                        comp_like = comp_fishlen_like,
                                                        ln_theta = ln_FishLen_theta,
                                                        ln_theta_agg = ln_FishLen_theta_agg,
                                                        corr_pars = FishLen_corr_pars,
                                                        corr_pars_agg = FishLen_corr_pars_agg,
                                                        comp_type = FishLenComps_Type,
                                                        n_sexes = n_sexes,
                                                        n_regions = n_regions,
                                                        n_cat = n_lens,
                                                        Obs = ObsFishLenComps,
                                                        pop_specific = FALSE,
                                                        age_or_len = 1)

              # Sample fishery lengths (population specific)
              sim_env$ObsFishLenComps_pop <- simulate_comps(r = r,
                                                            y = y,
                                                            seas = seas,
                                                            f = f,
                                                            sim = sim,
                                                            Exp = CAL,
                                                            ISS_pop = ISS_FishLenComps_pop,
                                                            AgeingError = NULL,
                                                            pop_comp_like = comp_fishlen_pop_like,
                                                            ln_pop_theta = ln_FishLen_pop_theta,
                                                            ln_pop_theta_agg = ln_FishLen_pop_theta_agg,
                                                            pop_corr_pars = FishLen_pop_corr_pars,
                                                            pop_corr_pars_agg = FishLen_pop_corr_pars_agg,
                                                            pop_comp_type = FishLenComps_pop_Type,
                                                            n_sexes = n_sexes,
                                                            n_regions = n_regions,
                                                            n_pop = n_pop,
                                                            n_cat = n_lens,
                                                            Obs = ObsFishLenComps_pop,
                                                            pop_specific = TRUE,
                                                            age_or_len = 1)

              # sample fishery conditional age-at-length. the joint of length and age is formed inside
              # the sampler from SizeAgeTrans and CAA, so no joint array is kept per replicate
              if(exists("do_fish_caal") && isTRUE(do_fish_caal)) {
                sim_env$ObsFish_caal <- simulate_caal(
                  r = r,
                  y = y,
                  f = f,
                  seas = seas,
                  sim = sim,
                  SizeAgeTrans = if(exists("SizeAgeTrans_fish") && !is.null(SizeAgeTrans_fish)) array(SizeAgeTrans_fish[,,,,,,,f,], dim = dim(SizeAgeTrans_fish)[-8]) else SizeAgeTrans,
                  AtAge = CAA,
                  ISS = ISS_Fish_caal,
                  AgeingError = array(AgeingError_fish[,,,f,], dim = dim(AgeingError)),
                  comp_like = comp_fish_caal_like,
                  ln_theta = ln_Fish_caal_theta,
                  ln_theta_agg = ln_Fish_caal_theta_agg,
                  comp_type = Fish_caal_Type,
                  n_sexes = n_sexes,
                  n_regions = n_regions,
                  n_lens = n_lens,
                  Obs = ObsFish_caal
                )
              } # end fishery caal

            } # end if size age transition if availiable

            # if there is discarding occuring
            if(!all(DAA[,r,y,seas,,,f,sim] == 0)) {

              # Discarded Compositions (Dynamic ISS based on feedback fishing mortality)
              if(exists("ISS_FishAgeComps_discard_fill") && isTRUE(ISS_FishAgeComps_discard_fill == "F_pattern") && isTRUE(run_feedback) && y >= feedback_start_yr + 1 && r == 1 && f == 1) {
                sim_env$ISS_FishAgeComps_discard[,1:y,seas,,,sim] <- predict_sim_fish_iss_fmort(ISS_FishComps = ISS_FishAgeComps_discard, Fmort = Fmort, y = y, sim = sim, seas = seas)
              }
              if(exists("ISS_FishAgeComps_pop_discard_fill") && isTRUE(ISS_FishAgeComps_pop_discard_fill == "F_pattern") && isTRUE(run_feedback) && y >= feedback_start_yr + 1 && r == 1 && f == 1) {
                for(p in 1:n_pop) sim_env$ISS_FishAgeComps_discard_pop[p,,1:y,seas,,,sim] <- predict_sim_fish_iss_fmort(
                  ISS_FishComps = array(ISS_FishAgeComps_discard_pop[p,,1:y,,,,sim],
                                                                                                                                              dim = c(n_regions, length(1:y), n_seas, n_sexes, n_fish_fleets, n_sims)),
                  Fmort = Fmort,
                  y = y,
                  sim = sim,
                  seas = seas
                )
              }

              # Length Compositions (Dynamic ISS based on feedback fishing mortality)
              if(exists("ISS_FishLenComps_discard_fill") && isTRUE(ISS_FishLenComps_discard_fill == "F_pattern") && isTRUE(run_feedback) && y >= feedback_start_yr + 1 && r == 1 && f == 1) {
                sim_env$ISS_FishLenComps_discard[,1:y,seas,,,sim] <- predict_sim_fish_iss_fmort(ISS_FishComps = ISS_FishLenComps_discard, Fmort = Fmort, y = y, sim = sim, seas = seas)
              }
              if(exists("ISS_FishLenComps_pop_discard_fill") && isTRUE(ISS_FishLenComps_pop_discard_fill == "F_pattern") && isTRUE(run_feedback) && y >= feedback_start_yr + 1 && r == 1 && f == 1) {
                for(p in 1:n_pop) sim_env$ISS_FishLenComps_discard_pop[p,,1:y,seas,,,sim] <- predict_sim_fish_iss_fmort(
                  ISS_FishComps = array(ISS_FishLenComps_discard_pop[p,,1:y,,,,sim],
                                                                                                                                              dim = c(n_regions, length(1:y), n_seas, n_sexes, n_fish_fleets, n_sims)),
                  Fmort = Fmort,
                  y = y,
                  sim = sim,
                  seas = seas
                )
              }

              # Sample fishery lengths (non-population specific, discard compositions)
              sim_env$ObsFishAgeComps_discard <- simulate_comps(r = r,
                                                                y = y,
                                                                seas = seas,
                                                                f = f,
                                                                sim = sim,
                                                                Exp = DAA,
                                                                ISS = ISS_FishAgeComps_discard,
                                                                AgeingError = array(AgeingError_fish[,,,f,], dim = dim(AgeingError)),
                                                                comp_like = comp_fishage_discard_like,
                                                                ln_theta = ln_FishAge_discard_theta,
                                                                ln_theta_agg = ln_FishAge_discard_theta_agg,
                                                                corr_pars = FishAge_discard_corr_pars,
                                                                corr_pars_agg = FishAge_discard_corr_pars_agg,
                                                                comp_type = FishAgeComps_discard_Type,
                                                                n_sexes = n_sexes,
                                                                n_regions = n_regions,
                                                                n_cat = n_ages,
                                                                Obs = ObsFishAgeComps_discard,
                                                                pop_specific = FALSE,
                                                                age_or_len = 0)

              # Sample fishery ages (population specific, discard compositions)
              sim_env$ObsFishAgeComps_discard_pop <- simulate_comps(r = r,
                                                                    y = y,
                                                                    seas = seas,
                                                                    f = f,
                                                                    sim = sim,
                                                                    Exp = DAA,
                                                                    ISS_pop = ISS_FishAgeComps_discard_pop,
                                                                    AgeingError = array(AgeingError_fish[,,,f,], dim = dim(AgeingError)),
                                                                    pop_comp_like = comp_fishage_discard_pop_like,
                                                                    ln_pop_theta = ln_FishAge_discard_pop_theta,
                                                                    ln_pop_theta_agg = ln_FishAge_discard_pop_theta_agg,
                                                                    pop_corr_pars = FishAge_discard_pop_corr_pars,
                                                                    pop_corr_pars_agg = FishAge_discard_pop_corr_pars_agg,
                                                                    pop_comp_type = FishAgeComps_discard_pop_Type,
                                                                    n_sexes = n_sexes,
                                                                    n_regions = n_regions,
                                                                    n_pop = n_pop,
                                                                    n_cat = n_ages,
                                                                    Obs = ObsFishAgeComps_discard_pop,
                                                                    pop_specific = TRUE,
                                                                    age_or_len = 0)

              # Sample fishery lengths (retained compositions)
              if((exists("SizeAgeTrans") && !is.null(SizeAgeTrans)) || (exists("SizeAgeTrans_fish") && !is.null(SizeAgeTrans_fish))) {
                # Sample fishery ages (non-population specific, discard compositions)
                sim_env$ObsFishLenComps_discard <- simulate_comps(r = r,
                                                                  y = y,
                                                                  seas = seas,
                                                                  f = f,
                                                                  sim = sim,
                                                                  Exp = DAL,
                                                                  ISS = ISS_FishLenComps_discard,
                                                                  AgeingError = NULL,
                                                                  comp_like = comp_fishlen_discard_like,
                                                                  ln_theta = ln_FishLen_discard_theta,
                                                                  ln_theta_agg = ln_FishLen_discard_theta_agg,
                                                                  corr_pars = FishLen_discard_corr_pars,
                                                                  corr_pars_agg = FishLen_discard_corr_pars_agg,
                                                                  comp_type = FishLenComps_discard_Type,
                                                                  n_sexes = n_sexes,
                                                                  n_regions = n_regions,
                                                                  n_cat = n_lens,
                                                                  Obs = ObsFishLenComps_discard,
                                                                  pop_specific = FALSE,
                                                                  age_or_len = 1)

                # Sample fishery lengths (population specific, discard compositions)
                sim_env$ObsFishLenComps_discard_pop <- simulate_comps(r = r,
                                                                      y = y,
                                                                      seas = seas,
                                                                      f = f,
                                                                      sim = sim,
                                                                      Exp = DAL,
                                                                      ISS_pop = ISS_FishLenComps_discard_pop,
                                                                      AgeingError = NULL,
                                                                      pop_comp_like = comp_fishlen_discard_pop_like,
                                                                      ln_pop_theta = ln_FishLen_discard_pop_theta,
                                                                      ln_pop_theta_agg = ln_FishLen_discard_pop_theta_agg,
                                                                      pop_corr_pars = FishLen_discard_pop_corr_pars,
                                                                      pop_corr_pars_agg = FishLen_discard_pop_corr_pars_agg,
                                                                      pop_comp_type = FishLenComps_discard_pop_Type,
                                                                      n_sexes = n_sexes,
                                                                      n_regions = n_regions,
                                                                      n_pop = n_pop,
                                                                      n_cat = n_lens,
                                                                      Obs = ObsFishLenComps_discard_pop,
                                                                      pop_specific = TRUE,
                                                                      age_or_len = 1)

              } # end if size age transition if availiable

            } # end if dmr > 0
          } # end if Fmort > 0

        } # end f loop
      } # end r loop
    } # end seas loop

    # fleets reporting once a year keep the year's total in season one, with the observation
    # error applied to that total rather than to each season
    if(n_seas > 1) {

      catch_agg <- collapse_seas_obs(TrueCatch, ObsCatch, exp(ln_sigmaC), Catch_seas_Type,
                                     rep(0, n_fish_fleets), y, sim, n_seas, n_regions, n_fish_fleets)
      sim_env$TrueCatch <- catch_agg$true
      sim_env$ObsCatch <- catch_agg$obs

      catch_pop_agg <- collapse_seas_obs(sim_env$TrueCatch_pop, sim_env$ObsCatch_pop, exp(ln_sigmaC_pop),
                                         Catch_pop_seas_Type, rep(0, n_fish_fleets), y, sim, n_seas,
                                         n_regions, n_fish_fleets, pop = TRUE)
      sim_env$TrueCatch_pop <- catch_pop_agg$true
      sim_env$ObsCatch_pop <- catch_pop_agg$obs

      fidx_agg <- collapse_seas_obs(sim_env$TrueFishIdx, sim_env$ObsFishIdx, ObsFishIdx_SE,
                                    FishIdx_seas_Type, FishIdx_LikeType, y, sim, n_seas,
                                    n_regions, n_fish_fleets)
      sim_env$TrueFishIdx <- fidx_agg$true
      sim_env$ObsFishIdx <- fidx_agg$obs

      fidx_pop_agg <- collapse_seas_obs(sim_env$TrueFishIdx_pop, sim_env$ObsFishIdx_pop, ObsFishIdx_pop_SE,
                                        FishIdx_pop_seas_Type, FishIdx_LikeType, y, sim, n_seas,
                                        n_regions, n_fish_fleets, pop = TRUE)
      sim_env$TrueFishIdx_pop <- fidx_pop_agg$true
      sim_env$ObsFishIdx_pop <- fidx_pop_agg$obs

      # compositions are redrawn from the year's numbers at age rather than from one season's
      if(any(FishAgeComps_seas_Type == 1)) {
        CAA_yr <- collapse_seas_at_age(CAA, FishAgeComps_seas_Type, y, sim, n_seas)
        for(f in which(FishAgeComps_seas_Type == 1)) {
          for(r in 1:n_regions) {
            sim_env$ObsFishAgeComps[r,y,-1,,,f,sim] <- 0
            sim_env$ObsFishAgeComps <- simulate_comps(r = r, y = y, f = f, seas = 1, sim = sim,
                                                      Exp = CAA_yr,
                                                      ISS = ISS_FishAgeComps,
                                                      AgeingError = array(AgeingError_fish[,,,f,], dim = dim(AgeingError)),
                                                      comp_like = comp_fishage_like,
                                                      ln_theta = ln_FishAge_theta,
                                                      ln_theta_agg = ln_FishAge_theta_agg,
                                                      corr_pars = FishAge_corr_pars,
                                                      corr_pars_agg = FishAge_corr_pars_agg,
                                                      comp_type = FishAgeComps_Type,
                                                      n_sexes = n_sexes,
                                                      n_regions = n_regions,
                                                      n_cat = n_ages,
                                                      Obs = sim_env$ObsFishAgeComps,
                                                      age_or_len = 0)
          } # end r loop
        } # end f loop
      } # end if aggregated fishery age compositions

    } # end if more than one season
  })

}


#' Collapse a year's seasonal observations into one annual observation
#'
#' A data source set to \code{"aggSeas"} reports once a year rather than once a
#' season, so the operating model has to draw one observation from the year's
#' total rather than eleven from its parts. The total is written into season one
#' by convention and the other seasons are left at zero, which is where the
#' estimation model's \code{Use} array should mark it.
#'
#' The true quantity is summed first and the observation error is applied once to
#' that sum, so the error is on the annual total the way the observation is
#' reported, rather than on each season and then added up.
#'
#' @param true_arr,obs_arr Arrays \code{[region, year, season, fleet, sim]}, or
#'   with a leading population dim, of the true quantity and the observation.
#' @param se_arr Standard deviation array matching \code{true_arr} without the
#'   simulation dim.
#' @param seas_agg Integer vector, one per fleet.
#' @param like_type Integer vector, one per fleet, passed to
#'   \code{\link{draw_index_obs}}.
#' @param y,sim Year and replicate being drawn.
#' @param n_seas,n_regions,n_fleets Dimension sizes.
#' @param pop Logical, whether the arrays have a leading population dim.
#'
#' @return A list with the updated \code{true} and \code{obs} arrays.
#'
#' @keywords internal
collapse_seas_obs <- function(true_arr, obs_arr, se_arr, seas_agg, like_type,
                              y, sim, n_seas, n_regions, n_fleets, pop = FALSE) {

  if(!any(seas_agg == 1) || n_seas == 1) return(list(true = true_arr, obs = obs_arr))

  for(f in which(seas_agg == 1)) {
    for(r in 1:n_regions) {

      if(pop) {
        year_total <- rowSums(matrix(true_arr[,r,y,,f,sim], ncol = n_seas)) # one total per population
        true_arr[,r,y,,f,sim] <- 0
        true_arr[,r,y,1,f,sim] <- year_total
        obs_arr[,r,y,,f,sim] <- 0
        obs_arr[,r,y,1,f,sim] <- draw_index_obs(year_total, se_arr[,r,y,1,f], like_type[f])

      } else {
        year_total <- sum(true_arr[r,y,,f,sim])
        true_arr[r,y,,f,sim] <- 0
        true_arr[r,y,1,f,sim] <- year_total
        obs_arr[r,y,,f,sim] <- 0
        obs_arr[r,y,1,f,sim] <- draw_index_obs(year_total, se_arr[r,y,1,f], like_type[f])
      }

    } # end r loop
  } # end f loop

  list(true = true_arr, obs = obs_arr)

} # end collapse_seas_obs

#' Sum an at-age or at-length array over seasons into season one
#'
#' The counterpart of \code{\link{collapse_seas_obs}} for compositions, which are
#' drawn from the numbers behind them rather than from a total. Fleets left at
#' \code{"spltSeas"} are returned untouched.
#'
#' @param arr Array \code{[pop, region, year, season, bin, sex, fleet, sim]} of
#'   predicted numbers.
#' @param seas_agg Integer vector, one per fleet.
#' @param y,sim Year and replicate being drawn.
#' @param n_seas Number of seasons.
#'
#' @return \code{arr} with each aggregated fleet's year summed into season one.
#'
#' @keywords internal
collapse_seas_at_age <- function(arr, seas_agg, y, sim, n_seas) {

  if(!any(seas_agg == 1) || n_seas == 1) return(arr)

  for(f in which(seas_agg == 1)) {
    year_total <- apply(arr[,,y,,,,f,sim, drop = FALSE], c(1,2,5,6), sum) # pop, region, bin, sex
    arr[,,y,,,,f,sim] <- 0
    arr[,,y,1,,,f,sim] <- year_total
  } # end f loop

  arr

} # end collapse_seas_at_age

#' Generate survey indices and compositions in a simulation
#'
#' Computes survey index-at-age (\code{SrvIAA}) for all populations using
#' the mid-survey abundance formula \eqn{N \cdot s \cdot e^{-t_{\text{srv}} Z}},
#' derives index-at-length (\code{SrvIAL}) when a size-age transition matrix
#' is available, generates observed survey indices (with lognormal error) as
#' abundance, biomass, or the recruitment deviations depending on
#' \code{srv_idx_type}, and draws age and
#' length composition samples via \code{\link{simulate_comps}}.
#'
#' @param y Integer. Year index.
#' @param sim Integer. Simulation replicate index.
#' @param sim_env Simulation environment created by \code{\link{Setup_sim_env}}.
#'   Modified in place. The following elements are updated:
#'   \describe{
#'     \item{\code{SrvIAA}}{Survey index-at-age for all populations.}
#'     \item{\code{SrvIAL}}{Survey index-at-length if \code{SizeAgeTrans} is present.}
#'     \item{\code{TrueSrvIdx}, \code{ObsSrvIdx}}{Aggregated survey index values.}
#'     \item{\code{TrueSrvIdx_pop}, \code{ObsSrvIdx_pop}}{Population-specific survey index values.}
#'     \item{\code{ObsSrvAgeComps}, \code{ObsSrvAgeComps_pop}}{Observed survey age compositions.}
#'     \item{\code{ObsSrvLenComps}, \code{ObsSrvLenComps_pop}}{Observed survey length compositions if \code{SizeAgeTrans} is available.}
#'   }
#'
#' @details This function loops over seasons, regions, and survey fleets for all
#' populations and replicates. It computes mid-period abundance, applies survey
#' selectivity, calculates true survey indices (abundance or biomass), applies
#' lognormal observation error, and simulates age and length composition samples.
#' Population-specific compositions are also generated if requested.
#'
#' @return \code{invisible(NULL)}. All modifications are performed by reference
#'   within \code{sim_env}.
#'
#' @keywords internal
#' @seealso \code{\link{simulate_comps}}, \code{\link{Setup_sim_env}}
generate_survey_comp_idx <- function(y, sim, sim_env) {

  sim_env$y   <- y
  sim_env$sim <- sim

  with(sim_env, {

    for(seas in 1:n_seas) {

      # numbers at survey timing. a survey index is a snapshot inside the season, so continuous
      # movement needs partial propagation under the combined movement-mortality generator
      if(move_timing == 2) {
        SrvN_sim <- array(0, dim = c(n_pop, n_regions, n_ages, n_sexes, n_srv_fleets))
        for(p in 1:n_pop) for(sf in 1:n_srv_fleets) for(a in 1:n_ages) for(s in 1:n_sexes) {
          SrvN_sim[p,,a,s,sf] <- survey_state(NAA[p,,y,seas,a,s,sim], NULL, ZAA[p,,y,seas,a,s,sim],
                                              Mrate[p,,,y,seas,a,s,sim], seasdur[seas],
                                              t_srv[,seas,sf], move_timing, expm_nsub = expm_nsub)
        }
      }

      for(r in 1:n_regions) {
        for(sf in 1:n_srv_fleets) {

          for(p in 1:n_pop) {
            # Survey Ages Indexed (midpoint year)
            if(move_timing == 2) {
              sim_env$SrvIAA[p,r,y,seas,,,sf,sim] <- SrvN_sim[p,r,,,sf] * srv_sel[p,r,y,seas,,,sf,sim]
            } else {
              sim_env$SrvIAA[p,r,y,seas,,,sf,sim] <- NAA[p,r,y,seas,,,sim] * srv_sel[p,r,y,seas,,,sf,sim] * exp(-t_srv[r,seas,sf] * ZAA[p,r,y,seas,,,sim])
            }
            if((exists("SizeAgeTrans") && !is.null(SizeAgeTrans)) || (exists("SizeAgeTrans_srv") && !is.null(SizeAgeTrans_srv))) for(s in 1:n_sexes) sim_env$SrvIAL[p,r,y,seas,,s,sf,sim] <- (if(exists("SizeAgeTrans_srv") && !is.null(SizeAgeTrans_srv)) SizeAgeTrans_srv[p,r,y,seas,,,s,sf,sim] else SizeAgeTrans[p,r,y,seas,,,s,sim]) %*% SrvIAA[p,r,y,seas,,s,sf,sim] # Survey index at length
          } # end p loop

          # Survey Index - Regional
          if(srv_idx_type[sf] == 0) sim_env$TrueSrvIdx[r,y,seas,sf,sim] <- srv_q[r,y,sf,sim] * sum(SrvIAA[,r,y,seas,,,sf,sim]) # True Survey Index (abundance)
          if(srv_idx_type[sf] == 1) sim_env$TrueSrvIdx[r,y,seas,sf,sim] <- srv_q[r,y,sf,sim] * sum(SrvIAA[,r,y,seas,,,sf,sim] * WAA_srv[,r,y,seas,,,sf,sim]) # True Survey Index (biomass)
          if(srv_idx_type[sf] == 2) sim_env$TrueSrvIdx[r,y,seas,sf,sim] <- srv_q[r,y,sf,sim] * sum(ln_RecDevs[,r,y,sim]) # True Survey Index (recruitment deviations)

          # observed index. an mvn fleet takes its scale from the covariance's factor decomposition
          # rather than the SE array, with one factor draw shared across its series per replicate
          sidx_like <- if(exists("SrvIdx_LikeType")) SrvIdx_LikeType[sf] else 0
          if(sidx_like == 2) {
            if(is.na(srv_idx_u[sf,sim])) sim_env$srv_idx_u[sf,sim] <- stats::rnorm(1)
            sidx_fac <- resolve_idx_factor(srv_idx_mvn[[sf]], r, y, seas)
            sim_env$ObsSrvIdx[r,y,seas,sf,sim] <- draw_index_obs(
              TrueSrvIdx[r,y,seas,sf,sim],
              NA,
              2,
              d = sidx_fac$d,
              lambda = sidx_fac$lambda,
              u = srv_idx_u[sf,sim]
            )
          } else {
            sim_env$ObsSrvIdx[r,y,seas,sf,sim] <- draw_index_obs(TrueSrvIdx[r,y,seas,sf,sim], ObsSrvIdx_SE[r,y,seas,sf], sidx_like)
          }

          # survey index at age
          if(exists("use_srv_idx_aa") && use_srv_idx_aa[sf] == 1) {
            saa_dim <- c(n_pop, n_regions, n_ages, n_sexes)
            saa_obs_dim <- c(n_regions, n_obs_ages, n_sexes) # use flags and errors sit on the observed ages
            store_at_age_cell(sim_env, "SrvIdxAA",
                              sim_at_age_cell(array(SrvIAA[,,y,seas,,,sf,sim], dim = saa_dim),
                                              array(0, dim = saa_dim),
                                              array(UseSrvIdxAA[,y,seas,,,sf], dim = saa_obs_dim),
                                              array(ObsSrvIdxAA_SE[,y,seas,,,sf], dim = saa_obs_dim),
                                              array(ln_sigmaSrvIdxAA[,,sf], dim = c(n_obs_ages, n_sexes)),
                                              SrvIdxAA_Type[y,sf], SrvIdxAA_LikeType[sf], SrvIdxAA_sigma_form[sf],
                                              FALSE, r,
                                              ageing_error = array(AgeingError_srv[y,,,sf,sim], dim = dim(AgeingError_srv)[2:3])),
                              r, y, seas, sf, sim)
          }

          # population-specific index. the covariance describes the regional series only, so an mvn
          # fleet's population data source keeps lognormal error, mirroring the estimation model
          if(srv_idx_type[sf] == 0) sim_env$TrueSrvIdx_pop[,r,y,seas,sf,sim] <- srv_q[r,y,sf,sim] * apply(SrvIAA[,r,y,seas,,,sf,sim, drop = FALSE], 1, sum) # True Survey Index (abundance)
          if(srv_idx_type[sf] == 1) sim_env$TrueSrvIdx_pop[,r,y,seas,sf,sim] <- srv_q[r,y,sf,sim] * apply(SrvIAA[,r,y,seas,,,sf,sim, drop = FALSE] * WAA_srv[,r,y,seas,,,sf,sim, drop = FALSE], 1, sum) # True Survey Index (biomass)
          sim_env$ObsSrvIdx_pop[,r,y,seas,sf,sim] <- draw_index_obs(TrueSrvIdx_pop[,r,y,seas,sf,sim], ObsSrvIdx_pop_SE[,r,y,seas,sf], if(sidx_like == 1) 1 else 0)

          # Survey Compositions
          # Sample survey ages
          sim_env$ObsSrvAgeComps <- simulate_comps(r = r,
                                                   y = y,
                                                   f = sf,
                                                   seas = seas,
                                                   sim = sim,
                                                   Exp = SrvIAA,
                                                   ISS = ISS_SrvAgeComps,
                                                   AgeingError = array(AgeingError_srv[,,,sf,], dim = dim(AgeingError)),
                                                   comp_like = comp_srvage_like,
                                                   ln_theta = ln_SrvAge_theta,
                                                   ln_theta_agg = ln_SrvAge_theta_agg,
                                                   corr_pars = SrvAge_corr_pars,
                                                   corr_pars_agg = SrvAge_corr_pars_agg,
                                                   comp_type = SrvAgeComps_Type,
                                                   n_sexes = n_sexes,
                                                   n_regions = n_regions,
                                                   n_cat = n_ages,
                                                   Obs = ObsSrvAgeComps,
                                                   age_or_len = 0)

          # Sample survey ages (population specific)
          sim_env$ObsSrvAgeComps_pop <- simulate_comps(r = r,
                                                       y = y,
                                                       seas = seas,
                                                       f = sf,
                                                       sim = sim,
                                                       Exp = SrvIAA,
                                                       ISS_pop = ISS_SrvAgeComps_pop,
                                                       AgeingError = array(AgeingError_srv[,,,sf,], dim = dim(AgeingError)),
                                                       pop_comp_like = comp_srvage_pop_like,
                                                       ln_pop_theta = ln_SrvAge_pop_theta,
                                                       ln_pop_theta_agg = ln_SrvAge_pop_theta_agg,
                                                       pop_corr_pars = SrvAge_pop_corr_pars,
                                                       pop_corr_pars_agg = SrvAge_pop_corr_pars_agg,
                                                       pop_comp_type = SrvAgeComps_pop_Type,
                                                       n_sexes = n_sexes,
                                                       n_regions = n_regions,
                                                       n_pop = n_pop,
                                                       n_cat = n_ages,
                                                       Obs = ObsSrvAgeComps_pop,
                                                       pop_specific = TRUE,
                                                       age_or_len = 0)

          # Sample survey lengths
          if((exists("SizeAgeTrans") && !is.null(SizeAgeTrans)) || (exists("SizeAgeTrans_srv") && !is.null(SizeAgeTrans_srv))) {
            sim_env$ObsSrvLenComps <- simulate_comps(r = r,
                                                     y = y,
                                                     f = sf,
                                                     seas = seas,
                                                     sim = sim,
                                                     Exp = SrvIAL,
                                                     ISS = ISS_SrvLenComps,
                                                     AgeingError = NULL,
                                                     comp_like = comp_srvlen_like,
                                                     ln_theta = ln_SrvLen_theta,
                                                     ln_theta_agg = ln_SrvLen_theta_agg,
                                                     corr_pars = SrvLen_corr_pars,
                                                     corr_pars_agg = SrvLen_corr_pars_agg,
                                                     comp_type = SrvLenComps_Type,
                                                     n_sexes = n_sexes,
                                                     n_regions = n_regions,
                                                     n_cat = n_lens,
                                                     Obs = ObsSrvLenComps,
                                                     age_or_len = 1)

            # Sample survey lengths (population specific)
            sim_env$ObsSrvLenComps_pop <- simulate_comps(r = r,
                                                         y = y,
                                                         seas = seas,
                                                         f = sf,
                                                         sim = sim,
                                                         Exp = SrvIAL,
                                                         ISS_pop = ISS_SrvLenComps_pop,
                                                         AgeingError = NULL,
                                                         pop_comp_like = comp_srvlen_pop_like,
                                                         ln_pop_theta = ln_SrvLen_pop_theta,
                                                         ln_pop_theta_agg = ln_SrvLen_pop_theta_agg,
                                                         pop_corr_pars = SrvLen_pop_corr_pars,
                                                         pop_corr_pars_agg = SrvLen_pop_corr_pars_agg,
                                                         pop_comp_type = SrvLenComps_pop_Type,
                                                         n_sexes = n_sexes,
                                                         n_regions = n_regions,
                                                         n_pop = n_pop,
                                                         n_cat = n_lens,
                                                         Obs = ObsSrvLenComps_pop,
                                                         pop_specific = TRUE,
                                                         age_or_len = 1)

            # Sample survey conditional age-at-length
            if(exists("do_srv_caal") && isTRUE(do_srv_caal)) {
              sim_env$ObsSrv_caal <- simulate_caal(
                r = r,
                y = y,
                f = sf,
                seas = seas,
                sim = sim,
                SizeAgeTrans = if(exists("SizeAgeTrans_srv") && !is.null(SizeAgeTrans_srv)) array(SizeAgeTrans_srv[,,,,,,,sf,], dim = dim(SizeAgeTrans_srv)[-8]) else SizeAgeTrans,
                AtAge = SrvIAA,
                ISS = ISS_Srv_caal,
                AgeingError = array(AgeingError_srv[,,,sf,], dim = dim(AgeingError)),
                comp_like = comp_srv_caal_like,
                ln_theta = ln_Srv_caal_theta,
                ln_theta_agg = ln_Srv_caal_theta_agg,
                comp_type = Srv_caal_Type,
                n_sexes = n_sexes,
                n_regions = n_regions,
                n_lens = n_lens,
                Obs = ObsSrv_caal
              )
            } # end survey caal

          } # end if size age transition if availiable

        } # end sf loop
      } # end r loop
    } # end seas loop

    # surveys reporting once a year keep the year's total in season one
    if(n_seas > 1) {

      sidx_agg <- collapse_seas_obs(sim_env$TrueSrvIdx, sim_env$ObsSrvIdx, ObsSrvIdx_SE,
                                    SrvIdx_seas_Type, SrvIdx_LikeType, y, sim, n_seas,
                                    n_regions, n_srv_fleets)
      sim_env$TrueSrvIdx <- sidx_agg$true
      sim_env$ObsSrvIdx <- sidx_agg$obs

      sidx_pop_agg <- collapse_seas_obs(sim_env$TrueSrvIdx_pop, sim_env$ObsSrvIdx_pop, ObsSrvIdx_pop_SE,
                                        SrvIdx_pop_seas_Type, SrvIdx_LikeType, y, sim, n_seas,
                                        n_regions, n_srv_fleets, pop = TRUE)
      sim_env$TrueSrvIdx_pop <- sidx_pop_agg$true
      sim_env$ObsSrvIdx_pop <- sidx_pop_agg$obs

      # compositions are redrawn from the year's numbers at age rather than from one season's
      if(any(SrvAgeComps_seas_Type == 1)) {
        SrvIAA_yr <- collapse_seas_at_age(SrvIAA, SrvAgeComps_seas_Type, y, sim, n_seas)
        for(sf in which(SrvAgeComps_seas_Type == 1)) {
          for(r in 1:n_regions) {
            sim_env$ObsSrvAgeComps[r,y,-1,,,sf,sim] <- 0
            sim_env$ObsSrvAgeComps <- simulate_comps(r = r, y = y, f = sf, seas = 1, sim = sim,
                                                     Exp = SrvIAA_yr,
                                                     ISS = ISS_SrvAgeComps,
                                                     AgeingError = array(AgeingError_srv[,,,sf,], dim = dim(AgeingError)),
                                                     comp_like = comp_srvage_like,
                                                     ln_theta = ln_SrvAge_theta,
                                                     ln_theta_agg = ln_SrvAge_theta_agg,
                                                     corr_pars = SrvAge_corr_pars,
                                                     corr_pars_agg = SrvAge_corr_pars_agg,
                                                     comp_type = SrvAgeComps_Type,
                                                     n_sexes = n_sexes,
                                                     n_regions = n_regions,
                                                     n_cat = n_ages,
                                                     Obs = sim_env$ObsSrvAgeComps,
                                                     age_or_len = 0)
          } # end r loop
        } # end sf loop
      } # end if aggregated survey age compositions

    } # end if more than one season

  })
}

#' Release conventional tags in the simulation
#'
#' Spreads each cohort's \code{n_tags}, or \code{n_tags_rel_input} when supplied,
#' across populations, ages and sexes in proportion to the selectivity-weighted
#' abundance \code{NAA_bef} of the release platform, rounding to integers.
#' \code{conv_fish_tag_attr} is then applied through
#' \code{\link{marginalize_conv_fish_tags}} to give the observation-level release
#' array \code{conv_tagged_fish_attr}, at the resolution the recapture likelihood
#' reads.
#'
#' Without \code{n_tags_rel_input}, a survey or fishery platform scales the
#' region's total against the selectivity-weighted global abundance, and a
#' population platform against the region's share of total abundance.
#'
#' @param y Integer. Year index.
#' @param sim Integer. Simulation replicate index.
#' @param sim_env Simulation environment from \code{\link{Setup_sim_env}},
#'   modified in place: \code{$conv_tagged_fish} and
#'   \code{$conv_tagged_fish_attr} for each cohort released in year \code{y}.
#'
#' @return \code{invisible(NULL)}; everything is modified by reference within
#'   \code{sim_env}.
#'
#' @keywords internal
release_conv_tags <- function(y, sim, sim_env) {

  sim_env$y   <- y
  sim_env$sim <- sim

  with(sim_env, {
    for(seas in 1:n_seas) {
      for(r in 1:n_regions) {

        # Get indices for tag cohorts in the current year and region
        tag_rel <- which(conv_tag_release_indicator[,1] == r & conv_tag_release_indicator[,2] == y & conv_tag_release_indicator[,3] == seas) # Get tag cohort (release event)

        # Release Tags if any events
        if(length(tag_rel) != 0) {

          # Tag Indexing
          tr <- conv_tag_release_indicator[tag_rel,1] # tag release region
          ty <- conv_tag_release_indicator[tag_rel,2] # tag release year
          tseas <- conv_tag_release_indicator[tag_rel,3] # tag release season
          tplat <- conv_tag_release_platform[tag_rel, 1] # get tagging platform information
          tplat_f <- as.numeric(conv_tag_release_platform[tag_rel, 2]) # get tagging fleet

          # Survey
          # Survey
          if(tplat[1] == 'survey') {
            NAA_slice  <- NAA_bef[, tr, ty, tseas, , , sim, drop = FALSE]
            dim(NAA_slice) <- c(n_pop, n_ages, n_sexes)
            NAA_sel <- array(0, dim = c(n_pop, n_ages, n_sexes))
            for(p in 1:n_pop) {
              sel_slice <- array(srv_sel[p, tr, ty, tseas, , , tplat_f, sim], dim = c(n_ages, n_sexes))
              NAA_sel[p,,] <- NAA_slice[p,,] * sel_slice
            }
            if(!exists("n_tags_rel_input")) {
              denom <- 0
              for(p in 1:n_pop) {
                NAA_denom_p <- array(NAA_bef[p, , ty, tseas, , , sim], dim = c(n_regions, n_ages, n_sexes))
                for(rr in 1:n_regions) {
                  sel_d <- array(srv_sel[p, rr, ty, tseas, , , tplat_f, sim], dim = c(n_ages, n_sexes))
                  denom <- denom + sum(NAA_denom_p[rr,,] * sel_d)
                }
              }
              n_tags_rel <- round(sum(NAA_sel) / denom * n_tags)
            } else {
              n_tags_rel <- n_tags_rel_input[tag_rel]
            }
            tmp_props <- NAA_sel / sum(NAA_sel)
          }

          # distribute tags for fishery
          if(tplat[1] == 'fishery') {
            NAA_slice  <- NAA_bef[, tr, ty, tseas, , , sim, drop = FALSE]
            dim(NAA_slice) <- c(n_pop, n_ages, n_sexes)
            NAA_sel <- array(0, dim = c(n_pop, n_ages, n_sexes))
            for(p in 1:n_pop) {
              sel_slice <- array(fish_sel[p, tr, ty, tseas, , , tplat_f, sim], dim = c(n_ages, n_sexes))
              NAA_sel[p,,] <- NAA_slice[p,,] * sel_slice
            }
            if(!exists("n_tags_rel_input")) {
              denom <- 0
              for(p in 1:n_pop) {
                NAA_denom_p <- array(NAA_bef[p, , ty, tseas, , , sim], dim = c(n_regions, n_ages, n_sexes))
                for(rr in 1:n_regions) {
                  sel_d <- array(fish_sel[p, rr, ty, tseas, , , tplat_f, sim], dim = c(n_ages, n_sexes))
                  denom <- denom + sum(NAA_denom_p[rr,,] * sel_d)
                }
              }
              n_tags_rel <- round(sum(NAA_sel) / denom * n_tags)
            } else {
              n_tags_rel <- n_tags_rel_input[tag_rel]
            }
            tmp_props <- NAA_sel / sum(NAA_sel)
          }

          # distribute tags by population
          if(tplat[1] == 'population') {
            if(!exists("n_tags_rel_input")) {
              n_tags_rel <- round(sum(NAA_bef[,tr,ty,tseas,,,sim]) / sum(NAA_bef[,,ty,tseas,,,sim]) * n_tags) # get region specific tags
            } else {
              n_tags_rel <- n_tags_rel_input[tag_rel] # use input tags by cohort if availiable
            }
            tmp_props <- NAA_bef[, tr, ty, tseas, , , sim] / sum(NAA_bef[, tr, ty, tseas, , , sim]) # get proportions by population, age, and sex
          }

          # multiply by tags in each region and distributa cross ages and sexes
          sim_env$conv_tagged_fish[tag_rel, , , , sim] <- array(round(tmp_props * n_tags_rel), dim = c(n_pop, n_ages, n_sexes))

          # marginalize across appropriate dimensions of what recaptures attributes are and report out
          sim_env$conv_tagged_fish_attr[tag_rel, , , , sim] <- marginalize_conv_fish_tags(conv_tagged_fish[tag_rel, , , , sim],
                                                                                              conv_fish_tag_attr[tag_rel], n_pop, n_ages, n_sexes)

        } # end if no tag releases

      } # end r loop
    } # end seas loop
  })
}

#' Generate conventional tag recaptures from fisheries in simulation
#'
#' For each tag cohort and recovery season in year \code{y}, advances the
#' available tagged fish through movement and mortality, applies Baranov's
#' equation for the predicted recaptures, and draws the observed ones through
#' \code{\link{simulate_conv_tag_fish_recaptures}}. Cohorts not yet released,
#' already at \code{conv_tag_max_liberty}, or released in a future year are
#' skipped.
#'
#' @param y Integer. Year index.
#' @param sim Integer. Simulation replicate index.
#' @param sim_env Simulation environment from \code{\link{Setup_sim_env}},
#'   modified in place. It gains \code{conv_tag_fish_avail}, the tagged fish
#'   available by age, region, season and fleet for each cohort;
#'   \code{pred_conv_tag_fish_recap} and \code{obs_conv_tag_fish_recap}, the
#'   predicted recaptures and the observed ones after sampling error;
#'   \code{conv_tag_fish_surv}, the survivors after mortality and movement; and
#'   \code{conv_tag_fish_reported}, the reporting-adjusted counts by fleet and
#'   region.
#'
#' @return \code{invisible(NULL)}; everything is modified by reference within
#'   \code{sim_env}.
#'
#' @details
#' Total fishing mortality entering \eqn{Z} splits into retained,
#' \eqn{F \cdot s_{\text{fish}} \cdot s_{\text{ret}}}, and dead discards,
#' \eqn{F \cdot s_{\text{fish}} \cdot (1 - s_{\text{ret}}) \cdot \text{dmr}}, as in
#' \code{\link{apply_pop_dy}}, and the Baranov numerator takes the retained
#' component alone, since tags come back from retained catch.
#'
#' At release, tags enter \code{conv_tag_fish_avail[1, rseas, tc, ...]} discounted
#' by the initial tag-induced mortality \code{ln_init_conv_tag_mort[tc]}, and when
#' \code{conv_tag_t_tagging[tc] < 1} total mortality is scaled by the fraction of
#' the season remaining, for that cell alone. Chronic shedding
#' \code{ln_conv_tag_shed[tc]} joins natural and fishing mortality in the total
#' rate. All three are vectors of length \code{n_tag_rel_events} indexed by
#' release event, so timing, initial mortality and shedding may differ across
#' cohorts. At the end of a season the survivors advance to the next season, or to
#' the next year's first season with plus group accumulation, and the reporting
#' rates in \code{conv_tag_fish_reporting} are applied by fleet and region.
#'
#' @keywords internal
generate_fishery_conv_tags_recap <- function(y, sim, sim_env) {

  sim_env$y   <- y
  sim_env$sim <- sim

  with(sim_env,{

    for(rseas in 1:n_seas) {

      # mortality is the same for every cohort at liberty in this year and season, so work it out
      # once here. the simulation dim is dropped so the arrays match the estimation model's
      tag_mort <- get_tag_mort(
        y = y,
        rseas = rseas,
        n_pop = n_pop,
        n_regions = n_regions,
        n_ages = n_ages,
        n_sexes = n_sexes,
        n_fish_fleets = n_fish_fleets,
        use_conv_fish_tagging = use_conv_fish_tagging,
        Fmort = array(Fmort[,,,,sim], dim = dim(Fmort)[1:4]),
        fish_sel = array(fish_sel[,,,,,,,sim], dim = dim(fish_sel)[1:7]),
        ret_sel = array(ret_sel[,,,,,,,sim], dim = dim(ret_sel)[1:7]),
        dmr = array(dmr[,,,,sim], dim = dim(dmr)[1:4]),
        natmort = array(natmort[,,,,,,sim], dim = dim(natmort)[1:6]),
        seasdur = seasdur
      )

      for(tc in 1:n_tag_rel_events) {

        # get indexing
        tr <- conv_tag_release_indicator[tc,1] # tag release region
        ty <- conv_tag_release_indicator[tc,2] # tag release year
        tseas <- conv_tag_release_indicator[tc,3] # tag release seasons

        # Skipping stuff if hasn't occurred yet, or if max liberty
        if(y < ty || (y == ty && rseas < tseas)) next
        ry <- y - ty + 1 # get tag liberty
        if(ry > conv_tag_max_liberty) next # skip if max liberty

        # Cohort specific containers, with in what earlier years at liberty
        # already recorded for this cohort
        avail_tc <- array(conv_tag_fish_avail[, , tc, , , , , sim],
                          dim = c(conv_tag_max_liberty + 1, n_seas, n_pop, n_regions, n_ages, n_sexes))
        recap_tc <- array(pred_conv_tag_fish_recap[, , tc, , , , , , sim],
                          dim = c(conv_tag_max_liberty, n_seas, n_pop, n_regions, n_ages, n_sexes, n_fish_fleets))

        # get fishing and natural mortality
        tmp_FAA <- tag_mort$FAA
        tmp_ret_FAA <- tag_mort$ret_FAA
        tmp_disc_DAA <- tag_mort$disc_DAA

        # get total mortality, adding this cohort's shedding rate
        tmp_ZAA <- tag_mort$Z_before_shed + (exp(ln_conv_tag_shed[tc]) * seasdur[rseas])

        # Fraction of this season the tag cohort is at liberty for
        tag_frac <- if(ry == 1 && rseas == tseas) conv_tag_t_tagging[tc] else 1
        tag_dur <- seasdur[rseas] * tag_frac

        # discount by tagging time when tagging is not at the start of the season. must match
        # get_tagging_observation_model(): scale F and Z together so F/Z stays a fraction
        if(tag_frac != 1) {
          tmp_ZAA      <- tmp_ZAA      * tag_frac
          tmp_FAA      <- tmp_FAA      * tag_frac
          tmp_ret_FAA  <- tmp_ret_FAA  * tag_frac
          tmp_disc_DAA <- tmp_disc_DAA * tag_frac
        }

        if(ry == 1 && rseas == tseas) {
          # Input tagged fish into available tags for recapture and adjust initial number of tagged fish for tag induced mortality (exponential mortality process)
          avail_tc[1, rseas, , tr, , ] <- array(conv_tagged_fish[tc, , , , sim] * exp(-exp(ln_init_conv_tag_mort[tc])), dim = c(n_pop, n_ages, n_sexes))
        }

        # get temporary survival value
        tmp_SAA <- exp(-tmp_ZAA)
        tag_moves <- (conv_tag_t_tagging[tc] == 1 || ry != 1 || rseas != tseas)

        # Move tagged fish around (skip only in first release year + tagging season when tagging occurs mid-season).
        # Under move_timing 1 and 2 movement is set by the transition operator below instead.
        if(move_timing == 0 && tag_moves) {
          for(p in 1:n_pop) {
            # Movement of tag cohorts
            if(do_recruits_move == 0) {
              for(a in 2:n_ages) for(s in 1:n_sexes) {
                avail_tc[ry, rseas, p, , a, s] <-
                  t(avail_tc[ry, rseas, p, , a, s]) %*%
                  Movement[p, , , y, rseas, a, s, sim]
              }
            } else { # if recruits move
              for(a in 1:n_ages) for(s in 1:n_sexes) {
                avail_tc[ry, rseas, p, , a, s] <-
                  t(avail_tc[ry, rseas, p, , a, s]) %*%
                  Movement[p, , , y, rseas, a, s, sim]
              } # end s loop
            } # end else
          } # end p loop
        } # end if

        # Post-season tag numbers, before the ageing shift
        if(move_timing == 0 || n_regions == 1) {
          tag_step <- array(avail_tc[ry, rseas, , , , ] * tmp_SAA[,,1,,],
                            dim = c(n_pop, n_regions, n_ages, n_sexes))
        } else {
          tag_step <- array(0, dim = c(n_pop, n_regions, n_ages, n_sexes))
          for(p in 1:n_pop) {
            for(a in 1:n_ages) {
              moves <- tag_moves && (do_recruits_move == 1 || a > 1)
              for(s in 1:n_sexes) {
                Mv <- if(moves) Movement[p,,,y,rseas,a,s,sim] else diag(n_regions)
                Qv <- if(moves) Mrate[p,,,y,rseas,a,s,sim] else matrix(0, n_regions, n_regions)
                tag_step[p,,a,s] <- advance_seas(avail_tc[ry,rseas,p,,a,s], Mv,
                                                 tmp_ZAA[p,,1,a,s], Qv, tag_dur, move_timing, expm_nsub = expm_nsub)
              } # end s loop
            } # end a loop
          } # end p loop
        }

        # Apply mortality and ageing to tagged fish
        if(rseas < n_seas) {

          # Season mortality within a given year, advance to next season same year/age
          avail_tc[ry, rseas + 1, , , , ] <- tag_step

        } else {

          # End of year mortality and age advancement (end of season)
          avail_tc[ry + 1, 1, , , 2:n_ages, ] <- tag_step[,,1:(n_ages - 1),]

          # Accumulate plus group
          avail_tc[ry + 1, 1, , , n_ages, ] <-
            avail_tc[ry + 1, 1, , , n_ages, ] + tag_step[,,n_ages,]
        }

        # # Apply Baranov's to get predicted recaptures
        # (add tiny epsilon to avoid 0/0 when tmp_ZAA == 0, e.g. conv_tag_t_tagging == 0 at release)
        for(f in 1:n_fish_fleets) {
          for(p in 1:n_pop) {
            if(move_timing == 2) {
              # Spatial Baranov: tags redistribute among regions while being caught, so
              # recaptures use the season-integrated tag abundance
              tag_int <- array(0, dim = c(n_regions, n_ages, n_sexes))
              for(a in 1:n_ages) {
                moves <- tag_moves && (do_recruits_move == 1 || a > 1)
                for(s in 1:n_sexes) {
                  Qv <- if(moves) Mrate[p,,,y,rseas,a,s,sim] else matrix(0, n_regions, n_regions)
                  tag_int[,a,s] <- integrate_seas_abundance(avail_tc[ry,rseas,p,,a,s],
                                                            tmp_ZAA[p,,1,a,s], Qv, tag_dur, expm_nsub = expm_nsub)
                } # end s loop
              } # end a loop
              # array() guards against R dropping a length-1 sex dimension from the F slice
              tmp_ret_FAA_slice <- array(tmp_ret_FAA[p,,1,,,f], dim = c(n_regions, n_ages, n_sexes))
              recap_tc[ry,rseas,p,,,,f] <- conv_tag_fish_reporting[,y,f,sim] *
                tmp_ret_FAA_slice * tag_int
            } else {
              recap_tc[ry,rseas,p,,,,f] <- conv_tag_fish_reporting[,y,f,sim] *
                (tmp_ret_FAA[p,,1,,,f] / (tmp_ZAA[p,,1,,] + 1e-10)) *
                avail_tc[ry,rseas,p,,,] *
                (1 - tmp_SAA[p,,1,,])
            }
          } # end p loop
        } # end f loop

        # Store this cohort's tags and predicted recaptures, which the recapture
        # draw below reads
        sim_env$conv_tag_fish_avail[, , tc, , , , , sim] <- avail_tc
        sim_env$pred_conv_tag_fish_recap[, , tc, , , , , , sim] <- recap_tc

        # Simulate Tag Recoveries
        sim_env$obs_conv_tag_fish_recap <- simulate_conv_tag_fish_recaptures(
          conv_fish_tag_like = conv_fish_tag_like,
          tag_recaptures_attr = conv_fish_tag_attr[tc],
          conv_tagged_fish = conv_tagged_fish,
          pred_conv_tag_fish_recap = pred_conv_tag_fish_recap,
          obs_conv_tag_fish_recap = obs_conv_tag_fish_recap,
          ln_conv_fish_tag_theta = ln_conv_fish_tag_theta,
          ry = ry,
          rseas = rseas,
          tc = tc,
          sim = sim,
          n_pop = n_pop,
          n_regions = n_regions,
          n_ages = n_ages,
          n_sexes = n_sexes,
          n_fish_fleets = n_fish_fleets
        )

      } # end tc loop
    } # end rseas loop
  })
}

# Sample Size Prediction ----------------------------------------------------

#' Predict fishery and discarded ISS under projected fishing mortality
#'
#' Scales fishery input sample sizes for the projection year \code{y} based
#' on the relationship between fishing mortality and historical ISS values.
#' For each region-sex-fleet cell, the minimum and maximum ISS from the
#' conditioning period (\code{1:(y-1)}) are identified from years with
#' positive, non-NA values, and the projected ISS is obtained by linear
#' interpolation between those bounds using the ratio of projected
#' \eqn{F_y} to the historical maximum \eqn{F} (capped at 1). If no valid
#' historical observations exist for a cell, ISS is set to zero. If
#' conditions for scaling are not met (e.g., maximum historical \eqn{F = 0}),
#' the mean historical ISS is used as a fallback. All prior years
#' (\code{1:(y-1)}) are reused unchanged from \code{ISS_FishComps}.
#'
#' @param ISS_FishComps Array of fishery ISS values
#'   \code{[n_regions × n_yrs × n_seas × n_sexes × n_fish_fleets × n_sims]}.
#' @param Fmort Array of fishing mortality rates
#'   \code{[n_regions × n_yrs × n_seas × n_fish_fleets × n_sims]}.
#' @param y Integer. Projection year index for which ISS is predicted.
#' @param sim Integer. Simulation replicate index.
#' @param seas Integer. Season index.
#'
#' @return Array \code{[n_regions × y × 1 × n_sexes × n_fish_fleets]} with
#'   historical ISS values filled in for years \code{1:(y-1)} and the
#'   predicted ISS in year \code{y}.
#'
#'
#' @keywords internal
predict_sim_fish_iss_fmort <- function(ISS_FishComps,
                                       Fmort,
                                       y,
                                       seas,
                                       sim
                                       ) {

  # dimensions
  dims <- dim(ISS_FishComps)
  n_regions <- dims[1]
  n_seas <- dims[3]
  n_sexes <- dims[4]
  n_fish_fleets <- dims[5]

  # extract temp vars
  tmp_iss <- ISS_FishComps[, 1:(y - 1), seas , , , sim, drop = FALSE]
  tmp_fmort <- Fmort[, 1:(y - 1), seas ,, sim, drop = FALSE]

  # container
  iss_container <- array(0, dim = c(n_regions, length(1:y), 1, n_sexes, n_fish_fleets))
  iss_container[, 1:(y - 1), seas, , ] <- ISS_FishComps[, 1:(y - 1), seas , , , sim] # fill in values back

  for(r in 1:n_regions) {
    for(s in 1:n_sexes) {
      for(f in 1:n_fish_fleets) {
        # get ISS and Fmort
        iss_vec <- tmp_iss[r, , 1, s, f, ]
        fmort_vec <- tmp_fmort[r, , 1, f, ]
        # remove zeros/NAs
        valid_idx <- which(iss_vec > 0 & !is.na(iss_vec) & !is.na(fmort_vec))
        if(length(valid_idx) > 0) {
          iss_valid <- iss_vec[valid_idx]
          fmort_valid <- fmort_vec[valid_idx]
          # min and max ISS from conditioning period
          min_iss <- min(iss_valid)
          max_iss <- max(iss_valid)
          # max Fmort from conditioning period
          max_fmort_hist <- max(fmort_valid)
          # new Fmort
          fmort_new <- Fmort[r, y, seas, f, sim]
          # scale ISS proportionally to Fmort relative to historical max
          if(max_fmort_hist > 0 && fmort_new >= 0) {
            # linear scaling between min and max ISS
            scaling_factor <- min(fmort_new / max_fmort_hist, 1)  # cap scaling at 1
            new_iss <- min_iss + scaling_factor * (max_iss - min_iss) # linear scaling
          } else {
            # use mean ISS if conditions not met ...
            new_iss <- mean(iss_valid)
          }
          iss_container[r, y, 1, s, f] <- new_iss
        } else {
          iss_container[r, y, 1, s, f] <- 0
        }
      } # end f loop
    } # end s loop
  } # end r loop

  return(iss_container)
}
