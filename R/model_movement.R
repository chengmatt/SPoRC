# Stage 2 of 3: objective function
#
# Turns movement parameters into the movement arrays the dynamics use, for both
# the unstructured and the continuous time parameterizations.
# get_movement_dp_design_matrix and get_ctmc_bound_form are shared with setup_movement.R.

#' Names for the CTMC generator bound forms
#'
#' @param x The \code{ctmc_diffusion_bounds} value: the code \code{0}, \code{1} or
#'   \code{2}, or the matching name \code{"none"}, \code{"softplus"} or
#'   \code{"upwind"}.
#'
#' @return Character scalar naming the form, or \code{NA_character_} if \code{x} is
#'   not one of the accepted values.
#'
#' @keywords internal
get_ctmc_bound_form <- function(x) {
  forms <- c("none", "softplus", "upwind") # codes 0, 1 and 2 in this order
  if(is.numeric(x) && length(x) == 1 && x %in% (seq_along(forms) - 1)) return(forms[x + 1])
  if(is.character(x) && length(x) == 1 && x %in% forms) return(x)
  NA_character_
}

#' Get Design Matrices for CTMC Movement
#'
#' Constructs the design matrices for the diffusion and preference components of a
#' Continuous Time Markov Chain (CTMC) movement model. These matrices are used
#' to parameterize movement rates in terms of covariates specified by formulas.
#'
#' @param data A \code{data.frame} containing the CTMC covariates. Must include all variables
#'   referenced in \code{diffusion_formula} and \code{preference_formula}.
#' @param preference_formula An R formula describing the linear predictor for movement
#'   preference (taxis). Variables must exist in \code{data}.
#' @param diffusion_formula An R formula describing the linear predictor for diffusion
#'   rates. Variables must exist in \code{data}.
#'
#' @return A \code{list} with the following components:
#' \describe{
#'   \item{\code{n_theta}}{Number of diffusion parameters (columns in \code{X_zk}).}
#'   \item{\code{n_gamma}}{Number of preference parameters (columns in \code{W_zk}).}
#'   \item{\code{X_zk}}{Diffusion design matrix constructed from \code{diffusion_formula} and \code{data}.}
#'   \item{\code{W_zk}}{Preference (taxis) design matrix constructed from \code{preference_formula} and \code{data}.}
#' }
#'
#' @keywords internal
get_movement_dp_design_matrix <- function(data,
                                          preference_formula,
                                          diffusion_formula
) {
  X_zk = stats::model.matrix(diffusion_formula, data) # diffusion design matrix
  W_zk = stats::model.matrix(preference_formula, data) # preference design matrix
  return(list(
    n_theta = ncol(X_zk),
    n_gamma = ncol(W_zk),
    X_zk = X_zk,
    W_zk = W_zk
  ))
}

#' Construct Movement Matrices for Unstructured or CTMC Movement
#'
#' Generates movement matrices for a population model based on either unstructured
#' multinomial logit movement or a Continuous Time Markov Chain (CTMC) formulation.
#' Also calculates a movement penalty if applicable. For CTMC movement, projection
#' years are supported: covariate lookups are capped at the last historical year
#' (\code{n_yrs}), so CTMC parameters are frozen at their final historical values
#' during projection unless \code{ctmc_move_dat} is extended with projection-year rows.
#'
#' @param move_type Integer flag indicating movement type: 0 = unstructured Markov,
#'   1 = CTMC movement.
#' @param do_recruits_move Integer flag: 0 = recruits (age 1) do not move, 1 = recruits move.
#' @param n_regions Number of spatial regions.
#' @param n_yrs Number of years in the observed (historical) data.
#' @param n_proj_yrs_devs Number of projected years. Extends the year dimension of the
#'   movement array beyond \code{n_yrs}; CTMC covariate lookups are capped at \code{n_yrs}
#'   unless \code{ctmc_move_dat} contains projection-year rows.
#' @param n_ages Number of age classes.
#' @param n_sexes Number of sexes.
#' @param n_seas Number of seasons.
#' @param n_pop Number of populations.
#' @param move_pars Array of movement parameters for unstructured (multinomial logit)
#'   movement, dimensioned \code{[pop, from_region, counter, year, seas, age, sex]}.
#'   \code{counter} indexes the \code{n_regions - 1} non-reference destinations.
#'   Ignored when \code{move_type == 1} or \code{use_fixed_movement == 1}.
#' @param move_devs Array of origin-destination movement deviations, dimensioned
#'   \code{[pop, origin_region, counter, year, seas, age, sex]}, where \code{counter}
#'   indexes adjacent non-diagonal destinations from each origin region. Applied as
#'   multiplicative offsets (\code{exp(move_devs)}) to off-diagonal diffusion rates
#'   in CTMC movement, or as additive offsets to logit-scale parameters in unstructured
#'   movement. Always indexed on the actual (possibly projected) year.
#' @param use_fixed_movement Integer flag: 0 = estimate movement, 1 = use fixed matrix.
#' @param Fixed_Movement Fixed movement matrix used when \code{use_fixed_movement == 1},
#'   dimensioned \code{[pop, from_region, to_region, year, seas, age, sex]}. Ignored otherwise.
#' @param ctmc_move_dat Data.frame of CTMC covariates. Required when \code{move_type == 1}.
#'   Must include columns \code{pop}, \code{regions}, \code{years}, \code{seas},
#'   \code{ages}, \code{sexes}, and any covariates referenced in \code{diffusion_formula}
#'   or \code{preference_formula}.
#' @param preference_formula R formula specifying preference (taxis) covariates
#'   for CTMC movement. Required when \code{move_type == 1}.
#' @param diffusion_formula R formula specifying diffusion covariates
#'   for CTMC movement. Required when \code{move_type == 1}.
#' @param log_move_diffusion_pars Log-scale diffusion parameters (\eqn{\theta_k}) for
#'   CTMC movement. Exponentiated and squared internally (\code{exp(2 * log_theta)}).
#'   Required when \code{move_type == 1}.
#' @param move_preference_pars Preference (taxis) parameters (\eqn{\gamma_k}) for CTMC
#'   movement, on the natural scale. Required when \code{move_type == 1}.
#' @param area_r Numeric vector of region areas (length \code{n_regions}) used to
#'   scale diffusion rates. Required when \code{move_type == 1}.
#' @param adjacency_mat Square \code{[n_regions x n_regions]} matrix defining
#'   connectivity among regions (1 = adjacent, 0 = not adjacent). Required when
#'   \code{move_type == 1}.
#' @param ctmc_diffusion_bounds How the off-diagonal generator entries are kept
#'   non-negative (a valid generator). Every form is evaluated on the adjacency
#'   edges only, so non-edges stay exactly zero. Taking \eqn{d} as the preference
#'   gradient \eqn{\gamma_i - \gamma_j} along the edge from \eqn{j} to \eqn{i} and
#'   \eqn{\theta_j} as the diffusion rate out of \eqn{j}:
#'   \describe{
#'     \item{\code{"none"} (or \code{0})}{\eqn{q = \theta_j + d}, unbounded. Valid
#'       only where diffusion outweighs taxis everywhere.}
#'     \item{\code{"softplus"} (or \code{1})}{softplus of \eqn{\theta_j + d} with
#'       width \code{ctmc_diffusion_eps}. Smooth, but an edge where taxis cancels
#'       diffusion carries a floor of \code{eps * log(2)}, so the width is a
#'       minimum exchange rate and not only a smoothing constant.}
#'     \item{\code{"upwind"} (or \code{2})}{discontinuous Galerkin / finite volume upwind flux,
#'       \eqn{q = \theta_j + \max(d, 0)}: diffusion is carried whole and only the
#'       down-gradient half of the taxis flux is added, so positivity never depends
#'       on the two cancelling. }
#'   }
#' @param ctmc_diffusion_eps Positive numeric width of the softplus used when
#'   \code{ctmc_diffusion_bounds} is \code{"softplus"}. Default 0.1.
#' @param seasdur Numeric vector of length \code{n_seas} giving season durations
#'   (summing to 1). Used to scale the CTMC generator when
#'   \code{ctmc_scale_by_seasdur == 1}. Defaults to \code{rep(1, n_seas)}, which
#'   reproduces the unscaled behavior.
#' @param ctmc_scale_by_seasdur Integer flag controlling the time units of the CTMC
#'   generator. \code{1} = treat \eqn{Q} as an annual rate and exponentiate
#'   \eqn{Q \cdot \mathrm{seasdur}[s]} for each season; \code{0} = treat \eqn{Q} as
#'   a per season rate and exponentiate it once per season irrespective of season
#'   duration. Only affects \code{move_type == 1} with \code{n_seas > 1}. Defaults
#'   to \code{0} here so that callers passing an unscaled generator get the
#'   arithmetic they expect; the user facing default is \code{1}, set by
#'   \code{Setup_Mod_Movement}.
#'
#' @param expm_nsub Integer controlling how the CTMC generator is exponentiated
#'   into movement fractions: \code{0} (default) uses \code{Matrix::expm}, a value
#'   \eqn{n \ge 1} uses the implicit backward Euler scheme \eqn{(I - Q\Delta/n)^{-n}},
#'   which is cheaper to differentiate but is a first-order approximation. Only read
#'   when \code{move_type == 1} and \code{use_fixed_movement == 0}. See
#'   \code{\link{mat_exp}}.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{Movement}}{Array of movement fractions dimensioned
#'     \code{[pop, from_region, to_region, year, seas, age, sex]}.
#'     Populated for all three movement configurations.}
#'   \item{\code{Mrate}}{Instantaneous rate matrix (generator \eqn{Q}) dimensioned
#'     \code{[pop, from_region, to_region, year, seas, age, sex]}. Populated only
#'     when \code{move_type == 1} and \code{use_fixed_movement == 0}; \code{NULL} otherwise.
#'     Stored \emph{unscaled} by season duration, so that consumers combining it with
#'     mortality must apply \code{seasdur[seas]} themselves (see
#'     \code{build_seas_operator}).}
#'   \item{\code{move_pen}}{Numeric movement penalty used for regularization.
#'     For CTMC movement, equal to \eqn{\sum_k \gamma_k^2}, a ridge on the
#'     preference coefficients applied once (not per stratum) that pins the
#'     otherwise unidentified level and spread. Zero for unstructured or fixed
#'     movement.}
#' }
#'
#' @details
#' Three movement configurations are supported:
#' \enumerate{
#'   \item Fixed movement (\code{use_fixed_movement == 1}): \code{Fixed_Movement} is used
#'     directly with no estimation.
#'   \item Unstructured multinomial logit movement (\code{move_type == 0}): movement
#'     fractions are estimated via a softmax transform of \code{move_pars + move_devs}.
#'   \item CTMC movement (\code{move_type == 1}): a generator matrix \eqn{Q = D + Z}
#'     is constructed from diffusion (\eqn{D}) and taxis (\eqn{Z}) components, then
#'     exponentiated via \code{\link{mat_exp}} (either \code{Matrix::expm} or the
#'     implicit solve, per \code{expm_nsub}) to obtain movement fractions. During
#'     projection years (\code{y > n_yrs}), covariate lookups are capped at \code{n_yrs}
#'     (i.e., CTMC base parameters are frozen at their last historical values), while
#'     \code{move_devs} continue to use the actual projected year index.
#' }
#'
#' @keywords internal
Get_Movement <- function(move_type,
                         do_recruits_move,
                         n_pop,
                         n_regions,
                         n_yrs,
                         n_proj_yrs_devs,
                         n_ages,
                         n_sexes,
                         n_seas,
                         move_pars,
                         move_devs,
                         use_fixed_movement,
                         Fixed_Movement = NULL,
                         ctmc_move_dat = NULL,
                         preference_formula = NULL,
                         diffusion_formula = NULL,
                         log_move_diffusion_pars,
                         move_preference_pars,
                         area_r,
                         adjacency_mat,
                         ctmc_diffusion_bounds,
                         ctmc_diffusion_eps = 0.1,
                         seasdur = rep(1, n_seas),
                         ctmc_scale_by_seasdur = 0,
                         expm_nsub = 0
) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  move_pen = 0 # initialize movement penalty if used
  Mrate = NULL # initialize for non-CTMC cases

  dims = list(pop = 1:n_pop,
              from = 1:n_regions,
              to = 1:n_regions,
              years = 1:(n_yrs + n_proj_yrs_devs),
              seas = 1:n_seas,
              ages = 1:n_ages,
              sexes = 1:n_sexes)

  # use fixed movement matrix
  if(use_fixed_movement == 1) {
    Movement = array(Fixed_Movement, dim = sapply(dims, length), dimnames = dims)
  } else if(move_type == 0) { # Unstructured markov movement

    Movement = array(0, dim = sapply(dims, length), dimnames = dims)
    ref_region = 1 # Set up reference region (always set at 0)
    n_move_yrs = n_yrs + n_proj_yrs_devs

    for(p in 1:n_pop) {
      for(r in 1:n_regions) {
        for(y in 1:n_move_yrs) {

          ypar = min(y, n_yrs) # projection years reuse the terminal year's parameters

          # movement out of this origin region in this year (container)
          move_pry = array(0, dim = c(n_regions, n_seas, n_ages, n_sexes))

          for(seas in 1:n_seas) {
            for(a in 1:n_ages) {
              for(s in 1:n_sexes) {

                # the reference region (1) sits at 0 and the remaining regions follow it
                move_tmp = c(0, move_pars[p,r,,ypar,seas,a,s] + move_devs[p,r,,y,seas,a,s])
                move_pry[,seas,a,s] = exp(move_tmp) / sum(exp(move_tmp)) # multinomial logit transform estimated movement

              } # end s loop
            } # end a loop
          } # end seas loop

          Movement[p,r,,y,,,] = move_pry

        } # end y loop
      } # end r loop
    } # end p loop

  } else if(move_type == 1) { # continuous markov chain movement with projection support

    # set up dimensions of movement matrix
    Mrate = Movement = array(0, dim = sapply(dims, length),  dimnames = dims)
    loop = expand.grid(dims[-(2:3)]) # get pop, year, age, and sexes to loop through
    if(do_recruits_move == 0) loop = loop[-which(loop$ages == 1),] # remove age 1, if recruits don't move

    # ctmc_move_dat holds one row per pop, region, year, season, age and sex, so its
    # rows can be found by position. The year axis is sized to whatever the covariates
    # carry, since ctmc_move_dat is allowed to hold projection year rows
    ctmc_key = sapply(c('pop','regions','years','seas','ages','sexes'), function(v) as.integer(ctmc_move_dat[,v])) # convert ctmc dataframe to matrix
    ctmc_row = array(NA_integer_, dim = pmax(c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes), apply(ctmc_key, 2, max))) # pmax to get projection year if there are any
    ctmc_row[ctmc_key] = seq_len(nrow(ctmc_move_dat))

    # setup design matrix
    design_mat = get_movement_dp_design_matrix(ctmc_move_dat, preference_formula, diffusion_formula)
    X_zk = design_mat$X_zk # diffusion
    W_zk = design_mat$W_zk # preference

    # diffusion rate from each region
    theta_k = exp(2 * log_move_diffusion_pars) # get diffusion parameter
    theta_z = (X_zk %*% theta_k)[,1] # multiply diffusion parameter by design matrix
    theta_z = theta_z/area_r[ctmc_move_dat[,'regions']]  # scale diffusion matrix by area

    # preference for each region. A formula with no terms (e.g. ~ 0) yields a zero-column
    # design matrix, which is the natural way to request pure diffusion with no taxis;
    # treat it as zero preference everywhere rather than indexing an empty object
    gamma_k = move_preference_pars # get preference parameters
    if(design_mat$n_gamma == 0) gamma_z = rep(0, nrow(ctmc_move_dat))
    else gamma_z = (W_zk %*% gamma_k)[,1] # multiply preference parameters by design matrix

    # ridge on the preference coefficients, applied once; pins both level and spread
    move_pen = move_pen + sum(gamma_k^2)

    # Generator edges, held as the (destination, origin) pairs the adjacency allows.
    # Every flow transform below runs on this value slot rather than on the full
    # n_regions x n_regions matrix, so a non-edge stays exactly zero by construction
    # rather than by cancellation, and the elementwise work scales with the number of
    # edges instead of with n_regions^2.
    bound_form = get_ctmc_bound_form(ctmc_diffusion_bounds)
    edge_ij = which(adjacency_mat == 1 & diag(1, n_regions) == 0, arr.ind = TRUE)
    edge_to = edge_ij[,1] # destination region
    edge_from = edge_ij[,2] # origin region
    edge_lin = edge_to + (edge_from - 1) * n_regions # linear index into [dest, origin]

    # move_devs counts the non-diagonal destinations of each origin in region order
    edge_dev = integer(length(edge_from))
    for(rr in 1:n_regions) {
      rr_edges = which(edge_from == rr)
      edge_dev[rr_edges[order(edge_to[rr_edges])]] = seq_along(rr_edges)
    } # end rr loop

    # Make instantaneous diffusion rate matrix
    for( index in seq_len(nrow(loop)) ){

      # stratum indices. Deviations always use the actual year
      pop_idx = loop$pop[index]
      y_idx = loop$years[index]
      seas_idx = loop$seas[index]
      a_idx = loop$ages[index]
      s_idx = loop$sexes[index]

      # Cap year spline look up parameters at n_yrs
      y_lookup = min(y_idx, n_yrs)

      # rows of ctmc_move_dat holding this stratum, one per region in region order
      which_index = ctmc_row[cbind(pop_idx, 1:n_regions, y_lookup, seas_idx, a_idx, s_idx)]

      # preference and diffusion for each strata, year, age, sex combination,
      # gathered onto the edges: d_e is the preference gradient along the edge and
      # t_e the diffusion rate out of its origin region
      pref_s = gamma_z[which_index] # get corresponding gammas
      theta_base = theta_z[which_index]

      d_e = pref_s[edge_to] - pref_s[edge_from]
      t_e = theta_base[edge_from]

      # Deviations scale the diffusion rate of their edge, before the flow transform.
      # move_devs is indexed as [origin_region, counter, year, seas, age, sex] where
      # counter goes through non-diagonal destinations for that origin
      for(e in seq_along(edge_lin)) {
        t_e[e] = t_e[e] * exp(move_devs[pop_idx, edge_from[e], edge_dev[e], y_idx, seas_idx, a_idx, s_idx])
      } # end e loop

      # edge flows, kept non-negative by the requested form
      if(bound_form == "none") q_e = t_e + d_e
      if(bound_form == "softplus") {
        u_e = t_e + d_e
        eps = ctmc_diffusion_eps # softplus width; softplus(0) = eps * log(2)
        q_e = (u_e + abs(u_e)) / 2 + eps * log(1 + exp(-abs(u_e) / eps))
      }
      # discontinuous Galerkin (upwind) flux: diffusion is carried whole and only the
      # down-gradient half of the taxis flux is added, so the two never cancel and a
      # down-gradient edge carries theta exactly. (d + abs(d)) / 2 is the positive part
      # written branch-free, since a branch would freeze at its tape-build value
      if(bound_form == "upwind") q_e = t_e + (d_e + abs(d_e)) / 2

      # scatter the edge flows back and conserve abundance
      Q_ss = adjacency_mat * 0
      Q_ss[edge_lin] = q_e
      diag(Q_ss) = -1 * Matrix::colSums(Q_ss) # diag to enforce sum to 0
      Q_ss = as(Q_ss, "sparseMatrix") # force sparse

      # Time units for turning the generator into fractions.
      dur = if(ctmc_scale_by_seasdur == 1) seasdur[seas_idx] else 1
      M_ss = mat_exp(Q_ss * dur, expm_nsub) # turn rate matrix into fractions

      # populate matrices
      Movement[pop_idx,,,y_idx,seas_idx,a_idx,s_idx] = t(M_ss)
      Mrate[pop_idx,,,y_idx,seas_idx,a_idx,s_idx] = t(as.matrix(Q_ss))

      # return penalty (Lagrange multiplier)

    } # end index loop
  }

  return(list(Movement = Movement, Mrate = Mrate, move_pen = move_pen))
}
