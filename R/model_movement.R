# Stage 2 of 3: objective function
#
# Turns movement parameters into the movement arrays the dynamics use, for both the unstructured and
# the continuous time parameterizations. Two helpers are shared with setup_movement.R.

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
#' Movement fractions under unstructured multinomial logit movement, a CTMC
#' generator, or a fixed matrix, plus the movement penalty. Under CTMC movement
#' the covariate lookups are capped at \code{n_yrs}, so the parameters are frozen
#' at their last historical values through the projection unless
#' \code{ctmc_move_dat} holds projection-year rows.
#'
#' @param move_type Integer: 0 unstructured Markov, 1 CTMC.
#' @param do_recruits_move Integer: 0 keeps recruits (age 1) in place, 1 moves them.
#' @param n_regions,n_ages,n_sexes,n_seas,n_pop Model dimensions.
#' @param n_yrs Number of observed (historical) years.
#' @param n_proj_yrs_devs Number of projected years, extending the year dim of the
#'   movement array beyond \code{n_yrs}. CTMC covariate lookups stay capped at
#'   \code{n_yrs} unless \code{ctmc_move_dat} holds projection-year rows.
#' @param move_pars Unstructured movement parameters \code{[pop, from_region,
#'   counter, year, seas, age, sex]}, \code{counter} indexing the
#'   \code{n_regions - 1} non-reference destinations. Ignored when
#'   \code{move_type == 1} or \code{use_fixed_movement == 1}.
#' @param move_devs Movement deviations \code{[pop, region, counter, year, seas,
#'   age, sex]}, always indexed on the actual (possibly projected) year. Under
#'   unstructured movement \code{counter} indexes the non-reference destinations
#'   and the deviations are additive on the logit scale. Under CTMC movement they
#'   sit on each region's preference instead, so \code{counter} has length one and
#'   \code{region} is the region whose preference is shifted.
#' @param use_fixed_movement Integer: 0 estimates movement, 1 uses
#'   \code{Fixed_Movement}.
#' @param Fixed_Movement Fixed movement matrix \code{[pop, from_region, to_region,
#'   year, seas, age, sex]}, read when \code{use_fixed_movement == 1}.
#' @param ctmc_move_dat Data frame of CTMC covariates, required when
#'   \code{move_type == 1}, with columns \code{pop}, \code{regions},
#'   \code{years}, \code{seas}, \code{ages}, \code{sexes} and any covariate the
#'   formulas name.
#' @param preference_formula,diffusion_formula Formulas for the preference (taxis)
#'   and diffusion covariates. Required when \code{move_type == 1}.
#' @param log_move_diffusion_pars Log-scale diffusion parameters
#'   (\eqn{\theta_k}), exponentiated and squared internally as
#'   \code{exp(2 * log_theta)}. Required when \code{move_type == 1}.
#' @param move_preference_pars Preference (taxis) parameters (\eqn{\gamma_k}) on
#'   the natural scale. Required when \code{move_type == 1}.
#' @param area_r Numeric vector \code{[n_regions]} of region areas, scaling the
#'   diffusion rates. Required when \code{move_type == 1}.
#' @param adjacency_mat Square \code{[n_regions x n_regions]} connectivity matrix,
#'   1 for adjacent and 0 for not. Required when \code{move_type == 1}.
#' @param ctmc_diffusion_bounds How the off-diagonal generator entries are kept
#'   non-negative. Every form is evaluated on the adjacency edges only, so
#'   non-edges stay exactly zero. With \eqn{d} the preference gradient
#'   \eqn{\gamma_i - \gamma_j} along the edge from \eqn{j} to \eqn{i} and
#'   \eqn{\theta_j} the diffusion rate out of \eqn{j}: \code{"none"} (or \code{0})
#'   is \eqn{q = \theta_j + d}, valid only where diffusion outweighs taxis
#'   everywhere; \code{"softplus"} (or \code{1}) is a softplus of that of width
#'   \code{ctmc_diffusion_eps}, smooth, but an edge where taxis cancels diffusion
#'   has a floor of \code{eps * log(2)}, so the width is a minimum exchange rate
#'   and not only a smoothing constant; \code{"upwind"} (or \code{2}) is the
#'   finite volume upwind flux \eqn{q = \theta_j + \max(d, 0)}, which keeps
#'   diffusion whole and adds only the down-gradient half of the taxis flux, so
#'   positivity never depends on the two cancelling.
#' @param ctmc_diffusion_eps Positive width of the softplus under
#'   \code{ctmc_diffusion_bounds = "softplus"}. Default 0.1.
#' @param seasdur Numeric vector \code{[n_seas]} of season durations summing to 1,
#'   used to scale the generator when \code{ctmc_scale_by_seasdur == 1}. Defaults
#'   to \code{rep(1, n_seas)}, the unscaled behavior.
#' @param ctmc_scale_by_seasdur Integer flag for the generator's time units.
#'   \code{1} treats \eqn{Q} as an annual rate and exponentiates
#'   \eqn{Q \cdot \mathrm{seasdur}[s]} each season; \code{0} treats it as a per
#'   season rate and exponentiates it once per season. Only matters under
#'   \code{move_type == 1} with \code{n_seas > 1}. Defaults to \code{0} here so a
#'   caller passing an unscaled generator gets the arithmetic it expects; the user
#'   facing default is \code{1}, set by \code{Setup_Mod_Movement}.
#' @param expm_nsub How the generator is exponentiated into movement fractions.
#'   \code{0} (default) uses \code{Matrix::expm}; \eqn{n \ge 1} uses the implicit
#'   backward Euler scheme \eqn{(I - Q\Delta/n)^{-n}}, cheaper to differentiate but
#'   first order. Read when \code{move_type == 1} and
#'   \code{use_fixed_movement == 0}. See \code{\link{mat_exp}}.
#'
#' @return A list with \code{Movement}, the movement fractions \code{[pop,
#'   from_region, to_region, year, seas, age, sex]}, populated under every
#'   configuration; \code{Mrate}, the generator \eqn{Q} on the same dims, populated
#'   only when \code{move_type == 1} and \code{use_fixed_movement == 0} and stored
#'   unscaled by season duration, so a consumer combining it with mortality applies
#'   \code{seasdur[seas]} itself (see \code{build_seas_operator}); and
#'   \code{move_pen}, the penalty, which under CTMC movement is
#'   \eqn{\sum_k \gamma_k^2}, a ridge on the preference coefficients applied once
#'   rather than per stratum that pins the otherwise unidentified level and spread,
#'   and is zero for unstructured or fixed movement.
#'
#' @details
#' Fixed movement (\code{use_fixed_movement == 1}) uses \code{Fixed_Movement}
#' directly. Unstructured movement (\code{move_type == 0}) takes a softmax of
#' \code{move_pars + move_devs}. CTMC movement (\code{move_type == 1}) builds a
#' generator \eqn{Q = D + Z} from its diffusion and taxis components and
#' exponentiates it through \code{\link{mat_exp}}.
#'
#' Under CTMC movement a deviation is added to its own region's preference,
#' \eqn{h_r + \epsilon_r}, before the gradient \eqn{d = h_i - h_j} is taken along
#' each edge. One deviation therefore moves every edge that touches that region,
#' raising the rates into it and lowering the rates out of it, and the perturbed
#' taxis is still the gradient of a surface. Only differences reach the generator,
#' so a constant added to every region's deviation leaves movement unchanged and
#' the level of the field is held only by the process error penalty.
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

    # deviations sit on a region's preference here, one per region, so a dev array built for the
    # unstructured model would silently hand every region its first destination's deviation
    if(dim(move_devs)[3] != 1) stop("move_devs has ", dim(move_devs)[3], " destinations, and CTMC movement holds its deviations on each region's preference, so that axis must be 1. Rebuild the inputs with Setup_Mod_Movement.")

    # set up dimensions of movement matrix
    Mrate = Movement = array(0, dim = sapply(dims, length),  dimnames = dims)
    loop = expand.grid(dims[-(2:3)]) # get pop, year, age, and sexes to loop through
    if(do_recruits_move == 0) loop = loop[-which(loop$ages == 1),] # remove age 1, if recruits don't move

    # ctmc_move_dat holds one row per pop, region, year, season, age and sex, so rows are found
    # by position. the year axis is sized to the covariates, which may include projection years
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
    theta_z = theta_z / area_r[ctmc_move_dat[,'regions']]  # scale diffusion matrix by area

    # preference for each region. a formula with no terms (~ 0) gives a zero-column design matrix,
    # which is how pure diffusion is requested, so treat it as zero preference everywhere
    gamma_k = move_preference_pars # get preference parameters
    if(design_mat$n_gamma == 0) gamma_z = rep(0, nrow(ctmc_move_dat))
    else gamma_z = (W_zk %*% gamma_k)[,1] # multiply preference parameters by design matrix

    # ridge on the preference coefficients, applied once; pins both level and spread
    move_pen = move_pen + sum(gamma_k^2)

    # generator edges, as the (destination, origin) pairs the adjacency allows. every flow transform
    # runs on this value slot, so a non-edge stays exactly zero and work scales with edge count
    bound_form = get_ctmc_bound_form(ctmc_diffusion_bounds)
    edge_ij = which(adjacency_mat == 1 & diag(1, n_regions) == 0, arr.ind = TRUE)
    edge_to = edge_ij[,1] # destination region
    edge_from = edge_ij[,2] # origin region
    edge_lin = edge_to + (edge_from - 1) * n_regions # linear index into [dest, origin]

    # Make instantaneous diffusion rate matrix
    for(index in seq_len(nrow(loop))){

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

      # deviations shift a region's own preference, so one moves every edge that touches that region
      pref_s = gamma_z[which_index] + move_devs[pop_idx,,1,y_idx,seas_idx,a_idx,s_idx] # get corresponding gammas
      theta_base = theta_z[which_index] # get corresponding thetas

      d_e = pref_s[edge_to] - pref_s[edge_from]
      t_e = theta_base[edge_from]

      # edge flows, kept non-negative
      if(bound_form == "none") q_e = t_e + d_e
      if(bound_form == "softplus") {
        u_e = t_e + d_e
        eps = ctmc_diffusion_eps # softplus width; softplus(0) = eps * log(2)
        q_e = (u_e + abs(u_e)) / 2 + eps * log(1 + exp(-abs(u_e) / eps))
      }

      # discontinuous Galerkin (upwind) flux
      if(bound_form == "upwind") q_e = t_e + (d_e + abs(d_e)) / 2

      # distribute the edge flows back and conserve abundance
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

    } # end index loop
  }

  return(list(Movement = Movement, Mrate = Mrate, move_pen = move_pen))
}
