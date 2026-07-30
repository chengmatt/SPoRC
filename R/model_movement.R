# Stage 2 of 3: objective function
#
# Turns movement parameters into the movement arrays the dynamics use, for both
# the unstructured and the continuous time parameterizations.
# get_movement_dp_design_matrix is shared with setup_movement.R.

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
  X_zk = model.matrix(diffusion_formula, data) # diffusion design matrix
  W_zk = model.matrix(preference_formula, data) # preference design matrix
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
#' @param ctmc_diffusion_bounds Integer flag: 1 = shift diffusion columns to ensure
#'   all off-diagonal generator matrix entries are non-negative (valid generator);
#'   0 = no bounds applied.
#' @param seasdur Numeric vector of length \code{n_seas} giving season durations
#'   (summing to 1). Used to scale the CTMC generator when
#'   \code{ctmc_scale_by_seasdur == 1}. Defaults to \code{rep(1, n_seas)}, which
#'   reproduces the unscaled behaviour.
#' @param ctmc_scale_by_seasdur Integer flag controlling the time units of the CTMC
#'   generator. \code{1} = treat \eqn{Q} as an annual rate and exponentiate
#'   \eqn{Q \cdot \mathrm{seasdur}[s]} for each season; \code{0} = treat \eqn{Q} as
#'   a per season rate and exponentiate it once per season irrespective of season
#'   duration. Only affects \code{move_type == 1} with \code{n_seas > 1}. Defaults
#'   to \code{0} here so that callers passing an unscaled generator get the
#'   arithmetic they expect; the user facing default is \code{1}, set by
#'   \code{Setup_Mod_Movement}.
#'
#'   Whether this flag changes the fit or only reparameterises it depends on
#'   whether the generator varies by season. With a season-agnostic generator it
#'   matters: the seasonal steps commute, so scaling on composes across the year to
#'   \eqn{\exp(Q)} regardless of \code{n_seas}, whereas scaling off composes to
#'   \eqn{\exp(n_{seas} Q)} and inflates movement as the seasonal time step shrinks.
#'   With a season-varying generator the scaling is absorbed --- \eqn{Q} is linear in
#'   \eqn{\theta}, so scaling \eqn{Q} by \code{seasdur[s]} equals shifting
#'   \code{log_move_diffusion_pars} by \eqn{\tfrac{1}{2}\log \mathrm{seasdur}[s]} ---
#'   but only if \emph{both} formulas carry a season term, since the flag scales
#'   diffusion and taxis together while only \eqn{\theta} can absorb it.
#'
#'   Note that even where the scaling is absorbed, it changes what the estimated
#'   diffusion means (per-season step vs annual rate, hence not comparable across
#'   seasons of unequal length), and the Dirichlet movement prior is evaluated on
#'   annual fractions \eqn{\exp(Q)} irrespective of this flag. Under
#'   \code{move_timing = 2} the flag governs only the reported \code{Movement}
#'   diagnostic: the dynamics take \code{Mrate} and apply \code{seasdur[seas]}
#'   themselves.
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
#'     For CTMC movement, equal to \eqn{\sum_{\text{strata}} (\sum_r \gamma_{z,r})^2},
#'     which penalizes large net preference values across regions. Zero for unstructured
#'     or fixed movement.}
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
#'     exponentiated via \code{Matrix::expm(Q)} to obtain movement fractions. During
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
                         seasdur = rep(1, n_seas),
                         ctmc_scale_by_seasdur = 0
) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  move_pen = 0 # initialize movement penalty if used
  Mrate = NULL # initialize for non-CTMC cases
  dims = list(pop = 1:n_pop, from = 1:n_regions, to = 1:n_regions, years = 1:(n_yrs + n_proj_yrs_devs), seas = 1:n_seas, ages = 1:n_ages, sexes = 1:n_sexes)

  # use fixed movement matrix
  if(use_fixed_movement == 1) {
    Movement = array(Fixed_Movement, dim = sapply(dims, length), dimnames = dims)
  } else if(move_type == 0) { # Unstructured markov movement

    Movement = array(0, dim = sapply(dims, length), dimnames = dims)
    ref_region = 1 # Set up reference region (always set at 0)

    for(p in 1:n_pop) {
      for(r in 1:n_regions) {
        for(y in 1:(n_yrs + n_proj_yrs_devs)) {
          for(seas in 1:n_seas) {
            for(a in 1:n_ages) {
              for(s in 1:n_sexes) {

                move_tmp = rep(0, n_regions) # temporary movement vector to store values
                counter = 1  # counter

                for(rr in 1:n_regions) {
                  if(rr != ref_region) {
                    # extract movement parameters
                    if(y <= n_yrs) tmp_move_pars = move_pars[p,r,counter,y,seas,a,s]
                    else tmp_move_pars = move_pars[p,r,counter,n_yrs,seas,a,s]
                    move_tmp[rr] = tmp_move_pars + move_devs[p,r,counter,y,seas,a,s]
                    counter = counter + 1
                  } # end if not reference region
                } # end rr loop

                Movement[p,r,,y,seas,a,s] = exp(move_tmp) / sum(exp(move_tmp)) # multinomial logit transform estimated movement

              } # end s loop
            } # end a loop
          } # end seas loop
        } # end y loop
      } # end r loop
    } # end p loop

  } else if(move_type == 1) { # continuous markov chain movement with projection support

    # set up dimensions of movement matrix
    Mrate = Movement = Taxis = Diffusion = array(0, dim = sapply(dims, length),  dimnames = dims)
    loop = expand.grid(dims[-(2:3)]) # get pop, year, age, and sexes to loop through
    if(do_recruits_move == 0) loop = loop[-which(loop$ages == 1),] # remove age 1, if recruits don't move

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

    # Make instantaneous diffusion rate matrix
    for( index in seq_len(nrow(loop)) ){

      # get pop, year, age, and sex specific indices for a given stratum combination
      which_rows = expand.grid(loop[index,"pop"], 1:n_regions, loop[index,"years"], loop[index,'seas'], loop[index,"ages"], loop[index,"sexes"] )
      which_rows$index = NA
      colnames(which_rows) = c("pop", "regions", "years", "seas", "ages", "sexes", "index" )

      # Cap year spline look up parameters at n_yrs
      y_lookup = min(loop[index,"years"], n_yrs)

      # Match the current stratum (region, year, seas, age, sex) to rows in ctmc_move_dat
      for( i2 in seq_len(nrow(which_rows)) ){
        which_rows$index[i2] = which((which_rows[i2,'pop'] == ctmc_move_dat[,'pop']) &
                                      (which_rows[i2,'regions'] == ctmc_move_dat[,'regions']) &
                                       y_lookup == ctmc_move_dat[,"years"] &
                                       (which_rows[i2,'seas'] == ctmc_move_dat[,'seas']) &
                                       (which_rows[i2,'ages'] == ctmc_move_dat[,'ages']) &
                                       (which_rows[i2,'sexes'] == ctmc_move_dat[,'sexes']) )
      }

      # preference for each strata, year, age, sex combination
      pref_s = gamma_z[which_rows$index] # get corresponding gammas
      Z_ss = adjacency_mat * outer( pref_s, pref_s, FUN = "-" )

      # base diffusion parameters for this stratum
      theta_base = theta_z[which_rows$index]

      # create base diffusion matrix (w/ corresponding thetas)
      D_ss = adjacency_mat %*% diag(theta_base, n_regions)

      # Add origin-destination deviations (always uses actual year, not y_lookup)
      pop_idx = loop$pop[index]
      y_idx = loop$years[index]
      seas_idx = loop$seas[index]
      a_idx = loop$ages[index]
      s_idx = loop$sexes[index]

      # Note: move_devs is indexed as [origin_region, counter, year, seas, age, sex]
      # where counter goes through non-diagonal destinations for that origin
      for(rr in 1:n_regions) {  # rr = origin (from)
        counter = 1  # Reset counter for each origin region
        for(r in 1:n_regions) {  # r = destination (to)
          # Only apply deviations to off-diagonal elements (actual transitions, not residency)
          if(adjacency_mat[r, rr] == 1 && r != rr) {  # if adjacent BUT NOT diagonal
            # Apply deviation: rr is origin, counter indexes non-diagonal destinations
            D_ss[r, rr] = D_ss[r, rr] * exp(move_devs[pop_idx, rr, counter, y_idx, seas_idx, a_idx, s_idx])
            counter = counter + 1  # Increment counter for next valid destination from rr
          } # end if
        } # end r (to)
      } # end rr (from)

      # apply diffusion bounds to ensure valid generator matrix
      if(ctmc_diffusion_bounds == 1) { # ensure D_ss(i,j) + Z_ss(i,j) > 0 for all i != j
        eps <- 0.1
        for(j in 1:n_regions) {
          zj = Z_ss[,j]
          mj = sum(zj) / n_regions
          minval = mj - eps * log(sum(exp((mj - zj) / eps)))  # smooth minimum of column j
          for(i in 1:n_regions) D_ss[i,j] = D_ss[i,j] - adjacency_mat[i,j] * minval
        } # end j loop
      }

      # conserve abundance
      diag(D_ss) = -1 * Matrix::colSums(D_ss) # diag to enforce sum to 0
      diag(Z_ss) = -1 * Matrix::colSums(Z_ss) # diag to enforce sum to 0
      D_ss = as(D_ss, "sparseMatrix") # force sparse

      Q_ss = D_ss + Z_ss # rate matrix

      # Time units for turning the generator into fractions.
      dur = if(ctmc_scale_by_seasdur == 1) seasdur[seas_idx] else 1
      M_ss = Matrix::expm( Q_ss * dur ) # turn rate matrix into fractions

      # populate matrices
      Movement[pop_idx,,,y_idx,seas_idx,a_idx,s_idx] = t(as.matrix(M_ss))
      Taxis[pop_idx,,,y_idx,seas_idx,a_idx,s_idx] = t(as.matrix(Z_ss))
      Diffusion[pop_idx,,,y_idx,seas_idx,a_idx,s_idx] = t(as.matrix(D_ss))
      Mrate[pop_idx,,,y_idx,seas_idx,a_idx,s_idx] = t(as.matrix(Q_ss))

      # return penalty (Lagrange multiplier)
      move_pen = move_pen + sum(pref_s)^2

    } # end index loop
  }

  return(list(Movement = Movement, Mrate = Mrate, move_pen = move_pen))
}
