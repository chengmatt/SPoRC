#' Get Design Matrices for CTMC Movement
#'
#' Constructs the design matrices for the diffusion and preference components of a
#' Continuous Time Markov Chain (CTMC) movement model. These matrices are used
#' to parameterize movement rates in terms of covariates specified by formulas.
#'
#' @param data A \code{data.frame} containing the CTMC covariates. Must include all variables
#'   referenced in \code{diffusion_formula} and \code{preference_formula}.
#' @param preference_formula An R formula describing the linear predictor for movement
#'   preference. Variables must exist in \code{data}.
#' @param diffusion_formula An R formula describing the linear predictor for diffusion
#'   rates. Variables must exist in \code{data}.
#'
#' @return A \code{list} with the following components:
#' \describe{
#'   \item{\code{n_theta}}{Number of diffusion parameters (number of columns in diffusion design matrix).}
#'   \item{\code{n_gamma}}{Number of preference parameters (number of columns in preference design matrix).}
#'   \item{\code{X_zk}}{Diffusion design matrix constructed from \code{diffusion_formula} and \code{data}.}
#'   \item{\code{W_zk}}{Preference design matrix constructed from \code{preference_formula} and \code{data}.}
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
#' years are supported: base parameters (preference/diffusion) can either be frozen
#' at the last historical year or extended via user-provided covariates in ctmc_move_dat.
#'
#' @param move_type Integer flag indicating movement type: 0 = unstructured Markov,
#'   1 = CTMC movement.
#' @param do_recruits_move Integer flag: 0 = recruits do not move, 1 = recruits move.
#' @param n_regions Number of spatial regions.
#' @param n_yrs Number of years in the observed data.
#' @param n_proj_yrs_devs Number of projected years for deviations.
#' @param n_ages Number of age classes.
#' @param n_sexes Number of sexes.
#' @param move_pars Array of movement parameters for unstructured movement.
#' @param move_devs Array of movement deviations (applies to both unstructured and CTMC movement).
#' @param use_fixed_movement Integer flag: 0 = estimate movement, 1 = use fixed matrix.
#' @param Fixed_Movement Optional fixed movement matrix.
#' @param ctmc_move_dat Data.frame with CTMC covariates used to build design matrices
#'   for diffusion and preference. Required columns (when \code{move_type == 1}) include
#'   \code{regions}, \code{years}, \code{ages}, and \code{sexes}, plus any covariates
#'   referenced in \code{diffusion_formula} and \code{preference_formula}.
#'   Can include projection years (years > n_yrs) with projected covariate values.
#'   Year effects in formulas (e.g., splines) are automatically capped at \code{n_yrs}
#'   to prevent extrapolation, while other covariates use their actual projected values.
#' @param preference_formula R formula specifying preference covariates for CTMC movement.
#' @param diffusion_formula R formula specifying diffusion covariates for CTMC movement.
#' @param log_move_diffusion_pars Log-transformed diffusion parameters for CTMC movement.
#' @param move_preference_pars Preference parameters for CTMC movement.
#' @param area_r Vector of areas for each region (used for scaling diffusion rates).
#' @param adjacency_mat Square adjacency matrix defining connectivity between regions for CTMC movement.
#' @param ctmc_diffusion_bounds Integer flag: 1 = apply diffusion bounds to generator matrix, 0 = no bounds.
#' @param n_seas Number of seasons
#' @param n_pop Number of populations
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{Movement}}{Array of movement fractions for each stratum (pop x from regions × to regions × years × seas, ages × sexes).}
#'   \item{\code{Mrate}}{Instantaneous movement rate matrix if CTMC movement is used (pop x from regions × to regions × years × seas x ages × sexes); otherwise NULL.}
#'   \item{\code{move_pen}}{Numeric value of movement penalty calculated from preference parameters (for CTMC only).}
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
                         ctmc_diffusion_bounds
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

    # preference for each region
    gamma_k = move_preference_pars # get preference parameters
    gamma_z = (W_zk %*% gamma_k)[,1] # multiply preference parameters by design matrix

    # Make instantaneous diffusion rate matrix
    for( index in seq_len(nrow(loop)) ){

      # get pop, year, age, and sex specific indices for a given stratum combination
      which_rows = expand.grid(1:n_pop, 1:n_regions, loop[index,"years"], loop[index,'seas'], loop[index,"ages"], loop[index,"sexes"] )
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
        for(j in 1:n_regions) {
          minval = min(Z_ss[,j])  # minimum value in column j
          for(i in 1:n_regions) D_ss[i,j] = D_ss[i,j] - adjacency_mat[i,j] * minval
        } # end j loop
      }

      # conserve abundance
      diag(D_ss) = -1 * Matrix::colSums(D_ss) # diag to enforce sum to 0
      diag(Z_ss) = -1 * Matrix::colSums(Z_ss) # diag to enforce sum to 0
      D_ss = as(D_ss, "sparseMatrix") # force sparse

      Q_ss = D_ss + Z_ss # rate matrix
      M_ss = Matrix::expm( Q_ss ) # turn rate matrix into fractions

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
