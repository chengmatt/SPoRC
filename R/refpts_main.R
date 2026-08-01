# Stage 3 of 3: post fit
#
# Reference point entry point. Get_Reference_Points chooses a solver based on the
# spatial structure and the requested quantity, then optim_ref_pts runs it as its
# own small RTMB problem. build_plus_group_T and solve_plus_group give the
# solvers a shared treatment of the plus group under movement.

#' Optimize Reference Point Models
#'
#' Constructs and optimizes an RTMB automatic differentiation objective function
#' for estimating fisheries reference points (e.g., SPR-based or Fmsy-based
#' biological reference points). After optimization, retrieves the model report
#' and standard deviation report from the best parameter estimates.
#'
#' @param model_name Function. An RTMB-compatible model function (e.g.,
#'   \code{global_SPR}, \code{global_BH_Fmsy}, \code{local_BH_Fmsy}) that
#'   defines the objective function for the reference point calculation.
#' @param data_list List. A named list of data inputs passed to \code{model_name}
#'   via \code{\link[RTMB]{MakeADFun}}. Should contain all quantities treated as
#'   fixed data within the reference point model (e.g., biological parameters,
#'   selectivity, movement matrices).
#' @param pars_list List. A named list of initial parameter values passed to
#'   \code{\link[RTMB]{MakeADFun}}. These are the parameters over which the
#'   objective function is optimized (e.g., \code{ln_Fmsy}, \code{ln_F_spr}).
#'
#' @returns An RTMB AD function object (list) with the following additional
#'   elements appended after optimization:
#'   \describe{
#'     \item{\code{$optim}}{Output from \code{\link[stats]{nlminb}}, including
#'       convergence code, final objective value, and optimized parameter
#'       estimates.}
#'     \item{\code{$rep}}{Named list of reported quantities from the model
#'       (e.g., equilibrium SSB, yield, reference point values), evaluated at
#'       the best parameter estimates via \code{$env$last.par.best}.}
#'     \item{\code{$sd_rep}}{Output from \code{\link[RTMB]{sdreport}}, containing
#'       standard errors and summary statistics for all estimated and derived
#'       quantities.}
#'   }
#'
#'
#' @keywords internal
optim_ref_pts <- function(model_name, data_list, pars_list) {
  tmp_obj <- RTMB::MakeADFun(cmb(model_name, data_list), parameters = pars_list, random = NULL, silent = TRUE)
  tmp_obj$optim <- stats::nlminb(tmp_obj$par, tmp_obj$fn, tmp_obj$gr, control = list(iter.max = 1e6, eval.max = 1e6, rel.tol = 1e-15))
  tmp_obj$rep <- tmp_obj$report(tmp_obj$env$last.par.best) # get report
  tmp_obj$sd_rep <- sdreport(tmp_obj)
  return(tmp_obj)
}

#' Build annual transition matrices for the plus-group analytical solution
#'
#' Constructs the four annual transition matrices needed to solve for the
#' equilibrium plus-group abundance analytically. Each matrix accumulates
#' survival and movement across all seasons for either the penultimate age or
#' the plus-group age, under either unfished or fished conditions.
#'
#' The equilibrium plus-group vector \eqn{N_+} satisfies
#' \deqn{N_+ = T_+ N_+ + T_{n-1} N_{n-1}}
#' which rearranges to \eqn{(I - T_+) N_+ = T_{n-1} N_{n-1}}, solved in
#' \code{\link{solve_plus_group}}.
#'
#' All arguments are sliced by the caller to remove the population dimension,
#' so this helper works identically for the single-population spatial case
#' (\code{global_BH_Fmsy}, \code{local_BH_Fmsy_sglpop}) and the
#' multi-population case (\code{global_SPR}, \code{local_BH_Fmsy_multipop}).
#'
#' @param M_penult Numeric vector \code{[n_regions]}. Natural mortality for
#'   the penultimate age class, used as an annual rate (scaled by
#'   \code{seasdur} internally).
#' @param M_plus Numeric vector \code{[n_regions]}. Natural mortality for the
#'   plus-group age class.
#' @param F_penult Numeric matrix \code{[n_regions, n_seas]}. Total fishing
#'   mortality per season for the penultimate age, already summed across
#'   fleets.
#' @param F_plus Numeric matrix \code{[n_regions, n_seas]}. Total fishing
#'   mortality per season for the plus-group age.
#' @param Mov_penult Numeric array \code{[n_regions, n_regions, n_seas]}.
#'   Movement transition matrices for the penultimate age. Entry
#'   \code{[r1, r2, s]} is the probability of moving from region \code{r1}
#'   to region \code{r2} in season \code{s}.
#' @param Mov_plus Numeric array \code{[n_regions, n_regions, n_seas]}.
#'   Movement transition matrices for the plus-group age.
#' @param n_regions Integer. Number of spatial regions.
#' @param n_seas Integer. Number of seasons.
#' @param seasdur Numeric vector \code{[n_seas]}. Fractional duration of each
#'   season (must sum to one).
#' @param Mrate_penult Numeric array \code{[n_regions, n_regions, n_seas]}.
#'   Instantaneous movement rate matrices (generator \eqn{Q}) for the penultimate
#'   age, in the same \code{[from, to]} convention as \code{Mov_penult}. Required
#'   when \code{move_timing = 2}, ignored otherwise.
#' @param Mrate_plus Numeric array \code{[n_regions, n_regions, n_seas]}.
#'   Instantaneous movement rate matrices for the plus-group age. Required when
#'   \code{move_timing = 2}, ignored otherwise.
#' @param move_timing Integer flag for movement/mortality sequencing:
#'   \code{0} = movement then mortality (default), \code{1} = mortality then
#'   movement, \code{2} = continuous. See \code{\link{build_seas_operator}}.
#'
#' @return A named list with four transition matrices, each of dimension
#'   \code{[n_regions, n_regions]}:
#'   \describe{
#'     \item{\code{T_penult_unfished}}{Annual transition for the penultimate
#'       age under unfished conditions.}
#'     \item{\code{T_plus_unfished}}{Annual transition for the plus-group age
#'       under unfished conditions.}
#'     \item{\code{T_penult_fished}}{Annual transition for the penultimate age
#'       under fished conditions.}
#'     \item{\code{T_plus_fished}}{Annual transition for the plus-group age
#'       under fished conditions.}
#'   }
#'
#' @keywords internal
build_plus_group_T <- function(M_penult, M_plus, F_penult, F_plus,
                               Mov_penult, Mov_plus,
                               n_regions, n_seas, seasdur,
                               Mrate_penult = NULL, Mrate_plus = NULL,
                               move_timing = 0) {

  # initialize transition matrices
  T_pu <- T_lu <- T_pf <- T_lf <- diag(n_regions)

  # build_seas_operator returns row convention [from, to]; the plus-group recursion
  # composes column-convention operators (applied as T %*% N), hence the transpose.
  seas_op <- function(Move, Z, Q, dur) t(build_seas_operator(Move, Z, Q, dur, move_timing))

  for (seas in seq_len(n_seas)) {
    Zu_p <- M_penult * seasdur[seas]                      # unfished total mortality
    Zu_l <- M_plus   * seasdur[seas]
    Zf_p <- M_penult * seasdur[seas] + F_penult[,seas]    # fished total mortality
    Zf_l <- M_plus   * seasdur[seas] + F_plus[,seas]
    Mp   <- Mov_penult[,, seas] # movement
    Ml   <- Mov_plus[,,   seas] # movement
    Qp   <- if(is.null(Mrate_penult)) NULL else Mrate_penult[,, seas] # instantaneous rates
    Ql   <- if(is.null(Mrate_plus))   NULL else Mrate_plus[,, seas]
    # transition matrices
    T_pu <- seas_op(Mp, Zu_p, Qp, seasdur[seas]) %*% T_pu
    T_lu <- seas_op(Ml, Zu_l, Ql, seasdur[seas]) %*% T_lu
    T_pf <- seas_op(Mp, Zf_p, Qp, seasdur[seas]) %*% T_pf
    T_lf <- seas_op(Ml, Zf_l, Ql, seasdur[seas]) %*% T_lf
  } # end seas loop


  list(T_penult_unfished = T_pu,
       T_plus_unfished = T_lu,
       T_penult_fished = T_pf,
       T_plus_fished = T_lf)
}


#' Solve for equilibrium plus-group numbers given transition matrices
#'
#' Given the annual transition matrices produced by
#' \code{\link{build_plus_group_T}} and the penultimate-age abundance vector,
#' solves the linear system \eqn{(I - T_+) N_+ = T_{n-1} N_{n-1}} for the
#' equilibrium plus-group abundance under both unfished and fished conditions.
#'
#' @param Ts Named list returned by \code{\link{build_plus_group_T}}.
#' @param N_penult_u Numeric vector \code{[n_regions]}. Unfished
#'   penultimate-age abundance (per-recruit) at the start of the year.
#' @param N_penult_f Numeric vector \code{[n_regions]}. Fished
#'   penultimate-age abundance (per-recruit) at the start of the year.
#' @param n_regions Integer. Number of spatial regions.
#'
#' @return A named list:
#'   \describe{
#'     \item{\code{unfished}}{Numeric vector \code{[n_regions]}. Equilibrium
#'       plus-group abundance per recruit under unfished conditions.}
#'     \item{\code{fished}}{Numeric vector \code{[n_regions]}. Equilibrium
#'       plus-group abundance per recruit under fished conditions.}
#'   }
#'
#' @keywords internal
solve_plus_group <- function(Ts, N_penult_u, N_penult_f, n_regions) {
  I <- diag(n_regions)
  list(
    unfished = solve(I - Ts$T_plus_unfished, Ts$T_penult_unfished %*% N_penult_u),
    fished   = solve(I - Ts$T_plus_fished,   Ts$T_penult_fished   %*% N_penult_f)
  )
}


#' Compute fishing and biological reference points from an assessment or simulation
#'
#' Wrapper that constructs the appropriate data list, calls the relevant
#' inner objective function via RTMB, and returns fishing and biological
#' reference points for use in projections or harvest control rules. Supports
#' single-region and spatially explicit multi-region models, with options for
#' SPR-based or Beverton-Holt MSY-based reference points.
#'
#' @param data List. SPoRC data object containing age structure, weight-at-age,
#'   maturity, natural mortality, seasons, and spatial configuration.
#' @param rep List. SPoRC report object from RTMB containing estimated or
#'   simulated quantities including \code{Fmort}, \code{fish_sel},
#'   \code{natmort}, \code{Rec}, \code{SSB}, \code{h_trans}, \code{R0},
#'   \code{rec_region_prop}, \code{rec_seas_prop}, \code{Movement}, and
#'   \code{stray_rate}.
#' @param SPR_x Numeric. Target spawning potential ratio fraction (e.g. 0.4).
#'   Required when \code{what} is \code{"SPR"}, \code{"independent_SPR"}, or
#'   \code{"global_SPR"}.
#' @param t_spawn Numeric. Fraction of the spawning season elapsed before
#'   spawning, used for the mid-season mortality correction. Default = 0.
#' @param sex_ratio_f Numeric array \code{[n_pop, n_regions]}. Female sex
#'   ratio at recruitment. Default = 0.5 everywhere.
#' @param calc_rec_st_yr Integer. First year included when computing mean
#'   historical recruitment for biological reference point scaling. Default = 1.
#' @param rec_age Integer. Recruitment lag in years, used to exclude the most
#'   recent years from the mean recruitment calculation. Default = 1.
#' @param type Character. Spatial structure of the model:
#'   \describe{
#'     \item{\code{"single_region"}}{No spatial movement; supports
#'       \code{"SPR"} and \code{"BH_MSY"}.}
#'     \item{\code{"multi_region"}}{Spatially explicit; supports
#'       \code{"independent_SPR"}, \code{"independent_BH_MSY"},
#'       \code{"global_SPR"}, \code{"global_BH_MSY"}, and
#'       \code{"local_BH_MSY"}.}
#'   }
#' @param what Character. Reference point method:
#'   \describe{
#'     \item{\code{"SPR"}}{Single-region \eqn{F_{SPR_x}}.}
#'     \item{\code{"BH_MSY"}}{Single-region Beverton-Holt \eqn{F_{MSY}}.}
#'     \item{\code{"independent_SPR"}}{Per-region \eqn{F_{SPR_x}} computed
#'       independently for each region without movement.}
#'     \item{\code{"independent_BH_MSY"}}{Per-region \eqn{F_{MSY}} computed
#'       independently for each region without movement.}
#'     \item{\code{"global_SPR"}}{Single shared \eqn{F_{SPR_x}} with
#'       movement, integrated across all regions.}
#'     \item{\code{"global_BH_MSY"}}{Single shared \eqn{F_{MSY}} with
#'       movement. Valid for single-population models only.}
#'     \item{\code{"local_BH_MSY"}}{Region-specific \eqn{F_{MSY}} values
#'       that jointly maximise total yield with movement. Valid for both
#'       single- and multi-population models.}
#'   }
#' @param n_avg_yrs Integer. Number of terminal years over which demographic
#'   rates (selectivity, natural mortality, weight, maturity, movement) are
#'   averaged before computing reference points. Default = 1.
#' @param local_bh_msy_newton_steps Integer. Number of Newton-Raphson
#'   iterations used to solve for equilibrium recruitment by origin region
#'   when \code{what = "local_BH_MSY"}. Increase if convergence is suspect.
#'   Default = 6.
#' @param is_discard_fleet Integer vector \code{[n_fish_fleets]}. Indicator
#'   for fleets whose catch should be excluded from landed yield when
#'   computing MSY-based reference points (0 = landing fleet, 1 = discard-only
#'   fleet). These fleets still contribute to total fishing mortality \code{Z}
#'   and affect population dynamics and spawning biomass. Only used by
#'   Beverton-Holt MSY methods (\code{"BH_MSY"}, \code{"independent_BH_MSY"},
#'   \code{"global_BH_MSY"}, \code{"local_BH_MSY"}); ignored for SPR-based
#'   methods. Default is all zeros (all fleets are landing fleets).
#'
#' @return A named list:
#'   \describe{
#'     \item{\code{f_ref_pt}}{Numeric vector \code{[n_regions]}. Fishing
#'       mortality reference point by region. All regions share the same value
#'       for global methods; regions have independent values for local or
#'       independent methods.}
#'     \item{\code{b_ref_pt}}{Numeric array \code{[n_pop, n_regions]}.
#'       Equilibrium spawning biomass at the reference point by population and
#'       region (\eqn{SBPR_F \times R_{eq}} or \eqn{SBPR_F \times \bar{R}}).}
#'     \item{\code{virgin_b_ref_pt}}{Numeric array \code{[n_pop, n_regions]}.
#'       Virgin (unfished) spawning biomass by population and region
#'       (\eqn{SBPR_0 \times R_0} or \eqn{SBPR_0 \times \bar{R}}).}
#'     \item{\code{pop_b_ref_pt}}{Numeric array \code{[n_pop, n_regions]}.
#'       Population-specific effective spawning biomass at the reference point,
#'       evaluated at each population's natal region and incorporating stray
#'       contributions from other populations.}
#'     \item{\code{virgin_pop_b_ref_pt}}{Numeric array \code{[n_pop, n_regions]}.
#'       Population-specific effective virgin spawning biomass, evaluated at
#'       each population's natal region.}
#'   }
#'
#' @importFrom stats nlminb
#' @import RTMB
#' @export Get_Reference_Points
#' @family Reference Points and Projections
Get_Reference_Points <- function(data,
                                 rep,
                                 SPR_x = NULL,
                                 t_spawn = 0,
                                 sex_ratio_f = array(0.5, dim = c(data$n_pop, data$n_regions)),
                                 calc_rec_st_yr = 1,
                                 rec_age = 1,
                                 type,
                                 what,
                                 n_avg_yrs = 1,
                                 local_bh_msy_newton_steps = 6,
                                 is_discard_fleet = array(0, dim = data$n_fish_fleets)
                                 ) {

  if(all(is_discard_fleet == 1)) stop("is_discard_fleet is all 1's! At least one fleet needs to be a retention fleet (0).")

  # Dimensions
  n_years <- length(data$years) # number of years
  n_ages <- length(data$ages) # number of ages
  n_pop <- data$n_pop # number of populations
  n_seas <- data$n_seas # number of populations
  n_regions <- data$n_regions # number of regions
  n_fish_fleets <- data$n_fish_fleets # number of fleets

  f_ref_pt <- vector()
  virgin_pop_b_ref_pt <- pop_b_ref_pt <- virgin_b_ref_pt <- b_ref_pt <- array(0, dim = c(n_pop, n_regions))

  # determine years to average over demographics
  n_yrs <- length(data$years)
  avg_yrs <- (n_yrs - n_avg_yrs + 1):n_yrs

  # precompute number of populations per natal region (used throughout for mass-balance stray scaling)
  n_pop_in_region <- rep(0, n_regions)
  for(p in 1:n_pop) n_pop_in_region[data$natal_region[p]] <- n_pop_in_region[data$natal_region[p]] + 1

  # movement / mortality sequencing; absent for input lists built before this option existed
  move_timing <- if(is.null(data$move_timing)) 0 else data$move_timing
  if(move_timing == 2 && is.null(rep$Mrate))
    stop("move_timing == 2 (continuous movement) requires the CTMC generator, but rep$Mrate is NULL. ",
         "Continuous movement reference points need move_type == 1 with use_fixed_movement == 0.")

  # Average the generator on the RATE scale when movement is continuous. Averaging the
  # exponentiated movement fractions instead would be wrong, since mean(expm(Q)) != expm(mean(Q)),
  # and the discrepancy grows with movement rate. Returns a zero placeholder of matching
  # shape when movement is discrete, so downstream data lists always carry the field.
  avg_Mrate <- function(pop_idx, out_dim) {
    if(move_timing != 2) return(array(0, dim = out_dim))
    array(apply(rep$Mrate[pop_idx,,,avg_yrs,,,1,drop = FALSE], c(1,2,3,5,6), mean), dim = out_dim)
  }

  if(type == "single_region") {

    if(n_regions > 1) stop("Single region reference points specified, but n_regions > 1!")

    data_list <- list() # set up data list

    # Seasonal stuff
    data_list$t_spawn <- t_spawn # specified mortality time up until spawning
    data_list$n_seas <- data$n_seas # number of seasons
    data_list$seasdur <- data$seasdur # seasonal duration
    data_list$spawn_seas <- data$spawn_seas # spawning season

    # Dimensions
    data_list$n_pop <- n_pop
    data_list$n_ages <- n_ages

    if(!what %in% c("SPR", "BH_MSY")) stop("what is not correctly specified! Should be SPR, BH_MSY for type = single_region")

    # setup shared data lists
    data_list$F_fract_flt <- array(rep$Fmort[1,n_years,,] / sum(rep$Fmort[1,n_years,,]), dim = c(data$n_seas, data$n_fish_fleets)) # fishing mortality fraction
    data_list$dmr <- array(rep$dmr[1,n_years,,], dim = c(data$n_seas, data$n_fish_fleets)) # discard mortality rate (fraction)
    fish_sel_avg <- apply(rep$fish_sel[, 1, avg_yrs, , , 1, , drop = FALSE], c(1, 4, 5, 7), mean)
    data_list$fish_sel <- array(fish_sel_avg, dim = c(n_pop, n_seas, n_ages, data$n_fish_fleets)) # get female selectivity for all fleets
    ret_sel_avg <- apply(rep$ret_sel[, 1, avg_yrs, , , 1, , drop = FALSE], c(1, 4, 5, 7), mean)
    data_list$ret_sel <- array(ret_sel_avg, dim = c(n_pop, n_seas, n_ages, data$n_fish_fleets)) # get female retention selectivity for all fleets
    natmort_avg <- apply(rep$natmort[,1,avg_yrs,,1,drop = FALSE], c(1,4), mean)
    data_list$natmort <- natmort_avg # get female natural mortality
    WAA_avg <- apply(data$WAA[,1,avg_yrs,,,1,drop = FALSE], c(1,4,5), mean)
    data_list$WAA <- array(WAA_avg, dim = c(n_pop, data$n_seas, n_ages)) # weight at age for females
    MatAA_avg <- apply(data$MatAA[,1,avg_yrs,,,1,drop = FALSE], c(1,4,5), mean)
    data_list$MatAA <- array(MatAA_avg, dim = c(n_pop, data$n_seas, n_ages)) # maturity at age for females

    # Other recruitment stuff
    data_list$sex_ratio_f <- sex_ratio_f # recruitment sex ratio
    data_list$stray_rate <- array(apply(rep$stray_rate[,avg_yrs, drop = FALSE], 1, mean), dim = n_pop) # stray rate
    data_list$rec_seas_prop <- array(rep$rec_seas_prop[,], dim = c(n_pop, data$n_seas)) # recruitment seasonal proportion
    data_list$natal_region <- data$natal_region
    data_list$n_pop_in_region <- n_pop_in_region

    if(what == 'SPR') {

      data_list$SPR_x <- SPR_x # SPR fraction
      par_list <- list() # set up parameter list
      par_list$log_F_x <- log(0.1) # F_x starting value

      # Make adfun object
      tmp_obj <- optim_ref_pts(single_region_SPR, data_list, par_list)

      # Output reference points
      f_ref_pt[1] <- tmp_obj$rep$F_x

      # Compute population specific reference points, by using stray rates
      mean_rec <- apply(rep$Rec[,1,calc_rec_st_yr:(n_years-rec_age),drop=FALSE], 1, mean)

      for(p2 in 1:n_pop) {
        r <- data$natal_region[p2]
        pop_b_ref_pt[p2,1]        <- tmp_obj$rep$SB[p2]  * mean_rec[p2]
        virgin_pop_b_ref_pt[p2,1] <- tmp_obj$rep$SB0[p2] * mean_rec[p2]
        for(p in 1:n_pop) {
          if(p != p2) {
            sc <- data_list$stray_rate[p] / n_pop_in_region[r]
            pop_b_ref_pt[p2,1]        <- pop_b_ref_pt[p2,1]        + sc * tmp_obj$rep$SB[p]  * mean_rec[p]
            virgin_pop_b_ref_pt[p2,1] <- virgin_pop_b_ref_pt[p2,1] + sc * tmp_obj$rep$SB0[p] * mean_rec[p]
          }
        }
      }

      # Compute global reference points (sum across populations)
      b_ref_pt[,1] <- tmp_obj$rep$SB * apply(rep$Rec[,1,calc_rec_st_yr:(n_years - rec_age), drop = F], 1, mean)
      virgin_b_ref_pt[,1] <- tmp_obj$rep$SB0 * apply(rep$Rec[,1,calc_rec_st_yr:(n_years - rec_age), drop = F], 1, mean)

    } # end SPR reference points

    if(what == 'BH_MSY') {

      # extract out beverton-holt parameters
      data_list$h <- array(rep$h_trans[,1], dim = n_pop) # steepness
      data_list$R0 <- array(rep$R0, dim = n_pop) # unfished recruitment
      data_list$is_discard_fleet <- is_discard_fleet # discarding

      par_list <- list() # set up parameter list
      par_list$log_Fmsy <- log(0.1) # Fmsy starting value

      # make adfun etc
      tmp_obj <- optim_ref_pts(single_region_BH_Fmsy, data_list, par_list)

      # Output reference points
      f_ref_pt[1] <- tmp_obj$rep$Fmsy

      # Accumulate biomass reference points
      for(p2 in 1:n_pop) {
        r <- data$natal_region[p2]
        pop_b_ref_pt[p2,1]        <- tmp_obj$rep$SB[p2]  * tmp_obj$rep$Req[p2]
        virgin_pop_b_ref_pt[p2,1] <- tmp_obj$rep$SB0[p2] * rep$R0[p2]
        for(p in 1:n_pop) {
          if(p != p2) {
            sc <- data_list$stray_rate[p] / n_pop_in_region[r]
            pop_b_ref_pt[p2,1]        <- pop_b_ref_pt[p2,1]        + sc * tmp_obj$rep$SB[p]  * tmp_obj$rep$Req[p]
            virgin_pop_b_ref_pt[p2,1] <- virgin_pop_b_ref_pt[p2,1] + sc * tmp_obj$rep$SB0[p] * data_list$R0[p]
          }
        }
      }

      b_ref_pt[,1] <- tmp_obj$rep$SB * tmp_obj$rep$Req
      virgin_b_ref_pt[,1] <- tmp_obj$rep$SB0 * data_list$R0

    }
  }

  if(type == 'multi_region') {

    if(!what %in% c("independent_SPR", "independent_BH_MSY", "global_SPR", "global_BH_MSY", "local_BH_MSY"))
      stop("what is not correctly specified! Should be independent_SPR, independent_BH_MSY, global_SPR, global_BH_MSY, local_BH_MSY for type = multi_region")

    data_list <- list() # set up data list

    # Seasonal stuff
    data_list$t_spawn <- t_spawn # specified mortality time up until spawning
    data_list$n_seas <- data$n_seas # number of seasons
    data_list$seasdur <- data$seasdur # seasonal duration
    data_list$spawn_seas <- data$spawn_seas # spawning season

    # Dimensions
    data_list$n_pop <- n_pop
    data_list$n_ages <- n_ages
    data_list$n_regions <- n_regions

    if(what %in% c("independent_SPR", 'independent_BH_MSY')) {

      tmp_obj <- list() # save optimized object as a list

      for(r in 1:data$n_regions) {

        # create shared data lists
        data_list$F_fract_flt <- array(rep$Fmort[r,n_years,,] / sum(rep$Fmort[r,n_years,,]), dim = c(data$n_seas, data$n_fish_fleets)) # get fleet F fraction to derive population level selectivity
        data_list$dmr <- array(rep$dmr[r,n_years,,], dim = c(data$n_seas, data$n_fish_fleets)) # get dmr
        fish_sel_avg <- apply(rep$fish_sel[, r, avg_yrs, , , 1, , drop = FALSE], c(1, 4, 5, 7), mean)
        data_list$fish_sel <- array(fish_sel_avg, dim = c(n_pop, n_seas, n_ages, data$n_fish_fleets)) # get female total selectivity for all fleets
        ret_sel_avg <- apply(rep$ret_sel[, r, avg_yrs, , , 1, , drop = FALSE], c(1, 4, 5, 7), mean)
        data_list$ret_sel <- array(ret_sel_avg, dim = c(n_pop, n_seas, n_ages, data$n_fish_fleets)) # get female retained selectivity for all fleets
        natmort_avg <- apply(rep$natmort[,r,avg_yrs,,1,drop = FALSE], c(1,4), mean)
        data_list$natmort <- array(natmort_avg, dim = c(n_pop, n_ages)) # get female natural mortality
        WAA_avg <- apply(data$WAA[,r,avg_yrs,,,1,drop = FALSE], c(1,4,5), mean)
        data_list$WAA <- array(WAA_avg, dim = c(n_pop, data$n_seas, n_ages)) # weight at age for females
        MatAA_avg <- apply(data$MatAA[,r,avg_yrs,,,1,drop = FALSE], c(1,4,5), mean)
        data_list$MatAA <- array(MatAA_avg, dim = c(n_pop, data$n_seas, n_ages)) # maturity at age for females
        data_list$sex_ratio_f <- array(sex_ratio_f[,r], dim = n_pop) # recruitment sex ratio
        data_list$stray_rate <- array(apply(rep$stray_rate[,avg_yrs, drop = FALSE], 1, mean), dim = n_pop) # stray rate
        data_list$rec_seas_prop <- array(rep$rec_seas_prop, dim = c(n_pop, data$n_seas)) # recruitment seasonal proportion
        data_list$natal_region <- data$natal_region
        data_list$n_pop_in_region <- n_pop_in_region

        if(what == 'independent_SPR') {

          data_list$SPR_x <- SPR_x # SPR fraction

          par_list <- list() # set up parameter list
          par_list$log_F_x <- log(0.1) # F_x starting value

          # Make adfun object
          tmp_obj[[r]] <- optim_ref_pts(single_region_SPR, data_list, par_list)

          # Output reference points
          f_ref_pt[r] <- tmp_obj[[r]]$rep$F_x

          # Compute population specific reference points, by using stray rates
          mean_rec <- apply(rep$Rec[,r,calc_rec_st_yr:(n_years-rec_age),drop=FALSE], 1, mean)
          if(n_pop > 1) {
            for(p2 in 1:n_pop) {
              rn <- data$natal_region[p2]
              pop_b_ref_pt[p2,r]        <- tmp_obj[[r]]$rep$SB[p2]  * mean_rec[p2]
              virgin_pop_b_ref_pt[p2,r] <- tmp_obj[[r]]$rep$SB0[p2] * mean_rec[p2]
              for(p in 1:n_pop) {
                if(p != p2) {
                  sc <- data_list$stray_rate[p] / n_pop_in_region[rn]
                  pop_b_ref_pt[p2,r]        <- pop_b_ref_pt[p2,r]        + sc * tmp_obj[[r]]$rep$SB[p]  * mean_rec[p]
                  virgin_pop_b_ref_pt[p2,r] <- virgin_pop_b_ref_pt[p2,r] + sc * tmp_obj[[r]]$rep$SB0[p] * mean_rec[p]
                } # end if
              } # end p loop
            } # end p2 loop
          }

          # Compute global reference points (sum across populations)
          b_ref_pt[,r] <- tmp_obj[[r]]$rep$SB * apply(rep$Rec[,r,calc_rec_st_yr:(n_years - rec_age), drop = F], 1, mean)
          virgin_b_ref_pt[,r] <- tmp_obj[[r]]$rep$SB0 * apply(rep$Rec[,r,calc_rec_st_yr:(n_years - rec_age), drop = F], 1, mean)

        } # independent SPR

        if(what == 'independent_BH_MSY') {

          # Beverton Holt parameters
          data_list$h <- array(rep$h_trans[,r], dim = n_pop) # steepness
          data_list$R0 <- array(rep$R0 * rep$rec_region_prop[,r], dim = n_pop) # unfished recruitment by region
          data_list$is_discard_fleet <- is_discard_fleet # discarding

          par_list <- list() # set up parameter list
          par_list$log_Fmsy <- log(0.1) # Fmsy starting value

          # optimize model
          tmp_obj[[r]] <- optim_ref_pts(single_region_BH_Fmsy, data_list, par_list)
          f_ref_pt[r] <- tmp_obj[[r]]$rep$Fmsy

          # get and accumulate biomass reference points
          if(n_pop > 1) {
            for(p2 in 1:n_pop) {
              rn <- data$natal_region[p2]
              pop_b_ref_pt[p2,r]        <- tmp_obj[[r]]$rep$SB[p2]  * tmp_obj[[r]]$rep$Req[p2]
              virgin_pop_b_ref_pt[p2,r] <- tmp_obj[[r]]$rep$SB0[p2] * rep$R0[p2]
              for(p in 1:n_pop) {
                if(p != p2) {
                  sc <- data_list$stray_rate[p] / n_pop_in_region[rn]
                  pop_b_ref_pt[p2,r]        <- pop_b_ref_pt[p2,r]        + sc * tmp_obj[[r]]$rep$SB[p]  * tmp_obj[[r]]$rep$Req[p]
                  virgin_pop_b_ref_pt[p2,r] <- virgin_pop_b_ref_pt[p2,r] + sc * tmp_obj[[r]]$rep$SB0[p] * data_list$R0[p]
                }
              }
            }
          }

          b_ref_pt[,r] <-  tmp_obj[[r]]$rep$SB * tmp_obj[[r]]$rep$Req
          virgin_b_ref_pt[,r] <-  tmp_obj[[r]]$rep$SB0 * data_list$R0

        } # independent BH MSY

      } # end r loop

      # sum up biomass reference points to pop-specific quantities
      if(n_pop == 1) {
        pop_b_ref_pt[1,1] = sum(b_ref_pt)
        virgin_pop_b_ref_pt[1,1] = sum(virgin_b_ref_pt)
      }

    } # end independent methods

    # Global SPR
    if(what == 'global_SPR') {

      # create data lists
      fratio <- array(0, dim = c(n_regions, data$n_seas, data$n_fish_fleets))
      terminal_F <- array(rep$Fmort[,n_years,,], dim = dim(fratio))
      for(r in 1:n_regions) for(seas in 1:data$n_seas) for(f in 1:data$n_fish_fleets) fratio[r,seas,f] <- terminal_F[r,seas,f] / sum(terminal_F[r,,])
      data_list$F_fract_flt <- fratio
      data_list$dmr <- array(rep$dmr[,n_years,,], dim = c(data$n_regions, data$n_seas, data$n_fish_fleets)) # get dmr
      fish_sel_avg <- apply(rep$fish_sel[, , avg_yrs, , , 1, , drop = FALSE], c(1, 2, 4, 5, 7), mean)
      data_list$fish_sel <- array(fish_sel_avg, dim = c(n_pop, n_regions, n_seas, n_ages, data$n_fish_fleets)) # get female total selectivity for all fleets
      ret_sel_avg <- apply(rep$ret_sel[, , avg_yrs, , , 1, , drop = FALSE], c(1, 2, 4, 5, 7), mean)
      data_list$ret_sel <- array(ret_sel_avg, dim = c(n_pop, n_regions, n_seas, n_ages, data$n_fish_fleets)) # get female retained selectivity for all fleets
      natmort_avg <- apply(rep$natmort[,,avg_yrs,,1,drop = FALSE], c(1,2,4), mean)
      data_list$natmort <- array(natmort_avg, dim = c(n_pop, n_regions, n_ages)) # get female natural mortality
      WAA_avg <- apply(data$WAA[,,avg_yrs,,,1,drop = FALSE], c(1, 2, 4, 5), mean)
      data_list$WAA <- array(WAA_avg, dim = c(n_pop, n_regions, data$n_seas, n_ages)) # weight at age for females
      MatAA_avg <- apply(data$MatAA[,,avg_yrs,,,1,drop = FALSE], c(1, 2, 4, 5), mean)
      data_list$MatAA <- array(MatAA_avg, dim = c(n_pop, n_regions, data$n_seas, n_ages)) # maturity at age for females
      Movement_avg <- apply(rep$Movement[,,,avg_yrs,,,1,drop = FALSE], c(1,2,3,5,6), mean)
      data_list$Movement <- array(Movement_avg, dim = c(n_pop, n_regions, n_regions, n_seas, n_ages)) # Movement
      data_list$move_timing <- move_timing
      data_list$seasdur <- data$seasdur
      data_list$Mrate <- avg_Mrate(seq_len(n_pop), c(n_pop, n_regions, n_regions, n_seas, n_ages)) # instantaneous rates

      # Recruitment options
      data_list$do_recruits_move <- data$do_recruits_move # whether recruits move
      data_list$sex_ratio_f <- sex_ratio_f # recruitment sex ratio
      data_list$rec_seas_prop <- array(rep$rec_seas_prop[,], dim = c(n_pop, data$n_seas)) # recruitment seasonal proportion
      data_list$stray_rate <- array(apply(rep$stray_rate[,avg_yrs, drop = FALSE], 1, mean), dim = n_pop) # stray rate
      sgl_seas_spawning_movement_avg <- apply(rep$sgl_seas_spawning_movement[,,,avg_yrs,,1,drop = FALSE], c(1,2,3,5), mean)
      data_list$sgl_seas_spawning_movement <- array(sgl_seas_spawning_movement_avg, dim = c(n_pop, n_regions, n_regions, n_ages)) # Movement
      data_list$natal_region <- data$natal_region
      data_list$n_pop_in_region <- n_pop_in_region

      data_list$SPR_x <- SPR_x # SPR fraction
      mean_rec <- apply(rep$Rec[,,calc_rec_st_yr:(n_years-rec_age),drop=FALSE], c(1,2), mean) # [n_pop, n_regions]
      total_mean_rec <- apply(mean_rec, 1, sum) # [n_pop] - total recruitment across regions
      data_list$rec_region_prop <- mean_rec / total_mean_rec # recruitment proportions

      par_list <- list() # set up parameter list
      par_list$log_F_x <- log(0.1) # F_x starting value

      # make adfn object
      tmp_obj <- optim_ref_pts(global_SPR, data_list, par_list)

      # output reference points
      f_ref_pt <- rep(tmp_obj$rep$F_x, n_regions)

      # Region-specific physical SSB
      for(r in 1:n_regions) {
        b_ref_pt[,r]        <- tmp_obj$rep$SB[,r]  * total_mean_rec
        virgin_b_ref_pt[,r] <- tmp_obj$rep$SB0[,r] * total_mean_rec
      }

      # Population-specific effective SSB at natal region
      if(n_pop > 1) {
        for(p2 in 1:n_pop) {
          r_natal <- data$natal_region[p2]
          pop_b_ref_pt[p2, r_natal]        <- tmp_obj$rep$SB[p2, r_natal]  * total_mean_rec[p2]
          virgin_pop_b_ref_pt[p2, r_natal] <- tmp_obj$rep$SB0[p2, r_natal] * total_mean_rec[p2]
          for(p in 1:n_pop) {
            if(p != p2) {
              sc <- data_list$stray_rate[p] / n_pop_in_region[r_natal]
              pop_b_ref_pt[p2, r_natal]        <- pop_b_ref_pt[p2, r_natal]        + sc * tmp_obj$rep$SB[p, r_natal]  * total_mean_rec[p]
              virgin_pop_b_ref_pt[p2, r_natal] <- virgin_pop_b_ref_pt[p2, r_natal] + sc * tmp_obj$rep$SB0[p, r_natal] * total_mean_rec[p]
            }
          }
        }
      } else {
        pop_b_ref_pt[1,1] = sum(b_ref_pt)
        virgin_pop_b_ref_pt[1,1] = sum(virgin_b_ref_pt)
      }

    } # end global SPR

    # Global BH MSY
    if(what == 'global_BH_MSY') {

      # Error out if invalid recruitment density dependent option
      if(n_pop > 1) stop("Invalid reference point option! When n_pop > 1 reference points must either be independent_SPR, independent_BH_MSY, global_SPR, or local_BH_MSY.")

      # create a data list
      fratio <- array(0, dim = c(n_regions, data$n_seas, data$n_fish_fleets))
      terminal_F <- array(rep$Fmort[,n_years,,], dim = dim(fratio))
      for(r in 1:n_regions) for(seas in 1:data$n_seas) for(f in 1:data$n_fish_fleets) fratio[r,seas,f] <- terminal_F[r,seas,f] / sum(terminal_F[r,,])
      data_list$F_fract_flt <- fratio
      data_list$dmr <- array(rep$dmr[,n_years,,], dim = c(data$n_regions, data$n_seas, data$n_fish_fleets)) # get dmr
      fish_sel_avg <- apply(rep$fish_sel[, , avg_yrs, , , 1, , drop = FALSE], c(1, 2, 4, 5, 7), mean)
      data_list$fish_sel <- array(fish_sel_avg, dim = c(n_pop, n_regions, n_seas, n_ages, data$n_fish_fleets)) # get female total selectivity for all fleets
      ret_sel_avg <- apply(rep$ret_sel[, , avg_yrs, , , 1, , drop = FALSE], c(1, 2, 4, 5, 7), mean)
      data_list$ret_sel <- array(ret_sel_avg, dim = c(n_pop, n_regions, n_seas, n_ages, data$n_fish_fleets)) # get female retained selectivity for all fleets
      natmort_avg <- apply(rep$natmort[1,,avg_yrs,,1,drop = FALSE], c(2,4), mean)
      data_list$natmort <- array(natmort_avg, dim = c(n_regions, n_ages)) # get female natural mortality
      WAA_avg <- apply(data$WAA[,,avg_yrs,,,1,drop = FALSE], c(2, 4, 5), mean)
      data_list$WAA <- array(WAA_avg, dim = c(n_regions, data$n_seas, n_ages)) # weight at age for females
      MatAA_avg <- apply(data$MatAA[,,avg_yrs,,,1,drop = FALSE], c(2, 4, 5), mean)
      data_list$MatAA <- array(MatAA_avg, dim = c(n_regions, data$n_seas, n_ages)) # maturity at age for females
      Movement_avg <- apply(rep$Movement[1,,,avg_yrs,,,1,drop = FALSE], c(1,2,3,5,6), mean)
      data_list$Movement <- array(Movement_avg, dim = c(n_regions, n_regions, n_seas, n_ages)) # Movement
      data_list$move_timing <- move_timing
      data_list$seasdur <- data$seasdur
      data_list$Mrate <- avg_Mrate(1, c(n_regions, n_regions, n_seas, n_ages)) # instantaneous rates
      data_list$is_discard_fleet <- is_discard_fleet

      # Recruitment options
      data_list$do_recruits_move <- data$do_recruits_move # whether recruits move
      data_list$rec_region_prop <- rep$rec_region_prop[1,] # recruitment proportions
      data_list$sex_ratio_f <- sex_ratio_f # recruitment sex ratio to use
      data_list$rec_seas_prop <- rep$rec_seas_prop[1,] # seasonal recruitment
      data_list$h <- mean(rep$h_trans[1,]) # steepness
      data_list$R0 <- rep$R0[1]  # unfished recruitment

      par_list <- list() # set up parameter list
      par_list$log_Fmsy <- log(0.1) # Fmsy starting value

      # make adfn object
      tmp_obj <- optim_ref_pts(global_BH_Fmsy, data_list, par_list)

      # Output reference points
      f_ref_pt <- rep(tmp_obj$rep$Fmsy, n_regions)
      b_ref_pt[1,] <- apply(tmp_obj$rep$SB_age[2,,,drop = F], 2, sum) * tmp_obj$rep$Req
      pop_b_ref_pt[1,1] <- sum(b_ref_pt)
      virgin_b_ref_pt[1,] <- apply(tmp_obj$rep$SB_age[1,,,drop = F], 2, sum) * data_list$R0
      virgin_pop_b_ref_pt[1,1]  <- sum(virgin_b_ref_pt)

    } # end global BH MSY

    if(what == 'local_BH_MSY') {

      # create data list
      fratio <- array(0, dim = c(n_regions, data$n_seas, data$n_fish_fleets))
      terminal_F <- array(rep$Fmort[,n_years,,], dim = dim(fratio))
      for(r in 1:n_regions) for(seas in 1:data$n_seas) for(f in 1:data$n_fish_fleets) fratio[r,seas,f] <- terminal_F[r,seas,f] / sum(terminal_F[r,,])
      data_list$F_fract_flt <- fratio
      data_list$dmr <- array(rep$dmr[,n_years,,], dim = c(data$n_regions, data$n_seas, data$n_fish_fleets)) # get dmr
      fish_sel_avg <- apply(rep$fish_sel[, , avg_yrs, , , 1, , drop = FALSE], c(1, 2, 4, 5, 7), mean)
      data_list$fish_sel <- array(fish_sel_avg, dim = c(n_pop, n_regions, n_seas, n_ages, data$n_fish_fleets)) # get female total selectivity for all fleets
      ret_sel_avg <- apply(rep$ret_sel[, , avg_yrs, , , 1, , drop = FALSE], c(1, 2, 4, 5, 7), mean)
      data_list$ret_sel <- array(ret_sel_avg, dim = c(n_pop, n_regions, n_seas, n_ages, data$n_fish_fleets)) # get female retained selectivity for all fleets
      natmort_avg <- apply(rep$natmort[,,avg_yrs,,1,drop = FALSE], c(1,2,4), mean)
      data_list$natmort <- array(natmort_avg, dim = c(if(n_pop > 1) n_pop else NULL, n_regions, n_ages)) # get female natural mortality
      WAA_avg <- apply(data$WAA[,,avg_yrs,,,1,drop = FALSE], c(1, 2, 4, 5), mean)
      data_list$WAA <- array(WAA_avg, dim = c(if(n_pop > 1) n_pop else NULL, n_regions, data$n_seas, n_ages)) # weight at age for females
      MatAA_avg <- apply(data$MatAA[,,avg_yrs,,,1,drop = FALSE], c(1, 2, 4, 5), mean)
      data_list$MatAA <- array(MatAA_avg, dim = c(if(n_pop > 1) n_pop else NULL, n_regions, data$n_seas, n_ages)) # maturity at age for females
      Movement_avg <- apply(rep$Movement[,,,avg_yrs,,,1,drop = FALSE], c(1,2,3,5,6), mean)
      data_list$Movement <- array(Movement_avg, dim = c(if(n_pop > 1) n_pop else NULL, n_regions, n_regions, n_seas, n_ages)) # Movement
      data_list$move_timing <- move_timing
      data_list$seasdur <- data$seasdur
      data_list$Mrate <- avg_Mrate(seq_len(n_pop), c(if(n_pop > 1) n_pop else NULL, n_regions, n_regions, n_seas, n_ages)) # instantaneous rates
      data_list$is_discard_fleet <- is_discard_fleet

      # Recruitment options
      data_list$do_recruits_move <- data$do_recruits_move # whether recruits move
      data_list$rec_region_prop <- array(rep$rec_region_prop, dim = c(if(n_pop > 1) n_pop else NULL, n_regions)) # recruitment proportions
      data_list$sex_ratio_f <- array(sex_ratio_f, dim = c(if(n_pop > 1) n_pop else NULL, n_regions)) # recruitment sex ratio to use
      data_list$rec_seas_prop <- array(rep$rec_seas_prop, dim = c(if(n_pop > 1) n_pop else NULL, data$n_seas)) # seasonal recruitment
      data_list$h <- array(rep$h_trans, dim = c(if(n_pop > 1) n_pop else NULL, n_regions)) # steepness
      data_list$R0 <- rep$R0  # unfished recruitment
      data_list$stray_rate <- array(apply(rep$stray_rate[,avg_yrs, drop = FALSE], 1, mean), dim = data$n_pop) # stray rate
      data_list$newton_steps <- local_bh_msy_newton_steps # number of newton steps to take
      data_list$natal_region <- data$natal_region
      data_list$n_pop_in_region <- n_pop_in_region
      sgl_seas_spawning_movement_avg <- apply(rep$sgl_seas_spawning_movement[,,,avg_yrs,,1,drop = FALSE], c(1,2,3,5), mean)
      data_list$sgl_seas_spawning_movement <- array(sgl_seas_spawning_movement_avg, dim = c(n_pop, n_regions, n_regions, n_ages)) # Movement

      par_list <- list() # set up parameter list
      par_list$log_Fmsy <- rep(log(0.1), n_regions) # Fmsy starting value

      # Make adfun object
      tmp_obj <- optim_ref_pts(if(n_pop == 1) local_BH_Fmsy_sglpop else local_BH_Fmsy_multipop, data_list, par_list)

      # Output reference points
      f_ref_pt <- tmp_obj$rep$Fmsy

      # multi-population reference points
      if(n_pop > 1) {
        for(r in 1:n_regions) {
          b_ref_pt[,r]        <- tmp_obj$rep$SB[,r]  * tmp_obj$rep$Req_o
          virgin_b_ref_pt[,r] <- tmp_obj$rep$SB0[,r] * data_list$R0
        }

        # accumulate stray rates (divided by n_pop_in_region for mass balance)
        for(p2 in 1:n_pop) {
          r_natal <- data$natal_region[p2]
          pop_b_ref_pt[p2, r_natal]        <- tmp_obj$rep$SB[p2, r_natal]  * tmp_obj$rep$Req_o[p2]
          virgin_pop_b_ref_pt[p2, r_natal] <- tmp_obj$rep$SB0[p2, r_natal] * data_list$R0[p2]

          for(p in 1:n_pop) {
            if(p != p2) {
              sc <- data_list$stray_rate[p] / n_pop_in_region[r_natal]
              pop_b_ref_pt[p2, r_natal]        <- pop_b_ref_pt[p2, r_natal]        + sc * tmp_obj$rep$SB[p, r_natal]  * tmp_obj$rep$Req_o[p]
              virgin_pop_b_ref_pt[p2, r_natal] <- virgin_pop_b_ref_pt[p2, r_natal] + sc * tmp_obj$rep$SB0[p, r_natal] * data_list$R0[p]
            }
          }
        }
      } else { # single population reference points
        for(r in 1:n_regions) {
          b_ref_pt[1, r]        <- sum(tmp_obj$rep$SB_fished_mat[, r]   * tmp_obj$rep$Req_o)
          virgin_b_ref_pt[1, r] <- sum(tmp_obj$rep$SB_unfished_mat[, r] * data_list$R0 * data_list$rec_region_prop)
        }
        pop_b_ref_pt[1, 1]        <- sum(b_ref_pt)
        virgin_pop_b_ref_pt[1, 1] <- sum(virgin_b_ref_pt)
      }
      # see if Newton Raphson calcs for equil rec converged
      if(sum(tmp_obj$rep$iter_vec) > 1e-10) warning("Calculations for equilibrium recruits from origin regions might not have converged! Try increasing local_bh_msy_newton_steps or be wary of these values!")
      if(sum(tmp_obj$rep$Fmsy) == sum(exp(par_list$log_Fmsy))) warning("It is unlikely this converged. Starting values of log Fmsy have not changed (specified at log (0.1).")

    }

  } # end multi region

  return(list(f_ref_pt = f_ref_pt,
              b_ref_pt = b_ref_pt,
              virgin_b_ref_pt = virgin_b_ref_pt,
              pop_b_ref_pt = pop_b_ref_pt,
              virgin_pop_b_ref_pt = virgin_pop_b_ref_pt,
              obj = tmp_obj))

}
