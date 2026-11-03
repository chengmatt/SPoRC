# Stage 3 of 3: post fit
#
# Constraints on the fishing mortality a control rule or input gives the projection: a catch cap and a
# stability band, built as functions that Do_Population_Projection applies through F_constraint.

#' Turn Model Quantities Into Plain Numbers
#'
#' Inside the model, RTMB records each calculation so it can differentiate it, and a
#' recorded (AD) quantity cannot be compared with \code{if}. A harvest rule or a constraint
#' that switches on stock status or catch reads a plain copy through this function and
#' compares that instead. The copy has no derivative, so the outcome of the comparison is
#' fixed, while the quantity itself is unchanged and still differentiated wherever the rule
#' uses it in arithmetic. RTMB records the model once, when it is built, so the record keeps
#' whichever side of the comparison was true then. With a projection set up by
#' \code{\link{Setup_Mod_Projection}}, \code{\link{fit_model}} builds the model again at the
#' estimates, so that side is the estimated stock's; in a model of your own, build it again
#' at the estimates yourself. Lists are read element by element, and plain numbers are
#' returned unchanged.
#'
#' @param x A model quantity (a number, vector or array), a list of them, or any plain
#'   R object.
#'
#' @return The same shape, names and dimnames, with plain numbers in place of model
#'   quantities.
#' @export ad_value
#' @examples
#' ad_value(list(a = 1:3, b = NULL)) # unchanged outside the model
ad_value <- function(x) {

  # do list objects
  if(is.list(x)) {
    out <- x
    for(i in seq_along(x)) if(!is.null(x[[i]])) out[[i]] <- ad_value(x[[i]]) # NULL elements are kept, not dropped
    return(out)
  }

  # do regular objects
  if(!inherits(x, "advector")) return(x)
  get_values <- utils::getFromNamespace("getValues", "RTMB")
  values <- get_values(x)
  if(is.null(dim(x))) names(values) <- names(x) else {
    dim(values) <- dim(x)
    dimnames(values) <- dimnames(x)
  }

  return(values)
}


#' Cap The Catch A Projection Takes
#'
#' Builds a constraint for \code{F_constraint} in \code{\link{Do_Population_Projection}}.
#' Where the retained catch the rule's F would take in a year is above the cap, F is
#' lowered until the catch equals it; elsewhere F is left as the rule set it. The
#' comparison is made on the plain numbers (\code{\link{ad_value}}), so the constraint
#' also runs inside the model, which keeps the side of the comparison that was true when
#' it was built. \code{\link{fit_model}} builds it at the estimates.
#'
#' @param cap Numeric. The cap on retained catch in biomass, in the units of
#'   \code{proj_Catch}. Under \code{scope = "region"}, one value for every region and
#'   year, a vector \code{[n_regions]}, or a matrix \code{[n_regions, n_proj_yrs]}
#'   indexed by the year the catch is taken, as \code{catch_input} is. Under
#'   \code{scope = "global"}, one value or a vector \code{[n_proj_yrs]}.
#' @param scope \code{"region"} (default) caps each region's own catch. \code{"global"}
#'   caps the catch summed over regions and shares the capped total across regions in
#'   proportion to the catch the rule's F would have taken, so the rule's regional
#'   pattern is kept. Pairs with \code{fmort_opt = "HCR_global"}, but runs under any rule.
#'
#' @return A function \code{(F, state)} for \code{F_constraint}.
#' @export catch_cap
#' @family Reference Points and Projections
#' @examples
#' \dontrun{
#' # a 25 kt cap on every region and year
#' Do_Population_Projection(..., F_constraint = catch_cap(25))
#' # a cap on the total catch, compared against one global stock status
#' Do_Population_Projection(..., fmort_opt = "HCR_global",
#'                          F_constraint = catch_cap(60, scope = "global"))
#' }
catch_cap <- function(cap, scope = "region") {

  # some preliminary checks
  if(!scope %in% c("region", "global")) stop("scope should be 'region' or 'global'.")
  if(!is.numeric(cap) || any(!is.finite(cap)) || any(cap < 0)) stop("cap should be finite, non-negative catch in biomass.")
  if(scope == "global" && is.matrix(cap)) stop("Under scope = 'global', cap is one value or a vector over projection years, not a matrix by region.")

  # nested function within for closure - catch_cap returns a function that has an environmnet w/ memory of cap and scope
  function(F, state) {

    "[<-" <- RTMB::ADoverload("[<-")

    # figure out some dimensioning stuff
    n_regions <- length(F)

    # get catch
    catch <- state$catch_at_F(F) # retained catch this year at the incoming F, by region
    catch_values <- ad_value(catch) # decisions are taken on the values

    # one total, shared across regions as the incoming F would have taken it
    if(scope == "global") {
      if(length(cap) > 1 && state$y > length(cap)) stop(paste0("cap has ", length(cap), " years, but the projection reached year ", state$y, "."))
      cap_y <- if(length(cap) == 1) cap else cap[state$y]
      if(sum(catch_values) <= cap_y) return(F)
      return(state$F_at_catch(catch * cap_y / sum(catch), F)) # the shares stay on the tape
    }

    # one cap per region, solved jointly; lowering one region moves fish into its neighbors, so
    # re-check until no further region is over
    cap_y <- constraint_year_value(cap, state$y, n_regions, "cap")
    over <- catch_values > cap_y
    if(!any(over)) return(F)
    F_out <- F
    for(pass in 1:n_regions) {
      target <- rep(NA, n_regions)
      target[over] <- cap_y[over] # regions over the cap are solved to it, the rest keep F
      F_out <- state$F_at_catch(target, F)
      catch_values <- ad_value(state$catch_at_F(F_out))
      newly_over <- !over & catch_values > cap_y * (1 + 1e-4) # the solver's own tolerance
      if(!any(newly_over)) break
      over <- over | newly_over
    } # end pass loop

    return(F_out)
  }
}


#' Keep A Projection's Catch Within some percentage of Last Year's
#'
#' Builds a constraint for \code{F_constraint} in \code{\link{Do_Population_Projection}}.
#' Where the retained catch the rule's F would take in a year is more than \code{down}
#' below or \code{up} above last year's catch, F is lowered or raised until the catch sits
#' at the nearer edge of that band; elsewhere F is left as the rule set it. The
#' comparison is made on the plain numbers (\code{\link{ad_value}}), so the constraint
#' also runs inside the model, which keeps the side of the comparison that was true when
#' it was built. \code{\link{fit_model}} builds it at the estimates.
#'
#' @param down Fraction between 0 and 1. The most the catch may fall below last year's.
#'   Default 0.2, a catch no lower than 80 percent of last year's.
#' @param up Non-negative fraction. The most the catch may rise above last year's.
#'   Default 0.2. \code{Inf} leaves the rise unbounded, so only the floor applies.
#' @param scope \code{"region"} (default) holds each region's own catch to the band of
#'   its own last year's catch. \code{"global"} holds the catch summed over regions to
#'   the band of last year's total, and shares the allowed total across regions in
#'   proportion to the catch the rule's F would have taken.
#'
#' @details Last year's catch is read from the projection, so in projection year 2 it is
#'   the terminal year's replayed catch. A region (or, under \code{"global"}, a
#'   projection) that took no catch last year has no catch to hold to, and is left at the
#'   rule's F. Applied after \code{catch_cap} in a list, a floor can raise the catch back
#'   above the cap, so put the cap last when nothing is to exceed it.
#'
#' @return A function \code{(F, state)} for \code{F_constraint}.
#' @export catch_stability
#' @family Reference Points and Projections
#' @examples
#' \dontrun{
#' # advice may change by at most 15 percent a year, and never exceed 25 kt
#' Do_Population_Projection(..., F_constraint = list(catch_stability(down = 0.15, up = 0.15),
#'                                                   catch_cap(25)))
#' }
catch_stability <- function(down = 0.2, up = 0.2, scope = "region") {

  # preliminary checks
  if(!scope %in% c("region", "global")) stop("scope should be 'region' or 'global'.")
  if(!is.numeric(down) || length(down) != 1 || !is.finite(down) || down < 0 || down > 1) stop("down should be a fraction between 0 and 1, the most the catch may fall below last year's.")
  if(!is.numeric(up) || length(up) != 1 || is.na(up) || up < 0) stop("up should be a non-negative fraction, the most the catch may rise above last year's, or Inf for no limit on the rise.")

  function(F, state) {

    "[<-" <- RTMB::ADoverload("[<-")
    n_regions <- length(F)
    catch <- state$catch_at_F(F) # retained catch this year at the incoming F, by region
    last <- state$last_catch # retained catch last year, by region
    catch_values <- ad_value(catch) # decisions are taken on the values
    last_values <- ad_value(last)

    # one total, kept to the band of last year's total
    if(scope == "global") {
      total <- sum(catch_values)
      last_total <- sum(last_values)
      if(last_total <= 0) return(F) # nothing taken last year, so nothing to hold to
      if(total < (1 - down) * last_total) return(state$F_at_catch(catch * ((1 - down) * sum(last)) / sum(catch), F)) # raised to the floor
      if(total > (1 + up) * last_total) return(state$F_at_catch(catch * ((1 + up) * sum(last)) / sum(catch), F)) # lowered to the ceiling
      return(F)
    }

    # one band per region, open where last year's catch was zero; solved jointly, and re-checked
    # since moving one region's F moves fish into its neighbors
    lower <- (1 - down) * last
    upper <- (1 + up) * last
    lower_values <- (1 - down) * last_values
    upper_values <- (1 + up) * last_values
    outside <- last_values > 0 & (catch_values < lower_values | catch_values > upper_values)
    if(!any(outside)) return(F)

    F_out <- F
    for(pass in 1:n_regions) {
      target <- as.numeric(rep(NA, n_regions)) # numeric, since a model quantity written into a logical NA comes back as complex numbers
      for(r in which(outside)) { # the nearer edge of the band
        if(abs(catch_values[r] - lower_values[r]) < abs(catch_values[r] - upper_values[r])) target[r] <- lower[r]
        else target[r] <- upper[r]
      } # end r loop
      F_out <- state$F_at_catch(target, F) # figure out what f is at a given catch
      catch_values <- ad_value(state$catch_at_F(F_out)) # figure out what catch it actually gives
      newly_outside <- last_values > 0 & !outside & (catch_values < lower_values * (1 - 1e-4) | catch_values > upper_values * (1 + 1e-4)) # check tolerance
      if(!any(newly_outside)) break
      outside <- outside | newly_outside
    } # end pass loop

    return(F_out)
  }
}


#' The Value A Constraint Reads In One Projection Year
#'
#' @param x One value, a vector \code{[n_regions]}, or a matrix \code{[n_regions, n_proj_yrs]}.
#' @param y Integer. Projection year.
#' @param n_regions Integer. Number of regions.
#' @param name Character. What \code{x} is, for error messages.
#' @return Numeric \code{[n_regions]}.
#' @keywords internal
#' @noRd
constraint_year_value <- function(x, y, n_regions, name) {
  if(is.matrix(x)) {
    if(nrow(x) != n_regions) stop(paste0(name, " has ", nrow(x), " rows, but the projection has ", n_regions, " regions."))
    if(y > ncol(x)) stop(paste0(name, " has ", ncol(x), " years, but the projection reached year ", y, "."))
    return(as.numeric(x[,y]))
  }
  if(length(x) == 1) return(rep(x, n_regions))
  if(length(x) == n_regions) return(as.numeric(x))
  stop(paste0(name, " should be one value, a vector [n_regions], or a matrix [n_regions, n_proj_yrs], not length ", length(x), " for ", n_regions, " regions."))
}


#' Apply The Constraints To One Projection Year's Fishing Mortality
#'
#' Builds what a constraint reads: the rule's F, last year's catch and spawning biomass,
#' this year's starting numbers, the catch this year's F would take, and the F that takes
#' a given catch. Then applies each constraint in turn, each reading the F the one before
#' it returned. The catch is solved on the values and the solved F then takes two Newton
#' steps, on the tape when the projection is taped, so it has the implicit derivative.
#'
#' @param y Integer. Projection year the F applies to, 2 or later.
#' @param F_rule Numeric \code{[n_regions]}. Annual F from the rule or input.
#' @param F_constraint List of constraint functions \code{(F, state)}.
#' @param seas_share Numeric matrix \code{[n_regions, n_seas]}. Terminal year seasonal
#'   shares of annual F.
#' @param state,tmp_rec,proj_args Projection state and inputs, as
#'   \code{solve_proj_year_F} takes them.
#' @param catch_f_max,catch_tol,catch_max_iter Solver settings.
#'
#' @return Numeric \code{[n_regions]}. The constrained annual F.
#' @keywords internal
#' @noRd
constrain_proj_F <- function(y, F_rule, F_constraint, seas_share, state, tmp_rec, proj_args,
                             catch_f_max, catch_tol, catch_max_iter) {

  "[<-" <- RTMB::ADoverload("[<-")
  n_pop <- proj_args$n_pop
  n_regions <- proj_args$n_regions
  n_seas <- proj_args$n_seas

  # the solve runs on non-ad values and an AD state gets its derivative back from the newton steps below
  state_values <- ad_value(state)
  tmp_rec_values <- ad_value(tmp_rec)
  proj_args_values <- ad_value(proj_args)
  seas_share_values <- ad_value(seas_share)

  # retained catch in biomass this year at an annual F by region, summed over populations, seasons and fleets
  catch_at_F <- function(F) {
    if(length(F) != n_regions) stop(paste0("catch_at_F takes an F per region, length ", n_regions, ", not ", length(F), "."))
    F_y <- array(F * seas_share, dim = c(n_regions, n_seas)) # annual F split over seasons at the terminal year shares
    return(rowSums(proj_catch_at_F(F_y, y, state, tmp_rec, proj_args)))
  }

  # the annual F by region that takes target where it is given, with the other regions kept at F_keep
  F_at_catch <- function(target, F_keep) {
    "[<-" <- RTMB::ADoverload("[<-")
    target_values <- ad_value(target)
    F_keep_values <- ad_value(F_keep)
    if(length(target) != n_regions || length(F_keep) != n_regions) stop(paste0("F_at_catch takes a target and an F per region, length ", n_regions, "."))
    if(any(!is.na(target_values) & (!is.finite(target_values) | target_values < 0))) stop("F_at_catch targets should be NA, to keep F, or finite non-negative catch.")
    solved <- which(!is.na(target_values))
    F_out <- F_keep
    if(length(solved) == 0) return(F_out)

    # regions without a target fish at F_keep; regions with one get their F from the solve
    F_base <- array(F_keep_values * seas_share_values, dim = c(n_regions, n_seas))
    F_base[solved,] <- 0
    target_all <- rep(0, n_regions)
    target_all[solved] <- target[solved] # on the tape when the target is, as a global share is
    target_all_values <- ad_value(target_all)

    # solved on the values
    sol <- solve_proj_F_catch(y = y,
                              target = target_all_values,
                              seas_profile = seas_share_values,
                              F_base = F_base,
                              target_seas = NULL,
                              state = state_values,
                              tmp_rec = tmp_rec_values,
                              proj_args = proj_args_values,
                              f_start = F_keep_values,
                              catch_f_max = catch_f_max,
                              catch_tol = catch_tol,
                              catch_max_iter = catch_max_iter,
                              label = paste0("F_constraint in projection year ", y))

    # newton steps from the solved F with the other regions at F_keep as the model has them, so the F moves
    # with the estimates
    F_base_model <- array(F_keep * seas_share, dim = c(n_regions, n_seas))
    F_base_model[solved,] <- 0
    F_reg <- refine_proj_F_catch(sol = sol, y = y, target = target_all, seas_profile = seas_share, F_base = F_base_model,
                                 target_seas = NULL, state = state, tmp_rec = tmp_rec, proj_args = proj_args)
    F_out[solved] <- F_reg[solved]

    return(F_out)
  }

  # last year's retained catch and spawning biomass, this year's starting numbers
  last_catch <- rep(0, n_regions)
  for(r in 1:n_regions) last_catch[r] <- sum(state$proj_Catch[,r,y - 1,,])
  constraint_state <- list(y = y,
                           F_rule = F_rule,
                           last_catch = last_catch,
                           SSB = array(state$proj_SSB[,,y - 1], dim = c(n_pop, n_regions)),
                           NAA = array(state$proj_NAA[,,y,1,,], dim = c(n_pop, n_regions, proj_args$n_ages, proj_args$n_sexes)),
                           catch_at_F = catch_at_F,
                           F_at_catch = F_at_catch)

  # each constraint reads the F the one before it returned; checked on the values
  F_out <- F_rule
  for(i in seq_along(F_constraint)) {
    F_out <- F_constraint[[i]](F_out, constraint_state)
    F_out_values <- ad_value(F_out)
    if(!is.numeric(F_out_values) || length(F_out_values) != n_regions || any(!is.finite(F_out_values)) || any(F_out_values < 0))
      stop(paste0("F_constraint[[", i, "]] should return finite, non-negative F of length ", n_regions, " in projection year ", y, "."))
    if(!inherits(F_out, "advector")) F_out <- as.numeric(F_out)
  } # end i loop

  return(F_out)
} # end constrain_proj_F
