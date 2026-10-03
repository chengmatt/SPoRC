# Stage 1 of 3: model setup
#
# Catchability deviations: each year's deviation from a fleet's block catchability, under iid,
# random walk, ar1 or a dsem. Used by both the fishery and the survey setups.

#' Deviation arrays a dsem may set from its arrows alone
#'
#' A dsem series with an sd of zero is a function of its covariates rather than
#' a deviation with a density, so something has to write it into the model. The
#' density sits just before the observation models, so catchability is still
#' ahead of it; growth, movement, recruitment and the numbers at age are all
#' consumed before it, so a zero sd on one of those is refused. This is the one
#' list, read by that refusal and by the objective.
#'
#' @return Character vector of parameter array names.
#'
#' @keywords internal
#' @export
q_dev_par_names <- function() c("ln_fish_q_devs", "ln_srv_q_devs")

#' Map catchability deviations, their sigma and their correlation
#'
#' A deviation is estimated where its fleet has a continuous catchability model
#' and that region and fleet hold index observations. Sigma and rho follow the
#' sharing strings over region and fleet, blanked wherever the deviations are.
#'
#' @param input_list List with \code{data}, \code{par} and \code{map}.
#' @param q_model Integer vector \code{[n_fleets]} of process error codes.
#' @param sigma_q_spec Sharing string for the deviation sigma, one of
#'   \code{"est_all"}, \code{"est_shared_r"}, \code{"est_shared_f"},
#'   \code{"est_shared_r_f"} or \code{"fix"}.
#' @param q_rho_spec Sharing string for the ar1 correlation, same options.
#' @param prefix \code{"fish"} or \code{"srv"}.
#' @param fleet_field Name of the fleet count in \code{data}.
#' @param use_field Name of the index use array in \code{data}.
#' @param q_re_years List \code{[n_fleets]}, each element the model years (as
#'   indices into \code{data$years}) that fleet's deviations are estimated
#'   over, or \code{NULL} for every year. \code{NULL} (default) gives every
#'   fleet every year. Projection years are always estimated. A year left out
#'   holds its deviation fixed at its starting value (zero unless set by
#'   \code{...}), which a random walk or ar1 then steps through like any
#'   other fixed cell, except at a fleet's first estimated year: see
#'   \code{\link{Get_q_dev_penalty}}.
#'
#' @return \code{input_list} with the three maps set.
#'
#' @keywords internal
do_q_devs_mapping <- function(input_list, q_model, sigma_q_spec, q_rho_spec, prefix, fleet_field, use_field, q_re_years = NULL) {

  n_regions <- input_list$data$n_regions
  n_fleets <- input_list$data[[fleet_field]]
  devs_name <- paste0("ln_", prefix, "_q_devs")
  use_arr <- input_list$data[[use_field]]
  use_pop_arr <- input_list$data[[paste0(use_field, "_pop")]]

  # only estimate devs for region with data, and only over a fleet's own active years
  map_devs <- input_list$par[[devs_name]]
  map_devs[] <- NA
  n_yr <- dim(map_devs)[2] # fit years plus any projection years, which stay active regardless of q_re_years
  n_yr_fit <- length(input_list$data$years)
  proj_years <- seq_len(n_yr)[-seq_len(n_yr_fit)]
  dev_counter <- 0
  for(f in 1:n_fleets) {
    if(q_model[f] == 1) next # "none", so this fleet has no deviations at all
    active_years <- if(is.null(q_re_years) || is.null(q_re_years[[f]])) seq_len(n_yr) else sort(union(q_re_years[[f]], proj_years))
    for(r in 1:n_regions) {
      has_data <- sum(use_arr[r,,,f]) > 0 || (!is.null(use_pop_arr) && sum(use_pop_arr[,r,,,f]) > 0)
      if(!has_data) next
      map_devs[r,active_years,f] <- dev_counter + seq_len(length(active_years))
      dev_counter <- dev_counter + length(active_years)
    } # end r loop
  } # end f loop

  input_list$map[[devs_name]] <- factor(map_devs)

  # the mirror the penalty reads, refreshed by sync_dev_map_data and blanked where a dsem takes over
  input_list$data[[paste0("map_", devs_name)]] <- array(as.numeric(factor(map_devs)), dim = dim(map_devs))

  dims <- c(region = n_regions, fleet = n_fleets)
  abbrev <- c(r = "region", f = "fleet")
  sigma_name <- paste0("ln_sigma_", prefix, "_q")
  rho_name <- paste0(prefix, "_q_rho")

  input_list$map[[sigma_name]] <- build_shared_spec_map(dims = dims, spec = sigma_q_spec, dim_abbrev = abbrev)
  input_list$map[[rho_name]] <- build_shared_spec_map(dims = dims, spec = q_rho_spec, dim_abbrev = abbrev)

  # a fleet with no deviations reads no sigma, and only an ar1 reads a correlation
  sigma_map <- as.integer(input_list$map[[sigma_name]])
  rho_map <- as.integer(input_list$map[[rho_name]])
  dim(sigma_map) <- dims
  dim(rho_map) <- dims

  has_devs <- apply(!is.na(map_devs), c(1, 3), any) # region by fleet, TRUE where a year series is estimated

  for(f in 1:n_fleets) {

    if(q_model[f] %in% c(1, 5)) sigma_map[,f] <- NA # none, or a dsem whose arrows hold the sd instead (none or dsem)
    if(q_model[f] != 4) rho_map[,f] <- NA # if not ar1

    # a region with no index data has no series for a sigma or a correlation to describe, and the penalty skips it
    sigma_map[!has_devs[,f], f] <- NA
    rho_map[!has_devs[,f], f] <- NA

  } # end f loop

  input_list$map[[sigma_name]] <- factor(renumber_map_levels(sigma_map))
  input_list$map[[rho_name]] <- factor(renumber_map_levels(rho_map))

  return(input_list)

} # end function

#' Renumber a map so its levels run from one with no gaps
#'
#' Blanking cells of a shared map leaves holes in the level numbering, which
#' \code{factor} keeps as unused levels.
#'
#' @param x Integer array holding map levels and \code{NA}.
#'
#' @return The same array with its levels renumbered.
#'
#' @keywords internal
renumber_map_levels <- function(x) {
  keep <- !is.na(x)
  if(!any(keep)) return(x)
  x[keep] <- match(x[keep], sort(unique(x[keep])))
  return(x)
}

#' Set up catchability deviations for a fishery or survey
#'
#' Builds the deviation array, its sigma and correlation, and refuses the
#' combinations that cannot be identified. Called by
#' \code{\link{Setup_Mod_Fishsel_and_Q}} and
#' \code{\link{Setup_Mod_Srvsel_and_Q}}.
#'
#' @param input_list List with \code{data}, \code{par} and \code{map}. The
#'   catchability blocks, the index use arrays and \code{q_type} must already
#'   be set.
#' @param q_model Character vector \code{[n_fleets]}, one of \code{"none"}
#'   (default), \code{"iid"}, \code{"rw"}, \code{"ar1"} or \code{"dsem"}.
#' @param sigma_q_spec Sharing string for the deviation sigma over region and
#'   fleet. Default \code{"est_all"}.
#' @param q_rho_spec Sharing string for the ar1 correlation. Default
#'   \code{"est_all"}.
#' @param q_rw_init_sigma Standard deviation of the first estimated year of a
#'   random walk. \code{NA} (default) starts the walk at zero under its own
#'   sigma, which keeps the block catchability as the level of the series.
#' @param q_type Character vector \code{[n_fleets]} of \code{"est"},
#'   \code{"arith"} or \code{"geo"}, as already validated by the caller.
#' @param prefix \code{"fish"} or \code{"srv"}.
#' @param fleet_field Name of the fleet count in \code{data}.
#' @param use_field Name of the index use array in \code{data}.
#' @param fleet_label Label used in messages.
#' @param starting_values Named list of starting values.
#' @param q_re_years As in \code{\link{do_q_devs_mapping}}.
#'
#' @return \code{input_list} with the deviation parameters, their maps and
#'   \code{<prefix>_q_model} in \code{data}.
#'
#' @keywords internal
setup_q_devs <- function(input_list,
                         q_model,
                         sigma_q_spec,
                         q_rho_spec,
                         q_rw_init_sigma,
                         q_type,
                         prefix,
                         fleet_field,
                         use_field,
                         fleet_label,
                         starting_values,
                         q_re_years = NULL) {

  n_regions <- input_list$data$n_regions
  n_fleets <- input_list$data[[fleet_field]]
  forms <- c("none", "iid", "rw", "ar1", "dsem")
  if(is.null(q_model)) q_model <- rep("none", n_fleets)
  check_fleet_spec_length(q_model, n_fleets, paste0(prefix, "_q_model"))
  if(!all(q_model %in% forms)) stop(prefix, "_q_model should be one of: ", paste(forms, collapse = ", "), ".")
  q_model_val <- match(q_model, forms)
  blocks_arr <- input_list$data[[paste0(prefix, "_q_blocks")]]

  # q_re_years: a list of one integer vector (or NULL) per fleet, indexing data$years
  n_yr_fit <- length(input_list$data$years)
  if(!is.null(q_re_years)) {
    if(!is.list(q_re_years)) stop(prefix, "_q_re_years must be a list of one element per fleet (or NULL).")
    check_fleet_spec_length(q_re_years, n_fleets, paste0(prefix, "_q_re_years"))
    for(f in 1:n_fleets) {
      yrs_f <- q_re_years[[f]]
      if(is.null(yrs_f)) next
      if(!all(yrs_f %in% seq_len(n_yr_fit)) || length(yrs_f) == 0)
        stop(prefix, "_q_re_years[[", f, "]] must index the model years (1 to ", n_yr_fit, ").")
      if(q_model_val[f] == 1) collect_message(prefix, "_q_re_years[[", f, "]] is ignored: ", fleet_label, " ", f, " has ", prefix, "_q_model = 'none'.")
    } # end f loop
  }

  for(f in 1:n_fleets) {
    if(q_model_val[f] == 1) next

    # dont allow analtical q
    if(q_type[f] != "est") stop(prefix, "_q_type is '", q_type[f], "' for ", fleet_label, " ", f, ", which solves catchability from the observations, ",
                                "and ", prefix, "_q_model is '", q_model[f], "'. The solve would absorb the deviations. Use ", prefix,
                                "_q_type = 'est' for that fleet, or ", prefix, "_q_model = 'none'.")

    # make block q and qdevs together invalid
    n_blks <- length(unique(as.vector(blocks_arr[,,f])))
    if(n_blks > 1) stop(prefix, "_q_blocks gives ", fleet_label, " ", f, " ", n_blks, " catchability blocks while ", prefix,
                        "_q_model is '", q_model[f], "'. Both describe time variation in catchability, so use one or the other.")

  } # end f loop

  for(f in 1:n_fleets) collect_message(paste0("Catchability deviations for ", fleet_label, " ", f, " are: ", q_model[f]))

  # a dsem fleet gets its density from the arrows, so its own sigma is fixed and not penalized
  if(any(q_model_val == 5)) {
    input_list$data$dsem_declared <- union(input_list$data$dsem_declared, paste0(prefix, "_q"))
    collect_message(prefix, "_q_model = 'dsem': the catchability deviations' density comes from Setup_Mod_DSEM, and ln_sigma_", prefix, "_q stays at its start.")
  }

  # Parameters --------------------------------------------------------------
  n_dev_yrs <- length(input_list$data$years) + input_list$data$n_proj_yrs_devs
  devs_name <- paste0("ln_", prefix, "_q_devs")
  sigma_name <- paste0("ln_sigma_", prefix, "_q")
  rho_name <- paste0(prefix, "_q_rho")

  input_list$par[[devs_name]] <- array(0, dim = c(n_regions, n_dev_yrs, n_fleets))
  input_list$par[[devs_name]] <- use_starting_value(input_list$par[[devs_name]], starting_values, devs_name)
  input_list$par[[sigma_name]] <- array(log(0.1), dim = c(n_regions, n_fleets))
  input_list$par[[sigma_name]] <- use_starting_value(input_list$par[[sigma_name]], starting_values, sigma_name)
  input_list$par[[rho_name]] <- array(0, dim = c(n_regions, n_fleets))
  input_list$par[[rho_name]] <- use_starting_value(input_list$par[[rho_name]], starting_values, rho_name)

  # Data --------------------------------------------------------------------
  input_list$data[[paste0(prefix, "_q_model")]] <- q_model_val
  input_list$data[[paste0(prefix, "_q_rw_init_sigma")]] <- q_rw_init_sigma

  input_list <- do_q_devs_mapping(input_list, q_model_val, sigma_q_spec, q_rho_spec, prefix, fleet_field, use_field, q_re_years)

  return(input_list)

} # end function
