# Stage 3 of 3: simulation
#
# Catchability deviations in the operating model

#' Set up catchability deviations for the operating model
#'
#' Gives each fleet a process error on annual catchability, matching
#' \code{srv_q_model} and \code{fish_q_model} in the estimation setup.
#' Catchability in replicate \eqn{i} becomes
#' \eqn{q_{r,y,f} \exp(\epsilon_{r,y,f,i})}, where \eqn{q} is the mean
#' \code{Setup_Sim_Fishery} and \code{Setup_Sim_Survey} supplied.
#'
#' Conditioning years read their deviations from the fit and draw nothing, so a
#' random walk or ar1 in the projection continues from the fit's last
#' deviation rather than restarting at zero.
#'
#' @param sim_list List of operating model inputs, with \code{fish_q} and
#'   \code{srv_q} already set.
#' @param fish_q_model,srv_q_model Character vectors, one per fleet, of
#'   \code{"none"} (default), \code{"iid"}, \code{"rw"}, \code{"ar1"} or
#'   \code{"dsem"}. A \code{"dsem"} fleet draws nothing here and takes its
#'   series from \code{\link{Setup_Sim_DSEM}}.
#' @param sigma_fish_q,sigma_srv_q Numeric arrays
#'   \code{[n_regions, n_fleets]} of deviation standard deviations on the log
#'   scale, or a single value used for every region and fleet. Default zero.
#' @param fish_q_rho,srv_q_rho Numeric arrays \code{[n_regions, n_fleets]} of
#'   ar1 correlations on the natural scale, read only under \code{"ar1"}.
#'   Default zero.
#'
#' @return \code{sim_list} with the deviation settings stored.
#'
#' @export
Setup_Sim_q_devs <- function(sim_list,
                             fish_q_model = NULL,
                             srv_q_model = NULL,
                             sigma_fish_q = 0,
                             sigma_srv_q = 0,
                             fish_q_rho = 0,
                             srv_q_rho = 0) {

  forms <- c("none", "iid", "rw", "ar1", "dsem") # catchability forms

  for(prefix in c("fish", "srv")) {

    # figure out if fishery sor survey
    n_fleets <- sim_list[[paste0("n_", prefix, "_fleets")]]
    model <- if(prefix == "fish") fish_q_model else srv_q_model
    sigma <- if(prefix == "fish") sigma_fish_q else sigma_srv_q
    rho <- if(prefix == "fish") fish_q_rho else srv_q_rho

    # if left null no model form (none)
    if(is.null(model)) model <- rep("none", n_fleets)
    if(length(model) != n_fleets) stop(prefix, "_q_model should hold one value per fleet (", n_fleets, ").")
    if(!all(model %in% forms)) stop(prefix, "_q_model should be one of: ", paste(forms, collapse = ", "), ".")
    if(any(rho < -1) || any(rho >= 1)) stop(prefix, "_q_rho must sit inside (-1, 1). An ar1 at one or beyond has no stationary variance.")

    # put stuff into sim list
    dims <- c(sim_list$n_regions, n_fleets)
    sim_list[[paste0(prefix, "_q_model")]] <- match(model, forms)
    sim_list[[paste0("sigma_", prefix, "_q")]] <- array(sigma, dim = dims)
    sim_list[[paste0(prefix, "_q_rho")]] <- array(rho, dim = dims)

  } # end prefix loop

  # conditioning years read catchability from the actual q alues, so a process error there would never be drawn
  n_cond <- if(is.null(sim_list$n_cond_yrs)) 0L else as.integer(sim_list$n_cond_yrs)
  if(isTRUE(n_cond >= sim_list$n_yrs) && any(c(sim_list$fish_q_model, sim_list$srv_q_model) > 1))
    stop("Every one of the ", sim_list$n_yrs, " operating model years is a conditioning year, where catchability ",
         "is read from the fit, so a catchability process error would never be drawn. Run the operating model ",
         "past the conditioning years, or leave the q models at 'none'.")

  return(sim_list)

} # end function

#' Split a fit's reported catchability into its mean and its deviations
#'
#' The reported catchability is \eqn{\exp(\ln q_{blk} + \epsilon)}. The mean is
#' the level a process error walks from, the deviations are what the
#' conditioning years reproduce, and an analytically solved fleet reads no
#' \code{ln_q} so its reported value is already its mean.
#'
#' @param rep_q Reported catchability \code{[region, year, fleet]}.
#' @param ln_q Log block catchability \code{[region, block, fleet]}.
#' @param q_blocks Block index \code{[region, year, fleet]}.
#' @param q_type Integer vector \code{[fleet]}, zero where catchability is
#'   estimated and non-zero where it is solved from the index.
#'
#' @return List of \code{q_mean} and \code{devs}, both
#'   \code{[region, year, fleet]}.
#'
#' @keywords internal
split_reported_q <- function(rep_q, ln_q, q_blocks, q_type = NULL) {

  n_regions <- dim(rep_q)[1]
  n_yrs <- dim(rep_q)[2]
  n_fleets <- dim(rep_q)[3]
  if(is.null(q_type)) q_type <- rep(0, n_fleets)

  q_mean <- rep_q # an analytic fleet is its own mean

  for(f in 1:n_fleets) {
    if(q_type[f] != 0) next
    for(r in 1:n_regions) q_mean[r,,f] <- exp(ln_q[r,q_blocks[r,,f],f])
  } # end f loop

  devs <- log(rep_q) - log(q_mean)
  devs[!is.finite(devs)] <- 0 # a fleet kept at zero catchability has no deviation to read

  return(list(q_mean = q_mean, devs = devs))

} # end function

#' Dont allow a catchability series the operating model cannot draw
#'
#' A dsem series given an sd of zero is derived: it has no innovation of its
#' own, the objective works it out from its covariates, and the fitted
#' deviation parameter stays at zero. \code{\link{draw_dsem_sim}} draws the
#' unknown cells from the field's precision, which a series with no variance
#' has no finite entry in, so every cell it has to draw comes back \code{NaN}.
#'
#' Catchability is the only process a derived series is allowed on, since
#' \code{\link{Setup_Mod_DSEM}} refuses one everywhere else, and the
#' conditioning years are read from the fit rather than drawn. So the one case
#' with cells left to draw is a derived catchability series on an operating
#' model that runs past the conditioning period, and that is what this refuses.
#'
#' @param sim_list Simulation list, or the environment
#'   \code{\link{Setup_sim_env}} built from one.
#'
#' @return Nothing. Called for its refusal.
#'
#' @keywords internal
check_q_dsem_drawable <- function(sim_list) {

  if(is.null(sim_list$dsem_model) || !any(sim_list$dsem_model$derived)) return(invisible(NULL))
  n_cond <- if(is.null(sim_list$n_cond_yrs)) 0L else as.integer(sim_list$n_cond_yrs)

  for(s in seq_along(sim_list$dsem_link_par)) {

    if(!sim_list$dsem_link_par[s] %in% q_dev_par_names()) next
    if(!isTRUE(sim_list$dsem_model$derived[sim_list$dsem_link_col[s]])) next
    n_yrs <- sim_list$n_yrs
    if(is.null(n_yrs) || n_yrs <= n_cond) next

    n_draw <- n_yrs - n_cond
    stop(sim_list$dsem_model$variables[sim_list$dsem_link_col[s]], " has its sd fixed at zero, so it is a derived ",
         "series with no innovation to draw from, and the operating model runs ", n_draw, " year",
         if(n_draw > 1) "s" else "", " past the ", n_cond, " conditioning years. Those years would come back NaN ",
         "rather than a catchability. Give the series an sd line, a small one if the effect is meant to stay a ",
         "regression, or run the operating model over the conditioning years alone, where the series is read from ",
         "the fit rather than drawn.")

  } # end s loop

  return(invisible(NULL))

} # end function

#' Draw one replicate's catchability deviations and apply them
#'
#' Called at the start of each replicate. Conditioning years keep the fit's own
#' deviations, a cell a dsem drew keeps that value, and every other year is this
#' replicate's own innovation, so a random walk or ar1 continues from the fit's
#' last deviation rather than restarting at zero.
#'
#' A fleet a dsem wrote for is drawn whatever its own \code{<prefix>_q_model}
#' says, so the series reaches catchability rather than an array nothing reads.
#'
#' @param sim Replicate index.
#' @param sim_env Simulation environment.
#'
#' @return Nothing. \code{sim_env$fish_q} and \code{sim_env$srv_q} are scaled by
#'   the deviations and the deviations themselves are stored.
#'
#' @keywords internal
draw_sim_q_devs <- function(sim, sim_env) {

  for(prefix in c("fish", "srv")) {

    q_model <- sim_env[[paste0(prefix, "_q_model")]]
    devs_name <- paste0("ln_", prefix, "_q_devs")
    drawn <- sim_env$dsem_drawn[[devs_name]]
    devs_in <- sim_env[[devs_name]]
    if(is.null(q_model) || (all(q_model == 1) && !any(drawn) && all(devs_in == 0))) next

    q_name <- paste0(prefix, "_q")
    base_name <- paste0(prefix, "_q_base")

    # the mean catchability before any deviation, kept so a replicate does not scale the last one
    if(is.null(sim_env[[base_name]])) sim_env[[base_name]] <- sim_env[[q_name]]

    n_regions <- sim_env$n_regions
    n_yrs <- dim(sim_env[[q_name]])[2]
    n_fleets <- dim(sim_env[[q_name]])[3]
    sigma <- sim_env[[paste0("sigma_", prefix, "_q")]]
    rho <- sim_env[[paste0(prefix, "_q_rho")]]
    n_cond <- if(is.null(sim_env$n_cond_yrs)) 0L else min(as.integer(sim_env$n_cond_yrs), n_yrs)

    for(f in 1:n_fleets) {

      fleet_drawn <- !is.null(drawn) && any(drawn[,,f])
      if(q_model[f] == 1 && !fleet_drawn && all(devs_in[,,f,] == 0)) next

      for(r in 1:n_regions) {

        devs <- devs_in[r,,f,sim] # the conditioning years already hold the fit's own deviations

        for(y in seq_len(n_yrs)) {
          if(y <= n_cond) next # read from the fit rather than drawn
          if(fleet_drawn && isTRUE(drawn[r,y,f])) next # a dsem wrote this cell
          if(q_model[f] %in% c(1, 5)) { devs[y] <- 0; next } # no innovation of its own to draw
          dev_mu <- if(q_model[f] == 2 || y == 1) 0 else if(q_model[f] == 3) devs[y - 1] else rho[r,f] * devs[y - 1]
          dev_sd <- if(q_model[f] == 4 && y == 1) sigma[r,f] / sqrt(1 - rho[r,f]^2) else sigma[r,f] # an ar1 starts at its stationary spread
          devs[y] <- stats::rnorm(1, dev_mu, dev_sd)
        } # end y loop

        sim_env[[devs_name]][r,,f,sim] <- devs
        sim_env[[q_name]][r,,f,sim] <- sim_env[[base_name]][r,,f,sim] * exp(devs)

      } # end r loop
    } # end f loop
  } # end prefix loop

  return(invisible(NULL))

} # end function
