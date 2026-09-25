# Stage 3 of 3: post fit
#
# Jitter analysis: refit from perturbed starting values to check the optimizer is finding the same
# minimum.

#' Draw jittered starting values for one jitter iteration
#'
#' \code{nlminb()} only ever sees the fixed effects, so a jittered random effect
#' has to reach the model another way. TMB starts the inner Laplace solve from
#' \code{obj$env$last.par}, so the draws are returned separately and written
#' there by the caller.
#'
#' @param obj An RTMB model object from \code{\link[RTMB]{MakeADFun}}.
#' @param par_vec Starting values, either the fixed-effect vector
#'   (\code{length(obj$par)}) or the joint fixed and random vector
#'   (\code{length(obj$env$par)}), or \code{NULL} for the model's own start.
#' @param sd Standard deviation of the additive normal draws.
#' @param jitter_random Whether the random effects are perturbed as well.
#'
#' @return A list with \code{fixed}, the jittered fixed-effect vector, and
#'   \code{random}, the random-effect starting values in \code{obj$env$random}
#'   order, jittered when \code{jitter_random = TRUE} and taken straight from
#'   \code{par_vec} otherwise. Length zero when the model has no random effects.
#' @keywords internal
#' @importFrom stats rnorm
jitter_start_values <- function(obj, par_vec, sd, jitter_random) {

  n_fixed <- length(obj$par)
  n_joint <- length(obj$env$par)
  rand_idx <- obj$env$random # NULL when the model has no random effects

  # a joint vector splits on the random indices; a fixed-length one leaves the
  # random effects at whatever start the model was built with
  if(is.null(par_vec)) {
    fixed_start <- obj$par
    rand_start <- if(is.null(rand_idx)) numeric(0) else obj$env$par[rand_idx]
  } else if(length(par_vec) == n_fixed) {
    fixed_start <- par_vec
    rand_start <- if(is.null(rand_idx)) numeric(0) else obj$env$par[rand_idx]
  } else if(!is.null(rand_idx) && length(par_vec) == n_joint) {
    fixed_start <- par_vec[-rand_idx]
    rand_start <- par_vec[rand_idx]
  } else {
    stop("par_vec has length ", length(par_vec), ", but this model has ", n_fixed,
         " fixed effects",
         if(is.null(rand_idx)) "" else paste0(" and ", n_joint, " joint fixed and random values"),
         ". Pass a vector of one of those lengths.")
  }

  fixed <- fixed_start + stats::rnorm(length(fixed_start), 0, sd)

  # the random effects are integrated out, so these only move the starting point of
  # the inner solve; they are returned jittered or not so the caller always seeds them
  random <- if(is.null(rand_idx)) numeric(0)
            else if(isTRUE(jitter_random)) rand_start + stats::rnorm(length(rand_start), 0, sd)
            else rand_start

  list(fixed = fixed, random = random)
}

#' Run Jitter Analysis for Model Diagnostics
#'
#' Tests how sensitive the optimization is to its starting values by perturbing
#' the parameter vector with additive normal noise, refitting, and recording the
#' resulting series and diagnostics. Each iteration perturbs the fixed effects,
#' and the random effects too under \code{jitter_random = TRUE}, optimizes with
#' \code{stats::nlminb()}, optionally takes extra Newton steps, and extracts the
#' reported quantities. Iterations run sequentially or in parallel through
#' \code{future}.
#'
#' @param data List of model data for the \code{RTMB} objective function.
#' @param parameters Named list of parameters for \code{RTMB::MakeADFun()}.
#' @param mapping Named list of parameter mappings.
#' @param random Character vector of random-effect parameters.
#' @param sd Standard deviation of the additive normal noise.
#' @param n_jitter Number of jittered runs.
#' @param n_newton_loops Extra Newton steps after \code{nlminb()} converges.
#'   Default 0.
#' @param do_par Logical, whether the iterations run in parallel.
#' @param n_cores Parallel workers used when \code{do_par = TRUE}.
#' @param par_vec Optional starting values to jitter, either the fixed-effect
#'   vector (\code{length(obj$par)}, such as \code{fit$optim$par}) or the joint
#'   fixed and random vector (\code{length(obj$env$par)}, such as
#'   \code{fit$env$last.par.best}). Any other length is an error, and \code{NULL}
#'   uses the model's own start.
#' @param jitter_random Logical, whether the random effects are perturbed
#'   alongside the fixed effects. Only the fixed effects are searched by
#'   \code{nlminb()}, so the random draws move the starting point of the inner
#'   Laplace solve and check whether it settles on the same modes. Either way the
#'   inner solve starts from the random values in \code{par_vec}, or from the
#'   model's own start when \code{par_vec} holds none. Default \code{FALSE}.
#'
#' @return A data frame of the jitter results: the spawning stock biomass and
#'   recruitment series of each run, with its jitter index, whether the Hessian
#'   was positive definite, the joint negative log-likelihood, and the maximum
#'   absolute gradient of the fixed effects.
#'
#' @import RTMB
#' @import future
#' @import future.apply
#' @import progressr
#' @import dplyr
#' @importFrom reshape2 melt
#' @importFrom stats rnorm nlminb optimHess
#'
#' @family Model Diagnostics
#' @export do_jitter
do_jitter <- function(data,
                      parameters,
                      mapping,
                      random = NULL,
                      sd,
                      n_jitter,
                      n_newton_loops = 0,
                      do_par,
                      n_cores,
                      par_vec = NULL,
                      jitter_random = FALSE
                      ) {

  jitter_all <- data.frame()

  obj <- RTMB::MakeADFun(
    cmb(SPoRC_rtmb, data),
    parameters = parameters,
    map = mapping,
    random = random,
    silent = TRUE
  )

  if(do_par == FALSE) {

    par_start <- obj$env$par # the model's own start, before any iteration moves the tape

    for(i in 1:n_jitter) {

      # one tape serves every iteration, so its state goes back to the model start each
      # time; value.best only ever falls, and a stale one repeats the earlier fit's report
      obj$env$last.par <- par_start
      obj$env$last.par.best <- par_start
      obj$env$value.best <- Inf

      # jitter starting values (additive normal draws)
      jit <- jitter_start_values(obj = obj, par_vec = par_vec, sd = sd, jitter_random = jitter_random)

      # the inner Laplace solve starts from last.par, so the random start goes there,
      # whether it came from par_vec, the model's own start, or a draw around either
      if(length(jit$random) > 0) {
        obj$env$last.par[obj$env$random] <- jit$random
        obj$env$last.par.best[obj$env$random] <- jit$random
      }

      # Now, optimize the function
      optim <- stats::nlminb(jit$fixed,
                             obj$fn,
                             obj$gr,
                             control = list(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-15))

      # newton steps
      tryCatch(
        expr =
                 for(j in 1:n_newton_loops) {
                   g = as.numeric(obj$gr(optim$par))
                   h = optimHess(optim$par, fn = obj$fn, gr = obj$gr)
                   optim$par = optim$par - solve(h,g)
                   optim$objective = obj$fn(optim$par)
                 },
        error = function(e) {e},
        warning = function(w) {w}
      )

      obj$rep <- obj$report(obj$env$last.par.best) # Get report
      obj$sd_rep <- RTMB::sdreport(obj) # Get sd report

      # put jitter results into a dataframe
      jitter_ts_df <- reshape2::melt(obj$rep$SSB) %>%
        dplyr::rename(Pop = Var1, Region = Var2, Year = Var3) %>%
        dplyr::mutate(Type = 'SSB') %>%
        dplyr::bind_rows(reshape2::melt(obj$rep$Rec) %>%
                           dplyr::rename(Pop = Var1, Region = Var2, Year = Var3) %>%
                           dplyr::mutate(Type = 'Recruitment')) %>%
        dplyr::mutate(jitter = i,
                      Hessian = obj$sd_rep$pdHess,
                      jnLL = obj$rep$jnLL,
                      Max_Gradient = max(abs(obj$sd_rep$gradient.fixed)))

      jitter_all <- rbind(jitter_all, jitter_ts_df) # bind dataframes

    } # end i loop
  } # don't parrallelize

  if(do_par == TRUE) {

    future::plan(future::multisession, workers = n_cores) # set up cores

    progressr::with_progress({

      p <- progressr::progressor(along = 1:n_jitter) # progress bar

      jitter_all <- future.apply::future_lapply(1:n_jitter, function(i) {

        # make obj
        obj <- RTMB::MakeADFun(
          cmb(SPoRC_rtmb, data),
          parameters = parameters,
          map = mapping,
          random = random,
          silent = TRUE
        )

        # Jitter starting values
        jit <- jitter_start_values(obj = obj, par_vec = par_vec, sd = sd, jitter_random = jitter_random)

        # the inner Laplace solve starts from last.par, so the random draws go there
        if(length(jit$random) > 0) {
          obj$env$last.par[obj$env$random] <- jit$random
          obj$env$last.par.best[obj$env$random] <- jit$random
        }

        # Optimize function
        optim <- stats::nlminb(jit$fixed,
                               obj$fn,
                               obj$gr,
                               control = list(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-15))

        # Newton steps
        tryCatch({
          for (j in 1:n_newton_loops) {
            g <- as.numeric(obj$gr(optim$par))
            h <- optimHess(optim$par, fn = obj$fn, gr = obj$gr)
            optim$par <- optim$par - solve(h, g)
            optim$objective <- obj$fn(optim$par)
          }
        }, error = function(e) e, warning = function(w) w)

        # get reports
        obj$rep <- obj$report(obj$env$last.par.best)
        obj$sd_rep <- RTMB::sdreport(obj)

        # put jitter results into a dataframe
        jitter_ts_df <- reshape2::melt(obj$rep$SSB) %>%
          dplyr::rename(Pop = Var1, Region = Var2, Year = Var3) %>%
          dplyr::mutate(Type = 'SSB') %>%
          dplyr::bind_rows(reshape2::melt(obj$rep$Rec) %>%
                             dplyr::rename(Pop = Var1, Region = Var2, Year = Var3) %>%
                             dplyr::mutate(Type = 'Recruitment')) %>%
          dplyr::mutate(jitter = i,
                        Hessian = obj$sd_rep$pdHess,
                        jnLL = obj$rep$jnLL,
                        Max_Gradient = max(abs(obj$sd_rep$gradient.fixed)))

        p() # update progress

        jitter_ts_df

      }, future.seed = TRUE) %>%
        bind_rows() # bine rows to combine results

      future::plan(future::sequential)  # Reset

    })

  } # end dor par

  return(jitter_all)
}
