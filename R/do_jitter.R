#' Run Jitter Analysis for Model Diagnostics
#'
#' Performs a jitter analysis to evaluate sensitivity of model optimization
#' to starting parameter values. The function repeatedly perturbs the
#' parameter vector with additive normal noise, refits the model, and records
#' resulting time series and diagnostic metrics.
#'
#' Each jitter iteration:
#' \itemize{
#'   \item Perturbs the starting parameter vector with random normal noise.
#'   \item Optimizes the objective function using \code{stats::nlminb()}.
#'   \item Optionally performs additional Newton steps to refine the solution.
#'   \item Extracts reported quantities (e.g., spawning biomass and recruitment)
#'   and diagnostic statistics.
#' }
#'
#' The analysis can be executed sequentially or in parallel using the
#' \code{future} framework.
#'
#' @param data A list of model data used to construct the \code{RTMB} objective
#'   function.
#' @param parameters A named list of model parameters used to initialize
#'   \code{RTMB::MakeADFun()}.
#' @param mapping A named list defining parameter mappings for
#'   \code{RTMB::MakeADFun()}.
#' @param random Character vector specifying random-effect parameters.
#' @param sd Numeric value specifying the standard deviation of the additive
#'   normal noise used to jitter parameters.
#' @param n_jitter Integer specifying the number of jittered optimization runs.
#' @param n_newton_loops Integer specifying the number of additional Newton
#'   optimization steps performed after \code{nlminb()} convergence.
#' @param do_par Logical indicating whether jitter iterations should be
#'   executed in parallel.
#' @param n_cores Integer specifying the number of parallel workers to use
#'   when \code{do_par = TRUE}.
#' @param par_vec Optional numeric vector of parameter values used as the
#'   starting point for jittering. If \code{NULL}, the model's default starting
#'   parameter vector is jittered. If provided, jittering is applied to this
#'   vector (for example, the maximum likelihood estimates).
#'
#' @return A \code{data.frame} containing jitter iteration results. The output
#' includes time series of spawning stock biomass (SSB) and recruitment,
#' along with diagnostic information for each jitter run, including:
#' \itemize{
#'   \item Jitter index
#'   \item Whether the Hessian is positive definite
#'   \item Joint negative log-likelihood
#'   \item Maximum absolute gradient of fixed effects
#' }
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
                      n_newton_loops,
                      do_par,
                      n_cores,
                      par_vec = NULL
                      ) {

  jitter_all <- data.frame()

  obj <- RTMB::MakeADFun(cmb(SPoRC_rtmb, data),
                         parameters = parameters,
                         map = mapping,
                         random = random,
                         silent = TRUE)

  if(do_par == FALSE) {
    for(i in 1:n_jitter) {

      # jitter original parameters (additive normal draws)
      if(is.null(par_vec)) jitter_pars <- obj$par + stats::rnorm(length(obj$par), 0, sd) # if not using mle parameters
      else jitter_pars <- par_vec + stats::rnorm(length(obj$par), 0, sd) # if using mle parameters

      # Now, optimize the function
      optim <- stats::nlminb(jitter_pars,
                             obj$fn,
                             obj$gr,
                             control = list(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-15))

      # newton steps
      try_improve <- tryCatch(expr =
                                for(j in 1:n_newton_loops) {
                                  g = as.numeric(obj$gr(optim$par))
                                  h = optimHess(optim$par, fn = obj$fn, gr = obj$gr)
                                  optim$par = optim$par - solve(h,g)
                                  optim$objective = obj$fn(optim$par)
                                }
                              , error = function(e){e}, warning = function(w){w})

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
        obj <- RTMB::MakeADFun(cmb(SPoRC_rtmb, data), parameters = parameters,  map = mapping, random = random, silent = TRUE)

        # Jitter original parameters
        if(is.null(par_vec)) jitter_pars <- obj$par + stats::rnorm(length(obj$par), 0, sd) # if not using mle parameters
        else jitter_pars <- par_vec + stats::rnorm(length(obj$par), 0, sd) # if using mle parameters

        # Optimize function
        optim <- stats::nlminb(jitter_pars,
                               obj$fn,
                               obj$gr,
                               control = list(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-15))

        # Newton steps
        try_improve <- tryCatch({
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

      }, future.seed = TRUE) %>% bind_rows() # bine rows to combine results

      future::plan(future::sequential)  # Reset

    })

  } # end dor par

  return(jitter_all)
}
