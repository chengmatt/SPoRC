#' Run Likelihood Profile
#'
#' Profiles the joint negative log-likelihood and all individual likelihood
#' components across a range of fixed values for a single parameter. Supports
#' both sequential and parallel execution.
#'
#' @param data Data list from the fitted model.
#' @param parameters Parameter list from the fitted model.
#' @param mapping Mapping list from the fitted model.
#' @param random Character vector of random effects to estimate. Default
#'   \code{NULL}.
#' @param what Character string. Name of the parameter to profile.
#' @param idx Vector pointing to the
#'   specific elements to fix when \code{parameters[[what]]} is an array.
#'   \code{NULL} for scalar parameters.
#' @param min_val Numeric. Minimum value of the profile range.
#' @param max_val Numeric. Maximum value of the profile range.
#' @param inc Numeric. Increment between profile values. Default \code{0.05}.
#' @param do_par Logical. Whether to use parallel processing. Default
#'   \code{FALSE}.
#' @param n_cores Integer. Number of parallel workers. If \code{NULL}
#'   (default), \code{parallel::detectCores() - 1} is used.
#'
#' @return A named list containing one dataframe per likelihood component,
#'   each with a \code{prof_val} column indicating the profiled parameter
#'   value, plus \code{agg_nLL} which aggregates all components across their
#'   respective dimensions. Components include:
#'   \describe{
#'     \item{Scalar penalties and priors}{\code{jnLL_df}, \code{rec_nLL_df},
#'       \code{M_nLL_df}, \code{sel_nLL_df}, \code{rec_prop_nLL_df},
#'       \code{Movement_nLL_df}, \code{h_nLL_df}, \code{TagRep_nLL_df},
#'       \code{Fmort_nLL_df}, \code{fish_q_nLL_df}, \code{srv_q_nLL_df}.}
#'     \item{Pooled data likelihoods}{\code{Catch_nLL_df}
#'       \code{[Region × Year × Seas × Fleet]},
#'       \code{FishIdx_nLL_df} and \code{SrvIdx_nLL_df}
#'       \code{[Region × Year × Seas × Fleet]},
#'       \code{FishAge_nLL_df}, \code{FishLen_nLL_df},
#'       \code{SrvAge_nLL_df}, \code{SrvLen_nLL_df}
#'       \code{[Region × Year × Seas × Sex × Fleet]},
#'       \code{conv_fish_tag_nLL_df}
#'       \code{[Recap_Year × Recap_Seas × Tag_Cohort × Region × Fleet]}.}
#'     \item{Population-specific data likelihoods}{\code{Catch_pop_nLL_df},
#'       \code{FishIdx_pop_nLL_df}, \code{SrvIdx_pop_nLL_df}
#'       \code{[Pop × Region × Year × Seas × Fleet]},
#'       \code{FishAge_pop_nLL_df}, \code{FishLen_pop_nLL_df},
#'       \code{SrvAge_pop_nLL_df}, \code{SrvLen_pop_nLL_df}
#'       \code{[Pop × Region × Year × Seas × Sex × Fleet]}.}
#'   }
#'
#' @import dplyr
#' @import RTMB
#' @importFrom reshape2 melt
#' @importFrom stats rnorm nlminb
#' @importFrom future plan multisession
#' @importFrom future.apply future_lapply
#' @importFrom progressr with_progress progressor
#' @importFrom parallel detectCores
#' @export do_likelihood_profile
#' @family Model Diagnostics
do_likelihood_profile <- function(data,
                                  parameters,
                                  mapping,
                                  random = NULL,
                                  what,
                                  idx = NULL,
                                  min_val,
                                  max_val,
                                  inc = 0.05,
                                  do_par = FALSE,
                                  n_cores = NULL
) {

  if(min_val > max_val) {
    stop("`min_val` is greater than `max_val`. This likely occurred because you are profiling a log-transformed parameter. Try swapping the values: use the current `min_val` as `max_val`, and vice versa.")
  }

  # create values to profile across
  vals <- seq(min_val, max_val, inc)

  # Create objects to store values
  jnLL <- matrix(NA, nrow = length(vals), ncol = 1, dimnames = list(vals, NULL))
  rec_nLL <- matrix(NA, nrow = length(vals), ncol = 1, dimnames = list(vals, NULL))
  sel_nLL <- matrix(NA, nrow = length(vals), ncol = 1, dimnames = list(vals, NULL))
  M_nLL <- matrix(NA, nrow = length(vals), ncol = 1, dimnames = list(vals, NULL))
  rec_prop_nLL <- matrix(NA, nrow = length(vals), ncol = 1, dimnames = list(vals, NULL))
  h_nLL <- matrix(NA, nrow = length(vals), ncol = 1, dimnames = list(vals, NULL))
  Movement_nLL <- matrix(NA, nrow = length(vals), ncol = 1, dimnames = list(vals, NULL))
  TagRep_nLL <- matrix(NA, nrow = length(vals), ncol = 1, dimnames = list(vals, NULL))
  Fmort_nLL <- matrix(NA, nrow = length(vals), ncol = 1, dimnames = list(vals, NULL))
  dmr_nLL <- matrix(NA, nrow = length(vals), ncol = 1, dimnames = list(vals, NULL))
  conv_fish_tag_nLL <- data.frame()
  Catch_nLL <- data.frame()
  Discard_nLL <- data.frame()
  Discard_pop_nLL_df <- data.frame()
  FishAge_nLL <- data.frame()
  FishAgeComps_discard_nLL <- data.frame()
  SrvAge_nLL <- data.frame()
  SrvLen_nLL <- data.frame()
  FishLen_nLL <- data.frame()
  FishLenComps_discard_nLL <- data.frame()
  FishIdx_nLL <- data.frame()
  SrvIdx_nLL <- data.frame()
  Catch_pop_nLL <- data.frame()
  Discard_pop_nLL <- data.frame()
  FishAge_pop_nLL <- data.frame()
  FishLen_pop_nLL <- data.frame()
  FishAge_discard_pop_nLL <- data.frame()
  FishLen_discard_pop_nLL <- data.frame()
  SrvAge_pop_nLL <- data.frame()
  SrvLen_pop_nLL <- data.frame()
  FishIdx_pop_nLL <- data.frame()
  SrvIdx_pop_nLL <- data.frame()
  fish_q_nLL <- matrix(NA, nrow = length(vals), ncol = 1, dimnames = list(vals, NULL))
  srv_q_nLL  <- matrix(NA, nrow = length(vals), ncol = 1, dimnames = list(vals, NULL))

  # If there is more than one value in this parameter
  if(do_par == FALSE) {
    for(j in 1:length(vals)) {
      if(!is.null(dim(parameters[[what]]))) {
        # Input fixed values for all indices
        for(k in 1:length(idx)) {
          parameters[[what]] <- do.call(`[<-`, c(list(parameters[[what]]), idx[k], list(vals[j])))
        }
        # Build map with all target indices set to NA at once
        map_parameter <- parameters[[what]]
        for(k in 1:length(idx)) {
          map_parameter <- do.call(`[<-`, c(list(map_parameter), idx[k], list(NA)))
        }
        # Renumber non-NA positions
        counter <- 1
        non_na <- which(!is.na(map_parameter))
        for(i in 1:length(non_na)) {
          map_parameter[non_na[i]] <- counter
          counter <- counter + 1
        }
        mapping[[what]] <- factor(map_parameter)
      } else { # else, there is only one value in this parameter
        parameters[[what]] <- vals[j]
        mapping[[what]] <- factor(NA)
      }

      # make adfun
      SPoRC_rtmb_model <- RTMB::MakeADFun(cmb(SPoRC_rtmb, data), parameters = parameters, map = mapping, random = random, silent = TRUE)

      # Within loop
      tryCatch({
        SPoRC_optim <- stats::nlminb(SPoRC_rtmb_model$par, SPoRC_rtmb_model$fn, SPoRC_rtmb_model$gr,
                                     control = list(iter.max = 1e6, eval.max = 1e6, rel.tol = 1e-15))

        # Get report
        report <- SPoRC_rtmb_model$report(SPoRC_rtmb_model$env$last.par.best)

        # Store values and save (note, some need to save wt*nLL, because of how nlL are combined in the jnLL for the model)
        jnLL[j,1] <- report$jnLL
        rec_nLL[j,1] <- sum(data$Wt_Rec * report$Init_Rec_nLL, data$Wt_Rec * report$Rec_nLL)
        M_nLL[j,1] <- report$M_nLL
        sel_nLL[j,1] <- report$sel_nLL
        rec_prop_nLL[j,1] <- report$rec_prop_nLL
        Movement_nLL[j,1] <- report$Movement_nLL
        h_nLL[j,1] <- report$h_nLL
        TagRep_nLL[j,1] <- report$TagRep_nLL
        Fmort_nLL[j,1] <- sum(data$Wt_F * report$Fmort_nLL)
        dmr_nLL[j,1] <- sum(data$Wt_D * report$dmr_nLL)
        conv_fish_tag_nLL <- rbind(conv_fish_tag_nLL, reshape2::melt(data$Wt_Tagging * report$conv_fish_tag_nLL) %>% dplyr::mutate(prof_val = vals[j]))
        Catch_nLL <- rbind(Catch_nLL, reshape2::melt(data$Wt_Catch * report$Catch_nLL) %>% dplyr::mutate(prof_val = vals[j]))
        Discard_nLL <- rbind(Discard_nLL, reshape2::melt(data$Wt_Discard * report$Discard_nLL) %>% dplyr::mutate(prof_val = vals[j]))
        FishAge_nLL <- rbind(FishAge_nLL, reshape2::melt(report$FishAgeComps_nLL) %>% dplyr::mutate(prof_val = vals[j]))
        FishAgeComps_discard_nLL <- rbind(FishAgeComps_discard_nLL, reshape2::melt(report$FishAgeComps_discard_nLL) %>% dplyr::mutate(prof_val = vals[j]))
        SrvAge_nLL <- rbind(SrvAge_nLL, reshape2::melt(report$SrvAgeComps_nLL) %>% dplyr::mutate(prof_val = vals[j]))
        SrvLen_nLL <- rbind(SrvLen_nLL, reshape2::melt(report$SrvLenComps_nLL) %>% dplyr::mutate(prof_val = vals[j]))
        FishLen_nLL <- rbind(FishLen_nLL, reshape2::melt(report$FishLenComps_nLL) %>% dplyr::mutate(prof_val = vals[j]))
        FishLenComps_discard_nLL <- rbind(FishLenComps_discard_nLL, reshape2::melt(report$FishLenComps_discard_nLL) %>% dplyr::mutate(prof_val = vals[j]))
        FishIdx_nLL <- rbind(FishIdx_nLL, reshape2::melt(data$Wt_FishIdx * report$FishIdx_nLL) %>% dplyr::mutate(prof_val = vals[j]))
        SrvIdx_nLL <- rbind(SrvIdx_nLL, reshape2::melt(data$Wt_SrvIdx * report$SrvIdx_nLL) %>% dplyr::mutate(prof_val = vals[j]))
        fish_q_nLL[j,1] <- report$fish_q_nLL
        srv_q_nLL[j,1]  <- report$srv_q_nLL
        Catch_pop_nLL   <- rbind(Catch_pop_nLL,   reshape2::melt(data$Wt_Catch_pop   * report$Catch_pop_nLL) %>% dplyr::mutate(prof_val = vals[j]))
        Discard_pop_nLL   <- rbind(Discard_pop_nLL,   reshape2::melt(data$Wt_Discard_pop   * report$Discard_pop_nLL) %>% dplyr::mutate(prof_val = vals[j]))
        FishIdx_pop_nLL <- rbind(FishIdx_pop_nLL, reshape2::melt(data$Wt_FishIdx_pop * report$FishIdx_pop_nLL) %>% dplyr::mutate(prof_val = vals[j]))
        SrvIdx_pop_nLL  <- rbind(SrvIdx_pop_nLL,  reshape2::melt(data$Wt_SrvIdx_pop  * report$SrvIdx_pop_nLL)  %>% dplyr::mutate(prof_val = vals[j]))
        FishAge_pop_nLL <- rbind(FishAge_pop_nLL, reshape2::melt(report$FishAgeComps_pop_nLL) %>% dplyr::mutate(prof_val = vals[j]))
        FishLen_pop_nLL <- rbind(FishLen_pop_nLL, reshape2::melt(report$FishLenComps_pop_nLL) %>% dplyr::mutate(prof_val = vals[j]))
        FishAge_discard_pop_nLL <- rbind(FishAge_discard_pop_nLL, reshape2::melt(report$FishAgeComps_discard_pop_nLL) %>% dplyr::mutate(prof_val = vals[j]))
        FishLen_discard_pop_nLL <- rbind(FishLen_discard_pop_nLL, reshape2::melt(report$FishLenComps_discard_pop_nLL) %>% dplyr::mutate(prof_val = vals[j]))
        SrvAge_pop_nLL  <- rbind(SrvAge_pop_nLL,  reshape2::melt(report$SrvAgeComps_pop_nLL) %>% dplyr::mutate(prof_val = vals[j]))
        SrvLen_pop_nLL  <- rbind(SrvLen_pop_nLL,  reshape2::melt(report$SrvLenComps_pop_nLL) %>% dplyr::mutate(prof_val = vals[j]))

        print(paste("Likelihood profile is at:", round(j / length(vals) * 100, 2), "%"))

      }, error = function(e) {
        message("Failed to optimize: ", e$message)
      })

    } # end j loop
  }

  if(do_par == TRUE) {

    # Set up parallel processing
    future::plan(future::multisession, workers = n_cores)

    progressr::with_progress({

      p <- progressr::progressor(along = 1:length(vals))

      # Run parallel processing
      profile_results <- future.apply::future_lapply(1:length(vals), function(j) {

        # Create local copies to avoid conflicts
        local_data <- data
        local_parameters <- parameters
        local_mapping <- mapping

        # Initialize result list
        result <- list(
          j = j,
          prof_val = vals[j],
          success = FALSE,
          jnLL = NA,
          rec_nLL = NA,
          M_nLL = NA,
          sel_nLL = NA,
          rec_prop_nLL = NA,
          Movement_nLL = NA,
          h_nLL = NA,
          TagRep_nLL = NA,
          Fmort_nLL = NA,
          dmr_nLL = NA,
          conv_fish_tag_nLL = data.frame(),
          Catch_nLL = data.frame(),
          Discard_nLL = data.frame(),
          FishAge_nLL = data.frame(),
          FishAgeComps_discard_nLL = data.frame(),
          SrvAge_nLL = data.frame(),
          SrvLen_nLL = data.frame(),
          FishLen_nLL = data.frame(),
          FishLenComps_discard_nLL = data.frame(),
          FishIdx_nLL = data.frame(),
          SrvIdx_nLL = data.frame(),
          fish_q_nLL = NA,
          srv_q_nLL  = NA,
          Catch_pop_nLL   = data.frame(),
          Discard_pop_nLL = data.frame(),
          FishIdx_pop_nLL = data.frame(),
          SrvIdx_pop_nLL  = data.frame(),
          FishAge_pop_nLL = data.frame(),
          FishLen_pop_nLL = data.frame(),
          FishAge_discard_pop_nLL = data.frame(),
          FishLen_discard_pop_nLL = data.frame(),
          SrvAge_pop_nLL  = data.frame(),
          SrvLen_pop_nLL  = data.frame()
        )

        if(!is.null(dim(parameters[[what]]))) {
          # Input fixed values for all indices
          for(k in 1:length(idx)) {
            parameters[[what]] <- do.call(`[<-`, c(list(parameters[[what]]), idx[k], list(vals[j])))
          }
          # Build map with all target indices set to NA at once
          map_parameter <- parameters[[what]]
          for(k in 1:length(idx)) {
            map_parameter <- do.call(`[<-`, c(list(map_parameter), idx[k], list(NA)))
          }
          # Renumber non-NA positions
          counter <- 1
          non_na <- which(!is.na(map_parameter))
          for(i in 1:length(non_na)) {
            map_parameter[non_na[i]] <- counter
            counter <- counter + 1
          }
          mapping[[what]] <- factor(map_parameter)
        } else { # else, there is only one value in this parameter
          local_parameters[[what]] <- vals[j]
          local_mapping[[what]] <- factor(NA)
        }

        # make adfun
        tryCatch({
          SPoRC_rtmb_model <- RTMB::MakeADFun(cmb(SPoRC_rtmb, local_data), parameters = local_parameters, map = local_mapping, random = random, silent = TRUE)

          SPoRC_optim <- stats::nlminb(SPoRC_rtmb_model$par, SPoRC_rtmb_model$fn, SPoRC_rtmb_model$gr,
                                       control = list(iter.max = 1e6, eval.max = 1e6, rel.tol = 1e-15))

          # Get report
          report <- SPoRC_rtmb_model$report(SPoRC_rtmb_model$env$last.par.best)

          # Store values and save
          result$jnLL <- report$jnLL
          result$rec_nLL <- sum(data$Wt_Rec * report$Init_Rec_nLL, data$Wt_Rec * report$Rec_nLL)
          result$M_nLL <- report$M_nLL
          result$sel_nLL <- report$sel_nLL
          result$rec_prop_nLL <- report$rec_prop_nLL
          result$Movement_nLL <- report$Movement_nLL
          result$h_nLL <- report$h_nLL
          result$TagRep_nLL <- report$TagRep_nLL
          result$Fmort_nLL <- sum(data$Wt_F * report$Fmort_nLL)
          result$dmr_nLL <- sum(data$Wt_D * report$dmr_nLL)
          result$conv_fish_tag_nLL <- reshape2::melt(data$Wt_Tagging * report$conv_fish_tag_nLL) %>% dplyr::mutate(prof_val = vals[j])
          result$Catch_nLL <- reshape2::melt(data$Wt_Catch * report$Catch_nLL) %>% dplyr::mutate(prof_val = vals[j])
          result$Discard_nLL <- reshape2::melt(data$Wt_Discard * report$Discard_nLL) %>% dplyr::mutate(prof_val = vals[j])
          result$FishAge_nLL <- reshape2::melt(report$FishAgeComps_nLL) %>% dplyr::mutate(prof_val = vals[j])
          result$FishAgeComps_discard_nLL <- reshape2::melt(report$FishAgeComps_discard_nLL) %>% dplyr::mutate(prof_val = vals[j])
          result$SrvAge_nLL <- reshape2::melt(report$SrvAgeComps_nLL) %>% dplyr::mutate(prof_val = vals[j])
          result$SrvLen_nLL <- reshape2::melt(report$SrvLenComps_nLL) %>% dplyr::mutate(prof_val = vals[j])
          result$FishLen_nLL <- reshape2::melt(report$FishLenComps_nLL) %>% dplyr::mutate(prof_val = vals[j])
          result$FishLenComps_discard_nLL <- reshape2::melt(report$FishLenComps_discard_nLL) %>% dplyr::mutate(prof_val = vals[j])
          result$FishIdx_nLL <- reshape2::melt(data$Wt_FishIdx * report$FishIdx_nLL) %>% dplyr::mutate(prof_val = vals[j])
          result$SrvIdx_nLL <- reshape2::melt(data$Wt_SrvIdx * report$SrvIdx_nLL) %>% dplyr::mutate(prof_val = vals[j])
          result$fish_q_nLL    <- report$fish_q_nLL
          result$srv_q_nLL     <- report$srv_q_nLL
          result$Catch_pop_nLL   <- reshape2::melt(data$Wt_Catch_pop   * report$Catch_pop_nLL)   %>% dplyr::mutate(prof_val = vals[j])
          result$Discard_pop_nLL   <- reshape2::melt(data$Wt_Discard_pop   * report$Discard_pop_nLL)   %>% dplyr::mutate(prof_val = vals[j])
          result$FishIdx_pop_nLL <- reshape2::melt(data$Wt_FishIdx_pop * report$FishIdx_pop_nLL) %>% dplyr::mutate(prof_val = vals[j])
          result$SrvIdx_pop_nLL  <- reshape2::melt(data$Wt_SrvIdx_pop  * report$SrvIdx_pop_nLL)  %>% dplyr::mutate(prof_val = vals[j])
          result$FishAge_pop_nLL <- reshape2::melt(report$FishAgeComps_pop_nLL) %>% dplyr::mutate(prof_val = vals[j])
          result$FishLen_pop_nLL <- reshape2::melt(report$FishLenComps_pop_nLL) %>% dplyr::mutate(prof_val = vals[j])
          result$FishAge_discard_pop_nLL <- reshape2::melt(report$FishAgeComps_discard_pop_nLL) %>% dplyr::mutate(prof_val = vals[j])
          result$FishLen_discard_pop_nLL <- reshape2::melt(report$FishLenComps_discard_pop_nLL) %>% dplyr::mutate(prof_val = vals[j])
          result$SrvAge_pop_nLL  <- reshape2::melt(report$SrvAgeComps_pop_nLL) %>% dplyr::mutate(prof_val = vals[j])
          result$SrvLen_pop_nLL  <- reshape2::melt(report$SrvLenComps_pop_nLL) %>% dplyr::mutate(prof_val = vals[j])
          result$success <- TRUE

        }, error = function(e) {
          result$error_msg <- e$message
        })

        # Update progress
        p(sprintf("j=%g", j))

        return(result)

      }, future.seed = TRUE)

    })

    # Process parallel results back into the original matrices/data.frames
    for(res in profile_results) {
      if(res$success) {
        j <- res$j
        jnLL[j,1] <- res$jnLL
        rec_nLL[j,1] <- res$rec_nLL
        M_nLL[j,1] <- res$M_nLL
        sel_nLL[j,1] <- res$sel_nLL
        rec_prop_nLL[j,1] <- res$rec_prop_nLL
        Movement_nLL[j,1] <- res$Movement_nLL
        h_nLL[j,1] <- res$h_nLL
        TagRep_nLL[j,1] <- res$TagRep_nLL
        Fmort_nLL[j,1] <- res$Fmort_nLL
        dmr_nLL[j,1] <- res$dmr_nLL
        conv_fish_tag_nLL <- rbind(conv_fish_tag_nLL, res$conv_fish_tag_nLL)
        Catch_nLL <- rbind(Catch_nLL, res$Catch_nLL)
        Discard_nLL <- rbind(Discard_nLL, res$Discard_nLL)
        FishAge_nLL <- rbind(FishAge_nLL, res$FishAge_nLL)
        FishAgeComps_discard_nLL <- rbind(FishAgeComps_discard_nLL, res$FishAgeComps_discard_nLL)
        SrvAge_nLL <- rbind(SrvAge_nLL, res$SrvAge_nLL)
        SrvLen_nLL <- rbind(SrvLen_nLL, res$SrvLen_nLL)
        FishLen_nLL <- rbind(FishLen_nLL, res$FishLen_nLL)
        FishLenComps_discard_nLL <- rbind(FishLenComps_discard_nLL, res$FishLenComps_discard_nLL)
        FishIdx_nLL <- rbind(FishIdx_nLL, res$FishIdx_nLL)
        SrvIdx_nLL <- rbind(SrvIdx_nLL, res$SrvIdx_nLL)
        fish_q_nLL[j,1]  <- res$fish_q_nLL
        srv_q_nLL[j,1]   <- res$srv_q_nLL
        Catch_pop_nLL   <- rbind(Catch_pop_nLL,   res$Catch_pop_nLL)
        Discard_pop_nLL   <- rbind(Discard_pop_nLL,   res$Discard_pop_nLL)
        FishIdx_pop_nLL <- rbind(FishIdx_pop_nLL, res$FishIdx_pop_nLL)
        SrvIdx_pop_nLL  <- rbind(SrvIdx_pop_nLL,  res$SrvIdx_pop_nLL)
        FishAge_pop_nLL <- rbind(FishAge_pop_nLL, res$FishAge_pop_nLL)
        FishLen_pop_nLL <- rbind(FishLen_pop_nLL, res$FishLen_pop_nLL)
        FishAge_discard_pop_nLL <- rbind(FishAge_discard_pop_nLL, res$FishAge_discard_pop_nLL)
        FishLen_discard_pop_nLL <- rbind(FishLen_discard_pop_nLL, res$FishLen_discard_pop_nLL)
        SrvAge_pop_nLL  <- rbind(SrvAge_pop_nLL,  res$SrvAge_pop_nLL)
        SrvLen_pop_nLL  <- rbind(SrvLen_pop_nLL,  res$SrvLen_pop_nLL)

      } else {
        message("Failed iteration ", res$j, " (value = ", res$prof_val, "): ", res$error_msg)
      }
    }

    # Close parallel workers
    future::plan(future::sequential)
    message("Parallel processing completed.")
  }


  # Doing some residual munging into the correct format
  jnLL_df <- reshape2::melt(jnLL) %>%
    dplyr::select(-Var2) %>%
    dplyr::rename(prof_val = Var1) %>%
    dplyr::mutate(type = 'jnLL')
  rec_nLL_df <- reshape2::melt(rec_nLL) %>%
    dplyr::select(-Var2) %>%
    dplyr::rename(prof_val = Var1) %>%
    dplyr::mutate(type = 'RecPen')
  M_nLL_df <- reshape2::melt(M_nLL) %>%
    dplyr::select(-Var2) %>%
    dplyr::rename(prof_val = Var1) %>%
    dplyr::mutate(type = 'M Prior')
  rec_prop_nLL_df <- reshape2::melt(rec_prop_nLL) %>%
    dplyr::select(-Var2) %>%
    dplyr::rename(prof_val = Var1) %>%
    dplyr::mutate(type = 'Recruitment Prop Prior')
  sel_nLL_df <- reshape2::melt(sel_nLL) %>%
    dplyr::select(-Var2) %>%
    dplyr::rename(prof_val = Var1) %>%
    dplyr::mutate(type = 'Selex Pen')
  Movement_nLL_df <- reshape2::melt(Movement_nLL) %>%
    dplyr::select(-Var2) %>%
    dplyr::rename(prof_val = Var1) %>%
    dplyr::mutate(type = 'Move Prior')
  h_nLL_df <- reshape2::melt(h_nLL) %>% dplyr::select(-Var2) %>% dplyr::rename(prof_val = Var1) %>% dplyr::mutate(type = 'h Prior')
  TagRep_nLL_df <- reshape2::melt(TagRep_nLL) %>%
    dplyr::select(-Var2) %>%
    dplyr::rename(prof_val = Var1) %>%
    dplyr::mutate(type = 'TagRep Prior')
  Fmort_nLL_df <- reshape2::melt(Fmort_nLL) %>%
    dplyr::select(-Var2) %>%
    dplyr::rename(prof_val = Var1) %>%
    dplyr::mutate(type = 'FmortPen')
  dmr_nLL_df <- reshape2::melt(dmr_nLL) %>%
    dplyr::select(-Var2) %>%
    dplyr::rename(prof_val = Var1) %>%
    dplyr::mutate(type = 'dmrPen')
  Catch_nLL_df <- Catch_nLL %>%
    dplyr::rename(Region = Var1, Year = Var2, Seas = Var3, Fleet = Var4) %>%
    dplyr::mutate(type = 'Catch')
  Discard_nLL_df <- Discard_nLL %>%
    dplyr::rename(Region = Var1, Year = Var2, Seas = Var3, Fleet = Var4) %>%
    dplyr::mutate(type = 'Discard')
  FishAge_nLL_df <- FishAge_nLL %>%
    dplyr::rename(Region = Var1, Year = Var2, Seas = Var3, Sex = Var4, Fleet = Var5) %>%
    dplyr::mutate(type = 'FishAge')
  FishAgeComps_discard_nLL_df <- FishAgeComps_discard_nLL %>%
    dplyr::rename(Region = Var1, Year = Var2, Seas = Var3, Sex = Var4, Fleet = Var5) %>%
    dplyr::mutate(type = 'FishAgeDiscard')
  SrvAge_nLL_df <- SrvAge_nLL %>%
    dplyr::rename(Region = Var1, Year = Var2, Seas = Var3, Sex = Var4, Fleet = Var5) %>%
    dplyr::mutate(type = 'SrvAge')
  FishLen_nLL_df <- FishLen_nLL %>%
    dplyr::rename(Region = Var1, Year = Var2, Seas = Var3, Sex = Var4, Fleet = Var5) %>%
    dplyr::mutate(type = 'FishLen')
  FishLenComps_discard_nLL_df <- FishLenComps_discard_nLL %>%
    dplyr::rename(Region = Var1, Year = Var2, Seas = Var3, Sex = Var4, Fleet = Var5) %>%
    dplyr::mutate(type = 'FishLenDiscard')
  SrvLen_nLL_df <- SrvLen_nLL %>%
    dplyr::rename(Region = Var1, Year = Var2, Seas = Var3, Sex = Var4, Fleet = Var5) %>%
    dplyr::mutate(type = 'SrvLen')
  FishIdx_nLL_df <- FishIdx_nLL %>%
    dplyr::rename(Region = Var1, Year = Var2, Seas = Var3, Fleet = Var4) %>%
    dplyr::mutate(type = 'FishIdx')
  SrvIdx_nLL_df <- SrvIdx_nLL %>%
    dplyr::rename(Region = Var1, Year = Var2, Seas = Var3, Fleet = Var4) %>%
    dplyr::mutate(type = 'SrvIdx')
  conv_fish_tag_nLL_df <- conv_fish_tag_nLL %>%
    dplyr::rename(Recap_Year = Var1, Recap_Seas = Var2, Tag_Cohort = Var3, Region = Var4, Fleet = Var5) %>%
    dplyr::mutate(type = 'Tagging')
  fish_q_nLL_df <- reshape2::melt(fish_q_nLL) %>%
    dplyr::select(-Var2) %>% dplyr::rename(prof_val = Var1) %>% dplyr::mutate(type = 'FishQ Prior')
  srv_q_nLL_df <- reshape2::melt(srv_q_nLL) %>%
    dplyr::select(-Var2) %>% dplyr::rename(prof_val = Var1) %>% dplyr::mutate(type = 'SrvQ Prior')
  Catch_pop_nLL_df <- Catch_pop_nLL %>%
    dplyr::rename(Pop = Var1, Region = Var2, Year = Var3, Seas = Var4, Fleet = Var5) %>%
    dplyr::mutate(type = 'CatchPop')
  Discard_pop_nLL_df <- Discard_pop_nLL %>%
    dplyr::rename(Pop = Var1, Region = Var2, Year = Var3, Seas = Var4, Fleet = Var5) %>%
    dplyr::mutate(type = 'DiscardPop')
  FishIdx_pop_nLL_df <- FishIdx_pop_nLL %>%
    dplyr::rename(Pop = Var1, Region = Var2, Year = Var3, Seas = Var4, Fleet = Var5) %>%
    dplyr::mutate(type = 'FishIdxPop')
  SrvIdx_pop_nLL_df <- SrvIdx_pop_nLL %>%
    dplyr::rename(Pop = Var1, Region = Var2, Year = Var3, Seas = Var4, Fleet = Var5) %>%
    dplyr::mutate(type = 'SrvIdxPop')
  FishAge_pop_nLL_df <- FishAge_pop_nLL %>%
    dplyr::rename(Pop = Var1, Region = Var2, Year = Var3, Seas = Var4, Sex = Var5, Fleet = Var6) %>%
    dplyr::mutate(type = 'FishAgePop')
  FishAge_discard_pop_nLL_df <- FishAge_discard_pop_nLL %>%
    dplyr::rename(Pop = Var1, Region = Var2, Year = Var3, Seas = Var4, Sex = Var5, Fleet = Var6) %>%
    dplyr::mutate(type = 'FishAgeDiscardPop')
  FishLen_pop_nLL_df <- FishLen_pop_nLL %>%
    dplyr::rename(Pop = Var1, Region = Var2, Year = Var3, Seas = Var4, Sex = Var5, Fleet = Var6) %>%
    dplyr::mutate(type = 'FishLenPop')
  FishLen_discard_pop_nLL_df <- FishLen_discard_pop_nLL %>%
    dplyr::rename(Pop = Var1, Region = Var2, Year = Var3, Seas = Var4, Sex = Var5, Fleet = Var6) %>%
    dplyr::mutate(type = 'FishLenDiscardPop')
  SrvAge_pop_nLL_df <- SrvAge_pop_nLL %>%
    dplyr::rename(Pop = Var1, Region = Var2, Year = Var3, Seas = Var4, Sex = Var5, Fleet = Var6) %>%
    dplyr::mutate(type = 'SrvAgePop')
  SrvLen_pop_nLL_df <- SrvLen_pop_nLL %>%
    dplyr::rename(Pop = Var1, Region = Var2, Year = Var3, Seas = Var4, Sex = Var5, Fleet = Var6) %>%
    dplyr::mutate(type = 'SrvLenPop')

  # Get likelihoods aggregated across all dimensions
  agg_nLL <- rbind(jnLL_df, rec_nLL_df, M_nLL_df, rec_prop_nLL_df, Movement_nLL_df, h_nLL_df,
                   TagRep_nLL_df,Fmort_nLL_df, dmr_nLL_df, sel_nLL_df,
                   Catch_nLL_df %>% dplyr::group_by(prof_val, type) %>%
                     dplyr::summarize(value = sum(value)),
                   Discard_nLL_df %>% dplyr::group_by(prof_val, type) %>%
                     dplyr::summarize(value = sum(value)),
                   conv_fish_tag_nLL_df %>% dplyr::group_by(prof_val, type) %>%
                     dplyr::summarize(value = sum(value, na.rm = T)),
                   FishAge_nLL_df %>% dplyr::group_by(prof_val, type) %>%
                     dplyr::summarize(value = sum(value, na.rm = T)),
                   FishAgeComps_discard_nLL_df %>% dplyr::group_by(prof_val, type) %>%
                     dplyr::summarize(value = sum(value, na.rm = T)),
                   SrvAge_nLL_df %>% dplyr::group_by(prof_val, type) %>%
                     dplyr::summarize(value = sum(value, na.rm = T)),
                   FishLen_nLL_df %>% dplyr::group_by(prof_val, type) %>%
                     dplyr::summarize(value = sum(value, na.rm = T)),
                   FishLenComps_discard_nLL_df %>% dplyr::group_by(prof_val, type) %>%
                     dplyr::summarize(value = sum(value, na.rm = T)),
                   SrvLen_nLL_df %>% dplyr::group_by(prof_val, type) %>%
                     dplyr::summarize(value = sum(value, na.rm = T)),
                   FishIdx_nLL_df %>% dplyr::group_by(prof_val, type) %>%
                     dplyr::summarize(value = sum(value, na.rm = T)),
                   SrvIdx_nLL_df %>% dplyr::group_by(prof_val, type) %>%
                     dplyr::summarize(value = sum(value, na.rm = T)),
                   fish_q_nLL_df,
                   srv_q_nLL_df,
                   Catch_pop_nLL_df   %>% dplyr::group_by(prof_val, type) %>% dplyr::summarize(value = sum(value, na.rm = T)),
                   Discard_pop_nLL_df   %>% dplyr::group_by(prof_val, type) %>% dplyr::summarize(value = sum(value, na.rm = T)),
                   FishIdx_pop_nLL_df %>% dplyr::group_by(prof_val, type) %>% dplyr::summarize(value = sum(value, na.rm = T)),
                   SrvIdx_pop_nLL_df  %>% dplyr::group_by(prof_val, type) %>% dplyr::summarize(value = sum(value, na.rm = T)),
                   FishAge_pop_nLL_df %>% dplyr::group_by(prof_val, type) %>% dplyr::summarize(value = sum(value, na.rm = T)),
                   FishLen_pop_nLL_df %>% dplyr::group_by(prof_val, type) %>% dplyr::summarize(value = sum(value, na.rm = T)),
                   FishAge_discard_pop_nLL_df %>% dplyr::group_by(prof_val, type) %>% dplyr::summarize(value = sum(value, na.rm = T)),
                   FishLen_discard_pop_nLL_df %>% dplyr::group_by(prof_val, type) %>% dplyr::summarize(value = sum(value, na.rm = T)),
                   SrvAge_pop_nLL_df  %>% dplyr::group_by(prof_val, type) %>% dplyr::summarize(value = sum(value, na.rm = T)),
                   SrvLen_pop_nLL_df  %>% dplyr::group_by(prof_val, type) %>% dplyr::summarize(value = sum(value, na.rm = T))
  )

  profile_list <- list(jnLL_df = jnLL_df,
                       rec_nLL_df = rec_nLL_df,
                       M_nLL_df = M_nLL_df,
                       sel_nLL_df = sel_nLL_df,
                       rec_prop_nLL_df = rec_prop_nLL_df,
                       Movement_nLL_df = Movement_nLL_df,
                       h_nLL_df = h_nLL_df,
                       TagRep_nLL_df = TagRep_nLL_df,
                       Fmort_nLL_df = Fmort_nLL_df,
                       dmr_nLL_df = dmr_nLL_df,
                       Catch_nLL_df = Catch_nLL_df,
                       Discard_nLL_df = Discard_nLL_df,
                       conv_fish_tag_nLL_df = conv_fish_tag_nLL_df,
                       FishAge_nLL_df = FishAge_nLL_df,
                       FishAgeComps_discard_nLL_df = FishAgeComps_discard_nLL_df,
                       SrvAge_nLL_df = SrvAge_nLL_df,
                       FishLen_nLL_df = FishLen_nLL_df,
                       FishLenComps_discard_nLL_df = FishLenComps_discard_nLL_df,
                       SrvLen_nLL_df = SrvLen_nLL_df,
                       FishIdx_nLL_df = FishIdx_nLL_df,
                       SrvIdx_nLL_df = SrvIdx_nLL_df,
                       agg_nLL = agg_nLL,
                       fish_q_nLL_df    = fish_q_nLL_df,
                       srv_q_nLL_df     = srv_q_nLL_df,
                       Catch_pop_nLL_df   = Catch_pop_nLL_df,
                       Discard_pop_nLL_df   = Discard_pop_nLL_df,
                       FishIdx_pop_nLL_df = FishIdx_pop_nLL_df,
                       SrvIdx_pop_nLL_df  = SrvIdx_pop_nLL_df,
                       FishAge_pop_nLL_df = FishAge_pop_nLL_df,
                       FishLen_pop_nLL_df = FishLen_pop_nLL_df,
                       FishAge_discard_pop_nLL_df = FishAge_discard_pop_nLL_df,
                       FishLen_discard_pop_nLL_df = FishLen_discard_pop_nLL_df,
                       SrvAge_pop_nLL_df  = SrvAge_pop_nLL_df,
                       SrvLen_pop_nLL_df  = SrvLen_pop_nLL_df
  )

  return(profile_list)
}
