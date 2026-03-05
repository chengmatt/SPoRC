#' Run Likelihood Profile
#'
#' @param data data list from model
#' @param parameters parameter list from model
#' @param mapping mapping list from model
#' @param random character vector of random effects to estimate
#' @param what parameter name we want to profile
#' @param idx Index for an parameter array, pointing to the value we want to map off (index is relative to a flattened array)
#' @param min_val minimum value of profile
#' @param max_val maximum value of profile
#' @param inc increment value between min and max value
#' @param do_par logical, whether to use parallel processing (default FALSE)
#' @param n_cores integer, number of cores to use for parallel processing (default is detectCores() - 1)
#'
#' @import dplyr
#' @import RTMB
#' @importFrom reshape2 melt
#' @importFrom stats rnorm nlminb
#' @importFrom future plan multisession
#' @importFrom future.apply future_lapply
#' @importFrom progressr with_progress progressor
#' @importFrom parallel detectCores
#' @returns Returns a list of likelihood profiled values for each data component with their respective dimensions (e.g., likelihood profiles by fleet, region, year, etc.) as well likelihood profiles for each data component, aggregated across all their respective dimensions.
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
  Rec_prop_nLL <- matrix(NA, nrow = length(vals), ncol = 1, dimnames = list(vals, NULL))
  h_nLL <- matrix(NA, nrow = length(vals), ncol = 1, dimnames = list(vals, NULL))
  Movement_nLL <- matrix(NA, nrow = length(vals), ncol = 1, dimnames = list(vals, NULL))
  TagRep_nLL <- matrix(NA, nrow = length(vals), ncol = 1, dimnames = list(vals, NULL))
  Fmort_nLL <- matrix(NA, nrow = length(vals), ncol = 1, dimnames = list(vals, NULL))
  Conv_Fish_Tag_nLL <- data.frame()
  Catch_nLL <- data.frame()
  FishAge_nLL <- data.frame()
  SrvAge_nLL <- data.frame()
  SrvLen_nLL <- data.frame()
  FishLen_nLL <- data.frame()
  FishIdx_nLL <- data.frame()
  SrvIdx_nLL <- data.frame()

  # If there is more than one value in this parameter
  if(do_par == FALSE) {
    for(j in 1:length(vals)) {
      if(!is.null(dim(parameters[[what]]))) {

        # Input fixed value into parameter list
        parameters[[what]] <- do.call(`[<-`, c(list(parameters[[what]]), idx, list(vals[j])))

        # Now, figure out which parameter to map off
        counter <- 1
        map_parameter <- do.call(`[<-`, c(list(parameters[[what]]), idx, list(NA))) # input NA
        non_na <- which(!is.na(map_parameter)) # figure out non NA positions and fill with unique numbers
        for(i in 1:length(non_na)) {
          map_parameter[non_na[i]] <- counter
          counter <- counter + 1 # update counter
        } # end i
        mapping[[what]] <- factor(map_parameter) # input NA into mapping list

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

        # Store values and save
        jnLL[j,1] <- report$jnLL
        rec_nLL[j,1] <- sum(data$Wt_Rec * report$Init_Rec_nLL, data$Wt_Rec * report$Rec_nLL)
        M_nLL[j,1] <- report$M_nLL
        sel_nLL[j,1] <- report$sel_nLL
        Rec_prop_nLL[j,1] <- report$Rec_prop_nLL
        Movement_nLL[j,1] <- report$Movement_nLL
        h_nLL[j,1] <- report$h_nLL
        TagRep_nLL[j,1] <- report$TagRep_nLL
        Fmort_nLL[j,1] <- sum(data$Wt_F * report$Fmort_nLL)
        Conv_Fish_Tag_nLL <- rbind(Conv_Fish_Tag_nLL, reshape2::melt(data$Wt_Tagging * report$Conv_Fish_Tag_nLL) %>% dplyr::mutate(prof_val = vals[j]))
        Catch_nLL <- rbind(Catch_nLL, reshape2::melt(data$Wt_Catch * report$Catch_nLL) %>% dplyr::mutate(prof_val = vals[j]))
        FishAge_nLL <- rbind(FishAge_nLL, reshape2::melt(report$FishAgeComps_nLL) %>% dplyr::mutate(prof_val = vals[j]))
        SrvAge_nLL <- rbind(SrvAge_nLL, reshape2::melt(report$SrvAgeComps_nLL) %>% dplyr::mutate(prof_val = vals[j]))
        SrvLen_nLL <- rbind(SrvLen_nLL, reshape2::melt(report$SrvLenComps_nLL) %>% dplyr::mutate(prof_val = vals[j]))
        FishLen_nLL <- rbind(FishLen_nLL, reshape2::melt(report$FishLenComps_nLL) %>% dplyr::mutate(prof_val = vals[j]))
        FishIdx_nLL <- rbind(FishIdx_nLL, reshape2::melt(data$Wt_FishIdx * report$FishIdx_nLL) %>% dplyr::mutate(prof_val = vals[j]))
        SrvIdx_nLL <- rbind(SrvIdx_nLL, reshape2::melt(data$Wt_SrvIdx * report$SrvIdx_nLL) %>% dplyr::mutate(prof_val = vals[j]))

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
          Rec_prop_nLL = NA,
          Movement_nLL = NA,
          h_nLL = NA,
          TagRep_nLL = NA,
          Fmort_nLL = NA,
          Conv_Fish_Tag_nLL = data.frame(),
          Catch_nLL = data.frame(),
          FishAge_nLL = data.frame(),
          SrvAge_nLL = data.frame(),
          SrvLen_nLL = data.frame(),
          FishLen_nLL = data.frame(),
          FishIdx_nLL = data.frame(),
          SrvIdx_nLL = data.frame()
        )

        if(!is.null(dim(local_parameters[[what]]))) {

          # Input fixed value into parameter list
          local_parameters[[what]] <- do.call(`[<-`, c(list(local_parameters[[what]]), idx, list(vals[j])))

          # Now, figure out which parameter to map off
          counter <- 1
          map_parameter <- do.call(`[<-`, c(list(local_parameters[[what]]), idx, list(NA))) # input NA
          non_na <- which(!is.na(map_parameter)) # figure out non NA positions and fill with unique numbers
          for(i in 1:length(non_na)) {
            map_parameter[non_na[i]] <- counter
            counter <- counter + 1 # update counter
          } # end i
          local_mapping[[what]] <- factor(map_parameter) # input NA into mapping list

        } else { # else, there is only one value in this parameter
          local_parameters[[what]] <- vals[j]
          local_mapping[[what]] <- factor(NA)
        }

        # make adfun
        tryCatch({
          SPoRC_rtmb_model <- RTMB::MakeADFun(cmb(SPoRC_rtmb, data), parameters = local_parameters, map = local_mapping, random = random, silent = TRUE)

          SPoRC_optim <- stats::nlminb(SPoRC_rtmb_model$par, SPoRC_rtmb_model$fn, SPoRC_rtmb_model$gr,
                                       control = list(iter.max = 1e6, eval.max = 1e6, rel.tol = 1e-15))

          # Get report
          report <- SPoRC_rtmb_model$report(SPoRC_rtmb_model$env$last.par.best)

          # Store values and save
          result$jnLL <- report$jnLL
          result$rec_nLL <- sum(data$Wt_Rec * report$Init_Rec_nLL, data$Wt_Rec * report$Rec_nLL)
          result$M_nLL <- report$M_nLL
          result$sel_nLL <- report$sel_nLL
          result$Rec_prop_nLL <- report$Rec_prop_nLL
          result$Movement_nLL <- report$Movement_nLL
          result$h_nLL <- report$h_nLL
          result$TagRep_nLL <- report$TagRep_nLL
          result$Fmort_nLL <- sum(data$Wt_F * report$Fmort_nLL)
          result$Conv_Fish_Tag_nLL <- reshape2::melt(data$Wt_Tagging * report$Conv_Fish_Tag_nLL) %>% dplyr::mutate(prof_val = vals[j])
          result$Catch_nLL <- reshape2::melt(data$Wt_Catch * report$Catch_nLL) %>% dplyr::mutate(prof_val = vals[j])
          result$FishAge_nLL <- reshape2::melt(report$FishAgeComps_nLL) %>% dplyr::mutate(prof_val = vals[j])
          result$SrvAge_nLL <- reshape2::melt(report$SrvAgeComps_nLL) %>% dplyr::mutate(prof_val = vals[j])
          result$SrvLen_nLL <- reshape2::melt(report$SrvLenComps_nLL) %>% dplyr::mutate(prof_val = vals[j])
          result$FishLen_nLL <- reshape2::melt(report$FishLenComps_nLL) %>% dplyr::mutate(prof_val = vals[j])
          result$FishIdx_nLL <- reshape2::melt(data$Wt_FishIdx * report$FishIdx_nLL) %>% dplyr::mutate(prof_val = vals[j])
          result$SrvIdx_nLL <- reshape2::melt(data$Wt_SrvIdx * report$SrvIdx_nLL) %>% dplyr::mutate(prof_val = vals[j])
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
        Rec_prop_nLL[j,1] <- res$Rec_prop_nLL
        Movement_nLL[j,1] <- res$Movement_nLL
        h_nLL[j,1] <- res$h_nLL
        TagRep_nLL[j,1] <- res$TagRep_nLL
        Fmort_nLL[j,1] <- res$Fmort_nLL
        Conv_Fish_Tag_nLL <- rbind(Conv_Fish_Tag_nLL, res$Conv_Fish_Tag_nLL)
        Catch_nLL <- rbind(Catch_nLL, res$Catch_nLL)
        FishAge_nLL <- rbind(FishAge_nLL, res$FishAge_nLL)
        SrvAge_nLL <- rbind(SrvAge_nLL, res$SrvAge_nLL)
        SrvLen_nLL <- rbind(SrvLen_nLL, res$SrvLen_nLL)
        FishLen_nLL <- rbind(FishLen_nLL, res$FishLen_nLL)
        FishIdx_nLL <- rbind(FishIdx_nLL, res$FishIdx_nLL)
        SrvIdx_nLL <- rbind(SrvIdx_nLL, res$SrvIdx_nLL)
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
  Rec_prop_nLL_df <- reshape2::melt(Rec_prop_nLL) %>%
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
  Catch_nLL_df <- Catch_nLL %>%
    dplyr::rename(Region = Var1, Year = Var2, Seas = Var3, Fleet = Var4) %>%
    dplyr::mutate(type = 'Catch')
  FishAge_nLL_df <- FishAge_nLL %>%
    dplyr::rename(Region = Var1, Year = Var2, Seas = Var3, Sex = Var4, Fleet = Var5) %>%
    dplyr::mutate(type = 'FishAge')
  SrvAge_nLL_df <- SrvAge_nLL %>%
    dplyr::rename(Region = Var1, Year = Var2, Seas = Var3, Sex = Var4, Fleet = Var5) %>%
    dplyr::mutate(type = 'SrvAge')
  FishLen_nLL_df <- FishLen_nLL %>%
    dplyr::rename(Region = Var1, Year = Var2, Seas = Var3, Sex = Var4, Fleet = Var5) %>%
    dplyr::mutate(type = 'FishLen')
  SrvLen_nLL_df <- SrvLen_nLL %>%
    dplyr::rename(Region = Var1, Year = Var2, Seas = Var3, Sex = Var4, Fleet = Var5) %>%
    dplyr::mutate(type = 'SrvLen')
  FishIdx_nLL_df <- FishIdx_nLL %>%
    dplyr::rename(Region = Var1, Year = Var2, Seas = Var3, Fleet = Var4) %>%
    dplyr::mutate(type = 'FishIdx')
  SrvIdx_nLL_df <- SrvIdx_nLL %>%
    dplyr::rename(Region = Var1, Year = Var2, Seas = Var3, Fleet = Var4) %>%
    dplyr::mutate(type = 'SrvIdx')
  Conv_Fish_Tag_nLL_df <- Conv_Fish_Tag_nLL %>%
    dplyr::rename(Recap_Year = Var1, Recap_Seas = Var2, Tag_Cohort = Var3, Region = Var4, Age = Var5, Sex = Var6) %>%
    dplyr::mutate(type = 'Tagging')

  # Get likelihoods aggregated across all dimensions
  agg_nLL <- rbind(jnLL_df, rec_nLL_df, M_nLL_df, Rec_prop_nLL_df, Movement_nLL_df, h_nLL_df,
                   TagRep_nLL_df,Fmort_nLL_df, sel_nLL_df,
                   Catch_nLL_df %>% dplyr::group_by(prof_val, type) %>%
                     dplyr::summarize(value = sum(value)),
                   Conv_Fish_Tag_nLL_df %>% dplyr::group_by(prof_val, type) %>%
                     dplyr::summarize(value = sum(value, na.rm = T)),
                   FishAge_nLL_df %>% dplyr::group_by(prof_val, type) %>%
                     dplyr::summarize(value = sum(value, na.rm = T)),
                   SrvAge_nLL_df %>% dplyr::group_by(prof_val, type) %>%
                     dplyr::summarize(value = sum(value, na.rm = T)),
                   FishLen_nLL_df %>% dplyr::group_by(prof_val, type) %>%
                     dplyr::summarize(value = sum(value, na.rm = T)),
                   SrvLen_nLL_df %>% dplyr::group_by(prof_val, type) %>%
                     dplyr::summarize(value = sum(value, na.rm = T)),
                   FishIdx_nLL_df %>% dplyr::group_by(prof_val, type) %>%
                     dplyr::summarize(value = sum(value, na.rm = T)),
                   SrvIdx_nLL_df %>% dplyr::group_by(prof_val, type) %>%
                     dplyr::summarize(value = sum(value, na.rm = T))
  )

  profile_list <- list(jnLL_df = jnLL_df,
                       rec_nLL_df = rec_nLL_df,
                       M_nLL_df = M_nLL_df,
                       sel_nLL_df = sel_nLL_df,
                       Rec_prop_nLL_df = Rec_prop_nLL_df,
                       Movement_nLL_df = Movement_nLL_df,
                       h_nLL_df = h_nLL_df,
                       TagRep_nLL_df = TagRep_nLL_df,
                       Fmort_nLL_df = Fmort_nLL_df,
                       Catch_nLL_df = Catch_nLL_df,
                       Conv_Fish_Tag_nLL_df = Conv_Fish_Tag_nLL_df,
                       FishAge_nLL_df = FishAge_nLL_df,
                       SrvAge_nLL_df = SrvAge_nLL_df,
                       FishLen_nLL_df = FishLen_nLL_df,
                       SrvLen_nLL_df = SrvLen_nLL_df,
                       FishIdx_nLL_df = FishIdx_nLL_df,
                       SrvIdx_nLL_df = SrvIdx_nLL_df,
                       agg_nLL = agg_nLL
  )

  return(profile_list)
}





