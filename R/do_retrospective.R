#' Truncate Model Inputs for Retrospective Diagnostics
#'
#' Internal helper used by \code{do_retrospective()} to truncate model inputs
#' when conducting retrospective diagnostics. The function removes the last
#' \code{j} years from the terminal portion of the time series and updates all
#' associated data objects, parameter arrays, and parameter mappings so that
#' their dimensions remain internally consistent.
#'
#' Specifically, the function adjusts the model \code{data}, \code{parameters},
#' and \code{mapping} lists used by the RTMB model by:
#' \itemize{
#'   \item Truncating the \code{years} vector.
#'   \item Removing terminal years from observations (catch, indices, and
#'   composition data).
#'   \item Truncating time-varying parameter arrays (e.g., recruitment
#'   deviations, fishing mortality deviations, selectivity deviations,
#'   movement parameters).
#'   \item Updating parameter mappings to match the truncated parameter
#'   dimensions.
#'   \item Adjusting block structures and auxiliary objects that depend on
#'   the number of modeled years.
#' }
#'
#' The resulting objects can be passed directly to the model to fit a
#' retrospective peel.
#'
#' @param j Integer specifying the number of terminal years to remove from the
#'   dataset. A value of \code{0} returns the full dataset with no truncation.
#' @param data List containing model data supplied to the RTMB model.
#' @param parameters List containing model parameters supplied to the RTMB
#'   model.
#' @param mapping List defining parameter mappings used during estimation.
#'
#' @returns A list containing truncated versions of the RTMB inputs:
#' \itemize{
#'   \item \code{retro_data} – Modified data list with terminal years removed.
#'   \item \code{retro_parameters} – Parameter list truncated to match the
#'   shortened time series.
#'   \item \code{retro_mapping} – Mapping list updated to match truncated
#'   parameter dimensions.
#' }
#'
#' @keywords internal
truncate_yr <- function(j,
                        data,
                        parameters,
                        mapping) {

  # set up retro data, parameters, and mapping
  retro_data <- data
  retro_parameters <- parameters
  retro_mapping <- mapping

  # Years
  retro_data$years <- data$years[1:(length(data$years) - j)] # remove j years from years vector
  if(!is.na(sum(retro_data$bias_year))) retro_data$bias_year[3:4] <- data$bias_year[3:4] - j # remove j years from bias correction vector (only applied to the full bias and descending limb)

# Recruitment -------------------------------------------------------------

  # Recruitment devs
  retro_parameters$ln_RecDevs <- parameters$ln_RecDevs[,,1:(dim(parameters$ln_RecDevs)[3] - j), drop = FALSE] # Recruitment deviations
  if(any(names(retro_mapping) == 'ln_RecDevs')) retro_mapping$ln_RecDevs <- factor(array(mapping$ln_RecDevs, dim = dim(parameters$ln_RecDevs))[,,1:(dim(parameters$ln_RecDevs)[3] - j), drop = FALSE]) # modify mapping if we have recruitment map

# Fishery -----------------------------------------------------------------

  # Catch, Fishery Index, and Compositions
  retro_data$ObsCatch <- data$ObsCatch[,1:(length(data$years) - j),,,drop = FALSE]
  retro_data$ObsFishIdx <- data$ObsFishIdx[,1:(length(data$years) - j),,,drop = FALSE]
  retro_data$ObsFishIdx_SE <- data$ObsFishIdx_SE[,1:(length(data$years) - j),,,drop = FALSE]
  retro_data$ObsFishAgeComps <- data$ObsFishAgeComps[,1:(length(data$years) - j),,,,,drop = FALSE]
  retro_data$ObsFishLenComps <- data$ObsFishLenComps[,1:(length(data$years) - j),,,,,drop = FALSE]

  # Fishery mortality devs
  retro_parameters$ln_F_devs <- parameters$ln_F_devs[,1:(length(data$years) - j),,,drop = FALSE] # modify F dev parameters
  retro_mapping$ln_F_devs <- factor(array(mapping$ln_F_devs, dim = dim(parameters$ln_F_devs))[,1:(length(data$years) - j),,,drop = FALSE]) # modify map

  # Fishery selectivity deviations
  retro_parameters$ln_fishsel_devs <- parameters$ln_fishsel_devs[,1:(length(data$years) - j),,,,drop = FALSE] # modify parameter length
  retro_mapping$ln_fishsel_devs <- factor(array(mapping$ln_fishsel_devs, dim = dim(parameters$ln_fishsel_devs))[,1:(length(data$years) - j),,,,drop = FALSE]) # modify map
  retro_data$map_ln_fishsel_devs <- data$map_ln_fishsel_devs[,1:(length(data$years) - j),,,,drop = FALSE]

  # Fishery selectivity and catchability blocks
  retro_data$fish_q_blocks <- data$fish_q_blocks[,1:(length(data$years) - j),, drop = FALSE]
  retro_data$fish_sel_blocks <- data$fish_sel_blocks[,1:(length(data$years) - j),, drop = FALSE]

  # Adjust fishery parameter blocks
  retro_parameters$ln_fish_q <- parameters$ln_fish_q[,1:max(retro_data$fish_q_blocks),,drop = FALSE]
  retro_parameters$ln_fish_fixed_sel_pars <- parameters$ln_fish_fixed_sel_pars[,,1:max(retro_data$fish_sel_blocks),,,drop = FALSE]

  # Adjust fishery mapping
  retro_mapping$ln_fish_q <- factor(array(mapping$ln_fish_q, dim = dim(parameters$ln_fish_q))[,1:max(retro_data$fish_q_blocks),,drop = FALSE])
  retro_mapping$ln_fish_fixed_sel_pars <- factor(array(mapping$ln_fish_fixed_sel_pars, dim = dim(parameters$ln_fish_fixed_sel_pars))[,,1:max(retro_data$fish_sel_blocks),,,drop = FALSE])

# Survey ------------------------------------------------------------------

  # Survey index and compositions
  retro_data$ObsSrvIdx <- data$ObsSrvIdx[,1:(length(data$years) - j),,,drop = FALSE]
  retro_data$ObsSrvIdx_SE <- data$ObsSrvIdx_SE[,1:(length(data$years) - j),,,drop = FALSE]
  retro_data$ObsSrvAgeComps <- data$ObsSrvAgeComps[,1:(length(data$years) - j),,,,,drop = FALSE]
  retro_data$ObsSrvLenComps <- data$ObsSrvLenComps[,1:(length(data$years) - j),,,,,drop = FALSE]

  # Survey selectivity deviations
  retro_parameters$ln_srvsel_devs <- parameters$ln_srvsel_devs[,1:(length(data$years) - j),,,,drop = FALSE] # Survey selectivity deviations
  retro_mapping$ln_srvsel_devs <- factor(array(mapping$ln_srvsel_devs, dim = dim(parameters$ln_srvsel_devs))[,1:(length(data$years) - j),,,,drop = FALSE]) # modify map
  retro_data$map_ln_srvsel_devs <- data$map_ln_srvsel_devs[,1:(length(data$years) - j),,,,drop = FALSE]

  # Survey selectivity and catchability blocks
  retro_data$srv_q_blocks <- data$srv_q_blocks[,1:(length(data$years) - j),, drop = FALSE]
  retro_data$srv_sel_blocks <- data$srv_sel_blocks[,1:(length(data$years) - j),, drop = FALSE]

  # Adjust survey parameter blocks
  retro_parameters$ln_srv_q <- parameters$ln_srv_q[,1:max(retro_data$srv_q_blocks),,drop = FALSE]
  retro_parameters$ln_srv_fixed_sel_pars <- parameters$ln_srv_fixed_sel_pars[,,1:max(retro_data$srv_sel_blocks),,,drop = FALSE]

  # Adjust survey mapping
  retro_mapping$ln_srv_q <- factor(array(mapping$ln_srv_q, dim = dim(parameters$ln_srv_q))[,1:max(retro_data$srv_q_blocks),,drop = FALSE])
  retro_mapping$ln_srv_fixed_sel_pars <- factor(array(mapping$ln_srv_fixed_sel_pars, dim = dim(parameters$ln_srv_fixed_sel_pars))[,,1:max(retro_data$srv_sel_blocks),,,drop = FALSE])

# Movement ----------------------------------------------------------------

  if(data$n_regions > 1) {
    # Movement stuff
    retro_parameters$move_pars <- parameters$move_pars[,,,1:(length(data$years) - j),,,,drop = FALSE]
    retro_parameters$move_devs <- parameters$move_devs[,,,1:(length(data$years) - j),,,,drop = FALSE]
    retro_mapping$move_pars <- factor(array(mapping$move_pars, dim = dim(parameters$move_pars))[,,,1:(length(data$years) - j),,,,drop = FALSE])
    retro_mapping$move_devs <- factor(array(mapping$move_devs, dim = dim(parameters$move_devs))[,,,1:(length(data$years) - j),,,,drop = FALSE])
    retro_data$Fixed_Movement <- data$Fixed_Movement[,,,1:(length(data$years) - j),,,,drop = FALSE]
    retro_data$sgl_seas_spawning_movement <- data$sgl_seas_spawning_movement[,,,1:(length(data$years) - j),,,drop = FALSE]
    retro_data$ctmc_move_dat <- retro_data$ctmc_move_dat[which(data$ctmc_move_dat$years %in% 1:(length(data$years) - j)),]
    retro_data$map_move_devs <- retro_data$map_move_devs[,,,1:(length(data$years) - j),,,,drop = FALSE]
  }

# Tagging -----------------------------------------------------------------

  if(any(data$use_conv_fish_tagging == 1)) {
    # Tag reporting
    retro_data$conv_tag_fish_reporting_blocks <- data$conv_tag_fish_reporting_blocks[,1:(length(data$years) - j),, drop = FALSE]
    if(!is.na(sum(retro_data$conv_tag_fish_reporting_blocks))) {
      retro_parameters$conv_tag_fish_reporting_pars <- parameters$conv_tag_fish_reporting_pars[,1:max(retro_data$conv_tag_fish_reporting_blocks),,drop = FALSE]
      retro_mapping$conv_tag_fish_reporting_pars <- factor(array(mapping$conv_tag_fish_reporting_pars, dim = dim(parameters$conv_tag_fish_reporting_pars))[,1:max(retro_data$conv_tag_fish_reporting_blocks),,drop = FALSE])
    }

    # Tag cohort stuff
    Tag_Release_Ind <- as.matrix(data$conv_tag_release_indicator)
    retro_data$conv_tag_release_indicator <- as.matrix(Tag_Release_Ind[which(Tag_Release_Ind[,2] %in% 1:(length(data$years) - j)), ])
    retro_data$conv_tag_release_platform <- data$conv_tag_release_platform[1:nrow(retro_data$conv_tag_release_indicator),]
    retro_data$n_conv_tag_cohorts <- nrow(retro_data$conv_tag_release_indicator)
    retro_data$conv_tagged_fish <- data$conv_tagged_fish[1:nrow(retro_data$conv_tag_release_indicator),,,,drop = FALSE] # remove data (not necessary, but helps with computational cost if using tagging)
    retro_data$obs_conv_tag_fish_recap <- data$obs_conv_tag_fish_recap[,,1:nrow(retro_data$conv_tag_release_indicator),,,,,,drop = FALSE] # remove data (not necessary, but helps with computational cost)
  }


# Data Weights, Composition Stuff, and use indicators ------------------------------------------------------------

  # Data weights and composition stuff
  retro_data$Wt_FishAgeComps <- data$Wt_FishAgeComps[,1:(length(data$years) - j),,,,drop = FALSE]
  retro_data$Wt_SrvAgeComps <- data$Wt_SrvAgeComps[,1:(length(data$years) - j),,,,drop = FALSE]
  retro_data$Wt_FishLenComps <- data$Wt_FishLenComps[,1:(length(data$years) - j),,,,drop = FALSE]
  retro_data$Wt_SrvLenComps <- data$Wt_SrvLenComps[,1:(length(data$years) - j),,,,drop = FALSE]
  retro_data$FishAgeComps_Type <- data$FishAgeComps_Type[1:(length(data$years) - j),,drop = FALSE]
  retro_data$FishLenComps_Type <- data$FishLenComps_Type[1:(length(data$years) - j),,drop = FALSE]
  retro_data$SrvLenComps_Type <- data$SrvLenComps_Type[1:(length(data$years) - j),,drop = FALSE]
  retro_data$SrvAgeComps_Type <- data$SrvAgeComps_Type[1:(length(data$years) - j),,drop = FALSE]
  if(length(dim(data$Wt_Catch)) == 4) retro_data$Wt_Catch <- data$Wt_Catch[,1:(length(data$years) - j),,,drop = FALSE] # Catch is dim = 4, b/c can accept scalar or array

  # data use indicators
  retro_data$UseFishAgeComps <- data$UseFishAgeComps[,1:(length(data$years) - j),,,drop = FALSE]
  retro_data$UseFishIdx <- data$UseFishIdx[,1:(length(data$years) - j),,,drop = FALSE]
  retro_data$UseCatch <- data$UseCatch[,1:(length(data$years) - j),,,drop = FALSE]
  retro_data$UseFishLenComps <- data$UseFishLenComps[,1:(length(data$years) - j),,,drop = FALSE]
  retro_data$UseSrvAgeComps <- data$UseSrvAgeComps[,1:(length(data$years) - j),,,drop = FALSE]
  retro_data$UseSrvIdx <- data$UseSrvIdx[,1:(length(data$years) - j),,,drop = FALSE]
  retro_data$UseSrvLenComps <- data$UseSrvLenComps[,1:(length(data$years) - j),,,drop = FALSE]


  return(list(retro_data = retro_data,
              retro_parameters = retro_parameters,
              retro_mapping = retro_mapping))
}



#' Run Retrospective Diagnostics for RTMB Models
#'
#' Conducts retrospective analyses by sequentially removing terminal years
#' ("peels") from the dataset and refitting the model. For each peel, the
#' function truncates the model inputs, optionally applies data lags and
#' Francis composition reweighting, fits the model, and extracts estimates
#' of spawning stock biomass (SSB) and recruitment.
#'
#' Retrospective analyses are commonly used to evaluate the stability of
#' model estimates through time and to diagnose potential model misspecification.
#'
#' @param n_retro Integer specifying the number of retrospective peels to perform.
#'   A value of \code{n_retro = 0} fits the model using the full dataset only.
#' @param data List containing the data supplied to the RTMB model.
#' @param parameters List containing the model parameters.
#' @param mapping List defining parameter mappings used during estimation.
#' @param random Character vector identifying random-effect parameters in the model.
#'   Default is \code{NULL}.
#' @param do_par Logical indicating whether retrospective peels should be run
#'   in parallel. Default is \code{FALSE}.
#' @param n_cores Integer specifying the number of cores to use when
#'   \code{do_par = TRUE}.
#' @param newton_loops Integer specifying the number of Newton optimization
#'   loops used during model fitting. Default is \code{3}.
#' @param do_francis Logical indicating whether Francis composition
#'   reweighting should be applied within each retrospective peel.
#'   Default is \code{FALSE}.
#' @param n_francis_iter Integer specifying the number of Francis reweighting
#'   iterations. Required if \code{do_francis = TRUE}.
#' @param nlminb_control List of control arguments passed to \code{stats::nlminb}
#'   during model fitting. Default is
#'   \code{list(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-15)}.
#' @param do_sdrep Logical indicating whether standard errors should be
#'   calculated using \code{RTMB::sdreport}. Default is \code{FALSE}.
#'
#' @param fishidx_datalag Integer array specifying lags applied to fishery
#'   index data \eqn{[region \times fleet]}. Default is zeros.
#' @param fishage_datalag Integer array specifying lags applied to fishery
#'   age-composition data \eqn{[region \times fleet]}. Default is zeros.
#' @param fishlen_datalag Integer array specifying lags applied to fishery
#'   length-composition data \eqn{[region \times fleet]}. Default is zeros.
#' @param srvidx_datalag Integer array specifying lags applied to survey
#'   index data \eqn{[region \times fleet]}. Default is zeros.
#' @param srvage_datalag Integer array specifying lags applied to survey
#'   age-composition data \eqn{[region \times fleet]}. Default is zeros.
#' @param srvlen_datalag Integer array specifying lags applied to survey
#'   length-composition data \eqn{[region \times fleet]}. Default is zeros.
#' @param conv_tag_datalag Integer specifying the lag applied to conventional
#'   tagging data. Default is \code{0}.
#'
#' @return A long-format \code{data.frame} containing retrospective estimates
#'   of spawning stock biomass and recruitment. Columns include:
#'   \itemize{
#'     \item \code{Pop} – Population index.
#'     \item \code{Region} – Region index.
#'     \item \code{Year} – Model year.
#'     \item \code{Type} – Quantity reported ("SSB" or "Recruitment").
#'     \item \code{peel} – Retrospective peel number (0 = full data, 1 = one-year peel, etc.).
#'     \item \code{value} – Estimated value of the quantity.
#'     \item \code{pdHess} – Logical indicator of positive-definite Hessian
#'     (returned when \code{do_sdrep = TRUE}).
#'     \item \code{max_grad} – Maximum absolute gradient of fixed effects
#'     (returned when \code{do_sdrep = TRUE}).
#'   }
#'
#' @export do_retrospective
#' @family Model Diagnostics
#'
#' @import RTMB
#' @import dplyr
#' @import future.apply
#' @import future
#' @import progressr
#' @importFrom reshape2 melt
#' @importFrom stats nlminb optimHess
do_retrospective <- function(n_retro,
                             data,
                             parameters,
                             mapping,
                             random = NULL,
                             do_par,
                             n_cores,
                             newton_loops = 3,
                             do_francis = FALSE,
                             n_francis_iter = NULL,
                             nlminb_control = list(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-15),
                             do_sdrep = FALSE,
                             fishidx_datalag = array(0, dim = c(data$n_regions, data$n_fish_fleets)),
                             fishage_datalag = array(0, dim = c(data$n_regions, data$n_fish_fleets)),
                             fishlen_datalag = array(0, dim = c(data$n_regions, data$n_fish_fleets)),
                             srvidx_datalag = array(0, dim = c(data$n_regions, data$n_srv_fleets)),
                             srvage_datalag = array(0, dim = c(data$n_regions, data$n_srv_fleets)),
                             srvlen_datalag = array(0, dim = c(data$n_regions, data$n_srv_fleets)),
                             conv_tag_datalag = 0
                             ) {

  # Loop through retrospective (no parrallelization)
  if(do_par == FALSE) {

    retro_all <- data.frame()

    for(j in 0:n_retro) {

      # truncate data
      init <- truncate_yr(j = j, data = data, parameters = parameters, mapping = mapping)
      init$retro_data$conv_tagged_fish
      SPoRC_rtmb(pars = init$retro_parameters, data = init$retro_data)

      # Fishery Data Lags
      start_col <- length(init$retro_data$years) # get start index
      for(f in 1:data$n_fish_fleets) {
        for(r in 1:data$n_regions) {
          # get lag indices
          fishage_tmp_lag <- fishage_datalag[r,f]
          fishlen_tmp_lag <- fishlen_datalag[r,f]
          fishidx_tmp_lag <- fishidx_datalag[r,f]
          # fishery ages
          if(fishage_tmp_lag > 0) {
            fishage_end_col <- max(start_col - fishage_tmp_lag + 1, 1) # get end index
            init$retro_data$UseFishAgeComps[r, start_col:fishage_end_col,, f] <- 0 # input 0 to lag incoming data into assessment
          }
          # fishery lengths
          if(fishlen_tmp_lag > 0) {
            fishlen_end_col <- max(start_col - fishlen_tmp_lag + 1, 1) # get end index
            init$retro_data$UseFishLenComps[r, start_col:fishlen_end_col,, f] <- 0 # input 0 to lag incoming data into assessment
          }
          # fishery index
          if(fishidx_tmp_lag > 0) {
            fishidx_end_col <- max(start_col - fishidx_tmp_lag + 1, 1) # get end index
            init$retro_data$UseFishIdx[r, start_col:fishidx_end_col,, f] <- 0 # input 0 to lag incoming data into assessment
          }
        } # end r loop
      } # end f loop

      # Survey Data Lags
      for(f in 1:data$n_srv_fleets) {
        for(r in 1:data$n_regions) {
          # get lag indices
          srvage_tmp_lag <- srvage_datalag[r,f]
          srvlen_tmp_lag <- srvlen_datalag[r,f]
          srvidx_tmp_lag <- srvidx_datalag[r,f]
          # survey ages
          if(srvage_tmp_lag > 0) {
            srvage_end_col <- max(start_col - srvage_tmp_lag + 1, 1) # get end index
            init$retro_data$UseSrvAgeComps[r, start_col:srvage_end_col,, f] <- 0 # input 0 to lag incoming data into assessment
          }
          # survey lengths
          if(srvlen_tmp_lag > 0) {
            srvlen_end_col <- max(start_col - srvlen_tmp_lag + 1, 1) # get end index
            init$retro_data$UseSrvLenComps[r, start_col:srvlen_end_col,, f] <- 0 # input 0 to lag incoming data into assessment
          }
          # survey index
          if(srvidx_tmp_lag > 0) {
            srvidx_end_col <- max(start_col - srvidx_tmp_lag + 1, 1) # get end index
            init$retro_data$UseSrvIdx[r, start_col:srvidx_end_col,, f] <- 0 # input 0 to lag incoming data into assessment
          }
        } # end r loop
      } # end f loop

      # Tagging Data Lags
      if(conv_tag_datalag > 0) {
        Tag_Release_Ind <- as.matrix(init$retro_data$conv_tag_release_indicator) # get tag release indicator
        tag_end_col <- max(start_col - conv_tag_datalag + 1, 1) # get end index
        init$retro_data$conv_tag_release_indicator <- as.matrix(Tag_Release_Ind[-which(Tag_Release_Ind[,2] %in% start_col:tag_end_col), ]) # remove tag data when lagged
        init$retro_data$n_conv_tag_cohorts <- nrow(init$retro_data$conv_tag_release_indicator)
        init$retro_data$conv_tag_release_platform <- init$retro_data$conv_tag_release_platform[1:nrow(init$retro_data$conv_tag_release_indicator),]
        init$retro_data$conv_tagged_fish <- init$retro_data$conv_tagged_fish[1:nrow(init$retro_data$conv_tag_release_indicator),,,,drop = FALSE] # remove data (not necessary, but helps with computational cost if using tagging)
        init$retro_data$obs_conv_tag_fish_recap <- init$retro_data$obs_conv_tag_fish_recap[,,1:nrow(init$retro_data$conv_tag_release_indicator),,,,,,drop = FALSE] # remove data (not necessary, but helps with computational cost)
      }

      if(do_francis == FALSE) { # don't do francis within retrospective loop

        # run model
        SPoRC_rtmb_model <- fit_model(
          data = init$retro_data,
          parameters = init$retro_parameters,
          mapping = init$retro_mapping,
          random = random,
          newton_loops = newton_loops,
          silent = TRUE
        )

        rep <- SPoRC_rtmb_model$rep # extract report

      } else {

        SPoRC_rtmb_model_francis <- run_francis(data = init$retro_data,
                                                parameters = init$retro_parameters,
                                                mapping = init$retro_mapping,
                                                random = random,
                                                n_francis_iter = n_francis_iter,
                                                newton_loops = newton_loops
                                                )

        SPoRC_rtmb_model <- SPoRC_rtmb_model_francis$obj # extract obj
        rep <- SPoRC_rtmb_model_francis$obj$rep # extract report

      } # end else

      # get ssb and recruitment
      retro_tmp <- reshape2::melt(rep$SSB) %>%
        dplyr::rename(Pop = Var1, Region = Var2, Year = Var3) %>%
        dplyr::mutate(Type = "SSB") %>%
        bind_rows(reshape2::melt(rep$Rec) %>%
                    dplyr::rename(Pop = Var1, Region = Var2, Year = Var3) %>%
                    dplyr::mutate(Type = "Recruitment")) %>%
        dplyr::mutate(peel = j)

      if(do_sdrep == TRUE) {
        sdrep <- RTMB::sdreport(SPoRC_rtmb_model) # get sdreport

        # input info about pdHess and gradients
        retro_tmp <- retro_tmp %>%
          dplyr::mutate(pdHess = sdrep$pdHess,
                        max_grad = max(abs(sdrep$gradient.fixed)))
      }

      retro_all <- rbind(retro_all, retro_tmp) # bind all rows

    } # end j
  } # iterative loop


  # Parrallelize Retrospective Loop
  if(do_par == TRUE) {

    future::plan(future::multisession, workers = n_cores) # set up cores

    progressr::with_progress({

      p <- progressr::progressor(along = 0:n_retro) # progress bar

      retro_all <- future.apply::future_lapply(0:n_retro, function(j) {

        # truncate data
        init <- truncate_yr(j = j, data = data, parameters = parameters, mapping = mapping)

        # Fishery Data Lags
        start_col <- length(init$retro_data$years) # get start index
        for(f in 1:data$n_fish_fleets) {
          for(r in 1:data$n_regions) {
            # get lag indices
            fishage_tmp_lag <- fishage_datalag[r,f]
            fishlen_tmp_lag <- fishlen_datalag[r,f]
            fishidx_tmp_lag <- fishidx_datalag[r,f]
            # fishery ages
            if(fishage_tmp_lag > 0) {
              fishage_end_col <- max(start_col - fishage_tmp_lag + 1, 1) # get end index
              init$retro_data$UseFishAgeComps[r, start_col:fishage_end_col,, f] <- 0 # input 0 to lag incoming data into assessment
            }
            # fishery lengths
            if(fishlen_tmp_lag > 0) {
              fishlen_end_col <- max(start_col - fishlen_tmp_lag + 1, 1) # get end index
              init$retro_data$UseFishLenComps[r, start_col:fishlen_end_col,, f] <- 0 # input 0 to lag incoming data into assessment
            }
            # fishery index
            if(fishidx_tmp_lag > 0) {
              fishidx_end_col <- max(start_col - fishidx_tmp_lag + 1, 1) # get end index
              init$retro_data$UseFishIdx[r, start_col:fishidx_end_col,, f] <- 0 # input 0 to lag incoming data into assessment
            }
          } # end r loop
        } # end f loop

        # Survey Data Lags
        for(f in 1:data$n_srv_fleets) {
          for(r in 1:data$n_regions) {
            # get lag indices
            srvage_tmp_lag <- srvage_datalag[r,f]
            srvlen_tmp_lag <- srvlen_datalag[r,f]
            srvidx_tmp_lag <- srvidx_datalag[r,f]
            # survey ages
            if(srvage_tmp_lag > 0) {
              srvage_end_col <- max(start_col - srvage_tmp_lag + 1, 1) # get end index
              init$retro_data$UsesrvAgeComps[r, start_col:srvage_end_col,, f] <- 0 # input 0 to lag incoming data into assessment
            }
            # survey lengths
            if(srvlen_tmp_lag > 0) {
              srvlen_end_col <- max(start_col - srvlen_tmp_lag + 1, 1) # get end index
              init$retro_data$UsesrvLenComps[r, start_col:srvlen_end_col,, f] <- 0 # input 0 to lag incoming data into assessment
            }
            # survey index
            if(srvidx_tmp_lag > 0) {
              srvidx_end_col <- max(start_col - srvidx_tmp_lag + 1, 1) # get end index
              init$retro_data$UsesrvIdx[r, start_col:srvidx_end_col,, f] <- 0 # input 0 to lag incoming data into assessment
            }
          } # end r loop
        } # end f loop

        # Tagging Data Lags
        if(conv_tag_datalag > 0) {
          Tag_Release_Ind <- as.matrix(init$retro_data$conv_tag_release_indicator) # get tag release indicator
          tag_end_col <- max(start_col - conv_tag_datalag + 1, 1) # get end index
          init$retro_data$conv_tag_release_indicator <- as.matrix(Tag_Release_Ind[-which(Tag_Release_Ind[,2] %in% start_col:tag_end_col), ]) # remove tag data when lagged
          init$retro_data$n_conv_tag_cohorts <- nrow(init$retro_data$conv_tag_release_indicator)
          init$retro_data$conv_tag_release_platform <- init$retro_data$conv_tag_release_platform[1:nrow(init$retro_data$conv_tag_release_indicator),]
          init$retro_data$conv_tagged_fish <- init$retro_data$conv_tagged_fish[1:nrow(init$retro_data$conv_tag_release_indicator),,,,drop = FALSE] # remove data (not necessary, but helps with computational cost if using tagging)
          init$retro_data$obs_conv_tag_fish_recap <- init$retro_data$obs_conv_tag_fish_recap[,,1:nrow(init$retro_data$conv_tag_release_indicator),,,,,,drop = FALSE] # remove data (not necessary, but helps with computational cost)
        }

        if(do_francis == FALSE) { # don't do francis within retrospective loop

          # run model
          SPoRC_rtmb_model <- fit_model(
            data = init$retro_data,
            parameters = init$retro_parameters,
            mapping = init$retro_mapping,
            random = random,
            newton_loops = newton_loops,
            silent = TRUE
          )

          rep <- SPoRC_rtmb_model$rep # extract report

        } else {

          SPoRC_rtmb_model_francis <- run_francis(data = init$retro_data,
                                                  parameters = init$retro_parameters,
                                                  mapping = init$retro_mapping,
                                                  random = random,
                                                  n_francis_iter = n_francis_iter,
                                                  newton_loops = newton_loops
          )

          SPoRC_rtmb_model <- SPoRC_rtmb_model_francis$obj # extract obj
          rep <- SPoRC_rtmb_model_francis$obj$rep # extract report

        } # end else

        retro_tmp <- reshape2::melt(rep$SSB) %>%
          dplyr::rename(Pop = Var1, Region = Var2, Year = Var3) %>%
          dplyr::mutate(Type = "SSB") %>%
          bind_rows(reshape2::melt(rep$Rec) %>%
                      dplyr::rename(Pop = Var1, Region = Var2, Year = Var3) %>%
                      dplyr::mutate(Type = "Recruitment")) %>%
          dplyr::mutate(peel = j)

        if(do_sdrep == TRUE) {
          sdrep <- RTMB::sdreport(SPoRC_rtmb_model) # get sdreport

          # input info about pdHess and gradients
          retro_tmp <- retro_tmp %>%
            dplyr::mutate(pdHess = sdrep$pdHess,
                          max_grad = max(abs(sdrep$gradient.fixed)))
        }

        # retro_all <- rbind(retro_all, retro_tmp) # bind all rows

        p() # update progress

        retro_tmp

      }, future.seed = TRUE) %>% bind_rows() # bine rows to combine results

      future::plan(future::sequential)  # Reset

    })
  } # do parrallelization for retrospective loop

  return(retro_all)
} # end function

#' Derive relative difference from terminal year from a retrospective analysis.
#'
#' @param retro_data Dataframe outputted from do_retrospective function
#'
#' @returns Returns a data frame with relative difference of SSB and recruitment from the terminal year
#' @export get_retrospective_relative_difference
#' @family Model Diagnostics
#' @import dplyr
#' @importFrom tidyr pivot_longer pivot_wider
get_retrospective_relative_difference <- function(retro_data) {

  unique_peels <- length(unique(retro_data$peel)) - 1 # get unique peels

  # Get the terminal year assessment
  terminal <- retro_data %>% dplyr::filter(peel == 0)

  # Get peels
  peels <- retro_data %>% filter(peel != 0) %>%
    tidyr::pivot_wider(names_from = peel, values_from = value, id_cols = c('Pop', 'Region', "Year", "Type"))

  # Summarize relative difference
  allret <- terminal %>%
    dplyr::left_join(peels, by = c("Pop", "Region", "Year", "Type")) %>%
    dplyr::mutate(across(as.character(1:unique_peels), ~ (.x - value) / value, .names = "{.col}"))

  # Pivot longer
  allret <- allret %>%
    dplyr::select(Pop, Region, Year, Type, as.character(1:unique_peels)) %>%
    tidyr::pivot_longer(cols = as.character(1:unique_peels), names_to = "peel", values_to = "rd") %>%
    dplyr::mutate(Pop = paste("Pop", Pop), Region = paste("Region", Region))

  return(allret)
}
