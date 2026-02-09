#' Get Time Series Plots
#'
#' @param rep List of n_models of `SPoRC` report lists
#' @param sd_rep List of n_models of `SPoRC` sdreport lists
#' @param model_names Vector of model names
#' @param do_ci Boolean for whether confidence intervals are plotted
#'
#' @returns Plots of spawning biomass, dynamic b0, total biomass, recruitment, and fishing mortality time-series across models
#' @export get_ts_plot
#' @family Plotting
#' @examples
#' \dontrun{
#'   get_ts_plot(list(rep1, rep2), list(sd_rep1, sd_rep2), c("Model1", "Model2"), do_ci = TRUE)
#' }
get_ts_plot <- function(rep,
                        sd_rep,
                        model_names,
                        do_ci = TRUE
                         ) {

  biom_rec_df <- data.frame() # empty dataframe to bind

  for(i in 1:length(rep)) {

    # Spawning Stock Biomass
    ssb_plot_df <- reshape2::melt(rep[[i]]$SSB) %>%
      dplyr::rename(Region = Var1, Year = Var2) %>%
      dplyr::bind_cols(se = sd_rep[[i]]$sd[names(sd_rep[[i]]$value) == "log(SSB)"]) %>%
      dplyr::mutate(lwr = exp(log(value) - 1.96 * se),
                    upr = exp(log(value) + 1.96 * se),
                    Region = paste("Region", Region),
                    Type = 'SSB',
                    Model = model_names[i])

    # Dynamic Unfished Spawning Stock Biomass
    ssb0_plot_df <- reshape2::melt(rep[[i]]$Dynamic_SSB0) %>%
      dplyr::rename(Region = Var1, Year = Var2) %>%
      dplyr::bind_cols(se = sd_rep[[i]]$sd[names(sd_rep[[i]]$value) == "log(Dynamic_SSB0)"]) %>%
      dplyr::mutate(lwr = exp(log(value) - 1.96 * se),
                    upr = exp(log(value) + 1.96 * se),
                    Region = paste("Region", Region),
                    Type = 'Dynamic SSB0',
                    Model = model_names[i])

    # Total Biomass
    totbiom_plot_df <- reshape2::melt(rep[[i]]$Total_Biom) %>%
      dplyr::rename(Region = Var1, Year = Var2) %>%
      dplyr::bind_cols(se = sd_rep[[i]]$sd[names(sd_rep[[i]]$value) == "log(Total_Biom)"]) %>%
      dplyr::mutate(lwr = exp(log(value) - 1.96 * se),
                    upr = exp(log(value) + 1.96 * se),
                    Region = paste("Region", Region),
                    Type = 'Total Biom',
                    Model = model_names[i])

    # Recruitment
    rec_plot_df <- reshape2::melt(rep[[i]]$Rec) %>%
      dplyr::rename(Region = Var1, Year = Var2) %>%
      dplyr::bind_cols(se = sd_rep[[i]]$sd[names(sd_rep[[i]]$value) == "log(Rec)"]) %>%
      dplyr::mutate(lwr = exp(log(value) - 1.96 * se),
                    upr = exp(log(value) + 1.96 * se),
                    Region = paste("Region", Region),
                    Type = 'Recruitment',
                    Model = model_names[i])

    # Fishing Mortality
    f_plot_df <- reshape2::melt(rep[[i]]$Fmort) %>%
      dplyr::rename(Region = Var1, Year = Var2, Seas = Var3, Type = Var4) %>%
      dplyr::mutate(Region = paste("Region", Region),
                    Type = paste("Seas", Seas, "Fleet", Type, "F"),
                    lwr = NA,
                    upr = NA,
                    se = NA,
                    Model = model_names[i]) %>%
      dplyr::select(-Seas)

    # bind together
    biom_rec_df <- rbind(ssb_plot_df, totbiom_plot_df, rec_plot_df, f_plot_df, biom_rec_df, ssb0_plot_df)
  }

  # Plot combined time series
  comb_ts_plot <- ggplot2::ggplot(biom_rec_df,
                             ggplot2::aes(x = Year, y = value, ymin = lwr, ymax = upr, color = factor(Model), fill = factor(Model))) +
    ggplot2::geom_line(lwd = 0.9) +
    ggplot2::facet_grid(Type~Region, scales = 'free') +
    ggplot2::labs(x = 'Year', y = 'Value', color = 'Model', fill = 'Model') +
    ggplot2::coord_cartesian(ylim = c(0,NA)) +
    theme_sablefish()

  # fishing mortality
  f_ts_plot <- ggplot2::ggplot(biom_rec_df %>% dplyr::filter(stringr::str_detect(Type, 'Fleet')),
                               ggplot2::aes(x = Year, y = value, color = factor(Model))) +
    ggplot2::geom_line(lwd = 0.9) +
    ggplot2::facet_grid(Type~Region, scales = 'free') +
    ggplot2::labs(x = 'Year', y = 'Value', color = 'Model', fill = 'Model') +
    ggplot2::coord_cartesian(ylim = c(0,NA)) +
    theme_sablefish()

  # recruitment
  rec_ts_plot <- ggplot2::ggplot(biom_rec_df %>% dplyr::filter(Type == 'Recruitment'),
                                 ggplot2::aes(x = Year, y = value, ymin = lwr, ymax = upr, color = factor(Model), fill = factor(Model))) +
    ggplot2::geom_line(lwd = 0.9) +
    ggplot2::facet_grid(Type~Region, scales = 'free') +
    ggplot2::labs(x = 'Year', y = 'Recruitment', color = 'Model', fill = 'Model') +
    ggplot2::coord_cartesian(ylim = c(0,NA)) +
    theme_sablefish()

  # ssb
  ssb_ts_plot <- ggplot2::ggplot(biom_rec_df %>% dplyr::filter(Type == 'SSB'),
                               ggplot2::aes(x = Year, y = value, ymin = lwr, ymax = upr, color = factor(Model), fill = factor(Model))) +
    ggplot2::geom_line(lwd = 0.9) +
    ggplot2::facet_grid(Type~Region, scales = 'free') +
    ggplot2::labs(x = 'Year', y = 'Spawning Stock Biomass', color = 'Model', fill = 'Model') +
    ggplot2::coord_cartesian(ylim = c(0,NA)) +
    theme_sablefish()

  # total biomass
  total_biom_plot <- ggplot2::ggplot(biom_rec_df %>% dplyr::filter(Type == 'Total Biom'),
                      ggplot2::aes(x = Year, y = value, ymin = lwr, ymax = upr, color = factor(Model), fill = factor(Model))) +
    ggplot2::geom_line(lwd = 0.9) +
    ggplot2::facet_grid(Type~Region, scales = 'free') +
    ggplot2::labs(x = 'Year', y = 'Total Biomass', color = 'Model', fill = 'Model') +
    ggplot2::coord_cartesian(ylim = c(0,NA)) +
    theme_sablefish()

  # dynamic b0
  ssb0_plot <- ggplot2::ggplot(biom_rec_df %>% dplyr::filter(Type == 'Dynamic SSB0'),
                               ggplot2::aes(x = Year, y = value, ymin = lwr, ymax = upr, color = factor(Model), fill = factor(Model))) +
    ggplot2::geom_line(lwd = 0.9) +
    ggplot2::facet_grid(Type~Region, scales = 'free') +
    ggplot2::labs(x = 'Year', y = 'Unfished Spawning Stock Biomass', color = 'Model', fill = 'Model') +
    ggplot2::coord_cartesian(ylim = c(0,NA)) +
    theme_sablefish()

  # dynamic b0 and SSB
  ssb_ssb0_plot <- ggplot2::ggplot(biom_rec_df %>% dplyr::filter(Type %in% c('Dynamic SSB0', "SSB")),
                               ggplot2::aes(x = Year, y = value, color = factor(Model), lty = Type)) +
    ggplot2::geom_line(lwd = 0.9) +
    ggplot2::facet_wrap(~Region, scales = 'free') +
    ggplot2::labs(x = 'Year', y = 'Spawning Stock Biomass', color = 'Model', fill = 'Model') +
    ggplot2::coord_cartesian(ylim = c(0,NA)) +
    theme_sablefish()

  # if using confidence intervals
  if(do_ci) {
    comb_ts_plot <- comb_ts_plot + ggplot2::geom_ribbon(alpha = 0.3, color = NA)
    rec_ts_plot <- rec_ts_plot + ggplot2::geom_ribbon(alpha = 0.3, color = NA)
    ssb_ts_plot <- ssb_ts_plot + ggplot2::geom_ribbon(alpha = 0.3, color = NA)
    total_biom_plot <- total_biom_plot + ggplot2::geom_ribbon(alpha = 0.3, color = NA)
    ssb0_plot <- ssb0_plot + ggplot2::geom_ribbon(alpha = 0.3, color = NA)
  }

  return(list(comb_ts_plot, f_ts_plot, rec_ts_plot, ssb_ts_plot, total_biom_plot, ssb0_plot, ssb_ssb0_plot))
}

#' Get Fishery and Survey Selectivity Plots
#'
#' @param rep List of n_models of `SPoRC` report lists
#' @param model_names Vector of model names
#' @param Selex_Type Character vector specifying whether to output age or length-based selectivity (age, length)
#' @param year_indx Year index for which selectivity year to plot (can be an integer or vector)
#'
#' @returns Plots of terminal year fishery and survey selectivity by fleet, region, and sex across models
#' @export get_selex_plot
#' @family Plotting
#' @examples
#' \dontrun{
#' get_selex_plot(list(rep1, rep2), c("Model1", "Model2"), year_indx = c(1:30))
#' }
get_selex_plot <- function(rep, model_names, Selex_Type = 'age', year_indx = NULL) {

  fishsel_plot_df <- data.frame()
  srvsel_plot_df <- data.frame()

  for (i in seq_along(rep)) {
    # Get appropriate selectivity components based on type
    fish_sel <- switch(Selex_Type,
                       age = rep[[i]]$fish_sel,
                       length = rep[[i]]$fish_sel_l)

    srv_sel <- switch(Selex_Type,
                      age = rep[[i]]$srv_sel,
                      length = rep[[i]]$srv_sel_l)

    # Define helper for reshaping and annotating
    reshape_and_annotate <- function(df, model) {
      reshape2::melt(df) %>%
        dplyr::rename(Region = Var1, Year = Var2, Bin = Var3, Sex = Var4, Fleet = Var5) %>%
        dplyr::mutate(Region = paste("Region", Region),
                      Fleet  = paste("Fleet", Fleet),
                      Sex    = paste("Sex", Sex),
                      Model  = model)
    }

    fishsel_plot_df <- rbind(fishsel_plot_df, reshape_and_annotate(fish_sel, model_names[i]))
    srvsel_plot_df  <- rbind(srvsel_plot_df, reshape_and_annotate(srv_sel,  model_names[i]))
  }

  if(is.null(year_indx)) year_indx <- max(fishsel_plot_df$Year)

  # fishery selectivity plot
  fish_sel_plot <- ggplot2::ggplot(
    fishsel_plot_df %>%
      dplyr::filter(Year %in% c(year_indx)),
    ggplot2::aes(x = Bin, y = value, group = interaction(Year, Model), color = Year)
  ) +
    ggplot2::geom_line(aes(linetype = factor(Model)), linewidth = 1.3, alpha = 0.8) +
    ggplot2::facet_grid(Sex ~ Region + Fleet) +
    ggplot2::scale_color_viridis_c() +
    ggplot2::labs(x = 'Bin', y = 'Fishery selectivity', linetype = 'Model', color = 'Year') +
    ggplot2::coord_cartesian(ylim = c(0, NA)) +
    theme_sablefish() +
    ggplot2::theme(legend.key.width = unit(2, "lines"))

  # survey selectivity plot
  srv_sel_plot <- ggplot2::ggplot(
    srvsel_plot_df %>%
      dplyr::filter(Year %in% c(year_indx)),
    ggplot2::aes(x = Bin, y = value, group = interaction(Year, Model), color = Year)
  ) +
    ggplot2::geom_line(aes(linetype = factor(Model)), linewidth = 1.3, alpha = 0.8) +
    ggplot2::facet_grid(Sex ~ Region + Fleet) +
    ggplot2::scale_color_viridis_c() +
    ggplot2::labs(x = 'Bin', y = 'Survey selectivity', linetype = 'Model', color = 'Year') +
    ggplot2::coord_cartesian(ylim = c(0, NA)) +
    theme_sablefish() +
    ggplot2::theme(legend.key.width = unit(2, "lines"))

  return(list(fish_sel_plot, srv_sel_plot))
}

#' Get Plots of Biological Quantities
#'
#' @param data List of n_models of `SPoRC` data lists
#' @param rep List of n_models of `SPoRC` report lists
#' @param model_names Vector of model names
#'
#' @returns A list of plots for terminal year and terminal season movement, natural mortality, weight-at-age, and maturity at age across models
#' @export get_biological_plot
#' @family Plotting
#' @examples
#' \dontrun{
#' get_biological_plot(list(data1, data2), list(rep1, rep2), c("Model1", "Model2"))
#' }
get_biological_plot <- function(data,
                                rep,
                                model_names) {

  move_plot_df <- data.frame()
  natmort_plot_df <- data.frame()
  waa_plot_df <- data.frame()
  mataa_plot_df <- data.frame()

  for(i in 1:length(rep)) {

    # Movement
    move_plot_tmp_df <- reshape2::melt(rep[[i]]$Movement) %>%
      dplyr::rename(Region_From = from, Region_To = to, Seas = seas, Year = years, Age = ages, Sex = sexes) %>%
      {if (data[[i]]$do_recruits_move == 0) filter(., Age != min(data[[i]]$ages)) else .} %>%
      dplyr::mutate(Region_From = paste("From Region", Region_From),
                    Region_To = paste("To Region", Region_To),
                    Sex = paste("Sex", Sex),
                    Model = model_names[i]
      )

    # Natural Mortality
    natmort_plot_tmp_df <- reshape2::melt(rep[[i]]$natmort) %>%
      dplyr::rename(Region = Var1, Year = Var2, Age = Var3, Sex = Var4) %>%
      dplyr::mutate(Region = paste("Region", Region),
                    Sex = paste("Sex", Sex),
                    Model = model_names[i]
      )

    # Spawning Weight-at-age
    waa_plot_tmp_df <- reshape2::melt(data[[i]]$WAA) %>%
      dplyr::rename(Region = Var1, Year = Var2, Seas = Var3, Age = Var4, Sex = Var5) %>%
      dplyr::mutate(Region = paste("Region", Region),
                    Sex = paste("Sex", Sex),
                    Model = model_names[i]
      )

    # Maturity at age
    mataa_plot_tmp_df <- reshape2::melt(data[[i]]$MatAA) %>%
      dplyr::rename(Region = Var1, Year = Var2, Seas = Var3, Age = Var4, Sex = Var5) %>%
      dplyr::mutate(Region = paste("Region", Region),
                    Sex = paste("Sex", Sex),
                    Model = model_names[i]
      )

    # bind all
    move_plot_df <- rbind(move_plot_df, move_plot_tmp_df)
    natmort_plot_df <- rbind(natmort_plot_df, natmort_plot_tmp_df)
    waa_plot_df <- rbind(waa_plot_df, waa_plot_tmp_df)
    mataa_plot_df <- rbind(mataa_plot_df, mataa_plot_tmp_df)
  }

  # Movement plot
  move_plot <- ggplot(move_plot_df %>% dplyr::filter(Year == max(move_plot_df$Year), Seas == max(move_plot_df$Seas)),
                      ggplot2::aes(x = Age, y = value, color = factor(Model), lty = Sex)) +
    ggplot2::geom_line(lwd = 1) +
    ggplot2::facet_grid(Region_To~Region_From) +
    ggplot2::labs(x = 'Age', y = 'Movement Probabilities', color = 'Model', lty = 'Sex') +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    theme_sablefish() +
    ggplot2::theme(legend.key.width = unit(2, "lines"))

  # Natural mortality plot
  natmort_plot <- ggplot2::ggplot(natmort_plot_df %>%
                                    dplyr::filter(Year == max(natmort_plot_df$Year), Seas == max(natmort_plot$Seas)),
                                  ggplot2::aes(x = Age, y = value, color = factor(Model))) +
    ggplot2::geom_line(lwd = 2) +
    ggplot2::facet_grid(Region~Sex) +
    ggplot2::labs(x = 'Age', y = 'Natural Mortality', color = 'Model') +
    ggplot2::coord_cartesian(ylim = c(0, NA)) +
    theme_sablefish() +
    ggplot2::theme(legend.key.width = unit(2, "lines"))

  # Weight at age plot
  waa_plot <- ggplot2::ggplot(waa_plot_df %>%
                                dplyr::filter(Year == max(waa_plot_df$Year), Seas == max(waa_plot$Seas)),
                              ggplot2::aes(x = Age, y = value, color = factor(Model))) +
    ggplot2::geom_line(lwd = 2) +
    ggplot2::facet_grid(Region~Sex) +
    ggplot2::labs(x = 'Age', y = 'Spawning Weight at Age', color = 'Model') +
    ggplot2::coord_cartesian(ylim = c(0, NA)) +
    theme_sablefish() +
    ggplot2::theme(legend.key.width = unit(2, "lines"))

  # Maturity plot
  mataa_plot <- ggplot2::ggplot(mataa_plot_df %>%
                                  dplyr::filter(Year == max(mataa_plot_df$Year), Seas == max(mataa_plot$Seas)),
                                ggplot2::aes(x = Age, y = value, color = factor(Model))) +
    ggplot2::geom_line(lwd = 2) +
    ggplot2::facet_grid(Region~Sex) +
    ggplot2::labs(x = 'Age', y = 'Maturity at Age', color = 'Model') +
    ggplot2::coord_cartesian(ylim = c(0, NA)) +
    theme_sablefish() +
    ggplot2::theme(legend.key.width = unit(2, "lines"))

  return(list(move_plot, natmort_plot, waa_plot, mataa_plot))

}

#' Get Data Fitted to Plot
#'
#' @param data List of n_models of `SPoRC` data lists
#' @param model_names Character vector of model names
#'
#' @returns A plot of data that were fitted to across models
#' @export get_data_fitted_plot
#' @family Plotting
#' @examples
#' \dontrun{
#' get_data_fitted_plot(list(data1, data2), c("Model1", "Model2"))
#' }
get_data_fitted_plot <- function(data,
                                 model_names
                                 ) {

  data_plot_all_df <- data.frame()
  for(i in 1:length(data)) {

    # Get tag release indicator
    if(data[[i]]$UseTagging == 1) {
      use_tag_indicator <- array(0, dim = c(max(data[[i]]$tag_release_indicator[,1]), max(data[[i]]$tag_release_indicator[,2]), max(data[[i]]$tag_release_indicator[,3])))
      use_tag_indicator[data[[i]]$tag_release_indicator[,1],data[[i]]$tag_release_indicator[,2],data[[i]]$tag_release_indicator[,3]] <- 1
    }

    # Bind all data indicators together
    data_plot_df <- reshape2::melt(data[[i]]$UseSrvLenComps) %>% # Survey lengths
      dplyr::rename(Region = Var1, Year = Var2, Seas = Var3, Fleet = Var4) %>%
      dplyr::mutate(Type = paste('Survey Lengths', "Seas", Seas, "Fleet", Fleet)) %>%

      dplyr::bind_rows(
        # survey ages
        reshape2::melt(data[[i]]$UseSrvAgeComps) %>%
          dplyr::rename(Region = Var1, Year = Var2, Seas = Var3, Fleet = Var4) %>%
          dplyr::mutate(Type = paste('Survey Ages', "Seas", Seas, "Fleet", Fleet)),

        # fishery lengths
        reshape2::melt(data[[i]]$UseFishLenComps) %>%
          dplyr::rename(Region = Var1, Year = Var2, Seas = Var3, Fleet = Var4) %>%
          dplyr::mutate(Type = paste('Fishery Lengths', "Seas", Seas, "Fleet", Fleet)),

        # fishery ages
        reshape2::melt(data[[i]]$UseFishAgeComps) %>%
          dplyr::rename(Region = Var1, Year = Var2, Seas = Var3, Fleet = Var4) %>%
          dplyr::mutate(Type = paste('Fishery Ages', "Seas", Seas, "Fleet", Fleet)),

        # fishery catches
        reshape2::melt(data[[i]]$UseCatch) %>%
          dplyr::rename(Region = Var1, Year = Var2, Seas = Var3, Fleet = Var4) %>%
          dplyr::mutate(Type = paste('Fishery Catch', "Seas", Seas, "Fleet", Fleet)),

        # fishery indices
        reshape2::melt(data[[i]]$UseFishIdx) %>%
          dplyr::rename(Region = Var1, Year = Var2, Seas = Var3, Fleet = Var4) %>%
          dplyr::mutate(Type = paste('Fishery Index', "Seas", Seas, "Fleet", Fleet)),

        # survey indices
        reshape2::melt(data[[i]]$UseSrvIdx) %>%
          dplyr::rename(Region = Var1, Year = Var2, Seas = Var3, Fleet = Var4) %>%
          dplyr::mutate(Type = paste('Survey Index', "Seas", Seas, "Fleet", Fleet))
      )

    # Add tagging if used
    if (data[[i]]$UseTagging == 1) {
      tag_df <- reshape2::melt(use_tag_indicator) %>%
        dplyr::rename(Region = Var1, Year = Var2, Seas = Var3) %>%
        dplyr::mutate(Type = paste('Tagging', "Seas", Seas), Fleet = NA)
      data_plot_df <- dplyr::bind_rows(data_plot_df, tag_df)
    }

    # remove data not fitted to
    data_plot_df <- data_plot_df %>%
      dplyr::filter(value != 0) %>%
      dplyr::mutate(Region = paste("Region", Region),
                    Model = model_names[i])

    data_plot_all_df <- rbind(data_plot_all_df, data_plot_df)
  }

  data_plot <- ggplot2::ggplot(data_plot_all_df,
                               ggplot2::aes(x = Year, y = Type, fill = Type)) +
    ggplot2::geom_point(size = 3, pch = 21, color = 'black', alpha = 0.8) +
    ggplot2::facet_grid(Model~Region) +
    theme_sablefish() +
    ggplot2::theme(legend.position = 'none') +
    ggplot2::labs(x = 'Year', y = '')

  return(data_plot)
}

#' Get plot of negative log likelihood values
#'
#' @param rep List of n_models of `SPoRC` report lists
#' @param model_names Vector of model names
#' @param data List of data from `SPoRC`
#'
#' @returns Plot and tables of negative log likelihood values across models
#' @export get_nLL_plot
#' @family Model Diagnostics
#' @examples
#' \dontrun{
#' get_nLL_plot(list(data1, data2), list(rep1, rep2), c("Model1", "Model2"))
#' }
get_nLL_plot <- function(data,
                         rep,
                         model_names
                         ) {


  nLL_all_df <- data.frame() # empty dataframe
  for(i in 1:length(rep)) {

    # Negative log likelihoods
    nLL_df <- data.frame(

      # nLL values
      value = c(safe_extract(rep[[i]], "jnLL"),
                safe_extract(rep[[i]], "h_nLL"),
                safe_extract(rep[[i]], "M_nLL"),
                safe_extract(rep[[i]], "Rec_prop_nLL"),
                sum(data[[i]]$Wt_Rec * safe_extract(rep[[i]], "Rec_nLL")),
                safe_extract(rep[[i]], "sel_nLL"),
                sum(data[[i]]$Wt_Tagging * safe_extract(rep[[i]], "Tag_nLL")),
                sum(data[[i]]$Wt_Catch * safe_extract(rep[[i]], "Catch_nLL")),
                sum(data[[i]]$Wt_F * safe_extract(rep[[i]], "Fmort_nLL")),
                safe_extract(rep[[i]], "srv_q_nLL"),
                safe_extract(rep[[i]], "fish_q_nLL"),
                sum(data[[i]]$Wt_SrvIdx * safe_extract(rep[[i]], "SrvIdx_nLL")),
                safe_extract(rep[[i]], "TagRep_nLL"),
                sum(data[[i]]$Wt_FishIdx * safe_extract(rep[[i]], "FishIdx_nLL")),
                sum(data[[i]]$Wt_Rec * safe_extract(rep[[i]], "Init_Rec_nLL")),
                safe_extract(rep[[i]], "Movement_nLL"),
                sum(safe_extract(rep[[i]], "SrvAgeComps_nLL")),
                sum(safe_extract(rep[[i]], "FishAgeComps_nLL")),
                sum(safe_extract(rep[[i]], "SrvLenComps_nLL")),
                sum(safe_extract(rep[[i]], "FishLenComps_nLL"))),

      # nLL names
      name = c("jnLL", "Steepness Prior", "M Prior", "Recruitment Prop Prior", "Recruitment Penalty",
               "Selectivity Penalty", "Tag nLL", "Catch nLL", "Fishing Mortality Penalty",
               "Survey Q Prior", "Fishery Q Prior", "Survey Index nLL", "Tag Reporting Prior",
               "Fishery Index nLL", "Initial Age Penalty", "Movement Prior",
               "Survey Age nLL", "Fishery Age nLL", "Survey Length nLL", "Fishery Length nLL"),

      # nLL types
      type = c('jnLL', 'Prior', "Prior", "Prior", "Penalty", "Penalty", "Tagging",
               "Catch", "Penalty", "Prior", "Prior", "Index",
               "Prior", "Index", "Penalty", "Prior", "Age",
               "Age", "Length", "Length"),

      Model = model_names[i]
    )

    nLL_all_df <- rbind(nLL_all_df, nLL_df)
  }

  nLL_plot <- ggplot2::ggplot(nLL_all_df %>% dplyr::filter(value != 0),
                              ggplot2::aes(x = name, y = value, fill = type)) +
    ggplot2::geom_col() +
    ggplot2::facet_wrap(~Model) +
    theme_sablefish() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 0.5)) +
    ggplot2::labs(x = 'Likelihood Component', y = 'Likelihood', fill = 'Type')


  # get nLL table
  nLL_table <- nLL_all_df %>%
    dplyr::select(-type) %>%
    dplyr::filter(value != 0) %>%  # remove zeros
    tidyr::pivot_wider(names_from = name, values_from = value)
  table_plot <- grid::grid.grabExpr(gridExtra::grid.table(nLL_table))
  table_plot1 <- cowplot::ggdraw() + cowplot::draw_grob(table_plot)

  return(list(nLL_plot, table_plot1))
}

#' Get Index Fits Plot
#'
#' @param data List of n_models of `SPoRC` data lists
#' @param rep List of n_models of `SPoRC` report lists
#' @param model_names Vector of model names
#'
#' @returns A plot of fitted values to various indices across models
#' @export get_idx_fits_plot
#' @family Model Diagnostics
#' @examples
#' \dontrun{
#' get_idx_fits_plot(list(data1, data2), list(rep1, rep2), c("Model1", "Model2"))
#' }
get_idx_fits_plot <- function(data,
                              rep,
                              model_names
                              ) {

  idx_fits_all <- data.frame()
  # get index fits data
  for(i in 1:length(rep)) {
    idx_fits <- get_idx_fits(data = data[[i]], rep = rep[[i]], year_labs = data[[i]]$years) %>%
      dplyr::mutate(Model = model_names[i])
    idx_fits_all <- rbind(idx_fits_all, idx_fits) # bind
  }

  # Plot index fits
  idx_fit_plot <- ggplot2::ggplot() +
    ggplot2::geom_line(idx_fits_all %>% dplyr::filter(obs != 0),
                       mapping = ggplot2::aes(x = Year, y = value, color = factor(Model)), lwd = 1.3) +
    ggplot2::geom_pointrange(idx_fits_all %>% dplyr::filter(obs != 0),
                             mapping = ggplot2::aes(x = Year, y = obs, ymin = lci, ymax = uci), color = 'black') +
    ggplot2::labs(x = "Year", y = 'Index', color = 'Model') +
    theme_sablefish() +
    ggplot2::coord_cartesian(ylim = c(0,NA)) +
    ggplot2::facet_grid(Category~Region, scales = 'free_y')

  return(idx_fit_plot)
}

#' Title Get Catch Fits Plot
#'
#' @param data List of n_models of `SPoRC` data lists
#' @param rep List of n_models of `SPoRC` report lists
#' @param model_names Vector of model names
#'
#' @returns A plot of fitted values to various catch time series across models
#' @export get_catch_fits_plot
#' @family Model Diagnostics
#' @examples
#' \dontrun{
#' get_catch_fits_plot(list(data1, data2), list(rep1, rep2), c("Model1", "Model2"))
#' }
get_catch_fits_plot <- function(data,
                                rep,
                                model_names
                                ) {

  catch_fits_all <- data.frame()
  # get catch fits data
  for(i in 1:length(rep)) {
    # Get catch fits
    catch_fits <- reshape2::melt(rep[[i]]$PredCatch) %>%
      dplyr::rename(Region = Var1, Year = Var2, Seas = Var3, Fleet = Var4) %>%
      dplyr::left_join(
        reshape2::melt(data[[i]]$ObsCatch) %>%
          dplyr::left_join(reshape2::melt(exp(rep[[i]]$ln_sigmaC) / data[[i]]$Wt_Catch) %>%
                             dplyr::rename(se = value),
                           by = c("Var1", "Var2", "Var3", "Var4")) %>%
          dplyr::rename(Region = Var1, Year = Var2, Seas = Var3, Fleet = Var4, obs = value),
        by = c("Region", "Year", "Seas", "Fleet")
      ) %>%
      dplyr::mutate(Model = model_names[i],
                    Region = paste('Region', Region),
                    Seas = paste("Seas", Seas),
                    Fleet = paste('Fleet', Fleet),
                    Seas_Fleet = paste(Seas, Fleet))
    catch_fits_all <- rbind(catch_fits_all, catch_fits) # bind
  }

  # Plot catch fits
  catch_fit_plot <- ggplot2::ggplot() +
    ggplot2::geom_line(catch_fits_all %>% dplyr::filter(obs != 0),
                       mapping = ggplot2::aes(x = Year, y = value, color = factor(Model)), lwd = 1.3) +
    ggplot2::geom_pointrange(catch_fits_all %>% dplyr::filter(obs != 0),
                             mapping = ggplot2::aes(x = Year, y = obs, ymin = exp(log(obs) - 1.96 * se),
                                                    ymax = exp(log(obs) + 1.96 * se)), color = 'black') +
    ggplot2::labs(x = "Year", y = 'Catch', color = 'Model') +
    theme_sablefish() +
    ggplot2::coord_cartesian(ylim = c(0,NA)) +
    ggplot2::facet_grid(Seas_Fleet~Region, scales = 'free_y')

  return(catch_fit_plot)
}



#' Get Retrospective Plot
#'
#' @param retro_output Dataframe generated from do_retrospective
#' @param Rec_Age Age in which recruitment occurs
#'
#' @returns A retrospective plot of recruitment and SSB in relative and absolute scales, as well as a retrospective plot of recruitment by cohort (squid plot)
#' @export get_retrospective_plot
#' @family Model Diagnostics
#' @examples
#' \dontrun{
#' # do retrospective
#' retro <- do_retrospective(n_retro = 7, # number of retro peels to run
#' data = data, # rtmb data
#' parameters = parameters, # rtmb parameters
#' mapping = mapping, # rtmb mapping
#' random = NULL, # if random effects are used
#' do_par = TRUE, # whether or not to parralleize
#' n_cores = 7, # if parallel, number of cores to use
#' do_francis = F, # if we want tod o Francis
#' n_francis_iter = NULL # Number of francis iterations to do
#' )
#' get_retrospective_plot(retro, Rec_Age = 2)
#' }
get_retrospective_plot <- function(retro_output, Rec_Age) {

  # Get relative differences
  ret_df <- get_retrospective_relative_difference(retro_output)

  # Get mohns rho (mean relative difference for a given terminal year peel to terminal year estimates)
  mohns_rho <- ret_df %>%
    dplyr::filter(peel == (max(Year) - Year)) %>%
    dplyr::group_by(Type, Region) %>%
    dplyr::summarize(rho = mean(rd))

  # get retrospective plot
  retro_plot <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, lty = 2, lwd = 1.3) +
    ggplot2::geom_line(ret_df,
                       mapping = ggplot2::aes(x = Year, y = rd, group = as.numeric(peel), color = as.numeric(peel)), lwd = 1.3) +
    ggplot2::geom_point(ret_df %>% dplyr::filter(peel == max(Year) - Year),
                        mapping = ggplot2::aes(x = Year, y = rd, group = as.numeric(peel), fill = as.numeric(peel)),
                        pch = 21, size = 6) +
    ggplot2::geom_text(mohns_rho, mapping = aes(x = -Inf, y = Inf, label = paste("Mohns Rho:", round(rho, 4))),
                       hjust = -0.3, vjust = 3, size = 5) +
    ggplot2::guides(color = ggplot2::guide_colourbar(barwidth = 15, barheight = 1.3)) +
    ggplot2::labs(x = 'Year', y = 'Relative Difference from Terminal Year', color = 'Retrospective Year', fill = 'Retrospective Year') +
    ggplot2::scale_color_viridis_c() +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::facet_grid(Region~Type, scales = 'free') +
    theme_sablefish() +
    ggplot2::theme(legend.position = 'top')

  # get absolute retro plot
  abs_retro_plot <- ggplot2::ggplot() +
    ggplot2::geom_line(retro_output %>% filter(peel != 0), mapping = ggplot2::aes(x = Year, y = value, group = peel, color = peel), lwd = 1) +
    ggplot2::geom_line(retro_output %>% filter(peel == 0), mapping = ggplot2::aes(x = Year, y = value), lty = 2, lwd = 1) +
    ggplot2::scale_color_viridis_c() +
    ggplot2::coord_cartesian(ylim = c(0, NA)) +
    ggplot2::facet_grid(Type~Region, scales = 'free_y') +
    theme_sablefish() +
    ggplot2::labs(x = 'Year', y = 'Value', color = 'Peel')

  # get squid plot
  squid_plot <- retro_output %>%
    dplyr::mutate(Year = Year, terminal = max(retro_output$Year) - peel, cohort = Year - Rec_Age, years_est = terminal-Year) %>%
    dplyr::filter(Type == 'Recruitment', cohort %in% seq(max(retro_output$Year) - 10, max(retro_output$Year), 1), terminal != Year) %>%
    ggplot2::ggplot(ggplot2::aes(x = years_est - 1, y = value, group = Year, color = factor(cohort))) +
    ggplot2::geom_line(lwd = 1.3) +
    ggplot2::geom_point(size = 4) +
    ggplot2::facet_wrap(~Region, scales = 'free') +
    ggplot2::theme_bw(base_size = 15) +
    ggplot2::labs(x = 'Years since cohort was last estimated', y = 'Recruitment (millions)', color = 'Cohort')

  return(list(retro_plot, abs_retro_plot, squid_plot))
}
#' Plotting function for all basic quantities
#'
#' @param data List of n_models of `SPoRC` data lists
#' @param rep List of n_models of `SPoRC` report lists
#' @param sd_rep List of n_models of sd report lists from `SPoRC`
#' @param out_path Path to the output directory. Users only need to specify the path.
#' @param model_names Character vector of model names
#'
#' @returns A series of plots compared across models outputted as a pdf in the specified directory
#' @export plot_all_basic
#' @family Plotting
#'
#' @examples
#' \dontrun{
#' plot_all_basic(
#'   data = list(data1, data2),
#'   rep = list(rep1, rep2),
#'   sd_rep = list(sd_rep1, sd_rep2),
#'   model_names = c("Model1", "Model2"),
#'   out_path = here::here()
#' )
#' }
plot_all_basic <- function(data,
                           rep,
                           sd_rep,
                           model_names,
                           out_path) {

  pdf(here::here(out_path, "plot_results.pdf"), width = 25, height = 13)
  print(get_biological_plot(data = data, rep = rep, model_names = model_names))
  print(get_data_fitted_plot(data = data, model_names = model_names))
  print(get_ts_plot(rep = rep, sd_rep = sd_rep, model_names = model_names))
  print(get_selex_plot(rep = rep, model_names = model_names))
  print(get_nLL_plot(data = data, rep = rep, model_names = model_names))
  dev.off()

}
#' Generate Key Projection Quantities and Table Plot
#'
#' Calculates biological and fishery reference points and performs population projections to estimate terminal spawning biomass, catch advice, and reference point values by model and region. Also returns a formatted table plot of key quantities.
#'
#' @param data A list of model input data objects, one for each model (i.e., a list of SPoRC-formatted data lists). Each element should contain information on regions, years, ages, fleets, and biological inputs (e.g., weight-at-age, maturity, mortality).
#' @param rep A list of model output objects, one for each model (i.e., a list of SPoRC-formatted report lists). Each element must include recruitment, selectivity, mortality, and numbers-at-age.
#' @param reference_points_opt A named list specifying options for reference point calculations. See \code{\link{Get_Reference_Points}} for more details. Must include:
#' \describe{
#'   \item{SPR_x}{Spawning potential ratio (e.g., 0.4) for calculating F reference points. May be \code{NULL} if using \code{bh_rec}.}
#'   \item{t_spwn}{Fraction of year when spawning occurs (e.g., 0.5).}
#'   \item{sex_ratio_f}{Proportion of recruits that are female.}
#'   \item{calc_rec_st_yr}{Start year for averaging recruitment.}
#'   \item{rec_age}{Recruitment age.}
#'   \item{type}{Reference point calculation method (e.g., "multi_region").}
#'   \item{what}{Type of output requested from the reference point function.}
#' }
#' @param proj_model_opt A named list of projection settings. See \code{\link{Do_Population_Projection}} for details. Must include:
#' \describe{
#'   \item{n_proj_yrs}{Number of years to project forward.}
#'   \item{HCR_function}{Harvest control rule function to use.}
#'   \item{recruitment_opt}{Recruitment assumption (e.g., "mean_rec", "bh_rec", "inv_gauss").}
#'   \item{fmort_opt}{Fishing mortality assumption (e.g., "input", "HCR").}
#'   \item{n_avg_yrs}{Number of years to average over for projection inputs.}
#' }
#' @param model_names A character vector of model identifiers (e.g., c("Base", "Alt1", "Alt2")), one for each element in \code{data} and \code{rep}.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{\code{[[1]]}}{A data.frame of key quantities by model and region, including terminal SSB, catch advice, and reference points.}
#'   \item{\code{[[2]]}}{A cowplot table plot (ggdraw object) of the same key quantities.}
#' }
#'
#' @details
#' This function checks input list completeness, calculates reference points using \code{Get_Reference_Points()}, performs population projections with \code{Do_Population_Projection()}, and assembles both tabular and visual summaries.
#'
#' If \code{recruitment_opt} is set to "inv_gauss", a warning is issued since only a single simulation will be run. This is typically not appropriate and an alternative assumption is recommended.
#'
#' @seealso \code{\link{Get_Reference_Points}}, \code{\link{Do_Population_Projection}}
#'
#' @examples
#' \dontrun{
#' reference_points_opt <- list(SPR_x = 0.4,
#'                              t_spwn = 0,
#'                              sex_ratio_f = 0.5,
#'                              calc_rec_st_yr = 20,
#'                              rec_age = 2,
#'                              type = "multi_region",
#'                              what = "global_SPR")
#'
#' proj_model_opt <- list(
#'   n_proj_yrs = 2,
#'   n_avg_yrs = 1,
#'   HCR_function = function(x, frp, brp, alpha = 0.05) {
#'     stock_status <- x / brp
#'     if (stock_status >= 1) f <- frp
#'     if (stock_status > alpha && stock_status < 1) f <- frp * (stock_status - alpha) / (1 - alpha)
#'     if (stock_status < alpha) f <- 0
#'     return(f)
#'   },
#'   recruitment_opt = "mean_rec",
#'   fmort_opt = "HCR"
#' )
#'
#' out <- get_key_quants(list(mlt_rg_sable_data),
#'                       list(mlt_rg_sable_rep),
#'                       reference_points_opt,
#'                       proj_model_opt,
#'                       "Model 1")
#' out[[1]]  # key quantities data.frame
#' out[[2]]  # table plot
#' }
#'
#' @export get_key_quants
#' @family Reference Points and Projections
get_key_quants <- function(data,
                           rep,
                           reference_points_opt,
                           proj_model_opt,
                           model_names
) {

  # required elements for reference points opt
  required <- c("SPR_x", "t_spwn", "sex_ratio_f", "calc_rec_st_yr", "rec_age", "type", "what")
  missing <- setdiff(required, names(reference_points_opt))

  # check to see if reference points opt has all of the required elements
  if (length(missing) > 0) {
    stop("`reference_points_opt` is missing the following required elements: ",
         paste(missing, collapse = ", "),
         ". It should include all of the following: ",
         paste(required, collapse = ", "), ".")
  }

  # required elements for catch projections opt
  required <- c("n_proj_yrs", "HCR_function", "recruitment_opt", "fmort_opt", "n_avg_yrs")
  missing <- setdiff(required, names(proj_model_opt))

  # check to see if reference points opt has all of the required elements
  if (length(missing) > 0) {
    stop("`proj_model_opt` is missing the following required elements: ",
         paste(missing, collapse = ", "),
         ". It should include all of the following: ",
         paste(required, collapse = ", "), ".")
  }

  if(proj_model_opt$recruitment_opt == 'inv_gauss') {
    warning("Recruitment during the projection period is set to 'inv_gauss', but only a single simulation will be run. This is likely inappropriate. Consider using an alternative recruitment option such as zero_rec, mean_rec, or bh_rec.")
  }

  # containers
  ref_pts <- list()
  key_quants_df <- data.frame()

  for(i in 1:length(rep)) {

    # get reference points
    tmp_ref_pts <- Get_Reference_Points(data = data[[i]],
                                        rep = rep[[i]],
                                        SPR_x = reference_points_opt$SPR_x,
                                        t_spwn = reference_points_opt$t_spwn,
                                        sex_ratio_f = reference_points_opt$sex_ratio_f,
                                        calc_rec_st_yr = reference_points_opt$calc_rec_st_yr,
                                        rec_age = reference_points_opt$rec_age,
                                        type = reference_points_opt$type,
                                        what = reference_points_opt$what
    )

    # input into list
    ref_pts[[i]] <- tmp_ref_pts

    # do population project to get catch advice
    n_proj_yrs <- proj_model_opt$n_proj_yrs # number of projection years
    t_spawn <- reference_points_opt$t_spwn # spawn timing

    # terminal estimates
    terminal_NAA <-  array(rep[[i]]$NAA[,length(data[[i]]$years),,], dim = c(data[[i]]$n_regions, length(data[[i]]$ages), data[[i]]$n_sexes)) # terminal NAA
    terminal_NAA0 <-  array(rep[[i]]$NAA0[,length(data[[i]]$years),,], dim = c(data[[i]]$n_regions, length(data[[i]]$ages), data[[i]]$n_sexes)) # terminal NAA
    terminal_F <- array(rep[[i]]$Fmort[,length(data[[i]]$years),], dim = c(data[[i]]$n_regions, data[[i]]$n_fish_fleets)) # terminal F
    recruitment <- array(rep[[i]]$Rec[,reference_points_opt$calc_rec_st_yr:(length(data[[i]]$years) - reference_points_opt$rec_age)],
                         dim = c(data[[i]]$n_regions, length(reference_points_opt$calc_rec_st_yr:(length(data[[i]]$years) - reference_points_opt$rec_age)))) # recruitment

    # demographics
    # determine years to average over demogrphaics
    n_avg_yrs <- proj_model_opt$n_avg_yrs
    n_proj_yrs <- proj_model_opt$n_proj_yrs
    n_yrs <- length(data[[i]]$years)
    avg_yrs <- (n_yrs - n_avg_yrs + 1):n_yrs

    # spawning weight-at-age
    WAA_avg <- apply(data[[i]]$WAA[,avg_yrs,,,drop = FALSE], c(1, 3, 4), mean)
    WAA <- array(rep(WAA_avg, each = n_proj_yrs), dim = c(data[[i]]$n_regions, n_proj_yrs, length(data[[i]]$ages), data[[i]]$n_sexes))

    WAA_fish_avg <- apply(data[[i]]$WAA_fish[,avg_yrs,,,,drop = FALSE], c(1, 3, 4, 5), mean)
    WAA_fish <- array(rep(WAA_fish_avg, each = n_proj_yrs), dim = c(data[[i]]$n_regions, n_proj_yrs, length(data[[i]]$ages), data[[i]]$n_sexes, data[[i]]$n_fish_fleets))

    # maturity at age
    MatAA_avg <- apply(data[[i]]$MatAA[,avg_yrs,,,drop = FALSE], c(1, 3, 4), mean)
    MatAA <- array(rep(MatAA_avg, each = proj_model_opt$n_proj_yrs), dim = c(data[[i]]$n_regions, proj_model_opt$n_proj_yrs, length(data[[i]]$ages), data[[i]]$n_sexes))

    # natural mortality
    natmort_avg <- apply(rep[[i]]$natmort[,avg_yrs,,,drop = FALSE], c(1, 3, 4), mean)
    natmort <- array(rep(natmort_avg, each = proj_model_opt$n_proj_yrs), dim = c(data[[i]]$n_regions, proj_model_opt$n_proj_yrs, length(data[[i]]$ages), data[[i]]$n_sexes))

    # fishery selectivity
    fish_sel_avg <- apply(rep[[i]]$fish_sel[,avg_yrs,,,,drop = FALSE], c(1, 3, 4, 5), mean)
    fish_sel <- array(rep(fish_sel_avg, each = proj_model_opt$n_proj_yrs), dim = c(data[[i]]$n_regions, proj_model_opt$n_proj_yrs, length(data[[i]]$ages), data[[i]]$n_sexes, data[[i]]$n_fish_fleets))

    # movement
    Movement_avg <- apply(rep[[i]]$Movement[,,avg_yrs,,,drop = FALSE], c(1,2,4,5), mean) # movement
    Movement <- aperm(abind::abind(replicate(n_proj_yrs, Movement_avg, simplify = FALSE), along = 5), perm = c(1,2,5,3,4)) # movement

    # Sex ratio
    sexratio_avg <- apply(rep[[i]]$sexratio[,avg_yrs,,drop = FALSE], c(1,3), mean)
    sexratio <- array(dim = c(data[[i]]$n_regions, proj_model_opt$n_proj_yrs, data[[i]]$n_sexes)) # empty array
    for(yr in 1:proj_model_opt$n_proj_yrs) sexratio[, yr, ] <- sexratio_avg # populate empty array

    # Now, set up inputs for reference points
    f_ref_pt = array(tmp_ref_pts$f_ref_pt, dim = c(data[[i]]$n_regions, n_proj_yrs))
    b_ref_pt = array(tmp_ref_pts$b_ref_pt, dim = c(data[[i]]$n_regions, n_proj_yrs))

    # Set up beverton-holt options if using beverton holt for projection
    if(proj_model_opt$recruitment_opt == 'bh_rec') {
      bh_rec_opt <- list(
        recruitment_dd = data[[i]]$rec_dd,
        rec_lag = data[[i]]$rec_lag,
        R0 = rep[[i]]$R0,
        Rec_Prop = rep[[i]]$Rec_trans_prop,
        h = rep[[i]]$h_trans,
        SSB = rep[[i]]$SSB,
        # Using first year for demographics of computing unfished SSB
        WAA = data[[i]]$WAA[,1,,,drop = FALSE],
        MatAA = data[[i]]$MatAA[,1,,,drop = FALSE],
        natmort = rep[[i]]$natmort[,1,,,drop = FALSE]
      )
    } else {
      bh_rec_opt <- NULL
    }

    # do population projection
    out_proj <- Do_Population_Projection(n_proj_yrs = n_proj_yrs, # Number of projection years
                                         n_regions = data[[i]]$n_regions, # number of regions
                                         n_ages = length(data[[i]]$ages), # number of ages
                                         n_sexes = data[[i]]$n_sexes, # number of sexes
                                         sexratio = sexratio, # sex ratio for recruitment
                                         n_fish_fleets = data[[i]]$n_fish_fleets, # number of fishery fleets
                                         do_recruits_move = data[[i]]$do_recruits_move, # whether recruits move
                                         recruitment = recruitment, # recruitment values to use for mean recruitment
                                         terminal_NAA = terminal_NAA, # terminal numbers at age
                                         terminal_NAA0 = terminal_NAA0, # terminal unfished numbers at age
                                         terminal_F = terminal_F, # terminal F
                                         natmort = natmort, # natural mortality values to use in projection
                                         WAA = WAA, # spawning weight at age values to use in projection
                                         WAA_fish = WAA_fish, # fishery weight at age to use in projection
                                         MatAA = MatAA, # maturity at age values to use in projection
                                         fish_sel = fish_sel, # fishery selectivity values to use in projection
                                         Movement = Movement, # movement values
                                         f_ref_pt = f_ref_pt, # fishery reference points
                                         b_ref_pt = b_ref_pt, # biological reference points
                                         HCR_function = proj_model_opt$HCR_function, # HCR function
                                         recruitment_opt = proj_model_opt$recruitment_opt, # recruitment assumption
                                         fmort_opt = proj_model_opt$fmort_opt, # Fishing mortality in projection years (whether input or HCR)
                                         t_spawn = reference_points_opt$t_spwn, # Spawn timing
                                         bh_rec_opt = bh_rec_opt # beverton holt projection options
    )

    # extract out quantities and store
    key_quants_tmp <- data.frame(Model = model_names[i],
                                 Region = 1:data[[i]]$n_regions,
                                 Terminal_SSB = round(out_proj$proj_SSB[,1], 5),
                                 Terminal_SSB0 = round(out_proj$proj_Dynamic_SSB0[,1], 5),
                                 Terminal_F = rowSums(terminal_F),
                                 Catch_Advice = round(apply(out_proj$proj_Catch[,2,,drop = FALSE], c(1,2), sum), 5), # sum across fleets
                                 B_Ref_Pt = round(tmp_ref_pts$b_ref_pt, 5),
                                 F_Ref_Pt = round(tmp_ref_pts$f_ref_pt, 5),
                                 B_over_B_Ref = round(out_proj$proj_SSB[,1] / tmp_ref_pts$b_ref_pt, 5),
                                 B_over_DynB_Ref = round(out_proj$proj_SSB[,1] / out_proj$proj_Dynamic_SSB0[,1], 5),
                                 F_over_F_Ref = round(rowSums(terminal_F) / tmp_ref_pts$f_ref_pt, 5)
    )

    key_quants_df <- rbind(key_quants_df, key_quants_tmp)

  } # end i

  # output table
  table_plot <- grid::grid.grabExpr(gridExtra::grid.table(key_quants_df))
  table_plot1 <- cowplot::ggdraw() + cowplot::draw_grob(table_plot)

  return(list(key_quants_df, table_plot1))

} # end function
