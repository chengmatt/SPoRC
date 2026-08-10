# Stage 3 of 3: post fit
#
# Plots and summary tables from a fitted model. Reads through the extraction
# layer in diag_fits.R and diag_retrospective.R rather than digging into the
# fitted object directly. theme_sablefish is the shared ggplot2 theme.

#' Get Time Series Plots
#'
#' Generates a suite of time series plots for key population dynamics quantities
#' (spawning stock biomass, dynamic unfished SSB, total biomass, recruitment, and
#' fishing mortality) across one or more SPoRC model runs. Optionally overlays
#' approximate 95% confidence intervals derived from the delta method via the
#' sdreport.
#'
#' @param rep List of length \code{n_models}, where each element is a SPoRC
#'   report list (i.e. the output of \code{obj$report()} after optimisation).
#' @param sd_rep List of length \code{n_models}, where each element is a SPoRC
#'   sdreport list (i.e. the output of \code{sdreport(obj)}). Used to extract
#'   delta-method standard errors on log-scale quantities (\code{log_SSB},
#'   \code{log_Dynamic_SSB0}, \code{log_Total_Biom}, \code{log_Rec}).
#' @param model_names Character vector of length \code{n_models} giving display
#'   names for each model run. Used as legend labels across all plots.
#' @param do_ci Logical. If \code{TRUE} (default), approximate 95% confidence
#'   ribbons are added to SSB, total biomass, dynamic SSB0, and recruitment
#'   plots. Ribbons are computed on the natural scale as
#'   \code{exp(log(value) ± 1.96 * se)}.
#'
#' @return A named list of seven \code{ggplot} objects:
#' \describe{
#'   \item{[[1]] comb_ts_plot}{All quantities faceted by Type × Region.}
#'   \item{[[2]] f_ts_plot}{Fishing mortality by fleet and season, faceted by Region.}
#'   \item{[[3]] rec_ts_plot}{Recruitment by population, faceted by Region.}
#'   \item{[[4]] ssb_ts_plot}{Spawning stock biomass by population, faceted by Region.
#'     Dynamic SSB0 is excluded from this panel.}
#'   \item{[[5]] total_biom_plot}{Total biomass by population, faceted by Region.}
#'   \item{[[6]] ssb0_plot}{Dynamic unfished SSB (B0) by population, faceted by Region.}
#'   \item{[[7]] ssb_ssb0_plot}{SSB and dynamic SSB0 overlaid on the same panel,
#'     distinguished by linetype, faceted by Region. Useful for visualising depletion.}
#' }
#'
#' @export get_ts_plot
#' @family Plotting
#'
#' @examples
#' \dontrun{
#'   plots <- get_ts_plot(
#'     rep        = list(rep1, rep2),
#'     sd_rep     = list(sdrep1, sdrep2),
#'     model_names = c("Base", "Sensitivity"),
#'     do_ci      = TRUE
#'   )
#'   plots[[4]] # SSB time series
#'   plots[[7]] # SSB vs dynamic B0
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
      dplyr::rename(Pop = Var1, Region = Var2, Year = Var3) %>%
      dplyr::bind_cols(se = sd_rep[[i]]$sd[names(sd_rep[[i]]$value) == "log_SSB"]) %>%
      dplyr::mutate(lwr = exp(log(value) - 1.96 * se),
                    upr = exp(log(value) + 1.96 * se),
                    Region = paste("Region", Region),
                    Pop = paste("Population", Pop),
                    Type = paste(Pop, 'SSB'),
                    Model = model_names[i])

    # Dynamic Unfished Spawning Stock Biomass
    ssb0_plot_df <- reshape2::melt(rep[[i]]$Dynamic_SSB0) %>%
      dplyr::rename(Pop = Var1, Region = Var2, Year = Var3) %>%
      dplyr::bind_cols(se = sd_rep[[i]]$sd[names(sd_rep[[i]]$value) == "log_Dynamic_SSB0"]) %>%
      dplyr::mutate(lwr = exp(log(value) - 1.96 * se),
                    upr = exp(log(value) + 1.96 * se),
                    Region = paste("Region", Region),
                    Pop = paste("Population", Pop),
                    Type = paste(Pop, 'Dynamic SSB0'),
                    Model = model_names[i])

    # Total Biomass
    totbiom_plot_df <- reshape2::melt(rep[[i]]$Total_Biom) %>%
      dplyr::rename(Pop = Var1, Region = Var2, Year = Var3) %>%
      dplyr::bind_cols(se = sd_rep[[i]]$sd[names(sd_rep[[i]]$value) == "log_Total_Biom"]) %>%
      dplyr::mutate(lwr = exp(log(value) - 1.96 * se),
                    upr = exp(log(value) + 1.96 * se),
                    Region = paste("Region", Region),
                    Pop = paste("Population", Pop),
                    Type = paste(Pop, 'Total Biom'),
                    Model = model_names[i])

    # Recruitment
    rec_plot_df <- reshape2::melt(rep[[i]]$Rec) %>%
      dplyr::rename(Pop = Var1, Region = Var2, Year = Var3) %>%
      dplyr::bind_cols(se = sd_rep[[i]]$sd[names(sd_rep[[i]]$value) == "log_Rec"]) %>%
      dplyr::mutate(lwr = exp(log(value) - 1.96 * se),
                    upr = exp(log(value) + 1.96 * se),
                    Region = paste("Region", Region),
                    Pop = paste("Population", Pop),
                    Type = paste(Pop, 'Recruitment'),
                    Model = model_names[i])

    # Fishing Mortality
    f_plot_df <- reshape2::melt(rep[[i]]$Fmort) %>%
      dplyr::rename(Region = Var1, Year = Var2, Seas = Var3, Type = Var4) %>%
      dplyr::mutate(Region = paste("Region", Region),
                    Type = paste("Seas", Seas, "Fleet", Type, "F"),
                    Pop = NA,
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
    ggplot2::labs(x = 'Year', y = 'Value', color = 'Model') +
    ggplot2::coord_cartesian(ylim = c(0,NA)) +
    theme_sablefish()

  # recruitment
  rec_ts_plot <- ggplot2::ggplot(biom_rec_df %>% dplyr::filter(str_detect(Type, 'Recruitment')),
                                 ggplot2::aes(x = Year, y = value, ymin = lwr, ymax = upr, color = factor(Model), fill = factor(Model))) +
    ggplot2::geom_line(lwd = 0.9) +
    ggplot2::facet_grid(Type~Region, scales = 'free') +
    ggplot2::labs(x = 'Year', y = 'Recruitment', color = 'Model', fill = 'Model') +
    ggplot2::coord_cartesian(ylim = c(0,NA)) +
    theme_sablefish()

  # ssb
  ssb_ts_plot <- ggplot2::ggplot(biom_rec_df %>% dplyr::filter(str_detect(Type, 'SSB') & !str_detect(Type, 'Dynamic')),
                               ggplot2::aes(x = Year, y = value, ymin = lwr, ymax = upr, color = factor(Model), fill = factor(Model))) +
    ggplot2::geom_line(lwd = 0.9) +
    ggplot2::facet_grid(Type~Region, scales = 'free') +
    ggplot2::labs(x = 'Year', y = 'Spawning Stock Biomass', color = 'Model', fill = 'Model') +
    ggplot2::coord_cartesian(ylim = c(0,NA)) +
    theme_sablefish()

  # total biomass
  total_biom_plot <- ggplot2::ggplot(biom_rec_df %>% dplyr::filter(str_detect(Type, 'Total')),
                      ggplot2::aes(x = Year, y = value, ymin = lwr, ymax = upr, color = factor(Model), fill = factor(Model))) +
    ggplot2::geom_line(lwd = 0.9) +
    ggplot2::facet_grid(Type~Region, scales = 'free') +
    ggplot2::labs(x = 'Year', y = 'Total Biomass', color = 'Model', fill = 'Model') +
    ggplot2::coord_cartesian(ylim = c(0,NA)) +
    theme_sablefish()

  # dynamic b0
  ssb0_plot <- ggplot2::ggplot(biom_rec_df %>% dplyr::filter(str_detect(Type, 'Dynamic')),
                               ggplot2::aes(x = Year, y = value, ymin = lwr, ymax = upr, color = factor(Model), fill = factor(Model))) +
    ggplot2::geom_line(lwd = 0.9) +
    ggplot2::facet_grid(Type~Region, scales = 'free') +
    ggplot2::labs(x = 'Year', y = 'Unfished Spawning Stock Biomass', color = 'Model', fill = 'Model') +
    ggplot2::coord_cartesian(ylim = c(0,NA)) +
    theme_sablefish()

  # dynamic b0 and SSB
  ssb_ssb0_plot <- ggplot2::ggplot(biom_rec_df %>% dplyr::filter(str_detect(Type, 'SSB')),
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
#' Generates age- or length-based selectivity plots for all fishery and survey
#' fleets across one or more SPoRC model runs. Selectivity curves are faceted by
#' sex, region, and fleet, and coloured by year when multiple years are requested,
#' allowing visualisation of time-varying selectivity alongside model comparisons.
#'
#' @param rep List of length \code{n_models}, where each element is a SPoRC
#'   report list (i.e. the output of \code{obj$report()} after optimisation).
#'   Must contain \code{fish_sel} and \code{srv_sel} (age-based) and/or
#'   \code{fish_sel_l} and \code{srv_sel_l} (length-based) depending on
#'   \code{Selex_Type}.
#' @param model_names Character vector of length \code{n_models} giving display
#'   names for each model run. Used as the linetype legend label.
#' @param Selex_Type Character string specifying whether to plot age- or
#'   length-based selectivity. One of \code{"age"} (default) or \code{"length"}.
#'   Must match the \code{Selex_Type} used when fitting the model (i.e.
#'   \code{fish_sel_l} / \code{srv_sel_l} are only populated when
#'   \code{Selex_Type = 1} in the model).
#' @param year_indx Integer or integer vector of year indices to include in the
#'   plot. When multiple years are supplied, curves are coloured by year on a
#'   continuous viridis scale, making time-variation in selectivity visible.
#'   If \code{NULL} (default), only the terminal year is plotted.
#'
#' @return A list of two \code{ggplot} objects:
#' \describe{
#'   \item{[[1]] fish_sel_plot}{Fishery selectivity curves faceted by
#'     Sex × (Region + Fleet). Lines are coloured by year and distinguished
#'     by linetype across models.}
#'   \item{[[2]] srv_sel_plot}{Survey selectivity curves with identical
#'     faceting and aesthetic structure as the fishery plot.}
#' }
#'
#' @export get_selex_plot
#' @family Plotting
#'
#' @examples
#' \dontrun{
#'   # Terminal year only
#'   plots <- get_selex_plot(list(rep1, rep2), c("Base", "Sensitivity"))
#'
#'   # All years to visualise time-varying selectivity
#'   plots <- get_selex_plot(list(rep1), "Base", Selex_Type = "age", year_indx = 1:40)
#'
#'   plots[[1]] # fishery selectivity
#'   plots[[2]] # survey selectivity
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
    reshape_and_annotate_age <- function(df, model) {
      reshape2::melt(df) %>%
        dplyr::rename(Pop = Var1, Region = Var2, Year = Var3, Seas = Var4, Bin = Var5, Sex = Var6, Fleet = Var7) %>%
        dplyr::mutate(Pop = paste('Pop', Pop),
                      Seas = paste('Seas', Seas),
                      Region = paste("Region", Region),
                      Fleet  = paste("Fleet", Fleet),
                      Sex    = paste("Sex", Sex),
                      Model  = model)
    }

    reshape_and_annotate_len <- function(df, model) {
      reshape2::melt(df) %>%
        dplyr::rename(Region = Var1, Year = Var2, Bin = Var3, Sex = Var4, Fleet = Var5) %>%
        dplyr::mutate(Region = paste("Region", Region),
                      Fleet  = paste("Fleet", Fleet),
                      Sex    = paste("Sex", Sex),
                      Model  = model)
    }

    if(Selex_Type == 'age') {
      fishsel_plot_df <- rbind(fishsel_plot_df, reshape_and_annotate_age(fish_sel, model_names[i]))
      srvsel_plot_df  <- rbind(srvsel_plot_df, reshape_and_annotate_age(srv_sel,  model_names[i]))
    } else {
      fishsel_plot_df <- rbind(fishsel_plot_df, reshape_and_annotate_len(fish_sel, model_names[i]))
      srvsel_plot_df  <- rbind(srvsel_plot_df, reshape_and_annotate_len(srv_sel,  model_names[i]))
    }
  }

  if(is.null(year_indx)) year_indx <- max(fishsel_plot_df$Year)

  # fishery selectivity plot
  # Set up facet based on selectivity type
  if(Selex_Type == 'age') facet <- ggplot2::facet_grid(Pop + Sex ~ Region + Fleet + Seas)
  if(Selex_Type == 'length') facet <- ggplot2::facet_grid(Sex ~ Region + Fleet)

  fish_sel_plot <- ggplot2::ggplot(
    fishsel_plot_df %>%
      dplyr::filter(Year %in% c(year_indx)),
    ggplot2::aes(x = Bin, y = value, group = interaction(Year, Model), color = Year)
  ) +
    ggplot2::geom_line(aes(linetype = factor(Model)), linewidth = 1.3, alpha = 0.8) +
    facet +
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
    facet +
    ggplot2::scale_color_viridis_c() +
    ggplot2::labs(x = 'Bin', y = 'Survey selectivity', linetype = 'Model', color = 'Year') +
    ggplot2::coord_cartesian(ylim = c(0, NA)) +
    theme_sablefish() +
    ggplot2::theme(legend.key.width = unit(2, "lines"))

  return(list(fish_sel_plot, srv_sel_plot))
}

#' Get Plots of Biological Quantities
#'
#' Generates plots of movement probabilities, natural mortality, weight-at-age,
#' and maturity-at-age across one or more SPoRC model runs. All plots use a
#' user-specified year index; movement plots additionally show all seasons and
#' are returned as a separate plot per population.
#'
#' @param data List of length \code{n_models}, where each element is a SPoRC
#'   data list. Used to extract \code{WAA}, \code{MatAA}, \code{ages},
#'   \code{n_pop}, and \code{do_recruits_move}.
#' @param rep List of length \code{n_models}, where each element is a SPoRC
#'   report list (i.e. the output of \code{obj$report()} after optimisation).
#'   Must contain \code{Movement} and \code{natmort}.
#' @param model_names Character vector of length \code{n_models} giving display
#'   names for each model run. Used as colour legend labels across all plots.
#'
#' @return A list of four elements:
#' \describe{
#'   \item{[[1]] move_plot}{A list of \code{n_pop} \code{ggplot} objects, one
#'     per population. Each plot shows age-specific movement probabilities for
#'     all seasons at \code{year_indx}, faceted by Region_To × (Region_From +
#'     Season), with lines coloured by model and distinguished by sex linetype.
#'     If \code{do_recruits_move = 0}, age-1 fish are excluded. Y-axis is fixed
#'     to [0, 1].}
#'   \item{[[2]] natmort_plot}{Natural mortality at age at \code{year_indx},
#'     faceted by Region × Sex. Lines coloured by model and distinguished by
#'     population linetype.}
#'   \item{[[3]] waa_plot}{Spawning weight-at-age at \code{year_indx}, faceted
#'     by Region × (Sex + Season). Lines coloured by model and distinguished by
#'     population linetype.}
#'   \item{[[4]] mataa_plot}{Maturity-at-age at \code{year_indx}, faceted by
#'     Region × (Sex + Season). Lines coloured by model and distinguished by
#'     population linetype.}
#' }
#'
#' @export get_biological_plot
#' @family Plotting
#'
#' @examples
#' \dontrun{
#'   plots <- get_biological_plot(
#'     data        = list(data1, data2),
#'     rep         = list(rep1, rep2),
#'     model_names = c("Base", "Sensitivity"),
#'     year_indx   = 40
#'   )
#'
#'   plots[[1]][[1]]  # movement plot for population 1
#'   plots[[1]][[2]]  # movement plot for population 2, etc
#'   plots[[2]]       # natural mortality
#'   plots[[3]]       # weight-at-age
#'   plots[[4]]       # maturity-at-age
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
      dplyr::rename(Pop = pop, Region_From = from, Region_To = to, Seas = seas, Year = years, Age = ages, Sex = sexes) %>%
      {if (data[[i]]$do_recruits_move == 0) filter(., Age != min(data[[i]]$ages)) else .} %>%
      dplyr::mutate(Region_From = paste("From Region", Region_From),
                    Region_To = paste("To Region", Region_To),
                    Sex = paste("Sex", Sex),
                    Seas = paste("Seas", Seas),
                    Pop = paste('Population', Pop),
                    Model = model_names[i]
      )

    # Natural Mortality
    natmort_plot_tmp_df <- reshape2::melt(rep[[i]]$natmort) %>%
      dplyr::rename(Pop = Var1, Region = Var2, Year = Var3, Age = Var4, Sex = Var5) %>%
      dplyr::mutate(Region = paste("Region", Region),
                    Sex = paste("Sex", Sex),
                    Pop = paste('Population', Pop),
                    Model = model_names[i]
      )

    # Spawning Weight-at-age
    waa_plot_tmp_df <- reshape2::melt(data[[i]]$WAA) %>%
      dplyr::rename(Pop = Var1, Region = Var2, Year = Var3, Seas = Var4, Age = Var5, Sex = Var6) %>%
      dplyr::mutate(Region = paste("Region", Region),
                    Sex = paste("Sex", Sex),
                    Seas = paste("Seas", Seas),
                    Pop = paste('Population', Pop),
                    Model = model_names[i]
      )

    # Maturity at age
    mataa_plot_tmp_df <- reshape2::melt(data[[i]]$MatAA) %>%
      dplyr::rename(Pop = Var1, Region = Var2, Year = Var3, Seas = Var4, Age = Var5, Sex = Var6) %>%
      dplyr::mutate(Region = paste("Region", Region),
                    Sex = paste("Sex", Sex),
                    Seas = paste("Seas", Seas),
                    Pop = paste('Population', Pop),
                    Model = model_names[i]
      )

    # bind all
    move_plot_df <- rbind(move_plot_df, move_plot_tmp_df)
    natmort_plot_df <- rbind(natmort_plot_df, natmort_plot_tmp_df)
    waa_plot_df <- rbind(waa_plot_df, waa_plot_tmp_df)
    mataa_plot_df <- rbind(mataa_plot_df, mataa_plot_tmp_df)
  }

  # Movement plot
  move_plot <- list()
  for(i in 1:data[[1]]$n_pop) {
    move_plot[[i]] <- ggplot(move_plot_df %>% dplyr::filter(Year == max(move_plot_df$Year),
                                                            Pop == paste("Population", i)),
                             ggplot2::aes(x = Age, y = value, color = factor(Model), lty = Sex)) +
      ggplot2::geom_line(lwd = 1) +
      ggplot2::facet_grid(Region_To ~ Region_From + Seas) +
      ggplot2::labs(x = 'Age', y = 'Movement Probabilities', color = 'Model', lty = 'Sex', title = paste("Population", i)) +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      theme_sablefish() +
      ggplot2::theme(legend.key.width = unit(2, "lines"))
  }

  # Natural mortality plot
  natmort_plot <- ggplot2::ggplot(natmort_plot_df %>%
                                    dplyr::filter(Year == max(natmort_plot_df$Year)),
                                  ggplot2::aes(x = Age, y = value, color = factor(Model), lty = factor(Pop))) +
    ggplot2::geom_line(lwd = 2) +
    ggplot2::facet_grid(Region~Sex) +
    ggplot2::labs(x = 'Age', y = 'Natural Mortality', color = 'Model', lty=  'Population') +
    ggplot2::coord_cartesian(ylim = c(0, NA)) +
    theme_sablefish() +
    ggplot2::theme(legend.key.width = unit(2, "lines"))

  # Weight at age plot
  waa_plot <- ggplot2::ggplot(waa_plot_df %>%
                                dplyr::filter(Year == max(waa_plot_df$Year)),
                              ggplot2::aes(x = Age, y = value, color = factor(Model), lty = Pop)) +
    ggplot2::geom_line(lwd = 2) +
    ggplot2::facet_grid(Region~Sex + Seas) +
    ggplot2::labs(x = 'Age', y = 'Spawning Weight at Age', color = 'Model', lty = 'Population') +
    ggplot2::coord_cartesian(ylim = c(0, NA)) +
    theme_sablefish() +
    ggplot2::theme(legend.key.width = unit(2, "lines"))

  # Maturity plot
  mataa_plot <- ggplot2::ggplot(mataa_plot_df %>%
                                  dplyr::filter(Year == max(mataa_plot_df$Year)),
                                ggplot2::aes(x = Age, y = value, color = factor(Model), lty = Pop)) +
    ggplot2::geom_line(lwd = 2) +
    ggplot2::facet_grid(Region~Sex+Seas) +
    ggplot2::labs(x = 'Age', y = 'Maturity at Age', color = 'Model', lty = 'Population') +
    ggplot2::coord_cartesian(ylim = c(0, NA)) +
    theme_sablefish() +
    ggplot2::theme(legend.key.width = unit(2, "lines"))

  return(list(move_plot, natmort_plot, waa_plot, mataa_plot))

}

#' Get Data Fitted to Plot
#'
#' Produces a dot-plot showing which data types were fitted to in each year,
#' region, and model. Each point represents a year in which a given data source
#' was active (i.e. the corresponding \code{Use*} indicator equals 1). This is
#' useful for quickly auditing data availability and comparing data structures
#' across model configurations.
#'
#' @param data List of length \code{n_models}, where each element is a SPoRC
#'   data list. Pooled \code{Use*} indicator arrays are extracted with
#'   dimensions \code{[n_regions × n_years × n_seas × n_fleets]}:
#'   \code{UseSrvLenComps}, \code{UseSrvAgeComps}, \code{UseFishLenComps},
#'   \code{UseFishAgeComps}, \code{UseCatch}, \code{UseFishIdx},
#'   \code{UseSrvIdx}. Population-specific \code{Use*_pop} indicator arrays
#'   with dimensions \code{[n_pop × n_regions × n_years × n_seas × n_fleets]}
#'   are included when any element equals 1:
#'   \code{UseFishAgeComps_pop}, \code{UseFishLenComps_pop},
#'   \code{UseFishIdx_pop}, \code{UseSrvAgeComps_pop},
#'   \code{UseSrvLenComps_pop}, \code{UseSrvIdx_pop}. If
#'   \code{use_conv_fish_tagging} contains any 1s,
#'   \code{conv_tag_release_indicator} is also used to construct a tagging
#'   activity indicator array.
#' @param model_names Character vector of length \code{n_models} giving display
#'   names for each model run. Used as row facet labels.
#'
#' @return A single \code{ggplot} object: a dot-plot with Year on the x-axis
#'   and data type (labelled by source, population where applicable, season,
#'   and fleet) on the y-axis, faceted by Model × Region. Points appear only
#'   in years where the corresponding \code{Use*} indicator is 1. The legend
#'   is suppressed; data types are distinguished by y-axis position and fill
#'   colour.
#'
#' @export get_data_fitted_plot
#' @family Plotting
#' @examples
#' \dontrun{
#'   get_data_fitted_plot(
#'     data        = list(data1, data2),
#'     model_names = c("Base", "Sensitivity")
#'   )
#' }
get_data_fitted_plot <- function(data,
                                 model_names
                                 ) {

  data_plot_all_df <- data.frame()
  for(i in 1:length(data)) {

    # Get tag release indicator
    if(any(data[[i]]$use_conv_fish_tagging == 1)) {
      use_tag_indicator <- array(0, dim = c(max(data[[i]]$conv_tag_release_indicator[,1]), max(data[[i]]$conv_tag_release_indicator[,2]), max(data[[i]]$conv_tag_release_indicator[,3]) ))
      use_tag_indicator[data[[i]]$conv_tag_release_indicator[,1], data[[i]]$conv_tag_release_indicator[,2], data[[i]]$conv_tag_release_indicator[,3]] <- 1
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
    if (any(data[[i]]$use_conv_fish_tagging == 1)) {
      tag_df <- reshape2::melt(use_tag_indicator) %>%
        dplyr::rename(Region = Var1, Year = Var2, Seas = Var3) %>%
        dplyr::mutate(Type = paste('Tagging', "Seas", Seas), Fleet = NA)
      data_plot_df <- dplyr::bind_rows(data_plot_df, tag_df)
    }

    # Population-specific data indicators
    pop_plot_df <- data.frame()

    if(any(data[[i]]$UseFishAgeComps_pop == 1)) {
      pop_plot_df <- dplyr::bind_rows(pop_plot_df,
                                      reshape2::melt(data[[i]]$UseFishAgeComps_pop) %>%
                                        dplyr::rename(Pop = Var1, Region = Var2, Year = Var3, Seas = Var4, Fleet = Var5) %>%
                                        dplyr::mutate(Type = paste('Pop Fishery Ages', "Seas", Seas, "Fleet", Fleet)))
    }
    if(any(data[[i]]$UseFishLenComps_pop == 1)) {
      pop_plot_df <- dplyr::bind_rows(pop_plot_df,
                                      reshape2::melt(data[[i]]$UseFishLenComps_pop) %>%
                                        dplyr::rename(Pop = Var1, Region = Var2, Year = Var3, Seas = Var4, Fleet = Var5) %>%
                                        dplyr::mutate(Type = paste('Pop Fishery Lengths', "Seas", Seas, "Fleet", Fleet)))
    }
    if(any(data[[i]]$UseFishIdx_pop == 1)) {
      pop_plot_df <- dplyr::bind_rows(pop_plot_df,
                                      reshape2::melt(data[[i]]$UseFishIdx_pop) %>%
                                        dplyr::rename(Pop = Var1, Region = Var2, Year = Var3, Seas = Var4, Fleet = Var5) %>%
                                        dplyr::mutate(Type = paste('Pop Fishery Index', "Seas", Seas, "Fleet", Fleet)))
    }
    if(any(data[[i]]$UseSrvAgeComps_pop == 1)) {
      pop_plot_df <- dplyr::bind_rows(pop_plot_df,
                                      reshape2::melt(data[[i]]$UseSrvAgeComps_pop) %>%
                                        dplyr::rename(Pop = Var1, Region = Var2, Year = Var3, Seas = Var4, Fleet = Var5) %>%
                                        dplyr::mutate(Type = paste('Pop Survey Ages', "Seas", Seas, "Fleet", Fleet)))
    }
    if(any(data[[i]]$UseSrvLenComps_pop == 1)) {
      pop_plot_df <- dplyr::bind_rows(pop_plot_df,
                                      reshape2::melt(data[[i]]$UseSrvLenComps_pop) %>%
                                        dplyr::rename(Pop = Var1, Region = Var2, Year = Var3, Seas = Var4, Fleet = Var5) %>%
                                        dplyr::mutate(Type = paste('Pop Survey Lengths', "Seas", Seas, "Fleet", Fleet)))
    }
    if(any(data[[i]]$UseSrvIdx_pop == 1)) {
      pop_plot_df <- dplyr::bind_rows(pop_plot_df,
                                      reshape2::melt(data[[i]]$UseSrvIdx_pop) %>%
                                        dplyr::rename(Pop = Var1, Region = Var2, Year = Var3, Seas = Var4, Fleet = Var5) %>%
                                        dplyr::mutate(Type = paste('Pop Survey Index', "Seas", Seas, "Fleet", Fleet)))
    }

    if(nrow(pop_plot_df) > 0) {
      pop_plot_df <- pop_plot_df %>% dplyr::select(-Pop)
      data_plot_df <- dplyr::bind_rows(data_plot_df, pop_plot_df)
    }

    # Remove data not fitted to
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

#' Get Plot of Negative Log Likelihood Values
#'
#' Extracts, weights, and visualises all negative log-likelihood (nLL)
#' components from one or more SPoRC model runs. Likelihood weights stored in
#' the data list are applied to the relevant components (catch, indices,
#' recruitment, tagging, fishing mortality) before plotting, so reported values
#' reflect the actual contribution of each component to the joint nLL. Components
#' with a value of zero are silently excluded from both the plot and table.
#'
#' @param data List of length \code{n_models}, where each element is a SPoRC
#'   data list. Likelihood weights are extracted from \code{Wt_Catch},
#'   \code{Wt_F}, \code{Wt_Rec}, \code{Wt_Tagging}, \code{Wt_SrvIdx}, and
#'   \code{Wt_FishIdx}.
#' @param rep List of length \code{n_models}, where each element is a SPoRC
#'   report list (i.e. the output of \code{obj$report()} after optimisation).
#'   The following nLL components are extracted: \code{jnLL}, \code{h_nLL},
#'   \code{M_nLL}, \code{rec_region_prop_nLL}, \code{Rec_nLL},
#'   \code{Init_Rec_nLL}, \code{sel_nLL}, \code{conv_fish_tag_nLL},
#'   \code{Catch_nLL}, \code{Fmort_nLL}, \code{srv_q_nLL}, \code{fish_q_nLL},
#'   \code{SrvIdx_nLL}, \code{FishIdx_nLL}, \code{TagRep_nLL},
#'   \code{Movement_nLL}, \code{SrvAgeComps_nLL}, \code{FishAgeComps_nLL},
#'   \code{SrvLenComps_nLL}, \code{FishLenComps_nLL}, and the
#'   population-specific counterparts \code{Catch_pop_nLL},
#'   \code{FishIdx_pop_nLL}, \code{SrvIdx_pop_nLL},
#'   \code{FishAgeComps_pop_nLL}, \code{FishLenComps_pop_nLL},
#'   \code{SrvAgeComps_pop_nLL}, \code{SrvLenComps_pop_nLL}. Missing
#'   components are handled via \code{safe_extract}.
#' @param model_names Character vector of length \code{n_models} giving display
#'   names for each model run. Used as facet labels on the bar chart.
#'
#' @return A list of two objects:
#' \describe{
#'   \item{[[1]] nLL_plot}{A stacked bar chart (\code{ggplot}) with one facet
#'     per model. Bars show the weighted nLL contribution of each component,
#'     coloured by component type (Prior, Penalty, Catch, Index, Age, Length,
#'     Tagging, jnLL). Components with value 0 are excluded.}
#'   \item{[[2]] table_plot}{A \code{ggdraw} table (via \code{gridExtra} and
#'     \code{cowplot}) showing weighted nLL values in wide format, with one
#'     column per non-zero component and one row per model.}
#' }
#'
#' @export get_nLL_plot
#' @family Model Diagnostics
#'
#' @examples
#' \dontrun{
#'   out <- get_nLL_plot(
#'     data        = list(data1, data2),
#'     rep         = list(rep1, rep2),
#'     model_names = c("Base", "Sensitivity")
#'   )
#'   out[[1]]  # bar chart
#'   out[[2]]  # table
#' }
get_nLL_plot <- function(data,
                         rep,
                         model_names
                         ) {


  nLL_all_df <- data.frame() # empty dataframe
  for(i in 1:length(rep)) {

    # nLL values
    nLL_df <- data.frame(

      value = c(safe_extract(rep[[i]], "jnLL"),
                safe_extract(rep[[i]], "h_nLL"),
                safe_extract(rep[[i]], "R0_nLL"),
                safe_extract(rep[[i]], "M_nLL"),
                safe_extract(rep[[i]], "rec_region_prop_nLL"),
                sum(data[[i]]$Wt_Rec * safe_extract(rep[[i]], "Rec_nLL")),
                safe_extract(rep[[i]], "sel_nLL"),
                sum(data[[i]]$Wt_Tagging * safe_extract(rep[[i]], "conv_fish_tag_nLL")),
                sum(data[[i]]$Wt_Catch * safe_extract(rep[[i]], "Catch_nLL")),
                sum(data[[i]]$Wt_F * safe_extract(rep[[i]], "Fmort_nLL")),
                safe_extract(rep[[i]], "srv_q_nLL"),
                safe_extract(rep[[i]], "fish_q_nLL"),
                sum(data[[i]]$Wt_SrvIdx * safe_extract(rep[[i]], "SrvIdx_nLL")),
                safe_extract(rep[[i]], "TagRep_nLL"),
                sum(data[[i]]$Wt_FishIdx * safe_extract(rep[[i]], "FishIdx_nLL")),
                sum(safe_extract(data[[i]], "Wt_Init_Rec") * safe_extract(rep[[i]], "Init_Rec_nLL")),
                safe_extract(rep[[i]], "Movement_nLL"),
                sum(safe_extract(rep[[i]], "SrvAgeComps_nLL")),
                sum(safe_extract(rep[[i]], "FishAgeComps_nLL")),
                sum(safe_extract(rep[[i]], "SrvLenComps_nLL")),
                sum(safe_extract(rep[[i]], "FishLenComps_nLL")),
                # population-specific
                sum(data[[i]]$Wt_Catch_pop    * safe_extract(rep[[i]], "Catch_pop_nLL")),
                sum(data[[i]]$Wt_FishIdx_pop   * safe_extract(rep[[i]], "FishIdx_pop_nLL")),
                sum(data[[i]]$Wt_SrvIdx_pop    * safe_extract(rep[[i]], "SrvIdx_pop_nLL")),
                sum(safe_extract(rep[[i]], "FishAgeComps_pop_nLL")),
                sum(safe_extract(rep[[i]], "FishLenComps_pop_nLL")),
                sum(safe_extract(rep[[i]], "SrvAgeComps_pop_nLL")),
                sum(safe_extract(rep[[i]], "SrvLenComps_pop_nLL"))),

      name = c("jnLL", "Steepness Prior", "R0 Prior", "M Prior", "Recruitment Prop Prior", "Recruitment Penalty",
               "Selectivity Penalty", "Conventional Tagging nLL", "Catch nLL", "Fishing Mortality Penalty",
               "Survey Q Prior", "Fishery Q Prior", "Survey Index nLL", "Tag Reporting Prior",
               "Fishery Index nLL", "Initial Age Penalty", "Movement Prior",
               "Survey Age nLL", "Fishery Age nLL", "Survey Length nLL", "Fishery Length nLL",
               "Pop Catch nLL", "Pop Fishery Index nLL", "Pop Survey Index nLL",
               "Pop Fishery Age nLL", "Pop Fishery Length nLL",
               "Pop Survey Age nLL", "Pop Survey Length nLL"),

      type = c('jnLL', 'Prior', "Prior", "Prior", "Prior", "Penalty", "Penalty", "Tagging",
               "Catch", "Penalty", "Prior", "Prior", "Index",
               "Prior", "Index", "Penalty", "Prior", "Age",
               "Age", "Length", "Length",
               "Catch", "Index", "Index",
               "Age", "Length", "Age", "Length"),

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
#' Plots observed survey and fishery indices alongside model-predicted values
#' for one or more SPoRC model runs. Observed values are shown as points with
#' approximate 95% confidence intervals; predicted trajectories are overlaid as
#' lines coloured by model. Years where the observed index is zero (i.e.
#' \code{Use*Idx = 0}) are excluded from both the points and lines.
#'
#' @param data List of length \code{n_models}, where each element is a SPoRC
#'   data list. Passed to \code{get_idx_fits} along with \code{rep[[i]]}; year
#'   labels are taken from \code{data[[i]]$years}.
#' @param rep List of length \code{n_models}, where each element is a SPoRC
#'   report list (i.e. the output of \code{obj$report()} after optimisation).
#'   Predicted index values are extracted internally via \code{get_idx_fits}.
#' @param model_names Character vector of length \code{n_models} giving display
#'   names for each model run. Used as the colour legend label on predicted
#'   trajectories.
#'
#' @return A single \code{ggplot} object. Observed indices are shown as
#'   \code{geom_pointrange} (black) with lower and upper confidence interval
#'   bounds from \code{get_idx_fits}. Predicted indices are shown as
#'   \code{geom_line} coloured by model.
#'
#' @export get_idx_fits_plot
#' @family Model Diagnostics
#'
#' @examples
#' \dontrun{
#'   get_idx_fits_plot(
#'     data        = list(data1, data2),
#'     rep         = list(rep1, rep2),
#'     model_names = c("Base", "Sensitivity")
#'   )
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

#' Get Catch and Discard Fits Plot
#'
#' Plots observed catch and discard time series alongside model-predicted
#' values for one or more SPoRC model runs, for both pooled (region-level)
#' and population-specific data streams.
#'
#' Pooled predictions are summed across populations before comparison with
#' observed data. Years where observed values are zero are excluded from both
#' observed and predicted layers to avoid issues in lognormal confidence
#' interval construction.
#'
#' @param data List of length \code{n_models}, where each element is a SPoRC
#'   data list. \code{ObsCatch}
#'   \code{[n_regions × n_yrs × n_seas × n_fish_fleets]} provides pooled observed
#'   catch, and \code{ObsDiscard} provides pooled observed discards.
#'   \code{Wt_Catch} and \code{Wt_Discard}, together with \code{ln_sigmaC},
#'   are used to reconstruct observation-error standard deviations for
#'   confidence interval construction.
#'
#'   When \code{UseCatch_pop} or \code{UseDiscard_pop} contain active elements,
#'   population-level data (\code{ObsCatch_pop}, \code{ObsDiscard_pop}) and
#'   corresponding weights are also used.
#'
#' @param rep List of length \code{n_models}, where each element is a SPoRC
#'   model report (output of \code{obj$report()} after optimisation).
#'   \code{PredCatch} and \code{PredDiscard}
#'   \code{[n_pop × n_regions × n_yrs × n_seas × n_fish_fleets]} are summed
#'   across populations for pooled trajectories and used directly for
#'   population-specific trajectories.
#'   \code{ln_sigmaC} and \code{ln_sigmaC_pop} are used for confidence interval
#'   construction.
#'
#' @param model_names Character vector of length \code{n_models} giving display
#'   names for each model run. Used in the legend for predicted trajectories.
#'
#' @return A list of \code{ggplot} objects:
#' \describe{
#'   \item{[[1]] catch_fit_rg_plot}{Pooled catch fits (region-level). Produced
#'   when \code{UseCatch == 1}. Observations shown as \code{geom_pointrange}
#'   with 95\% lognormal confidence intervals; predictions shown as lines
#'   colored by model. Faceted by (Season × Fleet) × Region with free y-scales.}
#'
#'   \item{[[2]] catch_fit_pop_plot}{Population-specific catch fits. Produced
#'   when \code{UseCatch_pop == 1}. Faceted by (Population × Season × Fleet) ×
#'   Region. Returns \code{NULL} if not used.}
#'
#'   \item{[[3]] discard_fit_rg_plot}{Pooled discard fits (region-level).
#'   Produced when \code{UseDiscard == 1}. Same structure as catch plots.}
#'
#'   \item{[[4]] discard_fit_pop_plot}{Population-specific discard fits.
#'   Produced when \code{UseDiscard_pop == 1}. Returns \code{NULL} if not used.}
#' }
#'
#' @details
#' This function produces diagnostic plots for both catch and discard data.
#' Observed and predicted time series are shown for each, with lognormal
#' confidence intervals derived from \code{ln_sigmaC} (or population-level
#' equivalents) and sampling weights.
#'
#' Observations equal to zero are excluded prior to plotting to avoid issues
#' in log-space confidence interval construction.
#'
#' Predicted catch and discard are aggregated across populations for pooled
#' diagnostics, while population-specific plots use unaggregated outputs.
#'
#' @export get_catch_fits_plot
#' @family Model Diagnostics
get_catch_fits_plot <- function(data,
                                rep,
                                model_names
                                ) {

  # Plot catch fits
  catch_fit_rg_plot <- if(any(data[[1]]$UseCatch == 1)) {

    catch_fits_rg_all <- data.frame()

    # get catch fits data
    for(i in 1:length(rep)) {

      # Get catch fits
      catch_fits <- reshape2::melt(rep[[i]]$PredCatch) %>%
        dplyr::rename(Pop = Var1, Region = Var2, Year = Var3, Seas = Var4, Fleet = Var5) %>%
        dplyr::group_by(Region, Year, Seas, Fleet) %>%
        dplyr::summarize(value = sum(value)) %>%
        dplyr::ungroup() %>%
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

      catch_fits_rg_all <- rbind(catch_fits_rg_all, catch_fits) # bind

    }

     ggplot2::ggplot() +
      ggplot2::geom_line(catch_fits_rg_all %>% dplyr::filter(obs != 0),
                         mapping = ggplot2::aes(x = Year, y = value, color = factor(Model)), lwd = 1.3) +
      ggplot2::geom_pointrange(catch_fits_rg_all %>% dplyr::filter(obs != 0),
                               mapping = ggplot2::aes(x = Year, y = obs, ymin = exp(log(obs) - 1.96 * se),
                                                      ymax = exp(log(obs) + 1.96 * se)), color = 'black') +
      ggplot2::labs(x = "Year", y = 'Catch', color = 'Model') +
      theme_sablefish() +
      ggplot2::coord_cartesian(ylim = c(0,NA)) +
      ggplot2::facet_grid(Seas_Fleet~Region, scales = 'free_y')
  } else NULL

  # Plot discards fits
  discard_fit_rg_plot <- if(any(data[[1]]$UseDiscard == 1)) {

    discard_fits_rg_all <- data.frame()

    # get discard fits data
    for(i in 1:length(rep)) {

      # Get discard fits
      discard_fits <- reshape2::melt(rep[[i]]$PredDiscard) %>%
        dplyr::rename(Pop = Var1, Region = Var2, Year = Var3, Seas = Var4, Fleet = Var5) %>%
        dplyr::group_by(Region, Year, Seas, Fleet) %>%
        dplyr::summarize(value = sum(value)) %>%
        dplyr::ungroup() %>%
        dplyr::left_join(
          reshape2::melt(data[[i]]$ObsDiscard) %>%
            dplyr::left_join(reshape2::melt(exp(rep[[i]]$ln_sigmaC) / data[[i]]$Wt_Discard) %>%
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

      discard_fits_rg_all <- rbind(discard_fits_rg_all, discard_fits) # bind

    }

    ggplot2::ggplot() +
      ggplot2::geom_line(discard_fits_rg_all %>% dplyr::filter(obs != 0),
                         mapping = ggplot2::aes(x = Year, y = value, color = factor(Model)), lwd = 1.3) +
      ggplot2::geom_pointrange(discard_fits_rg_all %>% dplyr::filter(obs != 0),
                               mapping = ggplot2::aes(x = Year, y = obs, ymin = exp(log(obs) - 1.96 * se),
                                                      ymax = exp(log(obs) + 1.96 * se)), color = 'black') +
      ggplot2::labs(x = "Year", y = 'Discard', color = 'Model') +
      theme_sablefish() +
      ggplot2::coord_cartesian(ylim = c(0,NA)) +
      ggplot2::facet_grid(Seas_Fleet~Region, scales = 'free_y')
  } else NULL

  catch_fit_pop_plot <- if(any(data[[1]]$UseCatch_pop == 1)) {

    catch_fits_pop_all <- data.frame()

    # get catch fits data
    for(i in 1:length(rep)) {

      # Get catch fits
      catch_fits <- reshape2::melt(rep[[i]]$PredCatch) %>%
        dplyr::rename(Pop = Var1, Region = Var2, Year = Var3, Seas = Var4, Fleet = Var5) %>%
        dplyr::left_join(
          reshape2::melt(data[[i]]$ObsCatch_pop) %>%
            dplyr::left_join(reshape2::melt(exp(rep[[i]]$ln_sigmaC_pop) / data[[i]]$Wt_Catch_pop) %>%
                               dplyr::rename(se = value),
                             by = c("Var1", "Var2", "Var3", "Var4", "Var5")) %>%
            dplyr::rename(Pop = Var1, Region = Var2, Year = Var3, Seas = Var4, Fleet = Var5, obs = value),
          by = c("Pop", "Region", "Year", "Seas", "Fleet")
        ) %>%
        dplyr::mutate(Model = model_names[i],
                      Pop = paste('Pop', Pop),
                      Region = paste('Region', Region),
                      Seas = paste("Seas", Seas),
                      Fleet = paste('Fleet', Fleet),
                      Pop_Seas_Fleet = paste(Pop, Seas, Fleet))

      catch_fits_pop_all <- rbind(catch_fits_pop_all, catch_fits) # bind
    }

    ggplot2::ggplot() +
      ggplot2::geom_line(catch_fits_pop_all %>% dplyr::filter(obs != 0),
                         mapping = ggplot2::aes(x = Year, y = value, color = factor(Model)), lwd = 1.3) +
      ggplot2::geom_pointrange(catch_fits_pop_all %>% dplyr::filter(obs != 0),
                               mapping = ggplot2::aes(x = Year, y = obs, ymin = exp(log(obs) - 1.96 * se),
                                                      ymax = exp(log(obs) + 1.96 * se)), color = 'black') +
      ggplot2::labs(x = "Year", y = 'Population-Specific Catch', color = 'Model') +
      theme_sablefish() +
      ggplot2::coord_cartesian(ylim = c(0,NA)) +
      ggplot2::facet_grid(Pop_Seas_Fleet~Region, scales = 'free_y')
  } else NULL

  discard_fit_pop_plot <- if(any(data[[1]]$UseDiscard_pop == 1)) {

    discard_fits_pop_all <- data.frame()

    # get discard fits data
    for(i in 1:length(rep)) {

      # Get discard fits
      discard_fits <- reshape2::melt(rep[[i]]$PredDiscard) %>%
        dplyr::rename(Pop = Var1, Region = Var2, Year = Var3, Seas = Var4, Fleet = Var5) %>%
        dplyr::left_join(
          reshape2::melt(data[[i]]$ObsDiscard_pop) %>%
            dplyr::left_join(reshape2::melt(exp(rep[[i]]$ln_sigmaC_pop) / data[[i]]$Wt_Discard_pop) %>%
                               dplyr::rename(se = value),
                             by = c("Var1", "Var2", "Var3", "Var4", "Var5")) %>%
            dplyr::rename(Pop = Var1, Region = Var2, Year = Var3, Seas = Var4, Fleet = Var5, obs = value),
          by = c("Pop", "Region", "Year", "Seas", "Fleet")
        ) %>%
        dplyr::mutate(Model = model_names[i],
                      Pop = paste('Pop', Pop),
                      Region = paste('Region', Region),
                      Seas = paste("Seas", Seas),
                      Fleet = paste('Fleet', Fleet),
                      Pop_Seas_Fleet = paste(Pop, Seas, Fleet))

      discard_fits_pop_all <- rbind(discard_fits_pop_all, discard_fits) # bind
    }

    ggplot2::ggplot() +
      ggplot2::geom_line(discard_fits_pop_all %>% dplyr::filter(obs != 0),
                         mapping = ggplot2::aes(x = Year, y = value, color = factor(Model)), lwd = 1.3) +
      ggplot2::geom_pointrange(discard_fits_pop_all %>% dplyr::filter(obs != 0),
                               mapping = ggplot2::aes(x = Year, y = obs, ymin = exp(log(obs) - 1.96 * se),
                                                      ymax = exp(log(obs) + 1.96 * se)), color = 'black') +
      ggplot2::labs(x = "Year", y = 'Population-Specific Discard', color = 'Model') +
      theme_sablefish() +
      ggplot2::coord_cartesian(ylim = c(0,NA)) +
      ggplot2::facet_grid(Pop_Seas_Fleet~Region, scales = 'free_y')
  } else NULL

  return(list(catch_fit_rg_plot, catch_fit_pop_plot, discard_fit_rg_plot, discard_fit_pop_plot))
}


#' Get Retrospective Plot
#'
#' Generates three retrospective diagnostic plots from a set of peeled model
#' runs: a relative-difference plot with Mohn's rho, an absolute-scale
#' trajectory plot, and a cohort-tracking squid plot for recruitment. Together
#' these plots diagnose systematic retrospective patterns in SSB and recruitment
#' estimates across sequential data removals.
#'
#' @param retro_output Data frame produced by \code{do_retrospective}, containing
#'   columns \code{Year}, \code{value}, \code{peel}, \code{Type} (e.g.
#'   \code{"SSB"}, \code{"Recruitment"}), \code{Pop}, and \code{Region}. Rows
#'   with \code{peel = 0} represent the full-data terminal run; rows with
#'   \code{peel > 0} represent sequential data removals. Rows with
#'   \code{value = 0} are excluded before plotting.
#' @param Rec_Age Integer giving the age at recruitment (e.g. \code{2} for
#'   age-2 recruitment). Used to convert model year to cohort year in the squid
#'   plot (\code{cohort = Year - Rec_Age}).
#'
#' @return A list of three \code{ggplot} objects:
#' \describe{
#'   \item{[[1]] retro_plot}{Relative-difference plot. Each line shows the
#'     proportional deviation of a peeled run from the terminal-year estimate
#'     at each year, computed via \code{get_retrospective_relative_difference}.
#'     Terminal-year points (where \code{peel = max(Year) - Year}) are
#'     highlighted. Mohn's rho (mean relative difference across terminal-year
#'     points) is annotated on each facet panel. Faceted by Region ×
#'     (Type + Population); lines and points coloured by peel year on a
#'     continuous viridis scale.}
#'   \item{[[2]] abs_retro_plot}{Absolute-scale trajectory plot. Peeled runs
#'     (\code{peel > 0}) are drawn as solid lines coloured by peel; the
#'     full-data run (\code{peel = 0}) is overlaid as a dashed black line.
#'     Faceted by Region × (Type + Population) with free y-scales.}
#'   \item{[[3]] squid_plot}{Cohort-tracking (squid) plot for recruitment.
#'     Restricted to the 10 most recent cohorts. X-axis shows years since the
#'     cohort was last estimated (\code{terminal - Year - 1}), y-axis shows
#'     recruitment value, and lines are grouped and coloured by cohort year.
#'     Faceted by Population × Region with free y-scales.}
#' }
#'
#' @export get_retrospective_plot
#' @family Model Diagnostics
#'
#' @examples
#' \dontrun{
#'   retro <- do_retrospective(
#'     n_retro        = 7,
#'     data           = data,
#'     parameters     = parameters,
#'     mapping        = mapping,
#'     random         = NULL,
#'     do_par         = TRUE,
#'     n_cores        = 7,
#'     do_francis     = FALSE,
#'     n_francis_iter = NULL
#'   )
#'   plots <- get_retrospective_plot(retro, Rec_Age = 2)
#'   plots[[1]]  # relative difference + Mohn's rho
#'   plots[[2]]  # absolute trajectories
#'   plots[[3]]  # squid plot
#' }
get_retrospective_plot <- function(retro_output, Rec_Age) {

  # Get relative differences
  ret_df <- get_retrospective_relative_difference(retro_output)
  ret_df <- ret_df %>% dplyr::filter(!is.nan(rd)) # remove 0s

  # some quick naming
  retro_output <- retro_output %>%
    dplyr::mutate(Pop = paste('Pop', Pop),
                  Region = paste('Region', Region)) %>%
    dplyr::filter(value != 0)

  # Get mohns rho (mean relative difference for a given terminal year peel to terminal year estimates)
  mohns_rho <- ret_df %>%
    dplyr::filter(peel == (max(Year) - Year)) %>%
    dplyr::group_by(Type, Pop, Region) %>%
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
    ggplot2::facet_grid(Region~Type + Pop, scales = 'free') +
    theme_sablefish() +
    ggplot2::theme(legend.position = 'top')

  # get absolute retro plot
  abs_retro_plot <- ggplot2::ggplot() +
    ggplot2::geom_line(retro_output %>% filter(peel != 0), mapping = ggplot2::aes(x = Year, y = value, group = peel, color = peel), lwd = 1) +
    ggplot2::geom_line(retro_output %>% filter(peel == 0), mapping = ggplot2::aes(x = Year, y = value), lty = 2, lwd = 1) +
    ggplot2::scale_color_viridis_c() +
    ggplot2::coord_cartesian(ylim = c(0, NA)) +
    ggplot2::facet_wrap(Region~Type+Pop, scales = 'free_y') +
    theme_sablefish() +
    ggplot2::labs(x = 'Year', y = 'Value', color = 'Peel')

  # get squid plot
  squid_plot <- retro_output %>%
    dplyr::mutate(Year = Year, terminal = max(retro_output$Year) - peel, cohort = Year - Rec_Age, years_est = terminal-Year) %>%
    dplyr::filter(Type == 'Recruitment', cohort %in% seq(max(retro_output$Year) - 10, max(retro_output$Year), 1), terminal != Year) %>%
    ggplot2::ggplot(ggplot2::aes(x = years_est - 1, y = value, group = Year, color = factor(cohort))) +
    ggplot2::geom_line(lwd = 1.3) +
    ggplot2::geom_point(size = 4) +
    ggplot2::facet_grid(Pop~Region, scales = 'free') +
    ggplot2::theme_bw(base_size = 15) +
    ggplot2::labs(x = 'Years since cohort was last estimated', y = 'Recruitment (millions)', color = 'Cohort')

  return(list(retro_plot, abs_retro_plot, squid_plot))
}

#' Plotting Function for All Basic Quantities
#'
#' Convenience wrapper that calls all core SPoRC plotting functions and writes
#' their output to a single PDF file. Equivalent to calling
#' \code{get_biological_plot}, \code{get_data_fitted_plot}, \code{get_ts_plot},
#' \code{get_selex_plot}, and \code{get_nLL_plot} in sequence and printing each
#' to the same device.
#'
#' @param data List of length \code{n_models}, where each element is a SPoRC
#'   data list. Passed to \code{get_biological_plot}, \code{get_data_fitted_plot},
#'   and \code{get_nLL_plot}.
#' @param rep List of length \code{n_models}, where each element is a SPoRC
#'   report list (i.e. the output of \code{obj$report()} after optimisation).
#'   Passed to \code{get_biological_plot}, \code{get_ts_plot},
#'   \code{get_selex_plot}, and \code{get_nLL_plot}.
#' @param sd_rep List of length \code{n_models}, where each element is a SPoRC
#'   sdreport list (i.e. the output of \code{sdreport(obj)}). Passed to
#'   \code{get_ts_plot} for delta-method confidence intervals.
#' @param model_names Character vector of length \code{n_models} giving display
#'   names for each model run. Passed to all sub-functions as legend and facet
#'   labels.
#' @param out_path Path to the output directory. The PDF is written to
#'   \code{here::here(out_path, "plot_results.pdf")} at 25 × 13 inches per page.
#'
#' @return Called for its side effect. Writes \code{plot_results.pdf} to
#'   \code{out_path} and returns \code{NULL} invisibly.
#'
#' @export plot_all_basic
#' @family Plotting
#'
#' @examples
#' \dontrun{
#'   plot_all_basic(
#'     data        = list(data1, data2),
#'     rep         = list(rep1, rep2),
#'     sd_rep      = list(sd_rep1, sd_rep2),
#'     model_names = c("Base", "Sensitivity"),
#'     out_path    = here::here()
#'   )
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
#' Calculates biological and fishery reference points and performs short-term
#' population projections to estimate terminal spawning biomass, catch advice,
#' and reference point ratios by model and region. Returns both a tidy
#' data frame and a formatted table plot of the assembled quantities.
#'
#' @note The quantities returned by this function are \strong{approximate and
#'   should not be treated as official catch advice}. This wrapper provides
#'   only a simplified projection interface; full projection capability,
#'   including stochastic recruitment, closed-loop feedback, and
#'   fleet-specific harvest control rules, requires calling
#'   \code{\link{Do_Population_Projection}} directly. Results here are
#'   intended for rapid model comparison and diagnostic screening only.
#'
#' @param data A list of length \code{n_models}, where each element is a
#'   SPoRC-formatted data list containing region, year, age, fleet, and
#'   biological inputs (e.g., weight-at-age, maturity, natural mortality).
#' @param rep A list of length \code{n_models}, where each element is a
#'   SPoRC-formatted report list (i.e., the output of \code{obj$report()}
#'   after optimisation). Each element must include recruitment, selectivity,
#'   fishing mortality, and numbers-at-age arrays.
#' @param reference_points_opt A named list of options passed to
#'   \code{\link{Get_Reference_Points}}. Required elements:
#'   \describe{
#'     \item{\code{SPR_x}}{Target spawning potential ratio (e.g., \code{0.4})
#'       for SPR-based F reference points. May be \code{NULL} when
#'       \code{type = "bh_msy"}.}
#'     \item{\code{t_spawn}}{Fraction of the year elapsed before spawning
#'       occurs (e.g., \code{0} for start-of-year, \code{0.5} for
#'       mid-year).}
#'     \item{\code{sex_ratio_f}}{Array of dimension
#'       \code{(n_pop, n_regions)} giving the proportion of recruits that
#'       are female.}
#'     \item{\code{calc_rec_st_yr}}{Index of the first model year to include
#'       when averaging recruitment for reference point calculations.}
#'     \item{\code{rec_age}}{Age at recruitment (used to lag the recruitment
#'       series relative to terminal year).}
#'     \item{\code{type}}{Reference point calculation method; e.g.,
#'       \code{"multi_region"} or \code{"single_region"}.}
#'     \item{\code{what}}{Output selector passed to
#'       \code{Get_Reference_Points}; e.g., \code{"global_SPR"} or
#'       \code{"local_MSY"}.}
#'   }
#' @param proj_model_opt A named list of projection settings passed to
#'   \code{\link{Do_Population_Projection}}. Required elements:
#'   \describe{
#'     \item{\code{n_proj_yrs}}{Number of years to project forward.}
#'     \item{\code{n_avg_yrs}}{Number of terminal model years over which
#'       demographic inputs (selectivity, weight-at-age, maturity, natural
#'       mortality, movement) are averaged before being held constant across
#'       the projection period.}
#'     \item{\code{HCR_function}}{A harvest control rule function with
#'       signature \code{function(x, frp, brp, ...)}, where \code{x} is
#'       current biomass, \code{frp} is the F reference point, and \code{brp}
#'       is the biomass reference point.}
#'     \item{\code{recruitment_opt}}{Recruitment assumption for projection
#'       years. One of \code{"mean_rec"}, \code{"bh_rec"},
#'       \code{"zero_rec"}, or \code{"inv_gauss"}. See Details.}
#'     \item{\code{fmort_opt}}{How fishing mortality is set during the
#'       projection. One of \code{"input"} (hold terminal F constant) or
#'       \code{"HCR"} (apply \code{HCR_function}).}
#'   }
#' @param model_names Character vector of length \code{n_models} giving
#'   display names for each model run (e.g., \code{c("Base", "Alt1")}).
#'
#' @return A list of length 2:
#'   \describe{
#'     \item{\code{[[1]]}}{A data frame of key quantities by model and
#'       region, with columns \code{Model}, \code{Region},
#'       \code{Terminal_SSB}, \code{Terminal_SSB0}, \code{Terminal_F},
#'       \code{Catch_Advice}, \code{B_Ref_Pt}, \code{F_Ref_Pt},
#'       \code{B_over_B_Ref}, \code{B_over_DynB_Ref}, and
#'       \code{F_over_F_Ref}.}
#'     \item{\code{[[2]]}}{A \code{cowplot} \code{ggdraw} object rendering
#'       the same quantities as a formatted table, suitable for inclusion in
#'       a PDF report.}
#'   }
#'
#' @details
#' For each model, the function: (1) computes reference points via
#' \code{Get_Reference_Points()}; (2) averages demographic inputs over the
#' last \code{n_avg_yrs} model years; (3) projects the population forward
#' \code{n_proj_yrs} years via \code{Do_Population_Projection()}; and (4)
#' extracts terminal SSB, dynamic unfished SSB, catch advice (year 2 of the
#' projection), and status ratios.
#'
#' When \code{recruitment_opt = "bh_rec"}, Beverton-Holt stock-recruit
#' parameters are passed to the projection via an internal \code{srr_opt}
#' list constructed from year-1 demographics to approximate unfished SSB. When
#' \code{recruitment_opt = "inv_gauss"}, a warning is issued because only a
#' single deterministic simulation is run; stochastic recruitment options
#' should be used within a full MSE loop rather than here.
#'
#' @seealso \code{\link{Get_Reference_Points}},
#'   \code{\link{Do_Population_Projection}}
#'
#' @export get_key_quants
#' @family Reference Points and Projections
#'
#' @examples
#' \dontrun{
#'   reference_points_opt <- list(
#'     SPR_x          = 0.4,
#'     t_spawn        = 0,
#'     sex_ratio_f    = array(0.5, dim = c(n_pop, n_regions)),
#'     calc_rec_st_yr = 20,
#'     rec_age        = 2,
#'     type           = "multi_region",
#'     what           = "global_SPR"
#'   )
#'
#'   proj_model_opt <- list(
#'     n_proj_yrs      = 2,
#'     n_avg_yrs       = 1,
#'     HCR_function    = function(x, frp, brp, alpha = 0.05) {
#'       stock_status <- x / brp
#'       if (stock_status >= 1)                            frp
#'       else if (stock_status > alpha) frp * (stock_status - alpha) / (1 - alpha)
#'       else                                              0
#'     },
#'     recruitment_opt = "mean_rec",
#'     fmort_opt       = "HCR"
#'   )
#'
#'   out <- get_key_quants(
#'     data                  = list(mlt_rg_sable_data),
#'     rep                   = list(mlt_rg_sable_rep),
#'     reference_points_opt  = reference_points_opt,
#'     proj_model_opt        = proj_model_opt,
#'     model_names           = "Model 1"
#'   )
#'
#'   out[[1]]  # key quantities data frame
#'   out[[2]]  # table plot
#' }
get_key_quants <- function(data,
                           rep,
                           reference_points_opt,
                           proj_model_opt,
                           model_names
) {

  # required elements for reference points opt
  required <- c("SPR_x", "t_spawn", "sex_ratio_f", "calc_rec_st_yr", "rec_age", "type", "what")
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
                                        t_spawn = reference_points_opt$t_spawn,
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
    t_spawn <- reference_points_opt$t_spawn # spawn timing

    # terminal estimates
    terminal_NAA <-  array(rep[[i]]$NAA[,,length(data[[i]]$years),,,], dim = c(data[[i]]$n_pop, data[[i]]$n_regions, data[[i]]$n_seas, length(data[[i]]$ages), data[[i]]$n_sexes)) # terminal NAA
    terminal_NAA0 <-  array(rep[[i]]$NAA0[,,length(data[[i]]$years),,,], dim = c(data[[i]]$n_pop, data[[i]]$n_regions, data[[i]]$n_seas, length(data[[i]]$ages), data[[i]]$n_sexes)) # terminal NAA
    terminal_F <- array(rep[[i]]$Fmort[,length(data[[i]]$years),,], dim = c(data[[i]]$n_regions, data[[i]]$n_seas, data[[i]]$n_fish_fleets)) # terminal F
    recruitment <- array(rep[[i]]$Rec[,,reference_points_opt$calc_rec_st_yr:(length(data[[i]]$years) - reference_points_opt$rec_age)],
                         dim = c(data[[i]]$n_pop, data[[i]]$n_regions, length(reference_points_opt$calc_rec_st_yr:(length(data[[i]]$years) - reference_points_opt$rec_age)))) # recruitment

    # demographics
    # determine years to average over demogrphaics
    n_avg_yrs <- proj_model_opt$n_avg_yrs
    n_proj_yrs <- proj_model_opt$n_proj_yrs
    n_yrs <- length(data[[i]]$years)
    avg_yrs <- (n_yrs - n_avg_yrs + 1):n_yrs

    # spawning weight-at-age
    WAA_avg <- apply(data[[i]]$WAA[,,avg_yrs,,,,drop = FALSE], c(1, 2, 4, 5, 6), mean)
    WAA <- aperm(
      array(rep(WAA_avg, times = n_proj_yrs),
            dim = c(data[[i]]$n_pop, data[[i]]$n_regions, data[[i]]$n_seas, length(data[[i]]$ages), data[[i]]$n_sexes, n_proj_yrs)),
      perm = c(1, 2, 6, 3, 4, 5))

    WAA_fish_avg <- apply(data[[i]]$WAA_fish[,,avg_yrs,,,,,drop = FALSE], c(1, 2, 4, 5, 6, 7), mean)
    WAA_fish <- aperm(
      array(rep(WAA_fish_avg, times = n_proj_yrs),
            dim = c(data[[i]]$n_pop, data[[i]]$n_regions, data[[i]]$n_seas, length(data[[i]]$ages), data[[i]]$n_sexes, data[[i]]$n_fish_fleets, n_proj_yrs)),
      perm = c(1, 2, 7, 3, 4, 5, 6))

    # maturity at age
    MatAA_avg <- apply(data[[i]]$MatAA[,,avg_yrs,,,,drop = FALSE], c(1, 2, 4, 5, 6), mean)
    MatAA <- aperm(
      array(rep(MatAA_avg, times = n_proj_yrs),
            dim = c(data[[i]]$n_pop, data[[i]]$n_regions, data[[i]]$n_seas, length(data[[i]]$ages), data[[i]]$n_sexes, n_proj_yrs)),
      perm = c(1, 2, 6, 3, 4, 5))

    # natural mortality
    natmort_avg <- apply(rep[[i]]$natmort[,,avg_yrs,,,drop = FALSE], c(1, 2, 4, 5), mean)
    natmort <- aperm(array(rep(natmort_avg, times = n_proj_yrs),
                           dim = c(data[[i]]$n_pop, data[[i]]$n_regions, length(data[[i]]$ages),
                                   data[[i]]$n_sexes, n_proj_yrs)), perm = c(1, 2, 5, 3, 4))
    # fishery selectivity
    fish_sel_avg <- apply(rep[[i]]$fish_sel[,avg_yrs,,,,drop = FALSE], c(1, 3, 4, 5), mean)
    fish_sel <- aperm(
      array(rep(fish_sel_avg, times = n_proj_yrs),
            dim = c(data[[i]]$n_regions, length(data[[i]]$ages), data[[i]]$n_sexes, data[[i]]$n_fish_fleets, n_proj_yrs)),
      perm = c(1, 5, 2, 3, 4))

    # movement
    Movement_avg <- apply(rep[[i]]$Movement[,,,avg_yrs,,,,drop = FALSE], c(1,2,3,5,6,7), mean) # movement
    Movement <- abind::abind(replicate(n_proj_yrs, Movement_avg, simplify = FALSE), along = 4)
    # Movement / mortality sequencing must follow the fitted model. Under continuous movement
    # the generator has to be averaged on the rate scale, since mean(expm(Q)) != expm(mean(Q)).
    proj_move_timing <- if(is.null(data[[i]]$move_timing)) 0 else data[[i]]$move_timing
    if(proj_move_timing == 2) {
      Mrate_avg <- apply(rep[[i]]$Mrate[,,,avg_yrs,,,,drop = FALSE], c(1,2,3,5,6,7), mean)
      proj_Mrate <- abind::abind(replicate(n_proj_yrs, Mrate_avg, simplify = FALSE), along = 4)
    } else proj_Mrate <- NULL
    sgl_seas_spawning_movement_avg <- apply(rep[[i]]$sgl_seas_spawning_movement[,,,avg_yrs,,1,drop = FALSE], c(1,2,3,5), mean)
    sgl_seas_spawning_movement <- array(sgl_seas_spawning_movement_avg, dim = c(data[[i]]$n_pop, data[[i]]$n_regions, data[[i]]$n_regions, data[[i]]$n_ages)) # Movement
    stray_rate <- array(apply(rep[[i]]$stray_rate[,avg_yrs, drop = FALSE], 1, mean), dim = c(data[[i]]$n_pop, n_proj_yrs))

    # Sex ratio
    sexratio_avg <- apply(rep[[i]]$sexratio[,,avg_yrs,,drop = FALSE], c(1,2,4), mean)
    sexratio <- array(dim = c(data[[i]]$n_pop, data[[i]]$n_regions, proj_model_opt$n_proj_yrs, data[[i]]$n_sexes)) # empty array
    for(yr in 1:proj_model_opt$n_proj_yrs) sexratio[,, yr, ] <- sexratio_avg # populate empty array

    # Now, set up inputs for reference points
    f_ref_pt = array(tmp_ref_pts$f_ref_pt, dim = c(data[[i]]$n_regions, n_proj_yrs))
    b_ref_pt = array(tmp_ref_pts$b_ref_pt, dim = c(data[[i]]$n_pop, data[[i]]$n_regions, n_proj_yrs))

    # Set up stock-recruit options if projecting under a stock-recruit curve.
    # Do_Population_Projection derives the curve from recruitment_opt itself, so
    # the same option list serves both forms.
    if(proj_model_opt$recruitment_opt %in% c('bh_rec', 'ricker_rec')) {
      # Reference year for the biologicals behind unfished spawning biomass per
      # recruit. Must match what the model was fitted with or S0, and therefore
      # the whole curve, shifts. Older data lists predate the option.
      sr_yr <- if(is.null(data[[i]]$SR_ref_yr)) 1 else data[[i]]$SR_ref_yr
      srr_opt <- list(
        do_recruits_move = data[[i]]$do_recruits_move,
        rec_dd = data[[i]]$rec_dd,
        rec_lag = data[[i]]$rec_lag,
        R0 = rep[[i]]$R0,
        rec_region_prop = rep[[i]]$rec_region_prop,
        h = rep[[i]]$h_trans,
        SSB = rep[[i]]$SSB,

        # Demographics for unfished SSB, taken at SR_ref_yr
        WAA = array(data[[i]]$WAA[,,sr_yr,,,1,drop = FALSE], dim = c(data[[i]]$n_pop, data[[i]]$n_regions, data[[i]]$n_seas, length(data[[i]]$ages)) ),
        MatAA = array(data[[i]]$MatAA[,,sr_yr,,,1,drop = FALSE], dim = c(data[[i]]$n_pop, data[[i]]$n_regions, data[[i]]$n_seas, length(data[[i]]$ages)) ) ,
        natmort = array(rep[[i]]$natmort[,,sr_yr,,1,drop = FALSE], dim = c(data[[i]]$n_pop, data[[i]]$n_regions, length(data[[i]]$ages) )),
        sgl_seas_spawning_movement = array(rep[[i]]$sgl_seas_spawning_movement[,,,sr_yr,,1], dim = c(data[[i]]$n_pop, data[[i]]$n_regions, data[[i]]$n_regions, length(data[[i]]$ages) )),
        Movement = array(rep[[i]]$Movement[,,,sr_yr,,,1], dim = c(data[[i]]$n_pop, data[[i]]$n_regions, data[[i]]$n_regions, data[[i]]$n_seas, length(data[[i]]$ages) )),
        # Instantaneous rates matched to Movement above, needed by the SPR machinery
        # behind Beverton-Holt when movement is continuous
        Mrate = if(proj_move_timing == 2) array(rep[[i]]$Mrate[,,,sr_yr,,,1], dim = c(data[[i]]$n_pop, data[[i]]$n_regions, data[[i]]$n_regions, data[[i]]$n_seas, length(data[[i]]$ages) )) else NULL,
        stray_rate = array(rep[[i]]$stray_rate[,1], dim = data[[i]]$n_pop),
        sex_ratio_f = array(if(data[[i]]$n_sexes == 1) 0.5 else rep[[i]]$sexratio[,,1,1], dim = c(data[[i]]$n_pop, data[[i]]$n_regions))
      )
    } else {
      srr_opt <- NULL
    }

    # do population projection
    out_proj <- Do_Population_Projection(n_proj_yrs = n_proj_yrs, # Number of projection years
                                         n_regions = data[[i]]$n_regions, # number of regions
                                         n_ages = length(data[[i]]$ages), # number of ages
                                         n_sexes = data[[i]]$n_sexes, # number of sexes
                                         n_pop = data[[i]]$n_pop, # number of populations
                                         n_seas = data[[i]]$n_seas, # number of seasons
                                         seasdur = data[[i]]$seasdur, # seasonal duration
                                         rec_seas_prop = rep[[i]]$rec_seas_prop, # recruitment seasonal proportion
                                         spawn_seas = data[[i]]$spawn_seas, # spawning season
                                         natal_region = data[[i]]$natal_region, # natal regions
                                         stray_rate = stray_rate, # stray rate
                                         t_spawn = data[[i]]$t_spawn, # spawn timing
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
                                         fmort_opt = 'Input', # Fishing mortality in projection years (whether input or HCR)
                                         srr_opt = srr_opt, # beverton holt projection options
                                         move_timing = proj_move_timing, # movement / mortality sequencing
                                         Mrate = proj_Mrate # instantaneous rates (continuous movement only)
    )

    # extract out quantities and store
    key_quants_tmp <- data.frame(Model = model_names[i],
                                 Region = 1:data[[i]]$n_regions,
                                 Terminal_SSB = round(apply(out_proj$proj_SSB[,,1, drop = FALSE], 2, sum), 5),
                                 Terminal_SSB0 = round(apply(out_proj$proj_Dynamic_SSB0[,,1, drop = FALSE], 2, sum), 5),
                                 Terminal_F = rowSums(terminal_F),
                                 Catch_Advice = round(apply(out_proj$proj_Catch[,,2,,,drop = FALSE], 2, sum), 5), # sum across fleets, season, and populations
                                 B_Ref_Pt = apply(tmp_ref_pts$b_ref_pt, 2, sum),
                                 F_Ref_Pt = round(tmp_ref_pts$f_ref_pt, 5),
                                 B_over_B_Ref = round(apply(out_proj$proj_SSB[,,1, drop = FALSE], 2, sum) / apply(tmp_ref_pts$b_ref_pt, 2, sum), 5),
                                 B_over_DynB_Ref = round(apply(out_proj$proj_SSB[,,1, drop = FALSE], 2, sum) / apply(out_proj$proj_Dynamic_SSB0[,,1, drop = FALSE], 2, sum), 5),
                                 F_over_F_Ref = round(rowSums(terminal_F) / tmp_ref_pts$f_ref_pt, 5)
    )

    key_quants_df <- rbind(key_quants_df, key_quants_tmp)

  } # end i

  # output table
  table_plot <- grid::grid.grabExpr(gridExtra::grid.table(key_quants_df))
  table_plot1 <- cowplot::ggdraw() + cowplot::draw_grob(table_plot)

  return(list(key_quants_df, table_plot1))

} # end function
#' ggplot2 theme for SPoRC plots
#'
#' A clean \code{theme_bw}-based ggplot2 theme with enlarged text elements
#' suitable for publication-quality SPoRC diagnostic and results figures.
#'
#' @return A \code{ggplot2} theme object.
#'
#' @import ggplot2
#' @export theme_sablefish
#' @family Plotting
theme_sablefish <- function() {
   theme_bw() +
    theme(legend.position = "top",
          strip.text = element_text(size = 17),
          title = element_text(size = 21, color = 'black'),
          axis.text = element_text(size = 15, color = "black"),
          axis.title = element_text(size = 17, color = 'black'),
          legend.text = element_text(size = 15, color = "black"),
          legend.title = element_text(size = 17, color = 'black'))
}
