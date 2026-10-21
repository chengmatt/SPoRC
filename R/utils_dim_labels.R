# Shared helpers
#
# Names the dims of a report, parameter list or data list (population, region, year, season, age, sex,
# fleet) on a copy, so the lists the model and every package function read by position stay as they were.

# Labeling ---------------------------------------------------------------------

#' Name the dimensions of a report, parameter list or data list
#'
#' Returns a copy of \code{x} with every array's dims named (\code{Pop},
#' \code{Region}, \code{Year}, \code{Seas}, \code{Age}, \code{Len}, \code{Sex},
#' \code{Fleet}, ...) and labeled with calendar years, model ages and length
#' bins. \code{x} can be the report (\code{obj$rep}), the parameter list
#' (\code{input_list$par} or \code{obj$env$parList()}) or the data list
#' (\code{input_list$data}), and \code{reshape2::melt()} on a labeled array
#' returns those columns by name. Package functions read these arrays by
#' position, so pass them the unlabeled lists.
#'
#' @param x Report, parameter or data list from a SPoRC model.
#' @param data The model's data list, \code{input_list$data}, which gives the
#'   years, ages and length bins. Defaults to \code{x}, for labeling the data
#'   list itself.
#' @param labels Optional named list of labels for \code{Pop}, \code{Region},
#'   \code{Seas}, \code{Sex}, \code{Fish_Fleet} or \code{Srv_Fleet}, one label per
#'   population, region and so on (e.g. \code{list(Sex = c("Female", "Male"))}).
#'   Unlabeled dims are numbered from one.
#'
#' @return \code{x}, with \code{dimnames} on its arrays and \code{names} on its
#'   vectors. Values are unchanged. Scalars, lists, vectors of indices and
#'   arrays of a shape this function does not describe are returned as they
#'   were. Dims that count blocks, parameters or selectivity bins
#'   (\code{Block}, \code{Sel_Par}, \code{Bin}, ...) keep any names they had
#'   and are otherwise numbered from one, as is any dim whose length differs
#'   from its label set,
#'   such as observed ages pooled differently from the model ages. Year dims
#'   past the last model year continue the calendar years, so \code{NAA} ends at
#'   the terminal year plus one.
#'
#' @examples
#' \dontrun{
#' labels <- list(Sex = c("Female", "Male"))
#' rep <- label_dims(obj$rep, input_list$data, labels)
#' par <- label_dims(input_list$par, input_list$data, labels)
#' dat <- label_dims(input_list$data, labels = labels)
#' rep$SSB[1, , "2020"]
#' ssb_df <- reshape2::melt(rep$SSB) # columns Pop, Region, Year, value
#' }
#'
#' @export label_dims
#' @family Utility
label_dims <- function(x,
                       data = x,
                       labels = list()
                       ) {

  if(is.null(data$years) || is.null(data$ages)) stop("data must be the model's data list, input_list$data, which gives the years and ages.")
  n_yrs <- length(data$years)

  # label sets a user can replace, and how many labels each needs
  n_labels <- c(Pop = data$n_pop,
                Region = data$n_regions,
                Seas = data$n_seas,
                Sex = data$n_sexes,
                Fish_Fleet = data$n_fish_fleets,
                Srv_Fleet = data$n_srv_fleets
                )

  # Input Validation ------------------------------------------------------------

  if(!is.list(labels)) stop("labels must be a named list, e.g. list(Region = c('BS', 'AI', 'WGOA'), Sex = c('Female', 'Male')).")
  bad_names <- setdiff(names(labels), names(n_labels))
  if(length(bad_names) > 0 || (length(labels) > 0 && is.null(names(labels)))) {
    stop("labels can only name ", paste(names(n_labels), collapse = ", "),
         ". Years, ages and length bins are read from data. Got: ", paste(bad_names, collapse = ", "))
  }
  for(dim_name in names(labels)) {
    if(length(labels[[dim_name]]) != n_labels[[dim_name]]) {
      stop("labels$", dim_name, " has ", length(labels[[dim_name]]), " entries but the model has ",
           n_labels[[dim_name]], ". Give one label per ", dim_name, ".")
    }
  } # end dim_name loop

  # Label Sets ------------------------------------------------------------------

  # numbered from one unless the user named them
  label_sets <- list()
  for(dim_name in names(n_labels)) {
    label_sets[[dim_name]] <- if(is.null(labels[[dim_name]])) seq_len(n_labels[[dim_name]]) else labels[[dim_name]]
  } # end dim_name loop

  label_sets$Year <- data$years # calendar years, continued below for longer year dims
  label_sets$Rec_Year <- data$years
  label_sets$Age <- data$ages
  label_sets$Len <- data$lens
  label_sets$Series <- data$dsem_var_names # NULL without a dsem, so numbered
  label_sets$Covariate <- data$dsem_var_names[data$dsem_cov_var_idx]

  # dims that read another dim's labels
  label_sets$Region_From <- label_sets$Region
  label_sets$Region_To <- label_sets$Region
  label_sets$Obs_Age <- label_sets$Age # observed ages are the model ages under identity ageing error
  label_sets$Obs_Len <- label_sets$Len
  label_sets$CAAL_Len <- label_sets$Len
  label_sets$Init_Age <- label_sets$Age[-1] # initial age deviations start at the second age

  # names the dims carry when the array has no other dim of that name
  dim_display <- c(Fish_Fleet = "Fleet", Srv_Fleet = "Fleet", Rec_Year = "Year", Init_Age = "Age",
                   Obs_Age = "Age", Obs_Len = "Len", CAAL_Len = "Len")

  # Name Dims -------------------------------------------------------------------

  for(x_name in intersect(names(x), names(dim_roles))) {

    # extract stuff out from dim and x
    x_dim <- x[[x_name]]
    roles <- dim_roles[[x_name]]
    extents <- if(is.null(dim(x_dim))) length(x_dim) else dim(x_dim)
    if(is.null(x_dim) || !is.atomic(x_dim) || length(extents) != length(roles)) next # a shape the table does not describe is left as it is
    old_labels <- if(is.null(dim(x_dim))) list(names(x_dim)) else dimnames(x_dim)

    # make dim labels container
    dim_labels <- vector("list", length(roles))
    for(i in seq_along(roles)) {

      dim_lab <- label_sets[[roles[i]]] # NULL for blocks, parameters and bins, so numbered

      # projection years and the terminal year plus one continue the calendar years
      if(roles[i] %in% c("Year", "Rec_Year") && extents[i] > n_yrs) dim_lab <- c(dim_lab, dim_lab[n_yrs] + seq_len(extents[i] - n_yrs))
      if(roles[i] == "Rec_Year" && extents[i] < n_yrs) dim_lab <- dim_lab[seq_len(extents[i])] # recruitment deviations start in the first year

      if(length(dim_lab) != extents[i] || anyNA(dim_lab)) dim_lab <- seq_len(extents[i]) # number it when the labels do not fit
      dim_labels[[i]] <- as.character(dim_lab)
      if(is.null(label_sets[[roles[i]]]) && !is.null(old_labels[[i]])) dim_labels[[i]] <- old_labels[[i]] # column names and the like are kept

    } # end i loop

    # Obs_Age reads as Age, unless the array has model ages too (e.g., the ageing error matrix)
    dim_names <- roles
    renamed <- roles %in% names(dim_display)
    dim_names[renamed] <- dim_display[roles[renamed]]
    clash <- dim_names %in% dim_names[duplicated(dim_names)] & renamed
    dim_names[clash] <- roles[clash]
    names(dim_labels) <- dim_names

    if(is.null(dim(x_dim))) names(x_dim) <- dim_labels[[1]] # vectors take names
    else dimnames(x_dim) <- dim_labels
    x[[x_name]] <- x_dim

  } # end x_name loop

  return(x)

} # end function


# Dim Roles --------------------------------------------------------------------

# what each dim of an array is, by name, across the report, parameter and data lists. a name in
# more than one list has the same dims in each. Bin, Block and the *_Par dims are positions
dim_roles <- list(

  # Report ---------------------------------------------------------------------

  # recruitment
  R0 = "Pop",
  rinit = "Pop",
  sr_R0 = "Pop",
  R0_yr = c("Pop", "Year"),
  sexratio = c("Pop", "Region", "Year", "Sex"),
  rec_region_prop = c("Pop", "Region"),
  rec_seas_prop = c("Pop", "Seas"),
  stray_rate = c("Pop", "Year"),
  h_trans = c("Pop", "Region"),
  SR_pred = c("Pop", "Region", "Year"),
  Rec = c("Pop", "Region", "Year"),
  ln_RecDevs = c("Pop", "Region", "Rec_Year"), # the first deviation is the first model year
  bias_ramp = "Rec_Year",
  init_bias_ramp = "Init_Age",

  # numbers at age and mortality
  NAA = c("Pop", "Region", "Year", "Seas", "Age", "Sex"), # years run to the terminal year plus one
  NAA0 = c("Pop", "Region", "Year", "Seas", "Age", "Sex"),
  NAA_bef = c("Pop", "Region", "Year", "Seas", "Age", "Sex"),
  NAA_aft = c("Pop", "Region", "Year", "Seas", "Age", "Sex"),
  NAA_int = c("Pop", "Region", "Year", "Seas", "Age", "Sex"),
  NAA_pred = c("Pop", "Region", "Year", "Seas", "Age", "Sex"),
  NAA_scalar = c("Pop", "Region", "Year", "Seas", "Age", "Sex"),
  ZAA = c("Pop", "Region", "Year", "Seas", "Age", "Sex"),
  natmort = c("Pop", "Region", "Year", "Seas", "Age", "Sex"),

  # movement
  Movement = c("Pop", "Region_From", "Region_To", "Year", "Seas", "Age", "Sex"),
  Mrate = c("Pop", "Region_From", "Region_To", "Year", "Seas", "Age", "Sex"),
  sgl_seas_spawning_movement = c("Pop", "Region_From", "Region_To", "Year", "Age", "Sex"),
  move_devs = c("Pop", "Region_From", "Move_Par", "Year", "Seas", "Age", "Sex"), # third dim is not a destination region

  # growth and weight
  SizeAgeTrans_fish = c("Pop", "Region", "Year", "Seas", "Len", "Age", "Sex", "Fish_Fleet"),
  SizeAgeTrans_srv = c("Pop", "Region", "Year", "Seas", "Len", "Age", "Sex", "Srv_Fleet"),
  SizeAgeTrans_spawn = c("Pop", "Region", "Year", "Len", "Age", "Sex"),
  mean_LAA_fish = c("Pop", "Region", "Year", "Seas", "Age", "Sex", "Fish_Fleet"),
  sd_LAA_fish = c("Pop", "Region", "Year", "Seas", "Age", "Sex", "Fish_Fleet"),
  mean_LAA_srv = c("Pop", "Region", "Year", "Seas", "Age", "Sex", "Srv_Fleet"),
  sd_LAA_srv = c("Pop", "Region", "Year", "Seas", "Age", "Sex", "Srv_Fleet"),
  mean_LAA_spawn = c("Pop", "Region", "Year", "Seas", "Age", "Sex"),
  sd_LAA_spawn = c("Pop", "Region", "Year", "Seas", "Age", "Sex"),
  L_beg = c("Pop", "Region", "Year", "Age", "Sex"),
  Linf = c("Pop", "Region", "Year", "Sex"),
  growth_pars_y = c("Pop", "Region", "Year", "Growth_Par", "Sex"),
  WAA = c("Pop", "Region", "Year", "Seas", "Age", "Sex"),
  WAA_fish = c("Pop", "Region", "Year", "Seas", "Age", "Sex", "Fish_Fleet"),
  WAA_srv = c("Pop", "Region", "Year", "Seas", "Age", "Sex", "Srv_Fleet"),

  # fishery
  init_F = c("Region", "Seas", "Fish_Fleet"),
  Fmort = c("Region", "Year", "Seas", "Fish_Fleet"),
  dmr = c("Region", "Year", "Seas", "Fish_Fleet"),
  tot_FAA = c("Pop", "Region", "Year", "Seas", "Age", "Sex", "Fish_Fleet"),
  ret_FAA = c("Pop", "Region", "Year", "Seas", "Age", "Sex", "Fish_Fleet"),
  disc_FAA = c("Pop", "Region", "Year", "Seas", "Age", "Sex", "Fish_Fleet"),
  CAA = c("Pop", "Region", "Year", "Seas", "Age", "Sex", "Fish_Fleet"),
  DAA = c("Pop", "Region", "Year", "Seas", "Age", "Sex", "Fish_Fleet"),
  CAL = c("Pop", "Region", "Year", "Seas", "Len", "Sex", "Fish_Fleet"),
  DAL = c("Pop", "Region", "Year", "Seas", "Len", "Sex", "Fish_Fleet"),
  Fish_caal = c("Pop", "Region", "Year", "Seas", "Len", "Age", "Sex", "Fish_Fleet"),
  Fish_caal_discard = c("Pop", "Region", "Year", "Seas", "Len", "Age", "Sex", "Fish_Fleet"),
  PredCatch = c("Pop", "Region", "Year", "Seas", "Fish_Fleet"),
  PredDiscard = c("Pop", "Region", "Year", "Seas", "Fish_Fleet"),
  PredFishIdx = c("Pop", "Region", "Year", "Seas", "Fish_Fleet"),
  fish_sel = c("Pop", "Region", "Year", "Seas", "Age", "Sex", "Fish_Fleet"),
  ret_sel = c("Pop", "Region", "Year", "Seas", "Age", "Sex", "Fish_Fleet"),
  fish_sel_l = c("Region", "Year", "Len", "Sex", "Fish_Fleet"),
  ret_sel_l = c("Region", "Year", "Len", "Sex", "Fish_Fleet"),
  fish_q = c("Region", "Year", "Fish_Fleet"),
  ln_fishsel_devs = c("Region", "Year", "Bin", "Sex", "Fish_Fleet"), # bins are ages, lengths or parameters by deviation form

  # survey
  PredSrvIdx = c("Pop", "Region", "Year", "Seas", "Srv_Fleet"),
  srv_sel = c("Pop", "Region", "Year", "Seas", "Age", "Sex", "Srv_Fleet"),
  srv_sel_l = c("Region", "Year", "Len", "Sex", "Srv_Fleet"),
  srv_q = c("Region", "Year", "Srv_Fleet"),
  SrvIAA = c("Pop", "Region", "Year", "Seas", "Age", "Sex", "Srv_Fleet"),
  SrvIAL = c("Pop", "Region", "Year", "Seas", "Len", "Sex", "Srv_Fleet"),
  Srv_caal = c("Pop", "Region", "Year", "Seas", "Len", "Age", "Sex", "Srv_Fleet"),
  ln_srvsel_devs = c("Region", "Year", "Bin", "Sex", "Srv_Fleet"),

  # observation error
  ln_sigmaC = c("Region", "Year", "Seas", "Fish_Fleet"),
  ln_sigmaC_pop = c("Pop", "Region", "Year", "Seas", "Fish_Fleet"),
  ln_sigmaFishIdx = "Fish_Fleet",
  ln_sigmaFishIdx_pop = "Fish_Fleet",
  ln_sigmaSrvIdx = "Srv_Fleet",
  ln_sigmaSrvIdx_pop = "Srv_Fleet",
  FishIdx_SD = c("Region", "Year", "Seas", "Fish_Fleet"),
  FishIdx_pop_SD = c("Pop", "Region", "Year", "Seas", "Fish_Fleet"),
  SrvIdx_SD = c("Region", "Year", "Seas", "Srv_Fleet"),
  SrvIdx_pop_SD = c("Pop", "Region", "Year", "Seas", "Srv_Fleet"),
  ln_sigmaCAA = c("Obs_Age", "Sex", "Fish_Fleet"),
  ln_sigmaSrvIdxAA = c("Obs_Age", "Sex", "Srv_Fleet"),
  trans_rho_catch = c("Region", "Sex", "Fish_Fleet"),
  trans_rho_discard = c("Region", "Sex", "Fish_Fleet"),
  trans_rho_srv_idx = c("Region", "Sex", "Srv_Fleet"),

  # at-age observations, on the observed ages
  PredCatchAA = c("Region", "Year", "Seas", "Obs_Age", "Sex", "Fish_Fleet"),
  PredDiscardAA = c("Region", "Year", "Seas", "Obs_Age", "Sex", "Fish_Fleet"),
  PredSrvIdxAA = c("Region", "Year", "Seas", "Obs_Age", "Sex", "Srv_Fleet"),
  PredCatchAA_pop = c("Pop", "Region", "Year", "Seas", "Obs_Age", "Sex", "Fish_Fleet"),
  PredDiscardAA_pop = c("Pop", "Region", "Year", "Seas", "Obs_Age", "Sex", "Fish_Fleet"),
  PredSrvIdxAA_pop = c("Pop", "Region", "Year", "Seas", "Obs_Age", "Sex", "Srv_Fleet"),
  CatchAA_nLL = c("Region", "Year", "Seas", "Obs_Age", "Sex", "Fish_Fleet"),
  DiscardAA_nLL = c("Region", "Year", "Seas", "Obs_Age", "Sex", "Fish_Fleet"),
  SrvIdxAA_nLL = c("Region", "Year", "Seas", "Obs_Age", "Sex", "Srv_Fleet"),
  CatchAA_pop_nLL = c("Pop", "Region", "Year", "Seas", "Obs_Age", "Sex", "Fish_Fleet"),
  DiscardAA_pop_nLL = c("Pop", "Region", "Year", "Seas", "Obs_Age", "Sex", "Fish_Fleet"),
  SrvIdxAA_pop_nLL = c("Pop", "Region", "Year", "Seas", "Obs_Age", "Sex", "Srv_Fleet"),

  # tagging
  pred_conv_tag_fish_recap = c("Liberty", "Seas", "Tag_Cohort", "Pop", "Region", "Age", "Sex", "Fish_Fleet"),
  conv_tag_fish_avail = c("Liberty", "Seas", "Tag_Cohort", "Pop", "Region", "Age", "Sex"), # liberty starts at release
  conv_tag_fish_reporting = c("Region", "Year", "Fish_Fleet"),

  # data likelihoods
  Catch_nLL = c("Region", "Year", "Seas", "Fish_Fleet"),
  Discard_nLL = c("Region", "Year", "Seas", "Fish_Fleet"),
  FishIdx_nLL = c("Region", "Year", "Seas", "Fish_Fleet"),
  SrvIdx_nLL = c("Region", "Year", "Seas", "Srv_Fleet"),
  Catch_pop_nLL = c("Pop", "Region", "Year", "Seas", "Fish_Fleet"),
  Discard_pop_nLL = c("Pop", "Region", "Year", "Seas", "Fish_Fleet"),
  FishIdx_pop_nLL = c("Pop", "Region", "Year", "Seas", "Fish_Fleet"),
  SrvIdx_pop_nLL = c("Pop", "Region", "Year", "Seas", "Srv_Fleet"),
  FishAgeComps_nLL = c("Region", "Year", "Seas", "Sex", "Fish_Fleet"),
  FishLenComps_nLL = c("Region", "Year", "Seas", "Sex", "Fish_Fleet"),
  FishAgeComps_discard_nLL = c("Region", "Year", "Seas", "Sex", "Fish_Fleet"),
  FishLenComps_discard_nLL = c("Region", "Year", "Seas", "Sex", "Fish_Fleet"),
  SrvAgeComps_nLL = c("Region", "Year", "Seas", "Sex", "Srv_Fleet"),
  SrvLenComps_nLL = c("Region", "Year", "Seas", "Sex", "Srv_Fleet"),
  FishAgeComps_pop_nLL = c("Pop", "Region", "Year", "Seas", "Sex", "Fish_Fleet"),
  FishLenComps_pop_nLL = c("Pop", "Region", "Year", "Seas", "Sex", "Fish_Fleet"),
  FishAgeComps_discard_pop_nLL = c("Pop", "Region", "Year", "Seas", "Sex", "Fish_Fleet"),
  FishLenComps_discard_pop_nLL = c("Pop", "Region", "Year", "Seas", "Sex", "Fish_Fleet"),
  SrvAgeComps_pop_nLL = c("Pop", "Region", "Year", "Seas", "Sex", "Srv_Fleet"),
  SrvLenComps_pop_nLL = c("Pop", "Region", "Year", "Seas", "Sex", "Srv_Fleet"),
  Fish_caal_nLL = c("Region", "Year", "Seas", "CAAL_Len", "Sex", "Fish_Fleet"),
  Srv_caal_nLL = c("Region", "Year", "Seas", "CAAL_Len", "Sex", "Srv_Fleet"),
  conv_fish_tag_nLL = c("Liberty", "Seas", "Tag_Cohort", "Region", "Fish_Fleet"),

  # priors and penalties
  Fmort_nLL = c("Region", "Year", "Seas", "Fish_Fleet"),
  dmr_nLL = c("Region", "Year", "Seas", "Fish_Fleet"),
  Rec_nLL = c("Pop", "Region", "Rec_Year"),
  Rec_level_nLL = c("Pop", "Region", "Year"),
  SR_pen_nLL = c("Pop", "Region", "Year"),
  Init_Rec_nLL = c("Pop", "Region", "Init_Age", "Sex"),
  Init_Sex_nLL = c("Pop", "Region", "Init_Age", "Sex"),

  # derived quantities
  Total_Biom = c("Pop", "Region", "Year"),
  SSB = c("Pop", "Region", "Year"),
  Dynamic_SSB0 = c("Pop", "Region", "Year"),
  eff_SSB = c("Pop", "Year"),
  Aggregated_SSB = "Year",
  Dynamic_Aggregated_SSB0 = "Year",

  # dsem
  dsem_x_grid = c("Year", "Series"), # grid rows start at the first model year
  dsem_margvar_grid = c("Year", "Series"),

  # Parameters -----------------------------------------------------------------

  # recruitment
  ln_global_R0 = c("Pop", "Block"),
  ln_sr_R0 = "Pop",
  ln_rinit = "Pop",
  rec_region_prop_pars = c("Pop", "Region_Par"),
  rec_seas_prop_pars = c("Pop", "Seas_Par"),
  steepness_h = c("Pop", "Region"),
  ln_InitDevs = c("Pop", "Region", "Init_Age", "Sex"),
  ln_sigmaR = c("Period", "Pop", "Region"),
  sexratio_pars = c("Pop", "Region", "Block"),
  RecDevs_rho = c("Pop", "Region"),
  stray_rate_pars = c("Pop", "Block"),

  # natural mortality and numbers at age
  ln_M = c("Pop_Block", "Region_Block", "Year_Block", "Seas_Block", "Age_Block", "Sex_Block"),
  ln_NAA = c("Pop", "Region", "Year", "Seas", "Age", "Sex"),
  ln_sigmaNAA = c("Pop_Block", "Region_Block", "Year_Block", "Seas_Block", "Age_Block", "Sex_Block"),
  NAA_pe_pars = c("Pop", "Region", "PE_Par", "Sex"),
  NAA_region_corr_pars = c("Pop", "Corr_Par", "Sex"),
  NAA_season_corr_pars = c("Pop", "Corr_Par", "Sex"),
  NAA_pop_corr_pars = "Corr_Par",

  # movement
  move_pars = c("Pop", "Region_From", "Move_Par", "Year", "Seas", "Age", "Sex"),
  move_pe_pars = c("Region_From", "Move_Par", "PE_Par"),
  move_pop_corr_pars = "Corr_Par",
  move_seas_corr_pars = "Corr_Par",

  # tagging
  ln_init_conv_tag_mort = "Tag_Cohort",
  ln_conv_tag_shed = "Tag_Cohort",
  conv_tag_fish_reporting_pars = c("Region", "Block", "Fish_Fleet"),

  # growth
  ln_growth_pars = c("Pop", "Region", "Sex", "Growth_Par"),
  ln_growth_devs = c("Pop", "Region", "Year", "Growth_Par", "Sex"),
  ln_growth_semipar_devs = c("Pop", "Region", "Year", "Age", "Sex"),
  growth_pe_pars = c("Pop", "Region", "PE_Par", "Sex", "PE_Set"),

  # selectivity and catchability
  ln_fishsel_bin_devs = c("Region", "Year", "Bin", "Sex", "Fish_Fleet"),
  fishsel_bin_devs_pe_pars = c("Region", "Bin", "Sex", "Fish_Fleet"),
  fish_fixed_sel_pars = c("Region", "Sel_Par", "Block", "Sex", "Fish_Fleet"),
  ln_fish_q = c("Region", "Block", "Fish_Fleet"),
  fishsel_pe_pars = c("Region", "PE_Par", "Sex", "Fish_Fleet"),
  ln_fishsel_sex_scale = c("Region", "Block", "Sex", "Fish_Fleet"),
  ln_fish_q_devs = c("Region", "Year", "Fish_Fleet"),
  ln_sigma_fish_q = c("Region", "Fish_Fleet"),
  fish_q_rho = c("Region", "Fish_Fleet"),
  ln_retsel_bin_devs = c("Region", "Year", "Bin", "Sex", "Fish_Fleet"),
  retsel_bin_devs_pe_pars = c("Region", "Bin", "Sex", "Fish_Fleet"),
  ret_fixed_sel_pars = c("Region", "Sel_Par", "Block", "Sex", "Fish_Fleet"),
  retsel_pe_pars = c("Region", "PE_Par", "Sex", "Fish_Fleet"),
  ln_retsel_devs = c("Region", "Year", "Bin", "Sex", "Fish_Fleet"),
  ln_retsel_sex_scale = c("Region", "Block", "Sex", "Fish_Fleet"),
  ln_srvsel_bin_devs = c("Region", "Year", "Bin", "Sex", "Srv_Fleet"),
  srvsel_bin_devs_pe_pars = c("Region", "Bin", "Sex", "Srv_Fleet"),
  srv_fixed_sel_pars = c("Region", "Sel_Par", "Block", "Sex", "Srv_Fleet"),
  ln_srv_q = c("Region", "Block", "Srv_Fleet"),
  srvsel_pe_pars = c("Region", "PE_Par", "Sex", "Srv_Fleet"),
  ln_srvsel_sex_scale = c("Region", "Block", "Sex", "Srv_Fleet"),
  ln_srv_q_devs = c("Region", "Year", "Srv_Fleet"),
  ln_sigma_srv_q = c("Region", "Srv_Fleet"),
  srv_q_rho = c("Region", "Srv_Fleet"),

  # fishing and discard mortality
  init_F_par = c("Region", "Seas", "Fish_Fleet"),
  ln_sigmaF = c("Region", "Seas", "Fish_Fleet"),
  Fdev_rho = c("Region", "Seas", "Fish_Fleet"),
  ln_F_mean = c("Region", "Seas", "Fish_Fleet"),
  ln_F_devs = c("Region", "Year", "Seas", "Fish_Fleet"),
  ln_sigma_dmr = c("Region", "Seas", "Fish_Fleet"),
  logit_dmr_mean = c("Region", "Seas", "Fish_Fleet"),
  logit_dmr_devs = c("Region", "Year", "Seas", "Fish_Fleet"),

  # observation error
  trans_rho_catch_year = c("Region", "Sex", "Fish_Fleet"),
  trans_rho_catch_us = c("Corr_Par", "Region", "Sex", "Fish_Fleet"),
  trans_rho_catch_pop = c("Pop", "Region", "Sex", "Fish_Fleet"),
  trans_rho_catch_pop_year = c("Pop", "Region", "Sex", "Fish_Fleet"),
  trans_rho_catch_pop_us = c("Corr_Par", "Pop", "Region", "Sex", "Fish_Fleet"),
  trans_rho_discard_year = c("Region", "Sex", "Fish_Fleet"),
  trans_rho_discard_us = c("Corr_Par", "Region", "Sex", "Fish_Fleet"),
  trans_rho_discard_pop = c("Pop", "Region", "Sex", "Fish_Fleet"),
  trans_rho_discard_pop_year = c("Pop", "Region", "Sex", "Fish_Fleet"),
  trans_rho_discard_pop_us = c("Corr_Par", "Pop", "Region", "Sex", "Fish_Fleet"),
  ln_sigmaCAA_pop = c("Pop", "Obs_Age", "Sex", "Fish_Fleet"),
  ln_sigmaDAA = c("Obs_Age", "Sex", "Fish_Fleet"),
  ln_sigmaDAA_pop = c("Pop", "Obs_Age", "Sex", "Fish_Fleet"),
  ln_sigmaD = c("Region", "Year", "Seas", "Fish_Fleet"),
  ln_sigmaD_pop = c("Pop", "Region", "Year", "Seas", "Fish_Fleet"),
  ln_FishAge_theta = c("Region", "Sex", "Fish_Fleet"),
  FishAge_corr_pars = c("Region", "Sex", "Fish_Fleet", "Corr_Par"),
  ln_FishAge_theta_agg = "Fish_Fleet",
  FishAge_corr_pars_agg = "Fish_Fleet",
  ln_FishLen_theta = c("Region", "Sex", "Fish_Fleet"),
  FishLen_corr_pars = c("Region", "Sex", "Fish_Fleet", "Corr_Par"),
  ln_FishLen_theta_agg = "Fish_Fleet",
  FishLen_corr_pars_agg = "Fish_Fleet",
  ln_FishAge_pop_theta = c("Pop", "Region", "Sex", "Fish_Fleet"),
  FishAge_pop_corr_pars = c("Pop", "Region", "Sex", "Fish_Fleet", "Corr_Par"),
  ln_FishAge_pop_theta_agg = c("Pop", "Fish_Fleet"),
  FishAge_pop_corr_pars_agg = c("Pop", "Fish_Fleet"),
  ln_FishLen_pop_theta = c("Pop", "Region", "Sex", "Fish_Fleet"),
  FishLen_pop_corr_pars = c("Pop", "Region", "Sex", "Fish_Fleet", "Corr_Par"),
  ln_FishLen_pop_theta_agg = c("Pop", "Fish_Fleet"),
  FishLen_pop_corr_pars_agg = c("Pop", "Fish_Fleet"),
  ln_FishAge_discard_theta = c("Region", "Sex", "Fish_Fleet"),
  FishAge_discard_corr_pars = c("Region", "Sex", "Fish_Fleet", "Corr_Par"),
  ln_FishAge_discard_theta_agg = "Fish_Fleet",
  FishAge_discard_corr_pars_agg = "Fish_Fleet",
  ln_FishLen_discard_theta = c("Region", "Sex", "Fish_Fleet"),
  FishLen_discard_corr_pars = c("Region", "Sex", "Fish_Fleet", "Corr_Par"),
  ln_FishLen_discard_theta_agg = "Fish_Fleet",
  FishLen_discard_corr_pars_agg = "Fish_Fleet",
  ln_FishAge_discard_pop_theta = c("Pop", "Region", "Sex", "Fish_Fleet"),
  FishAge_discard_pop_corr_pars = c("Pop", "Region", "Sex", "Fish_Fleet", "Corr_Par"),
  ln_FishAge_discard_pop_theta_agg = c("Pop", "Fish_Fleet"),
  FishAge_discard_pop_corr_pars_agg = c("Pop", "Fish_Fleet"),
  ln_FishLen_discard_pop_theta = c("Pop", "Region", "Sex", "Fish_Fleet"),
  FishLen_discard_pop_corr_pars = c("Pop", "Region", "Sex", "Fish_Fleet", "Corr_Par"),
  ln_FishLen_discard_pop_theta_agg = c("Pop", "Fish_Fleet"),
  FishLen_discard_pop_corr_pars_agg = c("Pop", "Fish_Fleet"),
  ln_Fish_caal_theta = c("Region", "Sex", "Fish_Fleet"),
  ln_Fish_caal_theta_agg = "Fish_Fleet",
  trans_rho_srv_idx_year = c("Region", "Sex", "Srv_Fleet"),
  trans_rho_srv_idx_us = c("Corr_Par", "Region", "Sex", "Srv_Fleet"),
  trans_rho_srv_idx_pop = c("Pop", "Region", "Sex", "Srv_Fleet"),
  trans_rho_srv_idx_pop_year = c("Pop", "Region", "Sex", "Srv_Fleet"),
  trans_rho_srv_idx_pop_us = c("Corr_Par", "Pop", "Region", "Sex", "Srv_Fleet"),
  ln_sigmaSrvIdxAA_pop = c("Pop", "Obs_Age", "Sex", "Srv_Fleet"),
  ln_SrvAge_theta = c("Region", "Sex", "Srv_Fleet"),
  SrvAge_corr_pars = c("Region", "Sex", "Srv_Fleet", "Corr_Par"),
  ln_SrvAge_theta_agg = "Srv_Fleet",
  SrvAge_corr_pars_agg = "Srv_Fleet",
  ln_SrvLen_theta = c("Region", "Sex", "Srv_Fleet"),
  SrvLen_corr_pars = c("Region", "Sex", "Srv_Fleet", "Corr_Par"),
  ln_SrvLen_theta_agg = "Srv_Fleet",
  SrvLen_corr_pars_agg = "Srv_Fleet",
  ln_SrvAge_pop_theta = c("Pop", "Region", "Sex", "Srv_Fleet"),
  SrvAge_pop_corr_pars = c("Pop", "Region", "Sex", "Srv_Fleet", "Corr_Par"),
  ln_SrvAge_pop_theta_agg = c("Pop", "Srv_Fleet"),
  SrvAge_pop_corr_pars_agg = c("Pop", "Srv_Fleet"),
  ln_SrvLen_pop_theta = c("Pop", "Region", "Sex", "Srv_Fleet"),
  SrvLen_pop_corr_pars = c("Pop", "Region", "Sex", "Srv_Fleet", "Corr_Par"),
  ln_SrvLen_pop_theta_agg = c("Pop", "Srv_Fleet"),
  SrvLen_pop_corr_pars_agg = c("Pop", "Srv_Fleet"),
  ln_Srv_caal_theta = c("Region", "Sex", "Srv_Fleet"),
  ln_Srv_caal_theta_agg = "Srv_Fleet",

  # dsem
  dsem_mu = "Series",
  dsem_x = c("Year", "Series"),
  dsem_delta0 = "Series",

  # Data -----------------------------------------------------------------------

  # dsem
  dsem_cov_obs = c("Year", "Covariate"),
  dsem_cov_fixed_sd = c("Year", "Covariate"),
  dsem_x_known = c("Year", "Series"),

  # weights
  Wt_Catch = c("Region", "Year", "Seas", "Fish_Fleet"),
  Wt_FishAgeComps = c("Region", "Year", "Seas", "Sex", "Fish_Fleet"),
  Wt_SrvAgeComps = c("Region", "Year", "Seas", "Sex", "Srv_Fleet"),
  Wt_FishLenComps = c("Region", "Year", "Seas", "Sex", "Fish_Fleet"),
  Wt_SrvLenComps = c("Region", "Year", "Seas", "Sex", "Srv_Fleet"),
  Wt_FishAgeComps_pop = c("Pop", "Region", "Year", "Seas", "Sex", "Fish_Fleet"),
  Wt_SrvAgeComps_pop = c("Pop", "Region", "Year", "Seas", "Sex", "Srv_Fleet"),
  Wt_FishLenComps_pop = c("Pop", "Region", "Year", "Seas", "Sex", "Fish_Fleet"),
  Wt_SrvLenComps_pop = c("Pop", "Region", "Year", "Seas", "Sex", "Srv_Fleet"),
  Wt_FishAgeComps_discard = c("Region", "Year", "Seas", "Sex", "Fish_Fleet"),
  Wt_FishLenComps_discard = c("Region", "Year", "Seas", "Sex", "Fish_Fleet"),
  Wt_FishAgeComps_discard_pop = c("Pop", "Region", "Year", "Seas", "Sex", "Fish_Fleet"),
  Wt_FishLenComps_discard_pop = c("Pop", "Region", "Year", "Seas", "Sex", "Fish_Fleet"),
  Wt_Fish_caal = c("Region", "Year", "Seas", "CAAL_Len", "Sex", "Fish_Fleet"),
  Wt_Srv_caal = c("Region", "Year", "Seas", "CAAL_Len", "Sex", "Srv_Fleet"),
  Wt_FishIdx = c("Region", "Year", "Seas", "Fish_Fleet"),
  Wt_Discard = c("Region", "Year", "Seas", "Fish_Fleet"),
  Wt_SrvIdx = c("Region", "Year", "Seas", "Srv_Fleet"),
  Wt_Catch_pop = c("Pop", "Region", "Year", "Seas", "Fish_Fleet"),
  Wt_FishIdx_pop = c("Pop", "Region", "Year", "Seas", "Fish_Fleet"),
  Wt_Discard_pop = c("Pop", "Region", "Year", "Seas", "Fish_Fleet"),
  Wt_SrvIdx_pop = c("Pop", "Region", "Year", "Seas", "Srv_Fleet"),
  Wt_Rec = c("Pop", "Region", "Rec_Year"),
  Wt_Init_Rec = c("Pop", "Region", "Init_Age", "Sex"),

  # recruitment
  natal_region = "Pop",
  seasdur = "Seas",
  rec_level_pen_yrs = "Year",
  sr_pen_yrs = "Year",
  init_F_prop = c("Region", "Seas", "Fish_Fleet"),
  fixed_stray_rate = c("Pop", "Year"),
  stray_rate_blocks = c("Pop", "Year"),
  sexratio_blocks = c("Pop", "Region", "Year"),
  fixed_rec_seas_prop = c("Pop", "Seas"),
  R0_blocks = c("Set", "Year", "Pop"),
  init_devs_pen_use = c("Pop", "Region", "Init_Age", "Sex"),
  map_ln_RecDevs = c("Pop", "Region", "Rec_Year"),
  map_ln_InitDevs = c("Pop", "Region", "Init_Age", "Sex"),

  # biology and growth
  MatAA = c("Pop", "Region", "Year", "Seas", "Age", "Sex"),
  AgeingError = c("Year", "Age", "Obs_Age"),
  AgeingError_fish = c("Year", "Age", "Obs_Age", "Fish_Fleet"),
  AgeingError_srv = c("Year", "Age", "Obs_Age", "Srv_Fleet"),
  SizeAgeTrans = c("Pop", "Region", "Year", "Seas", "Len", "Age", "Sex"),
  Fixed_natmort = c("Pop", "Region", "Year", "Seas", "Age", "Sex"),
  M_blocks = c("Pop", "Region", "Year", "Seas", "Age", "Sex"),
  map_ln_NAA = c("Pop", "Region", "Year", "Seas", "Age", "Sex"),
  naa_re_where = c("Pop", "Region"),
  naa_sigma_blocks = c("Pop", "Region", "Year", "Seas", "Age", "Sex"),
  growth_tv_dsem = "Growth_Par",
  growth_len_lower = "Len",
  wt_len_pars = c("Pop", "Region", "Sex", "Wt_Par"),
  growth_tv_model = "Growth_Par",
  growth_par_bounds = c("Growth_Par", "Bound"),
  map_ln_growth_devs = c("Pop", "Region", "Year", "Growth_Par", "Sex"),
  map_ln_growth_semipar_devs = c("Pop", "Region", "Year", "Age", "Sex"),
  LenBinMap = c("Len", "Obs_Len"),
  CAAL_LenBinMap = c("Len", "CAAL_Len"),

  # movement
  Fixed_Movement = c("Pop", "Region_From", "Region_To", "Year", "Seas", "Age", "Sex"),
  adjacency_mat = c("Region_From", "Region_To"),
  adjacency_collapsed = c("Region_From", "Move_Par"),
  area_r = "Region",
  map_move_devs = c("Pop", "Region_From", "Move_Par", "Year", "Seas", "Age", "Sex"),
  move_pop_block = "Pop",
  move_year_block = "Year",
  move_seas_block = "Seas",
  move_age_block = "Age",
  move_sex_block = "Sex",

  # tagging
  conv_tagged_fish = c("Tag_Cohort", "Pop", "Age", "Sex"),
  obs_conv_tag_fish_recap = c("Liberty", "Seas", "Tag_Cohort", "Pop", "Region", "Age", "Sex", "Fish_Fleet"),
  conv_tag_t_tagging = "Tag_Cohort",
  conv_tag_fish_reporting_blocks = c("Region", "Year", "Fish_Fleet"),
  conv_fish_tag_attr = "Tag_Cohort",
  conv_tag_release_indicator = c("Tag_Cohort", "Column"),
  conv_tag_release_platform = c("Tag_Cohort", "Column"),

  # selectivity and catchability
  fish_len_comp_sel = "Fish_Fleet",
  fish_waa_selected = "Fish_Fleet",
  srv_len_comp_sel = "Srv_Fleet",
  srv_waa_selected = "Srv_Fleet",
  cont_tv_fish_sel = c("Region", "Fish_Fleet"),
  fish_sel_blocks = c("Region", "Year", "Fish_Fleet"),
  fish_sel_model = c("Region", "Year", "Fish_Fleet"),
  fish_sel_bicubic_binnodes = c("Region", "Year", "Fish_Fleet"),
  fish_sel_bicubic_yrnodes = c("Region", "Year", "Fish_Fleet"),
  fish_sel_bicubic_selstyr = c("Region", "Year", "Fish_Fleet"),
  fish_sel_bicubic_nselbins = c("Region", "Year", "Fish_Fleet"),
  fish_q_type = "Fish_Fleet",
  fish_q_blocks = c("Region", "Year", "Fish_Fleet"),
  fishsel_pe_wt = "Fish_Fleet",
  fishsel_rw_init_sigma = "Fish_Fleet",
  fish_sel_bin_dev_bins = c("Bin", "Fish_Fleet"),
  cont_tv_fishsel_bin_devs = "Fish_Fleet",
  map_ln_fishsel_bin_devs = c("Region", "Year", "Bin", "Sex", "Fish_Fleet"),
  fish_sel_norm_bins = c("Bin", "Fish_Fleet"),
  use_fixed_fish_sel = "Fish_Fleet",
  fish_sel_bicubic_Wbin = c("Region", "Bin", "Node", "Block", "Fish_Fleet"),
  fish_sel_bicubic_Wyr = c("Region", "Year", "Node", "Block", "Fish_Fleet"),
  fishsel_sex_par_offset = "Fish_Fleet",
  fishsel_sex_scale_offset = "Fish_Fleet",
  fishsel_sex_apical_offset = "Fish_Fleet",
  fish_dbnrml_raw = c("Fish_Fleet", "Limb"),
  ret_dbnrml_raw = c("Fish_Fleet", "Limb"),
  fish_dbnrml_startbin = "Fish_Fleet",
  ret_dbnrml_startbin = "Fish_Fleet",
  fish_q_model = "Fish_Fleet",
  map_ln_fish_q_devs = c("Region", "Year", "Fish_Fleet"),
  map_ln_fishsel_devs = c("Region", "Year", "Bin", "Sex", "Fish_Fleet"),
  cont_tv_ret_sel = c("Region", "Fish_Fleet"),
  ret_sel_blocks = c("Region", "Year", "Fish_Fleet"),
  ret_sel_model = c("Region", "Year", "Fish_Fleet"),
  ret_sel_bicubic_binnodes = c("Region", "Year", "Fish_Fleet"),
  ret_sel_bicubic_yrnodes = c("Region", "Year", "Fish_Fleet"),
  ret_sel_bicubic_selstyr = c("Region", "Year", "Fish_Fleet"),
  ret_sel_bicubic_nselbins = c("Region", "Year", "Fish_Fleet"),
  use_fixed_ret_sel = "Fish_Fleet",
  ret_sel_input = c("Pop", "Region", "Year", "Seas", "Age", "Sex", "Fish_Fleet"),
  retsel_pe_wt = "Fish_Fleet",
  retsel_rw_init_sigma = "Fish_Fleet",
  ret_sel_bin_dev_bins = c("Bin", "Fish_Fleet"),
  cont_tv_retsel_bin_devs = "Fish_Fleet",
  map_ln_retsel_bin_devs = c("Region", "Year", "Bin", "Sex", "Fish_Fleet"),
  ret_sel_bicubic_Wbin = c("Region", "Bin", "Node", "Block", "Fish_Fleet"),
  ret_sel_bicubic_Wyr = c("Region", "Year", "Node", "Block", "Fish_Fleet"),
  retsel_sex_par_offset = "Fish_Fleet",
  retsel_sex_scale_offset = "Fish_Fleet",
  retsel_sex_apical_offset = "Fish_Fleet",
  map_ln_retsel_devs = c("Region", "Year", "Bin", "Sex", "Fish_Fleet"),
  cont_tv_srv_sel = c("Region", "Srv_Fleet"),
  srv_sel_blocks = c("Region", "Year", "Srv_Fleet"),
  srv_sel_model = c("Region", "Year", "Srv_Fleet"),
  srv_sel_bicubic_binnodes = c("Region", "Year", "Srv_Fleet"),
  srv_sel_bicubic_yrnodes = c("Region", "Year", "Srv_Fleet"),
  srv_sel_bicubic_selstyr = c("Region", "Year", "Srv_Fleet"),
  srv_sel_bicubic_nselbins = c("Region", "Year", "Srv_Fleet"),
  srv_q_blocks = c("Region", "Year", "Srv_Fleet"),
  srvsel_pe_wt = "Srv_Fleet",
  srvsel_rw_init_sigma = "Srv_Fleet",
  srv_sel_bin_dev_bins = c("Bin", "Srv_Fleet"),
  cont_tv_srvsel_bin_devs = "Srv_Fleet",
  map_ln_srvsel_bin_devs = c("Region", "Year", "Bin", "Sex", "Srv_Fleet"),
  srv_sel_norm_bins = c("Bin", "Srv_Fleet"),
  use_fixed_srv_sel = "Srv_Fleet",
  srv_sel_bicubic_Wbin = c("Region", "Bin", "Node", "Block", "Srv_Fleet"),
  srv_sel_bicubic_Wyr = c("Region", "Year", "Node", "Block", "Srv_Fleet"),
  srvsel_sex_par_offset = "Srv_Fleet",
  srvsel_sex_scale_offset = "Srv_Fleet",
  srvsel_sex_apical_offset = "Srv_Fleet",
  srv_dbnrml_raw = c("Srv_Fleet", "Limb"),
  srv_dbnrml_startbin = "Srv_Fleet",
  srv_q_type = "Srv_Fleet",
  srv_q_model = "Srv_Fleet",
  map_ln_srv_q_devs = c("Region", "Year", "Srv_Fleet"),
  map_ln_srvsel_devs = c("Region", "Year", "Bin", "Sex", "Srv_Fleet"),
  fishsel_bin_devs_rw_init_sigma = "Fish_Fleet",
  ret_sel_norm_bins = c("Bin", "Fish_Fleet"),
  retsel_bin_devs_rw_init_sigma = "Fish_Fleet",
  srvsel_bin_devs_rw_init_sigma = "Srv_Fleet",
  fish_sel_input = c("Pop", "Region", "Year", "Seas", "Age", "Sex", "Fish_Fleet"),
  srv_sel_input = c("Pop", "Region", "Year", "Seas", "Age", "Sex", "Srv_Fleet"),

  # catch and discards
  ObsCatchAA = c("Region", "Year", "Seas", "Obs_Age", "Sex", "Fish_Fleet"),
  UseCatchAA = c("Region", "Year", "Seas", "Obs_Age", "Sex", "Fish_Fleet"),
  ObsCatchAA_SE = c("Region", "Year", "Seas", "Obs_Age", "Sex", "Fish_Fleet"),
  ObsCatchAA_pop = c("Pop", "Region", "Year", "Seas", "Obs_Age", "Sex", "Fish_Fleet"),
  UseCatchAA_pop = c("Pop", "Region", "Year", "Seas", "Obs_Age", "Sex", "Fish_Fleet"),
  ObsCatchAA_pop_SE = c("Pop", "Region", "Year", "Seas", "Obs_Age", "Sex", "Fish_Fleet"),
  ObsDiscardAA = c("Region", "Year", "Seas", "Obs_Age", "Sex", "Fish_Fleet"),
  UseDiscardAA = c("Region", "Year", "Seas", "Obs_Age", "Sex", "Fish_Fleet"),
  ObsDiscardAA_SE = c("Region", "Year", "Seas", "Obs_Age", "Sex", "Fish_Fleet"),
  ObsDiscardAA_pop = c("Pop", "Region", "Year", "Seas", "Obs_Age", "Sex", "Fish_Fleet"),
  UseDiscardAA_pop = c("Pop", "Region", "Year", "Seas", "Obs_Age", "Sex", "Fish_Fleet"),
  ObsDiscardAA_pop_SE = c("Pop", "Region", "Year", "Seas", "Obs_Age", "Sex", "Fish_Fleet"),
  use_catch_aa = "Fish_Fleet",
  use_discard_aa = "Fish_Fleet",
  Catch_seas_Type = "Fish_Fleet",
  Catch_pop_seas_Type = "Fish_Fleet",
  Discard_seas_Type = "Fish_Fleet",
  Discard_pop_seas_Type = "Fish_Fleet",
  CatchAA_seas_Type = "Fish_Fleet",
  CatchAA_pop_seas_Type = "Fish_Fleet",
  DiscardAA_seas_Type = "Fish_Fleet",
  DiscardAA_pop_seas_Type = "Fish_Fleet",
  CatchAA_Type = c("Year", "Fish_Fleet"),
  CatchAA_pop_Type = c("Year", "Fish_Fleet"),
  DiscardAA_Type = c("Year", "Fish_Fleet"),
  DiscardAA_pop_Type = c("Year", "Fish_Fleet"),
  CatchAA_LikeType = "Fish_Fleet",
  CatchAA_sigma_form = "Fish_Fleet",
  CatchAA_pop_LikeType = "Fish_Fleet",
  CatchAA_pop_sigma_form = "Fish_Fleet",
  DiscardAA_LikeType = "Fish_Fleet",
  DiscardAA_sigma_form = "Fish_Fleet",
  DiscardAA_pop_LikeType = "Fish_Fleet",
  DiscardAA_pop_sigma_form = "Fish_Fleet",
  AgeObsCorr_catch = "Fish_Fleet",
  AgeObsCorr_catch_pop = "Fish_Fleet",
  AgeObsCorr_discard = "Fish_Fleet",
  AgeObsCorr_discard_pop = "Fish_Fleet",
  ObsCatch = c("Region", "Year", "Seas", "Fish_Fleet"),
  UseCatch = c("Region", "Year", "Seas", "Fish_Fleet"),
  UseCatch_pop = c("Pop", "Region", "Year", "Seas", "Fish_Fleet"),
  catch_units = "Fish_Fleet",
  UseDiscard = c("Region", "Year", "Seas", "Fish_Fleet"),
  UseDiscard_pop = c("Pop", "Region", "Year", "Seas", "Fish_Fleet"),
  discard_units = "Fish_Fleet",
  map_ln_F_devs = c("Region", "Year", "Seas", "Fish_Fleet"),
  map_logit_dmr_devs = c("Region", "Year", "Seas", "Fish_Fleet"),
  ObsDiscard = c("Region", "Year", "Seas", "Fish_Fleet"),
  ObsCatch_pop = c("Pop", "Region", "Year", "Seas", "Fish_Fleet"),
  ObsDiscard_pop = c("Pop", "Region", "Year", "Seas", "Fish_Fleet"),

  # fishery index and compositions
  ISS_FishAgeComps = c("Region", "Year", "Seas", "Sex", "Fish_Fleet"),
  ISS_FishLenComps = c("Region", "Year", "Seas", "Sex", "Fish_Fleet"),
  ISS_FishAgeComps_pop = c("Pop", "Region", "Year", "Seas", "Sex", "Fish_Fleet"),
  ISS_FishLenComps_pop = c("Pop", "Region", "Year", "Seas", "Sex", "Fish_Fleet"),
  ObsFishIdx = c("Region", "Year", "Seas", "Fish_Fleet"),
  ObsFishIdx_SE = c("Region", "Year", "Seas", "Fish_Fleet"),
  UseFishIdx = c("Region", "Year", "Seas", "Fish_Fleet"),
  UseFishIdx_pop = c("Pop", "Region", "Year", "Seas", "Fish_Fleet"),
  fish_idx_type = "Fish_Fleet",
  FishIdx_seas_Type = "Fish_Fleet",
  FishIdx_pop_seas_Type = "Fish_Fleet",
  FishAgeComps_seas_Type = "Fish_Fleet",
  FishAgeComps_pop_seas_Type = "Fish_Fleet",
  FishLenComps_seas_Type = "Fish_Fleet",
  FishLenComps_pop_seas_Type = "Fish_Fleet",
  FishAgeComps_discard_seas_Type = "Fish_Fleet",
  FishAgeComps_discard_pop_seas_Type = "Fish_Fleet",
  FishLenComps_discard_seas_Type = "Fish_Fleet",
  FishLenComps_discard_pop_seas_Type = "Fish_Fleet",
  FishIdx_LikeType = "Fish_Fleet",
  fish_idx_ages = c("Age", "Fish_Fleet"),
  FishAgeComps_bins = c("Obs_Age", "Fish_Fleet"),
  FishLenComps_bins = c("Obs_Len", "Fish_Fleet"),
  Fish_caal_bins = c("Obs_Age", "Fish_Fleet"),
  FishAgeComps_pop_bins = c("Obs_Age", "Fish_Fleet"),
  FishLenComps_pop_bins = c("Obs_Len", "Fish_Fleet"),
  t_fish = c("Region", "Seas", "Fish_Fleet"),
  ObsFishAgeComps = c("Region", "Year", "Seas", "Obs_Age", "Sex", "Fish_Fleet"),
  UseFishAgeComps = c("Region", "Year", "Seas", "Fish_Fleet"),
  ObsFishLenComps = c("Region", "Year", "Seas", "Obs_Len", "Sex", "Fish_Fleet"),
  UseFishLenComps = c("Region", "Year", "Seas", "Fish_Fleet"),
  UseFishAgeComps_pop = c("Pop", "Region", "Year", "Seas", "Fish_Fleet"),
  UseFishLenComps_pop = c("Pop", "Region", "Year", "Seas", "Fish_Fleet"),
  FishAgeComps_LikeType = "Fish_Fleet",
  FishLenComps_LikeType = "Fish_Fleet",
  FishAgeComps_pop_LikeType = "Fish_Fleet",
  FishLenComps_pop_LikeType = "Fish_Fleet",
  FishAgeComps_Type = c("Year", "Fish_Fleet"),
  FishLenComps_Type = c("Year", "Fish_Fleet"),
  FishAgeComps_pop_Type = c("Year", "Fish_Fleet"),
  FishLenComps_pop_Type = c("Year", "Fish_Fleet"),
  ISS_FishAgeComps_discard = c("Region", "Year", "Seas", "Sex", "Fish_Fleet"),
  ISS_FishLenComps_discard = c("Region", "Year", "Seas", "Sex", "Fish_Fleet"),
  ISS_FishAgeComps_discard_pop = c("Pop", "Region", "Year", "Seas", "Sex", "Fish_Fleet"),
  ISS_FishLenComps_discard_pop = c("Pop", "Region", "Year", "Seas", "Sex", "Fish_Fleet"),
  FishAgeComps_discard_bins = c("Obs_Age", "Fish_Fleet"),
  FishLenComps_discard_bins = c("Obs_Len", "Fish_Fleet"),
  FishAgeComps_discard_pop_bins = c("Obs_Age", "Fish_Fleet"),
  FishLenComps_discard_pop_bins = c("Obs_Len", "Fish_Fleet"),
  ObsFishAgeComps_discard = c("Region", "Year", "Seas", "Obs_Age", "Sex", "Fish_Fleet"),
  UseFishAgeComps_discard = c("Region", "Year", "Seas", "Fish_Fleet"),
  ObsFishLenComps_discard = c("Region", "Year", "Seas", "Obs_Len", "Sex", "Fish_Fleet"),
  UseFishLenComps_discard = c("Region", "Year", "Seas", "Fish_Fleet"),
  UseFishAgeComps_discard_pop = c("Pop", "Region", "Year", "Seas", "Fish_Fleet"),
  UseFishLenComps_discard_pop = c("Pop", "Region", "Year", "Seas", "Fish_Fleet"),
  FishAgeComps_discard_LikeType = "Fish_Fleet",
  FishLenComps_discard_LikeType = "Fish_Fleet",
  FishAgeComps_discard_pop_LikeType = "Fish_Fleet",
  FishLenComps_discard_pop_LikeType = "Fish_Fleet",
  FishAgeComps_discard_Type = c("Year", "Fish_Fleet"),
  FishLenComps_discard_Type = c("Year", "Fish_Fleet"),
  FishAgeComps_discard_pop_Type = c("Year", "Fish_Fleet"),
  FishLenComps_discard_pop_Type = c("Year", "Fish_Fleet"),
  ObsFish_caal = c("Region", "Year", "Seas", "CAAL_Len", "Obs_Age", "Sex", "Fish_Fleet"),
  UseFish_caal = c("Region", "Year", "Seas", "CAAL_Len", "Fish_Fleet"),
  ISS_Fish_caal = c("Region", "Year", "Seas", "CAAL_Len", "Sex", "Fish_Fleet"),
  Fish_caal_Type = c("Year", "Fish_Fleet"),
  Fish_caal_LikeType = "Fish_Fleet",
  ObsFishIdx_pop_SE = c("Pop", "Region", "Year", "Seas", "Fish_Fleet"),
  ObsFishIdx_pop = c("Pop", "Region", "Year", "Seas", "Fish_Fleet"),
  ObsFishAgeComps_pop = c("Pop", "Region", "Year", "Seas", "Obs_Age", "Sex", "Fish_Fleet"),
  ObsFishLenComps_pop = c("Pop", "Region", "Year", "Seas", "Obs_Len", "Sex", "Fish_Fleet"),
  ObsFishAgeComps_discard_pop = c("Pop", "Region", "Year", "Seas", "Obs_Age", "Sex", "Fish_Fleet"),
  ObsFishLenComps_discard_pop = c("Pop", "Region", "Year", "Seas", "Obs_Len", "Sex", "Fish_Fleet"),

  # survey index and compositions
  ISS_SrvAgeComps = c("Region", "Year", "Seas", "Sex", "Srv_Fleet"),
  ISS_SrvLenComps = c("Region", "Year", "Seas", "Sex", "Srv_Fleet"),
  ISS_SrvAgeComps_pop = c("Pop", "Region", "Year", "Seas", "Sex", "Srv_Fleet"),
  ISS_SrvLenComps_pop = c("Pop", "Region", "Year", "Seas", "Sex", "Srv_Fleet"),
  SrvIdx_seas_Type = "Srv_Fleet",
  SrvIdx_pop_seas_Type = "Srv_Fleet",
  SrvIdxAA_seas_Type = "Srv_Fleet",
  SrvIdxAA_pop_seas_Type = "Srv_Fleet",
  SrvAgeComps_seas_Type = "Srv_Fleet",
  SrvAgeComps_pop_seas_Type = "Srv_Fleet",
  SrvLenComps_seas_Type = "Srv_Fleet",
  SrvLenComps_pop_seas_Type = "Srv_Fleet",
  ObsSrvIdxAA = c("Region", "Year", "Seas", "Obs_Age", "Sex", "Srv_Fleet"),
  UseSrvIdxAA = c("Region", "Year", "Seas", "Obs_Age", "Sex", "Srv_Fleet"),
  ObsSrvIdxAA_SE = c("Region", "Year", "Seas", "Obs_Age", "Sex", "Srv_Fleet"),
  ObsSrvIdxAA_pop = c("Pop", "Region", "Year", "Seas", "Obs_Age", "Sex", "Srv_Fleet"),
  UseSrvIdxAA_pop = c("Pop", "Region", "Year", "Seas", "Obs_Age", "Sex", "Srv_Fleet"),
  ObsSrvIdxAA_pop_SE = c("Pop", "Region", "Year", "Seas", "Obs_Age", "Sex", "Srv_Fleet"),
  use_srv_idx_aa = "Srv_Fleet",
  SrvIdxAA_Type = c("Year", "Srv_Fleet"),
  SrvIdxAA_pop_Type = c("Year", "Srv_Fleet"),
  SrvIdxAA_LikeType = "Srv_Fleet",
  SrvIdxAA_sigma_form = "Srv_Fleet",
  SrvIdxAA_pop_LikeType = "Srv_Fleet",
  SrvIdxAA_pop_sigma_form = "Srv_Fleet",
  AgeObsCorr_srv_idx = "Srv_Fleet",
  AgeObsCorr_srv_idx_pop = "Srv_Fleet",
  ObsSrvIdx = c("Region", "Year", "Seas", "Srv_Fleet"),
  ObsSrvIdx_SE = c("Region", "Year", "Seas", "Srv_Fleet"),
  UseSrvIdx = c("Region", "Year", "Seas", "Srv_Fleet"),
  UseSrvIdx_pop = c("Pop", "Region", "Year", "Seas", "Srv_Fleet"),
  ObsSrvAgeComps = c("Region", "Year", "Seas", "Obs_Age", "Sex", "Srv_Fleet"),
  UseSrvAgeComps = c("Region", "Year", "Seas", "Srv_Fleet"),
  ObsSrvLenComps = c("Region", "Year", "Seas", "Obs_Len", "Sex", "Srv_Fleet"),
  UseSrvLenComps = c("Region", "Year", "Seas", "Srv_Fleet"),
  UseSrvAgeComps_pop = c("Pop", "Region", "Year", "Seas", "Srv_Fleet"),
  UseSrvLenComps_pop = c("Pop", "Region", "Year", "Seas", "Srv_Fleet"),
  SrvAgeComps_LikeType = "Srv_Fleet",
  SrvLenComps_LikeType = "Srv_Fleet",
  SrvAgeComps_pop_LikeType = "Srv_Fleet",
  SrvLenComps_pop_LikeType = "Srv_Fleet",
  SrvAgeComps_Type = c("Year", "Srv_Fleet"),
  SrvLenComps_Type = c("Year", "Srv_Fleet"),
  srv_idx_type = "Srv_Fleet",
  SrvAgeComps_pop_Type = c("Year", "Srv_Fleet"),
  SrvLenComps_pop_Type = c("Year", "Srv_Fleet"),
  SrvIdx_LikeType = "Srv_Fleet",
  srv_idx_ages = c("Age", "Srv_Fleet"),
  SrvAgeComps_bins = c("Obs_Age", "Srv_Fleet"),
  SrvLenComps_bins = c("Obs_Len", "Srv_Fleet"),
  Srv_caal_bins = c("Obs_Age", "Srv_Fleet"),
  SrvAgeComps_pop_bins = c("Obs_Age", "Srv_Fleet"),
  SrvLenComps_pop_bins = c("Obs_Len", "Srv_Fleet"),
  ObsSrv_caal = c("Region", "Year", "Seas", "CAAL_Len", "Obs_Age", "Sex", "Srv_Fleet"),
  UseSrv_caal = c("Region", "Year", "Seas", "CAAL_Len", "Srv_Fleet"),
  ISS_Srv_caal = c("Region", "Year", "Seas", "CAAL_Len", "Sex", "Srv_Fleet"),
  Srv_caal_Type = c("Year", "Srv_Fleet"),
  Srv_caal_LikeType = "Srv_Fleet",
  t_srv = c("Region", "Seas", "Srv_Fleet"),
  ObsSrvIdx_pop_SE = c("Pop", "Region", "Year", "Seas", "Srv_Fleet"),
  ObsSrvIdx_pop = c("Pop", "Region", "Year", "Seas", "Srv_Fleet"),
  ObsSrvAgeComps_pop = c("Pop", "Region", "Year", "Seas", "Obs_Age", "Sex", "Srv_Fleet"),
  ObsSrvLenComps_pop = c("Pop", "Region", "Year", "Seas", "Obs_Len", "Sex", "Srv_Fleet")

)
