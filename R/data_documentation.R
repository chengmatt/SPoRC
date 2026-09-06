# Package data Roxygen documentation for the example datasets shipped in data/. No code.

#' Sablefish data for single region case study
#'
#' A dataset containing the necessary elements for the Alaska sablefish case study.
#'
#' @format A list with multiple components needed for the single region sablefish model
#' @source 2024 Federal Alaska Sablefish Assessment
"sgl_rg_sable_data"

#' Sablefish data for multi region (5 area) case study
#'
#' A dataset containing the necessary elements for the Alaska sablefish spatial case study.
#'
#' @format A list with multiple components needed for the multi (5) region sablefish model
#' @source Cheng et al. 2025, Fish and Fisheries; Panmictic Panacea? Demonstrating Good Practices for Developing Spatial Stock Assessments through Application to Alaska Sablefish (Anoplopoma fimbria)
"mlt_rg_sable_data"

#' EBS Walleye Pollock data for single region case study
#'
#' A dataset containing the necessary elements for the EBS Walleye Pollock case study.
#'
#' @format A list with multiple components needed for the single region walleye pollock model
#' @source 2024 Federal EBS Walleye Pollock Assessment
"sgl_rg_ebswp_data"

#' Sablefish data for multi region (3 area) case study
#'
#' A dataset containing the necessary elements for the Alaska sablefish spatial case study.
#'
#' @format A list with multiple components needed for the multi (3) region sablefish model
#' @source Cheng et al. 2025, Fish and Fisheries; Panmictic Panacea? Demonstrating Good Practices for Developing Spatial Stock Assessments through Application to Alaska Sablefish (Anoplopoma fimbria)
"three_rg_sable_data"


#' Sablefish report for single region case study
#'
#' A report containing results for the Alaska sablefish case study.
#'
#' @format Report file from the single region sablefish case study
#' @source 2024 Federal Alaska Sablefish Assessment
"sgl_rg_sable_rep"

#' Sablefish report for 5 region case study
#'
#' A report containing results for the spatial Alaska sablefish case study.
#'
#' @format Report file from the 5 region sablefish case study
#' @source Cheng et al. 2025, Fish and Fisheries; Panmictic Panacea? Demonstrating Good Practices for Developing Spatial Stock Assessments through Application to Alaska Sablefish (Anoplopoma fimbria)
"mlt_rg_sable_rep"

#' Sablefish report for 3 region case study
#'
#' A report containing results for the spatial Alaska sablefish case study.
#'
#' @format Report file from the 3 region sablefish case study
#' @source Cheng et al. 2025, Fish and Fisheries; Panmictic Panacea? Demonstrating Good Practices for Developing Spatial Stock Assessments through Application to Alaska Sablefish (Anoplopoma fimbria)
"three_rg_sable_rep"

#' BSAI Blackspotted and Rougheye Rockfish data for single region case study
#'
#' A dataset containing the necessary elements for the BSAI blackspotted and rougheye rockfish case study,
#' including model inputs (dimensions, biologicals, catch, indices, and compositions) and the ADMB MLE
#' parameter estimates and derived quantities used for bridge verification.
#'
#' @format A list with components \code{inputs} (model dimensions and data), \code{mle} (ADMB maximum likelihood estimates),
#' and \code{admb} (ADMB derived quantities and likelihood components for bridge comparison)
"sgl_rg_rebs_data"

#' GOA Northern Rockfish data for single region case study
#'
#' A dataset containing the necessary elements for the GOA northern rockfish case study,
#' including model inputs (dimensions, biologicals, catch, indices, and compositions), the ADMB MLE
#' parameter estimates, and the ADMB derived quantities and likelihood components used for bridge verification.
#'
#' @format A list with model inputs at the top level, \code{mle} (ADMB maximum likelihood estimates),
#' and \code{admb} (ADMB derived quantities and likelihood components for bridge comparison)
#' @source 2024 Federal GOA Northern Rockfish Assessment
"sgl_rg_goa_nork_data"

#' BSAI Northern Rockfish data for single region case study
#'
#' A dataset containing the necessary elements for the BSAI northern rockfish case study,
#' including model inputs (dimensions, biologicals, catch, indices, and compositions), the ADMB MLE
#' parameter estimates, and the ADMB derived quantities and likelihood components used for bridge verification.
#'
#' @format A list with model inputs at the top level, \code{mle} (ADMB maximum likelihood estimates),
#' and \code{admb} (ADMB derived quantities and likelihood components for bridge comparison)
#' @source 2023 Federal BSAI Northern Rockfish Assessment
"sgl_rg_bsai_nork_data"

#' BSAI Pacific Ocean Perch data for single region case study
#'
#' A dataset containing the necessary elements for the BSAI Pacific ocean perch case study,
#' including model inputs (dimensions, biologicals, catch, two survey indices, and compositions), the ADMB MLE
#' parameter estimates, and the ADMB derived quantities and likelihood components used for bridge verification.
#'
#' @format A list with model inputs at the top level, \code{mle} (ADMB maximum likelihood estimates),
#' and \code{admb} (ADMB derived quantities and likelihood components for bridge comparison)
#' @source 2024 Federal BSAI Pacific Ocean Perch Assessment
"sgl_rg_bsai_pop_data"

#' BSAI Atka Mackerel data for single region case study
#'
#' A dataset containing the necessary elements for the BSAI Atka mackerel case study,
#' including model inputs (dimensions, biologicals, catch, survey index, and compositions), the ADMB MLE
#' parameter estimates, and the AMAK derived quantities and likelihood components used for bridge verification.
#'
#' @format A list with model inputs at the top level, \code{mle} (AMAK maximum likelihood estimates),
#' and \code{amak} (AMAK derived quantities and likelihood components for bridge comparison)
#' @source 2024 Federal BSAI Atka Mackerel Assessment, Model 16.0b
"sgl_rg_bsai_atka_data"

#' West Coast Sablefish data for single region case study
#'
#' A dataset containing the necessary elements for the 2025 West Coast sablefish case study,
#' including model inputs (dimensions, biologicals, catch, four survey indices, a recruitment
#' index, and sexed and unsexed compositions), the assessment's selectivity parameters and time
#' blocks, its maximum likelihood estimates, and its derived quantities and
#' likelihood components used for bridge verification.
#'
#' @format A list with model inputs at the top level, \code{sel_fish}/\code{sel_srv}/\code{sel_male}
#' (the assessment's double normal selectivity parameters, which of them it estimates, and which
#' blocks share one), \code{mle} (the assessment's maximum likelihood estimates, read from its
#' parameter file at twelve significant digits), and \code{ss3} (its derived
#' quantities and likelihood components for bridge comparison)
#' @source 2025 Pacific Fishery Management Council West Coast Sablefish Assessment
"sgl_rg_wc_sablefish_data"

#' GOA rex sole bridge data (Model 25.1)
#'
#' The 2025 Gulf of Alaska rex sole assessment (Model 25.1: two
#' areas with growth estimated separately in each, survey conditional
#' age-at-length, an ageing error matrix) parsed from its input files, with the
#' assessment's report quantities attached under \code{$ss3} and its maximum
#' likelihood estimates under \code{$mle}. Built by
#' \code{dev/rex_bridge/R/build_rex_data.R} from Carey McGilliard's goa_rex
#' repository. Used by \code{tests/testthat/helper-bridge_goa_rex.R} and the
#' rex sole regression test, which evaluates SPoRC at the assessment's own
#' estimate.
#'
#' @format A named list of observation arrays, biological inputs, the
#'   assessment's configuration, \code{$mle} (parameter values) and \code{$ss3} (report
#'   quantities: numbers at age, spawning biomass, recruitment, total biomass,
#'   catch, indices, growth tables, age-length keys, selectivity, one block of
#'   conditional age-at-length expected values, likelihood components).
#' @source \url{https://github.com/careymcgilliard/goa_rex}
"mlt_rg_goa_rex_data"


#' EBS Pacific cod bridge data (Model 24.1)
#'
#' The 2024 eastern Bering Sea Pacific cod assessment (Model 24.1: one area and
#' one sex, Richards growth with annual deviations on the
#' length at the young reference age and on K kept cohort by cohort,
#' length-based double normal selectivity with two fishery blocks and annual
#' deviations on the survey width, compositions recorded on coarser length bins
#' than the population has, two ageing error definitions) parsed from its
#' input files, with the assessment's report quantities attached under
#' \code{$ss3}, its year-by-year weight and fecundity at age under
#' \code{$wtatage}, and its maximum likelihood estimates under \code{$mle}.
#' Built by \code{dev/pcod_bridge/R/build_pcod_data.R} from the AFSC EBS_PCOD
#' repository. Used by \code{tests/testthat/helper-bridge_ebs_pcod.R} and the
#' Pacific cod regression test, which evaluates SPoRC at the assessment's own
#' estimate.
#'
#' @format A named list of observation arrays, biological inputs, the
#'   assessment's configuration, \code{$LenBinMap} (the population-to-data length bin map),
#'   \code{$MatAA} (maturity at age as the share of the weight at age that
#'   spawns), \code{$wtatage} (weight and fecundity at age by year),
#'   \code{$mle} (parameter values) and \code{$ss3} (report quantities:
#'   numbers at age, spawning biomass, recruitment, total biomass, catch, the
#'   index, growth parameters by year, selectivity at length and age, expected
#'   compositions, likelihood components and deviation penalties).
#' @source \url{https://github.com/afsc-assessments/EBS_PCOD}
"sgl_rg_ebs_pcod_data"

#' BSAI northern rock sole data for single region assessment case study
#'
#' A data list containing inputs for the 2024 BSAI Northern Rock Sole Assessment
#' (flatfish model, Model 24.2), together with the assessment's maximum
#' likelihood estimate and its reported output. The model is two sex with sex
#' specific natural mortality, free initial numbers at age estimated separately
#' for each sex, time varying logistic fishery selectivity with a male curve
#' offset, and a survey read as a July biomass index against January 1
#' compositions.
#'
#' @format Data list for single region northern rock sole assessment
#' @source McGilliard, C. R., Palsson, W., Haehn, R. 2024. Assessment of the northern rock sole stock in the Bering Sea and Aleutian Islands. North Pacific Fishery Management Council, Anchorage, AK.
"sgl_rg_bsai_nrs_data"

#' Dusky data for single region assessment case study
#'
#' A data list containing inputs for the 2024 GOA Dusky Rockfish Assessment
#'
#' @format Data list for single region dusky rockfish assessment
#' @source Omori, K. L., Williams, B. C., Hulson, P.-J., Ferriss, B. 2024. Assessment of the dusky rockfish stock in the Gulf of Alaska. North Pacific Fishery Management Council, Anchorage, AK.
"sgl_rg_dusky_data"

#' Dusky model outputs from single regino model
#'
#' A list containing inputs and outputs for the 2024 GOA Dusky Rockfish Assessment
#'
#' @format Data list for single region dusky rockfish assessment
#' @source Omori, K. L., Williams, B. C., Hulson, P.-J., Ferriss, B. 2024. Assessment of the dusky rockfish stock in the Gulf of Alaska. North Pacific Fishery Management Council, Anchorage, AK.
"dusky_rtmb_model"
