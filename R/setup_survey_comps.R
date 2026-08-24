# Stage 1 of 3: model setup
#
# Survey index and composition inputs. Setup_Mod_SrvIdx_and_Comps ingests the
# observed indices and the age and length compositions, and sets the
# composition likelihood choice along with its overdispersion and correlation
# parameters. Mirrors the fishery side but has no discard stream.

#' Set up observed survey indices and composition data
#'
#' Ingests observed survey index, age composition, and length composition data
#' (both pooled and population-specific) into \code{input_list$data},
#' initialises overdispersion and correlation starting values in
#' \code{input_list$par}, and constructs parameter maps via
#' \code{\link{do_comp_theta_mapping}} and \code{\link{do_comp_corr_pars_mapping}}
#' (called with \code{comp_prefix = "SrvAge"}/\code{"SrvLen"} and
#' \code{fleet_field = "n_srv_fleets"}). When \code{ISS_SrvAgeComps},
#' \code{ISS_SrvLenComps}, \code{ISS_SrvAgeComps_pop}, or
#' \code{ISS_SrvLenComps_pop} is \code{NULL}, input sample sizes are derived
#' automatically by summing observed composition counts across the appropriate
#' dimensions each year. Must be called after \code{\link{Setup_Mod_Dim}} and
#' before model compilation.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, \code{$map},
#'   and \code{$verbose} sublists, as returned by upstream setup functions.
#' @param ObsSrvIdx Observed survey index array
#'   \code{[n_regions × n_years × n_seas × n_srv_fleets]}.
#' @param ObsSrvIdx_SE Lognormal standard errors for \code{ObsSrvIdx}, same
#'   dimensions \code{[n_regions × n_years × n_seas × n_srv_fleets]}.
#' @param UseSrvIdx Binary indicator array
#'   \code{[n_regions × n_years × n_seas × n_srv_fleets]}. \code{1} = include
#'   in likelihood; \code{0} = exclude.
#' @param SrvLenComps_sel Character vector \code{[n_srv_fleets]}, whether a
#'   length-based selectivity is applied before or after the fish are spread
#'   over lengths. \code{"age"} (default) selects the index at age and spreads
#'   it afterwards; \code{"length"} spreads the numbers at each age over the
#'   key first and selects them length by length, so the survey sees the long
#'   fish of an age more often. The key is the survey's own, at \code{t_srv}.
#'   Requires length-based survey selectivity. Use \code{"length"} when
#'   selectivity is length based and the length compositions are what inform it.
#' @param srv_waa_selected Integer vector \code{[n_srv_fleets]} (0/1). With
#'   weight at age derived from growth and length-based selectivity, \code{1}
#'   makes a biomass index use the mean weight of the fish the survey sees at
#'   each age, \eqn{\sum_l P(l \mid a) s(l) w(l) / \sum_l P(l \mid a) s(l)},
#'   instead of the population mean weight at that age. The survey twin of
#'   \code{fish_waa_selected}. Only applies to an index in weight
#'   (\code{srv_idx_type = "biom"}).
#' @param srv_idx_ages Per-fleet selection of which ages contribute to the
#'   index total. Either a list with one element per survey fleet, where each
#'   element is a vector of ages or \code{NULL} for all ages, or an array
#'   \code{[n_ages x n_srv_fleets]} of 0/1 weights. Default \code{NULL} uses
#'   every age for every fleet. Restricting a fleet to a single age turns it
#'   into an index of that age alone, which is how an age-1 acoustic index is
#'   specified; the fleet's compositions are unaffected because the
#'   restriction applies to the index sum rather than to selectivity.
#' @param SrvIdx_LikeType Character vector \code{[n_srv_fleets]} giving the
#'   error structure of each survey index. Options are \code{"lognormal"}
#'   (default, the observation standard errors are on the log scale),
#'   \code{"normal"} (arithmetic scale), and \code{"mvn"} (multivariate
#'   normal on the arithmetic scale using a fixed covariance supplied through
#'   \code{SrvIdx_Cov}). One-step-ahead residuals are available only for
#'   lognormal fleets. A fleet's population-specific index stream follows the
#'   same choice for \code{"lognormal"} and \code{"normal"}, but stays
#'   lognormal under \code{"mvn"}, whose covariance describes the regional
#'   series only.
#' @param SrvIdx_Cov List with one element per survey fleet holding the fixed
#'   covariance matrix for fleets using \code{"mvn"}, and \code{NULL}
#'   otherwise. Each matrix must be square with one row per observation the
#'   fleet fits, ordered as the observations appear when scanning that fleet's
#'   \code{UseSrvIdx} slice in array order.
#' @param srv_idx_type Character vector \code{[n_srv_fleets]} specifying the
#'   index type per fleet. One of \code{"biom"} (biomass), \code{"abd"}
#'   (abundance), \code{"recdev"} (recruitment deviations), or \code{"none"}
#'   (no index for that fleet). Converted to integer codes (\code{1},
#'   \code{0}, \code{2}, \code{999}) before storage.
#'
#'
#'   A \code{"recdev"} fleet observes year class strength directly rather than
#'   any part of the population. Its predicted value is
#'   \code{q * (ln_RecDevs - mu)}, with \code{mu} the centre the recruitment
#'   penalty asserts for that year, so it measures the anomaly rather than the
#'   deviation as stored; under a bias ramp the two differ. Such a fleet reads
#'   no numbers at age, so its selectivity, survey timing and weight at age are
#'   unused and its compositions should be left off. It requires
#'   \code{SrvIdx_LikeType = "normal"}, since deviations are signed, and
#'   \code{RecDevs_pen_center = "fixed"} in \code{\link{Setup_Mod_Rec}}.
#' @param ObsSrvIdx_pop Observed population-specific survey index array
#'   \code{[n_pop × n_regions × n_years × n_seas × n_srv_fleets]}.
#' @param ObsSrvIdx_pop_SE Lognormal standard errors for \code{ObsSrvIdx_pop},
#'   same dimensions \code{[n_pop × n_regions × n_years × n_seas × n_srv_fleets]}.
#' @param UseSrvIdx_pop Binary indicator array
#'   \code{[n_pop × n_regions × n_years × n_seas × n_srv_fleets]}. \code{1} =
#'   include population-specific index in likelihood; \code{0} = exclude.
#'   Default: all zeros.
#' @param ObsSrvAgeComps Observed survey age compositions, array
#'   \code{[n_regions × n_years × n_seas × n_ages × n_sexes × n_srv_fleets]}.
#'   Values may be counts or proportions on a comparable scale.
#' @param UseSrvAgeComps Binary indicator array
#'   \code{[n_regions × n_years × n_seas × n_srv_fleets]}. \code{1} = fit age
#'   compositions; \code{0} = exclude.
#' @param ISS_SrvAgeComps Input sample sizes for survey age compositions, array
#'   \code{[n_regions × n_years × n_seas × n_sexes × n_srv_fleets]}, or
#'   \code{NULL} to derive automatically by summing \code{ObsSrvAgeComps}
#'   across the age dimension each year, respecting \code{SrvAgeComps_Type}.
#' @param ObsSrvLenComps Observed survey length compositions, array
#'   \code{[n_regions × n_years × n_seas × n_lens × n_sexes × n_srv_fleets]}.
#'   Only validated when \code{input_list$data$fit_lengths = 1} in \code{$data}.
#' @param UseSrvLenComps Binary indicator array
#'   \code{[n_regions × n_years × n_seas × n_srv_fleets]}. \code{1} = fit
#'   length compositions; \code{0} = exclude.
#' @param ISS_SrvLenComps Input sample sizes for survey length compositions,
#'   same structure as \code{ISS_SrvAgeComps}, or \code{NULL} for automatic
#'   derivation from \code{ObsSrvLenComps}.
#' @param SrvAgeComps_LikeType Character vector \code{[n_srv_fleets]}
#'   specifying the likelihood for survey age compositions. One of
#'   \code{"none"}, \code{"Multinomial"}, \code{"Dirichlet-Multinomial"},
#'   \code{"iid-Logistic-Normal"}, \code{"1d-Logistic-Normal"},
#'   \code{"2d-Logistic-Normal"}. Converted to integer codes
#'   (\code{999}, \code{0}-\code{4}) before storage.
#' @param SrvLenComps_LikeType Character vector \code{[n_srv_fleets]}
#'   specifying the likelihood for survey length compositions. Same options
#'   as \code{SrvAgeComps_LikeType}.
#' @param SrvAgeComps_Type Character vector defining the survey age composition
#'   structure per fleet and year range. Each element follows the format
#'   \code{"<type>_Year_<start>-<end>_Fleet_<fleet>"}. Use \code{"terminal"}
#'   in place of the end year to extend to the final model year. Valid types:
#'   \describe{
#'     \item{\code{"agg"}}{Aggregated across regions and sexes. Not compatible
#'       with \code{"2d-Logistic-Normal"}.}
#'     \item{\code{"spltRspltS"}}{Split by region and sex.}
#'     \item{\code{"spltRjntS"}}{Split by region, joint across sexes.}
#'     \item{\code{"none"}}{No composition data used.}
#'   }
#'   Parsed into a \code{[n_years × n_srv_fleets]} integer matrix before
#'   storage. An error is raised if any cell remains \code{NA} after parsing,
#'   indicating an incomplete year range specification.
#' @param SrvLenComps_Type Character vector defining the survey length
#'   composition structure. Same format and options as \code{SrvAgeComps_Type}.
#' @param ObsSrv_caal Observed conditional age-at-length array
#'   \code{[n_regions x n_years x n_seas x n_lens x n_ages x n_sexes x
#'   n_srv_fleets]}. A CAAL observation is the age composition of the fish aged
#'   from one length bin, so the age margin of each length row is what gets fit.
#'   \code{NULL} (default) for a model with no CAAL data.
#' @param UseSrv_caal Use flags \code{[n_regions x n_years x n_seas x n_lens x
#'   n_srv_fleets]}. Length bins with no aged fish carry a zero and are skipped.
#' @param ISS_Srv_caal Input sample sizes \code{[n_regions x n_years x n_seas x
#'   n_lens x n_sexes x n_srv_fleets]}. Summed from \code{ObsSrv_caal} when
#'   \code{NULL}.
#' @param Srv_caal_LikeType Character vector of length \code{n_srv_fleets}.
#'   One of \code{"none"}, \code{"Multinomial"} or
#'   \code{"Dirichlet-Multinomial"}. The logistic-normal families are not
#'   available for CAAL, since a single length bin's age sample is small and
#'   mostly zeros, which the additive log-ratio transform cannot handle.
#' @param Srv_caal_Type Composition type specification, using the same
#'   \code{"CompType_Year_x-y_Fleet_z"} convention as the marginal compositions.
#' @param SrvAgeComps_bins Which age bins each survey fleet's age composition is
#'   fitted over. Supply a list with one element per fleet, each a vector of age
#'   indices or \code{NULL} for all ages, or an \code{[n_ages x n_srv_fleets]}
#'   array of 0/1 weights. Both observed and expected compositions are
#'   restricted to the named bins and renormalized within them, so excluded bins
#'   are left out of the likelihood rather than being forced to be explained;
#'   this is how a fleet that only ages part of its age range is fitted. Indices
#'   refer to observed bins, that is after any ageing error has mapped model
#'   ages onto observed ones. Every fleet must retain at least one bin. Default
#'   \code{NULL}, which fits all ages for all fleets.
#' @param ObsSrvAgeComps_pop Observed population-specific survey age
#'   composition array
#'   \code{[n_pop × n_regions × n_years × n_seas × n_ages × n_sexes × n_srv_fleets]}.
#'   Required when any element of \code{UseSrvAgeComps_pop} is \code{1}.
#' @param UseSrvAgeComps_pop Binary indicator array
#'   \code{[n_pop × n_regions × n_years × n_seas × n_srv_fleets]}.
#'   \code{1} = fit population-specific age compositions; \code{0} = exclude.
#'   Default: all zeros.
#' @param ISS_SrvAgeComps_pop Input sample size array for population-specific
#'   survey age compositions
#'   \code{[n_pop × n_regions × n_years × n_seas × n_sexes × n_srv_fleets]}.
#'   If \code{NULL} (default), computed automatically by summing
#'   \code{ObsSrvAgeComps_pop} within each population-year-fleet-season-region
#'   cell according to \code{SrvAgeComps_pop_Type}.
#' @param ObsSrvLenComps_pop Observed population-specific survey length
#'   composition array
#'   \code{[n_pop × n_regions × n_years × n_seas × n_lens × n_sexes × n_srv_fleets]}.
#'   Required when \code{input_list$data$fit_lengths == 1} and any element of
#'   \code{UseSrvLenComps_pop} is \code{1}.
#' @param UseSrvLenComps_pop Binary indicator array
#'   \code{[n_pop × n_regions × n_years × n_seas × n_srv_fleets]}.
#'   \code{1} = fit population-specific length compositions; \code{0} = exclude.
#'   Default: all zeros.
#' @param ISS_SrvLenComps_pop Input sample size array for population-specific
#'   survey length compositions
#'   \code{[n_pop × n_regions × n_years × n_seas × n_sexes × n_srv_fleets]}.
#'   If \code{NULL} (default), derived automatically from
#'   \code{ObsSrvLenComps_pop}.
#' @param SrvAgeComps_pop_LikeType Character vector of length
#'   \code{n_srv_fleets} specifying the likelihood for population-specific
#'   survey age compositions. Same options as \code{SrvAgeComps_LikeType}.
#'   Default: \code{"none"} for all fleets.
#' @param SrvLenComps_pop_LikeType Character vector of length
#'   \code{n_srv_fleets} specifying the likelihood for population-specific
#'   survey length compositions. Same options as \code{SrvLenComps_LikeType}.
#'   Default: \code{"none"} for all fleets.
#' @param SrvAgeComps_pop_Type Character vector defining the composition
#'   structure for population-specific survey age compositions. Same format and
#'   options as \code{SrvAgeComps_Type}. Default: \code{"none"} for all fleets
#'   across all years.
#' @param SrvLenComps_pop_Type Character vector defining the composition
#'   structure for population-specific survey length compositions. Same format
#'   and options as \code{SrvLenComps_Type}. Default: \code{"none"} for all
#'   fleets across all years.
#' @param ... Optional named starting values for overdispersion and correlation
#'   parameters.
#'
#' @return The input \code{input_list} with survey data stored in \code{$data}
#'   (\code{ObsSrvIdx}, \code{ObsSrvIdx_SE}, \code{UseSrvIdx},
#'   \code{ObsSrvIdx_pop}, \code{ObsSrvIdx_pop_SE}, \code{UseSrvIdx_pop},
#'   \code{ObsSrvAgeComps}, \code{UseSrvAgeComps}, \code{ISS_SrvAgeComps},
#'   \code{ObsSrvLenComps}, \code{UseSrvLenComps}, \code{ISS_SrvLenComps},
#'   \code{ObsSrvAgeComps_pop}, \code{UseSrvAgeComps_pop},
#'   \code{ISS_SrvAgeComps_pop}, \code{ObsSrvLenComps_pop},
#'   \code{UseSrvLenComps_pop}, \code{ISS_SrvLenComps_pop},
#'   \code{SrvAgeComps_LikeType}, \code{SrvLenComps_LikeType},
#'   \code{SrvAgeComps_pop_LikeType}, \code{SrvLenComps_pop_LikeType},
#'   \code{SrvAgeComps_Type}, \code{SrvLenComps_Type},
#'   \code{SrvAgeComps_pop_Type}, \code{SrvLenComps_pop_Type},
#'   \code{srv_idx_type}); overdispersion and correlation starting values in
#'   \code{$par}; and factor maps in \code{$map} for all pooled and
#'   population-specific overdispersion and correlation parameter arrays.
#'
#' @export Setup_Mod_SrvIdx_and_Comps
#' @importFrom stringr str_detect
#' @family Model Setup
Setup_Mod_SrvIdx_and_Comps <- function(input_list,
                                       ObsSrvIdx,
                                       ObsSrvIdx_SE,
                                       UseSrvIdx,
                                       ObsSrvIdx_pop = NULL,
                                       ObsSrvIdx_pop_SE = NULL,
                                       UseSrvIdx_pop = array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_srv_fleets)),
                                       srv_idx_type,
                                       ObsSrvAgeComps,
                                       UseSrvAgeComps,
                                       ObsSrvLenComps,
                                       UseSrvLenComps,
                                       ISS_SrvAgeComps = NULL,
                                       ISS_SrvLenComps = NULL,
                                       SrvAgeComps_LikeType,
                                       SrvLenComps_LikeType,
                                       SrvAgeComps_Type,
                                       SrvLenComps_Type,
                                       ObsSrvAgeComps_pop = NULL,
                                       UseSrvAgeComps_pop = array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_srv_fleets)),
                                       ISS_SrvAgeComps_pop = NULL,
                                       ObsSrvLenComps_pop = NULL,
                                       UseSrvLenComps_pop = array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_srv_fleets)),
                                       ISS_SrvLenComps_pop = NULL,
                                       SrvAgeComps_pop_LikeType = rep("none", input_list$data$n_srv_fleets),
                                       SrvLenComps_pop_LikeType = rep("none", input_list$data$n_srv_fleets),
                                       SrvAgeComps_pop_Type = paste("none_Year_1-terminal_Fleet_", 1:input_list$data$n_srv_fleets, sep = ''),
                                       SrvLenComps_pop_Type = paste("none_Year_1-terminal_Fleet_", 1:input_list$data$n_srv_fleets, sep = ''),
                                       srv_idx_ages = NULL,
                                       SrvAgeComps_bins = NULL,
                                       SrvIdx_LikeType = rep("lognormal", input_list$data$n_srv_fleets),
                                       SrvLenComps_sel = rep("age", input_list$data$n_srv_fleets),
                                       srv_waa_selected = rep(0, input_list$data$n_srv_fleets),
                                       SrvIdx_Cov = NULL,

                                       # Conditional Age-at-Length
                                       ObsSrv_caal = NULL,
                                       UseSrv_caal = NULL,
                                       ISS_Srv_caal = NULL,
                                       Srv_caal_LikeType = rep("none", input_list$data$n_srv_fleets),
                                       Srv_caal_Type = paste("none_Year_1-terminal_Fleet_", 1:input_list$data$n_srv_fleets, sep = ''),
                                       ...
                                       ) {

  messages_list <<- character(0) # string to attach to for printing messages
  starting_values <- list(...)
  if(input_list$store_config) input_list$config$Setup_Mod_SrvIdx_and_Comps <- mget(names(formals()))[-1]

  # Input Validation --------------------------------------------------------

  # Survey Indices
  check_data_dimensions(ObsSrvIdx, n_regions = input_list$data$n_regions, n_seas = input_list$data$n_seas, n_years = length(input_list$data$years), n_srv_fleets = input_list$data$n_srv_fleets, what = 'ObsSrvIdx')
  check_data_dimensions(ObsSrvIdx_SE, n_regions = input_list$data$n_regions, n_seas = input_list$data$n_seas, n_years = length(input_list$data$years), n_srv_fleets = input_list$data$n_srv_fleets, what = 'ObsSrvIdx_SE')
  check_data_dimensions(UseSrvIdx, n_regions = input_list$data$n_regions, n_seas = input_list$data$n_seas, n_years = length(input_list$data$years), n_srv_fleets = input_list$data$n_srv_fleets, what = 'UseSrvIdx')
  if(any(UseSrvIdx_pop == 1)) {
    check_data_dimensions(ObsSrvIdx_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_srv_fleets = input_list$data$n_srv_fleets, what = 'ObsSrvIdx_pop')
    check_data_dimensions(ObsSrvIdx_pop_SE, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_srv_fleets = input_list$data$n_srv_fleets, what = 'ObsSrvIdx_pop_SE')
    check_data_dimensions(UseSrvIdx_pop, n_pop = input_list$data$n_pop,  n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_srv_fleets = input_list$data$n_srv_fleets, what = 'UseSrvIdx_pop')
  }
  if(!all(srv_idx_type %in% c("biom", "abd", "none", "recdev"))) stop("Invalid specification for srv_idx_type. Should be abd, biom, recdev, or none")

  # Survey compositions
  check_data_dimensions(ObsSrvAgeComps, n_regions = input_list$data$n_regions, n_seas = input_list$data$n_seas, n_years = length(input_list$data$years), n_sexes = input_list$data$n_sexes, n_srv_fleets = input_list$data$n_srv_fleets, what = 'ObsSrvAgeComps')
  check_data_dimensions(UseSrvAgeComps, n_regions = input_list$data$n_regions, n_seas = input_list$data$n_seas, n_years = length(input_list$data$years), n_srv_fleets = input_list$data$n_srv_fleets, what = 'UseSrvAgeComps')
  check_data_dimensions(UseSrvLenComps, n_regions = input_list$data$n_regions, n_seas = input_list$data$n_seas, n_years = length(input_list$data$years), n_srv_fleets = input_list$data$n_srv_fleets, what = 'UseSrvLenComps')
  if(input_list$data$fit_lengths == 1) check_data_dimensions(ObsSrvLenComps, n_seas = input_list$data$n_seas, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_lens = obs_len_bins(input_list), n_sexes = input_list$data$n_sexes, n_srv_fleets = input_list$data$n_srv_fleets, what = 'ObsSrvLenComps')
  if(!is.null(ISS_SrvAgeComps)) check_data_dimensions(ISS_SrvAgeComps, n_seas = input_list$data$n_seas, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_sexes = input_list$data$n_sexes, n_srv_fleets = input_list$data$n_srv_fleets, what = 'ISS_SrvAgeComps')
  if(!is.null(ISS_SrvLenComps)) check_data_dimensions(ISS_SrvLenComps, n_seas = input_list$data$n_seas, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_sexes = input_list$data$n_sexes, n_srv_fleets = input_list$data$n_srv_fleets, what = 'ISS_SrvLenComps')
  check_data_dimensions(SrvAgeComps_LikeType, n_srv_fleets = input_list$data$n_srv_fleets, what = 'SrvAgeComps_LikeType')
  check_data_dimensions(SrvLenComps_LikeType, n_srv_fleets = input_list$data$n_srv_fleets, what = 'SrvLenComps_LikeType')
  if(!all(SrvAgeComps_LikeType %in% c("none", "Multinomial", "Dirichlet-Multinomial", "iid-Logistic-Normal", "1d-Logistic-Normal", "2d-Logistic-Normal")))
    stop("Invalid specification for SrvAgeComps_LikeType Should be either none, Multinomial, Dirichlet-Multinomial, iid-Logistic-Normal, 1d-Logistic-Normal, 2d-Logistic-Normal")
  if(!all(SrvLenComps_LikeType %in% c("none", "Multinomial", "Dirichlet-Multinomial", "iid-Logistic-Normal", "1d-Logistic-Normal", "2d-Logistic-Normal")))
    stop("Invalid specification for SrvLenComps_LikeType Should be either none, Multinomial, Dirichlet-Multinomial, iid-Logistic-Normal, 1d-Logistic-Normal, 2d-Logistic-Normal")

  # Survey compositions (population-specific)
  if(any(UseSrvAgeComps_pop == 1)) check_data_dimensions(ObsSrvAgeComps_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_sexes = input_list$data$n_sexes, n_srv_fleets = input_list$data$n_srv_fleets, what = 'ObsSrvAgeComps_pop')
  check_data_dimensions(UseSrvAgeComps_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_srv_fleets = input_list$data$n_srv_fleets, what = 'UseSrvAgeComps_pop')
  check_data_dimensions(UseSrvLenComps_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_srv_fleets = input_list$data$n_srv_fleets, what = 'UseSrvLenComps_pop')
  if(input_list$data$fit_lengths == 1 && any(UseSrvLenComps_pop == 1)) check_data_dimensions(ObsSrvLenComps_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_seas = input_list$data$n_seas, n_years = length(input_list$data$years), n_lens = obs_len_bins(input_list), n_sexes = input_list$data$n_sexes, n_srv_fleets = input_list$data$n_srv_fleets, what = 'ObsSrvLenComps_pop')
  if(!is.null(ISS_SrvAgeComps_pop)) check_data_dimensions(ISS_SrvAgeComps_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_seas = input_list$data$n_seas, n_years = length(input_list$data$years), n_sexes = input_list$data$n_sexes, n_srv_fleets = input_list$data$n_srv_fleets, what = 'ISS_SrvAgeComps_pop')
  if(!is.null(ISS_SrvLenComps_pop)) check_data_dimensions(ISS_SrvLenComps_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_seas = input_list$data$n_seas, n_years = length(input_list$data$years), n_sexes = input_list$data$n_sexes, n_srv_fleets = input_list$data$n_srv_fleets, what = 'ISS_SrvLenComps_pop')
  check_data_dimensions(SrvAgeComps_pop_LikeType, n_srv_fleets = input_list$data$n_srv_fleets, what = 'SrvAgeComps_pop_LikeType')
  check_data_dimensions(SrvLenComps_pop_LikeType, n_srv_fleets = input_list$data$n_srv_fleets, what = 'SrvLenComps_pop_LikeType')
  if(!all(SrvAgeComps_pop_LikeType %in% c("none", "Multinomial", "Dirichlet-Multinomial", "iid-Logistic-Normal", "1d-Logistic-Normal", "2d-Logistic-Normal")))
    stop("Invalid specification for SrvAgeComps_pop_LikeType Should be either none, Multinomial, Dirichlet-Multinomial, iid-Logistic-Normal, 1d-Logistic-Normal, 2d-Logistic-Normal")
  if(!all(SrvLenComps_pop_LikeType %in% c("none", "Multinomial", "Dirichlet-Multinomial", "iid-Logistic-Normal", "1d-Logistic-Normal", "2d-Logistic-Normal")))
    stop("Invalid specification for SrvLenComps_pop_LikeType Should be either none, Multinomial, Dirichlet-Multinomial, iid-Logistic-Normal, 1d-Logistic-Normal, 2d-Logistic-Normal")

  # checking to make sure defaults are not applied
  if(any(UseSrvAgeComps_pop == 1)) {
    if(is.null(ObsSrvAgeComps_pop)) stop("ObsSrvAgeComps_pop is NULL, but UseSrvAgeComps_pop contains 1s!")
    if(any(str_detect(SrvAgeComps_pop_LikeType, "none"))) warning("SrvAgeComps_pop_LikeType has nones, but UseSrvAgeComps_pop contains 1s! Please verify!")
    if(any(str_detect(SrvAgeComps_pop_Type, "none"))) warning("SrvAgeComps_pop_Type has nones, but UseSrvAgeComps_pop contains 1s! Please verify!")
  }

  if(any(UseSrvLenComps_pop == 1)) {
    if(is.null(ObsSrvLenComps_pop)) stop("ObsSrvLenComps_pop is NULL, but UseSrvAgeComps_pop contains 1s!")
    if(any(str_detect(SrvLenComps_pop_LikeType, "none"))) warning("SrvLenComps_pop_LikeType has nones, but UseSrvAgeComps_pop contains 1s! Please verify!")
    if(any(str_detect(SrvLenComps_pop_Type, "none"))) warning("SrvLenComps_pop_Type has nones, but UseSrvAgeComps_pop contains 1s! Please verify!")
  }

  # Survey Index Options ----------------------------------------------------

  srv_idx_type_vals <- array(NA, dim = c( input_list$data$n_srv_fleets))
  for(f in 1:input_list$data$n_srv_fleets) {
    if(srv_idx_type[f] == 'biom') srv_idx_type_vals[f] <- 1 # biomass
    if(srv_idx_type[f] == 'abd') srv_idx_type_vals[f] <- 0 # abundance
    if(srv_idx_type[f] == 'recdev') {
      srv_idx_type_vals[f] <- 2 # recruitment deviations
      # a deviation is signed, so a lognormal cannot be used on it, and the anomaly is measured against the penalty's fixed centre
      if(SrvIdx_LikeType[f] != "normal") stop("srv_idx_type is 'recdev' for survey fleet ", f, ". Recruitment deviations are signed, so that fleet needs SrvIdx_LikeType = 'normal' rather than ", SrvIdx_LikeType[f], ".")
      if(!is.null(input_list$data$RecDevs_pen_center) && input_list$data$RecDevs_pen_center != 0) stop("srv_idx_type is 'recdev' for survey fleet ", f, ". It measures the deviation against the centre its penalty asserts, which is only defined when RecDevs_pen_center is 'fixed' in Setup_Mod_Rec.")
    }
    if(srv_idx_type[f] == 'none') srv_idx_type_vals[f] <- 999 # none
    collect_message(paste("Survey Index", "for survey fleet", f, "specified as:" , srv_idx_type[f]))
  } # end f loop


  # Survey Age Composition Options ------------------------------------------

  comp_srvage_like_vals <- vector()
  for(f in 1:input_list$data$n_srv_fleets) {
    if(SrvAgeComps_LikeType[f] == 'none') comp_srvage_like_vals <- c(comp_srvage_like_vals, 999)
    if(SrvAgeComps_LikeType[f] == "Multinomial") comp_srvage_like_vals <- c(comp_srvage_like_vals, 0)
    if(SrvAgeComps_LikeType[f] == "Dirichlet-Multinomial") comp_srvage_like_vals <- c(comp_srvage_like_vals, 1)
    if(SrvAgeComps_LikeType[f] == "iid-Logistic-Normal") comp_srvage_like_vals <- c(comp_srvage_like_vals, 2)
    if(SrvAgeComps_LikeType[f] == "1d-Logistic-Normal") comp_srvage_like_vals <- c(comp_srvage_like_vals, 3)
    if(SrvAgeComps_LikeType[f] == "2d-Logistic-Normal") comp_srvage_like_vals <- c(comp_srvage_like_vals, 4)
    collect_message(paste("Survey Age Composition Likelihoods", "for survey fleet", f, "specified as:" , SrvAgeComps_LikeType[f]))
  } # end f loop

  # Specifying composition type
  SrvAgeComps_Type_Mat <- array(NA, dim = c(length(input_list$data$years), input_list$data$n_srv_fleets))
  for(i in 1:length(SrvAgeComps_Type)) {

    # Extract out components from list
    tmp <- SrvAgeComps_Type[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    comps_type_tmp <- tmp_vec[1] # get composition type
    fleet <- as.numeric(tmp_vec[5]) # extract fleet index

    # Checking character string
    if(!comps_type_tmp %in% c("agg", "spltRspltS", "spltRjntS", 'none')) stop("SrvAgeComps_Type not specified correctly. Must be one of: agg, spltRspltS, spltRjntS, none")
    if(!fleet %in% c(1:input_list$data$n_srv_fleets)) stop("Invalid fleet specified for SrvAgeComps_Type This needs to be specified as CompType_Year_x-y_Fleet_x")

    # get year ranges
    if(!str_detect(tmp, "terminal")) { # if not terminal year
      year_range <- as.numeric(unlist(strsplit(tmp_vec[3], "-")))
      years <- year_range[1]:year_range[2] # get sequence of years
    } else { # if terminal year
      year_range <- unlist(strsplit(tmp_vec[3], '-'))[1] # get year range
      years <- as.numeric(year_range):length(input_list$data$years) # get sequence of years
    }

    # define composition types
    if(comps_type_tmp == "agg") {
      if(comp_srvage_like_vals[fleet] == 4) stop("Age composition likelihood specified as 2d-Logistic-Normal, but composition type is aggregated. This is not valid.")
      comps_type_val <- 0
    }
    if(comps_type_tmp == "spltRspltS") comps_type_val <- 1
    if(comps_type_tmp == "spltRjntS") comps_type_val <- 2
    if(comps_type_tmp == "none") comps_type_val <- 999

    # input into matrix
    SrvAgeComps_Type_Mat[years,fleet] <- comps_type_val
  } # end i

  if(any(is.na(SrvAgeComps_Type_Mat))) stop("SrvAgeComps_Type_Mat is returning an NA. Did you update the year range of SrvAgeComps_Type_Mat?")

  # Specifying composition likelihood for population-specific data
  comp_srvage_pop_like_vals <- vector()
  for(f in 1:input_list$data$n_srv_fleets) {
    if(SrvAgeComps_pop_LikeType[f] == 'none') comp_srvage_pop_like_vals <- c(comp_srvage_pop_like_vals, 999)
    if(SrvAgeComps_pop_LikeType[f] == "Multinomial") comp_srvage_pop_like_vals <- c(comp_srvage_pop_like_vals, 0)
    if(SrvAgeComps_pop_LikeType[f] == "Dirichlet-Multinomial") comp_srvage_pop_like_vals <- c(comp_srvage_pop_like_vals, 1)
    if(SrvAgeComps_pop_LikeType[f] == "iid-Logistic-Normal") comp_srvage_pop_like_vals <- c(comp_srvage_pop_like_vals, 2)
    if(SrvAgeComps_pop_LikeType[f] == "1d-Logistic-Normal") comp_srvage_pop_like_vals <- c(comp_srvage_pop_like_vals, 3)
    if(SrvAgeComps_pop_LikeType[f] == "2d-Logistic-Normal") comp_srvage_pop_like_vals <- c(comp_srvage_pop_like_vals, 4)
    collect_message(paste("Population Survey Age Composition Likelihoods", "for survey fleet", f, "specified as:" , SrvAgeComps_pop_LikeType[f]))
  } # end f loop

  # Specifying composition type
  SrvAgeComps_pop_Type_Mat <- array(NA, dim = c(length(input_list$data$years), input_list$data$n_srv_fleets))
  for(i in 1:length(SrvAgeComps_pop_Type)) {

    # Extract out components from list
    tmp <- SrvAgeComps_pop_Type[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    comps_type_tmp <- tmp_vec[1] # get composition type
    fleet <- as.numeric(tmp_vec[5]) # extract fleet index

    # Checking character string
    if(!comps_type_tmp %in% c("agg", "spltRspltS", "spltRjntS", 'none')) stop("SrvAgeComps_pop_Type not specified correctly. Must be one of: agg, spltRspltS, spltRjntS, none")
    if(!fleet %in% c(1:input_list$data$n_srv_fleets)) stop("Invalid fleet specified for SrvAgeComps_pop_Type. This needs to be specified as CompType_Year_x-y_Fleet_x")

    # get year ranges
    if(!str_detect(tmp, "terminal")) { # if not terminal year
      year_range <- as.numeric(unlist(strsplit(tmp_vec[3], "-")))
      years <- year_range[1]:year_range[2] # get sequence of years
    } else { # if terminal year
      year_range <- unlist(strsplit(tmp_vec[3], '-'))[1] # get year range
      years <- as.numeric(year_range):length(input_list$data$years) # get sequence of years
    }

    # Composition type
    # define composition types
    if(comps_type_tmp == "agg") {
      if(comp_srvage_pop_like_vals[fleet] == 4) stop("Population Age composition likelihood specified as 2d-Logistic-Normal, but composition type is aggregated. This is not valid.")
      comps_type_val <- 0
    }
    if(comps_type_tmp == "spltRspltS") comps_type_val <- 1
    if(comps_type_tmp == "spltRjntS") comps_type_val <- 2
    if(comps_type_tmp == "none") comps_type_val <- 999

    # input into matrix
    SrvAgeComps_pop_Type_Mat[years,fleet] <- comps_type_val
  } # end i

  if(any(is.na(SrvAgeComps_pop_Type_Mat))) stop("SrvAgeComps_pop_Type is returning an NA. Did you update the year range of SrvAgeComps_pop_Type?")

  # Survey Length Composition Options ---------------------------------------

  comp_srvlen_like_vals <- vector()
  for(f in 1:input_list$data$n_srv_fleets) {
    if(SrvLenComps_LikeType[f] == 'none') comp_srvlen_like_vals <- c(comp_srvlen_like_vals, 999)
    if(SrvLenComps_LikeType[f] == "Multinomial") comp_srvlen_like_vals <- c(comp_srvlen_like_vals, 0)
    if(SrvLenComps_LikeType[f] == "Dirichlet-Multinomial") comp_srvlen_like_vals <- c(comp_srvlen_like_vals, 1)
    if(SrvLenComps_LikeType[f] == "iid-Logistic-Normal") comp_srvlen_like_vals <- c(comp_srvlen_like_vals, 2)
    if(SrvLenComps_LikeType[f] == "1d-Logistic-Normal") comp_srvlen_like_vals <- c(comp_srvlen_like_vals, 3)
    if(SrvLenComps_LikeType[f] == "2d-Logistic-Normal") comp_srvlen_like_vals <- c(comp_srvlen_like_vals, 4)
    collect_message(paste("Survey Length Composition Likelihoods", "for survey fleet", f, "specified as:" , SrvLenComps_LikeType[f]))
  } # end f loop

  SrvLenComps_Type_Mat <- array(NA, dim = c(length(input_list$data$years), input_list$data$n_srv_fleets))
  for(i in 1:length(SrvLenComps_Type)) {

    # Extract out components from list
    tmp <- SrvLenComps_Type[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    comps_type_tmp <- tmp_vec[1] # get composition type
    fleet <- as.numeric(tmp_vec[5]) # extract fleet index

    # Checking character string
    if(!comps_type_tmp %in% c("agg", "spltRspltS", "spltRjntS", 'none')) stop("SrvLenComps_Type not specified correctly. Must be one of: agg, spltRspltS, spltRjntS, none")
    if(!fleet %in% c(1:input_list$data$n_srv_fleets)) stop("Invalid fleet specified for SrvLenComps_Type This needs to be specified as CompType_Year_x-y_Fleet_x")

    # get year ranges
    if(!str_detect(tmp, "terminal")) { # if not terminal year
      year_range <- as.numeric(unlist(strsplit(tmp_vec[3], "-")))
      years <- year_range[1]:year_range[2] # get sequence of years
    } else { # if terminal year
      year_range <- unlist(strsplit(tmp_vec[3], '-'))[1] # get year range
      years <- as.numeric(year_range):length(input_list$data$years) # get sequence of years
    }

    # define composition types
    if(comps_type_tmp == "agg") {
      if(comp_srvlen_like_vals[fleet] == 4) stop("Length composition likelihood specified as 2d-Logistic-Normal, but composition type is aggregated. This is not valid.")
      comps_type_val <- 0
    }
    if(comps_type_tmp == "spltRspltS") comps_type_val <- 1
    if(comps_type_tmp == "spltRjntS") comps_type_val <- 2
    if(comps_type_tmp == "none") comps_type_val <- 999

    # input into matrix
    SrvLenComps_Type_Mat[years,fleet] <- comps_type_val
  } # end i

  if(any(is.na(SrvLenComps_Type_Mat))) stop("SrvLenComps_Type_Mat is returning an NA. Did you update the year range of SrvLenComps_Type_Mat?")

  # Specifying composition likelihood for population-specific data
  comp_srvlen_pop_like_vals <- vector()
  for(f in 1:input_list$data$n_srv_fleets) {
    if(SrvLenComps_pop_LikeType[f] == 'none') comp_srvlen_pop_like_vals <- c(comp_srvlen_pop_like_vals, 999)
    if(SrvLenComps_pop_LikeType[f] == "Multinomial") comp_srvlen_pop_like_vals <- c(comp_srvlen_pop_like_vals, 0)
    if(SrvLenComps_pop_LikeType[f] == "Dirichlet-Multinomial") comp_srvlen_pop_like_vals <- c(comp_srvlen_pop_like_vals, 1)
    if(SrvLenComps_pop_LikeType[f] == "iid-Logistic-Normal") comp_srvlen_pop_like_vals <- c(comp_srvlen_pop_like_vals, 2)
    if(SrvLenComps_pop_LikeType[f] == "1d-Logistic-Normal") comp_srvlen_pop_like_vals <- c(comp_srvlen_pop_like_vals, 3)
    if(SrvLenComps_pop_LikeType[f] == "2d-Logistic-Normal") comp_srvlen_pop_like_vals <- c(comp_srvlen_pop_like_vals, 4)
    collect_message(paste("Population Survey Length Composition Likelihoods", "for survey fleet", f, "specified as:" , SrvLenComps_pop_LikeType[f]))
  } # end f loop

  # Specifying composition type
  SrvLenComps_pop_Type_Mat <- array(NA, dim = c(length(input_list$data$years), input_list$data$n_srv_fleets))
  for(i in 1:length(SrvLenComps_pop_Type)) {

    # Extract out components from list
    tmp <- SrvLenComps_pop_Type[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    comps_type_tmp <- tmp_vec[1] # get composition type
    fleet <- as.numeric(tmp_vec[5]) # extract fleet index

    # Checking character string
    if(!comps_type_tmp %in% c("agg", "spltRspltS", "spltRjntS", 'none')) stop("SrvLenComps_pop_Type not specified correctly. Must be one of: agg, spltRspltS, spltRjntS, none")
    if(!fleet %in% c(1:input_list$data$n_srv_fleets)) stop("Invalid fleet specified for SrvLenComps_pop_Type. This needs to be specified as CompType_Year_x-y_Fleet_x")

    # get year ranges
    if(!str_detect(tmp, "terminal")) { # if not terminal year
      year_range <- as.numeric(unlist(strsplit(tmp_vec[3], "-")))
      years <- year_range[1]:year_range[2] # get sequence of years
    } else { # if terminal year
      year_range <- unlist(strsplit(tmp_vec[3], '-'))[1] # get year range
      years <- as.numeric(year_range):length(input_list$data$years) # get sequence of years
    }

    # Composition type
    # define composition types
    if(comps_type_tmp == "agg") {
      if(comp_srvlen_pop_like_vals[fleet] == 4) stop("Population Len composition likelihood specified as 2d-Logistic-Normal, but composition type is aggregated. This is not valid.")
      comps_type_val <- 0
    }
    if(comps_type_tmp == "spltRspltS") comps_type_val <- 1
    if(comps_type_tmp == "spltRjntS") comps_type_val <- 2
    if(comps_type_tmp == "none") comps_type_val <- 999

    # input into matrix
    SrvLenComps_pop_Type_Mat[years,fleet] <- comps_type_val
  } # end i

  if(any(is.na(SrvLenComps_pop_Type_Mat))) stop("SrvLenComps_pop_Type is returning an NA. Did you update the year range of SrvLenComps_pop_Type?")

  # whether length selectivity is applied at length or through the size-age key.
  # Whether the selectivity is length based is only known once Setup_Mod_Srvsel_and_Q
  # runs, which checks this against it
  if(length(SrvLenComps_sel) != input_list$data$n_srv_fleets || !all(SrvLenComps_sel %in% c("age", "length"))) stop("SrvLenComps_sel must be one of age or length for each survey fleet")
  srv_len_comp_sel_vals <- as.numeric(SrvLenComps_sel == "length")
  for(sf in 1:input_list$data$n_srv_fleets) if(SrvLenComps_sel[sf] == "length") collect_message("Survey length compositions for fleet ", sf, " apply selectivity at length")

  # Survey Weight at Age Options ----------------------------------------------

  if(length(srv_waa_selected) != input_list$data$n_srv_fleets || !all(srv_waa_selected %in% c(0, 1))) stop("srv_waa_selected must be 0 or 1 for each survey fleet")
  if(any(srv_waa_selected == 1) && (is.null(input_list$data$derive_waa) || input_list$data$derive_waa != 1)) stop("srv_waa_selected = 1 needs waa_model = 'wt_len' in Setup_Mod_Biologicals")
  for(sf in which(srv_waa_selected == 1)) collect_message("Survey fleet ", sf, " takes its biomass index on the selection-weighted weight at age")

  # ISS Munging -------------------------------------------------------------

  # Survey Ages
  if(is.null(ISS_SrvAgeComps)) {
    collect_message("No ISS is specified for SrvAgeComps. ISS weighting is calculated by summing up values from ObsSrvAgeComps each year")
    ISS_SrvAgeComps <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_srv_fleets))
    for(y in 1:length(input_list$data$years)) {
      for(f in 1:input_list$data$n_srv_fleets) {
        for(seas in 1:input_list$data$n_seas) {
          # if aggregated across sexes and regions (0)
          if(SrvAgeComps_Type_Mat[y,f] == 0) ISS_SrvAgeComps[1,y,seas,1,f] <- sum(ObsSrvAgeComps[,y,seas,,,f])
          # if split by region and sex
          if(SrvAgeComps_Type_Mat[y,f] == 1) ISS_SrvAgeComps[,y,seas,,f] <- apply(ObsSrvAgeComps[,y,seas,,,f, drop = FALSE], c(1,4), sum)
          # if split by region, joint by sex
          if(SrvAgeComps_Type_Mat[y,f] == 2) ISS_SrvAgeComps[,y,seas,1,f] <- apply(ObsSrvAgeComps[,y,seas,,,f, drop = FALSE], 1, sum)
        } # end seas loop
      } # end f loop
    } # end y loop
  }

  # Survey Lengths
  if(is.null(ISS_SrvLenComps)) {
    collect_message("No ISS is specified for SrvLenComps. ISS weighting is calculated by summing up values from ObsSrvLenComps each year")
    ISS_SrvLenComps <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_srv_fleets))
    for(y in 1:length(input_list$data$years)) {
      for(f in 1:input_list$data$n_srv_fleets) {
        for(seas in 1:input_list$data$n_seas) {
          # if aggregated across sexes and regions (0)
          if(SrvLenComps_Type_Mat[y,f] == 0) ISS_SrvLenComps[1,y,seas,1,f] <- sum(ObsSrvLenComps[,y,seas,,,f])
          # if split by region and sex
          if(SrvLenComps_Type_Mat[y,f] == 1) ISS_SrvLenComps[,y,seas,,f] <- apply(ObsSrvLenComps[,y,seas,,,f, drop = FALSE], c(1,4), sum)
          # if split by region, joint by sex
          if(SrvLenComps_Type_Mat[y,f] == 2) ISS_SrvLenComps[,y,seas,1,f] <- apply(ObsSrvLenComps[,y,seas,,,f, drop = FALSE], 1, sum)
        } # end seas loop
      } # end f loop
    } # end y loop
  }

  # Srvery Ages
  if(is.null(ISS_SrvAgeComps_pop)) {
    collect_message("No ISS is specified for pop_SrvAgeComps. ISS weighting is calculated by summing up values from ObsSrvAgeComps_pop each year")
    ISS_SrvAgeComps_pop <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_srv_fleets))
    for(p in 1:input_list$data$n_pop) {
      for(y in 1:length(input_list$data$years)) {
        for(f in 1:input_list$data$n_srv_fleets) {
          for(seas in 1:input_list$data$n_seas) {
            # if aggregated across sexes and regions (0) or joint across sexes
            if(SrvAgeComps_pop_Type_Mat[y,f] == 0) ISS_SrvAgeComps_pop[p,1,y,seas,1,f] <- sum(ObsSrvAgeComps_pop[p,,y,seas,,,f])
            # if split by region and sex
            if(SrvAgeComps_pop_Type_Mat[y,f] == 1) ISS_SrvAgeComps_pop[p,,y,seas,,f] <- apply(ObsSrvAgeComps_pop[p,,y,seas,,,f, drop = FALSE], c(2,5), sum)
            # if split by region, joint by sex
            if(SrvAgeComps_pop_Type_Mat[y,f] == 2) ISS_SrvAgeComps_pop[p,,y,seas,1,f] <- apply(ObsSrvAgeComps_pop[p,,y,seas,,,f, drop = FALSE], 2, sum)
          } # end seas loop
        } # end f loop
      } # end y loop
    } # end p loop
  }

  # Srvery Lengths
  if(is.null(ISS_SrvLenComps_pop)) {
    collect_message("No ISS is specified for pop_SrvLenComps. ISS weighting is calculated by summing up values from ObsSrvLenComps_pop each year")
    ISS_SrvLenComps_pop <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_srv_fleets))
    for(p in 1:input_list$data$n_pop) {
      for(y in 1:length(input_list$data$years)) {
        for(f in 1:input_list$data$n_srv_fleets) {
          for(seas in 1:input_list$data$n_seas) {
            # if aggregated across sexes and regions (0)
            if(SrvLenComps_pop_Type_Mat[y,f] == 0) ISS_SrvLenComps_pop[p,1,y,seas,1,f] <- sum(ObsSrvLenComps_pop[p,,y,seas,,,f])
            # if split by region and sex
            if(SrvLenComps_pop_Type_Mat[y,f] == 1) ISS_SrvLenComps_pop[p,,y,seas,,f] <- apply(ObsSrvLenComps_pop[p,,y,seas,,,f, drop = FALSE], c(2,5), sum)
            # if split by region, joint by sex
            if(SrvLenComps_pop_Type_Mat[y,f] == 2) ISS_SrvLenComps_pop[p,,y,seas,1,f] <- apply(ObsSrvLenComps_pop[p,,y,seas,,,f, drop = FALSE], 2, sum)
          } # end seas loop
        } # end f loop
      } # end y loop
    } # end p loop
  }


  # Populate Data List ------------------------------------------------------

  input_list$data$ISS_SrvAgeComps <- ISS_SrvAgeComps
  input_list$data$ISS_SrvLenComps <- ISS_SrvLenComps
  input_list$data$ISS_SrvAgeComps_pop <- ISS_SrvAgeComps_pop
  input_list$data$ISS_SrvLenComps_pop <- ISS_SrvLenComps_pop
  input_list$data$ObsSrvIdx <- ObsSrvIdx
  input_list$data$ObsSrvIdx_SE <- ObsSrvIdx_SE
  input_list$data$UseSrvIdx <- UseSrvIdx
  input_list$data$ObsSrvIdx_pop <- ObsSrvIdx_pop
  input_list$data$ObsSrvIdx_pop_SE <- ObsSrvIdx_pop_SE
  input_list$data$UseSrvIdx_pop <- UseSrvIdx_pop
  input_list$data$ObsSrvAgeComps <- ObsSrvAgeComps
  input_list$data$UseSrvAgeComps <- UseSrvAgeComps
  input_list$data$ObsSrvLenComps <- ObsSrvLenComps
  input_list$data$UseSrvLenComps <- UseSrvLenComps
  input_list$data$ObsSrvAgeComps_pop <- ObsSrvAgeComps_pop
  input_list$data$UseSrvAgeComps_pop <- UseSrvAgeComps_pop
  input_list$data$ObsSrvLenComps_pop <- ObsSrvLenComps_pop
  input_list$data$UseSrvLenComps_pop <- UseSrvLenComps_pop
  input_list$data$SrvAgeComps_LikeType <- comp_srvage_like_vals
  input_list$data$SrvLenComps_LikeType <- comp_srvlen_like_vals
  input_list$data$SrvAgeComps_pop_LikeType <- comp_srvage_pop_like_vals
  input_list$data$SrvLenComps_pop_LikeType <- comp_srvlen_pop_like_vals
  input_list$data$SrvAgeComps_Type <- SrvAgeComps_Type_Mat
  input_list$data$SrvLenComps_Type <- SrvLenComps_Type_Mat
  input_list$data$srv_idx_type <- srv_idx_type_vals
  input_list$data$srv_len_comp_sel <- srv_len_comp_sel_vals
  input_list$data$srv_waa_selected <- srv_waa_selected
  input_list$data$SrvAgeComps_pop_Type <- SrvAgeComps_pop_Type_Mat
  input_list$data$SrvLenComps_pop_Type <- SrvLenComps_pop_Type_Mat

  # Index age selection and error structure ---------------------------------
  if(!all(SrvIdx_LikeType %in% c("lognormal", "normal", "mvn"))) stop("Invalid specification for SrvIdx_LikeType. Should be lognormal, normal, or mvn")
  if(length(SrvIdx_LikeType) != input_list$data$n_srv_fleets) stop("SrvIdx_LikeType is not length n_srv_fleets")

  srv_idx_like_vals <- convert_to_numeric(SrvIdx_LikeType, list(lognormal = 0, normal = 1, mvn = 2))
  srv_idx_ages_arr <- parse_idx_ages(srv_idx_ages, length(input_list$data$ages), input_list$data$n_srv_fleets, "srv_idx_ages")
  srv_idx_cov_parsed <- parse_idx_cov(SrvIdx_Cov, srv_idx_like_vals, UseSrvIdx, input_list$data$n_srv_fleets, "SrvIdx_Cov")

  for(f in 1:input_list$data$n_srv_fleets) {
    collect_message(paste("Survey Index likelihood for survey fleet", f, "specified as:", SrvIdx_LikeType[f]))
    if(sum(srv_idx_ages_arr[,f]) != length(input_list$data$ages)) {
      collect_message(paste("Survey Index for survey fleet", f, "is restricted to ages:", paste(which(srv_idx_ages_arr[,f] == 1), collapse = ", ")))
    }
  } # end f loop

  input_list$data$SrvIdx_LikeType <- srv_idx_like_vals
  input_list$data$srv_idx_ages <- srv_idx_ages_arr
  input_list$data$SrvAgeComps_bins <- parse_idx_ages(SrvAgeComps_bins, length(input_list$data$ages), input_list$data$n_srv_fleets, "SrvAgeComps_bins")
  input_list$data$SrvIdx_Cov <- srv_idx_cov_parsed

  # Populate Parameter List -------------------------------------------------

  # Dispersion parameters for the survey age comps
  if("ln_SrvAge_theta" %in% names(starting_values)) input_list$par$ln_SrvAge_theta <- starting_values$ln_SrvAge_theta
  else input_list$par$ln_SrvAge_theta <- array(0, dim = c(input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_srv_fleets))

  # logistic normal correlation parameters for survey age comps
  if("SrvAge_corr_pars" %in% names(starting_values)) input_list$par$SrvAge_corr_pars <- starting_values$SrvAge_corr_pars
  else input_list$par$SrvAge_corr_pars <- array(0.01, dim = c(input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_srv_fleets, 2))

  # aggregated
  if("ln_SrvAge_theta_agg" %in% names(starting_values)) input_list$par$ln_SrvAge_theta_agg <- starting_values$ln_SrvAge_theta_agg
  else input_list$par$ln_SrvAge_theta_agg <- array(0, dim = c(input_list$data$n_srv_fleets))

  # aggregated correlation parameters
  if("SrvAge_corr_pars_agg" %in% names(starting_values)) input_list$par$SrvAge_corr_pars_agg <- starting_values$SrvAge_corr_pars_agg
  else input_list$par$SrvAge_corr_pars_agg <- array(0.01, dim = c(input_list$data$n_srv_fleets))

  # Dispersion parameters for survey length comps
  if("ln_SrvLen_theta" %in% names(starting_values)) input_list$par$ln_SrvLen_theta <- starting_values$ln_SrvLen_theta
  else input_list$par$ln_SrvLen_theta <- array(0, dim = c(input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_srv_fleets))

  # logistic normal correlation parameters for survey length comps
  if("SrvLen_corr_pars" %in% names(starting_values)) input_list$par$SrvLen_corr_pars <- starting_values$SrvLen_corr_pars
  else input_list$par$SrvLen_corr_pars <- array(0.01, dim = c(input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_srv_fleets, 2))

  # aggregated
  if("ln_SrvLen_theta_agg" %in% names(starting_values)) input_list$par$ln_SrvLen_theta_agg <- starting_values$ln_SrvLen_theta_agg
  else input_list$par$ln_SrvLen_theta_agg <- array(0, dim = c(input_list$data$n_srv_fleets))

  if("SrvLen_corr_pars_agg" %in% names(starting_values)) input_list$par$SrvLen_corr_pars_agg <- starting_values$SrvLen_corr_pars_agg
  else input_list$par$SrvLen_corr_pars_agg <- array(0.01, dim = c(input_list$data$n_srv_fleets))

  # Dispersion parameters for the population survey age comps
  if("ln_SrvAge_pop_theta" %in% names(starting_values)) input_list$par$ln_SrvAge_pop_theta <- starting_values$ln_SrvAge_pop_theta
  else input_list$par$ln_SrvAge_pop_theta <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_srv_fleets))

  # logistic normal correlation parameters for population survey age comps
  if("SrvAge_pop_corr_pars" %in% names(starting_values)) input_list$par$SrvAge_pop_corr_pars <- starting_values$SrvAge_pop_corr_pars
  else input_list$par$SrvAge_pop_corr_pars <- array(0.01, dim = c(input_list$data$n_pop,input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_srv_fleets, 2))

  # aggregated population pars
  if("ln_SrvAge_pop_theta_agg" %in% names(starting_values)) input_list$par$ln_SrvAge_pop_theta_agg <- starting_values$ln_SrvAge_pop_theta_agg
  else input_list$par$ln_SrvAge_pop_theta_agg <- array(0, dim = c(input_list$data$n_pop,input_list$data$n_srv_fleets))

  # aggregated population correlation parameters
  if("SrvAge_pop_corr_pars_agg" %in% names(starting_values)) input_list$par$SrvAge_pop_corr_pars_agg <- starting_values$SrvAge_pop_corr_pars_agg
  else input_list$par$SrvAge_pop_corr_pars_agg <- array(0.01, dim = c(input_list$data$n_pop,input_list$data$n_srv_fleets))

  # Dispersion parameters for population survey length comps
  if("ln_SrvLen_pop_theta" %in% names(starting_values)) input_list$par$ln_SrvLen_pop_theta <- starting_values$ln_SrvLen_pop_theta
  else input_list$par$ln_SrvLen_pop_theta <- array(0, dim = c(input_list$data$n_pop,input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_srv_fleets))

  # logistic normal correlation parameters for population survey length comps
  if("SrvLen_pop_corr_pars" %in% names(starting_values)) input_list$par$SrvLen_pop_corr_pars <- starting_values$SrvLen_pop_corr_pars
  else input_list$par$SrvLen_pop_corr_pars <- array(0.01, dim = c(input_list$data$n_pop,input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_srv_fleets, 2))

  # aggregated population pars
  if("ln_SrvLen_pop_theta_agg" %in% names(starting_values)) input_list$par$ln_SrvLen_pop_theta_agg <- starting_values$ln_SrvLen_pop_theta_agg
  else input_list$par$ln_SrvLen_pop_theta_agg <- array(0, dim = c(input_list$data$n_pop,input_list$data$n_srv_fleets))

  if("SrvLen_pop_corr_pars_agg" %in% names(starting_values)) input_list$par$SrvLen_pop_corr_pars_agg <- starting_values$SrvLen_pop_corr_pars_agg
  else input_list$par$SrvLen_pop_corr_pars_agg <- array(0.01, dim = c(input_list$data$n_pop, input_list$data$n_srv_fleets))


  # Mapping Options ---------------------------------------------------------

  input_list <- do_comp_theta_mapping(input_list, comp_prefix = "SrvAge", fleet_field = "n_srv_fleets")
  input_list <- do_comp_theta_mapping(input_list, comp_prefix = "SrvLen", fleet_field = "n_srv_fleets")
  input_list <- do_comp_corr_pars_mapping(input_list, comp_prefix = "SrvAge", fleet_field = "n_srv_fleets")
  input_list <- do_comp_corr_pars_mapping(input_list, comp_prefix = "SrvLen", fleet_field = "n_srv_fleets")

  input_list <- do_comp_theta_mapping(input_list, comp_prefix = "SrvAge", has_pop = TRUE, fleet_field = "n_srv_fleets")
  input_list <- do_comp_theta_mapping(input_list, comp_prefix = "SrvLen", has_pop = TRUE, fleet_field = "n_srv_fleets")
  input_list <- do_comp_corr_pars_mapping(input_list, comp_prefix = "SrvAge", has_pop = TRUE, fleet_field = "n_srv_fleets")
  input_list <- do_comp_corr_pars_mapping(input_list, comp_prefix = "SrvLen", has_pop = TRUE, fleet_field = "n_srv_fleets")

  # Conditional Age-at-Length --------------------------------------------------
  input_list <- setup_caal_stream(
    input_list,
    ObsCAAL = ObsSrv_caal, UseCAAL = UseSrv_caal,
    ISS_CAAL = ISS_Srv_caal,
    CAAL_LikeType = Srv_caal_LikeType, CAAL_Type = Srv_caal_Type,
    fleet_type = "Srv"
  )

  # Print Messages ----------------------------------------------------------
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)
}
