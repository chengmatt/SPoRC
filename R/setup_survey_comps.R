# Stage 1 of 3: model setup
#
# Survey index and composition inputs. Setup_Mod_SrvIdx_and_Comps reads the indices and the age and length
# compositions and sets the composition likelihood. Mirrors the fishery side but has no discards.

#' Set up observed survey indices and composition data
#'
#' Sets the observed survey indices and compositions, pooled and
#' population-specific, with the overdispersion and correlation starting values and
#' the maps from \code{\link{do_comp_theta_mapping}} and
#' \code{\link{do_comp_corr_pars_mapping}}. A \code{NULL} \code{ISS_*} argument is
#' summed from the matching observed array each year, following that data source's
#' composition type. Call after \code{\link{Setup_Mod_Dim}}.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, \code{$map} and
#'   \code{$verbose}.
#' @param sigmaSrvIdx_spec,sigmaSrvIdx_pop_spec The estimated component of the
#'   aggregated and population-specific survey index observation error, one value
#'   per fleet. \code{"fix"} (default) uses the reported standard errors as they
#'   are. \code{"est_additive"} adds an estimated component to them,
#'   \code{"est_quadrature"} adds it in quadrature, and \code{"est_replace"}
#'   replaces them, as several ICES assessments do. An estimated component is
#'   confounded with a likelihood weight, since a weight on a normal likelihood is
#'   the same statement as dividing the variance, and \code{Setup_Mod_Weighting}
#'   warns when both are used. A fleet with a multivariate normal index takes its
#'   scale from the supplied covariance and cannot have one, which is an error.
#' @param sigmaSrvIdx_map,sigmaSrvIdx_pop_map Optional integer vectors
#'   \code{[n_srv_fleets]} of estimation groups for \code{ln_sigmaSrvIdx} and
#'   \code{ln_sigmaSrvIdx_pop}. Fleets sharing a value share a parameter and
#'   \code{NA} holds a fleet at its starting value. Defaults to one free parameter
#'   per fleet.
#' @param ObsSrvIdxAA Observed survey index at age \code{[n_regions, n_years,
#'   n_seas, n_obs_ages, n_sexes, n_srv_fleets]}, the ages being the columns of the
#'   fleet's ageing error matrix, through which the predicted index at each model
#'   age is read before it is compared. Supplying this fits the index at age
#'   directly, every age its own observation with its own catchability. The sex dim
#'   is required whatever the fleet reports: a data source summed over sexes has
#'   its observation in sex slot one. A fleet uses this or the aggregated index,
#'   never both.
#' @param UseSrvIdxAA Integer array shaped like \code{ObsSrvIdxAA}, \code{1} where
#'   an observation is fit.
#' @param ObsSrvIdxAA_pop,UseSrvIdxAA_pop Population-specific counterparts, with a
#'   leading population dim.
#' @param ObsSrvIdxAA_SE,ObsSrvIdxAA_pop_SE Reported standard errors shaped like
#'   their observation array, read only when \code{SrvIdxAA_sigma_form} asks for
#'   them.
#' @param sigmaSrvIdxAA_key,sigmaSrvIdxAA_pop_key Integer arrays \code{[n_obs_ages,
#'   n_sexes, n_srv_fleets]} coupling the index at age observation error, the key
#'   matrix ICES assessments use. Equal entries share a parameter and \code{NA}
#'   excludes one. The sex dim is required; a key coupling the sexes repeats its
#'   entries across them. The age shape of catchability is not set here: an index
#'   fit age by age puts it in selectivity through the \code{"nonparfree"} form.
#'   See \code{\link{Setup_Mod_Srvsel_and_Q}}.
#' @param sigmaSrvIdxAA_spec,sigmaSrvIdxAA_pop_spec \code{"est"} (default) or
#'   \code{"fix"}.
#' @param SrvIdxAA_Type,SrvIdxAA_pop_Type Which dims the fleet reports separately:
#'   \code{"agg"}, \code{"spltRaggS"} (default), \code{"aggRspltS"} or
#'   \code{"spltRspltS"}, or year and fleet specifications such as
#'   \code{"spltRaggS_Year_1-20_Fleet_1"}. See \code{\link{Setup_Mod_Catch_and_F}}.
#' @param SrvIdxAA_LikeType,SrvIdxAA_pop_LikeType \code{"lognormal"} (default) or
#'   \code{"normal"}, one setting for every fleet or one per fleet.
#' @param SrvIdxAA_sigma_form,SrvIdxAA_pop_sigma_form Where the observation error
#'   comes from: \code{"none"} (default), \code{"data"}, \code{"est_additive"} or
#'   \code{"est_quadrature"}.
#' @param AgeObsCorr_srv_idx,AgeObsCorr_srv_idx_pop Correlation across ages for the
#'   survey index at age: \code{"iid"} (default), \code{"1dar1"}, \code{"us"} or
#'   \code{"2dar1"}, one setting for every fleet or one per fleet. See
#'   \code{\link{Setup_Mod_Catch_and_F}}.
#' @param rho_srv_idx_spec,rho_srv_idx_pop_spec How the correlation parameters are
#'   shared over region, sex and fleet. \code{NULL} (default) gives one per fleet.
#'   See \code{\link{Setup_Mod_Catch_and_F}}.
#' @param ObsSrvIdx Observed survey index array \code{[n_regions × n_years × n_seas
#'   × n_srv_fleets]}.
#' @param ObsSrvIdx_SE Lognormal standard errors for \code{ObsSrvIdx}, same dims.
#' @param UseSrvIdx Binary array dimensioned like \code{ObsSrvIdx}, \code{1} to
#'   include the index in the likelihood.
#' @param SrvLenComps_sel Character vector \code{[n_srv_fleets]}, whether
#'   length-based selectivity applies before or after the fish are spread over
#'   lengths. \code{"age"} (default) selects the index at age and spreads it
#'   afterwards; \code{"length"} spreads first and selects length by length, so the
#'   survey sees the long fish of an age more often. The key is the survey's own at
#'   \code{t_srv}, and length-based survey selectivity is required.
#' @param srv_waa_selected Integer vector \code{[n_srv_fleets]} (0/1). With weight
#'   at age derived from growth and length-based selectivity, \code{1} makes a
#'   biomass index use the mean weight of the fish the survey sees at each age,
#'   \eqn{\sum_l P(l \mid a) s(l) w(l) / \sum_l P(l \mid a) s(l)}, rather than the
#'   population mean weight. The survey twin of \code{fish_waa_selected}, and only
#'   read for an index in weight.
#' @param srv_idx_ages Which ages contribute to each fleet's index total, either a
#'   list with one element per fleet (a vector of ages, or \code{NULL} for all) or
#'   an \code{[n_ages x n_srv_fleets]} array of 0/1 weights. \code{NULL} (default)
#'   uses every age. Restricting a fleet to one age makes it an index of that age
#'   alone, which is how an age-1 acoustic index is specified; the compositions are
#'   unaffected, since the restriction applies to the index sum.
#' @param SrvIdx_seas_Type,SrvIdx_pop_seas_Type,SrvIdxAA_seas_Type,SrvIdxAA_pop_seas_Type,SrvAgeComps_seas_Type,SrvAgeComps_pop_seas_Type,SrvLenComps_seas_Type,SrvLenComps_pop_seas_Type
#'   Whether a seasonal model reports this data source once a season or once a
#'   year, one value for every fleet or one per fleet. \code{"spltSeas"} (default)
#'   fits the observation against the prediction for the season it sits in;
#'   \code{"aggSeas"} sums the prediction over the year's seasons and fits one
#'   observation. Under \code{"aggSeas"} the observation stays in the season it was
#'   placed in, exactly one season per region and year may be on in the matching
#'   \code{Use} array, and the likelihood lands in that season. A survey measured at
#'   a point in time belongs in its own season with its own timing.
#' @param SrvIdx_LikeType Character vector \code{[n_srv_fleets]} of each index's
#'   error structure: \code{"lognormal"} (default, standard errors on the log
#'   scale), \code{"normal"} (arithmetic scale), or \code{"mvn"} (multivariate
#'   normal on the arithmetic scale with a fixed covariance from
#'   \code{SrvIdx_Cov}). One-step-ahead residuals are available for lognormal
#'   fleets only. A fleet's population-specific index follows the same choice for
#'   the first two and stays lognormal under \code{"mvn"}, whose covariance
#'   describes the regional series alone.
#' @param SrvIdx_Cov List with one element per fleet holding the fixed covariance
#'   for \code{"mvn"} fleets and \code{NULL} otherwise. Each matrix is square with
#'   one row per observation the fleet fits, ordered as they appear when scanning
#'   that fleet's \code{UseSrvIdx} slice in array order.
#' @param srv_idx_type Character vector \code{[n_srv_fleets]}: \code{"biom"},
#'   \code{"abd"}, \code{"recdev"} or \code{"none"}, stored as \code{1}, \code{0},
#'   \code{2} and \code{999}. A \code{"recdev"} fleet observes year class strength
#'   directly rather than any part of the population: its predicted value is
#'   \code{q * (ln_RecDevs - mu)}, with \code{mu} the center the recruitment
#'   penalty asserts for that year, so it measures the anomaly rather than the
#'   deviation as stored. It reads no numbers at age, so its selectivity, timing
#'   and weight at age are unused and its compositions should be left off. It
#'   requires \code{SrvIdx_LikeType = "normal"} and
#'   \code{RecDevs_pen_center = "fixed"} in \code{\link{Setup_Mod_Rec}}.
#' @param ObsSrvIdx_pop Observed population-specific index array \code{[n_pop ×
#'   n_regions × n_years × n_seas × n_srv_fleets]}.
#' @param ObsSrvIdx_pop_SE Lognormal standard errors for \code{ObsSrvIdx_pop}, same
#'   dims.
#' @param UseSrvIdx_pop Binary array dimensioned like \code{ObsSrvIdx_pop}. Default
#'   all zeros.
#' @param ObsSrvAgeComps Observed survey age compositions \code{[n_regions × n_years
#'   × n_seas × n_ages × n_sexes × n_srv_fleets]}, counts or proportions on a
#'   comparable scale.
#' @param UseSrvAgeComps Binary array \code{[n_regions × n_years × n_seas ×
#'   n_srv_fleets]}, \code{1} to fit the age compositions.
#' @param ISS_SrvAgeComps Input sample sizes \code{[n_regions × n_years × n_seas ×
#'   n_sexes × n_srv_fleets]}. \code{NULL} sums \code{ObsSrvAgeComps} over ages.
#' @param ObsSrvLenComps Observed survey length compositions \code{[n_regions ×
#'   n_years × n_seas × n_lens × n_sexes × n_srv_fleets]}. Only validated when
#'   \code{fit_lengths = 1}.
#' @param UseSrvLenComps Binary array \code{[n_regions × n_years × n_seas ×
#'   n_srv_fleets]}, \code{1} to fit the length compositions.
#' @param ISS_SrvLenComps Input sample sizes, structured as \code{ISS_SrvAgeComps}.
#'   \code{NULL} sums \code{ObsSrvLenComps}.
#' @param SrvAgeComps_LikeType Character vector \code{[n_srv_fleets]}:
#'   \code{"none"}, \code{"Multinomial"}, \code{"Dirichlet-Multinomial"}, the three
#'   logistic-normal forms or their three \code{-miss0} counterparts, stored as
#'   \code{999} and \code{0}-\code{7}.
#'
#'   The \code{miss0} forms drop the empty bins and renormalize the expected
#'   proportions over the bins that remain, rather than adding \code{addtocomp} to
#'   the zeros. Their standard deviation is divided by the square root of the input
#'   sample size, so the parameter is a per fish quantity and a year sampled harder
#'   is fit more tightly, and the change of variables from the log ratio is taken
#'   off so the result is a density on the composition itself. Their correlations
#'   run through the logistic function and are therefore positive, matching the
#'   autoregression they mirror. The \code{2d} form needs a composition joint
#'   across sexes. One-step-ahead residuals are not available for any of the three.
#' @param SrvLenComps_LikeType As \code{SrvAgeComps_LikeType}, for length
#'   compositions.
#' @param SrvAgeComps_Type Character vector of the composition structure per fleet
#'   and year range, each \code{"<type>_Year_<start>-<end>_Fleet_<fleet>"} with
#'   \code{"terminal"} allowed as the end year. Types are \code{"agg"} (aggregated
#'   across regions and sexes, not valid with \code{"2d-Logistic-Normal"}),
#'   \code{"spltRspltS"}, \code{"spltRjntS"} and \code{"none"}. Parsed into an
#'   \code{[n_years × n_srv_fleets]} integer matrix; a cell left \code{NA} means an
#'   incomplete year range and is an error.
#' @param SrvLenComps_Type As \code{SrvAgeComps_Type}, for length compositions.
#' @param ObsSrv_caal Observed conditional age-at-length array \code{[n_regions x
#'   n_years x n_seas x n_lens x n_ages x n_sexes x n_srv_fleets]}. An observation
#'   is the age composition of the fish aged from one length bin, so the age dim of
#'   each length row is what is fit. \code{NULL} (default) for no CAAL data.
#' @param UseSrv_caal Use flags \code{[n_regions x n_years x n_seas x n_lens x
#'   n_srv_fleets]}. Length bins with no aged fish take a zero and are skipped.
#' @param ISS_Srv_caal Input sample sizes \code{[n_regions x n_years x n_seas x
#'   n_lens x n_sexes x n_srv_fleets]}. Summed from \code{ObsSrv_caal} when
#'   \code{NULL}.
#' @param Srv_caal_LikeType Character vector \code{[n_srv_fleets]}: \code{"none"},
#'   \code{"Multinomial"} or \code{"Dirichlet-Multinomial"}. The logistic-normal
#'   families are not available for CAAL, since a length bin's age sample is small
#'   and mostly zeros.
#' @param Srv_caal_Type Composition type, in the same
#'   \code{"CompType_Year_x-y_Fleet_z"} vocabulary as the marginal compositions.
#' @param SrvAgeComps_bins Which age bins each fleet's age composition is fitted
#'   over, either a list with one element per fleet (bin indices, or \code{NULL}
#'   for all) or an \code{[n_obs_ages x n_srv_fleets]} array of 0/1 weights.
#'   Observed and expected are both restricted to the named bins and renormalized
#'   within them, so excluded bins leave the likelihood rather than being forced to
#'   be explained. Indices are observed bins, after any ageing error. For sex-joint
#'   comps the named bins are dropped from each sex's block, so the sex ratio
#'   becomes the ratio within the fitted bins. Every fleet must keep at least two
#'   bins. Default \code{NULL}, all bins.
#' @param SrvLenComps_bins Which length bins each fleet's length composition is
#'   fitted over, as \code{SrvAgeComps_bins}. Indices are observed length bins,
#'   after any \code{LenBinMap}.
#' @param Srv_caal_bins Which age bins each fleet's CAAL data are fitted over, as
#'   \code{SrvAgeComps_bins}, applied to every length bin's row of ages alike.
#' @param SrvAgeComps_pop_bins,SrvLenComps_pop_bins Which bins each fleet's
#'   population-specific age and length compositions are fitted over, as
#'   \code{SrvAgeComps_bins}.
#' @param ObsSrvAgeComps_pop Observed population-specific age composition array
#'   \code{[n_pop × n_regions × n_years × n_seas × n_ages × n_sexes ×
#'   n_srv_fleets]}. Required when any \code{UseSrvAgeComps_pop} is \code{1}.
#' @param UseSrvAgeComps_pop Binary array \code{[n_pop × n_regions × n_years ×
#'   n_seas × n_srv_fleets]}. Default all zeros.
#' @param ISS_SrvAgeComps_pop Input sample size array \code{[n_pop × n_regions ×
#'   n_years × n_seas × n_sexes × n_srv_fleets]}. \code{NULL} (default) sums
#'   \code{ObsSrvAgeComps_pop}.
#' @param ObsSrvLenComps_pop Observed population-specific length composition array
#'   \code{[n_pop × n_regions × n_years × n_seas × n_lens × n_sexes ×
#'   n_srv_fleets]}. Required when \code{fit_lengths == 1} and any
#'   \code{UseSrvLenComps_pop} is \code{1}.
#' @param UseSrvLenComps_pop Binary array \code{[n_pop × n_regions × n_years ×
#'   n_seas × n_srv_fleets]}. Default all zeros.
#' @param ISS_SrvLenComps_pop Input sample size array \code{[n_pop × n_regions ×
#'   n_years × n_seas × n_sexes × n_srv_fleets]}. \code{NULL} (default) sums
#'   \code{ObsSrvLenComps_pop}.
#' @param SrvAgeComps_pop_LikeType,SrvLenComps_pop_LikeType Character vectors
#'   \code{[n_srv_fleets]} for the population-specific compositions, with the same
#'   options as \code{SrvAgeComps_LikeType}. Default \code{"none"}.
#' @param SrvAgeComps_pop_Type,SrvLenComps_pop_Type Composition structure for the
#'   population-specific compositions, in the same format as
#'   \code{SrvAgeComps_Type}. Default \code{"none"} for every fleet and year.
#' @param ... Optional starting values for the overdispersion and correlation
#'   parameters.
#'
#' @return \code{input_list} with the survey observations, use flags, input sample
#'   sizes and integer-coded likelihood and composition types in \code{$data}, the
#'   overdispersion and correlation starting values in \code{$par}, and their
#'   factor maps, pooled and population-specific, in \code{$map}.
#'
#' @export Setup_Mod_SrvIdx_and_Comps
#' @importFrom stringr str_detect
#' @family Model Setup
Setup_Mod_SrvIdx_and_Comps <- function(input_list,
                                       ObsSrvIdx,
                                       ObsSrvIdx_SE,
                                       UseSrvIdx,
                                       ObsSrvIdxAA = NULL,
                                       UseSrvIdxAA = NULL,
                                       ObsSrvIdxAA_SE = NULL,
                                       ObsSrvIdxAA_pop = NULL,
                                       UseSrvIdxAA_pop = NULL,
                                       ObsSrvIdxAA_pop_SE = NULL,
                                       sigmaSrvIdxAA_key = NULL,
                                       sigmaSrvIdxAA_spec = "est",
                                       sigmaSrvIdxAA_pop_key = NULL,
                                       sigmaSrvIdxAA_pop_spec = "est",
                                       SrvIdxAA_Type = "spltRaggS",
                                       SrvIdxAA_pop_Type = "spltRaggS",
                                       SrvIdxAA_LikeType = "lognormal",
                                       SrvIdxAA_pop_LikeType = "lognormal",
                                       SrvIdxAA_sigma_form = "none",
                                       SrvIdxAA_pop_sigma_form = "none",
                                       AgeObsCorr_srv_idx = "iid",
                                       AgeObsCorr_srv_idx_pop = "iid",
                                       rho_srv_idx_spec = NULL,
                                       rho_srv_idx_pop_spec = NULL,
                                       sigmaSrvIdx_spec = "fix",
                                       sigmaSrvIdx_map = NULL,
                                       sigmaSrvIdx_pop_spec = "fix",
                                       sigmaSrvIdx_pop_map = NULL,
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
                                       SrvLenComps_bins = NULL,
                                       Srv_caal_bins = NULL,
                                       SrvAgeComps_pop_bins = NULL,
                                       SrvLenComps_pop_bins = NULL,
                                       SrvIdx_LikeType = rep("lognormal", input_list$data$n_srv_fleets),
                                       SrvIdx_seas_Type = NULL,
                                       SrvIdx_pop_seas_Type = NULL,
                                       SrvIdxAA_seas_Type = NULL,
                                       SrvIdxAA_pop_seas_Type = NULL,
                                       SrvAgeComps_seas_Type = NULL,
                                       SrvAgeComps_pop_seas_Type = NULL,
                                       SrvLenComps_seas_Type = NULL,
                                       SrvLenComps_pop_seas_Type = NULL,
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

  messages_list <<- character(0) # string to attach to for printing messages # nolint: object_usage_linter.
  starting_values <- list(...)
  if(input_list$store_config) input_list$config$Setup_Mod_SrvIdx_and_Comps <- mget(names(formals()))[-1]

  # Input Validation --------------------------------------------------------

  # Survey Indices
  check_data_dimensions(
    ObsSrvIdx,
    n_regions = input_list$data$n_regions,
    n_seas = input_list$data$n_seas,
    n_years = length(input_list$data$years),
    n_srv_fleets = input_list$data$n_srv_fleets,
    what = 'ObsSrvIdx'
  )
  check_data_dimensions(
    ObsSrvIdx_SE,
    n_regions = input_list$data$n_regions,
    n_seas = input_list$data$n_seas,
    n_years = length(input_list$data$years),
    n_srv_fleets = input_list$data$n_srv_fleets,
    what = 'ObsSrvIdx_SE'
  )
  check_data_dimensions(
    UseSrvIdx,
    n_regions = input_list$data$n_regions,
    n_seas = input_list$data$n_seas,
    n_years = length(input_list$data$years),
    n_srv_fleets = input_list$data$n_srv_fleets,
    what = 'UseSrvIdx'
  )
  if(any(UseSrvIdx_pop == 1)) {
    check_data_dimensions(
      ObsSrvIdx_pop,
      n_pop = input_list$data$n_pop,
      n_regions = input_list$data$n_regions,
      n_years = length(input_list$data$years),
      n_seas = input_list$data$n_seas,
      n_srv_fleets = input_list$data$n_srv_fleets,
      what = 'ObsSrvIdx_pop'
    )
    check_data_dimensions(
      ObsSrvIdx_pop_SE,
      n_pop = input_list$data$n_pop,
      n_regions = input_list$data$n_regions,
      n_years = length(input_list$data$years),
      n_seas = input_list$data$n_seas,
      n_srv_fleets = input_list$data$n_srv_fleets,
      what = 'ObsSrvIdx_pop_SE'
    )
    check_data_dimensions(
      UseSrvIdx_pop,
      n_pop = input_list$data$n_pop,
      n_regions = input_list$data$n_regions,
      n_years = length(input_list$data$years),
      n_seas = input_list$data$n_seas,
      n_srv_fleets = input_list$data$n_srv_fleets,
      what = 'UseSrvIdx_pop'
    )
  }
  if(!all(srv_idx_type %in% c("biom", "abd", "none", "recdev"))) stop("Invalid specification for srv_idx_type. Should be abd, biom, recdev, or none")

  # Survey compositions
  check_data_dimensions(
    ObsSrvAgeComps,
    n_regions = input_list$data$n_regions,
    n_seas = input_list$data$n_seas,
    n_years = length(input_list$data$years),
    n_sexes = input_list$data$n_sexes,
    n_srv_fleets = input_list$data$n_srv_fleets,
    what = 'ObsSrvAgeComps'
  )
  check_data_dimensions(
    UseSrvAgeComps,
    n_regions = input_list$data$n_regions,
    n_seas = input_list$data$n_seas,
    n_years = length(input_list$data$years),
    n_srv_fleets = input_list$data$n_srv_fleets,
    what = 'UseSrvAgeComps'
  )
  check_data_dimensions(
    UseSrvLenComps,
    n_regions = input_list$data$n_regions,
    n_seas = input_list$data$n_seas,
    n_years = length(input_list$data$years),
    n_srv_fleets = input_list$data$n_srv_fleets,
    what = 'UseSrvLenComps'
  )
  if(input_list$data$fit_lengths == 1) check_data_dimensions(
    ObsSrvLenComps,
    n_seas = input_list$data$n_seas,
    n_regions = input_list$data$n_regions,
    n_years = length(input_list$data$years),
    n_lens = obs_len_bins(input_list),
    n_sexes = input_list$data$n_sexes,
    n_srv_fleets = input_list$data$n_srv_fleets,
    what = 'ObsSrvLenComps'
  )
  if(!is.null(ISS_SrvAgeComps)) check_data_dimensions(
    ISS_SrvAgeComps,
    n_seas = input_list$data$n_seas,
    n_regions = input_list$data$n_regions,
    n_years = length(input_list$data$years),
    n_sexes = input_list$data$n_sexes,
    n_srv_fleets = input_list$data$n_srv_fleets,
    what = 'ISS_SrvAgeComps'
  )
  if(!is.null(ISS_SrvLenComps)) check_data_dimensions(
    ISS_SrvLenComps,
    n_seas = input_list$data$n_seas,
    n_regions = input_list$data$n_regions,
    n_years = length(input_list$data$years),
    n_sexes = input_list$data$n_sexes,
    n_srv_fleets = input_list$data$n_srv_fleets,
    what = 'ISS_SrvLenComps'
  )
  check_data_dimensions(SrvAgeComps_LikeType, n_srv_fleets = input_list$data$n_srv_fleets, what = 'SrvAgeComps_LikeType')
  check_data_dimensions(SrvLenComps_LikeType, n_srv_fleets = input_list$data$n_srv_fleets, what = 'SrvLenComps_LikeType')
  check_comp_like_type(SrvAgeComps_LikeType, "SrvAgeComps_LikeType")
  check_comp_like_type(SrvLenComps_LikeType, "SrvLenComps_LikeType")

  # Survey compositions (population-specific)
  if(any(UseSrvAgeComps_pop == 1)) check_data_dimensions(
    ObsSrvAgeComps_pop,
    n_pop = input_list$data$n_pop,
    n_regions = input_list$data$n_regions,
    n_years = length(input_list$data$years),
    n_seas = input_list$data$n_seas,
    n_sexes = input_list$data$n_sexes,
    n_srv_fleets = input_list$data$n_srv_fleets,
    what = 'ObsSrvAgeComps_pop'
  )
  check_data_dimensions(
    UseSrvAgeComps_pop,
    n_pop = input_list$data$n_pop,
    n_regions = input_list$data$n_regions,
    n_years = length(input_list$data$years),
    n_seas = input_list$data$n_seas,
    n_srv_fleets = input_list$data$n_srv_fleets,
    what = 'UseSrvAgeComps_pop'
  )
  check_data_dimensions(
    UseSrvLenComps_pop,
    n_pop = input_list$data$n_pop,
    n_regions = input_list$data$n_regions,
    n_years = length(input_list$data$years),
    n_seas = input_list$data$n_seas,
    n_srv_fleets = input_list$data$n_srv_fleets,
    what = 'UseSrvLenComps_pop'
  )
  if(input_list$data$fit_lengths == 1 && any(UseSrvLenComps_pop == 1)) check_data_dimensions(
    ObsSrvLenComps_pop,
    n_pop = input_list$data$n_pop,
    n_regions = input_list$data$n_regions,
    n_seas = input_list$data$n_seas,
    n_years = length(input_list$data$years),
    n_lens = obs_len_bins(input_list),
    n_sexes = input_list$data$n_sexes,
    n_srv_fleets = input_list$data$n_srv_fleets,
    what = 'ObsSrvLenComps_pop'
  )
  if(!is.null(ISS_SrvAgeComps_pop)) check_data_dimensions(
    ISS_SrvAgeComps_pop,
    n_pop = input_list$data$n_pop,
    n_regions = input_list$data$n_regions,
    n_seas = input_list$data$n_seas,
    n_years = length(input_list$data$years),
    n_sexes = input_list$data$n_sexes,
    n_srv_fleets = input_list$data$n_srv_fleets,
    what = 'ISS_SrvAgeComps_pop'
  )
  if(!is.null(ISS_SrvLenComps_pop)) check_data_dimensions(
    ISS_SrvLenComps_pop,
    n_pop = input_list$data$n_pop,
    n_regions = input_list$data$n_regions,
    n_seas = input_list$data$n_seas,
    n_years = length(input_list$data$years),
    n_sexes = input_list$data$n_sexes,
    n_srv_fleets = input_list$data$n_srv_fleets,
    what = 'ISS_SrvLenComps_pop'
  )
  check_data_dimensions(SrvAgeComps_pop_LikeType, n_srv_fleets = input_list$data$n_srv_fleets, what = 'SrvAgeComps_pop_LikeType')
  check_data_dimensions(SrvLenComps_pop_LikeType, n_srv_fleets = input_list$data$n_srv_fleets, what = 'SrvLenComps_pop_LikeType')
  check_comp_like_type(SrvAgeComps_pop_LikeType, "SrvAgeComps_pop_LikeType")
  check_comp_like_type(SrvLenComps_pop_LikeType, "SrvLenComps_pop_LikeType")

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

  srv_idx_type_vals <- array(NA, dim = c(input_list$data$n_srv_fleets))
  for(f in 1:input_list$data$n_srv_fleets) {
    if(srv_idx_type[f] == 'biom') srv_idx_type_vals[f] <- 1 # biomass
    if(srv_idx_type[f] == 'abd') srv_idx_type_vals[f] <- 0 # abundance
    if(srv_idx_type[f] == 'recdev') {
      srv_idx_type_vals[f] <- 2 # recruitment deviations
      # a deviation is signed, so a lognormal cannot be used on it, and the anomaly is measured against the penalty's fixed center
      if(SrvIdx_LikeType[f] != "normal") stop("srv_idx_type is 'recdev' for survey fleet ", f, ". Recruitment deviations are signed, so that fleet needs SrvIdx_LikeType = 'normal' rather than ", SrvIdx_LikeType[f], ".")
      if(!is.null(input_list$data$RecDevs_pen_center) && input_list$data$RecDevs_pen_center != 0) stop("srv_idx_type is 'recdev' for survey fleet ", f, ". It measures the deviation against the center its penalty asserts, which is only defined when RecDevs_pen_center is 'fixed' in Setup_Mod_Rec.")
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
    if(SrvAgeComps_LikeType[f] == "iid-Logistic-Normal-miss0") comp_srvage_like_vals <- c(comp_srvage_like_vals, 5)
    if(SrvAgeComps_LikeType[f] == "1d-Logistic-Normal-miss0") comp_srvage_like_vals <- c(comp_srvage_like_vals, 6)
    if(SrvAgeComps_LikeType[f] == "2d-Logistic-Normal-miss0") comp_srvage_like_vals <- c(comp_srvage_like_vals, 7)
    collect_message(paste("Survey Age Composition Likelihoods", "for survey fleet", f, "specified as:" , SrvAgeComps_LikeType[f]))
  } # end f loop

  check_miss0_osa(input_list, comp_srvage_like_vals, "SrvAgeComps_LikeType")

  # Specifying composition type
  SrvAgeComps_Type_Mat <- parse_year_fleet_spec(
    SrvAgeComps_Type, "SrvAgeComps_Type", input_list$data$n_srv_fleets, length(input_list$data$years),
    c(agg = 0, spltRspltS = 1, spltRjntS = 2, none = 999),
    check = function(value, fleet) {
      if(value %in% c("agg", "spltRspltS") && comp_srvage_like_vals[fleet] %in% c(4, 7))
        paste("The 2d logistic normal correlates bins with sexes, so it needs a",
              "composition joint across sexes. Use spltRjntS for this fleet.")
      else NULL
    })

  # Specifying composition likelihood for population-specific data
  comp_srvage_pop_like_vals <- vector()
  for(f in 1:input_list$data$n_srv_fleets) {
    if(SrvAgeComps_pop_LikeType[f] == 'none') comp_srvage_pop_like_vals <- c(comp_srvage_pop_like_vals, 999)
    if(SrvAgeComps_pop_LikeType[f] == "Multinomial") comp_srvage_pop_like_vals <- c(comp_srvage_pop_like_vals, 0)
    if(SrvAgeComps_pop_LikeType[f] == "Dirichlet-Multinomial") comp_srvage_pop_like_vals <- c(comp_srvage_pop_like_vals, 1)
    if(SrvAgeComps_pop_LikeType[f] == "iid-Logistic-Normal") comp_srvage_pop_like_vals <- c(comp_srvage_pop_like_vals, 2)
    if(SrvAgeComps_pop_LikeType[f] == "1d-Logistic-Normal") comp_srvage_pop_like_vals <- c(comp_srvage_pop_like_vals, 3)
    if(SrvAgeComps_pop_LikeType[f] == "2d-Logistic-Normal") comp_srvage_pop_like_vals <- c(comp_srvage_pop_like_vals, 4)
    if(SrvAgeComps_pop_LikeType[f] == "iid-Logistic-Normal-miss0") comp_srvage_pop_like_vals <- c(comp_srvage_pop_like_vals, 5)
    if(SrvAgeComps_pop_LikeType[f] == "1d-Logistic-Normal-miss0") comp_srvage_pop_like_vals <- c(comp_srvage_pop_like_vals, 6)
    if(SrvAgeComps_pop_LikeType[f] == "2d-Logistic-Normal-miss0") comp_srvage_pop_like_vals <- c(comp_srvage_pop_like_vals, 7)
    collect_message(paste("Population Survey Age Composition Likelihoods", "for survey fleet", f, "specified as:" , SrvAgeComps_pop_LikeType[f]))
  } # end f loop

  check_miss0_osa(input_list, comp_srvage_pop_like_vals, "SrvAgeComps_pop_LikeType")

  # Specifying composition type
  SrvAgeComps_pop_Type_Mat <- parse_year_fleet_spec(
    SrvAgeComps_pop_Type, "SrvAgeComps_pop_Type", input_list$data$n_srv_fleets, length(input_list$data$years),
    c(agg = 0, spltRspltS = 1, spltRjntS = 2, none = 999),
    check = function(value, fleet) {
      if(value %in% c("agg", "spltRspltS") && comp_srvage_pop_like_vals[fleet] %in% c(4, 7))
        paste("The 2d logistic normal correlates bins with sexes, so it needs a",
              "composition joint across sexes. Use spltRjntS for this fleet.")
      else NULL
    })

  # Survey Length Composition Options ---------------------------------------

  comp_srvlen_like_vals <- vector()
  for(f in 1:input_list$data$n_srv_fleets) {
    if(SrvLenComps_LikeType[f] == 'none') comp_srvlen_like_vals <- c(comp_srvlen_like_vals, 999)
    if(SrvLenComps_LikeType[f] == "Multinomial") comp_srvlen_like_vals <- c(comp_srvlen_like_vals, 0)
    if(SrvLenComps_LikeType[f] == "Dirichlet-Multinomial") comp_srvlen_like_vals <- c(comp_srvlen_like_vals, 1)
    if(SrvLenComps_LikeType[f] == "iid-Logistic-Normal") comp_srvlen_like_vals <- c(comp_srvlen_like_vals, 2)
    if(SrvLenComps_LikeType[f] == "1d-Logistic-Normal") comp_srvlen_like_vals <- c(comp_srvlen_like_vals, 3)
    if(SrvLenComps_LikeType[f] == "2d-Logistic-Normal") comp_srvlen_like_vals <- c(comp_srvlen_like_vals, 4)
    if(SrvLenComps_LikeType[f] == "iid-Logistic-Normal-miss0") comp_srvlen_like_vals <- c(comp_srvlen_like_vals, 5)
    if(SrvLenComps_LikeType[f] == "1d-Logistic-Normal-miss0") comp_srvlen_like_vals <- c(comp_srvlen_like_vals, 6)
    if(SrvLenComps_LikeType[f] == "2d-Logistic-Normal-miss0") comp_srvlen_like_vals <- c(comp_srvlen_like_vals, 7)
    collect_message(paste("Survey Length Composition Likelihoods", "for survey fleet", f, "specified as:" , SrvLenComps_LikeType[f]))
  } # end f loop

  check_miss0_osa(input_list, comp_srvlen_like_vals, "SrvLenComps_LikeType")

  SrvLenComps_Type_Mat <- parse_year_fleet_spec(
    SrvLenComps_Type, "SrvLenComps_Type", input_list$data$n_srv_fleets, length(input_list$data$years),
    c(agg = 0, spltRspltS = 1, spltRjntS = 2, none = 999),
    check = function(value, fleet) {
      if(value %in% c("agg", "spltRspltS") && comp_srvlen_like_vals[fleet] %in% c(4, 7))
        paste("The 2d logistic normal correlates bins with sexes, so it needs a",
              "composition joint across sexes. Use spltRjntS for this fleet.")
      else NULL
    })

  # Specifying composition likelihood for population-specific data
  comp_srvlen_pop_like_vals <- vector()
  for(f in 1:input_list$data$n_srv_fleets) {
    if(SrvLenComps_pop_LikeType[f] == 'none') comp_srvlen_pop_like_vals <- c(comp_srvlen_pop_like_vals, 999)
    if(SrvLenComps_pop_LikeType[f] == "Multinomial") comp_srvlen_pop_like_vals <- c(comp_srvlen_pop_like_vals, 0)
    if(SrvLenComps_pop_LikeType[f] == "Dirichlet-Multinomial") comp_srvlen_pop_like_vals <- c(comp_srvlen_pop_like_vals, 1)
    if(SrvLenComps_pop_LikeType[f] == "iid-Logistic-Normal") comp_srvlen_pop_like_vals <- c(comp_srvlen_pop_like_vals, 2)
    if(SrvLenComps_pop_LikeType[f] == "1d-Logistic-Normal") comp_srvlen_pop_like_vals <- c(comp_srvlen_pop_like_vals, 3)
    if(SrvLenComps_pop_LikeType[f] == "2d-Logistic-Normal") comp_srvlen_pop_like_vals <- c(comp_srvlen_pop_like_vals, 4)
    if(SrvLenComps_pop_LikeType[f] == "iid-Logistic-Normal-miss0") comp_srvlen_pop_like_vals <- c(comp_srvlen_pop_like_vals, 5)
    if(SrvLenComps_pop_LikeType[f] == "1d-Logistic-Normal-miss0") comp_srvlen_pop_like_vals <- c(comp_srvlen_pop_like_vals, 6)
    if(SrvLenComps_pop_LikeType[f] == "2d-Logistic-Normal-miss0") comp_srvlen_pop_like_vals <- c(comp_srvlen_pop_like_vals, 7)
    collect_message(paste("Population Survey Length Composition Likelihoods", "for survey fleet", f, "specified as:" , SrvLenComps_pop_LikeType[f]))
  } # end f loop

  check_miss0_osa(input_list, comp_srvlen_pop_like_vals, "SrvLenComps_pop_LikeType")

  # Specifying composition type
  SrvLenComps_pop_Type_Mat <- parse_year_fleet_spec(
    SrvLenComps_pop_Type, "SrvLenComps_pop_Type", input_list$data$n_srv_fleets, length(input_list$data$years),
    c(agg = 0, spltRspltS = 1, spltRjntS = 2, none = 999),
    check = function(value, fleet) {
      if(value %in% c("agg", "spltRspltS") && comp_srvlen_pop_like_vals[fleet] %in% c(4, 7))
        paste("The 2d logistic normal correlates bins with sexes, so it needs a",
              "composition joint across sexes. Use spltRjntS for this fleet.")
      else NULL
    })

  # whether length selectivity is applied at length or through the size-age key. only known to be
  # length based once the selectivity setup runs, which checks this against it
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
    for(y in seq_along(input_list$data$years)) {
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
    for(y in seq_along(input_list$data$years)) {
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
      for(y in seq_along(input_list$data$years)) {
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
      for(y in seq_along(input_list$data$years)) {
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

  # whether each data source is fit once a season or once a year as a season total
  n_srv <- input_list$data$n_srv_fleets
  input_list$data$SrvIdx_seas_Type <- parse_seas_agg_spec(SrvIdx_seas_Type, "SrvIdx_seas_Type", n_srv)
  input_list$data$SrvIdx_pop_seas_Type <- parse_seas_agg_spec(SrvIdx_pop_seas_Type, "SrvIdx_pop_seas_Type", n_srv)
  input_list$data$SrvIdxAA_seas_Type <- parse_seas_agg_spec(SrvIdxAA_seas_Type, "SrvIdxAA_seas_Type", n_srv)
  input_list$data$SrvIdxAA_pop_seas_Type <- parse_seas_agg_spec(SrvIdxAA_pop_seas_Type, "SrvIdxAA_pop_seas_Type", n_srv)
  input_list$data$SrvAgeComps_seas_Type <- parse_seas_agg_spec(SrvAgeComps_seas_Type, "SrvAgeComps_seas_Type", n_srv)
  input_list$data$SrvAgeComps_pop_seas_Type <- parse_seas_agg_spec(SrvAgeComps_pop_seas_Type, "SrvAgeComps_pop_seas_Type", n_srv)
  input_list$data$SrvLenComps_seas_Type <- parse_seas_agg_spec(SrvLenComps_seas_Type, "SrvLenComps_seas_Type", n_srv)
  input_list$data$SrvLenComps_pop_seas_Type <- parse_seas_agg_spec(SrvLenComps_pop_seas_Type, "SrvLenComps_pop_seas_Type", n_srv)

  # an aggregated observation is compared against the whole year, so only one season may be fit
  check_seas_agg_use(UseSrvIdx, input_list$data$SrvIdx_seas_Type, "UseSrvIdx")
  check_seas_agg_use(UseSrvIdx_pop, input_list$data$SrvIdx_pop_seas_Type, "UseSrvIdx_pop")
  check_seas_agg_use(UseSrvAgeComps, input_list$data$SrvAgeComps_seas_Type, "UseSrvAgeComps")
  check_seas_agg_use(UseSrvAgeComps_pop, input_list$data$SrvAgeComps_pop_seas_Type, "UseSrvAgeComps_pop")
  check_seas_agg_use(UseSrvLenComps, input_list$data$SrvLenComps_seas_Type, "UseSrvLenComps")
  check_seas_agg_use(UseSrvLenComps_pop, input_list$data$SrvLenComps_pop_seas_Type, "UseSrvLenComps_pop")

  for(sf in 1:n_srv) {
    if(input_list$data$SrvIdx_seas_Type[sf] == 1) collect_message("Survey index for survey fleet ", sf, " is fit as a season total")
    if(input_list$data$SrvAgeComps_seas_Type[sf] == 1) collect_message("Survey age compositions for survey fleet ", sf, " are fit as a season total")
    if(input_list$data$SrvLenComps_seas_Type[sf] == 1) collect_message("Survey length compositions for survey fleet ", sf, " are fit as a season total")
  } # end sf loop
  # Survey index at age. A fleet fits this or the aggregated index, never both.
  input_list <- do_at_age_data_setup(input_list, ObsSrvIdxAA, UseSrvIdxAA, ObsSrvIdxAA_SE,
                                     "SrvIdxAA", "n_srv_fleets")
  input_list <- do_at_age_data_setup(input_list, ObsSrvIdxAA_pop, UseSrvIdxAA_pop, ObsSrvIdxAA_pop_SE,
                                     "SrvIdxAA", "n_srv_fleets", pop = TRUE)

  use_srv_idx_aa <- rep(0, input_list$data$n_srv_fleets)
  for(sf in 1:input_list$data$n_srv_fleets) {
    if(any(input_list$data$UseSrvIdxAA[,,,,,sf] == 1) ||
       any(input_list$data$UseSrvIdxAA_pop[,,,,,,sf] == 1)) {
      use_srv_idx_aa[sf] <- 1
      if(any(UseSrvIdx[,,,sf] == 1)) {
        stop("Survey fleet ", sf, " has both an aggregated index and an index at age in use. ",
             "A fleet fits one or the other.")
      }
    }
  } # end sf loop
  input_list$data$use_srv_idx_aa <- use_srv_idx_aa

  check_seas_agg_use(input_list$data$UseSrvIdxAA, input_list$data$SrvIdxAA_seas_Type, "UseSrvIdxAA")
  check_seas_agg_use(input_list$data$UseSrvIdxAA_pop, input_list$data$SrvIdxAA_pop_seas_Type, "UseSrvIdxAA_pop")

  # at-age observation error for both data sources. Catchability at age is not set
  # here; it lives in selectivity, through the "nonparfree" form.
  input_list <- do_at_age_type_setup(input_list, SrvIdxAA_Type, "SrvIdxAA", "n_srv_fleets", "UseSrvIdxAA")
  input_list <- do_at_age_type_setup(input_list, SrvIdxAA_pop_Type, "SrvIdxAA", "n_srv_fleets", "UseSrvIdxAA_pop", pop = TRUE)
  input_list <- do_at_age_like_setup(input_list, SrvIdxAA_LikeType, SrvIdxAA_sigma_form, "SrvIdxAA", "n_srv_fleets")
  input_list <- do_at_age_like_setup(input_list, SrvIdxAA_pop_LikeType, SrvIdxAA_pop_sigma_form, "SrvIdxAA", "n_srv_fleets", pop = TRUE)

  input_list <- do_age_corr_setup(input_list, AgeObsCorr_srv_idx, "srv_idx", "n_srv_fleets",
                                  "UseSrvIdxAA", starting_values, rho_srv_idx_spec)
  input_list <- do_age_corr_setup(input_list, AgeObsCorr_srv_idx_pop, "srv_idx", "n_srv_fleets",
                                  "UseSrvIdxAA_pop", starting_values, rho_srv_idx_pop_spec, pop = TRUE)

  input_list <- do_key_mapping(input_list, sigmaSrvIdxAA_key,
                               at_age_sigma_spec(sigmaSrvIdxAA_spec, SrvIdxAA_sigma_form, any(use_srv_idx_aa == 1)),
                               "ln_sigmaSrvIdxAA", "n_srv_fleets", "UseSrvIdxAA", starting_values)
  input_list <- do_key_mapping(input_list, sigmaSrvIdxAA_pop_key,
                               at_age_sigma_spec(sigmaSrvIdxAA_pop_spec, SrvIdxAA_pop_sigma_form,
                                                 any(input_list$data$UseSrvIdxAA_pop == 1)),
                               "ln_sigmaSrvIdxAA_pop", "n_srv_fleets", "UseSrvIdxAA_pop",
                               starting_values, pop = TRUE)

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

  ## Index age selection and error structure --------------------------------
  if(!all(SrvIdx_LikeType %in% c("lognormal", "normal", "mvn"))) stop("Invalid specification for SrvIdx_LikeType. Should be lognormal, normal, or mvn")
  check_fleet_spec_length(SrvIdx_LikeType, input_list$data$n_srv_fleets, "SrvIdx_LikeType")

  srv_idx_like_vals <- convert_to_numeric(SrvIdx_LikeType, list(lognormal = 0, normal = 1, mvn = 2))
  srv_idx_ages_arr <- parse_bin_subset(srv_idx_ages, length(input_list$data$ages), input_list$data$n_srv_fleets, "srv_idx_ages")
  srv_idx_cov_parsed <- parse_idx_cov(SrvIdx_Cov, srv_idx_like_vals, UseSrvIdx, input_list$data$n_srv_fleets, "SrvIdx_Cov")

  for(f in 1:input_list$data$n_srv_fleets) {
    collect_message(paste("Survey Index likelihood for survey fleet", f, "specified as:", SrvIdx_LikeType[f]))
    if(sum(srv_idx_ages_arr[,f]) != length(input_list$data$ages)) {
      collect_message(paste("Survey Index for survey fleet", f, "is restricted to ages:", paste(which(srv_idx_ages_arr[,f] == 1), collapse = ", ")))
    }
  } # end f loop

  input_list$data$SrvIdx_LikeType <- srv_idx_like_vals
  input_list$data$srv_idx_ages <- srv_idx_ages_arr
  # composition bin restrictions. each data source is indexed on its own observed bins, so age data
  # sources are bounded by the observed age bins and length data sources by the length bins
  n_obs_age_bins <- obs_bin_count(input_list, ObsSrvAgeComps, 4, "age")
  n_obs_len_bins <- obs_bin_count(input_list, ObsSrvLenComps, 4, "len")
  # conditional age-at-length has its own observed age dimension, which need
  # not match the marginal age compositions, so it is measured off its own array
  n_obs_caal_bins <- obs_bin_count(input_list, ObsSrv_caal, 5, "age")
  # population-specific data sources likewise have their own arrays, with the
  # population dimension pushing the bins one place along
  n_obs_age_pop_bins <- obs_bin_count(input_list, ObsSrvAgeComps_pop, 5, "age")
  n_obs_len_pop_bins <- obs_bin_count(input_list, ObsSrvLenComps_pop, 5, "len")
  n_srv <- input_list$data$n_srv_fleets
  input_list$data$SrvAgeComps_bins <- check_comp_bins_min(parse_comp_bins(SrvAgeComps_bins, n_obs_age_bins, n_srv, "SrvAgeComps_bins"), comp_srvage_like_vals, "SrvAgeComps_bins")
  input_list$data$SrvLenComps_bins <- check_comp_bins_min(parse_comp_bins(SrvLenComps_bins, n_obs_len_bins, n_srv, "SrvLenComps_bins"), comp_srvlen_like_vals, "SrvLenComps_bins")
  input_list$data$Srv_caal_bins <- check_comp_bins_min(parse_comp_bins(Srv_caal_bins, n_obs_caal_bins, n_srv, "Srv_caal_bins"), ifelse(Srv_caal_LikeType == "none", 999, 0), "Srv_caal_bins")
  input_list$data$SrvAgeComps_pop_bins <- check_comp_bins_min(parse_comp_bins(SrvAgeComps_pop_bins, n_obs_age_pop_bins, n_srv, "SrvAgeComps_pop_bins"), comp_srvage_pop_like_vals, "SrvAgeComps_pop_bins")
  input_list$data$SrvLenComps_pop_bins <- check_comp_bins_min(parse_comp_bins(SrvLenComps_pop_bins, n_obs_len_pop_bins, n_srv, "SrvLenComps_pop_bins"), comp_srvlen_pop_like_vals, "SrvLenComps_pop_bins")

  # Reconcile the use flags with the restriction, so the fitting likelihood and
  # the residual routines agree on which blocks have data
  UseSrvAgeComps <- drop_empty_fitted_blocks(ObsSrvAgeComps, UseSrvAgeComps, input_list$data$SrvAgeComps_bins, 4, "SrvAgeComps")
  UseSrvLenComps <- drop_empty_fitted_blocks(ObsSrvLenComps, UseSrvLenComps, input_list$data$SrvLenComps_bins, 4, "SrvLenComps")
  UseSrvAgeComps_pop <- drop_empty_fitted_blocks(ObsSrvAgeComps_pop, UseSrvAgeComps_pop, input_list$data$SrvAgeComps_pop_bins, 5, "SrvAgeComps_pop")
  UseSrvLenComps_pop <- drop_empty_fitted_blocks(ObsSrvLenComps_pop, UseSrvLenComps_pop, input_list$data$SrvLenComps_pop_bins, 5, "SrvLenComps_pop")
  # written back over the copies stored before the restriction was known
  input_list$data$UseSrvAgeComps <- UseSrvAgeComps
  input_list$data$UseSrvLenComps <- UseSrvLenComps
  input_list$data$UseSrvAgeComps_pop <- UseSrvAgeComps_pop
  input_list$data$UseSrvLenComps_pop <- UseSrvLenComps_pop

  input_list$data$SrvIdx_Cov <- srv_idx_cov_parsed

  # Populate Parameter List -------------------------------------------------

  # Dispersion parameters for the survey age comps
  input_list$par$ln_SrvAge_theta <- array(0, dim = c(input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_srv_fleets))
  input_list$par$ln_SrvAge_theta <- use_starting_value(input_list$par$ln_SrvAge_theta, starting_values, "ln_SrvAge_theta")

  # logistic normal correlation parameters for survey age comps
  input_list$par$SrvAge_corr_pars <- array(0.01, dim = c(input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_srv_fleets, 2))
  input_list$par$SrvAge_corr_pars <- use_starting_value(input_list$par$SrvAge_corr_pars, starting_values, "SrvAge_corr_pars")

  # aggregated
  input_list$par$ln_SrvAge_theta_agg <- array(0, dim = c(input_list$data$n_srv_fleets))
  input_list$par$ln_SrvAge_theta_agg <- use_starting_value(input_list$par$ln_SrvAge_theta_agg, starting_values, "ln_SrvAge_theta_agg")

  # aggregated correlation parameters
  input_list$par$SrvAge_corr_pars_agg <- array(0.01, dim = c(input_list$data$n_srv_fleets))
  input_list$par$SrvAge_corr_pars_agg <- use_starting_value(input_list$par$SrvAge_corr_pars_agg, starting_values, "SrvAge_corr_pars_agg")

  # Dispersion parameters for survey length comps
  input_list$par$ln_SrvLen_theta <- array(0, dim = c(input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_srv_fleets))
  input_list$par$ln_SrvLen_theta <- use_starting_value(input_list$par$ln_SrvLen_theta, starting_values, "ln_SrvLen_theta")

  # logistic normal correlation parameters for survey length comps
  input_list$par$SrvLen_corr_pars <- array(0.01, dim = c(input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_srv_fleets, 2))
  input_list$par$SrvLen_corr_pars <- use_starting_value(input_list$par$SrvLen_corr_pars, starting_values, "SrvLen_corr_pars")

  # aggregated
  input_list$par$ln_SrvLen_theta_agg <- array(0, dim = c(input_list$data$n_srv_fleets))
  input_list$par$ln_SrvLen_theta_agg <- use_starting_value(input_list$par$ln_SrvLen_theta_agg, starting_values, "ln_SrvLen_theta_agg")

  input_list$par$SrvLen_corr_pars_agg <- array(0.01, dim = c(input_list$data$n_srv_fleets))
  input_list$par$SrvLen_corr_pars_agg <- use_starting_value(input_list$par$SrvLen_corr_pars_agg, starting_values, "SrvLen_corr_pars_agg")

  # Dispersion parameters for the population survey age comps
  input_list$par$ln_SrvAge_pop_theta <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_srv_fleets))
  input_list$par$ln_SrvAge_pop_theta <- use_starting_value(input_list$par$ln_SrvAge_pop_theta, starting_values, "ln_SrvAge_pop_theta")

  # logistic normal correlation parameters for population survey age comps
  input_list$par$SrvAge_pop_corr_pars <- array(0.01, dim = c(input_list$data$n_pop,input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_srv_fleets, 2))
  input_list$par$SrvAge_pop_corr_pars <- use_starting_value(input_list$par$SrvAge_pop_corr_pars, starting_values, "SrvAge_pop_corr_pars")

  # aggregated population pars
  input_list$par$ln_SrvAge_pop_theta_agg <- array(0, dim = c(input_list$data$n_pop,input_list$data$n_srv_fleets))
  input_list$par$ln_SrvAge_pop_theta_agg <- use_starting_value(input_list$par$ln_SrvAge_pop_theta_agg, starting_values, "ln_SrvAge_pop_theta_agg")

  # aggregated population correlation parameters
  input_list$par$SrvAge_pop_corr_pars_agg <- array(0.01, dim = c(input_list$data$n_pop,input_list$data$n_srv_fleets))
  input_list$par$SrvAge_pop_corr_pars_agg <- use_starting_value(input_list$par$SrvAge_pop_corr_pars_agg, starting_values, "SrvAge_pop_corr_pars_agg")

  # Dispersion parameters for population survey length comps
  input_list$par$ln_SrvLen_pop_theta <- array(0, dim = c(input_list$data$n_pop,input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_srv_fleets))
  input_list$par$ln_SrvLen_pop_theta <- use_starting_value(input_list$par$ln_SrvLen_pop_theta, starting_values, "ln_SrvLen_pop_theta")

  # logistic normal correlation parameters for population survey length comps
  input_list$par$SrvLen_pop_corr_pars <- array(0.01, dim = c(input_list$data$n_pop,input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_srv_fleets, 2))
  input_list$par$SrvLen_pop_corr_pars <- use_starting_value(input_list$par$SrvLen_pop_corr_pars, starting_values, "SrvLen_pop_corr_pars")

  # aggregated population pars
  input_list$par$ln_SrvLen_pop_theta_agg <- array(0, dim = c(input_list$data$n_pop,input_list$data$n_srv_fleets))
  input_list$par$ln_SrvLen_pop_theta_agg <- use_starting_value(input_list$par$ln_SrvLen_pop_theta_agg, starting_values, "ln_SrvLen_pop_theta_agg")

  input_list$par$SrvLen_pop_corr_pars_agg <- array(0.01, dim = c(input_list$data$n_pop, input_list$data$n_srv_fleets))
  input_list$par$SrvLen_pop_corr_pars_agg <- use_starting_value(input_list$par$SrvLen_pop_corr_pars_agg, starting_values, "SrvLen_pop_corr_pars_agg")


  # Survey index observation error stuff
  input_list$par$ln_sigmaSrvIdx <- rep(log(0.01), input_list$data$n_srv_fleets)
  input_list$par$ln_sigmaSrvIdx <- use_starting_value(input_list$par$ln_sigmaSrvIdx, starting_values, "ln_sigmaSrvIdx")

  input_list$par$ln_sigmaSrvIdx_pop <- rep(log(0.01), input_list$data$n_srv_fleets)
  input_list$par$ln_sigmaSrvIdx_pop <- use_starting_value(input_list$par$ln_sigmaSrvIdx_pop, starting_values, "ln_sigmaSrvIdx_pop")
  sigmaIdx_specs <- c("fix", "est_additive", "est_quadrature", "est_replace")
  for(spec_name in c("sigmaSrvIdx_spec", "sigmaSrvIdx_pop_spec")) {
    v <- get(spec_name)
    if(!v %in% sigmaIdx_specs) stop(spec_name, " is '", v, "', which is not recognized. Valid options: ",
                                    paste(sigmaIdx_specs, collapse = ", "))
  } # end spec_name loop

  sigmaIdx_forms <- list(fix = 0, est_additive = 1, est_quadrature = 2, est_replace = 3)
  input_list$data$sigmaSrvIdx_form <- convert_to_numeric(sigmaSrvIdx_spec, sigmaIdx_forms)
  input_list$data$sigmaSrvIdx_pop_form <- convert_to_numeric(sigmaSrvIdx_pop_spec, sigmaIdx_forms)

  # Check sigma is not estimated if using MVN index
  if(sigmaSrvIdx_spec != "fix" && any(srv_idx_like_vals == 2)) {
    stop("sigmaSrvIdx_spec is '", sigmaSrvIdx_spec, "' but survey fleet(s) ",
         paste(which(srv_idx_like_vals == 2), collapse = ", "),
         " use a multivariate normal index likelihood, which takes its scale from ",
         "SrvIdx_Cov and ignores the standard deviation. Fix the estimated sigma for ",
         "those fleets through sigmaSrvIdx_map, or use a lognormal or normal likelihood.")
  }

  # Mapping Options ---------------------------------------------------------

  input_list <- do_sigmaIdx_mapping(input_list, sigmaSrvIdx_spec, "n_srv_fleets",
                                    "ln_sigmaSrvIdx", sigmaSrvIdx_map)
  input_list <- do_sigmaIdx_mapping(input_list, sigmaSrvIdx_pop_spec, "n_srv_fleets",
                                    "ln_sigmaSrvIdx_pop", sigmaSrvIdx_pop_map)

  input_list <- do_comp_theta_mapping(input_list, comp_prefix = "SrvAge", fleet_field = "n_srv_fleets")
  input_list <- do_comp_theta_mapping(input_list, comp_prefix = "SrvLen", fleet_field = "n_srv_fleets")
  input_list <- do_comp_corr_pars_mapping(input_list, comp_prefix = "SrvAge", fleet_field = "n_srv_fleets")
  input_list <- do_comp_corr_pars_mapping(input_list, comp_prefix = "SrvLen", fleet_field = "n_srv_fleets")

  input_list <- do_comp_theta_mapping(input_list, comp_prefix = "SrvAge", has_pop = TRUE, fleet_field = "n_srv_fleets")
  input_list <- do_comp_theta_mapping(input_list, comp_prefix = "SrvLen", has_pop = TRUE, fleet_field = "n_srv_fleets")
  input_list <- do_comp_corr_pars_mapping(input_list, comp_prefix = "SrvAge", has_pop = TRUE, fleet_field = "n_srv_fleets")
  input_list <- do_comp_corr_pars_mapping(input_list, comp_prefix = "SrvLen", has_pop = TRUE, fleet_field = "n_srv_fleets")

  # Conditional Age-at-Length --------------------------------------------------
  input_list <- setup_caal_source(
    input_list,
    ObsCAAL = ObsSrv_caal,
    UseCAAL = UseSrv_caal,
    ISS_CAAL = ISS_Srv_caal,
    CAAL_LikeType = Srv_caal_LikeType,
    CAAL_Type = Srv_caal_Type,
    fleet_type = "Srv"
  )

  # Print Messages ----------------------------------------------------------
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)
}
