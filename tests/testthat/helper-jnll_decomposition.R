# Every term of the joint negative log likelihood, written out separately from the sum in
# model_objective.R, so a term dropped, double counted or weighted on the wrong side of sum() fails.
#
# Each entry names a REPORTed *_nLL object and the weight applied to it:
#
#   weight  name of the weight in input_list$data, or NA when the likelihood enters jnLL unweighted.
#           composition weights are applied inside the likelihood; priors have no weight at all
#   mode    "elementwise" for sum(Wt * nLL), "scalar" for Wt * sum(nLL). these agree on a scalar weight
#           and differ on an array. "at_age" collapses the age and sex dims the weight does not have
#   optional TRUE for a likelihood reported only when that data source is fit. every other term must
#           be reported by every model

jnLL_terms <- list(
  Catch_nLL                    = list(weight = "Wt_Catch",        mode = "elementwise"),
  Catch_pop_nLL                = list(weight = "Wt_Catch_pop",    mode = "elementwise"),
  Discard_nLL                  = list(weight = "Wt_Discard",      mode = "elementwise"),
  Discard_pop_nLL              = list(weight = "Wt_Discard_pop",  mode = "elementwise"),
  FishIdx_nLL                  = list(weight = "Wt_FishIdx",      mode = "elementwise"),
  FishIdx_pop_nLL              = list(weight = "Wt_FishIdx_pop",  mode = "elementwise"),
  SrvIdx_nLL                   = list(weight = "Wt_SrvIdx",       mode = "elementwise"),
  SrvIdx_pop_nLL               = list(weight = "Wt_SrvIdx_pop",   mode = "elementwise"),

  # Age-disaggregated data sources. Each is reported always but summed
  # under the same weight as its aggregated counterpart, and every entry is zero
  # unless that data source is fit, so they are elementwise like the aggregates.
  CatchAA_nLL                  = list(weight = "Wt_Catch",        mode = "at_age"),
  CatchAA_pop_nLL              = list(weight = "Wt_Catch_pop",    mode = "at_age"),
  DiscardAA_nLL                = list(weight = "Wt_Discard",      mode = "at_age"),
  DiscardAA_pop_nLL            = list(weight = "Wt_Discard_pop",  mode = "at_age"),
  SrvIdxAA_nLL                 = list(weight = "Wt_SrvIdx",       mode = "at_age"),
  SrvIdxAA_pop_nLL             = list(weight = "Wt_SrvIdx_pop",   mode = "at_age"),

  FishAgeComps_nLL             = list(weight = NA, mode = "scalar"),
  FishAgeComps_pop_nLL         = list(weight = NA, mode = "scalar"),
  FishLenComps_nLL             = list(weight = NA, mode = "scalar"),
  FishLenComps_pop_nLL         = list(weight = NA, mode = "scalar"),
  FishAgeComps_discard_nLL     = list(weight = NA, mode = "scalar"),
  FishAgeComps_discard_pop_nLL = list(weight = NA, mode = "scalar"),
  FishLenComps_discard_nLL     = list(weight = NA, mode = "scalar"),
  FishLenComps_discard_pop_nLL = list(weight = NA, mode = "scalar"),
  SrvAgeComps_nLL              = list(weight = NA, mode = "scalar"),
  SrvAgeComps_pop_nLL          = list(weight = NA, mode = "scalar"),
  SrvLenComps_nLL              = list(weight = NA, mode = "scalar"),
  SrvLenComps_pop_nLL          = list(weight = NA, mode = "scalar"),

  # the conditional age-at-length terms are reported only when that data source is
  # fit, so their absence is not a dropped likelihood
  Fish_caal_nLL                = list(weight = NA, mode = "scalar", optional = TRUE),
  Srv_caal_nLL                 = list(weight = NA, mode = "scalar", optional = TRUE),

  conv_fish_tag_nLL            = list(weight = "Wt_Tagging", mode = "scalar"),
  Fmort_nLL                    = list(weight = "Wt_F",       mode = "scalar"),
  dmr_nLL                      = list(weight = "Wt_D",       mode = "scalar"),
  Rec_nLL                      = list(weight = "Wt_Rec",     mode = "elementwise"),
  Init_Rec_nLL                 = list(weight = "Wt_Init_Rec", fallback = "Wt_Rec", mode = "elementwise"),
  Rec_level_nLL                = list(weight = NA,          mode = "scalar"),
  Init_Sex_nLL                 = list(weight = NA,          mode = "scalar"),
  SR_pen_nLL                   = list(weight = NA,          mode = "scalar"),

  NAA_state_nLL                = list(weight = NA, mode = "scalar"),
  sel_nLL                      = list(weight = NA, mode = "scalar"),
  growth_tv_nLL                = list(weight = NA, mode = "scalar"),
  growth_semipar_nLL           = list(weight = NA, mode = "scalar"),
  rinit_nLL                    = list(weight = NA, mode = "scalar"),
  M_nLL                        = list(weight = NA, mode = "scalar"),
  R0_nLL                       = list(weight = NA, mode = "scalar"),
  h_nLL                        = list(weight = NA, mode = "scalar"),
  Movement_nLL                 = list(weight = NA, mode = "scalar"),
  TagRep_nLL                   = list(weight = NA, mode = "scalar"),
  fish_q_nLL                   = list(weight = NA, mode = "scalar"),
  srv_q_nLL                    = list(weight = NA, mode = "scalar"),
  rec_prop_nLL                 = list(weight = NA, mode = "scalar")
)


#' Pull the report and data list out of whatever the caller passed
#'
#' Accepts a fitted object from fit_model (which attaches $rep and $data), or a
#' plain list(rep = , data = ) for models built with MakeADFun directly.
#'
#' @keywords internal
jnLL_rep_and_data <- function(model) {
  rep <- model$rep
  data <- model$data
  if(is.null(rep)) stop("model has no $rep; call fit_model() or pass list(rep = obj$report(...), data = data)")
  if(is.null(data)) stop("model has no $data; pass list(rep = ..., data = ...) if the object was built with MakeADFun directly")
  if(is.null(rep$jnLL)) stop("model$rep has no jnLL")
  list(rep = rep, data = data)
}


#' Per-component contributions to jnLL
#'
#' Recomputes each term of the jnLL sum from the reported likelihood components
#' and the weights in the data list, following \code{jnLL_terms}.
#'
#' @param model Fitted object from \code{fit_model}, or \code{list(rep =, data =)}.
#'
#' @return Data frame with one row per component: the component name, the weight
#'   applied, and that component's contribution to jnLL.
#'
#' @keywords internal
jnLL_contributions <- function(model) {
  bits <- jnLL_rep_and_data(model)
  rep <- bits$rep
  data <- bits$data

  rows <- lapply(names(jnLL_terms), function(term_name) {
    spec <- jnLL_terms[[term_name]]
    component <- rep[[term_name]]
    if(is.null(component)) return(NULL)

    # A weight the objective backfills for older input lists is absent from a
    # stored data list, and there the fallback weight is what actually controls
    # this component, so the effective NAME changes too and not just its value.
    wt_name <- spec$weight
    if(!is.na(wt_name[1]) && is.null(data[[wt_name]]) && !is.null(spec$fallback)) wt_name <- spec$fallback
    wt <- if(is.na(wt_name[1])) 1 else data[[wt_name]]
    if(is.null(wt)) stop("weight '", wt_name, "' listed in jnLL_terms for '", term_name, "' is absent from the data list")
    spec$weight <- wt_name

    contribution <- if(spec$mode == "at_age") {
      # An at-age component has age and sex dimensions the weight does not,
      # so both are summed within a cell before the cell's weight is applied,
      # which is what the objective does.
      nd <- length(dim(component))
      sum(wt * apply(component, seq_len(nd)[-c(nd - 2, nd - 1)], sum))
    } else if(spec$mode == "elementwise") sum(wt * component) else wt * sum(component)

    # Every term contributes one number to jnLL. More than one means the mode
    # recorded above is wrong for this weight's shape: "scalar" against an array
    # weight returns one value per weight element, which data.frame() would then
    # recycle into one row per element and count the component many times over.
    if(length(contribution) != 1)
      stop("'", term_name, "' resolved to ", length(contribution), " contributions rather than one. Its mode in jnLL_terms is '",
           spec$mode, "' but '", wt_name, "' has length ", length(wt),
           ", so check how ", term_name, " enters the jnLL sum in R/model_objective.R.")

    data.frame(component = term_name,
               weight = if(is.na(spec$weight[1])) NA_character_ else spec$weight,
               contribution = as.numeric(contribution),
               stringsAsFactors = FALSE)
  })

  do.call(rbind, rows)
}


#' Assert that jnLL equals the weighted sum of its reported components
#'
#' Two things are checked. First that the term list above still covers every
#' likelihood the model reports, so a newly added likelihood cannot go
#' unaccounted for. Second that summing those likelihoods with their weights
#' reproduces the reported jnLL.
#'
#' @param model Fitted object from \code{fit_model}, or \code{list(rep =, data =)}.
#' @param tolerance Relative tolerance for the reconstruction.
#' @param label Name for this test setup, used in failure messages.
#'
#' @keywords internal
expect_jnLL_decomposes <- function(model, tolerance = 1e-8, label = deparse(substitute(model))) {
  bits <- jnLL_rep_and_data(model)
  rep <- bits$rep

  reported <- grep("_nLL$", names(rep), value = TRUE)
  problems <- character()
  contributions <- NULL

  unaccounted <- setdiff(reported, names(jnLL_terms))
  if(length(unaccounted) > 0) {
    problems <- c(problems, sprintf(
      "SPoRC_rtmb reports likelihood(s) missing from jnLL_terms: %s. Add them to helper-jnll_decomposition.R and confirm they enter the jnLL sum in R/model_objective.R.",
      paste(unaccounted, collapse = ", ")))
  }

  optional <- names(jnLL_terms)[vapply(jnLL_terms, function(x) isTRUE(x$optional), logical(1))]
  stale <- setdiff(setdiff(names(jnLL_terms), reported), optional)
  if(length(stale) > 0) {
    problems <- c(problems, sprintf(
      "jnLL_terms lists likelihood(s) that SPoRC_rtmb no longer reports: %s. Either the REPORT was dropped or jnLL_terms is out of date.",
      paste(stale, collapse = ", ")))
  }

  # adding up a term list that no longer matches what the model reports says
  # nothing useful, so only sum once the two sides agree on the likelihoods
  if(length(problems) == 0) {
    contributions <- jnLL_contributions(model)
    reconstructed <- sum(contributions$contribution)
    reported_jnLL <- as.numeric(rep$jnLL)

    if(!isTRUE(all.equal(reconstructed, reported_jnLL, tolerance = tolerance))) {
      problems <- c(problems, sprintf(
        "jnLL does not equal the weighted sum of its reported components (reported %.10g, reconstructed %.10g, residual %g). Contributions:\n%s",
        reported_jnLL, reconstructed, reported_jnLL - reconstructed,
        paste(utils::capture.output(print(contributions[contributions$contribution != 0, ], row.names = FALSE)),
              collapse = "\n")))
    }
  }

  testthat::expect(
    length(problems) == 0,
    paste0(label, ": ", paste(problems, collapse = "\n"))
  )

  invisible(contributions)
}
