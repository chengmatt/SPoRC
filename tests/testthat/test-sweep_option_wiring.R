# Sweeps the options that choose a form rather than a sharing structure: the likelihood a data source
# is fit under, how its observations are split, what an index measures, how a catchability is computed.
#
# Two questions of each. Is it connected to anything, which a fixed jnLL cannot answer since an option read
# and never used looks exactly like one that is wired.
#
# And does it stay inside its own data source: these models are evaluated rather than fitted, so a survey
# option that moves a fishery likelihood is reading or writing something that is not its own.

# Which data source an option or a likelihood component belongs to, read off its name. Anything
# matching neither is shared, like recruitment or mortality, and is left out of the check.
wiring_stream_of <- function(x) {
  if(grepl("^Srv|^srv|Srv", x)) return("survey")
  if(grepl("^Fish|^fish|Fish|Catch|Discard|Fmort|dmr|conv_fish", x)) return("fishery")
  NA_character_
}

#' Every non-spec option that names its own legal values
#'
#' @keywords internal
wiring_catalog <- local({
  out <- list()
  for(stage in names(sweep_stage_slot)) {
    args <- names(formals(getExportedValue("SPoRC", stage)))
    args <- grep("LikeType$|_Type$|_type$|_form$|_units$", args, value = TRUE)
    for(a in args) {
      legal <- sweep_legal_specs(stage, a)
      if(length(legal) < 2) next
      out[[length(out) + 1]] <- list(
        stage = stage,
        arg = a,
        legal = legal,
        data_source = wiring_stream_of(a)
      )
    }
  }
  out
})

# Options this model cannot turn on, because they govern data it does not have: numbers at
# age, discards, conditional age-at-length, and the population-specific forms of each.
#
# Listing them keeps the gap visible rather than letting an inert option pass as tested.
wiring_unconfigured <- c(
  grep("AA_", vapply(wiring_catalog, function(x) x$arg, character(1)), value = TRUE),
  grep("_pop_|discard|caal", vapply(wiring_catalog, function(x) x$arg, character(1)),
       value = TRUE, ignore.case = TRUE),
  # the test setup fits ages rather than lengths
  "FishLenComps_LikeType", "FishLenComps_Type", "SrvLenComps_LikeType", "SrvLenComps_Type",
  # only the age setting builds without length data, so there is no second value to compare
  # against here. the length route is driven directly in its own test file instead
  "fish_selex_type", "ret_selex_type", "srv_selex_type"
)

#' Build one option value all the way through, ready to evaluate
#'
#' @keywords internal
wiring_dims <- list(
  n_regions = 2,
  n_sexes = 2,
  n_fish_fleets = 1,
  n_srv_fleets = 1,
  n_yrs = 8,
  n_ages = 5
)

wiring_build <- function(entry, value) {
  sweep_build_with(entry$stage, entry$arg,
                   sweep_format_value(entry$stage, entry$arg, value, wiring_dims),
                   dims = wiring_dims, full = TRUE)
}

#' Whether two option values produce a model that computes anything differently
#'
#' Writing a setting into the data list does not show it is connected, since a likelihood code
#' is stored whether or not the data exist. What counts is a change to what is estimated.
#'
#' @keywords internal
wiring_differs <- function(entry, a, b) {
  ia <- wiring_build(entry, a)
  ib <- wiring_build(entry, b)
  if(inherits(ia, "condition") || inherits(ib, "condition")) return(NA)
  d <- sweep_diff(sweep_structure(ia), sweep_structure(ib))
  if(length(d$map) > 0 || length(d$par) > 0) return(TRUE)

  ca <- wiring_contributions(entry, a)
  cb <- wiring_contributions(entry, b)
  if(is.null(ca) || is.null(cb)) return(NA)
  shared <- intersect(names(ca), names(cb))
  !isTRUE(all.equal(ca[shared], cb[shared], tolerance = 1e-12))
}

#' Per-component likelihood contributions for one option value
#'
#' @return Named numeric vector, or \code{NULL} if the value did not build.
#'
#' @keywords internal
wiring_contributions <- function(entry, value) {
  il <- wiring_build(entry, value)
  if(inherits(il, "condition")) return(NULL)
  fit <- tryCatch(fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE),
                  error = function(e) e)
  if(inherits(fit, "condition")) return(NULL)
  contrib <- jnLL_contributions(fit)
  stats::setNames(contrib$contribution, contrib$component)
}


test_that("the wiring sweep found options to sweep", {
  live <- Filter(function(x) !x$arg %in% wiring_unconfigured, wiring_catalog)
  expect_gt(length(wiring_catalog), 30)
  expect_gt(length(live), 5)
})


test_that("each option changes the model it is supposed to configure", {
  # An option whose every legal value builds an identical model is not connected
  # to anything the model reads.
  problems <- character()

  for(entry in wiring_catalog) {
    if(entry$arg %in% wiring_unconfigured) next
    verdicts <- vapply(entry$legal[-1], function(v) wiring_differs(entry, entry$legal[1], v), logical(1))

    if(all(is.na(verdicts))) {
      problems <- c(problems, sprintf("%s: no legal value could be built, so the option was never exercised",
                                      entry$arg))
    } else if(!any(verdicts %in% TRUE)) {
      problems <- c(problems, sprintf("%s: no legal value changes what the model estimates or evaluates",
                                      entry$arg))
    }
  }

  expect_equal(problems, character(0))
})


test_that("an option only moves the likelihood of its own data source", {
  # these models are evaluated at fixed parameters, so a survey setting has no route to a fishery
  # likelihood or the other way round, and anything crossing is reading the wrong array
  problems <- character()

  for(entry in wiring_catalog) {
    if(entry$arg %in% wiring_unconfigured || is.na(entry$data_source)) next
    other <- if(entry$data_source == "survey") "fishery" else "survey"

    base <- wiring_contributions(entry, entry$legal[1])
    if(is.null(base)) next

    for(v in entry$legal[-1]) {
      alt <- wiring_contributions(entry, v)
      if(is.null(alt)) next

      shared <- intersect(names(base), names(alt))
      foreign <- shared[vapply(shared, function(k) identical(wiring_stream_of(k), other), logical(1))]
      for(k in foreign) {
        if(!isTRUE(all.equal(base[[k]], alt[[k]], tolerance = 1e-12)))
          problems <- c(problems, sprintf("%s = '%s' moved %s from %.10g to %.10g",
                                          entry$arg, v, k, base[[k]], alt[[k]]))
      }
    }
  }

  expect_equal(problems, character(0))
})


test_that("the containment check can tell the data sources apart", {
  # if every component counted as shared, the check above would compare nothing, so both data
  # sources have to be represented among the likelihoods this model reports
  entry <- Filter(function(x) x$arg == "srv_idx_type", wiring_catalog)[[1]]
  contrib <- wiring_contributions(entry, "abd")

  data_sources <- vapply(names(contrib), wiring_stream_of, character(1))
  expect_true(any(data_sources == "survey", na.rm = TRUE))
  expect_true(any(data_sources == "fishery", na.rm = TRUE))
})


test_that("the list of options the sweep cannot reach is still accurate", {
  # An option listed as unreachable that now builds a live model should come off
  # the list so the checks above start covering it.
  became_reachable <- character()

  for(entry in wiring_catalog) {
    if(!entry$arg %in% wiring_unconfigured) next
    for(v in entry$legal[-1]) {
      if(isTRUE(wiring_differs(entry, entry$legal[1], v))) {
        became_reachable <- c(became_reachable, entry$arg)
        break
      }
    }
  }

  expect_equal(unique(became_reachable), character(0))
})
