# Sweeps every "<par>_spec" argument in the Setup_Mod_* API across the values it says are legal, and
# checks the map it builds against what its name claims.
#
# Written against the spec string rather than a stored map, so a spec added later is swept the moment it
# lands. Every check is a property of the map alone, so the sweep never builds an AD tape.

# Abbreviation to the dimension it names, as used across the est_shared_ vocabulary.
spec_abbrev_dims <- c(
  p = "pop",
  pop = "pop",
  r = "region",
  y = "year",
  seas = "season",
  s = "sex",
  f = "fleet",
  b = "block",
  x = "sex"
)

#' Dims of a parameter block the map is constant along
#'
#' A map shared over a dimension takes the same value for every cell that differs
#' only along that dimension. Reshaping the map to the parameter's own dim and
#' asking which dims it is constant along recovers what the spec actually did,
#' without this test needing to know the layout of each block.
#'
#' @param m Map factor for one parameter block.
#' @param d Dimensions of that block.
#'
#' @return Logical vector, one per dim, \code{NA} for dims of extent one
#'   where sharing cannot be distinguished from not sharing.
#'
#' @keywords internal
map_shared_dims <- function(m, d) {
  v <- as.character(m)
  if(length(v) != prod(d)) return(NULL)
  arr <- array(v, dim = d)
  vapply(seq_along(d), function(k) {
    if(d[k] < 2) return(NA)
    # constant along dim k when every slice equals the first
    first <- as.vector(apply(arr, seq_along(d)[-k], function(x) x[1]))
    all(vapply(seq_len(d[k]), function(i) {
      idx <- rep(list(bquote()), length(d)); idx[[k]] <- i
      identical(as.vector(do.call(`[`, c(list(arr), idx, list(drop = FALSE)))), first)
    }, logical(1)))
  }, logical(1))
}

#' Estimated-parameter count and shared dims for every block a spec touched
#'
#' @keywords internal
spec_block_profile <- function(il, blocks) {
  stats::setNames(lapply(blocks, function(b) {
    m <- il$map[[b]]
    p <- il$par[[b]]
    if(is.null(m)) return(list(n_free = length(p), shared = NULL))
    lv <- as.character(m)
    list(n_free = length(unique(lv[!is.na(lv)])),
         shared = map_shared_dims(m, if(is.null(dim(p))) length(p) else dim(p)))
  }), blocks)
}

#' Build one spec value in the configuration that makes it live
#'
#' @keywords internal
spec_build <- function(entry, value) {
  d <- sweep_live_dims(entry$arg)
  sweep_build_with(
    entry$stage,
    entry$arg,
    value,
    dims = d,
    extra = sweep_live_config(entry$arg),
    other = sweep_live_other(entry$arg)
  )
}

#' Whether a build outcome is the model legitimately declining this value
#'
#' @keywords internal
spec_declined <- function(res) {
  sweep_is_identifiability_refusal(res) || sweep_is_config_refusal(res)
}

# Specs the sweep cannot currently put into a live configuration. Listing them
# rather than skipping them silently keeps the gap visible: a spec added here
# needs a reason, and a spec that becomes reachable should be removed so the
# checks above start covering it.
sweep_unconfigured <- c(
  # the population-specific data sources need their own _pop data arrays
  # and a natal-homing configuration the sweep test setup does not yet build
  "rho_catch_pop_spec", "rho_discard_pop_spec", "rho_srv_idx_pop_spec", "sigmaCAA_pop_spec", "sigmaDAA_pop_spec",
  "sigmaSrvIdxAA_pop_spec",
  "sigmaFishIdx_pop_spec", "sigmaSrvIdx_pop_spec"
)


# One entry per spec argument, discovered rather than listed, so the sweep covers
# whatever the API currently offers.
sweep_spec_catalog <- local({
  out <- list()
  for(stage in names(sweep_stage_slot)) {
    args <- grep("_spec$", names(formals(getExportedValue("SPoRC", stage))), value = TRUE)
    for(a in args) {
      legal <- sweep_legal_specs(stage, a)
      if(length(legal) < 2) next
      out[[length(out) + 1]] <- list(stage = stage, arg = a, legal = legal)
    }
  }
  out
})


test_that("the sweep found spec arguments to sweep", {
  # a catalog that silently empties turns every test below into a no-op, so the
  # discovery itself is asserted
  expect_gt(length(sweep_spec_catalog), 15)
})


test_that("every legal spec value builds", {
  problems <- character()

  for(entry in sweep_spec_catalog) {
    for(v in entry$legal) {
      res <- spec_build(entry, v)
      if(!inherits(res, "condition") || spec_declined(res)) next
      problems <- c(problems, sprintf("%s = '%s': %s", entry$arg, v,
                                      substr(gsub("\n", " ", conditionMessage(res)), 1, 160)))
    }
  }

  expect_equal(problems, character(0))
})


test_that("each spec value produces a different estimation structure", {
  # A spec whose values all build the same model is not wired to anything. Each
  # value is built in the configuration its parameter needs to exist at all, so a
  # spec reported here is inert where it is supposed to act, not merely inert
  # beside a feature that happens to be off.
  #
  # A spec sharing only over dimensions this test setup has one of is excluded: it
  # equals est_all by arithmetic rather than by any fault in the wiring.
  problems <- character()

  for(entry in sweep_spec_catalog) {
    if(entry$arg %in% sweep_unconfigured || !"est_all" %in% entry$legal) next
    base <- spec_build(entry, "est_all")
    if(inherits(base, "condition")) next
    base_sig <- sweep_signature(base)

    for(v in setdiff(entry$legal, "est_all")) {
      if(sweep_spec_is_degenerate(entry$stage, v, sweep_live_dims(entry$arg))) next
      alt <- spec_build(entry, v)
      if(inherits(alt, "condition")) next
      if(sweep_identical(base_sig, sweep_signature(alt)))
        problems <- c(problems, sprintf("%s: '%s' is indistinguishable from 'est_all'",
                                        entry$arg, v))
    }
  }

  expect_equal(problems, character(0))
})


test_that("'fix' estimates nothing that 'est_all' estimated", {
  problems <- character()

  for(entry in sweep_spec_catalog) {
    if(entry$arg %in% sweep_unconfigured) next
    if(!all(c("fix", "est_all") %in% entry$legal)) next
    est <- spec_build(entry, "est_all")
    fx <- spec_build(entry, "fix")
    if(inherits(est, "condition") || inherits(fx, "condition")) next

    changed <- sweep_diff(sweep_signature(est), sweep_signature(fx))$map
    for(b in changed) {
      m <- fx$map[[b]]
      if(!is.null(m) && any(!is.na(as.character(m))))
        problems <- c(problems, sprintf("%s: block '%s' still has %d free parameters under 'fix'",
                                        entry$arg, b, length(unique(stats::na.omit(as.character(m))))))
    }
  }

  expect_equal(problems, character(0))
})


test_that("sharing never increases the number of estimated parameters", {
  # est_all is the finest partition a spec offers, so no sharing spec may leave a
  # block with more free parameters than it has there.
  problems <- character()

  for(entry in sweep_spec_catalog) {
    if(entry$arg %in% sweep_unconfigured) next
    shared <- grep("^est_shared", entry$legal, value = TRUE)
    if(!"est_all" %in% entry$legal || length(shared) == 0) next
    est <- spec_build(entry, "est_all")
    if(inherits(est, "condition")) next

    for(v in shared) {
      alt <- spec_build(entry, v)
      if(inherits(alt, "condition")) next
      blocks <- sweep_diff(sweep_signature(est), sweep_signature(alt))$map
      pe <- spec_block_profile(est, blocks)
      pa <- spec_block_profile(alt, blocks)
      for(b in blocks) {
        if(pa[[b]]$n_free > pe[[b]]$n_free)
          problems <- c(problems, sprintf("%s: '%s' leaves block '%s' with %d free parameters, more than est_all's %d",
                                          entry$arg, v, b, pa[[b]]$n_free, pe[[b]]$n_free))
      }
    }
  }

  expect_equal(problems, character(0))
})


test_that("two different sharing specs do not collapse the same dim", {
  # est_shared_r and est_shared_s name different dimensions, so on a model whose
  # dimensions all differ they must collapse different dims of the parameter.
  # Identical collapse means at least one of them is wired to the wrong dim,
  # which no pinned map detects as long as the pin was taken from the same wiring.
  problems <- character()

  for(entry in sweep_spec_catalog) {
    if(entry$arg %in% sweep_unconfigured) next
    single <- grep("^est_shared_[a-z]+$", entry$legal, value = TRUE)
    single <- single[sub("^est_shared_", "", single) %in% names(spec_abbrev_dims)]
    single <- single[!vapply(single, function(v) sweep_spec_is_degenerate(entry$stage, v, sweep_live_dims(entry$arg)), logical(1))]
    if(length(single) < 2 || !"est_all" %in% entry$legal) next

    est <- spec_build(entry, "est_all")
    if(inherits(est, "condition")) next

    seen <- list()
    for(v in single) {
      alt <- spec_build(entry, v)
      if(inherits(alt, "condition")) next
      blocks <- sweep_diff(sweep_signature(est), sweep_signature(alt))$map
      prof <- spec_block_profile(alt, blocks)
      # signature of which dims this spec collapsed, per block
      seen[[v]] <- paste(vapply(blocks, function(b)
        paste(b, paste(which(prof[[b]]$shared %in% TRUE), collapse = "/"), sep = ":"),
        character(1)), collapse = " ")
    }

    for(i in seq_along(seen)) for(j in seq_along(seen)) {
      if(j <= i) next
      dims_i <- sub("^est_shared_", "", names(seen)[i])
      dims_j <- sub("^est_shared_", "", names(seen)[j])
      # 's' and 'x' both abbreviate sex, so specs naming the same dimension are
      # expected to coincide
      if(spec_abbrev_dims[[dims_i]] == spec_abbrev_dims[[dims_j]]) next
      if(nzchar(seen[[i]]) && identical(seen[[i]], seen[[j]]))
        problems <- c(problems, sprintf("%s: '%s' and '%s' collapse the same dims (%s)",
                                        entry$arg, names(seen)[i], names(seen)[j], seen[[i]]))
    }
  }

  expect_equal(problems, character(0))
})


test_that("a spec sharing over more dimensions estimates no more than one sharing over fewer", {
  # est_shared_r_s shares everything est_shared_r shares and more, so it cannot
  # leave more free parameters. A spec that parses its dimension list wrongly
  # (dropping a component, or splitting on the wrong separator) shows up here.
  problems <- character()

  for(entry in sweep_spec_catalog) {
    if(entry$arg %in% sweep_unconfigured) next
    shared <- grep("^est_shared_", entry$legal, value = TRUE)
    parts <- lapply(shared, function(v) strsplit(sub("^est_shared_", "", v), "_")[[1]])
    names(parts) <- shared
    if(length(shared) < 2) next

    built <- list()
    for(v in shared) {
      b <- spec_build(entry, v)
      if(!inherits(b, "condition")) built[[v]] <- b
    }
    if(length(built) < 2) next
    ref <- sweep_signature(built[[1]])

    for(coarse in names(built)) for(fine in names(built)) {
      if(coarse == fine || !all(parts[[fine]] %in% parts[[coarse]])) next
      blocks <- union(sweep_diff(ref, sweep_signature(built[[coarse]]))$map,
                      sweep_diff(ref, sweep_signature(built[[fine]]))$map)
      if(length(blocks) == 0) next
      pc <- spec_block_profile(built[[coarse]], blocks)
      pf <- spec_block_profile(built[[fine]], blocks)
      for(b in blocks) {
        if(pc[[b]]$n_free > pf[[b]]$n_free)
          problems <- c(problems, sprintf("%s: '%s' shares over a superset of '%s' but leaves more free parameters in '%s' (%d vs %d)",
                                          entry$arg, coarse, fine, b, pc[[b]]$n_free, pf[[b]]$n_free))
      }
    }
  }

  expect_equal(problems, character(0))
})


test_that("every map is the same length as the parameter block it maps", {
  # RTMB pairs a map factor with a parameter block by position, so a map built at
  # the wrong length either errors inside MakeADFun or, when it happens to divide
  # evenly, recycles and ties together parameters that were meant to be separate.
  problems <- character()

  for(entry in sweep_spec_catalog) {
    for(v in entry$legal) {
      il <- spec_build(entry, v)
      if(inherits(il, "condition")) next
      for(b in names(il$map)) {
        if(is.null(il$par[[b]])) {
          problems <- c(problems, sprintf("map '%s' has no matching parameter block", b))
          next
        }
        if(length(il$map[[b]]) != length(il$par[[b]]))
          problems <- c(problems, sprintf("map '%s' has length %d against a parameter of length %d",
                                          b, length(il$map[[b]]), length(il$par[[b]])))
      }
    }
  }

  expect_equal(unique(problems), character(0))
})




test_that("the list of specs the sweep cannot reach is still accurate", {
  # A spec listed as unreachable that now builds a live model should come off the
  # list, so the checks above begin covering it rather than skipping it forever.
  became_reachable <- character()

  for(entry in sweep_spec_catalog) {
    if(!entry$arg %in% sweep_unconfigured || !"est_all" %in% entry$legal) next
    base <- spec_build(entry, "est_all")
    if(inherits(base, "condition")) next
    base_sig <- sweep_signature(base)
    live <- FALSE
    for(v in setdiff(entry$legal, "est_all")) {
      if(sweep_spec_is_degenerate(entry$stage, v, sweep_live_dims(entry$arg))) next
      alt <- spec_build(entry, v)
      if(inherits(alt, "condition")) next
      if(!sweep_identical(base_sig, sweep_signature(alt))) { live <- TRUE; break }
    }
    if(live) became_reachable <- c(became_reachable, entry$arg)
  }

  expect_equal(became_reachable, character(0))
})
