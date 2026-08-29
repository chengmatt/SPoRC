# Stage 1 of 3: model setup
#
# The generated half of the options documentation. vignettes/t_model_options.Rmd
# is written by hand and says why a setting exists and which to reach for; this
# reads the package itself and says what every setting is. Hand-written coverage
# drifts the moment an argument is added, so completeness is generated rather
# than maintained.

#' The setup stages, in the order a model is built
#'
#' @keywords internal
setup_stage_order <- function() {
  c("Setup_Mod_Dim", "Setup_Mod_Rec", "Setup_Mod_Biologicals", "Setup_Mod_Movement",
    "Setup_Mod_Tagging", "Setup_Mod_Catch_and_F", "Setup_Mod_FishIdx_and_Comps",
    "Setup_Mod_SrvIdx_and_Comps", "Setup_Mod_Fishsel_and_Q", "Setup_Mod_Srvsel_and_Q",
    "Setup_Mod_Weighting")
}

#' The package's Rd database, wherever it is being read from
#'
#' Three places this gets called from and each needs a different route. An
#' installed package answers to its name. A source tree answers to its root,
#' which is not the working directory when the caller is a vignette or a test, so
#' the root is walked up to. And \code{pkgload::load_all} shadows the installed
#' help without building its index, which is why the name route is tried and
#' allowed to fail rather than relied on.
#'
#' @return A named list of parsed Rd, empty when none can be found.
#'
#' @keywords internal
rd_database <- function() {
  db <- tryCatch(tools::Rd_db("SPoRC"), error = function(e) NULL)
  if(!is.null(db) && length(db) > 0) return(db)

  dir <- normalizePath(getwd(), mustWork = FALSE)
  for(i in seq_len(5)) {
    if(file.exists(file.path(dir, "DESCRIPTION")) && dir.exists(file.path(dir, "man"))) {
      db <- tryCatch(tools::Rd_db(dir = dir), error = function(e) NULL)
      if(!is.null(db)) return(db)
    }
    parent <- dirname(dir)
    if(identical(parent, dir)) break
    dir <- parent
  }
  list()
}

#' Argument descriptions from a function's own help page
#'
#' Reads the \code{\\arguments} section of the installed (or source) Rd, so the
#' reference and \code{?Setup_Mod_Rec} cannot disagree. One Rd item may document
#' several arguments at once (\code{"a,b,c"}), and each gets the shared text.
#'
#' @param topic Function name.
#' @param db Rd database, from \code{\link[tools]{Rd_db}}.
#'
#' @return Named character vector of descriptions, empty when the topic is absent.
#'
#' @keywords internal
rd_argument_text <- function(topic, db) {
  key <- paste0(topic, ".Rd")
  if(!key %in% names(db)) return(character(0))
  rd <- db[[key]]

  tags <- vapply(rd, function(x) {
    t <- attr(x, "Rd_tag"); if(is.null(t)) NA_character_ else t
  }, character(1))
  i <- which(tags == "\\arguments")
  if(length(i) == 0) return(character(0))

  flat <- function(x) paste(rapply(x, as.character, how = "unlist"), collapse = "")
  items <- Filter(function(x) identical(attr(x, "Rd_tag"), "\\item"), rd[[i[1]]])

  out <- character(0)
  for(it in items) {
    if(length(it) < 2) next
    names_i <- trimws(strsplit(trimws(flat(it[[1]])), ",")[[1]])
    text_i <- gsub("[[:space:]]+", " ", trimws(flat(it[[2]])))
    for(nm in names_i[nzchar(names_i)]) out[nm] <- text_i
  }
  out
}

#' Every argument the setup stages accept
#'
#' Assembled from \code{formals()} and the package's own Rd, so it covers the API
#' as it currently stands rather than as it stood when someone last wrote it
#' down. Regenerating the vignette regenerates this.
#'
#' @param stages Function names to document. Defaults to the eleven setup stages.
#' @param guide Path to the hand-written options guide, checked so the reference
#'   can report which settings it also discusses. \code{NULL} skips the check.
#'
#' @return A data frame with one row per argument: \code{stage}, \code{argument},
#'   \code{default}, \code{description}, and \code{in_guide}.
#'
#' @examples
#' \dontrun{
#' ref <- option_reference()
#' subset(ref, stage == "Tagging")
#' }
#'
#' @export
option_reference <- function(stages = setup_stage_order(), guide = NULL) {

  db <- rd_database()

  guide_text <- if(!is.null(guide) && file.exists(guide)) {
    paste(readLines(guide, warn = FALSE), collapse = "\n")
  } else NA_character_

  rows <- lapply(stages, function(s) {
    f <- tryCatch(formals(get(s, envir = asNamespace("SPoRC"))), error = function(e) NULL)
    if(is.null(f)) return(NULL)
    desc <- rd_argument_text(s, db)
    args <- setdiff(names(f), c("input_list", "...", "verbose"))
    if(length(args) == 0) return(NULL)

    default <- vapply(args, function(a) {
      d <- tryCatch(paste(deparse(f[[a]]), collapse = " "), error = function(e) "")
      if(!nzchar(d)) "required" else d
    }, character(1))

    data.frame(
      stage = sub("^Setup_Mod_", "", s),
      argument = args,
      default = unname(default),
      description = unname(ifelse(args %in% names(desc), desc[args], "")),
      in_guide = if(is.na(guide_text)) NA else
        vapply(args, function(a) grepl(a, guide_text, fixed = TRUE), logical(1)),
      stringsAsFactors = FALSE, row.names = NULL)
  })

  out <- do.call(rbind, rows)
  out[order(match(out$stage, sub("^Setup_Mod_", "", stages))), ]
}
