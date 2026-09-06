# North Sea sandeel area 1r, the smsR bridge. 1983-2021, ages 0-4, two half-year seasons, one region,
# one sex.
#
# The sandeel bundle sits in dev/dev_data/, outside an installed package, so this sources
# make_north_sea_sandeel_bridge_figs.R and stops before it fits. Skips when dev/ is not there.

#' Path to the sandeel bridge script, or NULL when the dev tree is absent
#'
#' @keywords internal
sandeel_script_path <- function() {
  # testthat runs from tests/testthat, R CMD check from a copy of it, and a
  # devtools::load_all session from the package root
  for(up in c("../..", "../../..", ".")) {
    p <- file.path(up, "dev", "make_sporc_obj_figs", "make_north_sea_sandeel_bridge_figs.R")
    d <- file.path(up, "dev", "dev_data", "north_sea_sandeel_1r.rds")
    if(file.exists(p) && file.exists(d)) return(normalizePath(p))
  } # end up loop
  NULL
}

#' Build the sandeel bridge seeded at smsR's maximum likelihood estimate
#'
#' Evaluates the case study script down to the point where it has an input list
#' seeded at smsR's estimate and a report from it, and no further. Everything
#' after that point in the script fits the model and draws figures, which this
#' test does not need.
#'
#' @return A list with the seeded report, the smsR reference series, and the
#'   relative difference in fishing mortality at age over every cell.
#'
#' @keywords internal
build_sandeel_bridge <- function() {

  path <- sandeel_script_path()
  if(is.null(path)) return(NULL)

  src <- readLines(path)

  # everything from the seeded stage's own report onward belongs to the figures
  # and the free fit
  stop_at <- grep('^message\\("seeded at smsR', src)
  if(!length(stop_at)) stop("The sandeel script no longer marks the end of its seeded stage, so this test needs updating")
  src <- src[seq_len(stop_at[1] - 1)]

  # the script sets up its own session and reads its data through here(); the
  # test already has the package loaded and resolves the bundle from the script
  root <- dirname(dirname(dirname(path)))
  drop <- grepl("^library\\(|^devtools::load_all\\(|^source\\(here\\(", src)
  src <- src[!drop]
  src <- sub('readRDS\\(here\\("dev", "dev_data", "north_sea_sandeel_1r.rds"\\)\\)',
             sprintf('readRDS(%s)', deparse(file.path(root, "dev", "dev_data", "north_sea_sandeel_1r.rds"))),
             src)

  env <- new.env(parent = globalenv())
  suppressWarnings(suppressMessages(
    eval(parse(text = paste(src, collapse = "\n")), envir = env)))

  list(
    rep = env$seed$rep,
    ref = env$ref,
    ref_F = env$ref_F,
    ref_N = env$ref_N,
    rel_F = env$rel,
    n_yrs = env$n_yrs,
    natmort = env$natmort,
    M_sms = env$M_sms,
    seasdur = env$seasdur,
    input_list = env$input_list
  )
}
