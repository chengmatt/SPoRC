# Purpose: run one vignette's chunks in a fresh session and report the first failure
# Creator: Matthew LH. Cheng
# Date Created: 9/8/26
#
# Usage: Rscript dev/check_vignettes/run_one.R <vignette.Rmd> <seconds> <result file>

args <- commandArgs(trailingOnly = TRUE)
vignette_path <- normalizePath(args[1])
time_cap <- if(length(args) > 1) as.numeric(args[2]) else Inf
result_path <- if(length(args) > 2) args[3] else ""
if(nzchar(result_path)) result_path <- file.path(normalizePath(dirname(result_path)), basename(result_path))

root <- here::here()
source(file.path(root, "dev", "check_vignettes", "helper-chunks.R"))
suppressMessages(devtools::load_all(root, helpers = FALSE, quiet = TRUE))

# the code runs from vignettes/, the directory knitr would put it in, so a relative
# path to figures/ resolves the same way it does at build
script_path <- tempfile(fileext = ".R")
n_chunks <- write_vignette_script(vignette_path, script_path)
script_lines <- readLines(script_path, warn = FALSE)
setwd(file.path(root, "vignettes"))
grDevices::pdf(NULL) # so a plotting chunk does not leave an Rplots.pdf behind

# a chunk that shows how to install a package must never install it from a check run
installers <- c("install.packages", "update.packages", "remove.packages", "pak", "pkg_install",
                "install_github", "install_git", "install_local", "install_version", "install")

is_install <- function(expr) {

  if(!is.call(expr)) return(FALSE)
  parts <- as.list(expr)
  fn <- parts[[1]]

  name <- if(is.name(fn)) as.character(fn)
          else if(is.call(fn) && identical(as.character(fn[[1]]), "::")) as.character(fn[[3]])
          else ""
  if(name %in% installers) return(TRUE)

  for(part in parts[-1]) if(!missing(part) && is.call(part) && is_install(part)) return(TRUE)
  FALSE

} # end is_install

status <- "ok"
detail <- ""
rmd_line <- NA_integer_
failed_call <- ""
n_skipped <- 0L
warnings_seen <- character(0)
started <- Sys.time()

exprs <- tryCatch(parse(script_path, keep.source = TRUE), error = function(e) e)

if(inherits(exprs, "error")) {
  status <- "parse"
  detail <- gsub("[\r\n]+", " ", conditionMessage(exprs))
  exprs <- expression()
} # end parse failure

srcs <- attr(exprs, "srcref")
env <- new.env(parent = globalenv())

# a warning must not abort the expression: tryCatch(warning = ) would leave every
# later chunk without the objects it needs and report failures that are not there
for(i in seq_along(exprs)) {

  if(is_install(exprs[[i]])) {
    n_skipped <- n_skipped + 1L
    next
  } # end install call

  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  setTimeLimit(elapsed = max(1, time_cap - elapsed), transient = TRUE)

  err <- tryCatch(
    withCallingHandlers(
      {eval(exprs[[i]], envir = env); NULL},
      warning = function(w) {
        warnings_seen <<- c(warnings_seen, gsub("[\r\n]+", " ", conditionMessage(w)))
        invokeRestart("muffleWarning")
      },
      message = function(m) invokeRestart("muffleMessage")
    ),
    error = function(e) conditionMessage(e)
  )

  if(!is.null(err)) {

    # the cap firing inside nlminb aborts the R call to obj$fn, and nlminb reports that as
    # "NA/NaN gradient evaluation" rather than a time limit. anything raised at the cap is
    # a timeout whatever it says, or a slow model reads as a broken one
    at_the_cap <- is.finite(time_cap) &&
                  as.numeric(difftime(Sys.time(), started, units = "secs")) >= 0.98 * time_cap

    status <- if(at_the_cap || grepl("reached elapsed time limit|reached CPU time limit", err)) "timeout" else "error"
    detail <- gsub("[\r\n]+", " ", err)
    if(at_the_cap) detail <- paste("reached the cap, message may be the optimizer's:", detail)
    rmd_line <- chunk_line_of(script_lines, as.integer(srcs[[i]][1]))
    failed_call <- gsub("[\r\n]+", " ", paste(deparse(exprs[[i]]), collapse = " "))
    break
  } # end failure

} # end i loop

setTimeLimit(elapsed = Inf)
minutes <- round(as.numeric(difftime(Sys.time(), started, units = "mins")), 2)

# the distinct warnings go to the log, the count goes in the row. a fit that stops on a
# NaN function value warns where a NaN gradient errors, so a muffled warning hides one
distinct_warnings <- table(warnings_seen)

if(length(distinct_warnings) > 0) {
  cat("\n=== warnings ===\n")
  for(k in order(-distinct_warnings)) {
    cat(distinct_warnings[[k]], "x ", substr(names(distinct_warnings)[k], 1, 200), "\n", sep = "")
  } # end k loop
} # end warnings

row <- sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
               basename(vignette_path), status, minutes,
               if(n_skipped > 0) paste0(n_chunks, " (", n_skipped, " install calls skipped)") else n_chunks,
               if(is.na(rmd_line)) "" else rmd_line,
               substr(failed_call, 1, 160), substr(detail, 1, 400), length(distinct_warnings))

if(nzchar(result_path)) writeLines(row, result_path) else cat(row, "\n")
