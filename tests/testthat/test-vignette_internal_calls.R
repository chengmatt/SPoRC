# A vignette is knitted against the installed package, so an evaluated chunk that calls an internal function
# bare fails the build. This reads every vignette's evaluated chunks and names such calls, SPoRC::: excepted.

vignette_dir <- testthat::test_path("..", "..", "vignettes")

# every function defined in R/, and which of them roxygen exports
package_functions <- function() {
  defined <- exported <- character(0)
  for(f in list.files(testthat::test_path("..", "..", "R"), full.names = TRUE)) {
    l <- readLines(f, warn = FALSE)
    at <- which(grepl("^[A-Za-z_][A-Za-z0-9_.]* *(<-|=) *function", l))
    nm <- sub(" *(<-|=).*", "", l[at])
    ex <- vapply(at, function(k) any(grepl("@export", l[max(1, k - 40):k])), TRUE)
    defined <- c(defined, nm)
    exported <- c(exported, nm[ex])
  } # end f loop
  list(defined = unique(defined), exported = unique(exported))
}

# the code of a vignette's evaluated chunks, none when the vignette turns evaluation off for every chunk
evaluated_code <- function(path) {
  lines <- readLines(path, warn = FALSE)
  if(any(grepl("opts_chunk\\$set\\(.*eval *= *FALSE", lines))) return(character(0))
  starts <- grep("^```\\{r", lines)
  ends <- grep("^```$", lines)
  code <- character(0)
  for(s in starts) {
    e <- ends[ends > s][1]
    if(is.na(e) || grepl("eval *= *FALSE", lines[s])) next
    code <- c(code, lines[(s + 1):(e - 1)])
  } # end s loop
  code
}

test_that("evaluated vignette chunks call internal functions with the package prefix", {

  skip_if_not(dir.exists(vignette_dir))
  fns <- package_functions()
  internal <- setdiff(fns$defined, fns$exported)
  bare <- character(0)
  for(path in list.files(vignette_dir, "\\.Rmd$", full.names = TRUE)) {
    code <- evaluated_code(path)
    if(length(code) == 0) next
    own <- sub(" *(<-|=).*", "", code[grepl("^[A-Za-z_][A-Za-z0-9_.]* *(<-|=) *function", trimws(code))]) # defined in the vignette itself
    called <- unlist(regmatches(code, gregexpr("(?<![A-Za-z0-9_.:])[A-Za-z_][A-Za-z0-9_.]*(?= *\\()", code, perl = TRUE)))
    hit <- setdiff(intersect(unique(called), internal), trimws(own))
    if(length(hit) > 0) bare <- c(bare, paste0(basename(path), ": ", paste(hit, collapse = ", ")))
  } # end path loop
  expect_equal(bare, character(0), label = "internal functions called bare in evaluated vignette chunks")

})
