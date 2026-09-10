# Purpose: read the r chunks out of a vignette so they can be run outside knitr
# Creator: Matthew LH. Cheng
# Date Created: 9/8/26

# every ```{r} block in a vignette, in order. a chunk headed purl = FALSE is a
# documentation fragment that cannot run on its own, and is dropped by default
read_vignette_chunks <- function(path, keep_fragments = FALSE) {

  txt <- readLines(path, warn = FALSE)
  opens <- grep("^```[{][rR]", txt) # chunk headers
  chunks <- vector("list", length(opens))

  for(i in seq_along(opens)) {

    open <- opens[i]
    below <- txt[(open + 1):length(txt)]
    after <- which(grepl("^```[ \t]*$", below))[1] # first bare fence below it, trailing spaces and all
    close <- if(is.na(after)) length(txt) + 1 else open + after
    body <- if(close > open + 1) txt[(open + 1):(close - 1)] else character(0)

    chunks[[i]] <- list(
      line = open, # header line in the Rmd
      header = txt[open],
      fragment = grepl("purl *= *(FALSE|F)\\b", txt[open]),
      code = body
    )

  } # end i loop

  if(!keep_fragments) chunks <- chunks[!vapply(chunks, function(x) x$fragment, logical(1))]
  chunks

} # end read_vignette_chunks

# write one vignette's chunks to a plain R script. each chunk opens with a #@ marker
# giving its line in the Rmd, so a failure is reported against the source and not the copy
write_vignette_script <- function(path, dest) {

  chunks <- read_vignette_chunks(path)
  code <- character(0)

  for(chunk in chunks) code <- c(code, paste0("#@ ", chunk$line), chunk$code, "")

  writeLines(code, dest)
  length(chunks)

} # end write_vignette_script

# the #@ marker above a line of the extracted script, so an error reports an Rmd line
chunk_line_of <- function(script_lines, line) {

  marks <- grep("^#@ ", script_lines)
  above <- marks[marks <= line]
  if(length(above) == 0) return(NA_integer_)
  as.integer(sub("^#@ ", "", script_lines[max(above)]))

} # end chunk_line_of
