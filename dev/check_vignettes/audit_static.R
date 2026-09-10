# Purpose: static checks on vignette code and on the report table, for what running does not catch
# Creator: Matthew LH. Cheng
# Date Created: 9/8/26
#
# Usage: Rscript dev/check_vignettes/audit_static.R
#
#   Findings go to dev/scratch/vignette_logs/audit_static.csv. A renamed argument, a
#   rename() that labels the wrong dim, or a subscript short of a dim all run without
#   error, so only these checks see them.

library(here)
root <- here::here()
source(file.path(root, "dev", "check_vignettes", "helper-chunks.R"))
suppressMessages(devtools::load_all(root, helpers = FALSE, quiet = TRUE))

# Reference dimensions --------------------------------------------------------
# d_model_report.Rmd names each report object's dims in order, which is what tells a
# Var2 from a Var3. a small model evaluated once says how many dims each one really has
report_axes <- local({

  txt <- readLines(file.path(root, "vignettes", "d_model_report.Rmd"), warn = FALSE)
  rows <- grep("^[|]", txt, value = TRUE)
  rows <- rows[!grepl("^[|][-: ]+[|]", rows)] # the header rules
  axes <- list()

  for(row in rows) {

    cells <- trimws(strsplit(sub("^[|]", "", sub("[|]$", "", row)), "[|]")[[1]])
    if(length(cells) != 3) next

    name <- gsub("`", "", cells[1])
    dims <- gsub("`", "", cells[3])
    if(!grepl("×", dims) || !grepl("^[A-Za-z0-9_]+$", name)) next

    parts <- trimws(strsplit(dims, "×")[[1]])
    axes[[name]] <- gsub("[^a-z_]", "", tolower(parts)) # (n_years + 1) reads as n_years

  } # end row loop

  axes

})

source(file.path(root, "tests", "testthat", "helper-objective_setup.R"))
reference <- evaluate_input(objective_setup_input(n_lens = 12))
rep_ndim <- vapply(reference$rep, function(x) length(dim(x)), integer(1))

# what a rename() may call each axis
axis_labels <- list(
  n_pop = c("pop", "population"),
  n_regions = c("region", "reg"),
  n_years = c("year", "yr", "yrs"),
  n_yrs = c("year", "yr", "yrs"),
  n_seas = c("season", "seas"),
  n_seasons = c("season", "seas"),
  n_ages = c("age", "ages"),
  n_lens = c("len", "length", "lens"),
  n_sexes = c("sex", "sexes"),
  n_fish_fleets = c("fleet", "fishfleet"),
  n_srv_fleets = c("fleet", "srvfleet"),
  n_fleets = c("fleet")
)

exported <- getNamespaceExports("SPoRC")
findings <- list()

note <- function(vignette, line, check, detail) {
  findings[[length(findings) + 1]] <<- data.frame(vignette = vignette, line = line,
                                                  check = check, detail = detail)
} # end note

# every call in an expression, with the names its arguments were given
calls_in <- function(expr) {

  if(!is.call(expr)) return(list())
  parts <- as.list(expr)
  fn <- parts[[1]]

  # SPoRC::Setup_Mod_Dim() and Setup_Mod_Dim() resolve to the same name
  fn_name <- if(is.name(fn)) as.character(fn)
             else if(is.call(fn) && identical(as.character(fn[[1]]), "::")) as.character(fn[[3]])
             else NA_character_

  out <- list(list(fn = fn_name, args = names(parts)[-1]))
  for(part in parts[-1]) if(!missing(part) && is.call(part)) out <- c(out, calls_in(part))
  out

} # end calls_in

# commas that separate subscripts, so a call inside a subscript does not add to the count
top_level_commas <- function(txt) {

  depth <- 0
  n <- 0

  for(ch in strsplit(txt, "")[[1]]) {
    if(ch %in% c("(", "[")) depth <- depth + 1
    else if(ch %in% c(")", "]")) depth <- depth - 1
    else if(ch == "," && depth == 0) n <- n + 1
  } # end ch loop

  n

} # end top_level_commas

# the text between a bracket at position open and its match
bracketed <- function(txt, open) {

  depth <- 0
  chars <- strsplit(txt, "")[[1]]

  for(i in open:length(chars)) {
    if(chars[i] == "[") depth <- depth + 1
    else if(chars[i] == "]") {
      depth <- depth - 1
      if(depth == 0) return(paste(chars[(open + 1):(i - 1)], collapse = ""))
    }
  } # end i loop

  NA_character_

} # end bracketed

# The report table against the model ------------------------------------------
for(name in names(report_axes)) {

  if(!name %in% names(rep_ndim)) next
  documented <- length(report_axes[[name]])
  actual <- rep_ndim[[name]]
  if(documented != actual) {
    note("d_model_report.Rmd", NA, "report table",
         paste0("rep$", name, " has ", actual, " dims, the table names ", documented))
  } # end mismatch

} # end name loop

# Packaged report objects against the model -----------------------------------
# a saved report is a snapshot. once a quantity gains a dim, vignette code written
# against the snapshot runs on it and fails on a fresh fit of the same model
saved_reports <- grep("_rep$", utils::data(package = "SPoRC")$results[, "Item"], value = TRUE)

for(name in saved_reports) {

  saved <- get(name)
  for(object in intersect(names(saved), names(rep_ndim))) {

    n_saved <- length(dim(saved[[object]]))
    if(n_saved == 0 || rep_ndim[[object]] == 0) next

    if(n_saved != rep_ndim[[object]]) {
      note(paste0("data/", name), NA, "saved report",
           paste0(object, " is ", n_saved, " dims here, ", rep_ndim[[object]], " in the model"))
    } # end mismatch

  } # end object loop
} # end name loop

# The vignette chunks ---------------------------------------------------------
vignettes <- list.files(file.path(root, "vignettes"), pattern = "[.]Rmd$", full.names = TRUE)

for(path in vignettes) {

  vignette <- basename(path)
  chunks <- read_vignette_chunks(path, keep_fragments = TRUE)

  for(chunk in chunks) {

    code <- chunk$code
    line <- chunk$line

    # 1. the chunk parses
    exprs <- tryCatch(parse(text = paste(code, collapse = "\n")), error = function(e) e)
    if(inherits(exprs, "error")) {
      note(vignette, line, "parse", sub("[\r\n].*", "", conditionMessage(exprs)))
      next
    } # end parse failure

    # 2. named arguments against the function's own formals. a function taking ... accepts
    # anything, so there is nothing to check
    for(expr in exprs) {
      for(call in calls_in(expr)) {

        fn_name <- call$fn
        if(is.na(fn_name) || !fn_name %in% exported) next

        formal_names <- names(formals(get(fn_name, envir = asNamespace("SPoRC"))))
        if("..." %in% formal_names) next

        supplied <- call$args
        for(arg in setdiff(supplied[nzchar(supplied)], formal_names)) {
          note(vignette, line, "argument", paste0(fn_name, "(", arg, " = )"))
        } # end arg loop

      } # end call loop
    } # end expr loop

    # 3. a melt() of a report object, then a rename() naming its dims. VarN is the Nth dim,
    # so a label on the wrong VarN plots the wrong axis and never errors
    melt_call <- "melt\\( *([A-Za-z0-9_.]+\\$)*rep\\$([A-Za-z0-9_]+) *\\)"

    for(k in grep(melt_call, code)) {

      object <- sub(".*\\$([A-Za-z0-9_]+) *\\).*", "\\1",
                    regmatches(code[k], regexpr(melt_call, code[k])))
      if(!object %in% names(report_axes)) next

      axes <- report_axes[[object]]
      window <- paste(code[k:min(k + 5, length(code))], collapse = " ")
      renames <- regmatches(window, gregexpr("rename\\([^)]*\\)", window))[[1]]
      pairs <- unlist(regmatches(renames, gregexpr("[A-Za-z_.]+ *= *Var[0-9]+", renames)))

      for(pair in pairs) {

        label <- tolower(gsub("[^A-Za-z]", "", sub(" *=.*", "", pair)))
        index <- as.integer(sub(".*Var", "", pair))

        if(index > length(axes)) {
          note(vignette, line, "melt", paste0("rep$", object, " has ", length(axes),
                                              " dims, ", pair))
          next
        } # end past the end

        allowed <- axis_labels[[axes[index]]]
        if(!is.null(allowed) && !label %in% allowed) {
          note(vignette, line, "melt", paste0(pair, " on rep$", object, ", but dim ",
                                              index, " is ", axes[index]))
        } # end wrong label

      } # end pair loop
    } # end k loop

    # 4. subscripts on a report object. one subscript is flat indexing and legitimate
    for(k in seq_along(code)) {

      found <- gregexpr("[A-Za-z0-9_.]*rep\\$[A-Za-z0-9_]+\\[", code[k])[[1]]
      if(found[1] == -1) next
      widths <- attr(found, "match.length")

      for(j in seq_along(found)) {

        hit <- substring(code[k], found[j], found[j] + widths[j] - 1) # sgl_rg_sable_rep$NAA[
        holder <- sub("\\$.*", "", hit)
        object <- sub("\\[$", "", sub(".*\\$", "", hit))

        # a saved report is measured against itself, since it is a snapshot rather than the
        # current model. a live obj$rep is measured against the model
        n_dims <- if(holder == "rep") rep_ndim[object]
                  else if(exists(holder) && !is.null(get(holder)[[object]])) length(dim(get(holder)[[object]]))
                  else NA_integer_
        if(is.na(n_dims)) next

        inside <- bracketed(code[k], found[j] + widths[j] - 1) # the opening bracket itself
        if(is.na(inside)) next

        inside <- gsub(", *drop *= *[A-Za-z]+", "", inside)
        n_subscripts <- top_level_commas(inside) + 1

        if(n_subscripts > 1 && n_subscripts != n_dims) {
          note(vignette, line, "subscript", paste0(holder, "$", object, "[", inside, "] uses ",
                                                   n_subscripts, " of ", n_dims, " dims"))
        } # end wrong count

      } # end j loop
    } # end k loop

  } # end chunk loop
} # end path loop

# Report ----------------------------------------------------------------------
log_dir <- file.path(root, "dev", "scratch", "vignette_logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
out <- if(length(findings) > 0) do.call(rbind, findings) else
  data.frame(vignette = character(0), line = integer(0), check = character(0), detail = character(0))
write.csv(out, file.path(log_dir, "audit_static.csv"), row.names = FALSE)

message(nrow(out), " finding(s), written to ", file.path(log_dir, "audit_static.csv"))
if(nrow(out) > 0) print(out, row.names = FALSE)
