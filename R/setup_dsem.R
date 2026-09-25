# Stage 1 of 3: model setup
#
# Links any SPoRC deviation process to covariate time series using dsem arrow lag notation. Holds the
# arrow reader, the table of processes that can be linked, and the setup function itself.

#' Switch TMB's atomic sparse log determinant on or off inside RTMB
#'
#' A moderated arrow puts a random effect in the precision, and TMB's atomic sparse log
#' determinant has no second derivative, so the plain taped factorization has to be used
#' instead. \code{Setup_Mod_DSEM} switches it off for the session when the arrows need it.
#'
#' @param value 0 for the plain taped factorization, 1 for TMB's default.
#'
#' @return The previous value, invisibly, or \code{NA} if RTMB does not expose the flag.
#'
#' @keywords internal
set_dsem_logdet_atomic <- function(value) {

  config_symbol <- tryCatch(getNativeSymbolInfo("TMBconfig", PACKAGE = "RTMB"), error = function(e) NULL)

  if(is.null(config_symbol)) {
    warning("This version of RTMB does not expose TMB's settings, so a moderated dsem that needs the determinant will return NaN.")
    return(invisible(NA_integer_))
  }

  settings_env <- new.env()
  .Call(config_symbol, settings_env, 1L) # 1 reads the settings into the environment
  settings <- eapply(settings_env, as.integer)
  previous <- settings$tmbad.atomic_sparse_log_determinant

  settings$tmbad.atomic_sparse_log_determinant <- as.integer(value)
  .Call(config_symbol, list2env(lapply(settings, as.integer), new.env()), 2L) # 2 writes them back

  return(invisible(previous))

} # end function

#' Map level of every cell of a deviation array
#'
#' The map says which cells are estimated and which are fixed, and an array with no map entry
#' has every cell estimated. Reading it this way gives one vector either way, \code{NA} at a
#' fixed cell, which is what the dsem checks read.
#'
#' @param input_list List with \code{par} and \code{map}.
#' @param par_name Name of the deviation array.
#'
#' @return Integer vector, one level per cell of the array, \code{NA} where fixed.
#'
#' @keywords internal
dev_map_levels <- function(input_list,
                           par_name) {

  if(is.null(input_list$map[[par_name]])) return(seq_along(input_list$par[[par_name]])) # no map, so every cell is estimated

  return(as.integer(input_list$map[[par_name]]))

} # end function

#' Order series so that everything needed within a year comes first
#'
#' @param n_series Number of series.
#' @param edge_from,edge_to Series indices of every dependency inside one year:
#'   the same-year paths, and each moderated arrow's moderating series.
#'
#' @return List with \code{order}, \code{NULL} when the edges form a loop, and
#'   \code{stuck}, the series left in that loop.
#'
#' @keywords internal
peel_series_order <- function(n_series,
                              edge_from,
                              edge_to) {

  left <- seq_len(n_series) # series still to be ordered
  series_order <- integer(0)

  repeat {
    live <- edge_from %in% left & edge_to %in% left # an edge out of a placed series no longer blocks
    ready <- setdiff(left, edge_to[live]) # nothing left pointing at them
    if(length(ready) == 0 || length(left) == 0) break
    series_order <- c(series_order, ready)
    left <- setdiff(left, ready)
  } # end repeat

  if(length(left) > 0) return(list(order = NULL, stuck = left)) # a loop, so no order exists

  return(list(order = series_order, stuck = integer(0)))

} # end function

#' Read dsem arrow and lag notation
#'
#' One arrow per line. \code{"from -> to, lag, name, start"} is a path: series
#' \code{from} in year \code{t - lag} affects series \code{to} in year \code{t}.
#' A lag may be negative, which reads \code{from} in a later year rather than an
#' earlier one. On a time axis that says the future affects the present, so it is
#' for axes that are not time, such as ages or length bins.
#' \code{"a <-> a, 0, name, start"} is the innovation sd of series \code{a}, and
#' \code{"a <-> b, 0, name, start"} a covariance term between two innovations,
#' read as dsem reads it. A name of \code{NA} fixes the arrow at \code{start},
#' and a name used on several arrows is one shared parameter.
#'
#' @param dsem_arrows Arrow lines, one per element of a character vector or one
#'   per line of a string. \code{#} starts a comment.
#' @param variables Character vector of series names the arrows may use.
#' @param covs Character vector of series groups whose innovations are allowed
#'   to covary, each group written as one string, for example \code{"a, b"}.
#'   A covariance term is added for every pair in a group, as dsem's \code{covs}
#'   does. Series with no sd line of their own get one added, also as dsem does.
#' @param mod_var_logscale Whether a moderated sd or covariance is the
#'   exponential of its series. \code{FALSE} (default) reads it on the natural
#'   scale, as dsem does.
#' @param variance What an sd line means. \code{"conditional"} (default) makes
#'   it the innovation sd, so a series' spread is that plus whatever its paths
#'   add: \eqn{s / \sqrt{1 - \rho^2}} under a self path \eqn{\rho}, which increases
#'   without bound under a random walk. \code{"diagonal"} and \code{"marginal"}
#'   make it the marginal sd of the series, with the innovation sd solved for
#'   year by year (\code{\link{get_dsem_matrices}}). The two are the same
#'   without covariance lines. With them, \code{"diagonal"} solves as if the
#'   innovations were independent, so a cell lands near its sd line, and
#'   \code{"marginal"} keeps the innovation correlations the lines imply and
#'   lands every cell on it exactly. These are \code{dsem}'s
#'   \code{constant_variance} settings. Both forms need stationary paths and
#'   do not allow a series with an sd of zero or a moderated sd.
#'
#' @return List with \code{arrows} (one row per arrow: type, from, to, lag,
#'   name, start, par, from_idx, to_idx, mod_idx), \code{variables},
#'   \code{beta_names} (paths and covariances, natural scale),
#'   \code{ln_sd_names} (sds, log scale), \code{mod_var_logscale} and
#'   \code{series_order}, the order the series can be drawn in within a year
#'   (\code{NULL} when same-year paths form a loop), and \code{derived}, the
#'   series whose sd is fixed at zero, worked out from what points into them,
#'   and \code{variance}, what the sd lines mean.
#'
#' @export
read_dsem_arrows <- function(dsem_arrows,
                             variables,
                             covs = NULL,
                             mod_var_logscale = FALSE,
                             variance = "conditional") {

  if(!is.character(variance) || length(variance) != 1 || !variance %in% c("conditional", "diagonal", "marginal")) {
    stop("variance should be 'conditional' (an sd line is the innovation sd), or 'diagonal' or 'marginal' ",
         "(an sd line is the marginal sd of its series).")
  }

  # get number of arrows from dsem map
  dsem_arrows <- paste(dsem_arrows, collapse = "\n") # one string with a line per arrow, or one arrow per element
  arrow_lines <- trimws(sub("#.*$", "", strsplit(dsem_arrows, "\n")[[1]]))
  arrow_lines <- arrow_lines[arrow_lines != ""] # drop blank lines
  n_arrows <- length(arrow_lines)
  if(n_arrows == 0) stop("dsem_arrows has no arrows. Please provide a valid map!")

  # setup dataframe for arrow
  arrows <- data.frame(type = character(n_arrows), # path, sd or covarince
                       from = character(n_arrows), # series the arrow comes from
                       to = character(n_arrows), # series the arrow points at
                       lag = integer(n_arrows), # years between from, to
                       name = character(n_arrows), # parameter name, NA when fixed
                       start = numeric(n_arrows)) # fixed value or starting value when estimated

  # fill in arrow dataframe
  for(i in seq_len(n_arrows)) {

    fields <- trimws(strsplit(arrow_lines[i], ",")[[1]]) # split commas
    fields <- c(fields, rep(NA, 4 - length(fields))) # some cleaning up (should only have 4 entires in a given line)
    two_headed <- grepl("<->", fields[1]) # see if there are double headed arrows
    ends <- trimws(strsplit(fields[1], if(two_headed) "<->" else "->")[[1]]) # figure out from, to

    # some checks
    if(length(ends) != 2 || any(ends == "")) {
      stop("Could not read the arrow in '", arrow_lines[i], "'. ",
           "Use 'from -> to, lag, name, start' or 'a <-> a, 0, name, start'.")
    }

    if(!all(ends %in% variables)) {
      stop("'", arrow_lines[i], "' uses a series that is not defined. Series available: ",
           paste(variables, collapse = ", "), ".")
    }

    if(is.na(fields[2]) || is.na(suppressWarnings(as.integer(fields[2])))) {
      stop("'", arrow_lines[i], "' needs an integer lag in years as its second field.")
    }

    # fill in arrow dataframe
    arrows$type[i] <- if(!two_headed) "path" else if(ends[1] == ends[2]) "sd" else "cov"
    arrows$from[i] <- ends[1]
    arrows$to[i] <- ends[2]
    arrows$lag[i] <- as.integer(fields[2])
    arrows$name[i] <- if(is.na(fields[3]) || fields[3] %in% c("", "NA")) NA else fields[3]
    arrows$start[i] <- suppressWarnings(as.numeric(fields[4]))

  } # end i loop

  # fill in covariances if left out of the map, before the variances, which is dsem's order
  for(group in covs) {

    pair_vars <- trimws(strsplit(group, "[ ,]+")[[1]])

    if(!all(pair_vars %in% variables)) {
      stop("covs names a series that is not defined: ", paste(setdiff(pair_vars, variables), collapse = ", "), ".")
    }

    for(i1 in seq_along(pair_vars)) {
      for(i2 in seq_along(pair_vars)) {

        if(i2 <= i1) next
        if(any(arrows$type == "cov" & arrows$from == pair_vars[i1] & arrows$to == pair_vars[i2] & arrows$lag == 0)) next # already in the map

        arrows <- rbind(arrows, data.frame(type = "cov",
                                           from = pair_vars[i1],
                                           to = pair_vars[i2],
                                           lag = 0L,
                                           name = paste0("C[", pair_vars[i1], ",", pair_vars[i2], "]"),
                                           start = NA_real_))

      } # end i2 loop
    } # end i1 loop
  } # end group loop

  # fill in variances if not provided, as dsem's add.variances() does
  for(series in setdiff(variables, arrows$to[arrows$type == "sd"])) {
    arrows <- rbind(arrows, data.frame(type = "sd",
                                       from = series,
                                       to = series,
                                       lag = 0L,
                                       name = paste0("V[", series, "]"),
                                       start = NA_real_))
  } # end series loop

  # more checks on valid notation
  if(any(abs(arrows$lag) >= 1 & arrows$type != "path")) {
    stop("Two-headed arrows (sds and covariances) cannot be lagged. Give them a lag of 0.")
  }

  if(any(arrows$type == "path" & arrows$lag == 0 & arrows$from == arrows$to)) {
    stop("A series cannot affect itself in the same year. Use a lag of 1 or more for an ",
         "autoregressive path, for example 'env -> env, 1, rho_env'.")
  }

  if(any(is.na(arrows$name) & is.na(arrows$start))) {
    stop("A fixed arrow (name NA) needs a value in its fourth field, for example ",
         "'rec -> rec, 1, NA, 1' for a random walk.")
  }

  if(anyDuplicated(paste(arrows$type, arrows$from, arrows$to, arrows$lag))) {
    stop("The same arrow is given twice. Give each from, to and lag combination one line.")
  }

  if(anyDuplicated(arrows$to[arrows$type == "sd"]) || !setequal(arrows$to[arrows$type == "sd"], variables)) {
    stop("Give every series exactly one sd line ('a <-> a, 0, name, start'). Series: ",
         paste(variables, collapse = ", "), ".")
  }

  if(any(arrows$type == "sd" & is.na(arrows$name) & !(arrows$start >= 0))) stop("A fixed sd cannot be negative.")

  # figure out arrow indexing stuff inside path matrix
  arrows$from_idx <- match(arrows$from, variables) # column of the from series
  arrows$to_idx <- match(arrows$to, variables) # column of the to series
  arrows$mod_idx <- match(arrows$name, variables, nomatch = 0L) # match arrow name to a variable for a moderating / latent var

  # setup moderating var stuff using Kahn's algorithim
  same_yr <- arrows$type == "path" & arrows$lag == 0 # find arrows that impact same year
  moderated <- arrows$mod_idx > 0

  # ordering to make sure things happen in a coherent order... (e.g., env before recruitment, moderating var, then recruitment)
  path_peel <- peel_series_order(length(variables), arrows$from_idx[same_yr], arrows$to_idx[same_yr]) # same-year paths on their own
  full_peel <- peel_series_order(length(variables),
                                 c(arrows$from_idx[same_yr], arrows$mod_idx[moderated]),
                                 c(arrows$to_idx[same_yr], arrows$to_idx[moderated])) # and the moderating edges
  series_order <- full_peel$order # NULL under a same-year loop, which dgmrf handles but the draw cannot

  # a loop the same-year paths alone do not make is one a moderator closes, so a series would set the coefficient on the arrow that sets it
  if(is.null(series_order) && !is.null(path_peel$order)) {
    stop("A moderating series reads a value that the arrow it moderates sets in the same year, through: ",
         paste(variables[full_peel$stuck], collapse = ", "), ".")
  }

  # check par naming stuff
  est_sd <- !is.na(arrows$name) & arrows$mod_idx == 0 & arrows$type == "sd"
  est_beta <- !is.na(arrows$name) & arrows$mod_idx == 0 & arrows$type != "sd"

  if(any(arrows$name[est_sd] %in% arrows$name[est_beta])) {
    stop("A parameter name is used for both an sd and a path or covariance. sds are estimated ",
         "on the log scale, so give them their own names.")
  }

  # figure out parameter indexing stuff (par is 0 if fixed or moderated)
  beta_names <- unique(arrows$name[est_beta])
  ln_sd_names <- unique(arrows$name[est_sd])
  arrows$par <- 0L
  arrows$par[est_beta] <- match(arrows$name[est_beta], beta_names)
  arrows$par[est_sd] <- match(arrows$name[est_sd], ln_sd_names)

  # figure out if any derived vars (i.e., if non-linear quadratic forms etc / nonlinear func of other time-series)
  derived <- logical(length(variables))
  zero_sd <- arrows$type == "sd" & arrows$par == 0 & arrows$mod_idx == 0 & arrows$start == 0
  derived[arrows$to_idx[zero_sd]] <- TRUE

  if(any(derived)) {

    if(is.null(series_order)) {
      stop("A series with an sd of zero takes its value from the series pointing into it, ",
           "so same-year paths cannot form a loop.")
    }

    if(any(arrows$type == "cov")) {
      stop("A series with an sd of zero cannot be given together with a covariance arrow.")
    }

    no_path_in <- setdiff(which(derived), arrows$to_idx[arrows$type == "path"])

    if(length(no_path_in) > 0) {
      stop("These series have an sd of zero and no path into them, so nothing sets their value: ",
           paste(variables[no_path_in], collapse = ", "), ".")
    }

  } # end if any derived series

  # the diagonal form solves each cell's innovation variance from its sd line, so it needs every series to have
  # one number for that line: no derived series (a target of zero) and no sd read off the grid year by year
  if(variance != "conditional") {

    if(any(derived)) {
      stop("variance = '", variance, "' makes each sd line the marginal sd of its series, and a series ",
           "with an sd of zero has no variance to solve for. Give every series an sd, or use ",
           "variance = 'conditional'.")
    }

    if(any(arrows$mod_idx > 0 & arrows$type != "path")) {
      stop("variance = '", variance, "' needs one number per sd line, and a moderated sd changes ",
           "year by year. Use variance = 'conditional' with a moderated sd.")
    }

  } # end if an sd line is a marginal sd

  # return stuff
  return(list(arrows = arrows,
              variables = variables,
              beta_names = beta_names,
              ln_sd_names = ln_sd_names,
              mod_var_logscale = mod_var_logscale,
              series_order = series_order,
              derived = derived,
              variance = variance))

} # end function

#' Grid cells each dsem arrow lands on
#'
#' Cells of the year by series grid are numbered series by series, years within
#' a series, which is the order \code{as.vector()} lays out a \code{[year,
#' series]} matrix. A path from series \code{j} to series \code{k} at lag
#' \code{L} fills cell \code{(t, k)} of row and cell \code{(t - L, j)} of column
#' in \eqn{B}, for every year \code{t > L}.
#'
#' @param dsem_model Output of \code{read_dsem_arrows}.
#' @param n_grid_yrs Number of years in the grid.
#'
#' @return List with \code{n_cells}, \code{n_grid_yrs} and, for \eqn{I - B} and
#'   \eqn{\Gamma}, a numbered sparse template (\code{m}), the entry filling each
#'   stored slot (\code{slot_entry}), the arrow behind each entry
#'   (\code{entry_arrow}, 0 for the diagonal of \eqn{I - B}), the row and
#'   column of each entry (\code{entry_row}, \code{entry_col}) and its year
#'   (\code{entry_yr}), which a moderated arrow reads its value at.
#'   \code{det_is_one} says whether the cells can be ordered so that everything
#'   an arrow comes from is set before what it points to, and
#'   \code{needs_dense_logdet} whether dgmrf will hold a random effect in the
#'   precision.
#'
#' @keywords internal
get_dsem_cells <- function(dsem_model,
                           n_grid_yrs) {

  # get dsem stuff
  arrows <- dsem_model$arrows
  n_cells <- n_grid_yrs * length(dsem_model$variables)

  # I - B starts from its diagonal of ones, entries of Gamma are all arrows
  path_row <- path_col <- path_arrow <- integer(0)
  gamma_row <- gamma_col <- gamma_arrow <- integer(0)

  # set up path matrices
  for(i in seq_len(nrow(arrows))) {

    if(abs(arrows$lag[i]) >= n_grid_yrs) next # the lag reaches past the grid, so the arrow lands nowhere

    # figure out where to put path coefficients
    to_yrs <- if(arrows$lag[i] >= 0) (arrows$lag[i] + 1):n_grid_yrs else 1:(n_grid_yrs + arrows$lag[i]) # to years
    rows <- (arrows$to_idx[i] - 1) * n_grid_yrs + to_yrs # to cells
    cols <- (arrows$from_idx[i] - 1) * n_grid_yrs + to_yrs - arrows$lag[i] # from cells

    if(arrows$type[i] == "path") {
      path_row <- c(path_row, rows)
      path_col <- c(path_col, cols)
      path_arrow <- c(path_arrow, rep(i, length(rows)))
    } else { # not a path coefficient
      gamma_row <- c(gamma_row, rows)
      gamma_col <- c(gamma_col, cols)
      gamma_arrow <- c(gamma_arrow, rep(i, length(rows)))
    } # end if path or two-headed

  } # end i loop

  # setup matrices for dsem to overwrite in the AD tape
  IminusB_entry_row <- c(1:n_cells, path_row) # fill in 1:n_cells b/c I - B

  # get path matrix (labels where coefficients will go in tape)
  IminusB_m <- Matrix::sparseMatrix(i = IminusB_entry_row,
                                    j = c(1:n_cells, path_col),
                                    x = seq_along(IminusB_entry_row),
                                    dims = c(n_cells, n_cells))

  # get sd and covariance matrix (labels where coefficients will go in tape)
  Gamma_m <- Matrix::sparseMatrix(i = gamma_row,
                                  j = gamma_col,
                                  x = seq_along(gamma_row),
                                  dims = c(n_cells, n_cells))

  # figure out if we can fill grid one cell at a time
  left <- seq_len(n_cells)
  repeat {
    live <- path_col %in% left & path_row %in% left
    ready <- setdiff(left, path_row[live])
    if(length(ready) == 0 || length(left) == 0) break
    left <- setdiff(left, ready)
  } # end repeat

  det_is_one <- length(left) == 0 # det(I-B) = 1 implies there is an ordering that allows 1 cell at a time
  has_lead <- any(arrows$lag < 0)
  IminusB_entry_yr <- (IminusB_entry_row - 1) %% n_grid_yrs + 1 # year of a cell, which a moderated arrow reads at
  gamma_entry_yr <- (gamma_row - 1) %% n_grid_yrs + 1

  # a moderated arrow puts a random effect in the precision, and TMB's atomic sparse log determinant
  # has no second derivative. a derived series never reaches dgmrf, so it does not need the switch
  needs_dense_logdet <- any(arrows$mod_idx > 0) && !any(dsem_model$derived)

  # where each matrix's entries sit, which is what the objective fills on the tape
  IminusB <- list(m = IminusB_m,
                  slot_entry = as.integer(IminusB_m@x),
                  entry_arrow = c(rep(0L, n_cells), path_arrow), # 0 marks the diagonal of I - B
                  entry_row = IminusB_entry_row,
                  entry_col = c(1:n_cells, path_col),
                  entry_yr = IminusB_entry_yr)

  Gamma <- list(m = Gamma_m,
                slot_entry = as.integer(Gamma_m@x),
                entry_arrow = gamma_arrow,
                entry_row = gamma_row,
                entry_col = gamma_col,
                entry_yr = gamma_entry_yr)

  return(list(n_cells = n_cells,
              n_grid_yrs = n_grid_yrs,
              IminusB = IminusB,
              Gamma = Gamma,
              has_cov = any(arrows$type == "cov"),
              has_mod = any(arrows$mod_idx > 0),
              has_lead = has_lead,
              det_is_one = det_is_one,
              needs_dense_logdet = needs_dense_logdet))

} # end function

#' Processes a dsem can link to
#'
#' @return List, one entry per process: the \code{label} its series names start
#'   with, the parameter array, the operating model array a drawn series is
#'   written into (\code{sim_par}, the log state becomes the innovation
#'   \code{naa_eta_all}), every dim of the parameter array in order (\code{"Yr"}
#'   marks the year dim), the sigma its own penalty reads, and whether that
#'   penalty reads the array's map mirror, which is what the dsem switches on.
#'
#' @keywords internal
dsem_process_table <- function() {

  list(
    list(
      label = "rec",
      par = "ln_RecDevs",
      sim_par = "ln_RecDevs",
      dim_names = c("Pop", "Region", "Yr"),
      sigma_par = "ln_sigmaR",
      penalty_reads_map = TRUE,
      why_not = ""
    ),
    list(
      label = "growth",
      par = "ln_growth_devs",
      sim_par = "ln_growth_devs",
      dim_names = c("Pop", "Region", "Yr", "Par", "Sex"),
      sigma_par = NA_character_,
      penalty_reads_map = TRUE,
      why_not = ""
    ),
    list(
      label = "growth_semipar",
      par = "ln_growth_semipar_devs",
      sim_par = "ln_growth_semipar_devs",
      dim_names = c("Pop", "Region", "Yr", "Age", "Sex"),
      sigma_par = NA_character_,
      penalty_reads_map = TRUE,
      why_not = ""
    ),
    list(
      label = "NAA",
      par = "ln_NAA",
      sim_par = "naa_eta_all",
      dim_names = c("Pop", "Region", "Yr", "Seas", "Age", "Sex"),
      sigma_par = "ln_sigmaNAA",
      penalty_reads_map = TRUE,
      why_not = ""
    ),
    list(
      label = "move",
      par = "move_devs",
      sim_par = "move_devs",
      dim_names = c("Pop", "From", "To", "Yr", "Seas", "Age", "Sex"),
      sigma_par = NA_character_,
      penalty_reads_map = TRUE,
      why_not = ""
    ),
    list(
      label = "fish_q",
      par = "ln_fish_q_devs",
      sim_par = "ln_fish_q_devs",
      dim_names = c("Region", "Yr", "Fleet"),
      sigma_par = "ln_sigma_fish_q",
      penalty_reads_map = TRUE,
      why_not = ""
    ),
    list(
      label = "srv_q",
      par = "ln_srv_q_devs",
      sim_par = "ln_srv_q_devs",
      dim_names = c("Region", "Yr", "Fleet"),
      sigma_par = "ln_sigma_srv_q",
      penalty_reads_map = TRUE,
      why_not = ""
    )
  )

} # end function

#' Series a dsem can link, and where each one's cells sit
#'
#' One series per deviation array, kept at every index except the year:
#' \code{ln_growth_devs} of a two sex model gives one series per population,
#' region, parameter and sex. An arrow names the ones to link.
#'
#' @param input_list List with \code{data}, \code{par} and \code{map}.
#' @param dsem_processes Labels from \code{dsem_process_table} to offer.
#' @param arrow_text The arrow lines as one string, comments already removed.
#' @param n_grid_yrs Rows of the dsem grid.
#'
#' @return List with \code{offered} (every name the processes could give, with
#'   \code{offered_label}, \code{offered_par} and \code{offered_cell}) and,
#'   for the linked ones, \code{name}, \code{par}, \code{label}, \code{sigma_par},
#'   \code{penalty_reads_map}, \code{grid_row} (the grid rows that array
#'   reaches) and \code{cell} (where each of those rows sits in the flattened
#'   array).
#'
#' @keywords internal
get_dsem_link <- function(input_list,
                          dsem_processes,
                          arrow_text,
                          n_grid_yrs) {

  arrow_text <- paste(arrow_text, collapse = " ") # one string, whether the lines came joined or one per element
  process_table <- dsem_process_table()
  labels <- vapply(process_table, function(x) x$label, "")
  unknown <- setdiff(dsem_processes, labels)

  if(length(unknown) > 0) {
    stop("dsem_processes names processes that cannot be linked: ", paste(unknown, collapse = ", "),
         ". Available: ", paste(labels, collapse = ", "), ".")
  }

  offered <- offered_label <- offered_par <- link_name <- link_par <- link_label <- link_sigma <- link_why <- character(0)
  link_reads_map <- logical(0)
  link_yr_dim <- integer(0)
  link_idx <- list()
  link_row <- link_cell <- offered_cell <- offered_idx <- list()

  for(entry in process_table) {

    if(!entry$label %in% dsem_processes) next
    if(is.null(input_list$par[[entry$par]])) next # this model does not hold that process

    dims <- dim(input_list$par[[entry$par]])

    if(length(dims) != length(entry$dim_names)) {
      stop(entry$par, " has ", length(dims), " dims, and the process table expects ",
           length(entry$dim_names), " (", paste(entry$dim_names, collapse = ", "),
           "). Update dsem_process_table.")
    }

    if(any(dims == 0)) next # an empty dim, so this model does not use the process (dusky's move_devs has no destination)

    yr_dim <- match("Yr", entry$dim_names)
    idx_dims <- setdiff(seq_along(dims), yr_dim) # every dim the series name is kept at
    stride <- cumprod(c(1, dims[-length(dims)])) # step in the flattened array per dim

    # the numbers at age state runs over chosen years, seasons and ages only, and NAA_pred is not
    # defined outside them, so those are the only cells a series can be linked at
    level <- lapply(dims, seq_len)
    live_yrs <- seq_len(dims[yr_dim])

    if(entry$par == "ln_NAA") {
      live_yrs <- input_list$data$naa_re_yrs
      level[[match("Seas", entry$dim_names)]] <- input_list$data$naa_re_seas
      level[[match("Age", entry$dim_names)]] <- input_list$data$naa_re_ages
      if(length(live_yrs) == 0) next # the state is off, so nothing to link
    }

    # one series per combination of the dims that are not years, and one series when there are none
    combos <- if(length(idx_dims) == 0) matrix(integer(0), nrow = 1L) else as.matrix(expand.grid(level[idx_dims]))

    for(k in seq_len(nrow(combos))) {

      idx <- integer(length(dims))
      idx[idx_dims] <- combos[k,]

      # label when every index dim has one level, which is how rec is named
      name <- if(all(dims[idx_dims] == 1)) entry$label else paste0(entry$label, "_", paste(entry$dim_names[idx_dims], combos[k,], sep = "_", collapse = "_"))

      # rows past the array's own years stay latent, and the dsem forecasts them
      grid_row <- live_yrs[live_yrs <= n_grid_yrs] # the state's own years, which need not start at one
      cell <- integer(length(grid_row))

      for(y in seq_along(grid_row)) {
        idx[yr_dim] <- grid_row[y]
        cell[y] <- 1L + sum((idx - 1L) * stride)
      } # end y loop

      offered <- c(offered, name)
      offered_label <- c(offered_label, entry$label)
      offered_par <- c(offered_par, entry$par)
      offered_cell[[length(offered)]] <- cell
      offered_idx[[length(offered)]] <- idx

      if(!grepl(paste0("(^|[^A-Za-z0-9_.])", name, "([^A-Za-z0-9_.]|$)"), arrow_text)) next # no arrow names it

      link_name <- c(link_name, name)
      link_par <- c(link_par, entry$par)
      link_label <- c(link_label, entry$label)
      link_sigma <- c(link_sigma, entry$sigma_par)
      link_reads_map <- c(link_reads_map, entry$penalty_reads_map)
      link_yr_dim <- c(link_yr_dim, yr_dim)
      idx[yr_dim] <- 0L # a placeholder, written per year wherever the cells are built
      link_idx[[length(link_idx) + 1L]] <- idx
      link_why <- c(link_why, entry$why_not)
      link_row[[length(link_row) + 1L]] <- grid_row
      link_cell[[length(link_cell) + 1L]] <- cell

    } # end k loop

  } # end process loop

  return(list(offered = offered,
              offered_label = offered_label,
              offered_par = offered_par,
              offered_cell = offered_cell,
              offered_idx = offered_idx,
              name = link_name,
              par = link_par,
              label = link_label,
              sigma_par = link_sigma,
              penalty_reads_map = link_reads_map,
              why_not = link_why,
              yr_dim = link_yr_dim,
              idx = link_idx,
              grid_row = link_row,
              cell = link_cell))

} # end function

#' Shorten a dsem to fewer grid years
#'
#' A retrospective peel drops years from every array, and the dsem has to follow:
#' the arrows do not change, but the cells they land on do. A linked series' cells
#' move whenever its array has a dim after the year one, so they are rebuilt from
#' the index that series sits at rather than truncated.
#'
#' @param data Data list holding the dsem fields.
#' @param parameters Parameter list, with its deviation arrays already peeled.
#' @param mapping Map list.
#' @param n_grid_yrs Rows the grid should keep.
#'
#' @return List with the peeled \code{data}, \code{parameters} and \code{mapping}.
#'
#' @keywords internal
peel_dsem_years <- function(data,
                            parameters,
                            mapping,
                            n_grid_yrs) {

  if(is.null(data$dsem_model)) return(list(data = data, parameters = parameters, mapping = mapping)) # no dsem to peel

  keep_rows <- seq_len(n_grid_yrs)
  old_rows <- data$dsem_n_grid_yrs
  n_var <- length(data$dsem_var_names)

  # the grid is [year, series] in all of these, so the peel is the same slice
  data$dsem_n_grid_yrs <- n_grid_yrs
  data$dsem_cells <- get_dsem_cells(data$dsem_model, n_grid_yrs) # where each arrow lands on the shorter grid
  data$dsem_cov_obs <- data$dsem_cov_obs[keep_rows,,drop = FALSE]
  if(!is.null(data$dsem_cov_fixed_sd)) data$dsem_cov_fixed_sd <- data$dsem_cov_fixed_sd[keep_rows,,drop = FALSE]
  if(!is.null(data$dsem_x_known)) data$dsem_x_known <- data$dsem_x_known[keep_rows,,drop = FALSE]
  parameters$dsem_x <- parameters$dsem_x[keep_rows,,drop = FALSE]

  # the map is a flat factor, so it goes back to [year, series] before the slice
  if(!is.null(mapping$dsem_x)) {
    mapping$dsem_x <- factor(matrix(as.integer(mapping$dsem_x), old_rows, n_var)[keep_rows,,drop = FALSE])
  }

  # a linked cell is a position in an array that was peeled too, and it only stays put when the year
  # dim is last. ln_RecDevs is; ln_growth_devs is not, so rebuild rather than truncate
  for(s in seq_along(data$dsem_link_par)) {

    dims <- dim(parameters[[data$dsem_link_par[s]]])
    stride <- cumprod(c(1, dims[-length(dims)])) # step in the flattened array per dim
    yr_dim <- data$dsem_link_yr_dim[s]
    idx <- data$dsem_link_idx[[s]]

    grid_row <- data$dsem_link_row[[s]][data$dsem_link_row[[s]] <= n_grid_yrs] # the rows that series had
    cell <- integer(length(grid_row))

    for(y in seq_along(grid_row)) {
      idx[yr_dim] <- grid_row[y]
      cell[y] <- 1L + sum((idx - 1L) * stride)
    } # end y loop

    data$dsem_link_row[[s]] <- grid_row
    data$dsem_link_cell[[s]] <- cell

  } # end s loop

  return(list(data = data, parameters = parameters, mapping = mapping))

} # end function

#' Take the linked deviations out of their own penalty
#'
#' Every deviation penalty reads \code{data$map_<parameter>} rather than the map
#' itself, so zeroing out the linked cells there takes those deviations out of it.
#' Call from \code{sync_dev_map_data}, which rebuilds it before \code{MakeADFun}.
#'
#' @param data Data list holding \code{dsem_link_par} and \code{dsem_link_cell}.
#'
#' @return \code{data} with the linked cells blanked in each mirror.
#'
#' @keywords internal
apply_dsem_link_switch <- function(data) {

  if(is.null(data$dsem_link_par)) return(data) # no dsem, or no linked series

  for(s in seq_along(data$dsem_link_par)) {
    mirror <- paste0("map_", data$dsem_link_par[s])
    if(is.null(data[[mirror]])) next # no mirror, so no penalty to take these cells out of
    data[[mirror]][data$dsem_link_cell[[s]]] <- NA
  } # end s loop

  return(data)

} # end function

#' The sd line of each linked series
#'
#' A series' own sd is the two-headed arrow from it to itself at lag zero. Under
#' \code{RecDevs_model = "dsem"} that value stands in for \code{sigmaR}, so the
#' objective needs to know which arrow it is. A moderated sd changes by year and
#' cannot stand in, so that series is given a 0.
#'
#' @param dsem_model From \code{\link{read_dsem_arrows}}.
#' @param link_col Grid column of each linked series.
#'
#' @return Integer vector, one arrow row per linked series, 0 where none serves.
#'
#' @keywords internal
get_dsem_link_sd_arrow <- function(dsem_model,
                                   link_col) {

  arrows <- dsem_model$arrows
  sd_arrow <- integer(length(link_col))

  for(s in seq_along(link_col)) {
    i <- which(arrows$type == "sd" & arrows$to_idx == link_col[s] & arrows$lag == 0 & arrows$mod_idx == 0)
    if(length(i) == 1) sd_arrow[s] <- i # none when the sd is moderated, which cannot stand in for sigmaR
  } # end s loop

  return(sd_arrow)

} # end function

#' Covariate family and link codes
#'
#' The dsem package's codes and names. Families: fixed 0, gaussian 1 (normal
#' is the same), bernoulli 2 (binomial is the same), poisson 3, Gamma 4 (gamma
#' is the same), gaussian_fixed_sd 5, lognormal 6, tweedie 7. Links: identity
#' 0, log 1, logit 2, cloglog 3. \code{dsem_default_link} gives each family
#' the link its dsem constructor defaults to, except the Gamma, whose stats
#' default (inverse) dsem has no code for, so it gets the log.
#'
#' @param family Family code, for \code{dsem_default_link}.
#'
#' @return Named integer vector, or one link code.
#'
#' @keywords internal
dsem_family_codes <- function() {

  c(fixed = 0, normal = 1, gaussian = 1, bernoulli = 2, binomial = 2, poisson = 3,
    gamma = 4, Gamma = 4, gaussian_fixed_sd = 5, lognormal = 6, tweedie = 7)

} # end function

#' @rdname dsem_family_codes
#' @keywords internal
dsem_link_codes <- function() {

  c(identity = 0, log = 1, logit = 2, cloglog = 3)

} # end function

#' @rdname dsem_family_codes
#' @keywords internal
dsem_default_link <- function(family) {

  c(0, 0, 2, 1, 1, 0, 1, 1)[family + 1] # by family code 0 to 7

} # end function

#' Refuse covariate values outside their family's support
#'
#' @param y Observations with NA for missing years.
#' @param family Family code.
#' @param name Covariate name for the message.
#'
#' @return \code{invisible(NULL)}; stops on a value the family cannot produce.
#'
#' @keywords internal
check_dsem_cov_support <- function(y,
                                   family,
                                   name) {

  y <- y[!is.na(y)]

  support <- switch(as.character(family),
                    "2" = list(ok = all(y %in% c(0, 1)), text = "0 or 1"), # bernoulli
                    "3" = list(ok = all(y >= 0 & y == round(y)), text = "non-negative integers"), # poisson
                    "4" = , "6" = list(ok = all(y > 0), text = "positive"), # gamma and lognormal
                    "7" = list(ok = all(y >= 0), text = "non-negative"), # tweedie, zeros allowed
                    list(ok = TRUE, text = "")) # the rest take any value

  if(!support$ok) {
    stop("Covariate ", name, " has values outside its family's support; ",
         names(dsem_family_codes())[match(family, dsem_family_codes())], " needs ", support$text, ".")
  }

  return(invisible(NULL))

} # end function

#' Starting mean of a covariate on its link scale
#'
#' The observed mean through the link: itself under identity, its log under
#' the log link, its logit or complementary log-log under those (the mean kept
#' inside 0.02 to 0.98 first). A lognormal takes the mean log instead. Zero
#' when nothing is observed or the value is not finite.
#'
#' @param y Observations with NA for missing years.
#' @param family Family code.
#' @param link Link code.
#'
#' @return Scalar.
#'
#' @keywords internal
dsem_cov_link_mean <- function(y,
                               family,
                               link) {

  y <- y[!is.na(y)]
  if(length(y) == 0) return(0)

  mean_obs <- mean(y)
  inside <- min(max(mean_obs, 0.02), 0.98) # a logit or cloglog link needs a mean strictly inside 0 and 1

  # a lognormal's cell is the log median, so its mean log is the start
  if(family == 6 && link == 1) start <- mean(log(y))
  else start <- switch(as.character(link),
                       "1" = if(mean_obs > 0) log(mean_obs) else 0, # log
                       "2" = stats::qlogis(inside), # logit
                       "3" = log(-log(1 - inside)), # cloglog
                       mean_obs) # identity

  if(!is.finite(start)) start <- 0

  return(start)

} # end function

#' Starting log sd of a covariate's observation error
#'
#' Half the observed spread on each family's own scale: the sd for normal, the
#' CV for gamma, the sd of the log for lognormal. Zero (a dispersion of one)
#' for the tweedie and for the families without a spread parameter.
#'
#' @param y Observations with NA for missing years.
#' @param family Family code.
#'
#' @return Scalar log sd.
#'
#' @keywords internal
dsem_cov_sd_start <- function(y,
                              family) {

  y <- y[!is.na(y)]
  if(length(y) < 2) return(0)

  start <- switch(as.character(family),
                  "1" = log(stats::sd(y) / 2), # normal
                  "4" = log(stats::sd(y) / mean(y) / 2), # gamma, on its CV
                  "6" = log(stats::sd(log(y)) / 2), # lognormal, on the log scale
                  0)

  if(!is.finite(start)) start <- 0

  return(start)

} # end function


#' Set up a dynamic structural equation model on any deviation process
#'
#' Deviations and covariate series become the columns of one year by series grid,
#' linked by arrow and lag lines. A linked series takes its penalty from the dsem
#' density instead of SPoRC's own, through its map mirror. There is no separate
#' objective: \code{\link{SPoRC_rtmb}} evaluates the density whenever a dsem is set
#' up. Fit with \code{random = c(<linked arrays>, "dsem_x")}.
#'
#' A linked recruitment cell is a random effect and takes the full lognormal
#' correction whenever its own penalty would take one: its mean drops by half its
#' variance under the arrows, that variance given the covariate values the model is
#' handed (\code{\link{get_dsem_margvar}}), so \eqn{R_0} scales mean recruitment
#' with or without the arrows. A ramp at zero means none, and a nonzero ramp on a
#' linked year is refused.
#'
#' @param input_list List with \code{data}, \code{par} and \code{map}.
#' @param dsem_arrows Arrow lines, one per element of a character vector or one per
#'   line of a string. Series are named for their process and every index dim, for
#'   example \code{rec}, \code{rec_Pop_1_Region_2} or
#'   \code{NAA_Pop_1_Region_1_Seas_1_Age_3_Sex_1}, with the bare label used when
#'   every index dim has one level. Covariates are named by the columns of
#'   \code{dsem_data}. A name of \code{NA} with a value in the fourth field fixes
#'   that arrow, an sd on the natural scale; to fix a parameter that several arrows
#'   share, write the value on each of them.
#' @param dsem_data Data frame with a \code{year} column and one column per
#'   covariate, \code{NA} where a covariate is not observed. \code{NULL} for a dsem
#'   among deviation series alone.
#' @param dsem_processes Processes whose series the arrows may name, from
#'   \code{dsem_process_table}. \code{NULL} (default) takes the processes whose module
#'   declared \code{"dsem"}: \code{RecDevs_model} in \code{\link{Setup_Mod_Rec}};
#'   \code{NAA_re}, \code{growth_tv_model} and \code{growth_semipar} in
#'   \code{\link{Setup_Mod_Biologicals}}; \code{cont_vary_movement} in
#'   \code{\link{Setup_Mod_Movement}}; \code{fish_q_model} in
#'   \code{\link{Setup_Mod_Fishsel_and_Q}}; and \code{srv_q_model} in
#'   \code{\link{Setup_Mod_Srvsel_and_Q}}. It is \code{"rec"} when none did.
#'
#'   A declared process has to have every series with an estimated cell in the arrows,
#'   since its own penalty is off. A process linked without a declaration keeps its
#'   penalty on the cells left out, but its sigma cannot stay estimated once every cell
#'   is linked.
#'
#'   A linked series covers every year of its array, so a cell mapped off (\code{NA} in
#'   the map) inside one is not left out: the dsem reads it at its fixed starting value
#'   and still evaluates its innovation. A recruitment series with such a mix is refused;
#'   for the other processes the mix is reported.
#'
#'   Under CTMC movement the deviations are the year to year part of habitat preference,
#'   so a \code{preference_formula} term that varies over years is refused once a
#'   movement series is linked.
#' @param dsem_family Named character vector, one entry per covariate, the
#'   distribution of its observations given the grid cells. \code{"fixed"}
#'   (default) takes an observed year as the cell itself, known, and leaves a
#'   missing year latent. \code{"normal"} leaves every year latent, each
#'   observation normal about the cell with sd \code{exp(ln_dsem_obs_sd)}.
#'   \code{"bernoulli"} (or \code{"binomial"}) takes 0/1 observations with the cell
#'   the logit of the probability; \code{"poisson"} takes counts with the cell the
#'   log mean; \code{"gamma"} takes positive values with the cell the log mean and
#'   \code{exp(ln_dsem_obs_sd)} the CV, shape \eqn{1/CV^2}; \code{"lognormal"} takes
#'   positive values with the cell the log median and \code{exp(ln_dsem_obs_sd)} the
#'   log-scale sd; \code{"tweedie"} takes non-negative values including zeros with
#'   the cell the log mean, \code{exp(ln_dsem_obs_sd)} the dispersion and
#'   \code{logit_dsem_tweedie_p} the power as \eqn{1 + \mathrm{plogis}(\cdot)}, in
#'   (1, 2) and starting at 1.5; and \code{"gaussian_fixed_sd"} is normal about the
#'   cell with a known sd per observed year from \code{dsem_fixed_sd}. dsem's own
#'   names \code{"gaussian"} and \code{"Gamma"} are accepted. Under a family whose
#'   link is not the identity the arrows, the mean under \code{dsem_mu_spec} and the
#'   grid are all on the link scale, and every latent cell starts at the series
#'   mean. The sd and power parameters are mapped off for the families with none.
#' @param dsem_link Optional named character vector, one entry per covariate, the
#'   link from the cell to the observation's mean: \code{"identity"}, \code{"log"},
#'   \code{"logit"} or \code{"cloglog"}. Defaults to each family's own: identity for
#'   fixed, normal and the fixed-sd normal, logit for bernoulli, log for the rest. A
#'   link applies whatever the family, as dsem does, so a Poisson under the identity
#'   can be handed a negative mean.
#' @param dsem_fixed_sd Data frame with a \code{year} column and one column per
#'   \code{"gaussian_fixed_sd"} covariate, holding that year's known sd. Needed on
#'   every observed year of such a covariate.
#' @param dsem_mu_spec \code{"est"} (default) estimates every covariate's mean,
#'   \code{"fix"} holds them all at the observed mean, a character vector estimates
#'   only those named, and a named numeric vector fixes those covariates at the
#'   values given with the rest at the observed mean. A linked series' mean is
#'   always fixed at zero, since its process already sits under a level parameter.
#' @param covs Passed to \code{read_dsem_arrows}: series groups whose innovations
#'   may covary.
#' @param dsem_delta0_spec \code{"none"} (default), \code{"est"}, or the names of
#'   the series whose first year offset is estimated.
#' @param mod_var_logscale Passed to \code{read_dsem_arrows}: whether a moderated
#'   sd is the exponential of its series.
#' @param dsem_variance Passed to \code{read_dsem_arrows} as \code{variance}.
#'   \code{"conditional"} (default) reads each sd line as the innovation sd;
#'   \code{"diagonal"} and \code{"marginal"} read it as the series' marginal sd,
#'   solving the innovation sd so the series comes out at that spread, and differ
#'   only when covariance lines are present. Under either, a linked recruitment
#'   cell's correction is half its sd line squared, as under the iid penalty,
#'   whatever paths feed it. A random walk cannot be written under them, since its
#'   paths alone carry a cell past any fixed spread after year one, and setup
#'   checks the starting values and refuses.
#'
#' @return \code{input_list} with the dsem data, parameters (\code{dsem_beta},
#'   \code{ln_dsem_sd}, \code{dsem_mu}, \code{dsem_x}, \code{ln_dsem_obs_sd},
#'   \code{logit_dsem_tweedie_p}, \code{dsem_delta0}) and their map.
#'
#' @export
Setup_Mod_DSEM <- function(input_list,
                           dsem_arrows,
                           dsem_data,
                           dsem_processes = NULL,
                           dsem_family = NULL,
                           dsem_link = NULL,
                           dsem_fixed_sd = NULL,
                           dsem_mu_spec = "est",
                           covs = NULL,
                           dsem_delta0_spec = "none",
                           mod_var_logscale = FALSE,
                           dsem_variance = "conditional") {

  messages_list <<- character(0)

  # Options -----------------------------------------------------------------

  years <- input_list$data$years
  n_grid_yrs <- length(years) + input_list$data$n_proj_yrs_devs # model years plus any projected deviation years
  grid_years <- years[1] + 0:(n_grid_yrs - 1) # calendar year of each grid row

  declared <- if(is.null(input_list$data$dsem_declared)) character(0) else input_list$data$dsem_declared # processes whose module said "dsem"
  if(is.null(dsem_processes)) dsem_processes <- if(length(declared) > 0) declared else "rec"

  cov_names <- if(is.null(dsem_data)) character(0) else setdiff(names(dsem_data), "year") # none is allowed: arrows among deviation series alone
  n_cov <- length(cov_names)

  if(is.null(dsem_family)) dsem_family <- stats::setNames(rep("fixed", n_cov), cov_names)
  family_code <- unname(dsem_family_codes()[dsem_family[cov_names]]) # the dsem package's own codes
  link_code <- if(is.null(dsem_link)) dsem_default_link(family_code) else unname(dsem_link_codes()[dsem_link[cov_names]]) # and its links

  # Input Validation --------------------------------------------------------

  if(!is.character(dsem_arrows) || length(dsem_arrows) == 0) {
    stop("dsem_arrows should be arrow lines, one per element of a character vector or one per line ",
         "of a string, for example c('env -> rec, 0, b_env', 'rec <-> rec, 0, sd_rec').")
  }

  dsem_arrows <- paste(dsem_arrows, collapse = "\n")

  if(!is.null(dsem_data) && (!is.data.frame(dsem_data) || !"year" %in% names(dsem_data))) {
    stop("dsem_data should be a data frame with a 'year' column and one column per covariate, ",
         "or NULL for no covariate.")
  }

  if(!is.null(dsem_data) && n_cov == 0) {
    stop("dsem_data has no covariate columns besides 'year'. Pass NULL for a dsem with no covariate.")
  }

  if(anyDuplicated(dsem_data$year)) stop("dsem_data has duplicated years. Give one row per year.")

  if(any(!dsem_data$year %in% grid_years)) {
    stop("dsem_data years should fall within ", min(grid_years), " to ", max(grid_years),
         " (model years plus n_proj_yrs_devs).")
  }

  if(any(!grepl("^[A-Za-z][A-Za-z0-9_.]*$", cov_names))) {
    stop("Covariate names should start with a letter and hold only letters, numbers, '_' and '.'.")
  }

  if(any(is.na(family_code))) {
    stop("dsem_family should name every covariate with one of '",
         paste(names(dsem_family_codes()), collapse = "', '"), "'.")
  }

  if(any(is.na(link_code))) {
    stop("dsem_link should name every covariate with one of '",
         paste(names(dsem_link_codes()), collapse = "', '"), "'.")
  }

  if(any(family_code == 5) && (!is.data.frame(dsem_fixed_sd) || !"year" %in% names(dsem_fixed_sd))) {
    stop("A gaussian_fixed_sd covariate needs dsem_fixed_sd, a data frame with a 'year' column ",
         "and one column of known sds per such covariate.")
  }

  if(!is.null(dsem_fixed_sd) && !all(cov_names[family_code == 5] %in% names(dsem_fixed_sd))) {
    stop("dsem_fixed_sd is missing a column for: ",
         paste(setdiff(cov_names[family_code == 5], names(dsem_fixed_sd)), collapse = ", "), ".")
  }

  # which deviation series the arrows name, and where each one's cells sit
  arrow_text <- paste(sub("#.*$", "", strsplit(dsem_arrows, "\n")[[1]]), collapse = " ")
  link <- get_dsem_link(input_list, dsem_processes, arrow_text, n_grid_yrs)

  if(length(link$name) == 0) {
    stop("No deviation series appears in dsem_arrows. With dsem_processes = c('",
         paste(dsem_processes, collapse = "', '"), "') the series available are: ",
         paste(utils::head(link$offered, 20), collapse = ", "),
         if(length(link$offered) > 20) ", ..." else "", ".")
  }

  if(any(cov_names %in% link$offered)) {
    stop("A covariate column has the name of a deviation series: ",
         paste(intersect(cov_names, link$offered), collapse = ", "), ". Rename the covariate.")
  }

  # the switch blanks map_<parameter>, so a linked deviation is penalized twice unless that mirror
  # both exists and is read by the penalty. the table says which penalties read theirs
  no_mirror <- unique(link$par[!paste0("map_", link$par) %in% names(input_list$data)])

  if(length(no_mirror) > 0) {
    stop("These arrays have no map mirror in the data list, so their own penalty cannot be switched ",
         "off and a linked deviation would be penalized twice, once there and once by the dsem: ",
         paste(no_mirror, collapse = ", "), ".")
  }

  deaf <- which(!link$penalty_reads_map)

  if(length(deaf) > 0) {
    stop("The penalty on ", paste(unique(link$par[deaf]), collapse = ", "), " does not read its map ",
         "mirror, so a linked deviation would be penalized twice, once there and once by the dsem. ",
         unique(link$why_not[deaf])[1], ".")
  }

  if("move_devs" %in% link$par && isTRUE(input_list$data$use_fixed_movement == 1)) {
    stop("A movement series is linked, but movement is fixed (use_fixed_movement = 1), so the model ",
         "never reads move_devs and the dsem would describe deviations that change nothing.")
  }

  # under CTMC movement the deviations are the year to year part of preference, so a preference
  # covariate that varies over years writes that part again, unpenalized and confounded with it
  if("move_devs" %in% link$par) {

    yr_pref <- get_yr_varying_pref_terms(input_list)

    if(length(yr_pref) > 0) {
      stop("A movement series is linked, and these preference_formula terms vary over years: ",
           paste(utils::head(yr_pref, 10), collapse = ", "), if(length(yr_pref) > 10) ", ..." else "",
           ". The deviations already are the year to year part of preference, so the formula would ",
           "fit it a second time with no penalty on it. Keep preference_formula to terms that stay ",
           "the same across years (region, age, depth) and give the year varying covariate to ",
           "dsem_data with an arrow into the movement series.")
    }

  } # end if a movement series is linked

  # a module that said "dsem" handed over every estimated cell of its process, so each of its series with an
  # estimated cell has to be in the arrows, or those cells would have no density at all
  for(lab in intersect(declared, dsem_processes)) {

    required <- character(0)

    for(i in which(link$offered_label == lab)) {
      # growth declares parameter by parameter, so only the series of a parameter that said "dsem" need one
      if(lab == "growth" && !isTRUE(input_list$data$growth_tv_dsem[link$offered_idx[[i]][4]] == 1)) next
      map_levels <- dev_map_levels(input_list, link$offered_par[i])
      if(any(!is.na(map_levels[link$offered_cell[[i]]]))) required <- c(required, link$offered[i])
    } # end i loop

    left_out <- setdiff(required, link$name)

    if(length(left_out) > 0) {
      stop("The ", lab, " module declared 'dsem', so every one of its series needs an arrow, and these ",
           "have none: ", paste(utils::head(left_out, 10), collapse = ", "),
           if(length(left_out) > 10) ", ..." else "",
           ". Add an sd line for each, or take the declaration back.")
    }

  } # end lab loop

  # dont_pen_recdev_first leaves the first years out of the penalty entirely, which a linked series
  # cannot do: the dsem is their density, so those years would still be penalized by it
  if("ln_RecDevs" %in% link$par && isTRUE(input_list$data$dont_pen_recdev_first > 0)) {
    stop("dont_pen_recdev_first is ", input_list$data$dont_pen_recdev_first, " and recruitment is linked ",
         "to the dsem. The setting takes the first years out of SPoRC's recruitment penalty, but the ",
         "dsem density replaces that penalty and still reads every year, so the setting would do ",
         "nothing. Set dont_pen_recdev_first = 0, or leave recruitment out of the arrows.")
  }

  # a linked cell is a random effect under the arrows and takes the full lognormal correction or none, so a
  # ramp on a linked year has nothing to act on
  if("ln_RecDevs" %in% link$par && isTRUE(input_list$data$do_rec_bias_ramp == 1)) {

    ramp <- get_rec_bias_ramp(1, input_list$data$bias_year, dim(input_list$par$ln_RecDevs)[3], input_list$data$max_bias_ramp_fct)
    linked_yrs <- unique(unlist(link$grid_row[link$par == "ln_RecDevs"]))

    if(any(ramp[linked_yrs[linked_yrs <= length(ramp)]] != 0)) {
      stop("The bias ramp is nonzero on recruitment years the dsem describes. Those deviations are ",
           "random effects under the arrows and take the full lognormal correction or none, so the ramp ",
           "has nothing to act on there. Set do_rec_bias_ramp = 0 for the full correction, or move ",
           "bias_year past the linked years for none.")
    }

  } # end if the bias ramp is on

  # dont_est_recdev_last drops the terminal deviations from the array, but the grid runs to the last year,
  # so those rows would be dsem values the population never reads
  if("ln_RecDevs" %in% link$par) {

    n_dropped <- length(input_list$data$years) - (dim(input_list$par$ln_RecDevs)[3] - input_list$data$n_proj_yrs_devs)

    if(isTRUE(n_dropped > 0)) {
      stop("dont_est_recdev_last leaves the last ", n_dropped, " recruitment deviation(s) out of the ",
           "model, and recruitment is linked to the dsem, whose grid runs to the last year. Those rows ",
           "would be estimated by the dsem alone and read by nothing. Set dont_est_recdev_last = 0, or ",
           "leave recruitment out of the arrows.")
    }

  } # end if recruitment is linked

  # a linked series covers every year of its array, so a cell mapped off is not left out the way the iid
  # penalty leaves one out: it enters the series at its fixed starting value and its innovation is still
  # evaluated, anchoring the process to that value. a recruitment year is never structurally zero, so a mix
  # there is refused; the other processes map cells off by design (growth years before the data, say), so
  # the mix is only reported
  for(s in seq_along(link$name)) {

    map_levels <- dev_map_levels(input_list, link$par[s])
    n_off <- sum(is.na(map_levels[link$cell[[s]]]))
    if(n_off == 0 || n_off == length(link$cell[[s]])) next # all estimated, or all fixed, which is checked below

    mix_msg <- paste0(n_off, " of the ", length(link$cell[[s]]), " cells of ", link$name[s],
                      " are mapped off (NA in map$", link$par[s], "), and the series is linked to the ",
                      "dsem. A linked series covers every year, so those cells are not left out: the ",
                      "dsem reads them at their fixed starting values and still evaluates their ",
                      "innovations, anchoring the process there.")

    if(link$par[s] == "ln_RecDevs") stop(mix_msg, " Estimate every year of the series, or leave it out of the arrows.")
    collect_message(mix_msg)

  } # end s loop

  variables <- c(cov_names, link$name) # grid column order: covariates, then linked series
  n_var <- length(variables)
  link_col <- match(link$name, variables)
  dsem_model <- read_dsem_arrows(dsem_arrows, variables, covs = covs, mod_var_logscale = mod_var_logscale, variance = dsem_variance)

  # a named numeric vector fixes those covariates' means at the values given (a known reference level, or a
  # self test's truth), the rest at the observed mean; character forms say which means are estimated
  mu_given <- NULL

  if(is.numeric(dsem_mu_spec)) {

    if(is.null(names(dsem_mu_spec)) || !all(names(dsem_mu_spec) %in% cov_names)) {
      stop("A numeric dsem_mu_spec fixes covariate means at the values given and needs their names. ",
           "Covariates here: ", paste(cov_names, collapse = ", "), ".")
    }

    mu_given <- dsem_mu_spec
    dsem_mu_spec <- "fix"

  } # end if the means are given as numbers

  if(!is.character(dsem_mu_spec)) {
    stop("dsem_mu_spec should be 'est', 'fix', the names of the covariates whose means are estimated, ",
         "or a named numeric vector of means to fix at.")
  }

  mu_est <- if(identical(dsem_mu_spec, "est")) cov_names else if(identical(dsem_mu_spec, "fix")) character(0) else dsem_mu_spec

  if(!all(mu_est %in% cov_names)) {
    stop("dsem_mu_spec names series that are not covariates: ",
         paste(setdiff(mu_est, cov_names), collapse = ", "), ". Covariates here: ",
         paste(cov_names, collapse = ", "), ".")
  }

  if(!is.character(dsem_delta0_spec)) {
    stop("dsem_delta0_spec should be 'none', 'est', or the names of the series whose first year offset ",
         "is estimated.")
  }

  delta0_est <- if(identical(dsem_delta0_spec, "est")) variables else if(identical(dsem_delta0_spec, "none")) character(0) else dsem_delta0_spec

  if(!all(delta0_est %in% variables)) {
    stop("dsem_delta0_spec names series the arrows do not use: ",
         paste(setdiff(delta0_est, variables), collapse = ", "), ". Series here: ",
         paste(variables, collapse = ", "), ".")
  }

  # a linked cell sharing a map level with another cell would put one parameter under two densities,
  # and a series with no estimated cell at all has nothing for the dsem to describe. cells fixed by
  # the map inside an otherwise estimated series enter as known values, as dsem's fixed family does
  for(s in seq_along(link$name)) {

    map_levels <- dev_map_levels(input_list, link$par[s])
    level_count <- table(map_levels) # how many cells of the whole array sit on each level
    levels_here <- map_levels[link$cell[[s]]]

    if(all(is.na(levels_here))) {
      stop("Every deviation of ", link$name[s], " is fixed by the map, so there is nothing for the dsem ",
           "to describe. Estimate them, or leave the series out of the arrows.")
    }

    if(any(is.na(levels_here))) {
      collect_message(link$name[s], ": ", sum(is.na(levels_here)), " of ", length(levels_here),
                      " cells are fixed by the map and enter the dsem as known values.")
    }

    levels_here <- stats::na.omit(levels_here)

    if(any(level_count[as.character(levels_here)] > 1)) {
      stop("Deviations of ", link$name[s], " share map levels with other cells. Shared deviations are ",
           "not supported for a linked series.")
    }

  } # end s loop

  # Populate Data List ------------------------------------------------------

  cov_obs <- matrix(NA_real_, n_grid_yrs, n_cov, dimnames = list(grid_years, cov_names)) # [year, covariate], NA where unobserved
  for(k in seq_len(n_cov)) cov_obs[match(dsem_data$year, grid_years), k] <- dsem_data[[cov_names[k]]]

  input_list$data$dsem_model <- dsem_model # arrows read from the lines
  input_list$data$dsem_cells <- get_dsem_cells(dsem_model, n_grid_yrs) # where each arrow lands

  # a moderated arrow puts a random effect in the precision, and TMB's atomic sparse log determinant
  # has no second derivative. the flag is read when MakeADFun builds the tape, so set it here
  if(input_list$data$dsem_cells$needs_dense_logdet) {
    set_dsem_logdet_atomic(0)
    collect_message("A moderated arrow puts a random effect in the precision, so TMB's atomic sparse ",
                    "log determinant is switched off for the rest of the session. ",
                    "Put it back with set_dsem_logdet_atomic(1).")
  }

  input_list$data$dsem_n_grid_yrs <- n_grid_yrs
  input_list$data$dsem_var_names <- variables
  for(k in seq_len(n_cov)) check_dsem_cov_support(cov_obs[,k], family_code[k], cov_names[k]) # each family's support, on the observed years

  # a known sd per observed year for the fixed-sd normal, on the grid like the observations
  cov_fixed_sd <- matrix(NA_real_, n_grid_yrs, n_cov, dimnames = list(grid_years, cov_names))

  for(k in which(family_code == 5)) {

    cov_fixed_sd[match(dsem_fixed_sd$year, grid_years),k] <- dsem_fixed_sd[[cov_names[k]]]
    seen <- !is.na(cov_obs[,k]) # the years this covariate is observed in

    if(any(is.na(cov_fixed_sd[seen,k]) | cov_fixed_sd[seen,k] <= 0)) {
      stop("Covariate ", cov_names[k], " is gaussian_fixed_sd and needs a positive sd in dsem_fixed_sd ",
           "for every observed year.")
    }

  } # end k loop

  input_list$data$dsem_cov_obs <- cov_obs
  input_list$data$dsem_cov_link <- link_code
  input_list$data$dsem_cov_fixed_sd <- cov_fixed_sd
  input_list$data$dsem_cov_var_idx <- match(cov_names, variables) # grid column of each covariate
  input_list$data$dsem_cov_family <- family_code
  input_list$data$dsem_link_par <- link$par # parameter array each linked series reads
  input_list$data$dsem_link_col <- link_col # its grid column
  input_list$data$dsem_link_row <- link$grid_row # the grid rows that array reaches
  input_list$data$dsem_link_cell <- link$cell # where those rows sit in the flattened array
  input_list$data$dsem_link_idx <- link$idx # the index each series sits at, so a resize can rebuild the cells
  input_list$data$dsem_link_yr_dim <- link$yr_dim # which dim of that array is years
  input_list$data$dsem_link_sd_arrow <- get_dsem_link_sd_arrow(dsem_model, link_col) # each series' own sd line, 0 when moderated

  # a derived series has no innovation, so its value is whatever the arrows give it. leaving the deviation
  # parameter free would put a second, unpenalized copy of it behind the computed one
  fillable <- q_dev_par_names() # the arrays the objective works out before the observations

  for(s in which(dsem_model$derived[link_col])) {

    par_name <- link$par[s]

    if(!par_name %in% fillable) {
      stop(link$name[s], " has an sd of zero, which makes it a function of its covariates rather than ",
           "a deviation with a density. ", par_name, " is built before the dsem runs, so nothing would ",
           "write that function into it and the deviations would stay at zero. Give the series an sd ",
           "line, or use a zero sd only on a catchability series.")
    }

    map_par <- dev_map_levels(input_list, par_name)
    map_par[link$cell[[s]]] <- NA # the computed cells are no longer parameters
    input_list$map[[par_name]] <- factor(renumber_map_levels(array(map_par, dim = dim(input_list$par[[par_name]]))))
    collect_message(link$name[s], " has an sd of zero, so its deviations are worked out from the arrows and the parameter is fixed.")

  } # end s loop

  # a declared recruitment series hands its sd to sigmaR, which the initial age deviations read, so it needs one of its own
  if("rec" %in% input_list$data$dsem_declared) {

    rec_series <- which(link$par == "ln_RecDevs")

    if(any(dsem_model$derived[link_col[rec_series]])) {
      stop("A recruitment series declared 'dsem' has its sd fixed at zero, so nothing stands in for ",
           "sigmaR, which the initial age deviations read. Give it an sd line.")
    }

    for(s in rec_series) {
      if(input_list$data$dsem_link_sd_arrow[s] == 0) {
        collect_message(link$name[s], "'s sd is moderated, so it cannot stand in for sigmaR and ln_sigmaR's own value is read for it.")
      }
    } # end s loop

  } # end if recruitment is declared

  input_list$data$dsem_delta0_use <- as.integer(length(delta0_est) > 0)

  # Populate Parameter List -------------------------------------------------

  arrows <- dsem_model$arrows
  beta_start <- arrows$start[match(dsem_model$beta_names, arrows$name)]
  sd_start <- arrows$start[match(dsem_model$ln_sd_names, arrows$name)]

  series_mean <- rep(0, n_var) # a linked series sits under its process' own level parameter

  # each covariate starts at its observed mean through its link, which is the scale the arrows work on
  for(k in seq_len(n_cov)) {
    series_mean[match(cov_names[k], variables)] <- dsem_cov_link_mean(cov_obs[,k], family_code[k], link_code[k])
  } # end k loop

  if(!is.null(mu_given)) series_mean[match(names(mu_given), variables)] <- mu_given

  input_list$par$dsem_beta <- ifelse(is.na(beta_start), 0, beta_start) # paths start at no effect
  input_list$par$ln_dsem_sd <- log(ifelse(is.na(sd_start), 1, sd_start)) # sds start at one
  input_list$par$dsem_mu <- series_mean
  input_list$par$dsem_x <- matrix(rep(series_mean, each = n_grid_yrs), n_grid_yrs, n_var) # [year, series], one column per series

  for(k in seq_len(n_cov)) {
    seen <- !is.na(cov_obs[,k]) & family_code[k] %in% c(0, 1, 5) & link_code[k] == 0 # a link-scale covariate starts every cell at its mean instead
    input_list$par$dsem_x[seen,match(cov_names[k], variables)] <- cov_obs[seen,k]
  } # end k loop

  # observation error starts at half the observed spread, on each family's own scale
  input_list$par$ln_dsem_obs_sd <- rep(0, n_cov)
  for(k in seq_len(n_cov)) input_list$par$ln_dsem_obs_sd[k] <- dsem_cov_sd_start(cov_obs[,k], family_code[k])
  input_list$par$logit_dsem_tweedie_p <- rep(0, n_cov) # power 1.5
  input_list$par$dsem_delta0 <- rep(0, n_var)

  # under the diagonal form the innovation variance is what the sd line leaves after the paths, and a path can
  # leave nothing: a random walk's cells sit at or past any fixed spread from year two on. check at the start
  if(dsem_variance != "conditional") {

    start_parts <- get_dsem_matrices(input_list$par$dsem_beta, input_list$par$ln_dsem_sd, dsem_model, input_list$data$dsem_cells)

    # every cell's variance under the solved innovations. an innovation sd of zero or NaN makes the
    # precision singular, so the solve is guarded and a failure reads as no room left
    start_var <- tryCatch({
      Q_start <- Matrix::t(start_parts$IminusB) %*% start_parts$Vinv %*% start_parts$IminusB # precision at the starting values
      as.numeric(Matrix::diag(Matrix::solve(Q_start))) # its inverse diagonal is each cell's variance
    }, error = function(e) rep(NaN, n_grid_yrs * n_var))
    if(!is.null(start_parts$sd_cell)) start_var[!is.finite(as.numeric(start_parts$sd_cell)) | as.numeric(start_parts$sd_cell) <= 0] <- NaN

    no_room <- which(!is.finite(start_var))

    if(length(no_room) > 0) {
      where <- paste0(variables[(no_room - 1) %/% n_grid_yrs + 1], " in grid year ", (no_room - 1) %% n_grid_yrs + 1)
      stop("dsem_variance = '", dsem_variance, "' reads each sd line as the series' marginal sd, and at ",
           "the starting values the paths alone carry these cells to or past that spread, so no ",
           "innovation variance is left: ", paste(utils::head(where, 5), collapse = "; "),
           if(length(where) > 5) "; ..." else "", ". A random walk cannot be written under this form. ",
           "Use an autoregressive path with |rho| < 1, or dsem_variance = 'conditional'.")
    }

    # what each series' sd line asks for, which the marginal form has to land on
    is_sd <- arrows$type == "sd"
    est_line <- !is.na(arrows$name[is_sd]) # an estimated line reads its parameter, a fixed one its own value
    sd_line <- arrows$start[is_sd]
    sd_line[est_line] <- exp(input_list$par$ln_dsem_sd)[arrows$par[is_sd][est_line]]
    target_var <- rep(sd_line[match(variables, arrows$to[is_sd])]^2, each = n_grid_yrs)
    miss <- max(abs(start_var - target_var) / target_var)

    if(dsem_variance == "marginal" && miss > 1e-6) {
      stop("dsem_variance = 'marginal' did not land every cell on its sd line at the starting values ",
           "(largest relative miss ", signif(miss, 3), "). The covariance lines are too strong for the ",
           "fixed point to settle; use dsem_variance = 'diagonal' or 'conditional'.")
    }

    approx_note <- "" # the diagonal form only lands on the sd lines exactly without covariance lines
    if(dsem_variance == "diagonal" && any(arrows$type == "cov")) {
      approx_note <- paste0(" (with covariance lines the diagonal form is approximate: largest relative ",
                            "miss at the start ", signif(miss, 3), ")")
    }

    collect_message("dsem_variance = '", dsem_variance, "': each sd line is the marginal sd of its ",
                    "series, and the innovation sd is solved for year by year", approx_note, ".")

  } # end if an sd line is a marginal sd

  # Mapping Options ---------------------------------------------------------

  input_list$map$dsem_beta <- factor(seq_along(input_list$par$dsem_beta)) # a fixed arrow (NA name) has no parameter here
  input_list$map$ln_dsem_sd <- factor(seq_along(input_list$par$ln_dsem_sd))

  map_mu <- rep(NA_integer_, n_var) # a linked series' mean stays at zero
  map_mu[match(mu_est, variables)] <- seq_along(mu_est)
  input_list$map$dsem_mu <- factor(map_mu)

  # a cell is estimated here only when nothing else sets it: not a known covariate, not a linked deviation.
  # rows before a linked series starts sit at its mean, so the first live year starts the way dsem starts a series
  map_x <- array(1:(n_grid_yrs * n_var), dim = c(n_grid_yrs, n_var))
  for(k in seq_len(n_cov)) if(family_code[k] == 0) map_x[!is.na(cov_obs[,k]), match(cov_names[k], variables)] <- NA
  for(s in seq_along(link$name)) map_x[c(seq_len(min(link$grid_row[[s]]) - 1), link$grid_row[[s]]), link_col[s]] <- NA
  input_list$map$dsem_x <- factor(map_x)

  # the cells the model is handed, which a linked recruitment cell's correction conditions on: a covariate's
  # observed years whatever its family, and the rows before a linked series starts, which sit at its mean
  x_known <- matrix(FALSE, n_grid_yrs, n_var)
  for(k in seq_len(n_cov)) x_known[!is.na(cov_obs[,k]), match(cov_names[k], variables)] <- TRUE
  for(s in seq_along(link$name)) x_known[seq_len(min(link$grid_row[[s]]) - 1), link_col[s]] <- TRUE
  input_list$data$dsem_x_known <- x_known

  input_list$map$ln_dsem_obs_sd <- factor(ifelse(family_code %in% c(1, 4, 6, 7), seq_len(n_cov), NA)) # normal, gamma, lognormal, tweedie
  input_list$map$logit_dsem_tweedie_p <- factor(ifelse(family_code == 7, seq_len(n_cov), NA))

  map_delta0 <- rep(NA_integer_, n_var)
  map_delta0[match(delta0_est, variables)] <- seq_along(delta0_est)
  input_list$map$dsem_delta0 <- factor(map_delta0)

  # Print Messages ----------------------------------------------------------

  collect_message("DSEM series: ", paste(variables, collapse = ", "))
  collect_message("DSEM grid: ", n_grid_yrs, " years, ", n_var, " series, ", nrow(arrows), " arrows")
  collect_message("Deviations linked (their own penalty comes off in sync_dev_map_data): ", paste(link$name, collapse = ", "))

  for(s in seq_along(link$name)) {
    short <- setdiff(1:n_grid_yrs, link$grid_row[[s]]) # grid years this series does not reach
    if(length(short) > 0) {
      collect_message(link$name[s], " reaches ", length(link$grid_row[[s]]), " of ", n_grid_yrs,
                      " grid years, so ", grid_years[min(short)], " onward is forecast by the dsem. ",
                      "Fix those rows in map$dsem_x if that is not wanted.")
    }
  } # end s loop

  # every estimated cell of an array linked leaves its sigma read by nothing through that penalty; the module's
  # own declaration fixes the sigma, so an estimated one here is a setup that was never told
  declaration_of <- c(ln_RecDevs = "RecDevs_model = 'dsem' in Setup_Mod_Rec",
                      ln_NAA = "NAA_re = 'dsem' in Setup_Mod_Biologicals",
                      ln_growth_devs = "growth_tv_model = 'dsem' in Setup_Mod_Biologicals",
                      ln_growth_semipar_devs = "growth_semipar = 'dsem' in Setup_Mod_Biologicals",
                      move_devs = "cont_vary_movement = 'dsem_...' in Setup_Mod_Movement",
                      ln_fish_q_devs = "fish_q_model = 'dsem' in Setup_Mod_Fishsel_and_Q",
                      ln_srv_q_devs = "srv_q_model = 'dsem' in Setup_Mod_Srvsel_and_Q")

  for(par_name in unique(link$par)) {

    sigma_par <- link$sigma_par[match(par_name, link$par)]
    map_levels <- dev_map_levels(input_list, par_name)
    all_linked <- length(setdiff(which(!is.na(map_levels)), unlist(link$cell[link$par == par_name]))) == 0
    sigma_estimated <- !is.na(sigma_par) && !is.null(input_list$map[[sigma_par]]) && any(!is.na(input_list$map[[sigma_par]]))

    if(all_linked && sigma_estimated) {
      stop("Every estimated ", par_name, " cell is linked, so nothing reads ", sigma_par, " through that ",
           "penalty any more, yet it is still estimated. Declare ", declaration_of[[par_name]],
           ", which fixes it, or fix it yourself.")
    }

  } # end par_name loop

  collect_message("DSEM means estimated: ", if(length(mu_est) == 0) "none" else paste(mu_est, collapse = ", "))
  if(length(delta0_est) > 0) collect_message("DSEM first year offsets estimated: ", paste(delta0_est, collapse = ", "))
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)

} # end function

#' Helper that defines the names of the deviation series a dsem can link
#'
#' A series is one time series of deviations, at fixed indices of every dim but the year:
#' recruitment on three regions is three series, one per region, named
#' \code{rec_Pop_1_Region_1} through \code{rec_Pop_1_Region_3}. An arrow line names the
#' series it acts on, and this writes those names from array indices rather than by hand.
#'
#' @param input_list List with \code{data}, \code{par} and \code{map}, after the setup
#'   function for that process has run.
#' @param process Deviations to name: \code{"rec"}, \code{"growth"},
#'   \code{"growth_semipar"}, \code{"NAA"}, \code{"move"}, \code{"fish_q"} or \code{"srv_q"}.
#' @param ... Indices to keep, named by dim in any case: \code{pop}, \code{region},
#'   \code{from}, \code{to}, \code{seas}, \code{age}, \code{sex}, \code{par} or
#'   \code{fleet}, whichever of those dims the process has. They are array indices and not
#'   labels, so \code{age = 3} is the third model age. A dim left out keeps every level.
#' @param estimated Whether to leave out series the map fixes, which an arrow has no
#'   parameter to link. \code{TRUE} by default.
#'
#' @return Character vector of series names, in the order \code{Setup_Mod_DSEM} reads them,
#'   and empty when the model has no such series, as movement on one region.
#'
#' @details Under CTMC movement a deviation sits on a region's preference rather than on a
#'   pair of regions, so \code{to} has one level and \code{from} is the region itself.
#'
#' @seealso \code{\link{Setup_Mod_DSEM}}
#'
#' @export dsem_series
#'
#' @examples
#' \dontrun{
#'   dsem_series(input_list, "rec")                          # every recruitment series
#'   s <- dsem_series(input_list, "move", from = 1, to = 2)   # one exchange, every age
#'   paste0(s, " <-> ", s, ", 0, sd_move")                    # arrow lines sharing one sd
#' }
dsem_series <- function(input_list,
                        process,
                        ...,
                        estimated = TRUE) {

  # Input Validation --------------------------------------------------------
  process_table <- dsem_process_table() # one entry per process, naming its array and that array's dims
  process_labels <- vapply(process_table, function(x) x$label, "")
  if(!is.character(process) || length(process) != 1 || !process %in% process_labels) {
    stop("process should be one of: ", paste(process_labels, collapse = ", "), ".")
  }
  entry <- process_table[[match(process, process_labels)]]
  name_dims <- setdiff(entry$dim_names, "Yr") # the dims a series name spells out, in the order it writes them
  filters <- list(...) # one vector of indices per dim being filtered
  if(length(filters) > 0 && (is.null(names(filters)) || any(names(filters) == ""))) {
    stop("Name every filter by its dim: ", paste(tolower(name_dims), collapse = ", "), ".")
  }
  filter_dim <- match(tolower(names(filters)), tolower(name_dims)) # the dim each filter picks on
  if(any(is.na(filter_dim))) {
    stop("'", paste(names(filters)[is.na(filter_dim)], collapse = "', '"), "' is not a dim of the ",
         process, " series. Dims: ", paste(tolower(name_dims), collapse = ", "), ".")
  }

  if(anyDuplicated(filter_dim)) stop("A dim is filtered twice. Give each dim one vector of indices.")
  dev_array <- input_list$par[[entry$par]] # the deviations this process estimates
  if(is.null(dev_array)) return(character(0)) # the model does not have the process at all
  n_levels <- dim(dev_array)[match(name_dims, entry$dim_names)] # levels each named dim has

  for(k in seq_along(filters)) {
    filter_idx <- filters[[k]] # the indices asked for on this dim
    if(!is.numeric(filter_idx) || any(filter_idx != round(filter_idx)) ||
       any(filter_idx < 1) || any(filter_idx > n_levels[filter_dim[k]])) {
      stop(names(filters)[k], " should be whole indices from 1 to ", n_levels[filter_dim[k]],
           ", the levels that dim has in ", entry$par, ".")
    }

  } # end k loop

  # Every Series the Array Has ----------------------------------------------
  n_grid_yrs <- length(input_list$data$years) + input_list$data$n_proj_yrs_devs # rows of the dsem grid
  all_series <- get_dsem_link(input_list, process, "", n_grid_yrs) # empty arrow text, so nothing is linked yet
  series_names <- all_series$offered # one name per series, in the order Setup_Mod_DSEM reads them
  series_index <- all_series$offered_idx # the array index each name sits at, in the array's own dim order
  series_cells <- all_series$offered_cell # where a series' years sit in the flattened array
  map_levels <- dev_map_levels(input_list, entry$par) # NA at a cell the map fixes

  # Keep the Series Asked For -----------------------------------------------
  keep <- rep(TRUE, length(series_names))
  for(s in seq_along(series_names)) {
    series_idx <- series_index[[s]][match(name_dims, entry$dim_names)] # this series' index per named dim
    for(k in seq_along(filters)) {
      if(!series_idx[filter_dim[k]] %in% filters[[k]]) keep[s] <- FALSE # a filter leaves this index out
    } # end k loop
    # nothing for an arrow to link when the map fixes the series in every year
    if(estimated && keep[s] && all(is.na(map_levels[series_cells[[s]]]))) keep[s] <- FALSE
  } # end s loop

  return(series_names[keep])

} # end function
