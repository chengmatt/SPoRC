# Stage 2 of 3: objective function
#
# Builds the precision matrix for correlated deviations over age, year and cohort. Used by the
# process error likelihoods in model_priors_penalties.R.

#' Construct a 3D sparse precision matrix over ages, years, and cohorts
#'
#' Builds a sparse \eqn{(n_{\text{ages}} \times n_{\text{yrs}}) \times
#' (n_{\text{ages}} \times n_{\text{yrs}})} precision matrix \eqn{Q} for a
#' Gaussian Markov random field (GMRF) with simultaneous autoregressive (SAR)
#' structure across three biological dimensions: age, year, and cohort
#' (age-year diagonal). The matrix is constructed via the path-matrix
#' factorization \eqn{Q = (I - B)^\top \Omega^{-1} (I - B)}, where \eqn{B}
#' encodes the partial correlations and \eqn{\Omega} is a diagonal variance
#' matrix. Two variance parameterizations are supported: marginal (stationary)
#' and conditional (non-stationary).
#'
#' @param n_ages Integer. Number of age classes.
#' @param n_yrs Integer. Number of years.
#' @param pcorr_age Numeric. Partial correlation along the age dimension
#'   (i.e., between adjacent ages within the same year).
#' @param pcorr_year Numeric. Partial correlation along the year dimension
#'   (i.e., between adjacent years within the same age).
#' @param pcorr_cohort Numeric. Partial correlation along the cohort diagonal
#'   (i.e., between the \eqn{(a-1, y-1)} and \eqn{(a, y)} cell).
#' @param ln_var_value Numeric. Log of the target variance. Exponentiated
#'   internally to \eqn{\sigma^2 = \exp(\text{ln\_var\_value})}.
#' @param Var_Type Integer. Variance parameterization: \code{0} = marginal
#'   (stationary) variance, where diagonal elements of \eqn{\Omega} are
#'   solved recursively via the accumulator \eqn{(I - B)^{-1}} to achieve a
#'   constant marginal variance \eqn{\sigma^2} at every node (slower);
#'   \code{1} = conditional (non-stationary) variance, where all diagonal
#'   elements of \eqn{\Omega} are set to \eqn{\sigma^2} directly (faster).
#'
#' @return A sparse \code{Matrix::sparseMatrix} precision matrix \eqn{Q} of
#'   dimension \eqn{(n_{\text{ages}} \times n_{\text{yrs}}) \times
#'   (n_{\text{ages}} \times n_{\text{yrs}})}, compatible with
#'   \code{RTMB::dgmrf}.
#'
#' @importFrom Matrix sparseMatrix
#' @importFrom methods as
#' @keywords internal
Get_3d_precision <- function(n_ages, n_yrs, pcorr_age, pcorr_year, pcorr_cohort, ln_var_value, Var_Type) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  index = expand.grid(seq_len(n_ages), seq_len(n_yrs)) # create index combinations to loop through
  i = j = x = numeric(0) # initialize posiiton to fill in precision matrix
  var_value = exp(ln_var_value) # transform to normal space

  for(n in seq_len(nrow(index))){
    age = index[n,1] # get age index out of all index combinations
    year = index[n,2] # get year index out of all index combinations
    if(age > 1) {
      i = c(i, n)
      j = c(j, which(index[,1] == (age - 1) & index[,2] == year))
      x = c(x, pcorr_age) # link to the age-adjacent neighbor in the same year
    }
    if(year > 1) {
      i = c(i, n)
      j = c(j, which(index[,1] == age & index[,2] == (year - 1)))
      x = c(x, pcorr_year) # link to the year-adjacent neighbor at the same age
    }
    if(age > 1 && year > 1) {
      i = c(i, n)
      j = c(j, which(index[,1] == (age - 1) & index[,2] == (year - 1)))
      x = c(x, pcorr_cohort) # cohort correlation indexing
    }
  } # end n loop

  # create B path matrix
  B = matrix(0, nrow = n_ages * n_yrs, ncol = n_ages * n_yrs)
  B[cbind(i, j)] = x
  B = as(B, "sparseMatrix")

  # identity matrix
  I = as(diag(1, n_ages * n_yrs, n_ages * n_yrs), "sparseMatrix")

  # Solve Omega recursively for stationary variance (accumulator function)
  if(Var_Type == 0) {
    L = solve(I - B) # solve to get accumulator function for stationary variance
    d = rep(0, nrow(index))
    for(n in seq_len(nrow(index))){
      if(n == 1) {
        d[n] = var_value
      }else{
        cumvar = sum(L[n,seq_len(n - 1)] * d[seq_len(n - 1)] * L[n,seq_len(n - 1)])
        d[n] = (var_value - cumvar) / L[n,n]^2
      }
    } # end n loop
  } # end marginal variance (stationary variance)

  if(Var_Type == 1) d = var_value # conditional variance (non-stationary variance)

  # omega matrix
  Omega_inv = diag(1 / d, n_ages * n_yrs, n_ages * n_yrs)
  Q = as((I - Matrix::t(B)) %*% Omega_inv %*% (I - B), "sparseMatrix") # solve for precision

  return(Q)
}

#' Values of the arrows from the parameter vectors
#'
#' Each arrow line becomes one number here. A path or covariance line reads its
#' coefficient from \code{dsem_beta}, an sd line reads \code{exp(ln_dsem_sd)}, a line
#' written with \code{NA} for a name keeps the value it was given, and a line named
#' after a series is left at zero because its value changes year to year and is read
#' off the grid where the matrices are filled.
#'
#' The lines \code{env -> rec, 0, b} (\code{dsem_beta[1] = 0.4}),
#' \code{env <-> env, 0, s} (\code{ln_dsem_sd[1] = log(0.9)}) and
#' \code{rec <-> rec, 0, NA, 1} give the values 0.4, 0.9 and 1.
#'
#' @param dsem_beta Numeric vector, one entry per estimated path or covariance.
#' @param ln_dsem_sd Numeric vector, one entry per estimated sd line, log scale.
#' @param dsem_model Output of \code{\link{read_dsem_arrows}}.
#'
#' @return Numeric vector with one value per arrow, in the arrows' order.
#'
#' @keywords internal
get_dsem_arrow_values = function(dsem_beta,
                                 ln_dsem_sd,
                                 dsem_model) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  arrows = dsem_model$arrows
  arrow_value = rep(0, nrow(arrows)) # container

  # fill in paths, covariances and sds, one branch each so a fixed arrow is not read again as estimated
  for(i in 1:nrow(arrows)) {
    if(arrows$mod_idx[i] > 0) next # moderated, so its value changes year to year
    else if(arrows$par[i] == 0) arrow_value[i] = arrows$start[i] # fixed arrow
    else if(arrows$type[i] == "sd") arrow_value[i] = exp(ln_dsem_sd[arrows$par[i]]) # sd, estimated on the log scale
    else arrow_value[i] = dsem_beta[arrows$par[i]] # path coefficient or covariance term
  } # end i loop

  return(arrow_value)

} # end function

#' Get DSEM matrices
#'
#' Stack the grid into one single vector, years within series. Every path arrow puts
#' its coefficient into \eqn{B} at the cell it points to (row) and the cell it
#' reads (column, the same series \code{lag} rows earlier), so
#' \eqn{(I - B)(x - \mu)} turns the grid into its innovations. Every sd line
#' puts its value on the diagonal of \eqn{\Gamma} and every covariance line
#' off it, so \eqn{V = \Gamma^{\top}\Gamma} is the covariance of the innovations
#' in one year. The positions are worked out once, outside of the tape, by
#' \code{\link{get_dsem_cells}}; this function only writes the numbers into the
#' stored slots, which is what lets the matrices be built on the tape.
#'
#' A moderated arrow has no single value: its coefficient in year \eqn{t} is
#' the moderating series' value in that year, read from \code{x_grid}, and a
#' moderated sd is that value or its exponential (\code{mod_var_logscale}).
#' Also, under \code{variance = "diagonal"} or \code{"marginal"} (set in
#' \code{\link{read_dsem_arrows}}) an sd line is the marginal sd of its series
#' rather than the innovation sd.
#'
#' @param dsem_beta,ln_dsem_sd,dsem_model As in \code{\link{get_dsem_arrow_values}}.
#' @param dsem_cells Output of \code{\link{get_dsem_cells}}.
#' @param x_grid Matrix \code{[year, series]} of the grid, needed when an arrow
#'   is moderated.
#' @param need_Vinv Whether to form \eqn{V^{-1}}; the derived-series route only
#'   needs the sds.
#'
#' @return List with \code{IminusB} (sparse, \eqn{I - B}), \code{Vinv} (sparse
#'   \eqn{V^{-1}}, or \code{NULL} when not asked for) and \code{sd_cell} (each
#'   cell's innovation sd, or \code{NULL} when a covariance line couples them).
#'
#' @keywords internal
get_dsem_matrices <- function(dsem_beta,
                              ln_dsem_sd,
                              dsem_model,
                              dsem_cells,
                              x_grid = NULL,
                              need_Vinv = TRUE) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  arrows = dsem_model$arrows
  arrow_value = get_dsem_arrow_values(dsem_beta, ln_dsem_sd, dsem_model)

  if(dsem_cells$has_mod && is.null(x_grid)) {
    stop("An arrow is named after a series, so its value comes from x_grid. Pass x_grid to get_dsem_precision.")
  }

  # get I - B (has the path coeffs)
  IminusB = RTMB::AD(dsem_cells$IminusB$m)
  entry_value = c(1, -arrow_value)[dsem_cells$IminusB$entry_arrow + 1] # entry arrow 0 is the diagonal

  # a moderated path takes the moderating series' value in the year the arrow points to
  path_mod = which(c(0L, arrows$mod_idx)[dsem_cells$IminusB$entry_arrow + 1] > 0)

  if(length(path_mod) > 0) {
    mod_series = arrows$mod_idx[dsem_cells$IminusB$entry_arrow[path_mod]] # the series each one is moderated by
    path_mod_cell = (mod_series - 1) * dsem_cells$n_grid_yrs + dsem_cells$IminusB$entry_yr[path_mod] # figure out moderated path cells
    entry_value[path_mod] = -x_grid[path_mod_cell] # add moderated path stuff into entry values
  } # end if any moderated path

  IminusB@x = RTMB::AD(entry_value[dsem_cells$IminusB$slot_entry]) # input path coefficients into I - B

  # the sd and covariance lines fill Gamma
  gamma_value = arrow_value[dsem_cells$Gamma$entry_arrow]
  gamma_mod = which(arrows$mod_idx[dsem_cells$Gamma$entry_arrow] > 0)

  if(length(gamma_mod) > 0) {
    mod_series = arrows$mod_idx[dsem_cells$Gamma$entry_arrow[gamma_mod]]
    gamma_mod_cell = (mod_series - 1) * dsem_cells$n_grid_yrs + dsem_cells$Gamma$entry_yr[gamma_mod]
    gamma_value[gamma_mod] = if(dsem_model$mod_var_logscale) exp(x_grid[gamma_mod_cell]) else x_grid[gamma_mod_cell]
  } # end if any moderated sd or covariance

  # the diagonal / marginal form: each sd line is the marginal sd of its cell
  if(isTRUE(dsem_model$variance %in% c("diagonal", "marginal"))) {

    # setup stuff
    n_cells = dsem_cells$n_cells
    IminusB_dense = matrix(0, n_cells, n_cells)
    IminusB_dense[(dsem_cells$IminusB$entry_col - 1) * n_cells + dsem_cells$IminusB$entry_row] = entry_value
    A = solve(IminusB_dense) # response of every cell to every innovation - solve dense and add back into gamma

    is_sd = arrows$type[dsem_cells$Gamma$entry_arrow] == "sd"
    sd_row = dsem_cells$Gamma$entry_row[is_sd]
    target_var = rep(0, n_cells)
    target_var[sd_row] = gamma_value[is_sd]^2 # the sd line squared, the variance each cell should end up with

    if(!dsem_cells$has_cov || dsem_model$variance == "diagonal") {

      # with independent innovations each cell's variance is (A o A) times them, so invert that
      innovation_var = solve(A*A, target_var)
      gamma_value[is_sd] = sqrt(innovation_var[sd_row])

    } else {

      # marginal version
      Gamma0 = matrix(0, n_cells, n_cells)
      Gamma0[(dsem_cells$Gamma$entry_col - 1) * n_cells + dsem_cells$Gamma$entry_row] = gamma_value
      V0 = t(Gamma0) %*% Gamma0 # the innovation covariance the lines imply, before any scaling
      v0 = diag(V0)
      scale2 = solve(A*A, target_var) / v0 # start with no cross terms

      # each pass scales the innovations by d, reads off what that gives every cell, and takes the
      # cross terms out of the target so the next solve only has to cover a cell's own share
      for(iter in 1:15) {
        W = A * rep(sqrt(scale2), each = n_cells) # A D, column k scaled by d_k
        total_var = rowSums((W %*% V0) * W) # diag(A D V0 D t(A)), every cell's variance at this d
        own_var = as.vector((A*A) %*% (scale2 * v0)) # the part without cross terms
        scale2 = solve(A*A, target_var - (total_var - own_var)) / v0
      } # end iter loop

      gamma_value = gamma_value * sqrt(scale2)[dsem_cells$Gamma$entry_col]

    } # end if covariance lines under the marginal form
  } # end if a marginal variance form

  if(!dsem_cells$has_cov) { # no covariance

    sd_cell = rep(0, dsem_cells$n_cells)
    sd_cell[dsem_cells$Gamma$entry_row] = gamma_value
    Vinv = NULL

    if(need_Vinv) {
      Vinv = RTMB::AD(Matrix::sparseMatrix(i = 1:dsem_cells$n_cells, j = 1:dsem_cells$n_cells, x = 1))
      Vinv@x = RTMB::AD(1 / sd_cell^2) # unit diagonal
    }

  } else {

    # covariance terms -> V = t(Gamma) Gamma
    Gamma = RTMB::AD(dsem_cells$Gamma$m)
    Gamma@x = RTMB::AD(gamma_value[dsem_cells$Gamma$slot_entry])
    Vinv = solve(Matrix::t(Gamma) %*% Gamma) # V = t(Gamma) Gamma, as dsem reads a covariance line
    sd_cell = NULL # the innovations are not independent, so leave sd as NULL

  } # end if covariance arrows

  return(list(IminusB = IminusB, Vinv = Vinv, sd_cell = sd_cell))

} # end function

#' The precision of the whole grid
#'
#' \eqn{Q = (I - B)^{\top} V^{-1} (I - B)}, the inverse covariance of the
#' stacked grid: if the innovations have covariance \eqn{V} and the grid is
#' \eqn{x - \mu = (I - B)^{-1}\varepsilon}, then the grid's covariance is
#' \eqn{(I - B)^{-1} V (I - B)^{-\top}} and its inverse is this product. It is
#' the same form \code{\link{Get_3d_precision}} builds for the numbers at age
#' field. \code{RTMB::dgmrf} evaluates the density from it without forming the
#' covariance.
#'
#' @inheritParams get_dsem_matrices
#'
#' @return Sparse precision matrix over the \code{n_grid_yrs * n_series} cells.
#'
#' @keywords internal
get_dsem_precision <- function(dsem_beta,
                               ln_dsem_sd,
                               dsem_model,
                               dsem_cells,
                               x_grid = NULL) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  parts <- get_dsem_matrices(dsem_beta, ln_dsem_sd, dsem_model, dsem_cells, x_grid)
  Q <- Matrix::t(parts$IminusB) %*% parts$Vinv %*% parts$IminusB

  return(Q)

} # end function

#' The variance of each grid cell given the cells the model is handed
#'
#' A lognormal deviation with variance \eqn{v} has \eqn{E[\exp(x)] = \exp(\mu + v/2)},
#' so recruitment keeps its mean at \eqn{R_0} only if the cell's mean drops by
#' \eqn{v/2}. The \eqn{v} that does it is the cell's variance given what the model
#' is handed (the observed covariate values, and the rows before a series starts),
#' not the sd line's square: under a self path \eqn{\rho} a settled year has
#' \eqn{\sigma^2 / (1 - \rho^2)}, and every lagged path adds to it. For unknown
#' cells \eqn{U} and known cells \eqn{K}, \eqn{\mathrm{Var}(x_U \mid x_K) = (Q_{UU})^{-1}},
#' the inverse of the unknown block of the precision, so each cell's variance is a
#' diagonal entry of that inverse.
#'
#' A series with an sd of zero (derived) has no innovation, so its rows leave the
#' quadratic form; an unknown cell pointing into one is refused, since the derived
#' value would then be random while the model treats it as fixed.
#'
#' @inheritParams get_dsem_matrices
#' @param known_cell Logical over the \code{n_grid_yrs * n_series} cells (years
#'   within series), \code{TRUE} where the value is handed to the model.
#'
#' @return Matrix \code{[year, series]} of variances, zero on the known cells.
#'
#' @keywords internal
get_dsem_margvar = function(dsem_beta,
                            ln_dsem_sd,
                            x_grid,
                            dsem_model,
                            dsem_cells,
                            known_cell) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  n_grid_yrs = nrow(x_grid)
  derived_cell = rep(dsem_model$derived, each = n_grid_yrs)
  any_derived = any(derived_cell)
  parts = get_dsem_matrices(dsem_beta, ln_dsem_sd, dsem_model, dsem_cells, x_grid, need_Vinv = !any_derived)

  unknown = which(!known_cell & !derived_cell) # the cells with a variance to work out to make stationary
  stochastic = which(!derived_cell) # the rows with an innovation of their own
  if(length(unknown) == 0) return(matrix(0, n_grid_yrs, ncol(x_grid)))

  # a derived series has no innovation, so its rows leave the quadratic form and the
  # remaining ones get a diagonal V of their own sds
  if(any_derived) {

    if(is.null(parts$sd_cell)) {
      stop("A series with an sd of zero cannot be given together with a covariance arrow: ",
           "its innovation is zero, so the covariance of the innovations is singular.")
    }

    feeds = parts$IminusB[which(derived_cell), unknown, drop = FALSE] # arrows from an unknown cell into a derived series

    if(length(feeds@x) > 0) {
      stop("A series with an sd of zero reads a cell whose value is estimated, so its own values ",
           "are not fixed and the recruitment correction cannot condition on them. Give the series ",
           "an sd, or observe every cell it reads.")
    }

    Vinv_s = RTMB::AD(Matrix::sparseMatrix(i = seq_along(stochastic), j = seq_along(stochastic), x = 1))
    Vinv_s@x = RTMB::AD(1 / parts$sd_cell[stochastic]^2)

  } else Vinv_s = parts$Vinv

  IminusB_su = parts$IminusB[stochastic, unknown, drop = FALSE]
  Q_uu = Matrix::t(IminusB_su) %*% Vinv_s %*% IminusB_su # solve for precision of the unknown cells given the known cells / derived
  var_cell = rep(0, length(known_cell))
  var_cell[unknown] = diag(solve(Q_uu))

  return(matrix(var_cell, n_grid_yrs, ncol(x_grid)))

} # end function

#' Work out the time-series where the sd is fixed at zero
#'
#' A series with an sd of zero has no innovation, so its value each year is
#' exactly what the arrows pointing into it give: its own mean plus, for every
#' path into it, that arrow's coefficient times how far the series it reads
#' sits from its mean in the year it reads. Series are visited in an order that
#' sets each one before anything reads it (\code{series_order}), one year at a
#' time.
#'
#' @param x_grid Matrix \code{[year, series]}, the derived columns overwritten.
#' @param only Integer vector of series to fill, or \code{NULL} (default) for
#'   every derived series. Used to fill the catchability series before the
#'   observations are built, since the rest need the population first.
#' @param mu_grid Matrix \code{[year, series]} of series means.
#' @param arrow_value Output of \code{\link{get_dsem_arrow_values}}.
#' @param dsem_model Output of \code{\link{read_dsem_arrows}}.
#'
#' @return \code{x_grid} with every derived cell set.
#'
#' @keywords internal
fill_dsem_derived = function(x_grid,
                             mu_grid,
                             arrow_value,
                             dsem_model,
                             only = NULL) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  arrows = dsem_model$arrows
  if(!any(dsem_model$derived)) return(x_grid) # every series has an innovation, so nothing to figure out

  for(yr in 1:nrow(x_grid)) {
    for(j in dsem_model$series_order) { # an order that sets each series before anything reading it

      if(!dsem_model$derived[j]) next # estimated, so its value is already in x_grid
      if(!is.null(only) && !j %in% only) next # this call fills a named subset of the derived series
      value = mu_grid[yr,j] # this series' own mean, before any arrow

      for(i in which(arrows$type == "path" & arrows$to_idx == j)) {

        from_yr = yr - arrows$lag[i] # year the arrow reads, this same year at lag 0
        from_j = arrows$from_idx[i] # series the arrow reads
        if(from_yr < 1) next # before the first year, so this arrow adds nothing

        # the slope is the arrow's own number, or the moderating series' value this year
        slope = if(arrows$mod_idx[i] > 0) x_grid[yr,arrows$mod_idx[i]] else arrow_value[i]
        deviation = x_grid[from_yr,from_j] - mu_grid[from_yr,from_j] # how far that series sits off its mean
        value = value + slope * deviation

      } # end i loop

      x_grid[yr,j] = value

    } # end j loop
  } # end yr loop

  return(x_grid)

} # end function

#' Negative log density of the dsem grid
#'
#' The grid is a multivariate normal with mean \eqn{\mu} and precision
#' \eqn{Q = (I - B)^{\top} V^{-1} (I - B)}:
#' \deqn{-\ell = -\tfrac{1}{2}\log|Q| + \tfrac{1}{2}(x - \mu)^{\top} Q (x - \mu) + \tfrac{n}{2}\log 2\pi,}
#' evaluated by \code{RTMB::dgmrf}.
#'
#' A first-year offset \eqn{\delta_0} shifts year one of each series and is
#' propagated through the paths by \eqn{(I - B)^{-1}}, which is how an AR1
#' started off its mean decays back at \eqn{\rho^{t-1}}.
#'
#' A series whose sd is fixed at zero makes \eqn{V} singular, so \eqn{Q} does
#' not exist. Those cells are worked out by \code{\link{fill_dsem_derived}}
#' and dropped, and what is left is the sum of normal densities of the
#' remaining innovations.
#'
#' @param dsem_beta,ln_dsem_sd,dsem_model,dsem_cells,x_grid As in
#'   \code{\link{get_dsem_matrices}}.
#' @param mu_grid Matrix \code{[year, series]} of series means: a covariate's
#'   mean, zero for a deviation series, and the log deterministic prediction for
#'   a numbers at age series. The density reads \code{x_grid - mu_grid}, so for
#'   numbers at age that difference is the innovation, which is stored nowhere.
#' @param delta0 Optional numeric vector, one first-year offset per series.
#'
#' @return Scalar negative log density.
#'
#' @keywords internal
get_dsem_nLL = function(dsem_beta,
                        ln_dsem_sd,
                        x_grid,
                        mu_grid,
                        dsem_model,
                        dsem_cells,
                        delta0 = NULL) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  n_grid_yrs = nrow(x_grid)
  any_derived = any(dsem_model$derived)
  arrow_value = get_dsem_arrow_values(dsem_beta, ln_dsem_sd, dsem_model) # get arrow values

  # if we have derived vars
  if(any_derived) {

    if(!is.null(delta0)) {
      stop("A first year offset and a series with an sd of zero cannot be given together: the offset ",
           "would have to be worked out from the same paths that set the derived series.")
    }

    x_grid = fill_dsem_derived(x_grid, mu_grid, arrow_value, dsem_model) # an sd of zero means the value is not estimated

  } # end if any derived series

  # get parts for constructing Q
  parts <- get_dsem_matrices(dsem_beta, ln_dsem_sd, dsem_model, dsem_cells, x_grid, need_Vinv = !any_derived)

  # an offset in the first year of each series (via delta 0), propagated through all the paths
  if(!is.null(delta0)) {
    delta_cell <- rep(0, dsem_cells$n_cells)
    delta_cell[(seq_along(delta0) - 1) * n_grid_yrs + 1] <- delta0
    mu_grid <- mu_grid + matrix(as.vector(solve(parts$IminusB, delta_cell)), n_grid_yrs, ncol(x_grid))
  }

  if(!any_derived) {

    Q = Matrix::t(parts$IminusB) %*% parts$Vinv %*% parts$IminusB
    dsem_nLL = -RTMB::dgmrf(as.vector(x_grid), mu = as.vector(mu_grid), Q = Q, log = TRUE)

  } else {

    # checks; an sd of zero makes V singular, so Q does not exist and dgmrf cannot be used. dropping those cells leaves a unit determinant
    if(is.null(parts$sd_cell)) {
      stop("A series with an sd of zero cannot be given together with a covariance arrow: ",
           "its innovation is zero, so the covariance of the innovations is singular.")
    }

    if(!dsem_cells$det_is_one) {
      stop("A series with an sd of zero cannot be given together with a lead or a same-year loop: ",
           "its value is worked out one year at a time, so every arrow into it has to be set first.")
    }

    # fill_dsem_derived set those cells to exactly what their arrows predict, so their innovation is
    # zero as an identity rather than a residual, and a density on it would be meaningless
    has_density = !rep(dsem_model$derived, each = n_grid_yrs)
    innovation = parts$IminusB %*% as.vector(x_grid - mu_grid)
    dsem_nLL = -sum(RTMB::dnorm(as.vector(innovation)[has_density], 0, parts$sd_cell[has_density], TRUE))

  } # end if any series has an sd of zero

  return(dsem_nLL)

} # end function

#' Covariate observation density by family and link
#'
#' The observations of one covariate given its grid cells. The cell goes through the
#' link to the mean (identity, exp, inverse logit or inverse cloglog), and the
#' family's density is taken about that mean. Fixed (0) has no density. Normal (1) has
#' an estimated sd, gaussian_fixed_sd (5) a known sd per observation, bernoulli (2) a
#' coin flip at the mean, poisson (3) that mean, Gamma (4) shape \eqn{1/CV^2} and that
#' mean, lognormal (6) that mean as its median, tweedie (7) that mean with a
#' dispersion and a power in (1, 2).
#'
#' @param y Observed values, no NA.
#' @param x Grid cells for those years, on the link scale.
#' @param family Family code.
#' @param link Link code.
#' @param obs_sd Measurement sd (normal), CV (Gamma), sd of the log
#'   (lognormal) or dispersion (tweedie). Unused otherwise.
#' @param tweedie_p Tweedie power in (1, 2). Unused otherwise.
#' @param fixed_sd Known sd per observation for gaussian_fixed_sd. Unused
#'   otherwise.
#'
#' @return Scalar negative log likelihood.
#'
#' @keywords internal
get_dsem_obs_nLL = function(y,
                            x,
                            family,
                            link,
                            obs_sd,
                            tweedie_p,
                            fixed_sd = NULL) {

  # the mean of the observation, the cell through the link
  mu = if(link == 1) exp(x) # log
       else if(link == 2) 1 / (1 + exp(-x)) # logit
       else if(link == 3) 1 - exp(-exp(x)) # cloglog
       else x # identity

  if(family == 1) return(-sum(RTMB::dnorm(y, mu, obs_sd, TRUE))) # normal
  if(family == 2) return(-sum(RTMB::dbinom(y, 1, mu, TRUE))) # bernoulli
  if(family == 3) return(-sum(RTMB::dpois(y, mu, TRUE))) # poisson
  if(family == 4) return(-sum(RTMB::dgamma(y, shape = 1 / obs_sd^2, scale = mu * obs_sd^2, log = TRUE))) # gamma, obs_sd the CV
  if(family == 5) return(-sum(RTMB::dnorm(y, mu, fixed_sd, TRUE))) # normal with a known sd per observation
  if(family == 6) return(-sum(RTMB::dlnorm(y, log(mu), obs_sd, TRUE))) # lognormal
  if(family == 7) return(-sum(RTMB::dtweedie(y, mu, obs_sd, tweedie_p, TRUE))) # tweedie

  return(0) # fixed values have no density

} # end function
