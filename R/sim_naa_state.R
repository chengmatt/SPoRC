# State Setup ---------------------------------------------------------------

#' Specify the state-space numbers-at-age process for simulation
#'
#' Turns on process error in the numbers at age for a simulated population. The
#' operating model applies the same centered state the estimation model does:
#' the deterministic mortality and ageing step is computed, then the numbers are
#' multiplied by \eqn{\exp(\eta)} with \eqn{\eta} drawn from the covariance the
#' arguments here describe.
#'
#' Arguments mirror \code{Setup_Mod_Biologicals}'s state-space options so a
#' simulated population and the model fitted to it are written the same way, so
#' a self test compares like with like rather than one translation against
#' another.
#'
#' @param sim_list Simulation list from \code{\link{Setup_Sim_Dim}}.
#' @param NAA_re Character. \code{"none"} (default) leaves the numbers at age
#'   deterministic past recruitment and the initial age structure. Otherwise one
#'   of \code{"iid"}, \code{"1dar1_a"}, \code{"1dar1_y"}, \code{"2dar1"},
#'   \code{"3dcond"} or \code{"3dmarg"}.
#' @param sigmaNAA Numeric. Conditional standard deviation of the innovations,
#'   the same quantity \code{ln_sigmaNAA} holds in the estimation model. Under an
#'   autoregressive form the marginal standard deviation is larger by
#'   \eqn{1/\sqrt{1 - \rho^2}} per correlated dim.
#' @param rho_age,rho_year,rho_cohort Numeric correlations in \eqn{(-1, 1)} over
#'   the age, year and cohort dims. Only the ones the chosen form reads are
#'   used.
#' @param NAA_re_pop,NAA_re_region,NAA_re_sex,NAA_re_season Character,
#'   \code{"iid"} (default) or \code{"us"}, an unstructured correlation across
#'   that dim.
#' @param pop_corr,region_corr,sex_corr,season_corr Numeric vectors of length
#'   \eqn{n(n-1)/2} giving the correlations for those dims, ordered as the
#'   strict lower triangle is filled by column. A single value is recycled.
#' @param NAA_re_ages,NAA_re_years Ages and year indices the state covers.
#'   \code{NULL} (default) uses everything from the second onward.
#' @param naa_eta_input Array
#'   \code{[n_pop, n_regions, n_cond_yrs, n_seas, n_ages, n_sexes, n_sims]} of
#'   log-scale state innovations, or \code{NULL} (default) to draw every year.
#'   Supplying a fit's own innovations, \eqn{\ln N - \ln \hat{N}}, makes the
#'   operating model reproduce that fit's numbers at age rather than a fresh
#'   realization of the same process. The year extent says how many leading years
#'   are held: give the fitted years and any year beyond them is still drawn, as
#'   a closed loop running past the data
#' needs.
#' @param NAA_re_seasons Seasons the state covers. \code{"annual"} (default) puts
#'   a state at season one only, leaving the numbers deterministic between
#'   seasons; \code{"all"} puts one at the start of every season, and an integer
#'   vector selects specific seasons.
#'
#' @return \code{sim_list} with the state-space settings attached.
#'
#' @export Setup_Sim_NAA_state
#' @family Simulation Setup
Setup_Sim_NAA_state <- function(sim_list,
                                NAA_re = "none",
                                sigmaNAA = 0.3,
                                rho_age = 0,
                                rho_year = 0,
                                rho_cohort = 0,
                                NAA_re_pop = "iid",
                                NAA_re_region = "iid",
                                NAA_re_sex = "iid",
                                NAA_re_season = "iid",
                                pop_corr = 0,
                                region_corr = 0,
                                sex_corr = 0,
                                season_corr = 0,
                                NAA_re_ages = NULL,
                                NAA_re_years = NULL,
                                NAA_re_seasons = "annual",
                                naa_eta_input = NULL) {

  codes <- c(none = 0, iid = 1, `1dar1_a` = 2, `1dar1_y` = 3, `2dar1` = 4, `3dcond` = 5, `3dmarg` = 6)
  if(length(NAA_re) != 1 || !NAA_re %in% names(codes))
    stop("NAA_re is '", NAA_re, "'. Valid options: ", paste(unique(names(codes)), collapse = ", "))

  margin_codes <- c(iid = 0, us = 1)
  for(name in c("NAA_re_pop", "NAA_re_region", "NAA_re_sex", "NAA_re_season")) {
    v <- get(name)
    if(length(v) != 1 || !v %in% names(margin_codes))
      stop(name, " is '", v, "'. Valid options: iid, us")
  } # end name loop

  n_ages <- sim_list$n_ages
  n_yrs <- sim_list$n_yrs
  n_seas <- if(is.null(sim_list$n_seas)) 1 else sim_list$n_seas

  sim_list$NAA_re <- codes[[NAA_re]]
  sim_list$sigmaNAA <- sigmaNAA
  sim_list$naa_rho <- c(age = rho_age, year = rho_year, cohort = rho_cohort)
  sim_list$NAA_re_pop <- margin_codes[[NAA_re_pop]]
  sim_list$NAA_re_region <- margin_codes[[NAA_re_region]]
  sim_list$NAA_re_sex <- margin_codes[[NAA_re_sex]]
  sim_list$NAA_re_season <- margin_codes[[NAA_re_season]]
  sim_list$naa_pop_corr <- pop_corr
  sim_list$naa_region_corr <- region_corr
  sim_list$naa_sex_corr <- sex_corr
  sim_list$naa_season_corr <- season_corr

  # the state covers ages two and older in years two onward, as it does in the estimation model:
  # age one is recruitment and year one at older ages is the initial age structure
  sim_list$naa_re_ages <- if(is.null(NAA_re_ages)) 2:n_ages else NAA_re_ages
  sim_list$naa_re_yrs <- if(is.null(NAA_re_years)) 2:n_yrs else NAA_re_years

  # season one is the year boundary, so it alone reproduces the annual state
  sim_list$naa_re_seas <- if(identical(NAA_re_seasons, "annual")) 1 else
                          if(identical(NAA_re_seasons, "all")) seq_len(n_seas) else
                          sort(unique(as.integer(NAA_re_seasons)))
  if(!all(sim_list$naa_re_seas %in% seq_len(n_seas)))
    stop("NAA_re_seasons is read as season indices into 1:", n_seas, ", or the strings ",
         "\"annual\" and \"all\". It was: ", paste(NAA_re_seasons, collapse = ", "))

  # innovations supplied rather than drawn. however many years are given are kept and the rest of the
  # run still draws, which is what a closed loop past the data needs
  if(!is.null(naa_eta_input)) {
    d <- as.integer(dim(naa_eta_input))
    want <- c(sim_list$n_pop, sim_list$n_regions, NA, n_seas, n_ages, sim_list$n_sexes, sim_list$n_sims)
    if(length(d) != 7 || !identical(d[-3], as.integer(want[-3])))
      stop("naa_eta_input should be [", paste(replace(want, 3, "n_cond_yrs"), collapse = ", "), "], population, ",
           "region, the years to hold, season, age, sex and replicate. It was [", paste(d, collapse = ", "), "].")
    if(d[3] > n_yrs)
      stop("naa_eta_input holds ", d[3], " years but the operating model runs ", n_yrs, ".")
    sim_list$naa_eta_input <- naa_eta_input
  }

  sim_list
}


#' The fit's numbers-at-age process, as the operating model draws it
#'
#' The sd by cell, read through the fit's sigma blocks; one correlation per dim
#' of the age-year field, a per-cell parameter averaged before it is transformed;
#' and the unstructured correlations across populations, regions, sexes and
#' seasons in the lower-triangle order \code{draw_naa_innovations} reads.
#' \code{simulation_self_test} draws the years after \code{n_cond_yrs} from it
#' and \code{condition_closed_loop_simulations} the projection years.
#'
#' @param data Data list of the fit.
#' @param pars Parameter list at the fitted values.
#'
#' @return List of \code{sigmaNAA} \code{[n_pop, n_regions, n_yrs, n_seas,
#'   n_ages, n_sexes]} over the fitted years, \code{naa_rho}, the four margin
#'   switches \code{NAA_re_pop}, \code{NAA_re_region}, \code{NAA_re_sex} and
#'   \code{NAA_re_season}, and their correlations \code{naa_pop_corr},
#'   \code{naa_region_corr}, \code{naa_sex_corr} and \code{naa_season_corr}.
#'
#' @keywords internal
naa_process_from_fit <- function(data, pars) {

  # the sd each cell's innovation is drawn at, through the fit's blocks
  blocks <- data$naa_sigma_blocks[,,seq_along(data$years),,,,drop = FALSE]
  sigma <- array(exp(pars$ln_sigmaNAA)[as.vector(blocks)], dim = dim(blocks))

  # one correlation per dim of the age-year field, a per-cell parameter averaged first
  pe <- pars$NAA_pe_pars
  naa_rho <- if(is.null(pe)) c(age = 0, year = 0, cohort = 0)
             else c(age = rho_trans(mean(pe[,,1,])), year = rho_trans(mean(pe[,,2,])), cohort = rho_trans(mean(pe[,,3,])))

  # unstructured correlations, in the lower-triangle order the draw reads them in
  corr_from <- function(corr_pars, n) {
    if(is.null(corr_pars) || n < 2) return(0)
    C <- build_us_corr(as.vector(corr_pars), n)
    C[lower.tri(C)]
  }
  margin_on <- function(opt_name) if(is.null(data[[opt_name]])) 0 else data[[opt_name]]
  n_seas_re <- if(is.null(data$naa_re_seas)) 1 else length(data$naa_re_seas)

  list(sigmaNAA = sigma,
       naa_rho = naa_rho,
       NAA_re_pop = margin_on("NAA_re_pop"),
       NAA_re_region = margin_on("NAA_re_region"),
       NAA_re_sex = margin_on("NAA_re_sex"),
       NAA_re_season = margin_on("NAA_re_season"),
       naa_pop_corr = if(margin_on("NAA_re_pop") == 1) corr_from(pars$NAA_pop_corr_pars, data$n_pop) else 0,
       naa_region_corr = if(margin_on("NAA_re_region") == 1) corr_from(pars$NAA_region_corr_pars, data$n_regions) else 0,
       naa_sex_corr = if(margin_on("NAA_re_sex") == 1) corr_from(pars$NAA_sex_corr_pars, data$n_sexes) else 0,
       naa_season_corr = if(margin_on("NAA_re_season") == 1) corr_from(pars$NAA_season_corr_pars, n_seas_re) else 0)

} # end naa_process_from_fit

# Innovation Draws ----------------------------------------------------------

#' Apply a correlation factor along one dim of an array
#'
#' Colors an array of independent normals so that one dim has a given
#' correlation, the reverse of the whitening the penalty uses. Operating on plain
#' doubles rather than on the AD tape, so the dim can simply be permuted to the
#' front rather than being reached by index arithmetic.
#'
#' @param x Numeric array.
#' @param L Lower triangular factor of the dim's correlation matrix.
#' @param dim_idx Integer dim to apply it along.
#'
#' @return An array of the same shape.
#'
#' @keywords internal
color_naa_dim <- function(x, L, dim_idx) {
  d <- dim(x)
  perm <- c(dim_idx, setdiff(seq_along(d), dim_idx))
  xp <- aperm(x, perm)
  z <- L %*% matrix(xp, nrow = d[dim_idx])
  aperm(array(z, dim = d[perm]), order(perm))
}


#' Draw state-space numbers-at-age innovations
#'
#' Draws \eqn{\eta} for a whole replicate at once. Drawing year by year would only
#' work when the year dim is independent or Markov; a separable autoregression
#' or a three-dimensional field correlates the whole span, so the array is built
#' up front and applied as the year loop reaches each boundary.
#'
#' Correlation is imposed dim by dim on independent normals, applying each
#' dim's Cholesky factor in turn, which is the reverse of how the penalty
#' whitens them. The three-dimensional field is the exception: its cohort term
#' couples age and year, so those two dims are drawn together from the sparse
#' precision rather than separately.
#'
#' @param sim_env Simulation environment with the settings from
#'   \code{\link{Setup_Sim_NAA_state}} and the dimensions.
#'
#' @return Array \code{[pop, region, year, season, age, sex]} of innovations, zero
#'   outside the active ages, years and seasons.
#'
#' @keywords internal
draw_naa_innovations <- function(sim_env) {

  np <- sim_env$n_pop
  nr <- sim_env$n_regions
  ns <- sim_env$n_sexes
  ny <- sim_env$n_yrs
  na <- sim_env$n_ages
  nk <- sim_env$n_seas
  code <- sim_env$NAA_re
  rho <- sim_env$naa_rho

  eta <- array(stats::rnorm(np * nr * ny * nk * na * ns), dim = c(np, nr, ny, nk, na, ns))

  ar1_chol <- function(n, r) if(n == 1 || r == 0) diag(n) else t(chol(r^abs(outer(1:n, 1:n, "-"))))

  if(code %in% c(5, 6)) {
    # the cohort term couples age and year, so those dims come from the joint precision
    Q <- Get_3d_precision(na, ny, rho[["age"]], rho[["year"]], rho[["cohort"]], 0,
                          Var_Type = if(code == 5) 1 else 0)
    Lq <- Matrix::Cholesky(methods::as(Matrix::forceSymmetric(Q), "sparseMatrix"), LDL = FALSE)
    for(p in 1:np) for(r in 1:nr) for(k in 1:nk) for(s in 1:ns) {
      # Get_3d_precision numbers its nodes age fastest, matching a [age, year] layout
      z <- as.vector(t(array(eta[p,r,,k,,s], dim = c(ny, na))))
      eta[p,r,,k,,s] <- t(array(as.vector(Matrix::solve(Lq, z, system = "Lt")), dim = c(na, ny)))
    } # end p, r, k, s loop
  } else {
    if(code %in% c(3, 4)) eta <- color_naa_dim(eta, ar1_chol(ny, rho[["year"]]), 3)
    if(code %in% c(2, 4)) eta <- color_naa_dim(eta, ar1_chol(na, rho[["age"]]), 5)
  }

  # correlation across the remaining dims, each an unstructured factor
  us_chol <- function(cor_vals, n) {
    if(n == 1) return(diag(1))
    v <- rep(cor_vals, length.out = n * (n - 1) / 2)
    C <- diag(n)
    C[lower.tri(C)] <- v
    C[upper.tri(C)] <- t(C)[upper.tri(C)]
    t(chol(C))
  }
  if(isTRUE(sim_env$NAA_re_pop == 1)) eta <- color_naa_dim(eta, us_chol(sim_env$naa_pop_corr, np), 1)
  if(isTRUE(sim_env$NAA_re_region == 1)) eta <- color_naa_dim(eta, us_chol(sim_env$naa_region_corr, nr), 2)
  if(isTRUE(sim_env$NAA_re_sex == 1)) eta <- color_naa_dim(eta, us_chol(sim_env$naa_sex_corr, ns), 6)
  # the season correlation spans the active seasons, matching the estimation model's factor
  if(isTRUE(sim_env$NAA_re_season == 1)) {
    ks <- sim_env$naa_re_seas
    eta[,,,ks,,] <- color_naa_dim(array(eta[,,,ks,,], dim = c(np, nr, ny, length(ks), na, ns)),
                                     us_chol(sim_env$naa_season_corr, length(ks)), 4)
  }

  # sigmaNAA is the conditional standard deviation, not the marginal one, so rescale here to get unit variance
  scale <- sim_env$sigmaNAA
  if(code %in% c(2, 4)) scale <- scale / sqrt(1 - rho[["age"]]^2)
  if(code %in% c(3, 4)) scale <- scale / sqrt(1 - rho[["year"]]^2)
  eta <- eta * scale

  # do bias correction
  if(isTRUE(sim_env$bias_correct_pe == 2)) eta <- eta - 0.5 * scale^2 # only if bias correcting "all" processess

  # cells outside the active rectangle stay deterministic
  out <- array(0, dim = dim(eta))
  out[,,sim_env$naa_re_yrs, sim_env$naa_re_seas, sim_env$naa_re_ages,] <-
    eta[,,sim_env$naa_re_yrs, sim_env$naa_re_seas, sim_env$naa_re_ages,]
  out
}
