# Independent implementations of the equations in vignettes/c_model_equations.Rmd.
#
# These are written from the vignette, not from R/. That is the whole point: a
# reimplementation that borrowed the model's own helpers would agree with it
# whatever either of them computed, and would say nothing about whether the model
# computes what it is documented to. Written separately, the two agree only if
# both are right.
#
# They are deliberately slow and literal. Every loop is written out, nothing is
# vectorized, and each function carries the equation it implements so the two can
# be read side by side. Speed does not matter here; legibility against the
# vignette does.
#
# The fixtures these run against carry one population, one region, one season, one
# sex and one fleet, so a reported array reduces to a year-by-age matrix. That
# keeps the reference code close to the printed equations, which are written
# without those subscripts too.

#' Reduce a reported array to a year-by-age matrix
#'
#' @param arr Array with year and age as its only non-singleton margins.
#' @param n_yrs,n_ages Expected extents, checked so a silently reshaped array
#'   fails here rather than comparing the wrong cells.
#'
#' @keywords internal
oracle_ya <- function(arr, n_yrs, n_ages) {
  m <- drop(arr)
  if(!identical(dim(m), c(as.integer(n_yrs), as.integer(n_ages))))
    stop("expected a ", n_yrs, " by ", n_ages, " array after dropping singleton margins, got ",
         paste(dim(m), collapse = "x"))
  m
}

#' Reduce a reported array to a vector over years
#'
#' @keywords internal
oracle_y <- function(arr, n_yrs) {
  v <- as.numeric(drop(arr))
  if(length(v) != n_yrs) stop("expected ", n_yrs, " years, got ", length(v))
  v
}


# ---------------------------------------------------------------------------
# Process equations
# ---------------------------------------------------------------------------

#' Numbers at age carried forward
#'
#' From "Population Projection". With one season, individuals advance in age at
#' the end of the year:
#'
#'   N[y+1, a+1] = N[y, a] exp(-Z[y, a])                      for 1 <= a < a+
#'   N[y+1, a+]  = N[y+1, a+] + N[y, a+] exp(-Z[y, a+])
#'
#' the second line being the plus group accumulating survivors of both the last
#' true age and of itself.
#'
#' @param n1 Numbers at age in the first year.
#' @param rec Recruitment entering age one in years 2..n_yrs.
#' @param Z Year-by-age total mortality.
#'
#' @return A (n_yrs + 1) by n_ages matrix.
#'
#' @keywords internal
oracle_project_naa <- function(n1, rec, Z) {
  n_yrs <- nrow(Z); n_ages <- ncol(Z)
  N <- matrix(0, n_yrs + 1, n_ages)
  N[1, ] <- n1

  for(y in seq_len(n_yrs)) {
    for(a in seq_len(n_ages - 1)) {
      N[y + 1, a + 1] <- N[y, a] * exp(-Z[y, a])
    } # end a loop
    # the plus group takes survivors of the oldest true age and of itself
    N[y + 1, n_ages] <- N[y + 1, n_ages] + N[y, n_ages] * exp(-Z[y, n_ages])
    if(y < n_yrs) N[y + 1, 1] <- rec[y + 1]
  } # end y loop

  N
}

#' Baranov catch at age
#'
#' From "Fishery Observation Model":
#'
#'   C[y, a] = retF[y, a] / Z[y, a] * N[y, a] * (1 - exp(-Z[y, a]))
#'
#' @param retF,Z,N Year-by-age retained fishing mortality, total mortality and
#'   numbers at age.
#'
#' @keywords internal
oracle_baranov <- function(retF, Z, N) {
  out <- matrix(0, nrow(Z), ncol(Z))
  for(y in seq_len(nrow(Z))) {
    for(a in seq_len(ncol(Z))) {
      out[y, a] <- retF[y, a] / Z[y, a] * N[y, a] * (1 - exp(-Z[y, a]))
    } # end a loop
  } # end y loop
  out
}

#' Spawning stock biomass
#'
#' From "Spawning Biomass Timing". The population is carried a fraction t_spawn
#' into the spawning season before it is weighed:
#'
#'   N_spawn[y, a] = N[y, a] exp(-t_spawn Z[y, a])
#'   SSB[y]        = sum_a N_spawn[y, a] WAA[y, a] MatAA[y, a]
#'
#' A single-sex model carries both sexes in one set of numbers, so the sum is
#' halved to leave females. The vignette states this in the sentence after the
#' equation rather than in the equation, which is why it is a separate argument
#' here: the clause is easy to read past.
#'
#' @param sex_ratio Fraction of the population that spawns. One half for a
#'   single-sex model, one when sexes are tracked separately.
#'
#' @keywords internal
oracle_ssb <- function(N, Z, WAA, MatAA, t_spawn, sex_ratio = 0.5) {
  vapply(seq_len(nrow(Z)), function(y) {
    total <- 0
    for(a in seq_len(ncol(Z))) {
      n_spawn <- N[y, a] * exp(-t_spawn * Z[y, a])
      total <- total + n_spawn * WAA[y, a] * MatAA[y, a]
    } # end a loop
    total * sex_ratio
  }, numeric(1))
}

#' Total biomass
#'
#' The same propagated numbers as spawning biomass, weighed without maturity.
#'
#' @keywords internal
oracle_total_biomass <- function(N, Z, WAA, t_spawn) {
  # total biomass counts the whole population, so no sex fraction is applied
  oracle_ssb(N, Z, WAA, MatAA = matrix(1, nrow(WAA), ncol(WAA)), t_spawn = t_spawn,
             sex_ratio = 1)
}


# ---------------------------------------------------------------------------
# Selectivity
# ---------------------------------------------------------------------------

#' Selectivity forms
#'
#' From "Fishery and Survey Selectivity". Each is written on a generic bin b.
#'
#'   logist1:  1 / (1 + exp(-k (b - b50)))
#'   logist2:  1 / (1 + 19^((b50 - b) / b95))
#'   gamma:    p   = 0.5 (sqrt(bmax^2 + 4 delta^2) - bmax)
#'             sel = (b / bmax)^(bmax / p) exp((bmax - b) / p)
#'   power:    1 / b^phi
#'
#' Two things the vignette leaves out. It gives each formula without the order its
#' parameters arrive in, which the model fixes: logist1 takes (b50, slope), not
#' (slope, b50), and the order is as much part of the interface as the curve. And
#' it prints the gamma root as sqrt(bmax + 4 delta^2), without the square on bmax
#' that the code carries; unsquared the expression is not dimensionally consistent
#' and is not the standard reparameterized gamma, so the vignette is what is wrong
#' there.
#'
#' @param form One of "logist1", "logist2", "gamma", "power".
#' @param bins Numeric bin vector.
#' @param pars Parameters on their natural scale.
#'
#' @keywords internal
oracle_selex <- function(form, bins, pars) {
  switch(form,
    logist1 = 1 / (1 + exp(-pars[2] * (bins - pars[1]))),
    logist2 = 1 / (1 + 19^((pars[1] - bins) / pars[2])),
    gamma   = {
      bmax <- pars[1]; delta <- pars[2]
      p <- 0.5 * (sqrt(bmax^2 + 4 * delta^2) - bmax)
      (bins / bmax)^(bmax / p) * exp((bmax - bins) / p)
    },
    power   = 1 / bins^pars[1],
    stop("oracle_selex has no reference implementation for '", form, "'"))
}


# ---------------------------------------------------------------------------
# Likelihoods
# ---------------------------------------------------------------------------

#' Multinomial composition negative log likelihood
#'
#' From "Fishery and Survey Compositions", in the negative log form the vignette
#' gives:
#'
#'   -l = ESS sum_b (O_b + c 1_const) [ log(O_b + c) - log(E_b + c) ]
#'
#' with c the guard constant against log(0) and 1_const controlling whether it is
#' also added to the weights. A perfect fit contributes zero, which is what makes
#' this the offset form rather than the bare multinomial kernel.
#'
#' @param obs,pred Observed and expected proportions.
#' @param ess Effective sample size.
#' @param const Guard constant.
#' @param const_obs Whether the constant is added to the weighting proportions.
#'
#' @keywords internal
oracle_multinomial_nll <- function(obs, pred, ess, const = 0, const_obs = TRUE) {
  total <- 0
  for(b in seq_along(obs)) {
    weight <- obs[b] + if(const_obs) const else 0
    total <- total + weight * (log(obs[b] + const) - log(pred[b] + const))
  } # end b loop
  ess * total
}

#' Effective sample size under the Dirichlet-multinomial
#'
#' From "Fishery and Survey Compositions":
#'
#'   ESS = 1 / (1 + theta) + ISS theta / (1 + theta)
#'
#' @keywords internal
oracle_dm_ess <- function(iss, theta) {
  1 / (1 + theta) + iss * theta / (1 + theta)
}

#' Lognormal observation negative log likelihood
#'
#' From "Fishery Catches" and "Fishery and Survey Indices". The vignette states
#' the density; its negative log, dropping no terms, is
#'
#'   -l = log(sigma) + 0.5 log(2 pi) + (log(obs) - log(pred))^2 / (2 sigma^2)
#'
#' @keywords internal
oracle_lognormal_nll <- function(obs, pred, sigma, include_constant = TRUE) {
  quad <- (log(obs) - log(pred))^2 / (2 * sigma^2)
  quad + log(sigma) + if(include_constant) 0.5 * log(2 * pi) else 0
}
