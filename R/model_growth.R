# Stage 2 of 3: objective function
#
# Parametric growth. Builds mean length at age and its spread from von
# Bertalanffy parameters in the Schnute form (length at a young reference age,
# length at an old reference age, K) or from the Richards generalization of that
# curve, turns them into the size-age transition matrix that length compositions,
# conditional age-at-length and length-based selectivity read, and optionally
# into weight at age through a weight-length relationship and into maturity at
# age through a maturity-at-length curve. Growth is linear from L0 at age zero up
# to L1 at A1, the CV interpolates between the two reference lengths, the
# plus-group mean size is adjusted for fish older than the accumulator age, and
# the age-length key is a binned normal with the tails accumulated into the end
# bins. Any growth parameter can vary over time through deviations, and the
# resulting size at age can either be read off each year's curve or carried
# cohort by cohort, where every cohort grows by the increment the current year's
# parameters imply. The cohort form runs one year at a time from the year loop
# of the population dynamics, because the plus group blends the cohort entering
# it with the fish already there by their numbers.

#' Mean length and its spread at a set of real ages
#'
#' Schnute-form von Bertalanffy growth: \code{L1} is the mean length at age
#' \code{A1}, \code{L2} the mean length at age \code{A2}, and \code{K} the
#' growth rate, so
#' \deqn{L_\infty = L_1 + \frac{L_2 - L_1}{1 - e^{-K(A_2 - A_1)}}}
#' and \eqn{L(x) = L_\infty + (L_1 - L_\infty) e^{-K(x - A_1)}} for real ages
#' \eqn{x \ge A_1}. Below \code{A1} growth is linear from \code{L0} at age zero,
#' \eqn{L(x) = L_0 + (x / A_1)(L_1 - L_0)}, the linear phase.
#' \code{L2_asymptote = 1} reads \code{L2} as the asymptote directly, with no
#' second reference age to solve it from.
#'
#' With a Richards coefficient \code{rho} other than one the curve is the
#' Richards generalization, which applies the same form to the lengths raised to
#' that power,
#' \deqn{L(x)^\rho = L_\infty^\rho + (L_1^\rho - L_\infty^\rho) e^{-K(x - A_1)}}
#' with \eqn{L_\infty^\rho = L_1^\rho + (L_2^\rho - L_1^\rho) / (1 - e^{-K(A_2 - A_1)})}
#' when \code{A2} is a real age. \code{rho = 1} is the von Bertalanffy curve.
#'
#' The coefficient of variation is \code{CV1} below \code{A1}, \code{CV2} at and
#' above \code{A2}, and in between interpolates linearly on mean length
#' (\code{cv_type = 0}) or on age (\code{cv_type = 1}).
#' The spread is \code{CV * L} under \code{sd_type = 0} and the parameter itself
#' under \code{sd_type = 1}.
#'
#' @param x Numeric vector of real ages (data, not parameters).
#' @param L0 Length at age zero, the anchor of the linear phase.
#' @param L1,L2,K,CV1,CV2 Growth parameters, natural scale, possibly AD.
#' @param A1,A2 Reference ages for \code{L1} and \code{L2}. Ignored for the
#'   asymptote under \code{L2_asymptote}, though \code{A2} still bounds the CV
#'   interpolation.
#' @param cv_type Integer, 0 interpolate the CV on length, 1 scale by age.
#' @param sd_type Integer, 0 the CV parameters scale the mean, 1 they are SDs.
#' @param A2_cv Age at and above which \code{CV2} applies. Defaults to \code{A2}.
#'   Under \code{L2_asymptote} there is no second reference age, so
#'   \code{\link{Setup_Mod_Biologicals}} sets \code{A2} to the accumulator age
#'   and the interpolation runs to there.
#' @param rho Richards coefficient, natural scale, possibly AD. One (the default)
#'   is the von Bertalanffy curve.
#' @param cv_ref Optional vector of the coefficient of variation at each element
#'   of \code{x}, used in place of the one this curve implies. Holds the spread
#'   at age at a reference year's while the mean moves, which is the convention
#'   for a time-varying growth curve.
#' @param L2_asymptote Integer, 0 (default) solves \eqn{L_\infty} from
#'   \code{L1} and \code{L2} at their reference ages, 1 reads \code{L2} as
#'   \eqn{L_\infty} itself. Set from \code{growth_A2 = "Linf"} in
#'   \code{\link{Setup_Mod_Biologicals}}.
#'
#' @return List with \code{L} (mean length), \code{sd} (spread), \code{Linf} and
#'   \code{cv}.
#'
#' @keywords internal
#' @import RTMB
get_laa_curve = function(x, L0, L1, L2, K, CV1, CV2, A1, A2, cv_type = 0, sd_type = 0, A2_cv = NULL, rho = 1, cv_ref = NULL,
                         L2_asymptote = 0) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  if(is.null(A2_cv)) A2_cv = A2
  # the asymptote, on the power scale for the Richards form. Under
  # L2_asymptote the second reference length is the asymptote itself, so there
  # is no reference age to solve it from
  LinfR = if(L2_asymptote == 1) L2^rho else L1^rho + (L2^rho - L1^rho) / (1 - exp(-K * (A2 - A1)))
  Linf = LinfR^(1 / rho)

  # A1 of zero puts every age on the curve and L(0) = L1; the guard keeps x / A1 finite.
  A1_plus = if(A1 > 0) A1 else 1e-6
  below = as.numeric(x < A1)
  # ages below A1 are masked out, and are evaluated at A1 so the power stays on a positive base
  on_curve = (LinfR + (L1^rho - LinfR) * exp(-K * pmax(x - A1, 0)))^(1 / rho)
  L = (L0 + (x / A1_plus) * (L1 - L0)) * below + on_curve * (1 - below)

  above = as.numeric(x >= A2_cv)
  mid = 1 - below - above
  cv_mid = if(cv_type == 0) CV1 + (L - L1) * (CV2 - CV1) / (L2 - L1) else CV1 + (x - A1) * (CV2 - CV1) / (A2_cv - A1)
  cv = if(is.null(cv_ref)) CV1 * below + CV2 * above + cv_mid * mid else cv_ref
  sd = if(sd_type == 0) cv * L else cv

  return(list(L = L, sd = sd, Linf = Linf, cv = cv))
}

#' Grow a mean length forward by a fraction of a year
#'
#' The increment of the Richards (or von Bertalanffy, \code{rho = 1}) curve
#' over an elapsed time \code{e} from a current mean length \code{L}:
#' \deqn{L(t + e)^\rho = L_\infty^\rho + (L^\rho - L_\infty^\rho) e^{-K e}}
#' Applied to a length that sits on the curve it returns the curve's value
#' \code{e} later, so splitting a year into seasons changes nothing; applied to a
#' length carried from an earlier year's parameters it is how a cohort keeps its
#' own history.
#'
#' @param L Mean length(s) at the start, possibly AD.
#' @param e Elapsed time in years (data).
#' @param K,Linf,rho Growth rate, asymptote and Richards coefficient in effect.
#'
#' @return Mean length(s) after the increment.
#'
#' @keywords internal
grow_increment = function(L, e, K, Linf, rho = 1) {
  (Linf^rho + (L^rho - Linf^rho) * exp(-K * e))^(1 / rho)
}

#' Mean size of the plus group
#'
#' Adjust for fish older than the accumulator age: the plus group is a
#' mixture of ages from the accumulator age onward, their numbers decaying at
#' \eqn{Z = 0.2} per year, and their size growing linearly from the curve's
#' value at the accumulator age to \eqn{L_\infty} over a second lifetime,
#' \deqn{\bar L_+ = \frac{\sum_{a=0}^{n} e^{-0.2 a}\,[L_n + (a/n)(L_\infty - L_n)]}{\sum_{a=0}^{n} e^{-0.2 a}}}
#' with \eqn{n} the accumulator age.
#'
#' @param L_acc Mean length at the accumulator age from the curve.
#' @param Linf Asymptotic length.
#' @param n_acc The accumulator age.
#'
#' @return The adjusted plus-group mean length.
#'
#' @keywords internal
plus_group_size = function(L_acc, Linf, n_acc) {
  a = 0:n_acc
  w = exp(-0.2 * a)
  return(sum(w * (L_acc + (a / n_acc) * (Linf - L_acc))) / sum(w))
}

#' Binned age-length key
#'
#' \eqn{P(l \mid a)} on the length bins given by their lower edges: the normal (or lognormal) CDF is taken at every lower edge,
#' differenced, and the tails below the first edge and above the last are
#' accumulated into the end bins, so every column sums to one.
#'
#' @param len_lower Numeric vector of lower bin edges.
#' @param mu Mean length at age, possibly AD.
#' @param sd Spread of length at age, possibly AD.
#' @param dist Integer, 0 normal, 1 lognormal (\code{sd} on the log scale, mean
#'   corrected so the arithmetic mean stays \code{mu}).
#'
#' @return Matrix \code{[n_lens x n_ages]}.
#'
#' @keywords internal
#' @import RTMB
get_alk = function(len_lower, mu, sd, dist = 0) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  n_lens = length(len_lower)
  n_ages = length(mu)
  alk = matrix(0, nrow = n_lens, ncol = n_ages)
  for(a in 1:n_ages) {
    z = if(dist == 0) (len_lower - mu[a]) / sd[a] else (log(len_lower) - (log(mu[a]) - 0.5 * sd[a]^2)) / sd[a]
    cdf = RTMB::pnorm(z)
    col = c(cdf[2:n_lens], 1) - cdf # upper tail accumulates into the last bin
    col[1] = col[1] + cdf[1] # lower tail accumulates into the first bin
    alk[,a] = col
  } # end a loop
  return(alk)
}

#' Growth parameters in effect in one year
#'
#' Applies the time-varying deviations of one stratum to its base parameters.
#' Under the log link a parameter in year \eqn{y} is \eqn{P \exp(\delta_y)};
#' under the logit link it is kept inside its bounds,
#' \deqn{P_y = lo + (hi - lo)\,\mathrm{logit}^{-1}\!\left(\log\frac{P - lo}{hi - P} + \delta_y\right)}
#' A random walk's deviation array holds the walk's position, so both
#' structures read the same way here and differ only in their penalty.
#'
#' @param ln_pars Log-scale base parameters of the stratum, length
#'   \code{n_gpars} (five for the von Bertalanffy form, six with the Richards
#'   coefficient last).
#' @param ln_devs Matrix \code{[n_yrs x n_gpars]} of the stratum's deviations.
#' @param tv_model Integer vector \code{[n_gpars]}, 0 constant, 1 iid
#'   deviations, 2 random walk.
#' @param tv_link Integer, 0 log link, 1 logit link within \code{bounds}.
#' @param bounds Matrix \code{[n_gpars x 2]} of lower and upper bounds, read
#'   under the logit link.
#' @param y Year index.
#'
#' @return Natural-scale parameter vector for the year, with a seventh element
#'   \code{rho} of one when the base carries five parameters.
#'
#' @keywords internal
get_growth_pars_year = function(ln_pars, ln_devs, tv_model, tv_link, bounds, y) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  gp = exp(ln_pars)
  n_gpars = length(ln_pars)
  for(k in 1:n_gpars) {
    if(tv_model[k] == 0) next
    dev = ln_devs[y, k]
    if(tv_link == 0) gp[k] = gp[k] * exp(dev)
    else {
      lo = bounds[k, 1]; hi = bounds[k, 2]
      # the 1e-7 keeps the logit finite at a bound
      base_logit = log((gp[k] - lo + 1e-7) / (hi - gp[k] + 1e-7))
      gp[k] = lo + (hi - lo) / (1 + exp(-base_logit - dev))
    }
  } # end k loop
  if(n_gpars == 5) gp = c(gp, 1) # von Bertalanffy: Richards coefficient of one
  return(gp)
}

#' Selection-weighted weight at age
#'
#' The mean weight of the fish a length-selective gear takes at each age,
#' \eqn{\sum_l P(l \mid a) s(l) w(l) / \sum_l P(l \mid a) s(l)}. This is the
#' weight a catch in biomass is made of when selectivity acts on length.
#'
#' @param key Matrix \code{[n_lens x n_ages]}, \eqn{P(l \mid a)}.
#' @param sel_l Selectivity at length, length \code{n_lens}.
#' @param w_len Weight at the bin midpoints, length \code{n_lens}.
#'
#' @return Vector of length \code{n_ages}.
#'
#' @keywords internal
get_selected_waa = function(key, sel_l, w_len) {
  num = as.vector(t(key) %*% (sel_l * w_len))
  den = as.vector(t(key) %*% sel_l)
  return(num / (den + 1e-20))
}

#' Growth state at the start of a year
#'
#' Evaluates the curve a stratum starts from and the quantities that stay fixed
#' over the years: the mean length at the start of every age in the first year
#' (the plus group adjusted), the CV at age and timing when it is held at the
#' first year's sizes, and the asymptote.
#'
#' @keywords internal
growth_start_state = function(gp, ages, growth_A1, growth_A2, growth_L0, growth_cv_type, growth_sd_type,
                              growth_plus_group, growth_L2_asymptote = 0) {

  n_ages = length(ages); n_acc = max(ages)
  L1 = gp[1]; L2 = gp[2]; K = gp[3]; CV1 = gp[4]; CV2 = gp[5]; rho = gp[6]
  crv = get_laa_curve(x = ages, L0 = growth_L0, L1 = L1, L2 = L2, K = K, CV1 = CV1, CV2 = CV2,
                      A1 = growth_A1, A2 = growth_A2, cv_type = growth_cv_type, sd_type = growth_sd_type, rho = rho,
                      L2_asymptote = growth_L2_asymptote)
  L_beg = crv$L
  L_beg[n_ages] = if(growth_plus_group == 1) plus_group_size(crv$L[n_ages], crv$Linf, n_acc) else crv$L[n_ages]
  return(list(L_beg = L_beg, Linf = crv$Linf))
}

#' Mean length and spread of every age at one point in a year
#'
#' The size each integer age has reached at elapsed time \code{e} into year
#' \code{y}. Under curve growth (\code{cohort = 0}) every age is read from the
#' year's curve at its real age, with the plus group grown from its adjusted
#' size. Under cohort growth the ages that are propagated start from the
#' start-of-year state \code{L_beg} and grow by the year's increment; the ages
#' still in the linear phase take the length at \code{A1} their own cohort was
#' born with; and the first integer age past \code{A1} sits on the current
#' year's curve, which is where the propagation picks it up. The CV at age is
#' that of the year's own curve under curve growth and is held at the first
#' year's sizes under cohort growth.
#'
#' @param len_devs Optional vector of log deviations on mean length at age, one
#'   per age, or \code{NULL} for a purely parametric curve.
#' @keywords internal
growth_laa_at = function(e, gp, ages, growth_A1, growth_A2, growth_L0, growth_cv_type, growth_sd_type,
                         cohort, L_beg, L1_birth, cv_ref, a_prop, len_devs = NULL, growth_L2_asymptote = 0) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  n_ages = length(ages)
  L1 = gp[1]; L2 = gp[2]; K = gp[3]; CV1 = gp[4]; CV2 = gp[5]; rho = gp[6]
  x = ages + e
  crv = get_laa_curve(x = x, L0 = growth_L0, L1 = L1, L2 = L2, K = K, CV1 = CV1, CV2 = CV2,
                      A1 = growth_A1, A2 = growth_A2, cv_type = growth_cv_type, sd_type = growth_sd_type,
                      A2_cv = growth_A2, rho = rho, cv_ref = cv_ref, L2_asymptote = growth_L2_asymptote)
  Linf = crv$Linf
  L = crv$L

  if(cohort == 0) {
    # the plus group grows from its adjusted size rather than from the curve
    L[n_ages] = grow_increment(L_beg[n_ages], e, K, Linf, rho)
  } else {
    for(i in 1:n_ages) {
      a = ages[i]
      if(x[i] < growth_A1) L[i] = growth_L0 + (x[i] / (if(growth_A1 > 0) growth_A1 else 1e-6)) * (L1_birth[i] - growth_L0)
      else if(x[i] == growth_A1) L[i] = L1_birth[i]
      else if(a >= a_prop) L[i] = grow_increment(L_beg[i], e, K, Linf, rho)
      # else: past A1 but not yet propagated, the current year's curve already in L
    } # end i loop
  }

  # Semi-parametric growth: the parametric mean at age times a year-by-age
  # deviation surface. Applied after the curve and any cohort propagation, so
  # the deviations describe departures from paraemtric growth
  if(!is.null(len_devs)) L = L * exp(len_devs)

  # the spread follows the CV rule at the size reached, or the first year's CV at age
  if(is.null(cv_ref)) {
    above = as.numeric(x >= growth_A2)
    below = as.numeric(x < growth_A1)
    mid = 1 - above - below
    cv_mid = if(growth_cv_type == 0) CV1 + (L - L1) * (CV2 - CV1) / (L2 - L1) else crv$cv
    cv = CV1 * below + CV2 * above + cv_mid * mid
  } else cv = cv_ref
  sd = if(growth_sd_type == 0) cv * L else cv

  return(list(L = L, sd = sd, cv = cv, Linf = Linf))
}

#' Growth module
#'
#' Builds the size-age transition arrays and, when asked, weight at age and
#' maturity at age from the growth parameters, for every year. Every fleet gets
#' its own key and weight, read at that fleet's timing within the season: each
#' fishery fleet's at \code{t_fish}, each survey's at \code{t_srv}, and the
#' spawning weight at \code{t_spawn}.
#'
#' @section Time variation:
#' With no deviations the key is built once and broadcast over the years. With
#' deviations on any parameter and \code{growth_tv_type = 0} ("curve") every
#' year's sizes are read from that year's own curve. Under
#' \code{growth_tv_type = 1} ("cohort") the sizes are carried forward cohort by
#' cohort from \code{growth_cohort_styr} on, which needs the numbers at age of
#' each year to blend the plus group; that form is run one year at a time from
#' the population dynamics through \code{\link{Get_Growth_Year}}, and this
#' function only evaluates the years before the propagation starts, which all
#' sit on the first year's curve.
#'
#' @section Timing within the year:
#' A season starts at the cumulative
#' duration of the seasons before it , and a point inside
#' a season is that start plus the fraction of the season elapsed times the
#' season's duration , so every evaluation is anchored at
#' the season start rather than compounded from the previous evaluation. Mean
#' length is the curve read at the real age, the integer age plus that elapsed
#' fraction of a year, and the plus group grows from its adjusted size by the
#' curve's increment over the same elapsed time. Under constant growth
#' the increment over a season is exactly the curve's difference between the two
#' real ages, so splitting the year into seasons of any durations changes no mean
#' length. The spawning
#' weight is read at the spawning fraction of the season, and each fleet's key
#' and weight at that fleet's own timing, \code{t_fish} or \code{t_srv}, so a
#' composition and the weight behind an index are formed at the point in the
#' season the observation is taken. Fleets that share a timing share one
#' evaluation. Seasonal multipliers on \eqn{K} are not carried.
#'
#' @param ln_growth_pars Array \code{[pop, region, sex, n_gpars]} of log growth
#'   parameters in the order L1, L2, K, CV1, CV2 and, for the Richards form, rho.
#' @param ln_growth_devs Array \code{[pop, region, year, n_gpars, sex]} of
#'   time-varying deviations, or \code{NULL} for none.
#' @param growth_tv_model Integer vector \code{[n_gpars]}, 0 constant, 1 iid,
#'   2 random walk, per parameter.
#' @param growth_tv_link Integer, 0 log link, 1 logit link within
#'   \code{growth_par_bounds}.
#' @param growth_par_bounds Matrix \code{[n_gpars x 2]} of bounds for the logit
#'   link.
#' @param growth_tv_type Integer, 0 each year on its own curve, 1 cohort
#'   propagation.
#' @param ln_growth_semipar_devs Array \code{[pop, region, year, age, sex]} of
#'   log deviations on mean length at age, or \code{NULL} for none. Multiplies
#'   the parametric mean at age, so the curve stays the parametric part and the
#'   deviations carry departures from it; the spread at age follows the deviated
#'   mean, leaving the coefficient of variation at age alone.
#' @param growth_semipar Integer process error code for those deviations,
#'   \code{0} for none. Only its being nonzero is read here; the structure is
#'   scored in the objective.
#' @param growth_cohort_styr Year index the cohort propagation starts from.
#' @param growth_A1,growth_A2 Reference ages; see \code{\link{get_laa_curve}}.
#' @param growth_L0 Length at age zero.
#' @param growth_len_lower Lower edges of the length bins.
#' @param t_fish Array \code{[region, season, fish_fleet]} of fishery timings,
#'   the fraction of the season elapsed at which each fishery fleet's key and
#'   weight at age are read.
#' @param t_srv Array \code{[region, season, srv_fleet]} of survey timings, the
#'   fraction of the season elapsed when each survey is taken, at which that
#'   survey's key and weight at age are read.
#' @param growth_cv_type,growth_sd_type,growth_dist Integer switches passed to
#'   \code{\link{get_laa_curve}} and \code{\link{get_alk}}.
#' @param growth_plus_group Integer (0/1); apply \code{\link{plus_group_size}}.
#' @param derive_waa Integer (0/1); build weight at age from the key.
#' @param wt_len_pars Array \code{[pop, region, sex, 2]} of weight-length
#'   parameters \eqn{a, b} in \eqn{W = a L^b}, read when \code{derive_waa = 1}.
#' @param ages Numeric vector of model ages.
#' @param seasdur Season durations as fractions of a year.
#' @param spawn_seas,t_spawn Spawning season and fraction of it elapsed at spawning.
#' @param n_pop,n_regions,n_yrs,n_seas,n_sexes,n_fish_fleets,n_srv_fleets Dimensions.
#' @param years_eval Integer vector of the years to evaluate; the default is
#'   every year.
#'
#' @return List with \code{SizeAgeTrans_fish} and \code{SizeAgeTrans_srv}
#'   \code{[pop, region, year, season, len, age, sex, fleet]}, one key per fleet
#'   at that fleet's timing, the spawning key \code{SizeAgeTrans_spawn} \code{[pop, region, year, len, age, sex]},
#'   \code{mean_LAA_fish}, \code{sd_LAA_fish}, \code{mean_LAA_srv},
#'   \code{sd_LAA_srv} \code{[pop, region, year, season, age, sex, fleet]},
#'   \code{mean_LAA_spawn} and \code{sd_LAA_spawn} \code{[pop, region, year,
#'   season, age, sex]}, \code{L_beg} \code{[pop, region, year, age, sex]} (the
#'   start-of-year mean length), \code{Linf} \code{[pop, region, year, sex]},
#'   \code{growth_pars_y} \code{[pop, region, year, 6, sex]}, and when
#'   \code{derive_waa = 1} also \code{WAA}, \code{WAA_fish} and \code{WAA_srv}.
#'
#' @keywords internal
#' @import RTMB
Get_Growth = function(ln_growth_pars, growth_A1, growth_A2, growth_L0, growth_len_lower,
                      growth_cv_type, growth_sd_type, growth_dist, growth_plus_group, growth_L2_asymptote = 0,
                      derive_waa, wt_len_pars, ages, seasdur, spawn_seas, t_spawn,
                      n_pop, n_regions, n_yrs, n_seas, n_sexes, n_fish_fleets, n_srv_fleets,
                      t_fish, t_srv,
                      ln_growth_devs = NULL, growth_tv_model = NULL, growth_tv_link = 0, growth_par_bounds = NULL,
                      growth_tv_type = 0, growth_cohort_styr = 1, years_eval = NULL,
                      ln_growth_semipar_devs = NULL, growth_semipar = 0) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  n_ages = length(ages)
  n_lens = length(growth_len_lower)
  n_gpars = dim(ln_growth_pars)[4]
  if(is.null(growth_tv_model)) growth_tv_model = rep(0, n_gpars)
  tv_any = any(growth_tv_model > 0)
  semipar = growth_semipar > 0 && !is.null(ln_growth_semipar_devs)
  # under cohort growth this function only builds the years before the propagation starts
  if(is.null(years_eval)) years_eval = if(growth_tv_type == 1) seq_len(growth_cohort_styr) else 1:n_yrs
  if(length(years_eval) == 0) years_eval = 1

  out = growth_containers(n_pop, n_regions, n_yrs, n_seas, n_lens, n_ages, n_sexes, n_fish_fleets, n_srv_fleets, derive_waa)

  for(p in 1:n_pop) {
    for(r in 1:n_regions) {
      for(s in 1:n_sexes) {

        devs_prs = if(is.null(ln_growth_devs)) matrix(0, n_yrs, n_gpars) else matrix(ln_growth_devs[p,r,,,s], n_yrs, n_gpars)
        # parameters of the first year set the state every year starts from under
        # constant or curve growth, and the first year's CV under cohort growth
        gp_1 = get_growth_pars_year(ln_growth_pars[p,r,s,], devs_prs, growth_tv_model, growth_tv_link, growth_par_bounds, 1)
        start = growth_start_state(gp_1, ages, growth_A1, growth_A2, growth_L0, growth_cv_type, growth_sd_type, growth_plus_group, growth_L2_asymptote)
        w_mid = wt_len_pars[p,r,s,1] * growth_len_mid(growth_len_lower)^wt_len_pars[p,r,s,2] # weight at the bin midpoints

        # years that share one curve are evaluated once and propagated out; a
        # deviation surface on mean length at age makes every year its own
        constant = (!tv_any || growth_tv_type == 1) && !semipar
        for(y in years_eval) {
          gp_y = if(constant) gp_1 else get_growth_pars_year(ln_growth_pars[p,r,s,], devs_prs, growth_tv_model, growth_tv_link, growth_par_bounds, y)
          st_y = if(constant) start else growth_start_state(gp_y, ages, growth_A1, growth_A2, growth_L0, growth_cv_type, growth_sd_type, growth_plus_group, growth_L2_asymptote)
          fill_yrs = if(constant) years_eval else y
          out = growth_fill_year(out, p, r, s, fill_yrs, gp_y, st_y$L_beg, L1_birth = NULL, cv_ref_fn = NULL, a_prop = NULL, cohort = 0,
                                 ages, growth_A1, growth_A2, growth_L0, growth_len_lower, growth_cv_type, growth_sd_type, growth_dist,
                                 derive_waa, w_mid, seasdur, spawn_seas, t_spawn, n_seas, t_fish, t_srv,
                                 len_devs = if(semipar) ln_growth_semipar_devs[p,r,y,,s] else NULL,
                                 growth_L2_asymptote = growth_L2_asymptote)
          if(constant) break
        } # end y loop

      } # end s loop
    } # end r loop
  } # end p loop

  return(out)
}

#' Bin midpoints from lower edges, the last bin taking the width of the one before
#' @keywords internal
growth_len_mid = function(len_lower) {
  n = length(len_lower)
  widths = diff(len_lower)
  len_lower + c(widths, widths[n - 1]) / 2
}

#' Empty containers for the growth module's output
#' @keywords internal
growth_containers = function(n_pop, n_regions, n_yrs, n_seas, n_lens, n_ages, n_sexes, n_fish_fleets, n_srv_fleets, derive_waa) {
  out = list(
    SizeAgeTrans_fish = array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_lens, n_ages, n_sexes, n_fish_fleets)),
    SizeAgeTrans_srv = array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_lens, n_ages, n_sexes, n_srv_fleets)),
    SizeAgeTrans_spawn = array(0, dim = c(n_pop, n_regions, n_yrs, n_lens, n_ages, n_sexes)),
    mean_LAA_fish = array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_fish_fleets)),
    sd_LAA_fish = array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_fish_fleets)),
    mean_LAA_srv = array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_srv_fleets)),
    sd_LAA_srv = array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_srv_fleets)),
    mean_LAA_spawn = array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes)),
    sd_LAA_spawn = array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes)),
    L_beg = array(0, dim = c(n_pop, n_regions, n_yrs, n_ages, n_sexes)),
    Linf = array(0, dim = c(n_pop, n_regions, n_yrs, n_sexes)),
    growth_pars_y = array(0, dim = c(n_pop, n_regions, n_yrs, 6, n_sexes))
  )
  if(derive_waa == 1) {
    out$WAA = array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes))
    out$WAA_fish = array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_fish_fleets))
    out$WAA_srv = array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_srv_fleets))
  }
  out
}

#' Fill one stratum's growth output for a set of years from one start-of-year state
#'
#' Evaluates every timing any fleet, composition or the spawning reads at, builds
#' the keys and the derived weight and maturity, and writes them into the years
#' named, which all share the state handed in.
#'
#' @keywords internal
growth_fill_year = function(out, p, r, s, fill_yrs, gp, L_beg, L1_birth, cv_ref_fn, a_prop, cohort,
                            ages, growth_A1, growth_A2, growth_L0, growth_len_lower, growth_cv_type, growth_sd_type, growth_dist,
                            derive_waa, w_mid, seasdur, spawn_seas, t_spawn, n_seas, t_fish, t_srv, len_devs = NULL,
                            growth_L2_asymptote = 0) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  n_ages = length(ages)
  n_yr_fill = length(fill_yrs)
  cum_before = c(0, cumsum(seasdur))[1:n_seas] # fraction of the year elapsed at each season start
  rep_yrs = function(x) rep(x, each = n_yr_fill) # lays the years fastest for a [years, ...] block

  for(y in fill_yrs) { out$L_beg[p,r,y,,s] = L_beg; out$growth_pars_y[p,r,y,,s] = gp }

  for(seas in 1:n_seas) {

    # Every point in the season something is read at: each fishery fleet's key
    # and weight at its t_fish, each survey's at its t_srv, and the spawning
    # weight at t_spawn in the spawning season. Each distinct point is
    # evaluated once.
    t_fish_seas = as.vector(t_fish[r, seas, ]); t_srv_seas = as.vector(t_srv[r, seas, ])
    t_unique = unique(c(t_fish_seas, t_srv_seas, if(seas == spawn_seas) t_spawn))

    for(k in seq_along(t_unique)) {

      tk = t_unique[k]
      elapsed = cum_before[seas] + tk * seasdur[seas] # real age offset of every age class at this point
      cv_ref = if(is.null(cv_ref_fn)) NULL else cv_ref_fn(elapsed)
      laa = growth_laa_at(e = elapsed, gp = gp, ages = ages, growth_A1 = growth_A1, growth_A2 = growth_A2, growth_L0 = growth_L0,
                          growth_cv_type = growth_cv_type, growth_sd_type = growth_sd_type, cohort = cohort,
                          L_beg = L_beg, L1_birth = L1_birth, cv_ref = cv_ref, a_prop = a_prop, len_devs = len_devs,
                          growth_L2_asymptote = growth_L2_asymptote)
      L = laa$L; sd = laa$sd
      for(y in fill_yrs) out$Linf[p,r,y,s] = laa$Linf
      alk = get_alk(growth_len_lower, L, sd, dist = growth_dist) # n_lens x n_ages (get alk for sizeage)
      waa = if(derive_waa == 1) as.vector(t(alk) %*% w_mid) else NULL

      # get fishery
      for(f in which(t_fish_seas == tk)) {
        out$SizeAgeTrans_fish[p,r,fill_yrs,seas,,,s,f] = rep_yrs(alk)
        out$mean_LAA_fish[p,r,fill_yrs,seas,,s,f] = rep_yrs(L)
        out$sd_LAA_fish[p,r,fill_yrs,seas,,s,f] = rep_yrs(sd)
        if(derive_waa == 1) out$WAA_fish[p,r,fill_yrs,seas,,s,f] = rep_yrs(waa)
      } # end f loop

      # get survey
      for(sf in which(t_srv_seas == tk)) {
        out$SizeAgeTrans_srv[p,r,fill_yrs,seas,,,s,sf] = rep_yrs(alk)
        out$mean_LAA_srv[p,r,fill_yrs,seas,,s,sf] = rep_yrs(L)
        out$sd_LAA_srv[p,r,fill_yrs,seas,,s,sf] = rep_yrs(sd)
        if(derive_waa == 1) out$WAA_srv[p,r,fill_yrs,seas,,s,sf] = rep_yrs(waa)
      } # end sf loop

      # get spawning
      if(seas == spawn_seas && tk == t_spawn) {
        out$SizeAgeTrans_spawn[p,r,fill_yrs,,,s] = rep_yrs(alk)
        out$mean_LAA_spawn[p,r,fill_yrs,seas,,s] = rep_yrs(L)
        out$sd_LAA_spawn[p,r,fill_yrs,seas,,s] = rep_yrs(sd)
        if(derive_waa == 1) out$WAA[p,r,fill_yrs,seas,,s] = rep_yrs(waa)
      }

    } # end k loop
  } # end seas loop

  out
}

#' Cohort growth, one year at a time
#'
#' The year-loop companion of \code{\link{Get_Growth}} under cohort
#' propagation. For year \code{y} it builds every key, length, weight and
#' weight of the year from the start-of-year state, and carries the state to
#' the next year with the year's parameters: each propagated age grows by the
#' year's increment, the first propagated age is placed on the year's curve, the
#' ages in the linear phase keep the length at \code{A1} their cohort was born
#' with, and the plus group's size next year blends the cohort entering it with
#' the fish already there by their numbers at the start of this year,
#' \deqn{\bar L_{+,y+1} = \frac{(N_{n-1} + 0.01)\, g(L_{n-1}) + (N_{n} + 0.01)\, g(L_{+})}{N_{n-1} + N_n + 0.02}}
#' with \eqn{g} the year's increment and \eqn{N} the start-of-year numbers of
#' the stratum. The CV at age is held at the first year's sizes.
#'
#' @param growth The list \code{\link{Get_Growth}} returned, carrying the
#'   output containers and the years before the propagation started.
#' @param y Year index to evaluate.
#' @param NAA_y Numbers at age at the start of year \code{y}, array
#'   \code{[pop, region, age, sex]}, read for the plus-group blend.
#' @param ln_growth_pars,ln_growth_devs,growth_tv_model,growth_tv_link,growth_par_bounds
#'   As in \code{\link{Get_Growth}}.
#' @param growth_A1,growth_A2,growth_L0,growth_len_lower,growth_cv_type,growth_sd_type,growth_dist,growth_plus_group
#'   As in \code{\link{Get_Growth}}.
#' @param derive_waa,wt_len_pars,ages,seasdur,spawn_seas,t_spawn,n_pop,n_regions,n_seas,n_sexes
#'   As in \code{\link{Get_Growth}}.
#' @param t_fish,t_srv Timings as in \code{\link{Get_Growth}}.
#'
#' @return The \code{growth} list with year \code{y} filled and, when a year
#'   follows, \code{L_beg[,, y + 1,,]} set to the state the next year starts
#'   from.
#'
#' @keywords internal
#' @import RTMB
Get_Growth_Year = function(growth, y, NAA_y, ln_growth_pars, ln_growth_devs, growth_tv_model, growth_tv_link, growth_par_bounds,
                           growth_A1, growth_A2, growth_L0, growth_len_lower, growth_cv_type, growth_sd_type, growth_dist,
                           growth_plus_group, growth_L2_asymptote = 0, derive_waa, wt_len_pars, ages, seasdur, spawn_seas, t_spawn,
                           n_pop, n_regions, n_seas, n_sexes, t_fish, t_srv,
                           ln_growth_semipar_devs = NULL, growth_semipar = 0) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  n_ages = length(ages); n_acc = max(ages); n_yrs = dim(growth$L_beg)[3]
  n_gpars = dim(ln_growth_pars)[4]
  # the first integer age whose start-of-year size is carried forward; the age
  # before it sits on the current year's curve
  a_prop = ceiling(growth_A1) + 1
  cum_before = c(0, cumsum(seasdur))[1:n_seas]

  for(p in 1:n_pop) {
    for(r in 1:n_regions) {
      for(s in 1:n_sexes) {

        devs_prs = if(is.null(ln_growth_devs)) matrix(0, n_yrs, n_gpars) else matrix(ln_growth_devs[p,r,,,s], n_yrs, n_gpars)
        pars_at = function(yy) get_growth_pars_year(ln_growth_pars[p,r,s,], devs_prs, growth_tv_model, growth_tv_link, growth_par_bounds, yy)
        gp_1 = pars_at(1)
        gp_y = pars_at(y)
        w_mid = wt_len_pars[p,r,s,1] * growth_len_mid(growth_len_lower)^wt_len_pars[p,r,s,2]

        # the length at A1 each cohort was born with: its birth year's L1, or the
        # first year's for cohorts older than the model
        L1_birth = rep(0, n_ages)
        for(i in 1:n_ages) L1_birth[i] = pars_at(max(1, y - ages[i]))[1]

        # The spread at age is held at the first year's, evaluated from that
        # year's own curve and parameters at each timing, and then carried
        # unchanged while the mean moves.
        start = growth_start_state(gp_1, ages, growth_A1, growth_A2, growth_L0, growth_cv_type, growth_sd_type, growth_plus_group, growth_L2_asymptote)
        cv_ref_fn = function(elapsed) growth_laa_at(e = elapsed, gp = gp_1, ages = ages, growth_A1 = growth_A1, growth_A2 = growth_A2,
                                                    growth_L0 = growth_L0, growth_cv_type = growth_cv_type, growth_sd_type = growth_sd_type,
                                                    cohort = 0, L_beg = start$L_beg, L1_birth = NULL, cv_ref = NULL, a_prop = a_prop,
                                                    growth_L2_asymptote = growth_L2_asymptote)$cv

        L_beg = growth$L_beg[p,r,y,,s]
        growth = growth_fill_year(growth, p, r, s, y, gp_y, L_beg, L1_birth, cv_ref_fn, a_prop, cohort = 1,
                                  ages, growth_A1, growth_A2, growth_L0, growth_len_lower, growth_cv_type, growth_sd_type, growth_dist,
                                  derive_waa, w_mid, seasdur, spawn_seas, t_spawn, n_seas, t_fish, t_srv,
                                  len_devs = if(growth_semipar > 0 && !is.null(ln_growth_semipar_devs)) ln_growth_semipar_devs[p,r,y,,s] else NULL,
                                  growth_L2_asymptote = growth_L2_asymptote)

        # carry the state to next year with this year's increment
        if(y < n_yrs) {
          K = gp_y[3]; rho = gp_y[6]; Linf = growth$Linf[p,r,y,s]
          L_next = growth_laa_at(e = 0, gp = gp_y, ages = ages, growth_A1 = growth_A1, growth_A2 = growth_A2, growth_L0 = growth_L0,
                                 growth_cv_type = growth_cv_type, growth_sd_type = growth_sd_type, cohort = 0,
                                 L_beg = L_beg, L1_birth = NULL, cv_ref = NULL, a_prop = a_prop,
                                 growth_L2_asymptote = growth_L2_asymptote)$L # the year's curve at integer ages
          grown = grow_increment(L_beg, 1, K, Linf, rho) # every age one year later
          for(i in 1:n_ages) {
            a = ages[i]
            if(a < a_prop) next # rebuilt from the linear phase and the curve next year
            if(a == a_prop) { L_next[i] = L_next[i]; next } # the first propagated age starts on this year's curve
            if(a < n_acc) L_next[i] = grown[i - 1]
          } # end i loop
          # the plus group: the cohort entering it and the fish already there,
          # blended by their numbers at the start of this year
          N_in = NAA_y[p,r,n_ages - 1,s]; N_old = NAA_y[p,r,n_ages,s]
          L_next[n_ages] = ((N_in + 0.01) * grown[n_ages - 1] + (N_old + 0.01) * grown[n_ages]) / (N_in + N_old + 0.02)
          growth$L_beg[p,r,y + 1,,s] = L_next
        }

      } # end s loop
    } # end r loop
  } # end p loop

  growth
}


#' Assemble the growth module's arguments from the model's data and parameters
#'
#' The growth settings are the same for every call, whether the whole series is
#' built up front or one year at a time from inside the population loop, so they
#' are gathered once here rather than written out at each call site.
#'
#' @param env Environment holding the unpacked data and parameters, i.e. the
#'   \code{SPoRC_rtmb} frame after \code{RTMB::getAll}. Defaults to the caller.
#'
#' @return A named list of arguments for \code{\link{Get_Growth}} and
#'   \code{\link{Get_Growth_Year}}.
#'
#' @keywords internal
growth_args_from_model = function(env = parent.frame()) {
  nm = c("ln_growth_pars", "growth_A1", "growth_A2", "growth_L0", "growth_len_lower",
         "growth_cv_type", "growth_sd_type", "growth_dist", "growth_plus_group", "growth_L2_asymptote",
         "derive_waa", "wt_len_pars", "ages", "seasdur", "spawn_seas", "t_spawn",
         "n_pop", "n_regions", "n_seas", "n_sexes", "t_fish", "t_srv",
         "ln_growth_devs", "growth_tv_model", "growth_tv_link", "growth_par_bounds",
         "ln_growth_semipar_devs", "growth_semipar")
  stats::setNames(lapply(nm, function(x) get(x, envir = env)), nm)
}

#' Copy one year of the growth module's output into the model's own arrays
#'
#' Under cohort growth a year's size-age transition, lengths, weights and growth parameters are
#' only known once the population loop reaches that year, so they are written
#' into the arrays the rest of the model reads one year at a time.
#'
#' @param dest Named list of the model's growth arrays, as
#'   \code{\link{growth_take_year}} returns them.
#' @param g The growth module's output, from \code{\link{Get_Growth_Year}}.
#' @param y Year index to copy.
#' @param derive_waa Integer (0/1); whether the weight arrays are present.
#'
#' @return \code{dest} with year \code{y} filled.
#'
#' @keywords internal
#' @import RTMB
growth_take_year = function(dest, g, y, derive_waa) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  dest$SizeAgeTrans_fish[,,y,,,,,] = g$SizeAgeTrans_fish[,,y,,,,,]
  dest$SizeAgeTrans_srv[,,y,,,,,] = g$SizeAgeTrans_srv[,,y,,,,,]
  dest$SizeAgeTrans_spawn[,,y,,,] = g$SizeAgeTrans_spawn[,,y,,,]
  dest$mean_LAA_fish[,,y,,,,] = g$mean_LAA_fish[,,y,,,,]; dest$sd_LAA_fish[,,y,,,,] = g$sd_LAA_fish[,,y,,,,]
  dest$mean_LAA_srv[,,y,,,,] = g$mean_LAA_srv[,,y,,,,]; dest$sd_LAA_srv[,,y,,,,] = g$sd_LAA_srv[,,y,,,,]
  dest$mean_LAA_spawn[,,y,,,] = g$mean_LAA_spawn[,,y,,,]; dest$sd_LAA_spawn[,,y,,,] = g$sd_LAA_spawn[,,y,,,]
  dest$Linf[,,y,] = g$Linf[,,y,]; dest$L_beg[,,y,,] = g$L_beg[,,y,,]
  dest$growth_pars_y[,,y,,] = g$growth_pars_y[,,y,,]
  if(derive_waa == 1) {
    dest$WAA[,,y,,,] = g$WAA[,,y,,,]
    dest$WAA_fish[,,y,,,,] = g$WAA_fish[,,y,,,,]
    dest$WAA_srv[,,y,,,,] = g$WAA_srv[,,y,,,,]
  }
  return(dest)
}

#' Selection-weighted weight at age for one year
#'
#' A length-selective gear does not take fish evenly across an age, so the mean
#' weight of what it takes is the weight averaged over the fleet's key
#' re-weighted by selectivity, not the population mean. Applied only to the
#' fleets that ask for it. Used for fishery and survey fleets alike.
#'
#' @param WAA_fleet Array \code{[pop, region, year, season, age, sex, fleet]} of
#'   the fleet type's weight at age, returned with year \code{y} overwritten
#'   for the fleets named in \code{waa_selected}.
#' @param SizeAgeTrans_fleet Each fleet's size-age key.
#' @param sel_l Array \code{[region, year, len, sex, fleet]} of selectivity at
#'   length.
#' @param wt_len_pars Array \code{[pop, region, sex, 2]} of weight-length
#'   parameters.
#' @param len_mid Bin midpoints the weight-length relationship is read at.
#' @param waa_selected Integer vector \code{[fleet]} (0/1).
#' @param y Year index.
#' @param n_pop,n_regions,n_seas,n_sexes Dimensions.
#'
#' @return \code{WAA_fleet} with year \code{y} updated.
#'
#' @keywords internal
#' @import RTMB
growth_selected_waa_year = function(WAA_fleet, SizeAgeTrans_fleet, sel_l, wt_len_pars, len_mid,
                                    waa_selected, y, n_pop, n_regions, n_seas, n_sexes) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  for(f in which(waa_selected == 1)) {
    for(p in 1:n_pop) {
      for(r in 1:n_regions) {
        for(seas in 1:n_seas) {
          for(s in 1:n_sexes) {

            w_mid = wt_len_pars[p,r,s,1] * len_mid^wt_len_pars[p,r,s,2]
            WAA_fleet[p,r,y,seas,,s,f] = get_selected_waa(SizeAgeTrans_fleet[p,r,y,seas,,,s,f], sel_l[r,y,,s,f], w_mid)

          } # end s loop
        } # end seas loop
      } # end r loop
    } # end p loop
  } # end f loop

  return(WAA_fleet)
}
