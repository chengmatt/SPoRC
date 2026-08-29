# The analytic gradient against a finite difference of the objective, and the
# score against zero on data the model generated itself.
#
# Two things go wrong in an AD objective that no comparison of jnLL to a stored
# number will show. A branch taken on the value of a parameter, an assignment into
# an AD array that drops the tape, or a cached tape evaluated off-tape all leave
# the objective returning the right number and its derivative returning the wrong
# one. The optimizer then walks a surface the reported likelihood does not
# describe, converges somewhere, and reports success.
#
# A gradient that disagrees with a finite difference of the very function it is
# supposed to differentiate is the direct statement of that fault.

#' Finite-difference gradient of the objective in a few coordinates
#'
#' A central difference on the coordinates with the largest analytic gradient,
#' which are the ones an optimizer will actually move along and the ones where a
#' dropped tape shows up as a difference rather than as two small numbers.
#'
#' @param obj Object from \code{fit_model(do_optim = FALSE)}.
#' @param k Number of coordinates to check.
#' @param h Step size, applied relative to the parameter's own magnitude.
#'
#' @keywords internal
fd_gradient_check <- function(obj, k = 6, h = 1e-5) {
  p <- obj$par
  g <- as.numeric(obj$gr(p))
  idx <- utils::head(order(abs(g), decreasing = TRUE), k)

  fd <- vapply(idx, function(i) {
    step <- h * max(abs(p[i]), 1)
    up <- p; up[i] <- up[i] + step
    dn <- p; dn[i] <- dn[i] - step
    (obj$fn(up) - obj$fn(dn)) / (2 * step)
  }, numeric(1))

  list(analytic = g[idx], numeric = fd, index = idx)
}

#' Configurations to differentiate
#'
#' One per feature that changes the shape of the objective, so a branch that is
#' only reachable under one option is still differentiated.
#'
#' @keywords internal
gradient_configs <- function() {
  list(
    single_region = list(dims = list(n_regions = 1, n_sexes = 1, n_fish_fleets = 1)),
    multi_region  = list(dims = list(n_regions = 3, n_sexes = 1, n_fish_fleets = 1)),
    two_sex       = list(dims = list(n_regions = 1, n_sexes = 2, n_fish_fleets = 1)),
    multi_fleet   = list(dims = list(n_regions = 1, n_sexes = 1, n_fish_fleets = 3)),
    seasonal      = list(dims = list(n_regions = 2, n_sexes = 1, n_fish_fleets = 1, n_seas = 2)),
    beverton_holt = list(dims = list(n_regions = 1, n_sexes = 1, n_fish_fleets = 1),
                         rec = list(rec_model = "bh_rec", h_spec = "fix")),
    dirichlet_mn  = list(dims = list(n_regions = 1, n_sexes = 1, n_fish_fleets = 1),
                         comps = "Dirichlet-Multinomial")
  )
}

#' Build and tape one gradient configuration
#'
#' @keywords internal
gradient_obj <- function(cfg) {
  dims <- utils::modifyList(list(n_srv_fleets = 1, n_yrs = 10, n_ages = 6), cfg$dims)
  nr <- dims$n_regions; nx <- dims$n_sexes; nf <- dims$n_fish_fleets
  ns <- if(is.null(dims$n_seas)) 1 else dims$n_seas
  NY <- dims$n_yrs; NAG <- dims$n_ages
  like <- if(is.null(cfg$comps)) "Multinomial" else cfg$comps

  il <- sweep_input(
    dims = dims,
    rec = if(is.null(cfg$rec)) list() else cfg$rec,
    catch = list(ObsCatch = array(1e4 / (nr * nf * ns), dim = c(nr, NY, ns, nf)),
                 UseCatch = array(1, dim = c(nr, NY, ns, nf))),
    fishidx = list(
      ObsFishIdx = array(1e5, dim = c(nr, NY, ns, nf)),
      ObsFishIdx_SE = array(0.2, dim = c(nr, NY, ns, nf)),
      UseFishIdx = array(1, dim = c(nr, NY, ns, nf)), fish_idx_type = rep("biom", nf),
      ObsFishAgeComps = array(1 / NAG, dim = c(nr, NY, ns, NAG, nx, nf)),
      UseFishAgeComps = array(1, dim = c(nr, NY, ns, nf)),
      ISS_FishAgeComps = array(100, dim = c(nr, NY, ns, nx, nf)),
      FishAgeComps_LikeType = rep(like, nf),
      ObsFishLenComps = array(0, dim = c(nr, NY, ns, 1, nx, nf)),
      UseFishLenComps = array(0, dim = c(nr, NY, ns, nf)),
      ISS_FishLenComps = array(0, dim = c(nr, NY, ns, nx, nf)),
      FishLenComps_LikeType = rep("none", nf)),
    srvidx = list(
      ObsSrvIdx = array(1e5, dim = c(nr, NY, ns, 1)),
      ObsSrvIdx_SE = array(0.2, dim = c(nr, NY, ns, 1)),
      UseSrvIdx = array(1, dim = c(nr, NY, ns, 1)), srv_idx_type = "abd",
      ObsSrvAgeComps = array(1 / NAG, dim = c(nr, NY, ns, NAG, nx, 1)),
      UseSrvAgeComps = array(1, dim = c(nr, NY, ns, 1)),
      ISS_SrvAgeComps = array(100, dim = c(nr, NY, ns, nx, 1)),
      SrvAgeComps_LikeType = like,
      ObsSrvLenComps = array(0, dim = c(nr, NY, ns, 1, nx, 1)),
      UseSrvLenComps = array(0, dim = c(nr, NY, ns, 1)),
      ISS_SrvLenComps = array(0, dim = c(nr, NY, ns, nx, 1)),
      SrvLenComps_LikeType = "none"),
    wt = list(Wt_FishAgeComps = array(1, dim = c(nr, NY, ns, nx, nf)),
              Wt_FishLenComps = array(1, dim = c(nr, NY, ns, nx, nf)),
              Wt_SrvAgeComps = array(1, dim = c(nr, NY, ns, nx, 1)),
              Wt_SrvLenComps = array(1, dim = c(nr, NY, ns, nx, 1))))

  fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)
}


test_that("the analytic gradient agrees with a finite difference of the objective", {
  problems <- character()

  for(nm in names(gradient_configs())) {
    obj <- tryCatch(gradient_obj(gradient_configs()[[nm]]), error = function(e) e)
    if(inherits(obj, "condition")) {
      problems <- c(problems, sprintf("%s: did not build (%s)", nm,
                                      substr(conditionMessage(obj), 1, 120)))
      next
    }

    chk <- fd_gradient_check(obj)
    # relative to the analytic value, with an absolute floor so a coordinate whose
    # gradient is near zero is not judged on a ratio of two roundoff errors
    rel <- abs(chk$analytic - chk$numeric) / pmax(abs(chk$analytic), 1)
    bad <- which(rel > 1e-4)
    for(i in bad) {
      problems <- c(problems, sprintf("%s: coordinate %d has analytic gradient %.8g against finite difference %.8g",
                                      nm, chk$index[i], chk$analytic[i], chk$numeric[i]))
    }
  }

  expect_equal(problems, character(0))
})


test_that("every configuration has a finite objective and a finite gradient", {
  # An infinite or missing gradient is how a log of zero, a division by an empty
  # sum, or an out-of-range index reaches the optimizer, which then either stalls
  # or steps somewhere arbitrary.
  problems <- character()

  for(nm in names(gradient_configs())) {
    obj <- tryCatch(gradient_obj(gradient_configs()[[nm]]), error = function(e) e)
    if(inherits(obj, "condition")) next
    if(!is.finite(obj$fn(obj$par)))
      problems <- c(problems, sprintf("%s: objective is not finite", nm))
    if(!all(is.finite(obj$gr(obj$par))))
      problems <- c(problems, sprintf("%s: gradient has %d non-finite entries", nm,
                                      sum(!is.finite(obj$gr(obj$par)))))
  }

  expect_equal(problems, character(0))
})


test_that("the gradient check would notice a wrong derivative", {
  # The comparison above only means something if the tolerance is tight enough to
  # separate a correct gradient from an incorrect one. Perturbing the point the
  # finite difference is taken at has to break it.
  obj <- gradient_obj(gradient_configs()$single_region)
  p <- obj$par
  g <- as.numeric(obj$gr(p))
  i <- which.max(abs(g))

  step <- 1e-5 * max(abs(p[i]), 1)
  up <- p; up[i] <- up[i] + step
  dn <- p; dn[i] <- dn[i] - step
  fd_here <- (obj$fn(up) - obj$fn(dn)) / (2 * step)

  # the same difference taken a long way from the evaluation point is a gradient
  # of somewhere else, and must not pass the tolerance the test uses
  far <- p; far[i] <- far[i] + 0.5
  up2 <- far; up2[i] <- up2[i] + step
  dn2 <- far; dn2[i] <- dn2[i] - step
  fd_far <- (obj$fn(up2) - obj$fn(dn2)) / (2 * step)

  expect_lt(abs(g[i] - fd_here) / max(abs(g[i]), 1), 1e-4)
  expect_gt(abs(g[i] - fd_far) / max(abs(g[i]), 1), 1e-4)
})
