library(SPoRC)
library(testthat)

test_that("nonparlog is non-parametric on the log scale and centered within each year", {

  n_ages <- 8
  set.seed(51)
  pars <- rnorm(n_ages, 0, 0.5)
  ln_seldevs <- array(0, dim = c(1, 3, n_ages, 1, 1))
  ln_seldevs[1, 2, , 1, 1] <- rnorm(n_ages, 0, 0.3)
  ln_seldevs[1, 3, , 1, 1] <- rnorm(n_ages, 0, 0.3)

  sel <- function(y, tv = 2) Get_Selex(
    Selex_Model = 9,
    TimeVary_Model = tv,
    pars = pars,
    ln_seldevs = ln_seldevs,
    Region = 1,
    Year = y,
    Bin = 1:n_ages,
    Sex = 1
  )

  test_that("every year averages to one across bins", {
    for(y in 1:3) expect_equal(mean(sel(y)), 1, tolerance = 1e-12)
  })

  test_that("it is exp(par + dev) rescaled by that year's own mean", {
    raw <- exp(pars + ln_seldevs[1, 2, , 1, 1])
    expect_equal(sel(2), raw / mean(raw), tolerance = 1e-12)
  })

  test_that("the level of the parameters is not identified, only their differences", {
    shifted <- Get_Selex(
      Selex_Model = 9,
      TimeVary_Model = 2,
      pars = pars + 3.7,
      ln_seldevs = ln_seldevs,
      Region = 1,
      Year = 2,
      Bin = 1:n_ages,
      Sex = 1
    )
    expect_equal(shifted, sel(2), tolerance = 1e-12)
  })

  test_that("without time variation every year is the same curve", {
    expect_equal(sel(1, tv = 0), sel(3, tv = 0), tolerance = 1e-12)
  })

  test_that("it differs from the logit-scale non-parametric form", {
    logit_form <- Get_Selex(
      Selex_Model = 5,
      TimeVary_Model = 2,
      pars = pars,
      ln_seldevs = ln_seldevs,
      Region = 1,
      Year = 2,
      Bin = 1:n_ages,
      Sex = 1
    )
    expect_false(isTRUE(all.equal(as.numeric(sel(2)), as.numeric(logit_form / mean(logit_form)))))
  })

  test_that("it is differentiable through an RTMB tape", {
    f <- function(p) {
      "c" <- RTMB::ADoverload("c")
      "[<-" <- RTMB::ADoverload("[<-")
      s <- Get_Selex(
        Selex_Model = 9,
        TimeVary_Model = 0,
        pars = p,
        ln_seldevs = ln_seldevs,
        Region = 1,
        Year = 1,
        Bin = 1:3,
        Sex = 1
      )
      sum(s^2)
    }
    obj <- RTMB::MakeADFun(f, c(0.1, 0.4, -0.2), silent = TRUE)
    expect_no_error(obj$fn(obj$par))
    expect_no_error(obj$gr(obj$par))
  })

})

test_that("bin overrides replace named bins with their own free deviations", {

  n_ages <- 6; n_yrs <- 4
  ln_seldevs <- array(0, dim = c(1, n_yrs, n_ages, 1, 1))
  bin_devs <- array(0, dim = c(1, n_yrs, n_ages, 1, 1))
  set.seed(52)
  for(y in 1:n_yrs) bin_devs[1, y, , 1, 1] <- rnorm(n_ages, -2, 0.4)

  base <- function(bins = NULL) t(sapply(1:n_yrs, function(y)
    Get_Selex(
      Selex_Model = 0,
      TimeVary_Model = 0,
      pars = c(log(3), log(1.2)),
      ln_seldevs = ln_seldevs,
      Region = 1,
      Year = y,
      Bin = 1:n_ages,
      Sex = 1,
      bin_devs = bin_devs,
      bin_dev_bins = bins
    )))

  plain <- base(NULL)
  overridden <- base(c(1, 2))

  test_that("named bins take exp(bin_dev) exactly", {
    for(y in 1:n_yrs) expect_equal(overridden[y, 1:2], exp(bin_devs[1, y, 1:2, 1, 1]), tolerance = 1e-12)
  })

  test_that("bins that are not named keep the parametric value", {
    expect_equal(overridden[, 3:n_ages], plain[, 3:n_ages], tolerance = 1e-12)
  })

  test_that("no override at all leaves the curve untouched", {
    expect_equal(base(integer(0)), plain, tolerance = 1e-12)
  })

  test_that("an overridden bin varies by year even when the form does not", {
    expect_equal(plain[1, 1], plain[4, 1], tolerance = 1e-12)
    expect_false(isTRUE(all.equal(overridden[1, 1], overridden[4, 1])))
  })

  test_that("the override survives a form that standardizes internally", {
    with_ovr <- Get_Selex(
      Selex_Model = 9,
      TimeVary_Model = 0,
      pars = rep(0, n_ages),
      ln_seldevs = ln_seldevs,
      Region = 1,
      Year = 1,
      Bin = 1:n_ages,
      Sex = 1,
      bin_devs = bin_devs,
      bin_dev_bins = 1
    )
    expect_equal(with_ovr[1], exp(bin_devs[1, 1, 1, 1, 1]), tolerance = 1e-12)
  })

  test_that("it is differentiable through an RTMB tape", {
    f <- function(p) {
      "c" <- RTMB::ADoverload("c")
      "[<-" <- RTMB::ADoverload("[<-")
      bd <- array(0, dim = c(1, 1, 3, 1, 1))
      for(i in 1:3) bd[1, 1, i, 1, 1] <- p[i]
      s <- Get_Selex(
        Selex_Model = 0,
        TimeVary_Model = 0,
        pars = c(log(2), log(1)),
        ln_seldevs = array(0, dim = c(1, 1, 3, 1, 1)),
        Region = 1,
        Year = 1,
        Bin = 1:3,
        Sex = 1,
        bin_devs = bd,
        bin_dev_bins = c(1, 2)
      )
      sum(s^2)
    }
    obj <- RTMB::MakeADFun(f, c(-1, -0.5, 0.2), silent = TRUE)
    expect_no_error(obj$fn(obj$par))
    # only the two overridden bins are on the tape, so the third gradient is zero
    expect_equal(as.numeric(obj$gr(obj$par))[3], 0, tolerance = 1e-12)
  })

})

test_that("setup_sel_bin_devs builds the parameter, map and process error arrays", {

  input_list <- Setup_Mod_Dim(
    years = 1:12,
    ages = 1:7,
    lens = NA,
    n_regions = 1,
    n_sexes = 1,
    n_fish_fleets = 1,
    n_srv_fleets = 2,
    n_seas = 1,
    n_pop = 1,
    natal_region = 1,
    verbose = FALSE
  )
  input_list$data$n_proj_yrs_devs <- 0

  out <- SPoRC:::setup_sel_bin_devs(
    input_list,
    list(1, c(1, 2)),
    c("none", "rw"),
    prefix = "srv",
    n_fleets = 2,
    bins = 7
  )

  test_that("the deviation array spans regions, years, bins, sexes and fleets", {
    expect_equal(dim(out$par$ln_srvsel_bin_devs), c(1, 12, 7, 1, 2))
    expect_equal(dim(out$par$srvsel_bin_devs_pe_pars), c(1, 7, 1, 2))
  })

  test_that("only the named bins are estimated", {
    map <- out$data$map_ln_srvsel_bin_devs
    expect_equal(sum(!is.na(map[1, , , 1, 1])), 12)      # fleet 1 overrides one bin
    expect_equal(sum(!is.na(map[1, , , 1, 2])), 24)      # fleet 2 overrides two bins
    expect_true(all(is.na(map[1, , 3:7, 1, 1])))
  })

  test_that("a process error hyperparameter exists only where a structure was asked for", {
    pe <- array(as.numeric(levels(out$map$srvsel_bin_devs_pe_pars))[out$map$srvsel_bin_devs_pe_pars],
                dim = dim(out$par$srvsel_bin_devs_pe_pars))
    expect_true(all(is.na(pe[1, , 1, 1])))               # fleet 1 has no structure
    expect_equal(sum(!is.na(pe[1, , 1, 2])), 2)          # fleet 2 has one per overridden bin
  })

  test_that("the override bins are recorded as a 0/1 array", {
    expect_equal(out$data$srv_sel_bin_dev_bins[, 1], c(1, 0, 0, 0, 0, 0, 0))
    expect_equal(out$data$srv_sel_bin_dev_bins[, 2], c(1, 1, 0, 0, 0, 0, 0))
  })

  test_that("errors on an out-of-range bin, a wrong-length list, and an unknown structure", {
    expect_error(SPoRC:::setup_sel_bin_devs(input_list, list(1, 99), c("none", "rw"), "srv", 2, 7))
    expect_error(SPoRC:::setup_sel_bin_devs(input_list, list(1), c("none", "rw"), "srv", 2, 7))
    expect_error(SPoRC:::setup_sel_bin_devs(input_list, list(1, 2), c("none", "bogus"), "srv", 2, 7))
  })

})
