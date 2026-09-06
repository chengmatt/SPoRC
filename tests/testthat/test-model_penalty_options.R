library(SPoRC)
library(testthat)

test_that("get_selex_fixed_penalty centers a named set of selectivity parameters", {

  # [region, par, block, sex, fleet]
  pars <- array(0, dim = c(1, 12, 1, 1, 2))
  set.seed(41)
  pars[1,,1,1,1] <- rnorm(12, 0, 0.4)
  pars[1,,1,1,2] <- rnorm(12, 0.3, 0.4)

  test_that("one row reproduces wt * log(mean(exp(pars)))^2 over the named set", {
    tbl <- data.frame(region = 1, fleet = 1, block = 1, sex = 1, wt = 10)
    tbl$par <- list(1:7)
    expect_equal(SPoRC:::get_selex_fixed_penalty(tbl, pars),
                 10 * log(mean(exp(pars[1,1:7,1,1,1])))^2, tolerance = 1e-12)
  })

  test_that("the set is penalized jointly, not parameter by parameter", {
    joint <- data.frame(region = 1, fleet = 1, block = 1, sex = 1, wt = 1)
    joint$par <- list(1:4)
    separate <- data.frame(region = rep(1, 4), fleet = 1, block = 1, sex = 1, wt = 1)
    separate$par <- as.list(1:4)
    expect_false(isTRUE(all.equal(SPoRC:::get_selex_fixed_penalty(joint, pars),
                                  SPoRC:::get_selex_fixed_penalty(separate, pars))))
  })

  test_that("a set whose exponentials average to one contributes nothing", {
    flat <- array(0, dim = c(1, 3, 1, 1, 1)) # exp(0) = 1, so mean is 1 and log is 0
    tbl <- data.frame(region = 1, fleet = 1, block = 1, sex = 1, wt = 100)
    tbl$par <- list(1:3)
    expect_equal(SPoRC:::get_selex_fixed_penalty(tbl, flat), 0, tolerance = 1e-12)
  })

  test_that("rows accumulate and address the right fleet", {
    tbl <- data.frame(region = c(1, 1), fleet = c(1, 2), block = 1, sex = 1, wt = c(10, 3))
    tbl$par <- list(1:7, 1:5)
    expect_equal(SPoRC:::get_selex_fixed_penalty(tbl, pars),
                 10 * log(mean(exp(pars[1,1:7,1,1,1])))^2 + 3 * log(mean(exp(pars[1,1:5,1,1,2])))^2,
                 tolerance = 1e-12)
  })

  test_that("it is differentiable through an RTMB tape", {
    f <- function(p) {
      "c" <- RTMB::ADoverload("c")
      "[<-" <- RTMB::ADoverload("[<-")
      arr <- array(0, dim = c(1, 3, 1, 1, 1))
      for(i in 1:3) arr[1,i,1,1,1] <- p[i]
      tbl <- data.frame(region = 1, fleet = 1, block = 1, sex = 1, wt = 10)
      tbl$par <- list(1:3)
      SPoRC:::get_selex_fixed_penalty(tbl, arr)
    }
    obj <- RTMB::MakeADFun(f, c(0.2, -0.1, 0.4), silent = TRUE)
    expect_no_error(obj$fn(obj$par))
    expect_no_error(obj$gr(obj$par))
    fd <- sapply(1:3, function(i) {
      up <- dn <- obj$par; up[i] <- up[i] + 1e-6; dn[i] <- dn[i] - 1e-6
      (obj$fn(up) - obj$fn(dn)) / 2e-6
    })
    expect_equal(as.numeric(obj$gr(obj$par)), fd, tolerance = 1e-5)
  })

})

test_that("validate_selex_penalty checks the table only when the flag is on", {

  good <- data.frame(region = 1, fleet = 1, block = 1, sex = 1, wt = 10)
  good$par <- list(1:3)

  test_that("a flag of zero passes the table through untouched", {
    expect_null(SPoRC:::validate_selex_penalty(NULL, 0, "x"))
  })

  test_that("a plain integer par column is normalized to a list", {
    plain <- data.frame(region = 1, fleet = 1, block = 1, sex = 1, par = 2, wt = 1)
    out <- SPoRC:::validate_selex_penalty(plain, 1, "x")
    expect_true(is.list(out$par))
    expect_equal(out$par[[1]], 2)
  })

  test_that("errors on a missing table, missing columns, a negative weight, and an empty set", {
    expect_error(SPoRC:::validate_selex_penalty(NULL, 1, "x"))
    expect_error(SPoRC:::validate_selex_penalty(good[, c("region", "fleet")], 1, "x"))
    bad_wt <- good; bad_wt$wt <- -1
    expect_error(SPoRC:::validate_selex_penalty(bad_wt, 1, "x"))
    empty <- good; empty$par <- list(integer(0))
    expect_error(SPoRC:::validate_selex_penalty(empty, 1, "x"))
  })

})

test_that("a deviation-specific Wt_Rec excludes years from the recruitment penalty", {

  n_pop <- 1; n_regions <- 1; n_ages <- 5; n_dev <- 6
  ln_RecDevs <- array(0, dim = c(n_pop, n_regions, n_dev))
  set.seed(42)
  ln_RecDevs[1,1,] <- rnorm(n_dev, 0, 0.6)
  ln_InitDevs <- array(rnorm(n_ages - 1, 0, 0.5), dim = c(n_pop, n_regions, n_ages - 1))

  run <- function(map = NULL) {
    SPoRC:::get_recruitment_penalty(
      n_pop = n_pop,
      n_regions = n_regions,
      n_ages = n_ages,
      n_est_rec_devs = n_dev,
      rec_dd = 0,
      natal_region = 1,
      rec_region_prop_spec = 0,
      rec_region_prop = array(1, dim = c(n_pop, n_regions)),
      equil_init_age_strc = 2,
      ln_InitDevs = ln_InitDevs,
      init_age_devs_shared = NULL,
      ln_sigmaR = array(log(0.6), dim = c(2, n_pop, n_regions)),
      bias_ramp = rep(0, n_dev),
      sigmaR_switch = 1,
      ln_RecDevs = ln_RecDevs,
      sigmaR2_early = array(0.36, dim = c(n_pop, n_regions)),
      sigmaR2_late = array(0.36, dim = c(n_pop, n_regions)),
      do_rec_bias_ramp = 0,
      map_ln_RecDevs = map
    )
  }

  pen <- run()

  test_that("the penalty is per deviation, so a weight array selects years", {
    wt <- array(1, dim = dim(ln_RecDevs))
    wt[1,1,c(2, 5)] <- 0
    expect_equal(sum(wt * pen$Rec_nLL), sum(pen$Rec_nLL[1,1,-c(2,5)]), tolerance = 1e-12)
  })

  test_that("a scalar weight still reproduces the old whole-penalty scaling", {
    expect_equal(sum(3 * pen$Rec_nLL), 3 * sum(pen$Rec_nLL), tolerance = 1e-12)
  })

  test_that("zero weight and an unestimated deviation reach the same objective by different routes", {
    map <- array(seq_len(n_dev), dim = dim(ln_RecDevs))
    map[1,1,3] <- NA
    wt <- array(1, dim = dim(ln_RecDevs)); wt[1,1,3] <- 0
    expect_equal(sum(run(map)$Rec_nLL), sum(wt * pen$Rec_nLL), tolerance = 1e-12)
  })

  test_that("the initial age penalty is separately dimensioned from the recruitment one", {
    expect_equal(dim(pen$Init_Rec_nLL)[3], n_ages - 1)
    expect_equal(dim(pen$Rec_nLL)[3], n_dev)
    expect_false(identical(dim(pen$Init_Rec_nLL), dim(pen$Rec_nLL)))
  })

})

test_that("Setup_Mod_Weighting validates the recruitment weight shapes", {

  base_list <- Setup_Mod_Dim(
    years = 1:10,
    ages = 1:5,
    lens = NA,
    n_regions = 1,
    n_sexes = 1,
    n_fish_fleets = 1,
    n_srv_fleets = 1,
    n_seas = 1,
    n_pop = 1,
    natal_region = 1
  )
  base_list <- Setup_Mod_Rec(input_list = base_list, rec_model = "mean_rec")

  rec_dim <- dim(base_list$par$ln_RecDevs)
  init_dim <- dim(base_list$par$ln_InitDevs)

  test_that("a scalar Wt_Rec fills in Wt_Init_Rec", {
    out <- Setup_Mod_Weighting(input_list = base_list, Wt_Rec = 2)
    expect_equal(out$data$Wt_Rec, 2)
    expect_equal(out$data$Wt_Init_Rec, 2)
  })

  test_that("an array Wt_Rec is accepted at the deviation array's shape", {
    wt <- array(1, dim = rec_dim); wt[1,1,rec_dim[3]] <- 0
    out <- Setup_Mod_Weighting(
      input_list = base_list,
      Wt_Rec = wt,
      Wt_Init_Rec = array(1, dim = init_dim)
    )
    expect_equal(dim(out$data$Wt_Rec), rec_dim)
  })

  test_that("an array Wt_Rec without an explicit Wt_Init_Rec is refused, since the two differ in shape", {
    expect_error(Setup_Mod_Weighting(input_list = base_list, Wt_Rec = array(1, dim = rec_dim)))
  })

  test_that("a wrongly shaped weight is refused rather than silently recycled", {
    expect_error(Setup_Mod_Weighting(
      input_list = base_list,
      Wt_Rec = array(1, dim = c(1, 1, rec_dim[3] + 3)),
      Wt_Init_Rec = 1
    ))
    expect_error(Setup_Mod_Weighting(
      input_list = base_list,
      Wt_Rec = 1,
      Wt_Init_Rec = array(1, dim = c(1, 1, init_dim[3] + 2))
    ))
  })

})
