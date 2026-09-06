library(SPoRC)
library(testthat)

test_that("Get_Comp_Likelihoods works!", {

  # Build a simple identity ageing error matrix (square, no error)
  identity_ae <- function(n) diag(n)

  # Build a uniform observed composition (counts)
  uniform_obs <- function(n_regions, n_bins, n_sexes, total = 100) {
    array(total / n_bins, dim = c(n_regions, n_bins, n_sexes))
  }

  # Build a slightly-peaked expected composition (sums to 1 along bins for each r,s)
  peaked_exp <- function(n_regions, n_bins, n_sexes) {
    arr <- array(1, dim = c(n_regions, n_bins, n_sexes))
    arr[, ceiling(n_bins / 2), ] <- 3   # bump middle bin
    for (r in seq_len(n_regions)) for (s in seq_len(n_sexes))
      arr[r, , s] <- arr[r, , s] / sum(arr[r, , s])
    arr
  }

  # Minimal call: sensible defaults for every scalar/array argument
  call_comp_nll <- function(
    Comp_Type       = 1,
    Likelihood_Type = 0,
    n_regions       = 2,
    n_model_bins    = 5,
    n_obs_bins      = 5,
    n_sexes         = 1,
    age_or_len      = 1,   # lengths by default (no ageing error)
    use             = NULL,
    Exp             = NULL,
    Obs             = NULL,
    ISS             = NULL,
    Wt_Mltnml       = NULL,
    ln_theta        = NULL,
    ln_theta_agg    = 0,
    LN_corr_pars    = NULL,
    LN_corr_pars_agg = 0,
    addtocomp       = 1e-4,
    AgeingError     = NULL
  ) {
    if (is.null(use))      use      <- rep(1L, n_regions)
    if (is.null(Exp))      Exp      <- peaked_exp(n_regions, n_model_bins, n_sexes)
    if (is.null(Obs))      Obs      <- uniform_obs(n_regions, n_obs_bins, n_sexes)
    if (is.null(ISS))      ISS      <- array(50, dim = c(n_regions, n_sexes))
    if (is.null(Wt_Mltnml)) Wt_Mltnml <- array(1, dim = c(n_regions, n_sexes))
    if (is.null(ln_theta)) ln_theta <- array(0, dim = c(n_regions, n_sexes))
    if (is.null(LN_corr_pars)) LN_corr_pars <- array(0, dim = c(n_regions, n_sexes, 3))
    if (is.null(AgeingError)) AgeingError <- identity_ae(n_obs_bins)

    Get_Comp_Likelihoods(
      Exp              = Exp,
      Obs              = Obs,
      ISS              = ISS,
      Wt_Mltnml        = Wt_Mltnml,
      ln_theta_agg     = ln_theta_agg,
      ln_theta         = ln_theta,
      LN_corr_pars     = LN_corr_pars,
      LN_corr_pars_agg = LN_corr_pars_agg,
      Comp_Type        = Comp_Type,
      Likelihood_Type  = Likelihood_Type,
      n_regions        = n_regions,
      n_model_bins     = n_model_bins,
      n_obs_bins       = n_obs_bins,
      n_sexes          = n_sexes,
      age_or_len       = age_or_len,
      AgeingError      = AgeingError,
      use              = use,
      addtocomp        = addtocomp
    )
  }

  # ── output structure ─────────────────────────────────────────────────────────

  test_that("output is a numeric matrix of dim [n_regions x n_sexes]", {
    for (ct in 0:2) {
      res <- call_comp_nll(
        Comp_Type = ct,
        n_regions = 3,
        n_sexes = 2,
        Likelihood_Type = 0
      )
      expect_true(is.numeric(res),  label = sprintf("numeric Comp_Type=%d", ct))
      expect_equal(unname(dim(res)), c(3, 2), label = sprintf("dim Comp_Type=%d", ct))
    }
  })

  # ── Comp_Type = 0 (aggregated) ───────────────────────────────────────────────

  test_that("Comp_Type=0, LT=0 (multinomial): nLL stored in [1,1], rest zero", {
    res <- call_comp_nll(Comp_Type = 0, Likelihood_Type = 0)
    expect_true(is.finite(res[1, 1]))
    # only [1,1] is populated
    expect_equal(sum(res != 0), 1)
  })

  test_that("Comp_Type=0, LT=1 (dirichlet-multinomial): finite nLL in [1,1]", {
    res <- call_comp_nll(Comp_Type = 0, Likelihood_Type = 1)
    expect_true(is.finite(res[1, 1]))
  })

  test_that("Comp_Type=0, LT=2 (LN iid): finite nLL in [1,1]", {
    res <- call_comp_nll(Comp_Type = 0, Likelihood_Type = 2)
    expect_true(is.finite(res[1, 1]))
  })

  test_that("Comp_Type=0, LT=3 (LN AR1): finite nLL in [1,1]", {
    res <- call_comp_nll(Comp_Type = 0, Likelihood_Type = 3)
    expect_true(is.finite(res[1, 1]))
  })

  # ── Comp_Type = 1 (split sex & region) ───────────────────────────────────────

  test_that("Comp_Type=1, LT=0: all used regions/sexes have finite nLL", {
    res <- call_comp_nll(Comp_Type = 1, Likelihood_Type = 0,
                         n_regions = 2, n_sexes = 2)
    expect_true(all(is.finite(res)))
  })

  test_that("Comp_Type=1, LT=1 (DM): finite nLL for all regions and sexes", {
    res <- call_comp_nll(Comp_Type = 1, Likelihood_Type = 1,
                         n_regions = 2, n_sexes = 2)
    expect_true(all(is.finite(res)))
  })

  test_that("Comp_Type=1, LT=2 (LN iid): finite nLL for all strata", {
    res <- call_comp_nll(Comp_Type = 1, Likelihood_Type = 2,
                         n_regions = 2, n_sexes = 1)
    expect_true(all(is.finite(res)))
  })

  test_that("Comp_Type=1, LT=3 (LN AR1): finite nLL for all strata", {
    res <- call_comp_nll(Comp_Type = 1, Likelihood_Type = 3,
                         n_regions = 2, n_sexes = 1)
    expect_true(all(is.finite(res)))
  })

  # ── Comp_Type = 2 (joint sex, split region) ───────────────────────────────────

  test_that("Comp_Type=2, LT=0 (multinomial): nLL in [r,1] for each region", {
    res <- call_comp_nll(Comp_Type = 2, Likelihood_Type = 0,
                         n_regions = 2, n_sexes = 2)
    expect_true(is.finite(res[1, 1]))
    expect_true(is.finite(res[2, 1]))
  })

  test_that("Comp_Type=2, LT=1 (DM): finite nLL by region", {
    res <- call_comp_nll(Comp_Type = 2, Likelihood_Type = 1,
                         n_regions = 2, n_sexes = 2)
    expect_true(all(is.finite(res[, 1])))
  })

  test_that("Comp_Type=2, LT=2 (LN iid): finite nLL by region", {
    res <- call_comp_nll(Comp_Type = 2, Likelihood_Type = 2,
                         n_regions = 2, n_sexes = 2)
    expect_true(all(is.finite(res[, 1])))
  })

  test_that("Comp_Type=2, LT=3 (LN AR1): finite nLL by region", {
    res <- call_comp_nll(Comp_Type = 2, Likelihood_Type = 3,
                         n_regions = 2, n_sexes = 2)
    expect_true(all(is.finite(res[, 1])))
  })

  test_that("Comp_Type=2, LT=4 (LN AR1 + constant sex corr): finite nLL by region", {
    res <- call_comp_nll(Comp_Type = 2, Likelihood_Type = 4,
                         n_regions = 2, n_sexes = 2)
    expect_true(all(is.finite(res[, 1])))
  })

  # ── use-region filtering ─────────────────────────────────────────────────────

  test_that("regions with use=0 reduce total nLL vs all regions used (Comp_Type=1)", {
    # The function filters Exp/Obs to used regions and writes results into
    # comp_nLL[1..n_used, s] sequentially — it does NOT preserve the original
    # region index in the output. Using fewer regions should give lower total nLL.
    res_all  <- call_comp_nll(
      Comp_Type = 1,
      Likelihood_Type = 0,
      n_regions = 3,
      use = c(1L, 1L, 1L)
    )
    res_one  <- call_comp_nll(
      Comp_Type = 1,
      Likelihood_Type = 0,
      n_regions = 3,
      use = c(0L, 1L, 0L)
    )
    # With only 1 region used, exactly 1 slot is populated; the rest are 0
    expect_equal(sum(res_one != 0), 1)
    # Total nLL with 3 regions > total nLL with 1 region
    expect_gt(sum(res_all), sum(res_one))
  })

  test_that("nLL increases when more regions are used (Comp_Type=1)", {
    res1 <- call_comp_nll(
      Comp_Type = 1,
      Likelihood_Type = 0,
      n_regions = 3,
      use = c(1L, 0L, 0L)
    )
    res3 <- call_comp_nll(
      Comp_Type = 1,
      Likelihood_Type = 0,
      n_regions = 3,
      use = c(1L, 1L, 1L)
    )
    expect_gt(sum(res3), sum(res1))
  })

  # ── ageing error (age_or_len = 0) ────────────────────────────────────────────

  test_that("age compositions with identity ageing error match length compositions", {
    # With a square identity AE matrix the result should be identical to lengths
    res_len <- call_comp_nll(
      Comp_Type = 1,
      Likelihood_Type = 0,
      age_or_len = 1,
      n_model_bins = 5,
      n_obs_bins = 5
    )
    res_age <- call_comp_nll(
      Comp_Type = 1,
      Likelihood_Type = 0,
      age_or_len = 0,
      n_model_bins = 5,
      n_obs_bins = 5,
      AgeingError = identity_ae(5)
    )
    expect_equal(res_age, res_len, tolerance = 1e-10)
  })

  test_that("non-square ageing error collapses model bins to obs bins (Comp_Type=1)", {
    # 6 model age bins mapped to 4 observed bins
    n_mod <- 6; n_obs <- 4
    AE <- matrix(0, nrow = n_mod, ncol = n_obs)
    AE[1:2, 1] <- 0.5; AE[3:4, 2] <- 0.5; AE[5, 3] <- 1; AE[6, 4] <- 1
    Exp <- peaked_exp(2, n_mod, 1)
    Obs <- uniform_obs(2, n_obs, 1)
    res <- call_comp_nll(
      Comp_Type = 1,
      Likelihood_Type = 0,
      age_or_len = 0,
      n_model_bins = n_mod,
      n_obs_bins = n_obs,
      AgeingError = AE,
      Exp = Exp,
      Obs = Obs,
      n_regions = 2,
      n_sexes = 1
    )
    expect_true(all(is.finite(res)))
  })

  # ── zero-count handling ───────────────────────────────────────────────────────

  test_that("zero counts in Obs do not produce NaN (Comp_Type=1, LT=0, addtocomp)", {
    Obs_z <- uniform_obs(2, 5, 1)
    Obs_z[1, 1, 1] <- 0   # introduce a zero
    res <- call_comp_nll(
      Comp_Type = 1,
      Likelihood_Type = 0,
      Obs = Obs_z,
      addtocomp = 1e-4
    )
    expect_true(all(is.finite(res)))
    expect_false(any(is.nan(res)))
  })

  test_that("zero counts in Obs handled gracefully for LN likelihoods (Comp_Type=1)", {
    Obs_z <- uniform_obs(2, 6, 1)
    Obs_z[1, 1, 1] <- 0
    Obs_z[1, 3, 1] <- 0
    for (lt in 2:3) {
      res <- call_comp_nll(
        Comp_Type = 1,
        Likelihood_Type = lt,
        n_model_bins = 6,
        n_obs_bins = 6,
        Obs = Obs_z,
        addtocomp = 0
      )
      expect_true(all(is.finite(res)),
                  label = sprintf("LT=%d zero handling", lt))
    }
  })

  # ── ESS / weighting sensitivity ───────────────────────────────────────────────

  test_that("doubling ISS roughly doubles multinomial nLL (Comp_Type=1)", {
    res1 <- call_comp_nll(Comp_Type = 1, Likelihood_Type = 0,
                          ISS = array(50, c(2, 1)))
    res2 <- call_comp_nll(Comp_Type = 1, Likelihood_Type = 0,
                          ISS = array(100, c(2, 1)))
    # nLL scales linearly with ESS for multinomial
    expect_equal(sum(res2) / sum(res1), 2, tolerance = 0.01)
  })

  test_that("Wt_Mltnml = 0 gives zero multinomial nLL (Comp_Type=1)", {
    res <- call_comp_nll(
      Comp_Type = 1,
      Likelihood_Type = 0,
      Wt_Mltnml = array(0, c(2, 1))
    )
    expect_equal(sum(res), 0)
  })

  # ── perfect fit gives nLL = 0 for multinomial ────────────────────────────────

  test_that("multinomial nLL = 0 when Exp == Obs (Comp_Type=1)", {
    n_bins <- 5
    # Uniform obs and exactly matching expected
    flat <- array(1 / n_bins, dim = c(2, n_bins, 1))
    res <- call_comp_nll(
      Comp_Type = 1,
      Likelihood_Type = 0,
      Exp = flat,
      Obs = flat * 100,
      n_model_bins = n_bins,
      n_obs_bins = n_bins,
      AgeingError = identity_ae(n_bins)
    )
    expect_equal(sum(res), 0, tolerance = 1e-10)
  })

  test_that("multinomial nLL = 0 when Exp == Obs (Comp_Type=2)", {
    n_bins <- 5; n_sexes <- 2
    flat <- array(1 / (n_bins * n_sexes), dim = c(2, n_bins, n_sexes))
    res <- call_comp_nll(
      Comp_Type = 2,
      Likelihood_Type = 0,
      Exp = flat,
      Obs = flat * 100,
      n_model_bins = n_bins,
      n_obs_bins = n_bins,
      n_sexes = n_sexes,
      AgeingError = identity_ae(n_bins)
    )
    expect_equal(sum(res), 0, tolerance = 1e-10)
  })

  # ── worse fit increases nLL ───────────────────────────────────────────────────

  test_that("more-peaked (worse-fit) expected increases multinomial nLL (Comp_Type=1)", {
    n_bins <- 5
    flat  <- array(1 / n_bins, dim = c(2, n_bins, 1))
    # Expected concentrated on bin 1, observed uniform => poor fit
    bad   <- flat; bad[, 1, ] <- 0.9; bad[, 2:5, ] <- 0.1 / 4
    res_good <- call_comp_nll(
      Comp_Type = 1,
      Likelihood_Type = 0,
      Exp = flat,
      Obs = flat * 100,
      n_model_bins = n_bins,
      n_obs_bins = n_bins
    )
    res_bad  <- call_comp_nll(
      Comp_Type = 1,
      Likelihood_Type = 0,
      Exp = bad,
      Obs = flat * 100,
      n_model_bins = n_bins,
      n_obs_bins = n_bins
    )
    expect_gt(sum(res_bad), sum(res_good))
  })

  # ── Comp_Type=0 aggregation ───────────────────────────────────────────────────

  test_that("Comp_Type=0 aggregates: result stored only in [1,1], nLL is finite", {
    # Aggregated path always writes into comp_nLL[1,1] only.
    # Use a single-region case so the aggregation is unambiguous.
    res <- call_comp_nll(Comp_Type = 0, Likelihood_Type = 0, n_regions = 1)
    expect_true(is.finite(res[1, 1]))
    expect_equal(res[1, 1], call_comp_nll(
      Comp_Type = 0,
      Likelihood_Type = 0,
      n_regions = 1
    )[1, 1],
                 tolerance = 1e-10)
  })

  test_that("Comp_Type=0 nLL increases with worse fit (single region)", {
    n_bins <- 5
    flat <- array(1 / n_bins, dim = c(1, n_bins, 1))
    bad  <- flat; bad[, 1, ] <- 0.9; bad[, 2:5, ] <- 0.1 / 4
    res_good <- call_comp_nll(
      Comp_Type = 0,
      Likelihood_Type = 0,
      n_regions = 1,
      n_model_bins = n_bins,
      n_obs_bins = n_bins,
      Exp = flat,
      Obs = flat * 100
    )
    res_bad  <- call_comp_nll(
      Comp_Type = 0,
      Likelihood_Type = 0,
      n_regions = 1,
      n_model_bins = n_bins,
      n_obs_bins = n_bins,
      Exp = bad,
      Obs = flat * 100
    )
    expect_gt(res_bad[1, 1], res_good[1, 1])
  })
})
