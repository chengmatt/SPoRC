library(SPoRC)
library(testthat)

test_that("OSA composition likelihood pipeline works!", {

  identity_ae <- function(n) diag(n)
  fam_of <- function(lt) if (lt %in% c(0, 1)) "discrete" else "continuous"
  alr <- function(p) log(p[-length(p)]) - log(p[length(p)])

  # normalize a [region, bin, sex] array so each (r,s) slice sums to 1
  normalize_rs <- function(arr) {
    d <- dim(arr)
    for (r in seq_len(d[1])) for (s in seq_len(d[3])) {
      arr[r, , s] <- arr[r, , s] / sum(arr[r, , s])
    }
    arr
  }

  # Build a [n_regions, length(vec), n_sexes] array where every (r,s) slice
  # is exactly `vec`. NOTE: array(rep(vec, n_regions), dim = c(n_regions, ...))
  # does NOT do this -- R arrays fill column-major (region fastest), so
  # rep()-ing the vector and reshaping scrambles which bin goes to which
  # region. Build it explicitly instead.
  make_prop_array <- function(vec, n_regions, n_sexes) {
    arr <- array(0, dim = c(n_regions, length(vec), n_sexes))
    for (r in seq_len(n_regions)) for (s in seq_len(n_sexes)) arr[r, , s] <- vec
    arr
  }

  # Pack a *single* (year=1, season=1, fleet=1, pop=1) group and return the
  # tracked OSA vector, using the real pack_comp_osa() machinery so that the
  # unit tests below stay in lock-step with the packer's actual conventions.
  pack_single_group <- function(prop_or_counts, ISS, Wt, use,
                                n_obs_bins, n_sexes, comp_type, like_type,
                                addtocomp = 1e-4) {
    n_regions <- length(use)
    ObsArr <- array(prop_or_counts, dim = c(n_regions, 1, 1, n_obs_bins, n_sexes, 1))
    ISSArr <- array(ISS, dim = c(n_regions, 1, 1, n_sexes, 1))
    WtArr  <- array(Wt,  dim = c(n_regions, 1, 1, n_sexes, 1))
    UseArr <- array(use, dim = c(n_regions, 1, 1, 1))
    TypeMat <- matrix(comp_type, nrow = 1, ncol = 1)
    LikeTypeVec <- like_type

    pack_comp_osa(
      ObsArr = ObsArr, ISSArr = ISSArr, WtArr = WtArr, UseArr = UseArr,
      TypeMat = TypeMat, LikeTypeVec = LikeTypeVec,
      n_yrs = 1, n_seas = 1, n_fleets = 1, n_sexes = n_sexes,
      addtocomp = addtocomp, family = fam_of(like_type),
      pop = FALSE, n_pop = 1
    )
  }

  # Pack + evaluate Get_Comp_Likelihoods_OSA() for a single group, with
  # sensible defaults for every argument (mirrors call_comp_nll() in the
  # non-OSA test file).
  call_osa_nll <- function(
    Comp_Type, Likelihood_Type,
    n_regions    = 2,
    n_model_bins = 4,
    n_obs_bins   = 4,
    n_sexes      = 1,
    age_or_len   = 1,     # lengths by default (no ageing error needed)
    use          = NULL,
    Exp          = NULL,
    Obs_raw      = NULL,  # raw proportions/counts BEFORE packing
    ISS          = NULL,
    Wt           = NULL,
    ln_theta     = NULL,
    ln_theta_agg = 0,
    LN_corr_pars = NULL,
    LN_corr_pars_agg = 0,
    addtocomp    = 1e-4,
    AgeingError  = NULL
  ) {
    if (is.null(use)) use <- rep(1L, n_regions)
    if (is.null(Exp)) {
      Exp <- normalize_rs(array(1, dim = c(n_regions, n_model_bins, n_sexes)))
    }
    if (is.null(Obs_raw)) {
      Obs_raw <- array(100 / n_obs_bins, dim = c(n_regions, n_obs_bins, n_sexes))
    }
    if (is.null(ISS))      ISS      <- array(1000, dim = c(n_regions, n_sexes))
    if (is.null(Wt))       Wt       <- array(1,    dim = c(n_regions, n_sexes))
    if (is.null(ln_theta)) ln_theta <- array(0,    dim = c(n_regions, n_sexes))
    if (is.null(LN_corr_pars)) LN_corr_pars <- array(0, dim = c(n_regions, n_sexes, 3))
    if (is.null(AgeingError))  AgeingError  <- identity_ae(n_obs_bins)

    tracked <- pack_single_group(Obs_raw, ISS, Wt, use, n_obs_bins, n_sexes,
                                 Comp_Type, Likelihood_Type, addtocomp)

    Get_Comp_Likelihoods_OSA(
      Exp = Exp, Obs = tracked, ISS = ISS,
      ln_theta = ln_theta, ln_theta_agg = ln_theta_agg,
      LN_corr_pars = LN_corr_pars, LN_corr_pars_agg = LN_corr_pars_agg,
      Comp_Type = Comp_Type, Likelihood_Type = Likelihood_Type,
      n_regions = n_regions, n_model_bins = n_model_bins, n_obs_bins = n_obs_bins,
      n_sexes = n_sexes, age_or_len = age_or_len, AgeingError = AgeingError,
      use = use, addtocomp = addtocomp
    )
  }

  # Get_Comp_Likelihoods_OSA: structure
  test_that("OSA output is a numeric matrix of dim [n_regions x n_sexes]", {
    for (ct in 0:2) {
      lts <- if (ct == 2) 0:4 else 0:3
      for (lt in lts) {
        res <- call_osa_nll(Comp_Type = ct, Likelihood_Type = lt,
                            n_regions = 3, n_sexes = 2)
        expect_true(is.numeric(res),
                    label = sprintf("numeric Comp_Type=%d LT=%d", ct, lt))
        expect_equal(unname(dim(res)), c(3, 2),
                     label = sprintf("dim Comp_Type=%d LT=%d", ct, lt))
      }
    }
  })

  test_that("Comp_Type=0 (OSA): nLL stored only in [1,1]", {
    for (lt in 0:3) {
      res <- call_osa_nll(Comp_Type = 0, Likelihood_Type = lt, n_regions = 2, n_sexes = 2)
      expect_true(is.finite(res[1, 1]), label = sprintf("finite LT=%d", lt))
      expect_equal(sum(res != 0), 1, label = sprintf("only [1,1] populated LT=%d", lt))
    }
  })

  test_that("Comp_Type=1 (OSA): all used region/sex strata are finite", {
    for (lt in 0:3) {
      res <- call_osa_nll(Comp_Type = 1, Likelihood_Type = lt, n_regions = 2, n_sexes = 2)
      expect_true(all(is.finite(res)), label = sprintf("LT=%d", lt))
    }
  })

  test_that("Comp_Type=2 (OSA): nLL stored only in column 1, one entry per region", {
    for (lt in 0:4) {
      res <- call_osa_nll(Comp_Type = 2, Likelihood_Type = lt, n_regions = 2, n_sexes = 2)
      expect_true(all(is.finite(res[, 1])), label = sprintf("LT=%d", lt))
      expect_equal(sum(res[, 2] != 0), 0, label = sprintf("sex column 2 untouched LT=%d", lt))
    }
  })

  # Get_Comp_Likelihoods_OSA: fit-quality sanity checks
  test_that("multinomial (OSA) favors Exp matching the true generating proportions", {
    n_bins <- 5
    p_true <- c(0.05, 0.1, 0.6, 0.15, 0.1)
    p_bad  <- c(0.6, 0.1, 0.05, 0.1, 0.15)

    Exp_good <- normalize_rs(make_prop_array(p_true, 2, 1))
    Exp_bad  <- normalize_rs(make_prop_array(p_bad,  2, 1))
    Obs_raw  <- make_prop_array(p_true, 2, 1) * 100

    res_good <- call_osa_nll(Comp_Type = 1, Likelihood_Type = 0, n_regions = 2,
                             n_model_bins = n_bins, n_obs_bins = n_bins,
                             Exp = Exp_good, Obs_raw = Obs_raw)
    res_bad  <- call_osa_nll(Comp_Type = 1, Likelihood_Type = 0, n_regions = 2,
                             n_model_bins = n_bins, n_obs_bins = n_bins,
                             Exp = Exp_bad, Obs_raw = Obs_raw)

    expect_gt(sum(res_bad), sum(res_good))
  })

  test_that("Dirichlet-multinomial (OSA) favors Exp matching the true generating proportions", {
    n_bins <- 5
    p_true <- c(0.05, 0.1, 0.6, 0.15, 0.1)
    p_bad  <- c(0.6, 0.1, 0.05, 0.1, 0.15)

    Exp_good <- normalize_rs(make_prop_array(p_true, 2, 1))
    Exp_bad  <- normalize_rs(make_prop_array(p_bad,  2, 1))
    Obs_raw  <- make_prop_array(p_true, 2, 1) * 100

    res_good <- call_osa_nll(Comp_Type = 1, Likelihood_Type = 1, n_regions = 2,
                             n_model_bins = n_bins, n_obs_bins = n_bins,
                             Exp = Exp_good, Obs_raw = Obs_raw)
    res_bad  <- call_osa_nll(Comp_Type = 1, Likelihood_Type = 1, n_regions = 2,
                             n_model_bins = n_bins, n_obs_bins = n_bins,
                             Exp = Exp_bad, Obs_raw = Obs_raw)

    expect_gt(sum(res_bad), sum(res_good))
  })

  test_that("logistic-normal iid (OSA) favors Exp matching the true generating proportions", {
    n_bins <- 5
    p_true <- c(0.05, 0.1, 0.6, 0.15, 0.1)
    p_bad  <- c(0.6, 0.1, 0.05, 0.1, 0.15)

    Exp_good <- normalize_rs(make_prop_array(p_true, 2, 1))
    Exp_bad  <- normalize_rs(make_prop_array(p_bad,  2, 1))
    Obs_raw  <- make_prop_array(p_true, 2, 1)

    res_good <- call_osa_nll(Comp_Type = 1, Likelihood_Type = 2, n_regions = 2,
                             n_model_bins = n_bins, n_obs_bins = n_bins,
                             Exp = Exp_good, Obs_raw = Obs_raw)
    res_bad  <- call_osa_nll(Comp_Type = 1, Likelihood_Type = 2, n_regions = 2,
                             n_model_bins = n_bins, n_obs_bins = n_bins,
                             Exp = Exp_bad, Obs_raw = Obs_raw)

    expect_gt(sum(res_bad), sum(res_good))
  })

  test_that("logistic-normal AR1 (OSA) favors Exp matching the true generating proportions", {
    n_bins <- 5
    p_true <- c(0.05, 0.1, 0.6, 0.15, 0.1)
    p_bad  <- c(0.6, 0.1, 0.05, 0.1, 0.15)

    Exp_good <- normalize_rs(make_prop_array(p_true, 2, 1))
    Exp_bad  <- normalize_rs(make_prop_array(p_bad,  2, 1))
    Obs_raw  <- make_prop_array(p_true, 2, 1)

    res_good <- call_osa_nll(Comp_Type = 1, Likelihood_Type = 3, n_regions = 2,
                             n_model_bins = n_bins, n_obs_bins = n_bins,
                             Exp = Exp_good, Obs_raw = Obs_raw)
    res_bad  <- call_osa_nll(Comp_Type = 1, Likelihood_Type = 3, n_regions = 2,
                             n_model_bins = n_bins, n_obs_bins = n_bins,
                             Exp = Exp_bad, Obs_raw = Obs_raw)

    expect_gt(sum(res_bad), sum(res_good))
  })

  test_that("2D AR1 + constant sex correlation (OSA, Comp_Type=2, LT=4) is finite and fit-sensitive", {
    n_bins <- 4
    p_true <- c(0.1, 0.6, 0.2, 0.1)
    p_bad  <- c(0.6, 0.1, 0.1, 0.2)

    Exp_good <- normalize_rs(make_prop_array(p_true, 2, 2))
    Exp_bad  <- normalize_rs(make_prop_array(p_bad,  2, 2))
    Obs_raw  <- make_prop_array(p_true, 2, 2)

    res_good <- call_osa_nll(Comp_Type = 2, Likelihood_Type = 4, n_regions = 2, n_sexes = 2,
                             n_model_bins = n_bins, n_obs_bins = n_bins,
                             Exp = Exp_good, Obs_raw = Obs_raw)
    res_bad  <- call_osa_nll(Comp_Type = 2, Likelihood_Type = 4, n_regions = 2, n_sexes = 2,
                             n_model_bins = n_bins, n_obs_bins = n_bins,
                             Exp = Exp_bad, Obs_raw = Obs_raw)

    expect_true(all(is.finite(res_good[, 1])))
    expect_gt(sum(res_bad[, 1]), sum(res_good[, 1]))
  })

  test_that("zeros in the composition don't break LN (OSA) as long as addtocomp > 0", {
    n_bins <- 6
    Obs_raw <- array(1 / n_bins, dim = c(2, n_bins, 1))
    Obs_raw[1, 1, 1] <- 0
    Obs_raw[1, 3, 1] <- 0
    for (lt in 2:3) {
      res <- call_osa_nll(Comp_Type = 1, Likelihood_Type = lt, n_regions = 2,
                          n_model_bins = n_bins, n_obs_bins = n_bins,
                          Obs_raw = Obs_raw, addtocomp = 1e-4)
      expect_true(all(is.finite(res)), label = sprintf("LT=%d zero handling", lt))
    }
  })

  test_that("use-region filtering (OSA): fewer used regions -> fewer populated strata", {
    res_one <- call_osa_nll(Comp_Type = 1, Likelihood_Type = 0, n_regions = 3,
                            use = c(0L, 1L, 0L))
    expect_equal(sum(res_one != 0), 1)
  })

  # pack_comp_osa: shapes, NULL handling, numeric correctness
  test_that("pack_comp_osa returns NULL when no fleet matches the requested family", {
    n_regions <- 2; n_obs_bins <- 4; n_sexes <- 1
    ObsArr <- array(1, dim = c(n_regions, 1, 1, n_obs_bins, n_sexes, 1))
    ISSArr <- array(50, dim = c(n_regions, 1, 1, n_sexes, 1))
    WtArr  <- array(1,  dim = c(n_regions, 1, 1, n_sexes, 1))
    UseArr <- array(1L, dim = c(n_regions, 1, 1, 1))
    TypeMat <- matrix(1, nrow = 1, ncol = 1)
    LikeTypeVec <- 0   # multinomial only -> "discrete"

    out <- pack_comp_osa(ObsArr, ISSArr, WtArr, UseArr, TypeMat, LikeTypeVec,
                         n_yrs = 1, n_seas = 1, n_fleets = 1, n_sexes = n_sexes,
                         addtocomp = 1e-4, family = "continuous",
                         pop = FALSE, n_pop = 1)
    expect_null(out)
  })

  test_that("pack_comp_osa returns NULL when no region has use=1", {
    n_regions <- 2; n_obs_bins <- 4; n_sexes <- 1
    ObsArr <- array(1, dim = c(n_regions, 1, 1, n_obs_bins, n_sexes, 1))
    ISSArr <- array(50, dim = c(n_regions, 1, 1, n_sexes, 1))
    WtArr  <- array(1,  dim = c(n_regions, 1, 1, n_sexes, 1))
    UseArr <- array(0L, dim = c(n_regions, 1, 1, 1))
    TypeMat <- matrix(1, nrow = 1, ncol = 1)
    LikeTypeVec <- 0

    out <- pack_comp_osa(ObsArr, ISSArr, WtArr, UseArr, TypeMat, LikeTypeVec,
                         n_yrs = 1, n_seas = 1, n_fleets = 1, n_sexes = n_sexes,
                         addtocomp = 1e-4, family = "discrete",
                         pop = FALSE, n_pop = 1)
    expect_null(out)
  })

  test_that("pack_comp_osa produces the documented vector lengths (discrete)", {
    n_obs_bins <- 5; n_sexes <- 2; use <- c(1L, 1L, 0L)  # n_ru = 2
    n_ru <- sum(use)

    for (ct in 0:2) {
      g <- pack_single_group(
        prop_or_counts = array(1 / n_obs_bins, dim = c(3, n_obs_bins, n_sexes)),
        ISS = array(50, dim = c(3, n_sexes)), Wt = array(1, dim = c(3, n_sexes)),
        use = use, n_obs_bins = n_obs_bins, n_sexes = n_sexes,
        comp_type = ct, like_type = 0
      )
      expected_len <- if (ct == 0) n_obs_bins else n_ru * n_obs_bins * n_sexes
      expect_equal(length(g), expected_len, label = sprintf("discrete ct=%d", ct))
    }
  })

  test_that("pack_comp_osa produces the documented vector lengths (continuous / LN)", {
    n_obs_bins <- 5; n_sexes <- 2; use <- c(1L, 1L, 0L)
    n_ru <- sum(use)

    for (ct in 0:2) {
      g <- pack_single_group(
        prop_or_counts = array(1 / n_obs_bins, dim = c(3, n_obs_bins, n_sexes)),
        ISS = array(50, dim = c(3, n_sexes)), Wt = array(1, dim = c(3, n_sexes)),
        use = use, n_obs_bins = n_obs_bins, n_sexes = n_sexes,
        comp_type = ct, like_type = 2
      )
      expected_len <- if (ct == 0) n_obs_bins - 1
      else if (ct == 1) n_ru * (n_obs_bins - 1) * n_sexes
      else n_ru * (n_obs_bins * n_sexes - 1)
      expect_equal(length(g), expected_len, label = sprintf("continuous ct=%d", ct))
    }
  })

  test_that("pack_comp_osa (multinomial, ct=0) produces round(pr * ISS * Wt) counts", {
    n_obs_bins <- 4
    props <- c(0.1, 0.2, 0.3, 0.4)
    prop_arr <- make_prop_array(props, 2, 1)  # only region1/sex1 matters for ct=0
    iss <- 200; wt <- 0.5

    g <- pack_single_group(prop_arr, ISS = array(iss, dim = c(2, 1)),
                           Wt = array(wt, dim = c(2, 1)), use = c(1L, 1L),
                           n_obs_bins = n_obs_bins, n_sexes = 1,
                           comp_type = 0, like_type = 0, addtocomp = 1e-4)

    pr_expected <- (props + 1e-4) / sum(props + 1e-4)
    expected <- round(pr_expected * iss * wt)
    expect_equal(as.numeric(g), as.numeric(expected))
  })

  test_that("pack_comp_osa (Dirichlet-multinomial, ct=0) produces round(pr * ISS) counts (no Wt)", {
    n_obs_bins <- 4
    props <- c(0.1, 0.2, 0.3, 0.4)
    prop_arr <- make_prop_array(props, 2, 1)
    iss <- 150; wt <- 999  # Wt should be ignored for DM

    g <- pack_single_group(prop_arr, ISS = array(iss, dim = c(2, 1)),
                           Wt = array(wt, dim = c(2, 1)), use = c(1L, 1L),
                           n_obs_bins = n_obs_bins, n_sexes = 1,
                           comp_type = 0, like_type = 1, addtocomp = 1e-4)

    pr_expected <- (props + 1e-4) / sum(props + 1e-4)
    expected <- round(pr_expected * iss)
    expect_equal(as.numeric(g), as.numeric(expected))
  })

  test_that("pack_comp_osa (LN, ct=0) produces the ALR transform of the proportions", {
    n_obs_bins <- 4
    props <- c(0.1, 0.2, 0.3, 0.4)
    prop_arr <- make_prop_array(props, 2, 1)

    g <- pack_single_group(prop_arr, ISS = array(50, dim = c(2, 1)),
                           Wt = array(1, dim = c(2, 1)), use = c(1L, 1L),
                           n_obs_bins = n_obs_bins, n_sexes = 1,
                           comp_type = 0, like_type = 2, addtocomp = 0)

    pr_expected <- props / sum(props)
    expected <- alr(pr_expected)
    expect_equal(length(g), n_obs_bins - 1)
    expect_equal(as.numeric(g), as.numeric(expected), tolerance = 1e-8)
  })

  # eval_comp_osa: multi-fleet round trip, zero_init behavior
  test_that("eval_comp_osa round-trips pack_comp_osa output for a single discrete fleet", {
    n_regions <- 2; n_obs_bins <- 4; n_model_bins <- 4; n_sexes <- 1

    Exp_true <- normalize_rs(array(c(0.1, 0.2, 0.3, 0.4,
                                     0.1, 0.2, 0.3, 0.4),
                                   dim = c(n_regions, n_model_bins, n_sexes)))
    Obs_raw  <- array(c(0.1, 0.2, 0.3, 0.4,
                        0.1, 0.2, 0.3, 0.4), dim = c(n_regions, n_obs_bins, n_sexes))

    ObsArr <- array(Obs_raw, dim = c(n_regions, 1, 1, n_obs_bins, n_sexes, 1))
    ISSArr <- array(200, dim = c(n_regions, 1, 1, n_sexes, 1))
    WtArr  <- array(1,   dim = c(n_regions, 1, 1, n_sexes, 1))
    UseArr <- array(1L,  dim = c(n_regions, 1, 1, 1))
    TypeMat <- matrix(1, nrow = 1, ncol = 1)      # Comp_Type = 1 (split)
    LikeTypeVec <- 0                              # multinomial

    tracked <- pack_comp_osa(ObsArr, ISSArr, WtArr, UseArr, TypeMat, LikeTypeVec,
                             n_yrs = 1, n_seas = 1, n_fleets = 1, n_sexes = n_sexes,
                             addtocomp = 1e-4, family = "discrete",
                             pop = FALSE, n_pop = 1)

    ExpArrFn <- function(p, y, seas, f) Exp_true
    lnThetaArr    <- array(0, dim = c(n_regions, n_sexes, 1))
    lnThetaAggVec <- c(0)
    LNcorrArr     <- array(0, dim = c(n_regions, n_sexes, 1, 3))
    LNcorrAggVec  <- c(0)

    nLL_init <- array(0, dim = c(n_regions, 1, 1, n_sexes, 1))

    res <- eval_comp_osa(
      nLL_arr = nLL_init, tracked = tracked, ExpArrFn = ExpArrFn,
      UseArr = UseArr, TypeMat = TypeMat, LikeTypeVec = LikeTypeVec,
      ISSArr = ISSArr, lnThetaArr = lnThetaArr, lnThetaAggVec = lnThetaAggVec,
      LNcorrArr = LNcorrArr, LNcorrAggVec = LNcorrAggVec,
      n_regions = n_regions, n_yrs = 1, n_seas = 1, n_fleets = 1, n_sexes = n_sexes,
      n_model_bins = n_model_bins, n_obs_bins = n_obs_bins, age_or_len = 1,
      AgeingErrorFn = NULL, addtocomp = 1e-4,
      family = "discrete", zero_init = TRUE, pop = FALSE, n_pop = 1
    )

    # cross-check against a direct call of Get_Comp_Likelihoods_OSA on the
    # same tracked vector
    direct <- Get_Comp_Likelihoods_OSA(
      Exp = Exp_true, Obs = tracked,
      ISS = array(200, dim = c(n_regions, n_sexes)),
      ln_theta = array(0, dim = c(n_regions, n_sexes)), ln_theta_agg = 0,
      LN_corr_pars = array(0, dim = c(n_regions, n_sexes, 3)), LN_corr_pars_agg = 0,
      Comp_Type = 1, Likelihood_Type = 0,
      n_regions = n_regions, n_model_bins = n_model_bins, n_obs_bins = n_obs_bins,
      n_sexes = n_sexes, age_or_len = 1, AgeingError = identity_ae(n_obs_bins),
      use = rep(1L, n_regions), addtocomp = 1e-4
    )

    expect_equal(as.numeric(res[, 1, 1, , 1]), as.numeric(direct[, 1]), tolerance = 1e-8)
    expect_true(all(is.finite(res)))
  })

  test_that("eval_comp_osa preserves other fleets' slots when zero_init = FALSE, and wipes them when TRUE", {
    n_regions <- 2; n_obs_bins <- 4; n_model_bins <- 4; n_sexes <- 1; n_fleets <- 2

    # Fleet 1: multinomial (discrete), Comp_Type = 1
    # Fleet 2: logistic-normal iid (continuous), Comp_Type = 1
    LikeTypeVec <- c(0, 2)
    TypeMat <- matrix(c(1, 1), nrow = 1, ncol = 2)

    Obs_raw <- array(0.25, dim = c(n_regions, n_obs_bins, n_sexes))
    ObsArr <- array(rep(as.numeric(Obs_raw), n_fleets),
                    dim = c(n_regions, 1, 1, n_obs_bins, n_sexes, n_fleets))
    ISSArr <- array(100, dim = c(n_regions, 1, 1, n_sexes, n_fleets))
    WtArr  <- array(1,   dim = c(n_regions, 1, 1, n_sexes, n_fleets))
    UseArr <- array(1L,  dim = c(n_regions, 1, 1, n_fleets))

    tracked_discrete   <- pack_comp_osa(ObsArr, ISSArr, WtArr, UseArr, TypeMat, LikeTypeVec,
                                        n_yrs = 1, n_seas = 1, n_fleets = n_fleets, n_sexes = n_sexes,
                                        addtocomp = 1e-4, family = "discrete", pop = FALSE, n_pop = 1)
    tracked_continuous <- pack_comp_osa(ObsArr, ISSArr, WtArr, UseArr, TypeMat, LikeTypeVec,
                                        n_yrs = 1, n_seas = 1, n_fleets = n_fleets, n_sexes = n_sexes,
                                        addtocomp = 1e-4, family = "continuous", pop = FALSE, n_pop = 1)

    Exp_true <- normalize_rs(array(0.25, dim = c(n_regions, n_model_bins, n_sexes)))
    ExpArrFn <- function(p, y, seas, f) Exp_true

    lnThetaArr    <- array(0, dim = c(n_regions, n_sexes, n_fleets))
    lnThetaAggVec <- rep(0, n_fleets)
    LNcorrArr     <- array(0, dim = c(n_regions, n_sexes, n_fleets, 3))
    LNcorrAggVec  <- rep(0, n_fleets)

    nLL_init <- array(0, dim = c(n_regions, 1, 1, n_sexes, n_fleets))

    res1 <- eval_comp_osa(nLL_init, tracked_discrete, ExpArrFn, UseArr, TypeMat, LikeTypeVec,
                          ISSArr, lnThetaArr, lnThetaAggVec, LNcorrArr, LNcorrAggVec,
                          n_regions = n_regions, n_yrs = 1, n_seas = 1, n_fleets = n_fleets,
                          n_sexes = n_sexes, n_model_bins = n_model_bins, n_obs_bins = n_obs_bins,
                          age_or_len = 1, AgeingErrorFn = NULL, addtocomp = 1e-4,
                          family = "discrete", zero_init = TRUE, pop = FALSE, n_pop = 1)

    # fleet 1 (discrete) populated, fleet 2 (continuous) still zero
    expect_true(all(is.finite(res1[, 1, 1, , 1])))
    expect_true(all(res1[, 1, 1, , 1] != 0))
    expect_true(all(res1[, 1, 1, , 2] == 0))

    res2 <- eval_comp_osa(res1, tracked_continuous, ExpArrFn, UseArr, TypeMat, LikeTypeVec,
                          ISSArr, lnThetaArr, lnThetaAggVec, LNcorrArr, LNcorrAggVec,
                          n_regions = n_regions, n_yrs = 1, n_seas = 1, n_fleets = n_fleets,
                          n_sexes = n_sexes, n_model_bins = n_model_bins, n_obs_bins = n_obs_bins,
                          age_or_len = 1, AgeingErrorFn = NULL, addtocomp = 1e-4,
                          family = "continuous", zero_init = FALSE, pop = FALSE, n_pop = 1)

    # fleet 1's values from the discrete pass are preserved...
    expect_equal(as.numeric(res2[, 1, 1, , 1]), as.numeric(res1[, 1, 1, , 1]))
    # ...and fleet 2 is now populated
    expect_true(all(is.finite(res2[, 1, 1, , 2])))
    expect_true(all(res2[, 1, 1, , 2] != 0))

    # zero_init = TRUE wipes the previously-filled discrete slot
    res3 <- eval_comp_osa(res2, tracked_continuous, ExpArrFn, UseArr, TypeMat, LikeTypeVec,
                          ISSArr, lnThetaArr, lnThetaAggVec, LNcorrArr, LNcorrAggVec,
                          n_regions = n_regions, n_yrs = 1, n_seas = 1, n_fleets = n_fleets,
                          n_sexes = n_sexes, n_model_bins = n_model_bins, n_obs_bins = n_obs_bins,
                          age_or_len = 1, AgeingErrorFn = NULL, addtocomp = 1e-4,
                          family = "continuous", zero_init = TRUE, pop = FALSE, n_pop = 1)

    expect_true(all(res3[, 1, 1, , 1] == 0))       # fleet 1 wiped
    expect_true(all(res3[, 1, 1, , 2] != 0))       # fleet 2 still populated
  })

  test_that("eval_comp_osa returns nLL_arr unchanged when tracked is NULL", {
    n_regions <- 2; n_sexes <- 1
    nLL_init <- array(runif(n_regions * n_sexes), dim = c(n_regions, 1, 1, n_sexes, 1))
    UseArr <- array(1L, dim = c(n_regions, 1, 1, 1))
    TypeMat <- matrix(1, nrow = 1, ncol = 1)
    LikeTypeVec <- 0
    ISSArr <- array(100, dim = c(n_regions, 1, 1, n_sexes, 1))
    lnThetaArr <- array(0, dim = c(n_regions, n_sexes, 1))
    lnThetaAggVec <- c(0)
    LNcorrArr <- array(0, dim = c(n_regions, n_sexes, 1, 3))
    LNcorrAggVec <- c(0)

    res <- eval_comp_osa(nLL_init, tracked = NULL,
                         ExpArrFn = function(p, y, seas, f) stop("should not be called"),
                         UseArr = UseArr, TypeMat = TypeMat, LikeTypeVec = LikeTypeVec,
                         ISSArr = ISSArr, lnThetaArr = lnThetaArr, lnThetaAggVec = lnThetaAggVec,
                         LNcorrArr = LNcorrArr, LNcorrAggVec = LNcorrAggVec,
                         n_regions = n_regions, n_yrs = 1, n_seas = 1, n_fleets = 1,
                         n_sexes = n_sexes, n_model_bins = 4, n_obs_bins = 4, age_or_len = 1,
                         AgeingErrorFn = NULL, addtocomp = 1e-4,
                         family = "discrete", zero_init = TRUE, pop = FALSE, n_pop = 1)

    expect_identical(res, nLL_init)
  })
})
