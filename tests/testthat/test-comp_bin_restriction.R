library(SPoRC)
library(testthat)

# Composition bin restriction (the *_bins arguments). The invariant the whole
# feature rests on: fitting a stream over a subset of its observed bins must give
# exactly the likelihood you would get by handing in only those bins in the first
# place, renormalized within them. That has to hold for every composition type,
# not just the aggregated one, and for lengths as well as ages.

# Expected composition, mildly peaked so the bins are not interchangeable
peaked_exp <- function(n_regions, n_bins, n_sexes) {
  arr <- array(1, dim = c(n_regions, n_bins, n_sexes))
  for(r in seq_len(n_regions)) for(s in seq_len(n_sexes)) {
    arr[r,,s] <- seq(1, 3, length.out = n_bins) * (1 + 0.1 * r + 0.2 * s)
    arr[r,,s] <- arr[r,,s] / sum(arr[r,,s])
  } # end r and s loops
  arr
}

# Observed counts, deliberately uneven across bins and sexes
lumpy_obs <- function(n_regions, n_bins, n_sexes, total = 200) {
  arr <- array(0, dim = c(n_regions, n_bins, n_sexes))
  for(r in seq_len(n_regions)) for(s in seq_len(n_sexes)) {
    v <- (seq_len(n_bins) %% 4) + 1 + r + s
    arr[r,,s] <- total * v / sum(v)
  } # end r and s loops
  arr
}

comp_nll <- function(Exp, Obs, n_model_bins, n_obs_bins, Comp_Type, Likelihood_Type,
                     n_regions, n_sexes, age_or_len = 0, AgeingError = NULL,
                     comp_bins = NULL, use = NULL) {
  if(is.null(use)) use <- rep(1L, n_regions)
  if(is.null(AgeingError) && age_or_len == 0) AgeingError <- diag(1, n_model_bins)
  if(is.null(AgeingError)) AgeingError <- NA
  Get_Comp_Likelihoods(
    Exp = Exp, Obs = Obs,
    ISS = array(50, dim = c(n_regions, n_sexes)),
    Wt_Mltnml = array(1, dim = c(n_regions, n_sexes)),
    ln_theta = array(0, dim = c(n_regions, n_sexes)), ln_theta_agg = 0,
    LN_corr_pars = array(0.3, dim = c(n_regions, n_sexes, 3)), LN_corr_pars_agg = 0.3,
    Comp_Type = Comp_Type, Likelihood_Type = Likelihood_Type,
    n_regions = n_regions, n_model_bins = n_model_bins, n_obs_bins = n_obs_bins,
    n_sexes = n_sexes, age_or_len = age_or_len, AgeingError = AgeingError,
    use = use, addtocomp = 1e-4, comp_bins = comp_bins
  )
}

test_that("restricting to a contiguous bin range equals fitting only those bins", {
  n_bins <- 8
  keep <- 3:7

  for(ct in 0:2) {
    # The aggregated comps read Obs[1,,1] and reshape Exp as though bins ran
    # fastest, so they are a single-region construct and are exercised as one
    n_regions <- if(ct == 0) 1 else 2
    n_sexes <- 2
    Exp <- peaked_exp(n_regions, n_bins, n_sexes)
    Obs <- lumpy_obs(n_regions, n_bins, n_sexes)
    # LikeType 4 is a sex-joint family, so it only applies to Comp_Type 2.
    # Comp_Type 2 with LikeType 3 runs one AR1 across the whole [bin x sex]
    # stack, so its lags do not survive being re-read as a shorter stack; it
    # gets its own test below.
    lts <- if(ct == 2) c(0, 1, 2, 4) else 0:3
    for(lt in lts) {
      restricted <- comp_nll(Exp, Obs, n_bins, n_bins, ct, lt, n_regions, n_sexes,
                             comp_bins = keep)
      direct <- comp_nll(Exp[, keep, , drop = FALSE], Obs[, keep, , drop = FALSE],
                         length(keep), length(keep), ct, lt, n_regions, n_sexes)
      expect_equal(restricted, direct, tolerance = 1e-10,
                   info = sprintf("Comp_Type=%d Likelihood_Type=%d", ct, lt))
    } # end lt loop
  } # end ct loop
})

test_that("the sex-joint AR1 keeps true bin lags inside each sex block", {
  # Comp_Type 2 with LikeType 3 lays one AR1 over the [bin x sex] stack. Cutting
  # that covariance down to the fitted cells keeps the lag between two fitted
  # bins of the same sex equal to their true separation, which is what the
  # correlation is meant to describe. The lag across the sex boundary stretches
  # as a side effect, and is arbitrary either way.
  n_regions <- 1; n_sexes <- 2; n_bins <- 8
  keep <- 3:7
  Exp <- peaked_exp(n_regions, n_bins, n_sexes)
  Obs <- lumpy_obs(n_regions, n_bins, n_sexes)
  restricted <- comp_nll(Exp, Obs, n_bins, n_bins, 2, 3, n_regions, n_sexes, comp_bins = keep)
  expect_true(all(is.finite(restricted)))

  rho <- 2 / (1 + exp(-2 * 0.3)) - 1
  full <- get_AR1_CorrMat(n_bins * n_sexes, rho)
  joint <- as.vector(outer(keep, (seq_len(n_sexes) - 1) * n_bins, "+"))
  cut <- as.matrix(full)[joint, joint]
  # within sex 1, bins 3 and 4 stay one bin apart
  expect_equal(cut[1, 2], as.matrix(full)[3, 4], tolerance = 1e-12)
  # across the sex boundary the lag is longer than a naive short stack would give
  expect_lt(cut[length(keep), length(keep) + 1],
            as.matrix(get_AR1_CorrMat(length(keep) * n_sexes, rho))[length(keep), length(keep) + 1])
})

test_that("naming every bin is a no-op", {
  n_regions <- 2; n_sexes <- 2; n_bins <- 6
  Exp <- peaked_exp(n_regions, n_bins, n_sexes)
  Obs <- lumpy_obs(n_regions, n_bins, n_sexes)
  for(ct in 0:2) {
    lts <- if(ct == 2) 0:4 else 0:3
    for(lt in lts) {
      expect_equal(comp_nll(Exp, Obs, n_bins, n_bins, ct, lt, n_regions, n_sexes,
                            comp_bins = seq_len(n_bins)),
                   comp_nll(Exp, Obs, n_bins, n_bins, ct, lt, n_regions, n_sexes),
                   tolerance = 1e-12,
                   info = sprintf("Comp_Type=%d Likelihood_Type=%d", ct, lt))
    } # end lt loop
  } # end ct loop
})

test_that("length compositions restrict the same way ages do", {
  n_sexes <- 1; n_bins <- 10
  keep <- 4:9
  for(ct in 0:2) {
    n_regions <- if(ct == 0) 1 else 2   # aggregated comps are single-region, as above
    Exp <- peaked_exp(n_regions, n_bins, n_sexes)
    Obs <- lumpy_obs(n_regions, n_bins, n_sexes)
    for(lt in 0:2) {
      restricted <- comp_nll(Exp, Obs, n_bins, n_bins, ct, lt, n_regions, n_sexes,
                             age_or_len = 1, AgeingError = NA, comp_bins = keep)
      direct <- comp_nll(Exp[, keep, , drop = FALSE], Obs[, keep, , drop = FALSE],
                         length(keep), length(keep), ct, lt, n_regions, n_sexes,
                         age_or_len = 1, AgeingError = NA)
      expect_equal(restricted, direct, tolerance = 1e-10,
                   info = sprintf("lengths Comp_Type=%d Likelihood_Type=%d", ct, lt))
    } # end lt loop
  } # end ct loop
})

test_that("a length bin map and a bin restriction compose", {
  # 8 model bins collapsed onto 4 observed bins, then fitted over 2 of those 4
  n_regions <- 1; n_sexes <- 1; n_model <- 8; n_obs <- 4
  map <- matrix(0, nrow = n_model, ncol = n_obs)
  for(i in seq_len(n_model)) map[i, ceiling(i / 2)] <- 1
  Exp <- peaked_exp(n_regions, n_model, n_sexes)
  Obs <- lumpy_obs(n_regions, n_obs, n_sexes)
  keep <- 2:3

  restricted <- comp_nll(Exp, Obs, n_model, n_obs, 1, 0, n_regions, n_sexes,
                         age_or_len = 1, AgeingError = map, comp_bins = keep)
  # mapping by hand, then fitting the kept observed bins directly
  mapped <- array((Exp[1,,1] / sum(Exp[1,,1])) %*% map, dim = c(1, n_obs, 1))
  direct <- comp_nll(mapped[, keep, , drop = FALSE], Obs[, keep, , drop = FALSE],
                     length(keep), length(keep), 1, 0, n_regions, n_sexes,
                     age_or_len = 1, AgeingError = NA)
  expect_equal(restricted, direct, tolerance = 1e-10)
})

test_that("the sex-joint stack is restricted within each sex, not across the stack", {
  # Comp_Type 2 evaluates a [bin x sex] stack. Restricting must drop the named
  # bins from every sex's block, leaving n_fit_bins * n_sexes cells, rather than
  # slicing the flattened stack.
  n_regions <- 1; n_sexes <- 2; n_bins <- 6
  keep <- c(2, 5)
  Exp <- peaked_exp(n_regions, n_bins, n_sexes)
  Obs <- lumpy_obs(n_regions, n_bins, n_sexes)
  restricted <- comp_nll(Exp, Obs, n_bins, n_bins, 2, 0, n_regions, n_sexes,
                         comp_bins = keep)
  direct <- comp_nll(Exp[, keep, , drop = FALSE], Obs[, keep, , drop = FALSE],
                     length(keep), length(keep), 2, 0, n_regions, n_sexes)
  expect_equal(restricted, direct, tolerance = 1e-10)
  # and the sex ratio the joint comps carry is now the ratio within the kept bins
  expect_false(isTRUE(all.equal(
    as.numeric(restricted),
    as.numeric(comp_nll(Exp, Obs, n_bins, n_bins, 2, 0, n_regions, n_sexes))
  )))
})

test_that("AR1 lags are measured over the observed range, not the fitted one", {
  # With a gap in the kept bins, the covariance between the bins either side of
  # the gap must reflect the true lag, so a gapped restriction differs from
  # simply handing in the kept bins as if they were adjacent.
  n_regions <- 1; n_sexes <- 1; n_bins <- 8
  keep <- c(1, 2, 7, 8)
  Exp <- peaked_exp(n_regions, n_bins, n_sexes)
  Obs <- lumpy_obs(n_regions, n_bins, n_sexes)
  gapped <- comp_nll(Exp, Obs, n_bins, n_bins, 1, 3, n_regions, n_sexes, comp_bins = keep)
  as_adjacent <- comp_nll(Exp[, keep, , drop = FALSE], Obs[, keep, , drop = FALSE],
                          length(keep), length(keep), 1, 3, n_regions, n_sexes)
  expect_true(all(is.finite(gapped)))
  expect_false(isTRUE(all.equal(as.numeric(gapped), as.numeric(as_adjacent))))
})

test_that("parse_bin_subset accepts the list and array forms alike", {
  n_bins <- 10; n_fleets <- 3
  from_list <- parse_bin_subset(list(NULL, 2:6, 10), n_bins, n_fleets, "x")
  arr <- array(1, dim = c(n_bins, n_fleets))
  arr[, 2] <- 0; arr[2:6, 2] <- 1
  arr[, 3] <- 0; arr[10, 3] <- 1
  expect_equal(from_list, arr)
  expect_equal(parse_bin_subset(arr, n_bins, n_fleets, "x"), arr)
  expect_equal(parse_bin_subset(NULL, n_bins, n_fleets, "x"),
               array(1, dim = c(n_bins, n_fleets)))

  expect_error(parse_bin_subset(list(1:3, 2:6), n_bins, n_fleets, "x"), "length 2")
  expect_error(parse_bin_subset(list(NULL, NULL, 99), n_bins, n_fleets, "x"), "outside 1:10")
  expect_error(parse_bin_subset(list(NULL, NULL, integer(0)), n_bins, n_fleets, "x"), "no bins")
})

test_that("check_bin_map holds AgeingError and LenBinMap to the same rules", {
  ok <- diag(1, 6)
  expect_silent(check_bin_map(ok, 6, "AgeingError"))
  # a row of zeros drops a model bin from the observations, which is allowed
  shifted <- diag(1, 6)[, 2:6]
  expect_silent(check_bin_map(shifted, 6, "AgeingError"))
  # a many-to-one collapse, the length bin map's whole purpose
  collapse <- matrix(0, 6, 3)
  for(i in 1:6) collapse[i, ceiling(i / 2)] <- 1
  expect_silent(check_bin_map(collapse, 6, "LenBinMap"))

  # published ageing error matrices are rounded at source, so ordinary rounding
  # must not trip the row-sum rule for them. A length bin map is written by hand
  # rather than read off a rounded table, so it keeps the tight default it has
  # always been held to.
  rounded <- ok; rounded[1, 1] <- 0.997; rounded[2, 2] <- 1.002
  expect_silent(check_bin_map(rounded, 6, "AgeingError", strict = FALSE, tol = 0.05))
  expect_error(check_bin_map(rounded, 6, "LenBinMap"), "sum to neither")

  expect_error(check_bin_map(ok, 5, "AgeingError"), "one row per model bin")
  half <- ok; half[1, 1] <- 0.5
  # strict, which is what LenBinMap has always been
  expect_error(check_bin_map(half, 6, "LenBinMap"), "sum to neither")
  # non-strict, which is what AgeingError gets so an existing model still runs
  messages_list <<- character(0)
  expect_silent(check_bin_map(half, 6, "AgeingError", strict = FALSE))
  expect_true(any(grepl("sum to neither", messages_list)))

  empty_col <- cbind(diag(1, 6), 0)
  expect_error(check_bin_map(empty_col, 6, "LenBinMap"), "nothing mapped into them")
  messages_list <<- character(0)
  expect_silent(check_bin_map(empty_col, 6, "AgeingError", strict = FALSE))
  expect_true(any(grepl("nothing mapped into them", messages_list)))
  neg <- ok; neg[2, 3] <- -1; neg[2, 2] <- 2
  expect_error(check_bin_map(neg, 6, "AgeingError"), "negative")
})

# One-step-ahead residuals have to be packed and walked over the same restricted
# bins the likelihood fits, or the diagnostic silently disagrees with the fit.
# The packer and the evaluator keep a shared pointer, so a stride computed on the
# wrong bin count desynchronizes every group after the first.

test_that("the OSA packer and evaluator agree on the restricted bin count", {
  n_regions <- 2; n_sexes <- 2; n_obs_bins <- 8; n_fleets <- 2
  keep <- 3:7
  BinsArr <- array(1, dim = c(n_obs_bins, n_fleets))
  BinsArr[, 1] <- 0; BinsArr[keep, 1] <- 1   # fleet 1 restricted, fleet 2 not

  ObsArr <- array(0, dim = c(n_regions, 2, 1, n_obs_bins, n_sexes, n_fleets))
  for(r in 1:n_regions) for(y in 1:2) for(s in 1:n_sexes) for(f in 1:n_fleets) {
    v <- (seq_len(n_obs_bins) %% 3) + 1 + r + s + f
    ObsArr[r, y, 1, , s, f] <- 100 * v / sum(v)
  } # end loops
  ISSArr <- array(60, dim = c(n_regions, 2, 1, n_sexes, n_fleets))
  WtArr <- array(1, dim = c(n_regions, 2, 1, n_sexes, n_fleets))
  UseArr <- array(1, dim = c(n_regions, 2, 1, n_fleets))

  for(ct in 0:2) {
    for(lt in 0:3) {
      TypeMat <- matrix(ct, nrow = 2, ncol = n_fleets)
      fam <- if(lt %in% c(0, 1)) "discrete" else "continuous"
      vec <- pack_comp_osa(
        ObsArr = ObsArr, ISSArr = ISSArr, WtArr = WtArr, UseArr = UseArr,
        TypeMat = TypeMat, LikeTypeVec = rep(lt, n_fleets),
        n_yrs = 2, n_seas = 1, n_fleets = n_fleets, n_sexes = n_sexes,
        addtocomp = 1e-4, family = fam, BinsArr = BinsArr
      )

      # length the evaluator will walk, worked out independently of the packer
      block <- function(nb, n_ru) {
        if(lt %in% c(0, 1)) return(if(ct == 0) nb else n_ru * nb * n_sexes)
        if(ct == 0) return(nb - 1)
        if(ct == 1) return(n_ru * (nb - 1) * n_sexes)
        n_ru * (nb * n_sexes - 1)
      }
      expected_len <- 2 * (block(length(keep), n_regions) + block(n_obs_bins, n_regions))
      expect_equal(length(vec), expected_len,
                   info = sprintf("packed length Comp_Type=%d Likelihood_Type=%d", ct, lt))

      # the evaluator must consume exactly that vector, leaving no group adrift
      nLL <- array(0, dim = c(n_regions, 2, 1, n_sexes, n_fleets))
      out <- eval_comp_osa(
        nLL_arr = nLL, tracked = vec,
        ExpArrFn = function(p, y, seas, f) {
          e <- array(0, dim = c(n_regions, n_obs_bins, n_sexes))
          for(r in 1:n_regions) for(s in 1:n_sexes) {
            e[r,,s] <- seq(1, 2, length.out = n_obs_bins)
            e[r,,s] <- e[r,,s] / sum(e[r,,s])
          } # end loops
          e
        },
        UseArr = UseArr, TypeMat = TypeMat, LikeTypeVec = rep(lt, n_fleets),
        ISSArr = ISSArr, lnThetaArr = array(0, dim = c(n_regions, n_sexes, n_fleets)),
        lnThetaAggVec = rep(0, n_fleets),
        LNcorrArr = array(0.3, dim = c(n_regions, n_sexes, n_fleets, 3)),
        LNcorrAggVec = rep(0.3, n_fleets),
        n_regions = n_regions, n_yrs = 2, n_seas = 1, n_fleets = n_fleets,
        n_sexes = n_sexes, n_model_bins = n_obs_bins, n_obs_bins = n_obs_bins,
        age_or_len = 1, AgeingErrorFn = NULL, addtocomp = 1e-4,
        family = fam, BinsArr = BinsArr
      )
      expect_true(all(is.finite(out)),
                  info = sprintf("finite nLL Comp_Type=%d Likelihood_Type=%d", ct, lt))
    } # end lt loop
  } # end ct loop
})

test_that("OSA labels report true observed bin numbers under a restriction", {
  n_regions <- 1; n_sexes <- 1; n_obs_bins <- 6
  keep <- 2:5
  BinsArr <- array(0, dim = c(n_obs_bins, 1)); BinsArr[keep, 1] <- 1

  ObsArr <- array(0, dim = c(n_regions, 1, 1, n_obs_bins, n_sexes, 1))
  ObsArr[1, 1, 1, , 1, 1] <- 100 * (1:n_obs_bins) / sum(1:n_obs_bins)

  res <- pack_comp_osa(
    ObsArr = ObsArr, ISSArr = array(50, dim = c(n_regions, 1, 1, n_sexes, 1)),
    WtArr = array(1, dim = c(n_regions, 1, 1, n_sexes, 1)),
    UseArr = array(1, dim = c(n_regions, 1, 1, 1)),
    TypeMat = matrix(0, 1, 1), LikeTypeVec = 0,
    n_yrs = 1, n_seas = 1, n_fleets = 1, n_sexes = n_sexes,
    addtocomp = 1e-4, family = "discrete", return_labels = TRUE, BinsArr = BinsArr
  )
  expect_equal(nrow(res$labels), length(keep))
  # bins are named by their observed index, not renumbered 1..n_fit_bins, so a
  # restricted fleet's residuals stay comparable with an unrestricted one's
  expect_equal(res$labels$bin, keep)
  expect_true(res$labels$last_in_group[length(keep)])
  expect_false(any(res$labels$last_in_group[-length(keep)]))
})

test_that("an empty fitted block contributes nothing rather than NaN", {
  # restricting onto bins a region never sampled leaves that region's block empty.
  # Comp_Type 1 has always skipped such a block; Comp_Type 2 must too, since a
  # restriction can empty a block that the full composition filled.
  n_regions <- 2; n_sexes <- 2; n_bins <- 6
  keep <- 5:6
  Exp <- peaked_exp(n_regions, n_bins, n_sexes)
  Obs <- lumpy_obs(n_regions, n_bins, n_sexes)
  Obs[2, keep, ] <- 0            # region 2 has nothing in the bins being fit

  for(ct in 1:2) {
    for(lt in 0:2) {
      res <- comp_nll(Exp, Obs, n_bins, n_bins, ct, lt, n_regions, n_sexes, comp_bins = keep)
      expect_true(all(is.finite(res)),
                  info = sprintf("finite Comp_Type=%d Likelihood_Type=%d", ct, lt))
      expect_equal(unname(res[2, ]), rep(0, n_sexes), tolerance = 1e-12,
                   info = sprintf("empty block zeroed Comp_Type=%d Likelihood_Type=%d", ct, lt))
    } # end lt loop
  } # end ct loop
})

# Degenerate restrictions. A single fitted bin breaks the OSA machinery in three
# separate places, so it is refused at setup where the message can name the
# argument. And a bin array indexed on the wrong number of bins would let the
# packer and the evaluator disagree silently, which is the one failure this
# machinery cannot survive, so it is refused at the packer.

test_that("a restriction leaving fewer than two bins is refused at setup", {
  arr <- array(0, dim = c(10, 2)); arr[5, 1] <- 1; arr[, 2] <- 1
  expect_error(check_comp_bins_min(arr, c(0, 0), "FishAgeComps_bins"),
               "leaves fleet 1 with 1 fitted bin")
  # a fleet that is not fitted is skipped, since its bins are never read
  expect_silent(check_comp_bins_min(arr, c(999, 0), "FishAgeComps_bins"))
  ok <- array(0, dim = c(10, 1)); ok[4:5, 1] <- 1
  expect_silent(check_comp_bins_min(ok, 0, "FishAgeComps_bins"))
  # every family degenerates on one bin, not only the logistic-normals
  for(lt in c(0, 1, 2, 3, 4)) {
    expect_error(check_comp_bins_min(arr, c(lt, lt), "SrvAgeComps_bins"), "at least two bins")
  } # end lt loop
})

test_that("a zero-length slice request reads nothing rather than counting backwards", {
  # eval_comp_osa asked for tracked[k:(k + slice_length - 1)], which for
  # slice_length 0 counts DOWN and grabs the previous group's elements
  tracked <- 1:10
  k <- 5
  expect_equal(length(tracked[seq.int(from = k, length.out = 0)]), 0)
  expect_equal(tracked[seq.int(from = k, length.out = 3)], c(5L, 6L, 7L))
})

test_that("a bin array indexed on the wrong bin count is refused, not silently mismatched", {
  n_regions <- 1; n_sexes <- 1; n_obs_bins <- 6
  ObsArr <- array(0, dim = c(n_regions, 1, 1, n_obs_bins, n_sexes, 1))
  ObsArr[1, 1, 1, , 1, 1] <- 100 * (1:n_obs_bins) / sum(1:n_obs_bins)
  wrong <- array(1, dim = c(n_obs_bins + 3, 1))   # indexed on a different stream's bins
  wrong[1, 1] <- 0

  expect_error(
    pack_comp_osa(ObsArr = ObsArr, ISSArr = array(50, dim = c(n_regions, 1, 1, n_sexes, 1)),
                  WtArr = array(1, dim = c(n_regions, 1, 1, n_sexes, 1)),
                  UseArr = array(1, dim = c(n_regions, 1, 1, 1)),
                  TypeMat = matrix(0, 1, 1), LikeTypeVec = 0,
                  n_yrs = 1, n_seas = 1, n_fleets = 1, n_sexes = n_sexes,
                  addtocomp = 1e-4, family = "discrete", BinsArr = wrong),
    "rows but the observations carry")
})

test_that("a restriction that empties a block clears its use flag", {
  # the fitting likelihood already skips such a block; clearing the flag stops the
  # residual packer normalizing (0 + addtocomp) into a flat composition and fitting it
  n_regions <- 2; n_yrs <- 3; n_seas <- 1; n_bins <- 8; n_sexes <- 1; n_fleets <- 1
  obs <- array(0, dim = c(n_regions, n_yrs, n_seas, n_bins, n_sexes, n_fleets))
  for(r in 1:n_regions) for(y in 1:n_yrs) obs[r, y, 1, , 1, 1] <- 10
  obs[2, 2, 1, 5:8, 1, 1] <- 0            # region 2, year 2 has nothing above bin 4
  use <- array(1, dim = c(n_regions, n_yrs, n_seas, n_fleets))
  bins <- array(0, dim = c(n_bins, n_fleets)); bins[5:8, 1] <- 1

  messages_list <<- character(0)
  out <- drop_empty_fitted_blocks(obs, use, bins, 4, "FishAgeComps")
  expect_equal(out[2, 2, 1, 1], 0)                          # emptied block cleared
  expect_equal(sum(out), sum(use) - 1)                      # and only that one
  expect_true(any(grepl("use flag was cleared", messages_list)))

  # an unrestricted stream is returned untouched, with nothing reported
  messages_list <<- character(0)
  expect_equal(drop_empty_fitted_blocks(obs, use, array(1, dim = c(n_bins, n_fleets)), 4, "x"), use)
  expect_equal(length(messages_list), 0)

  # a bins array sized for a different stream is left for the packers to reject
  expect_equal(drop_empty_fitted_blocks(obs, use, array(1, dim = c(n_bins + 2, n_fleets)), 4, "x"), use)

  # a block of NAs counts as empty too, since that is what the likelihood's guard does
  obs_na <- obs
  obs_na[1, 3, 1, 5:8, 1, 1] <- NA
  out_na <- drop_empty_fitted_blocks(obs_na, use, bins, 4, "x")
  expect_equal(out_na[1, 3, 1, 1], 0)
  expect_equal(sum(out_na), sum(use) - 2)   # the NA block and the zero block, nothing else
})
