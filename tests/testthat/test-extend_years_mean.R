# extend_years(fill = "mean") builds the projection-year index standard errors and input sample sizes in
# condition_closed_loop_simulations(), where year is the second dim and sexes, fleets and sims follow it.

library(SPoRC)
library(testthat)

test_that("extend_years 'mean' keeps each sex and fleet mean in place when year is not the last dim", {
  # dims follow the closed loop ISS arrays: region, year, season, sex, fleet, sim
  n_regions <- 2
  n_hist <- 5 # conditioning years
  n_proj <- 3 # appended projection years
  n_seas <- 2
  n_sexes <- 2
  n_fleets <- 3
  n_sims <- 2
  a <- array(seq_len(n_regions * n_hist * n_seas * n_sexes * n_fleets * n_sims), # distinct value in every cell
             dim = c(n_regions, n_hist, n_seas, n_sexes, n_fleets, n_sims))
  a[1, 2, 1, 2, 3, 1] <- 0 # zeros and NA are years without data and drop out of the mean
  a[2, 4, 2, 1, 2, 2] <- NA
  a[1, , 2, 2, 1, 2] <- 0 # no valid years, so the mean is zero

  out <- SPoRC:::extend_years(a, n_years = n_proj, yr_dim = 2, fill = "mean")
  expect_equal(dim(out), c(n_regions, n_hist + n_proj, n_seas, n_sexes, n_fleets, n_sims))
  expect_equal(out[, 1:n_hist, , , , , drop = FALSE], a, ignore_attr = TRUE) # conditioning years unchanged

  # per-element mean over the conditioning years
  expected <- array(NA_real_, dim = c(n_regions, n_seas, n_sexes, n_fleets, n_sims))
  for(r in 1:n_regions) {
    for(k in 1:n_seas) {
      for(s in 1:n_sexes) {
        for(f in 1:n_fleets) {
          for(i in 1:n_sims) {
            x <- a[r, , k, s, f, i]
            x <- x[!is.na(x) & x != 0]
            expected[r, k, s, f, i] <- if(length(x) == 0) 0 else mean(x)
          } # end i loop
        } # end f loop
      } # end s loop
    } # end k loop
  } # end r loop
  expect_equal(expected[1, 2, 2, 1, 2], 0)
  expect_equal(length(unique(as.vector(expected))), length(expected)) # no two cells share a mean, so a swap cannot pass

  # every projection year holds the same per-element mean
  for(y in n_hist + 1:n_proj) expect_equal(out[, y, , , , ], expected, ignore_attr = TRUE)
})
