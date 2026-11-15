# get_comp_prop's aggregated compositions (Comp_Type 0) are the total over regions and sexes in each bin, as the
# likelihood fits them. Restrc_Comps flattened a region-first array by bin, which mixed regions and bins.

library(SPoRC)
library(testthat)

test_that("aggregated comps sum over regions and sexes bin by bin", {

  n_regions <- 3
  n_bins <- 4
  n_sexes <- 2

  # expected values that differ by region, bin and sex
  Exp <- array(0, dim = c(n_regions, 1, 1, n_bins, n_sexes, 1))
  for(r in 1:n_regions) for(b in 1:n_bins) for(s in 1:n_sexes) Exp[r,1,1,b,s,1] <- 10 * r + b^2 + s / 10
  Obs <- array(1, dim = dim(Exp))

  # the total over regions and sexes in each bin, as proportions
  agg_total <- apply(Exp, 4, sum)
  agg_prop <- agg_total / sum(agg_total)

  # lengths
  len_comps <- SPoRC:::Restrc_Comps(Exp = Exp, Obs = Obs, Comp_Type = 0, age_or_len = 1, AgeingError = NA)
  expect_equal(as.numeric(len_comps$Exp[1,,1]), agg_prop, tolerance = 1e-14)

  # ages, through ageing error
  ageing_error <- matrix(c(0.8, 0.2, 0.0, 0.0,
                           0.1, 0.8, 0.1, 0.0,
                           0.0, 0.1, 0.8, 0.1,
                           0.0, 0.0, 0.2, 0.8), nrow = n_bins, byrow = TRUE)
  age_comps <- SPoRC:::Restrc_Comps(Exp = Exp, Obs = Obs, Comp_Type = 0, age_or_len = 0, AgeingError = ageing_error)
  age_prop <- as.numeric(agg_prop %*% ageing_error)
  expect_equal(as.numeric(age_comps$Exp[1,,1]), age_prop / sum(age_prop), tolerance = 1e-14)

  # one region is unchanged
  one_region <- SPoRC:::Restrc_Comps(Exp = Exp[1,,,,,,drop = FALSE], Obs = Obs[1,,,,,,drop = FALSE], Comp_Type = 0, age_or_len = 1, AgeingError = NA)
  one_total <- apply(Exp[1,,,,,,drop = FALSE], 4, sum)
  expect_equal(as.numeric(one_region$Exp[1,,1]), one_total / sum(one_total), tolerance = 1e-14)
})
