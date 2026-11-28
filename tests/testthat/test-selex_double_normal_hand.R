# Double normal selectivity against the documented equations written out by hand (1e-10): both limbs anchored at
# the end bins, either limb left as its own curve, a plateau height below one, and a start bin above the first.

library(SPoRC)
library(testthat)

# Setup ----------------------------------------------------------------------

dbnrml_hand <- list(
  bins = 1:12,
  pars = c(4, 0, log(4), log(9), qlogis(0.1), qlogis(0.3)) # p1 to p6 before their transforms
)

dbnrml_model <- function(apical = 1, dbnrml_raw = c(0, 0), dbnrml_startbin = 1) {
  as.numeric(Get_Selex(
    Selex_Model = 4,
    TimeVary_Model = 0,
    pars = dbnrml_hand$pars,
    ln_seldevs = array(0, dim = c(1, 1, 6, 1, 1)),
    Region = 1,
    Year = 1,
    Bin = dbnrml_hand$bins,
    Sex = 1,
    apical = apical,
    dbnrml_raw = dbnrml_raw,
    dbnrml_startbin = dbnrml_startbin
  ))
}

# the curve one bin at a time, from the equations in the model equations vignette
dbnrml_calc <- function(apical = 1, dbnrml_raw = c(0, 0), dbnrml_startbin = 1) {

  bins <- dbnrml_hand$bins
  binwidth <- bins[2] - bins[1]
  p1 <- dbnrml_hand$pars[1] # bin where the plateau begins
  p2 <- p1 + binwidth + (0.99 * max(bins) - p1 - binwidth) / (1 + exp(-dbnrml_hand$pars[2])) # bin where it ends
  p3 <- exp(dbnrml_hand$pars[3]) # ascending width
  p4 <- exp(dbnrml_hand$pars[4]) # descending width
  p5 <- 1 / (1 + exp(-dbnrml_hand$pars[5])) # selectivity at the first bin
  p6 <- 1 / (1 + exp(-dbnrml_hand$pars[6])) # selectivity at the last bin

  asc_min <- exp(-(bins[dbnrml_startbin] - p1)^2 / p3) # ascending limb at the start bin
  desc_min <- exp(-(max(bins) - p2)^2 / p4) # descending limb at the last bin
  selex <- rep(0, length(bins))

  for(b in seq_along(bins)) {

    asc <- exp(-(bins[b] - p1)^2 / p3)
    desc <- exp(-(bins[b] - p2)^2 / p4)

    # a limb is either anchored at its end bin or left as its own curve built up to the plateau height
    asc_limb <- p5 + (apical - p5) * (asc - asc_min) / (1 - asc_min)
    if(dbnrml_raw[1] == 1) asc_limb <- apical * asc
    desc_limb <- apical + (p6 - apical) * (desc - 1) / (desc_min - 1)
    if(dbnrml_raw[2] == 1) desc_limb <- apical * desc

    join1 <- 1 / (1 + exp(-20 * (bins[b] - p1) / (1 + abs(bins[b] - p1))))
    join2 <- 1 / (1 + exp(-20 * (bins[b] - p2) / (1 + abs(bins[b] - p2))))
    selex[b] <- asc_limb * (1 - join1) + join1 * (apical * (1 - join2) + desc_limb * join2)

  } # end b loop

  # bins below the start bin fall off as the square of their bin relative to it
  if(dbnrml_startbin > 1) {
    for(b in 1:(dbnrml_startbin - 1)) selex[b] <- (bins[b] / bins[dbnrml_startbin])^2 * selex[dbnrml_startbin]
  }

  selex
}

# Tests ----------------------------------------------------------------------

test_that("the double normal matches the documented equations under every limb option", {

  for(dbnrml_raw in list(c(0, 0), c(1, 0), c(0, 1), c(1, 1))) {
    for(apical in c(1, 0.8)) {
      expect_equal(dbnrml_model(apical = apical, dbnrml_raw = dbnrml_raw), dbnrml_calc(apical = apical, dbnrml_raw = dbnrml_raw), tolerance = 1e-10,
                   label = sprintf("limbs left as their own curve (asc, desc) = (%d, %d), plateau height %g", dbnrml_raw[1], dbnrml_raw[2], apical))
    } # end apical loop
  } # end dbnrml_raw loop

  # a start bin above the first
  expect_equal(dbnrml_model(dbnrml_startbin = 3), dbnrml_calc(dbnrml_startbin = 3), tolerance = 1e-10)

})

test_that("anchored limbs end on the first and last bin selectivities, and a free descending limb does not", {

  anchored <- dbnrml_model()
  expect_equal(anchored[1], 0.1, tolerance = 1e-4) # p5, within what the joiner leaves at the first bin
  expect_equal(anchored[12], 0.3, tolerance = 1e-4) # p6
  expect_equal(max(anchored), 1, tolerance = 1e-4) # the plateau

  # left as its own curve the descending limb ends where the curve takes it, below the anchored value here
  free_desc <- dbnrml_model(dbnrml_raw = c(0, 1))
  expect_lt(free_desc[12], 0.25)
  expect_equal(free_desc[1:5], anchored[1:5], tolerance = 1e-4) # the ascending side is unchanged

  # the plateau height scales the top of the curve and leaves the anchors
  lower <- dbnrml_model(apical = 0.8)
  expect_equal(max(lower), 0.8, tolerance = 1e-4)
  expect_equal(lower[1], 0.1, tolerance = 1e-4)

})
