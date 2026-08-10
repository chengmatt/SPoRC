# Compute a model-agnostic selectivity smoothness / dome-shape penalty (Positive Scale)

Regularization penalty operating directly on a realized
selectivity-at-bin-at-year surface, rather than on any particular
parameterization's deviations. Because it only ever looks at the
resulting selectivity values, it applies uniformly to any selectivity
functional form and any fleet, called once per fleet from the
"Selectivity Smoothness Penalty" section of `SPoRC_rtmb.R`.

## Usage

``` r
Get_Selex_Smoothness_Penalty(
  sel_vals,
  wt_bin_curve = 0,
  wt_bin_diff = 0,
  wt_yr_diff = 0,
  wt_yr_curve = 0,
  wt_dome = 0,
  wt_mean_center = 0,
  normalize = TRUE,
  bin_range = NULL,
  yr_diff_ref = NULL
)
```

## Arguments

- sel_vals:

  Array of selectivity values dimensioned `[1, year, bin, sex, 1]`.
  Evaluated on the log scale internally.

- wt_bin_curve:

  Non-negative weight on the age/bin curvature penalty: the sum of
  squared second differences of log-selectivity across bins, within each
  year, normalized by the number of bins. Penalizes jagged (non-smooth)
  selectivity-at-age curves. `0` (default) disables this term. Requires
  at least 3 bins to have any effect.

- wt_bin_diff:

  Non-negative weight on the unconditional bin first-difference penalty:
  the sum of squared first differences of log-selectivity across bins,
  within each year, normalized by the number of bins. Unlike `wt_dome`
  (which only penalizes decreases), both increases and decreases
  contribute. Requires at least 2 bins to have any effect.

- wt_yr_diff:

  Non-negative weight on the inter-annual first-difference penalty: the
  sum of squared first differences of log-selectivity across years,
  within each bin, normalized by the number of years. Penalizes abrupt
  year-to-year jumps in selectivity-at-bin. `0` (default) disables this
  term. Requires at least 2 years to have any effect.

- wt_yr_curve:

  Non-negative weight on the inter-annual second-difference (smoothness)
  penalty: the sum of squared second differences of log-selectivity
  across years, within each bin, normalized by the number of years.
  Penalizes jagged (non-smooth) year-to-year selectivity trajectories.
  `0` (default) disables this term. Requires at least 3 years to have
  any effect.

- wt_dome:

  Non-negative weight on the dome-shape (non-monotonicity) penalty: for
  each year, penalizes any decrease in log-selectivity moving from one
  bin to the next (i.e. discourages, but does not forbid, dome shaped
  dynamics. `0` (default) disables this term.

- wt_mean_center:

  Non-negative weight on a per-year mean-centering (sum-to-zero)
  regularization: for each year, penalizes the squared mean of
  log-selectivity across bins. `0` (default) disables this term; set to
  `10000`.

- normalize:

  Logical. If `TRUE` (default), `wt_bin_curve` is divided by the number
  of bins the penalties act over and `wt_yr_diff`/`wt_yr_curve` are
  divided by the number of years. `SPoRC_rtmb.R` always calls this with
  `normalize = TRUE`.

- bin_range:

  Length-two vector giving the first and last bin the penalties act
  over, or `NULL` (default) for every bin. Restricting the range is how
  a shape penalty is confined to the older ages where a curve is
  expected to flatten, without constraining the ascending limb.

## Value

Numeric scalar: the positive log-likelihood contribution from the
requested penalty terms. Negated externally to form the negative
log-likelihood.

## Details

Every `wt_` argument accepts either a single number applied to all
years, or a vector with one value per year. A per-year vector lets a
penalty act only in the years where selectivity is allowed to change, or
act with a different strength in each year, which is how a random walk
with a year-specific standard deviation is expressed: set the year's
weight to `1 / (2 * sigma^2)` and pass `normalize = FALSE`. Years whose
weight is zero are skipped entirely.
