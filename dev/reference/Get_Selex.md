# Calculate Selectivity

Selectivity-at-bin under any of the parametric, semi-parametric and
non-parametric forms, constant or time-varying. Equations are in the
model equations vignette.

## Usage

``` r
Get_Selex(
  Selex_Model,
  TimeVary_Model,
  pars,
  ln_seldevs,
  Region,
  Year,
  Bin,
  Sex,
  Wbin_bicubic = NULL,
  Wyr_bicubic = NULL,
  n_bin_nodes_bicubic = NULL,
  n_yr_nodes_bicubic = NULL,
  bin_devs = NULL,
  bin_dev_bins = NULL,
  sel_norm_bins = NULL,
  n_sel_bins = NULL,
  apical = 1,
  dbnrml_raw = c(0, 0),
  dbnrml_startbin = 1
)
```

## Arguments

- Selex_Model:

  Integer selectivity form. `0` logistic (b50, slope), `1` gamma dome
  (bin-at-peak, curvature), `2` power \\1/\text{bin}^{p}\\, `3` logistic
  (b50, b95), `4` double normal with plateau and flexible tails (six
  parameters, \\p_1\\ the bin the plateau starts at on the bin scale,
  \\p_5\\ and \\p_6\\ the selectivity at the first and last bins), `6`
  and `7` the two logistics scaled by an asymptote \\\alpha \in (0,1)\\,
  `5` non-parametric on the logit scale through `plogis`, standardized
  over years and bins jointly, `9` non-parametric on the log scale with
  the level free and each year centered over `sel_norm_bins`, and `8` a
  bicubic spline over a bin-node by year-node grid, which takes no
  `TimeVary_Model` deviations.

- TimeVary_Model:

  Integer temporal structure. `0` none, `1` iid and `2` random walk
  deviations on the model parameters, `3` and `4` a 3D GMRF on the
  marginal or conditional variance and `5` a separable 2D AR(1), all
  three applied at bin level.

- pars:

  Numeric vector of selectivity parameters on the transformed scale, in
  the order the form expects: `c(ln_b50, ln_slope)` for model 0,
  `c(ln_bmax, ln_delta)` for 1, `c(ln_power)` for 2, `c(ln_b50, ln_b95)`
  for 3, `c(p1, ..., p6)` for 4, one logit value per bin for 5,
  `c(logit_alpha, ln_b50, ln_k)` for 6, `c(logit_alpha, ln_b50, ln_b95)`
  for 7, one log value per bin for 9, and for 8 the flattened bin-node
  by year-node grid, filled column-major into a `n_yr_nodes` by
  `n_bin_nodes` matrix.

- ln_seldevs:

  Array of log-scale deviations
  `[n_regions, n_years, n_parameters_or_bins, n_sexes, 1]`. Under
  `TimeVary_Model` 1 and 2 they apply to the parameters after
  exponentiation; under 3 to 5 they apply multiplicatively at bin level
  to the constructed curve; under `Selex_Model = 5` they act on the
  bin-level logit parameters before the transformation.

- Region:

  Integer region index.

- Year:

  Integer absolute year index, a row index into `Wyr_bicubic`. Read
  directly by `Selex_Model == 8`, and otherwise only to index
  `ln_seldevs`.

- Bin:

  Numeric vector of bins (ages or lengths).

- Sex:

  Integer sex index.

- Wbin_bicubic:

  Numeric `length(Bin) x n_bin_nodes` natural cubic spline weight matrix
  (see
  [`Get_Natural_Cubic_Spline_Weights`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Natural_Cubic_Spline_Weights.md))
  mapping bin-node log-selectivity onto `Bin`. Read under
  `Selex_Model == 8` only, and may be `NULL` otherwise. Zero padding in
  unused columns contributes nothing.

- Wyr_bicubic:

  Numeric `n_yrs_total x n_yr_nodes` weight matrix mapping year-node
  log-selectivity onto every absolute model year; row `Year` is used for
  this call. Read under `Selex_Model == 8` only. A single column of ones
  gives a time-invariant bin-only spline.

- n_bin_nodes_bicubic, n_yr_nodes_bicubic:

  Integer, this fleet and block's own node counts. Supply them whenever
  `Wbin_bicubic` or `Wyr_bicubic` may have been zero-padded wider than
  this block's grid, because another block shares the storage array: the
  padded width would misassign which flattened parameter lands in which
  cell. `NULL` (default) falls back to
  [`ncol()`](https://rdrr.io/r/base/nrow.html), which is right when no
  padding is possible.

- bin_devs:

  Array of log-scale bin-override deviations
  `[n_regions, n_years, n_bins, n_sexes, 1]`, or `NULL`. Supplies the
  value for every bin named in `bin_dev_bins`.

- bin_dev_bins:

  Integer vector of the bins whose selectivity is replaced by
  `exp(bin_devs[...])` rather than taken from the functional form, or
  `NULL` for none. Applied after everything else, including any
  standardization, so the named bins are governed by their own
  deviations while the rest of the curve keeps its parametric shape.

- sel_norm_bins:

  Integer vector of the bins the mean-one standardization averages over
  under `Selex_Model = 9`, or `NULL` for every bin. A gear whose
  catchability is defined against part of the age range standardizes
  over that part, which shifts the scale absorbed by q.

- n_sel_bins:

  Integer, or `NULL`/`0` for none. Bins beyond this one are kept at its
  computed value rather than evaluated through the functional form, the
  `NSelBins` plateau convention several assessments apply. Applied after
  the form and its parameter deviations, but before the bin-level
  semi-parametric deviations and the bin overrides.

- apical:

  Numeric, the height the double normal's limbs are built up to and the
  plateau sits at. `1` (default) is the ordinary curve; a sex with an
  apical offset takes `exp(ln_*sel_sex_scale)` here, which moves the
  middle of its curve and leaves the first and last bins where their own
  parameters put them. Ignored by every other form.

- dbnrml_raw:

  Integer vector of length two (0/1) for the double normal: whether the
  ascending and descending limbs are left as raw Gaussians, \\\exp(-(x -
  peak)^2 / width)\\ built up to the apical value, rather than anchored
  to `p5` and `p6` at the first and last bins. Default `c(0, 0)`.

- dbnrml_startbin:

  Integer, the bin the double normal's ascending limb is anchored at and
  built up from. `1` (default) anchors at the first bin. Anchor later
  when the compositions start above the population's first length bin:
  `p5` is then the selectivity at that bin, and every bin below takes
  \\(b / b\_{start})^2\\ times the selectivity there.

## Value

Numeric vector of selectivity at `Bin`, on the natural scale and not
normalized unless a downstream component says so.
