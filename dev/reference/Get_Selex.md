# Calculate Selectivity

Computes selectivity-at-bin using a suite of parametric,
semi-parametric, and non-parametric formulations. Supports constant,
time-varying, and fully flexible selectivity structures.

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
  n_yr_nodes_bicubic = NULL
)
```

## Arguments

- Selex_Model:

  Integer specifying the selectivity model:

  0

  :   Logistic (b50, slope): \\1 / (1 + \exp(-k(\text{bin} -
      b\_{50})))\\

  1

  :   Gamma-shaped dome (bin-at-peak \\b\_{\max}\\, curvature
      \\\delta\\).

  2

  :   Power function (monotonic decreasing): \\1 /
      \text{bin}^{\text{power}}\\.

  3

  :   Logistic (b50, b95 parameterization).

  4

  :   Double-normal dome with plateau and flexible tails (6 parameters).

  5

  :   Non-parametric selectivity: bin-level logit parameters mapped via
      `plogis`, optionally modified by time-varying deviations.

  6

  :   Logistic selectivity with asymptote: \\\alpha / (1 +
      \exp(-k(\text{bin} - b\_{50})))\\. Allows maximum selectivity
      \\\alpha \in (0,1)\\.

  7

  :   Logistic selectivity with asymptote (b50, b95 parameterization):
      \\\alpha / (1 + 19^{(b\_{50} - \text{bin})/b\_{95}})\\. Equivalent
      to Model 3 scaled by asymptote \\\alpha\\.

  8

  :   Bicubic spline over a bin-node x year-node grid (see
      `Wbin_bicubic`, `Wyr_bicubic`). One generalized form covers three
      cases depending on how the caller constructs the node grid and
      interpolation weights: a single smooth bin x year surface
      (bicubic), a time-invariant bin-only spline (`n_yr_nodes == 1`),
      or a bin-only spline re-fit independently within each of several
      year blocks (`n_yr_nodes == 1` within each of SPoRC's existing
      selectivity blocks). No `TimeVary_Model` deviation layering
      applies to this model.

- TimeVary_Model:

  Integer specifying temporal structure:

  0

  :   No time variation.

  1

  :   IID deviations applied multiplicatively to model parameters.

  2

  :   Random walk deviations applied multiplicatively to model
      parameters.

  3

  :   3D GMRF (marginal variance): deviations applied multiplicatively
      at bin level.

  4

  :   3D GMRF (conditional variance): deviations applied
      multiplicatively at bin level.

  5

  :   Separable 2D AR(1): deviations applied multiplicatively at bin
      level.

- pars:

  Numeric vector of log-scale selectivity parameters. Parameters are
  exponentiated or transformed depending on model specification:

  Model 0

  :   `c(ln_b50, ln_slope)`

  Model 1

  :   `c(ln_bmax, ln_delta)`

  Model 2

  :   `c(ln_power)`

  Model 3

  :   `c(ln_b50, ln_b95)`

  Model 4

  :   `c(p1, p2, p3, p4, p5, p6)`

  Model 5

  :   `c(logit_sel_1, ..., logit_sel_nbins)`

  Model 6

  :   `c(logit_alpha, ln_b50, ln_k)`

  Model 7

  :   `c(logit_alpha, ln_b50, ln_b95)`

  Model 8

  :   Flattened bin-node x year-node log-selectivity grid, length
      `ncol(Wyr_bicubic) * ncol(Wbin_bicubic)`, filled column-major into
      a `ncol(Wyr_bicubic)` (rows, year nodes) by `ncol(Wbin_bicubic)`
      (columns, bin nodes) matrix.

- ln_seldevs:

  Array of log-scale selectivity deviations with dimension
  `[n_regions, n_years, n_parameters_or_bins, n_sexes, 1]`.

  For `TimeVary_Model = 1-2`: deviations apply to selectivity parameters
  on the natural scale after exponentiation.

  For `TimeVary_Model = 3-5`: deviations apply multiplicatively at the
  bin level to the constructed selectivity curve.

  For `Selex_Model = 5`: deviations act directly on bin-level logit
  selectivity parameters prior to logistic transformation.

- Region:

  Integer region index.

- Year:

  Integer year index (absolute, i.e. a row index into `Wyr_bicubic`).
  Only used directly by `Selex_Model == 8`; otherwise only used to index
  `ln_seldevs`.

- Bin:

  Numeric vector of bins (ages or lengths).

- Sex:

  Integer sex index.

- Wbin_bicubic:

  Numeric `length(Bin) x n_bin_nodes` natural cubic spline weight matrix
  (see
  [`Get_Natural_Cubic_Spline_Weights`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Natural_Cubic_Spline_Weights.md)),
  mapping bin-node log-selectivity values onto `Bin`. Only used when
  `Selex_Model == 8`; ignored (may be `NULL`) otherwise. Zero padding in
  unused columns (e.g. when a shared parameter array is padded to a
  common width across fleets/blocks) contributes nothing, since it is
  multiplied through to zero.

- Wyr_bicubic:

  Numeric `n_yrs_total x n_yr_nodes` interpolation weight matrix mapping
  year-node log-selectivity values onto every absolute model year (rows
  beyond the fitted block are typically constructed to hold the boundary
  node constant). Row `Year` is used for this call. Only used when
  `Selex_Model == 8`; ignored (may be `NULL`) otherwise. A single-column
  matrix of all-1s (`n_yr_nodes == 1`) yields a time-invariant bin-only
  spline, since every year maps onto the same single node.

- n_bin_nodes_bicubic, n_yr_nodes_bicubic:

  Integer. This fleet/block's own true number of bin nodes / year nodes.
  Only used when `Selex_Model == 8`. **Must** be supplied whenever
  `Wbin_bicubic`/`Wyr_bicubic` may have been zero-padded wider than this
  specific block's own grid (e.g. because some *other* fleet/block
  shares the same padded storage array but uses a larger bicubic grid),
  `ncol(Wbin_bicubic)`/`ncol(Wyr_bicubic)` give the padded (shared)
  width, not this block's true node counts, and using the padded width
  to reshape `pars` would misassign which flattened parameter values
  land in which (bin-node, year-node) cell. Default `NULL` falls back to
  `ncol(Wbin_bicubic)`/`ncol(Wyr_bicubic)` which is correct whenever no
  padding mismatch is possible (a single bicubic block or fleet, or
  direct unit testing).

## Value

Numeric vector of selectivity values corresponding to `Bin`. Values are
on the natural scale and are not normalized unless specified in
downstream components.

## Details

For `TimeVary_Model = 0`, only the base parametric form is evaluated.

For `TimeVary_Model = 1-2`, deviations modify model parameters
multiplicatively on the natural scale after transformation.

For `TimeVary_Model = 3-5`, deviations act directly on the constructed
selectivity curve as multiplicative log-normal perturbations:
\$\$\text{selex} = \text{selex} \cdot \exp(\delta\_{r,y,b,s})\$\$.

For `Selex_Model = 5`, selectivity is fully non-parametric: bin-specific
logit parameters are optionally adjusted by time-varying deviations and
transformed via: \$\$\text{selex}\_b = \text{logit}^{-1}(\eta_b)\$\$.

Models 6 and 7 extend logistic selectivity by introducing an asymptote
parameter \\\alpha \in (0,1)\\ that allows selectivity to saturate below
full vulnerability.
