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
  Sex
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

- ln_seldevs:

  Array of log-scale selectivity deviations with dimension
  `[n_regions, n_years, n_parameters_or_bins, n_sexes, 1]`.

  For `TimeVary_Model = 1–2`: deviations apply to selectivity parameters
  on the natural scale after exponentiation.

  For `TimeVary_Model = 3–5`: deviations apply multiplicatively at the
  bin level to the constructed selectivity curve.

  For `Selex_Model = 5`: deviations act directly on bin-level logit
  selectivity parameters prior to logistic transformation.

- Region:

  Integer region index.

- Year:

  Integer year index.

- Bin:

  Numeric vector of bins (ages or lengths).

- Sex:

  Integer sex index.

## Value

Numeric vector of selectivity values corresponding to `Bin`. Values are
on the natural scale and are not normalized unless specified in
downstream components.

## Details

For `TimeVary_Model = 0`, only the base parametric form is evaluated.

For `TimeVary_Model = 1–2`, deviations modify model parameters
multiplicatively on the natural scale after transformation.

For `TimeVary_Model = 3–5`, deviations act directly on the constructed
selectivity curve as multiplicative log-normal perturbations:
\$\$\text{selex} = \text{selex} \cdot \exp(\delta\_{r,y,b,s})\$\$.

For `Selex_Model = 5`, selectivity is fully non-parametric: bin-specific
logit parameters are optionally adjusted by time-varying deviations and
transformed via: \$\$\text{selex}\_b = \text{logit}^{-1}(\eta_b)\$\$.

Models 6 and 7 extend logistic selectivity by introducing an asymptote
parameter \\\alpha \in (0,1)\\ that allows selectivity to saturate below
full vulnerability.
