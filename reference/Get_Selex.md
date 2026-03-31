# Calculate Selectivity

Computes selectivity using one of several parametric or semi-parametric
models. Supports both constant and time-varying selectivity, including
random effects and GMRF-based deviations.

## Usage

``` r
Get_Selex(
  Selex_Model,
  TimeVary_Model,
  ln_Pars,
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

  :   Logistic selectivity: uses b50 and slope parameters

  1

  :   Gamma-shaped (dome) selectivity: uses bin-at-peak and delta
      parameters

  2

  :   Power function selectivity: decreasing selectivity with bin

  3

  :   Logistic selectivity using b50 and b95

  4

  :   Double-normal (dome-shaped) selectivity with plateau and flexible
      tails

- TimeVary_Model:

  Integer specifying time variation structure:

  0

  :   No time variation (constant or blocked)

  1

  :   IID deviations

  2

  :   Random walk over time

  3

  :   3D AR1-GMRF marginal

  4

  :   3D AR1-GMRF conditional

  5

  :   2D AR1-GMRF

- ln_Pars:

  Vector of log-transformed selectivity parameters. Interpretation
  depends on \`Selex_Model\`.

- ln_seldevs:

  Array of selectivity deviations (may be log-scale), dimensioned as:
  \[n_regions, n_years, n_bins, n_sexes, 1\]. Used for time-varying or
  semi-parametric selectivity.

- Region:

  Integer index for region

- Year:

  Integer index for year

- Bin:

  Numeric vector of bins to compute selectivity for

- Sex:

  Integer index for sex

## Value

A numeric vector of selectivity values corresponding to the bins
specified in the model.

## Details

Selectivity parameters are transformed internally (typically using
\`exp()\` or logistic transformations) to ensure they remain in valid
ranges. Deviations (\`ln_seldevs\`) apply multiplicatively to these
transformed parameters when time-varying models are used. For
semi-parametric models (TimeVary_Model 3–5), deviations are applied
directly to the resulting selectivity curve.
