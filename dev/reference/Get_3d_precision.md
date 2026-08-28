# Construct a 3D sparse precision matrix over ages, years, and cohorts

Builds a sparse \\(n\_{\text{ages}} \times n\_{\text{yrs}}) \times
(n\_{\text{ages}} \times n\_{\text{yrs}})\\ precision matrix \\Q\\ for a
Gaussian Markov random field (GMRF) with simultaneous autoregressive
(SAR) structure across three biological dimensions: age, year, and
cohort (age-year diagonal). The matrix is constructed via the
path-matrix factorization \\Q = (I - B)^\top \Omega^{-1} (I - B)\\,
where \\B\\ encodes the partial correlations and \\\Omega\\ is a
diagonal variance matrix. Two variance parameterizations are supported:
marginal (stationary) and conditional (non-stationary).

## Usage

``` r
Get_3d_precision(
  n_ages,
  n_yrs,
  pcorr_age,
  pcorr_year,
  pcorr_cohort,
  ln_var_value,
  Var_Type
)
```

## Arguments

- n_ages:

  Integer. Number of age classes.

- n_yrs:

  Integer. Number of years.

- pcorr_age:

  Numeric. Partial correlation along the age dimension (i.e., between
  adjacent ages within the same year).

- pcorr_year:

  Numeric. Partial correlation along the year dimension (i.e., between
  adjacent years within the same age).

- pcorr_cohort:

  Numeric. Partial correlation along the cohort diagonal (i.e., between
  the \\(a-1, y-1)\\ and \\(a, y)\\ cell).

- ln_var_value:

  Numeric. Log of the target variance. Exponentiated internally to
  \\\sigma^2 = \exp(\text{ln\\var\\value})\\.

- Var_Type:

  Integer. Variance parameterization: `0` = marginal (stationary)
  variance, where diagonal elements of \\\Omega\\ are solved recursively
  via the accumulator \\(I - B)^{-1}\\ to achieve a constant marginal
  variance \\\sigma^2\\ at every node (slower); `1` = conditional
  (non-stationary) variance, where all diagonal elements of \\\Omega\\
  are set to \\\sigma^2\\ directly (faster).

## Value

A sparse
[`Matrix::sparseMatrix`](https://rdrr.io/pkg/Matrix/man/sparseMatrix.html)
precision matrix \\Q\\ of dimension \\(n\_{\text{ages}} \times
n\_{\text{yrs}}) \times (n\_{\text{ages}} \times n\_{\text{yrs}})\\,
compatible with
[`RTMB::dgmrf`](https://rdrr.io/pkg/RTMB/man/MVgauss.html).
