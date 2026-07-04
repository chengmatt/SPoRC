# Construct a constant (exchangeable) correlation matrix

Builds an \\n \times n\\ correlation matrix with 1 on the diagonal and
\\\rho\\ on all off-diagonal elements, corresponding to a compound
symmetry (exchangeable) covariance structure. Used in SPoRC to model
constant correlation across sexes in the 2D logistic-normal composition
likelihood (`comp_like = 4`).

## Usage

``` r
get_Constant_CorrMat(n, rho)
```

## Arguments

- n:

  Integer. Matrix dimension (typically `n_sexes`).

- rho:

  Numeric. Off-diagonal correlation in \\(-1, 1)\\. A value of 0
  produces the identity matrix; a value approaching 1 produces
  near-perfect correlation across sexes.

## Value

Numeric \\n \times n\\ correlation matrix.

## Examples

``` r
if (FALSE) { # \dontrun{
get_Constant_CorrMat(2, 0.5)
} # }
```
