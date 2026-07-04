# Construct an AR(1) correlation matrix

Builds an \\n \times n\\ correlation matrix whose \\(i,j)\\ element
equals \\\rho^{\|i-j\|}\\, corresponding to a stationary AR(1) process
with autocorrelation parameter \\\rho\\.

## Usage

``` r
get_AR1_CorrMat(n, rho)
```

## Arguments

- n:

  Integer. Matrix dimension (number of bins, ages, or lengths).

- rho:

  Numeric. AR(1) autocorrelation parameter in \\(-1, 1)\\. Values close
  to 1 produce strong positive correlation between adjacent bins; values
  close to 0 approach the identity matrix.

## Value

Numeric \\n \times n\\ correlation matrix.

## Examples

``` r
if (FALSE) { # \dontrun{
get_AR1_CorrMat(10, 0.5)
} # }
```
