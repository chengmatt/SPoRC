# Build an unstructured correlation matrix from unconstrained parameters

An unstructured correlation across ages is also somtimes desirable as a
check. It places no shape on how ages covary, at the cost of
\\n(n-1)/2\\ parameters.

## Usage

``` r
build_us_corr(pars, n)
```

## Arguments

- pars:

  Numeric vector of length \\n(n-1)/2\\, unconstrained.

- n:

  Number of ages.

## Value

An \\n \times n\\ correlation matrix.

## Details

The parameters fill the strict lower triangle of a matrix whose diagonal
is one. Normalizing each row to unit length makes that matrix a Cholesky
factor, so the product with its transpose is a correlation matrix with
ones on the diagonal and is positive definite for any parameter values.
