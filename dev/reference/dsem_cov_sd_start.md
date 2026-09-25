# Starting log sd of a covariate's observation error

Half the observed spread on each family's own scale: the sd for normal,
the CV for gamma, the sd of the log for lognormal. Zero (a dispersion of
one) for the tweedie and for the families without a spread parameter.

## Usage

``` r
dsem_cov_sd_start(y, family)
```

## Arguments

- y:

  Observations with NA for missing years.

- family:

  Family code.

## Value

Scalar log sd.
