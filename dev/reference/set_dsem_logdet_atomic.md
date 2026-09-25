# Switch TMB's atomic sparse log determinant on or off inside RTMB

A moderated arrow puts a random effect in the precision, and TMB's
atomic sparse log determinant has no second derivative, so the plain
taped factorization has to be used instead. `Setup_Mod_DSEM` switches it
off for the session when the arrows need it.

## Usage

``` r
set_dsem_logdet_atomic(value)
```

## Arguments

- value:

  0 for the plain taped factorization, 1 for TMB's default.

## Value

The previous value, invisibly, or `NA` if RTMB does not expose the flag.
