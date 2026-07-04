# Squeeze a probability onto the open unit interval

Maps values in the closed interval \\\[0, 1\]\\ onto the open interval
to avoid boundary evaluations (e.g. `log(0)` or division by zero) inside
the conditional composition likelihoods. Identical to TMB's
`convenience.hpp` `squeeze()` (with `eps` equal to machine epsilon), so
this implementation agrees with WHAM's OSA composition likelihood to the
last digit.

## Usage

``` r
osa_squeeze(u)
```

## Arguments

- u:

  Numeric or AD scalar/vector of probabilities in \\\[0, 1\]\\.

## Value

The input mapped into the open interval \\(eps, 1 - eps)\\.
