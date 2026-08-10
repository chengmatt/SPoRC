# Derivative of equilibrium recruitment with respect to spawning biomass

Supplies the Jacobian term for the Newton solves in the spatial cases.

## Usage

``` r
equil_rec_ssb_deriv(S, S0, R0, h, rec_model)
```

## Arguments

- S:

  Equilibrium spawning biomass.

- S0:

  Unfished spawning biomass.

- R0:

  Unfished recruitment.

- h:

  Steepness.

- rec_model:

  1 for Beverton-Holt, 2 for Ricker.
