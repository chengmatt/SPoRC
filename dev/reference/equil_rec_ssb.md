# Equilibrium recruitment from spawning biomass

The same two curves evaluated at a spawning biomass rather than a
per-recruit quantity, for the spatial solvers that iterate on
recruitment directly.

## Usage

``` r
equil_rec_ssb(S, S0, R0, h, rec_model)
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
