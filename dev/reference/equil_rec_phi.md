# Equilibrium recruitment from spawning biomass per recruit

Solves R = f(R \* phi) for R, where f is the stock-recruit curve. Both
forms pass through (S0, R0) and diverge away from it, so a model fitted
with one and given reference points from the other is silently wrong
rather than noisy.

## Usage

``` r
equil_rec_phi(phi, phi0, R0, h, rec_model)
```

## Arguments

- phi:

  Spawning biomass per recruit under fishing.

- phi0:

  Unfished spawning biomass per recruit.

- R0:

  Unfished recruitment.

- h:

  Steepness.

- rec_model:

  1 for Beverton-Holt, 2 for Ricker.
